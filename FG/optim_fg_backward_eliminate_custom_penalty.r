library(dplyr)
library(prodlim)    # Used by pec
library(survival)   # Used by pec
library(riskRegression) # Used by pec
library(pec)        # Wolbers C-index
library(caret)      # createFolds
library(readr)

set.seed(20250521)

############ MUST BE EDITED ########################################
out_base_name = "fg_penalized_BE_result"

# factor_vars 및 cont_vars 정의
factor_vars <- c('sex', 'diabetes', 'ckd_hd', # baseline
                 'hbv', 'hcv', 'image_lc', 'ald', 'masld', 'ldlt',  # cause of liver disease
                 'pre_tx', 'tace', 'rfa', 'resection', 'rtx'  # previous HCC treatment
)
cont_vars <- c('lr_5_ct', 'lr_tr_v_ct', 'lr_3_ct', 'lr_4_ct',
	       'lr_m_ct', 'lr_tr_nv_ct', 'lr_tr_e_ct',  # LIRADS count information
	       'log_lr_5_size_ct', 'log_lr_tr_v_size_whole_ct',
               'log_lr_tr_v_size_enhancing_ct', 'log_lr_3_size_ct',
               'log_lr_4_size_ct', 'log_lr_m_size_ct', 'log_lr_tr_nv_size_ct',
               'log_lr_tr_e_size_ct', # CT LIRADS size information
               'log_afp', 'log_pivka', # Tumor marker
               'alb', 'log_tbil', 'cr', 'log_pt_inr', 'na', 'log_crp', # lab test
               'age_lt', 'bmi' # demographic
)

# variables to be transformed log scale
vars_to_log <- c('afp', 'pivka', 'tbil', 'pt_inr', 'crp',
		 'lr_5_size_ct', 'lr_tr_v_size_whole_ct',
		 'lr_tr_v_size_enhancing_ct', 'lr_3_size_ct',
		 'lr_4_size_ct', 'lr_m_size_ct', 'lr_tr_nv_size_ct',
		 'lr_tr_e_size_ct'
)

########################################################################

# function for valid factor variable selection (기존 유지)
check_n_filter_factor_vars <- function(data, var_names) {
    valid_factor_vars <- c()
    for (fv in var_names) {
        if (fv %in% colnames(data)) {
            if (is.factor(data[[fv]])) {
                num_levels <- length(levels(data[[fv]]))
                if (num_levels < 2) {
                    cat(paste0("Warning: '", fv, "' is excluded due to single value.\n"))
                } else {
                    valid_factor_vars <- c(valid_factor_vars, fv)
                }
            } else {
                cat(paste0("Warning: '", fv, "' is not factor, and excluded.\n"))
            }
        } else {
            cat(paste0("Warning: '", fv, "' does not exist in data.\n"))
        }
    }
    return(valid_factor_vars)
}

# 데이터 로드 함수 (기존 유지)
load_data <- function() {
    data <- read.csv("../../data/internal_competing_risk_dataset.csv")
    
    # Apply log transformation
    for (var_nm in vars_to_log) {
        new_column_nm <- paste0("log_", var_nm)
        data[[new_column_nm]] <- log(data[[var_nm]] + 1)
    }
    
    # Convert to factors
    data[factor_vars] <- lapply(data[factor_vars], factor)    
    return(data)
}

# 5 fold stratification by fstatus (기존 유지)
create_cv_folds <- function(data, k=5) {
    set.seed(20250521)
    folds <- createFolds(data$fstatus, k=k, list=TRUE)
    return(folds)
}

# =============================================================================
# Penalty 계산 함수들
# =============================================================================

# 다양한 penalty 방식 계산
calculate_penalized_score <- function(c_index, sd_c_index, n_vars, penalty_method = "moderate") {

	# 여러 penalty 방법 중 선택
	penalty_configs <- list(
        # 보수적 (변수 개수에 민감)
        conservative = list(
            var_penalty = 0.008,     # 변수 1개당 0.008 감점
            sd_penalty = 0.5,        # SD가 클수록 더 큰 penalty
            min_improvement = 0.005   # 최소 개선 임계값
        ),
        
        # 중간 (균형잡힌 접근)
        moderate = list(
            var_penalty = 0.005,     # 변수 1개당 0.005 감점  
            sd_penalty = 0.3,        # SD penalty
            min_improvement = 0.003  # 최소 개선 임계값
        ),
        
        # 관대 (성능 우선)
        liberal = list(
            var_penalty = 0.003,     # 변수 1개당 0.003 감점
            sd_penalty = 0.1,        # 낮은 SD penalty
            min_improvement = 0.001  # 낮은 최소 개선 임계값
        )
    )
    
    config <- penalty_configs[[penalty_method]]
    
    # Penalty 계산
    var_penalty <- config$var_penalty * n_vars
    sd_penalty <- config$sd_penalty * sd_c_index
    total_penalty <- var_penalty + sd_penalty
    
    # Penalized score
    penalized_score <- c_index - total_penalty
    
    return(list(
        penalized_score = penalized_score,
        var_penalty = var_penalty,
        sd_penalty = sd_penalty,
        total_penalty = total_penalty,
        min_improvement_threshold = config$min_improvement
    ))
}

# 개선이 의미있는지 판단
is_meaningful_improvement <- function(old_c_index, new_c_index, old_sd, new_sd, 
                                    old_n_vars, new_n_vars, penalty_method = "moderate") {
    
    # 기존 모델과 새 모델의 penalized score
    old_penalty_info <- calculate_penalized_score(old_c_index, old_sd, old_n_vars, penalty_method)
    new_penalty_info <- calculate_penalized_score(new_c_index, new_sd, new_n_vars, penalty_method)
    
    old_score <- old_penalty_info$penalized_score
    new_score <- new_penalty_info$penalized_score
    
    improvement <- new_score - old_score
    min_threshold <- old_penalty_info$min_improvement_threshold
    
    is_meaningful <- improvement > min_threshold
    
    return(list(
        is_meaningful = is_meaningful,
        improvement = improvement,
        old_penalized = old_score,
        new_penalized = new_score,
        threshold = min_threshold,
        old_penalty_breakdown = old_penalty_info,
        new_penalty_breakdown = new_penalty_info
    ))
}

# =============================================================================
# Enhanced CV 평가 (기존 코드 간소화 버전)
# =============================================================================

evaluate_fg_cv <- function(data, folds, variables) {
    c_indices <- numeric(length(folds))
    formula_str <- paste("Hist(ftime, fstatus) ~", paste(variables, collapse=" + "))
    
    for(i in seq_along(folds)) {
        test_idx <- folds[[i]]
        train_data <- data[-test_idx, ]
        test_data <- data[test_idx, ]
        
        # Select only required columns
        train_data_for_fg <- train_data %>%
            select(ftime, fstatus, all_of(variables)) %>%
            filter(complete.cases(.))
            
        test_data_for_fg <- test_data %>%
            select(ftime, fstatus, all_of(variables)) %>%
            filter(complete.cases(.))
        
        tryCatch({
            fit <- suppressWarnings(FGR(as.formula(formula_str), 
                                      data=train_data_for_fg, 
                                      cause=1))
            cindex <- suppressMessages(pec::cindex(list(fg=fit), 
                                                 formula=as.formula(formula_str), 
                                                 data=test_data_for_fg,
                                                 cause=1))
            c_indices[i] <- cindex$AppCindex$fg
        }, error=function(e) {
            cat(sprintf("Fold %d error: %s\n", i, e$message))
            c_indices[i] <- NA
        })
    }
    
    return(list(
        mean_c_index = mean(c_indices, na.rm=TRUE),
        sd_c_index = sd(c_indices, na.rm=TRUE),
        c_indices = c_indices
    ))
}

# =============================================================================
# Penalized Backward Elimination 메인 함수
# =============================================================================

fg_penalized_backward_elimination <- function(data, k=5, penalty_method="moderate", max_iterations=50) {
    
    cat("\n=== Penalized Backward Elimination for Fine-Gray ===\n")
    cat(sprintf("Penalty method: %s\n", penalty_method))
    
    # 유효한 변수 선택
    val_factor_vars <- check_n_filter_factor_vars(data, factor_vars)
    all_candidate_vars <- c(cont_vars, val_factor_vars)
    available_vars <- all_candidate_vars[all_candidate_vars %in% colnames(data)]

    if (length(available_vars) == 0) {
        stop("No valid variables found for analysis.")
    }

    cat(sprintf("Starting penalized backward elimination with %d variables:\n", length(available_vars)))
    cat(paste(available_vars, collapse=", "), "\n")

    # CV folds 생성
    folds <- create_cv_folds(data, k)

    # 초기 baseline 성능 계산
    cat("\n=== Initial Baseline Performance ===\n")
    initial_result <- evaluate_fg_cv(data, folds, available_vars)
    initial_penalty_info <- calculate_penalized_score(
        initial_result$mean_c_index, 
        initial_result$sd_c_index, 
        length(available_vars), 
        penalty_method
    )

    cat(sprintf("Raw C-index: %.4f (±%.4f)\n", 
               initial_result$mean_c_index, initial_result$sd_c_index))
    cat(sprintf("Penalized Score: %.4f (penalty=%.4f)\n", 
               initial_penalty_info$penalized_score, initial_penalty_info$total_penalty))

    # 최적 결과 추적용
    best_overall <- list(
        mean_c_index = initial_result$mean_c_index,
        sd_c_index = initial_result$sd_c_index,
        penalized_score = initial_penalty_info$penalized_score,
        variables = available_vars,
        iteration = 0
    )

    # Backward elimination 시작
    current_vars <- available_vars
    current_result <- initial_result
    current_penalty_info <- initial_penalty_info
    iteration_results <- list()
    consecutive_no_improvement <- 0

    for(i in 1:max_iterations) {
        if (length(current_vars) <= 3) {  # 최소 3개 변수는 유지
            cat("\nStopping: Minimum variable count reached (3)\n")
            break
        }

        cat(sprintf("\n=== Iteration %d: Testing removal from %d variables ===\n",
                    i, length(current_vars)))
        
        cat(sprintf("Current penalized score: %.4f (raw C=%.4f, penalty=%.4f)\n",
                   current_penalty_info$penalized_score,
                   current_result$mean_c_index,
                   current_penalty_info$total_penalty))

        # 이번 iteration의 최고 성능 추적
        iteration_best <- list(
            penalized_score = -Inf,
            var_to_remove = NULL,
            remaining_vars = NULL,
            result = NULL,
            improvement_analysis = NULL
        )

        # 각 변수를 하나씩 제거해보며 테스트
        for(var_to_test in current_vars) {
            test_vars <- current_vars[current_vars != var_to_test]

            if (length(test_vars) == 0) next

            result <- evaluate_fg_cv(data, folds, test_vars)
            
            # 의미있는 개선인지 판단
            improvement_analysis <- is_meaningful_improvement(
                current_result$mean_c_index,
                result$mean_c_index,
                current_result$sd_c_index,
                result$sd_c_index,
                length(current_vars),
                length(test_vars),
                penalty_method
            )

            cat(sprintf("  Remove '%s': Raw C=%.4f, Penalized=%.4f, Improvement=%+.4f%s\n", 
                       var_to_test, 
                       result$mean_c_index,
                       improvement_analysis$new_penalized,
                       improvement_analysis$improvement,
                       ifelse(improvement_analysis$is_meaningful, " ✓", "")))

            # 이번 iteration 최고 성능 업데이트
            if(improvement_analysis$new_penalized > iteration_best$penalized_score) {
                iteration_best$penalized_score <- improvement_analysis$new_penalized
                iteration_best$var_to_remove <- var_to_test
                iteration_best$remaining_vars <- test_vars
                iteration_best$result <- result
                iteration_best$improvement_analysis <- improvement_analysis
            }

            # 전체 최고 성능 업데이트
            if(improvement_analysis$new_penalized > best_overall$penalized_score) {
                best_overall <- list(
                    mean_c_index = result$mean_c_index,
                    sd_c_index = result$sd_c_index,
                    penalized_score = improvement_analysis$new_penalized,
                    variables = test_vars,
                    iteration = i
                )
            }
        }

        # 의미있는 개선이 있는지 확인
        if(iteration_best$improvement_analysis$is_meaningful) {
            consecutive_no_improvement <- 0
            cat(sprintf("\n✓ Meaningful improvement found: %.4f > %.4f (threshold: %.4f)\n", 
                       iteration_best$improvement_analysis$improvement,
                       iteration_best$improvement_analysis$threshold,
                       iteration_best$improvement_analysis$threshold))
        } else {
            consecutive_no_improvement <- consecutive_no_improvement + 1
            cat(sprintf("\n✗ No meaningful improvement: %.4f <= %.4f (consecutive: %d)\n",
                       iteration_best$improvement_analysis$improvement,
                       iteration_best$improvement_analysis$threshold,
                       consecutive_no_improvement))
        }

        # 2번 연속 개선 없으면 조기 종료
        if(consecutive_no_improvement >= 2) {
            cat("\nEarly stopping: No meaningful improvement for 2 consecutive iterations\n")
            break
        }

        # 결과 저장 및 변수 업데이트 (의미있는 개선이 있을 때만)
        if(iteration_best$improvement_analysis$is_meaningful) {
            iteration_results[[i]] <- iteration_best
            current_vars <- iteration_best$remaining_vars
            current_result <- iteration_best$result
            current_penalty_info <- iteration_best$improvement_analysis$new_penalty_breakdown

            # Iteration 결과 출력
            cat(sprintf("\n--- Iteration %d Best ---\n", i))
            cat(sprintf("Raw C-index: %.4f (±%.4f)\n",
                       iteration_best$result$mean_c_index,
                       iteration_best$result$sd_c_index))
            cat(sprintf("Penalized Score: %.4f\n", iteration_best$penalized_score))
            cat(sprintf("Removed: %s\n", iteration_best$var_to_remove))
            cat(sprintf("Remaining (%d): %s\n",
                       length(iteration_best$remaining_vars),
                       paste(iteration_best$remaining_vars, collapse=", ")))
        } else {
            # 의미있는 개선이 없으면 반복 중단
            cat("No meaningful improvement found, stopping iterations\n")
            break
        }
    }

    # 최종 결과 출력
    cat("\n\nFINAL BEST RESULT (PENALIZED APPROACH):\n")
    cat(sprintf("Best Penalized Score: %.4f\n", best_overall$penalized_score))
    cat(sprintf("Raw C-index: %.4f (±%.4f)\n", 
               best_overall$mean_c_index, best_overall$sd_c_index))
    cat(sprintf("Found at iteration: %s\n",
               ifelse(best_overall$iteration == 0, "0 (baseline)", as.character(best_overall$iteration))))
    cat(sprintf("Variables (%d): %s\n",
               length(best_overall$variables),
               paste(best_overall$variables, collapse=", ")))

    # Penalty 분석
    final_penalty_info <- calculate_penalized_score(
        best_overall$mean_c_index,
        best_overall$sd_c_index,
        length(best_overall$variables),
        penalty_method
    )
    
    cat("\n--- Penalty Breakdown ---\n")
    cat(sprintf("Variable penalty: %.4f (%d vars × penalty rate)\n",
               final_penalty_info$var_penalty, length(best_overall$variables)))
    cat(sprintf("Variance penalty: %.4f (%.4f SD × penalty rate)\n",
               final_penalty_info$sd_penalty, best_overall$sd_c_index))
    cat(sprintf("Total penalty: %.4f\n", final_penalty_info$total_penalty))

    return(list(
        initial_baseline = list(
            mean_c_index = initial_result$mean_c_index,
            sd_c_index = initial_result$sd_c_index,
            variables = available_vars,
            penalized_score = initial_penalty_info$penalized_score
        ),
        iteration_results = iteration_results,
        best_overall = best_overall,
        penalty_method = penalty_method,
        final_penalty_info = final_penalty_info
    ))
}

# =============================================================================
# 결과 분석 및 저장 함수 (Enhanced)
# =============================================================================

analyze_and_save_results <- function(results, out_base_name) {
    current_time <- format(Sys.time(), "%Y%m%d_%H%M")
    base_filename <- paste(out_base_name, current_time, sep="_")
    
    # 텍스트 결과 파일
    txt_output_file <- paste0(base_filename, ".txt")
    
    sink(txt_output_file, append=FALSE)
    
    cat("=== Fine-Gray Penalized Backward Elimination Results ===\n\n")
    
    # Initial baseline
    initial <- results$initial_baseline
    cat("INITIAL BASELINE:\n")
    cat(sprintf("Raw C-index: %.4f (±%.4f)\n",
               initial$mean_c_index, initial$sd_c_index))
    cat(sprintf("Penalized Score: %.4f\n", initial$penalized_score))
    cat(sprintf("Variables (%d): %s\n\n", length(initial$variables),
               paste(initial$variables, collapse=", ")))
    
    # Final best result
    best <- results$best_overall
    cat("FINAL BEST MODEL (PENALIZED APPROACH):\n")
    cat(sprintf("Raw C-index: %.4f (±%.4f)\n", best$mean_c_index, best$sd_c_index))
    cat(sprintf("Penalized Score: %.4f\n", best$penalized_score))
    cat(sprintf("Found at iteration: %d\n", best$iteration))
    cat(sprintf("Variables (%d): %s\n\n", length(best$variables),
               paste(best$variables, collapse=", ")))
    
    # Penalty method and breakdown
    cat("PENALTY ANALYSIS:\n")
    cat(sprintf("Penalty method: %s\n", results$penalty_method))
    penalty_info <- results$final_penalty_info
    cat(sprintf("Variable penalty: %.4f\n", penalty_info$var_penalty))
    cat(sprintf("Variance penalty: %.4f\n", penalty_info$sd_penalty))
    cat(sprintf("Total penalty: %.4f\n\n", penalty_info$total_penalty))
    
    # Iteration details
    if (length(results$iteration_results) > 0) {
        cat("ITERATION DETAILS:\n")
        for (i in seq_along(results$iteration_results)) {
            iter_result <- results$iteration_results[[i]]
            cat(sprintf("Iteration %d:\n", i))
            cat(sprintf("  Raw C-index: %.4f (±%.4f)\n", 
                       iter_result$result$mean_c_index, iter_result$result$sd_c_index))
            cat(sprintf("  Penalized Score: %.4f\n", iter_result$penalized_score))
            cat(sprintf("  Removed variable: %s\n", iter_result$var_to_remove))
            cat(sprintf("  Remaining variables (%d): %s\n\n", 
                       length(iter_result$remaining_vars),
                       paste(iter_result$remaining_vars, collapse=", ")))
        }
    }
    
    # Comparison with non-penalized approach
    cat("COMPARISON WITH NON-PENALIZED APPROACH:\n")
    cat("Non-penalized would select based on raw C-index only.\n")
    cat("Penalized approach considers both performance and model complexity.\n")
    cat("This helps prevent overfitting and improves generalization.\n\n")
    
    # Recommendations
    cat("RECOMMENDATIONS:\n")
    cat("1. Validate selected model on external test set\n")
    cat("2. Consider clinical interpretability of selected variables\n")
    cat("3. Monitor for overfitting in deployment\n")
    cat("4. Consider ensemble methods for additional robustness\n")
    
    sink()
    
    # RDS 파일로 전체 결과 저장
    rds_file <- paste0(base_filename, ".rds")
    saveRDS(results, rds_file)
    
    cat(sprintf("Results saved to: %s\n", txt_output_file))
    cat(sprintf("Complete results saved to: %s\n", rds_file))
    
    return(list(
        txt_file = txt_output_file,
        rds_file = rds_file
    ))
}

# =============================================================================
# 메인 실행 코드
# =============================================================================

cat("\n=== Starting Penalized Fine-Gray Analysis ===\n")

# 1. 데이터 로드 및 전처리
cat("Loading and preprocessing data...\n")
data <- load_data()
cat(sprintf("Loaded data with %d observations and %d variables\n", nrow(data), ncol(data)))

# 2. Penalized backward elimination 실행
cat("\nStarting penalized backward selection...\n")
#results <- fg_penalized_backward_elimination(data, k=5, penalty_method="moderate")
#results <- fg_penalized_backward_elimination(data, k=5, penalty_method="conservative")
results <- fg_penalized_backward_elimination(data, k=5, penalty_method="liberal")

# 3. 결과 분석 및 저장
cat("\nAnalyzing and saving results...\n")
output_files <- analyze_and_save_results(results, out_base_name)

cat("\n=== Analysis Complete ===\n")
cat("Output files created:\n")
for (file in output_files) {
    if (!is.null(file)) cat(sprintf("  - %s\n", file))
}

cat("\n=== PENALTY METHOD EXPLANATION ===\n")
cat("The 'moderate' penalty method used:\n")
cat("- Variable penalty: 0.005 per variable\n")
cat("- Variance penalty: 0.3 × CV standard deviation\n")
cat("- Minimum improvement threshold: 0.003\n")
cat("\nThis prevents small, potentially random improvements\n")
cat("from adding complexity to the model.\n")

# 선택사항: 다른 penalty method들과 비교
cat("\nTo try different penalty methods, run:\n")
cat("results_conservative <- fg_penalized_backward_elimination(data, penalty_method='conservative')\n")
cat("results_liberal <- fg_penalized_backward_elimination(data, penalty_method='liberal')\n")
