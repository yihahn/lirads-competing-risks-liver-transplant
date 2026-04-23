library(dplyr)
library(prodlim)    # Used by pec
library(survival)   # Used by pec
library(riskRegression) # Used by pec
library(pec)        # Wolbers C-index
library(caret)      # createFolds
library(randomForestSRC) # For random forest competing risks models
library(readr)

set.seed(20250521)

############ MUST BE EDITED ########################################
out_base_name = "rsf_unified_cv_BE_result"
in_base_name = "../data/stratified_kfold_splits/fold"

# factor_vars 및 cont_vars 정의
factor_vars <- c('sex', 'diabetes', 'ckd_hd', # baseline
                 'hbv', 'hcv', 'image_lc', 'ald', 'masld', 'ldlt',  # cause of liver disease
                 'pre_tx', 'tace', 'rfa', 'resection', 'rtx',  # previous HCC treatment
		 'lr_5_ct', 'lr_tr_v_ct', 'lr_3_ct', 'lr_4_ct',
		 'lr_m_ct', 'lr_tr_nv_ct', 'lr_tr_e_ct'  # LIRADS count information
)
cont_vars <- c('log_lr_5_size_ct', 'log_lr_tr_v_size_whole_ct',
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

# 파라미터 그리드 정의
param_grid <- expand.grid(
    ntree = c(500, 1000),
    nodesize = c(5, 10, 15)
)

# function for valid factor variable selection (keep this as is)
check_n_filter_factor_vars <- function(data, var_names) {
    valid_factor_vars <- c()
    for (fv in var_names) {
        # 1. Check the variables exist in data
        if (fv %in% colnames(data)) {
            # 2. Check the variables is factor
            if (is.factor(data[[fv]])) {
                # 3. Check unique values of the variable
                num_levels <- length(levels(data[[fv]]))
                if (num_levels < 2) {
                    cat(paste0("Warning: '", fv,
                               "' is excluded due to single value.\n"))
                } else {
                    valid_factor_vars <- c(valid_factor_vars, fv)
                }
            } else {
                cat(paste0("Warning: '", fv,
                           "' is not factor, and excluded in RF-SRC.\n"))
            }
        } else {
            cat(paste0("Warning: '", fv, "' does not exist in data. and excluded in RF-SRC.\n"))
        }
    }
    return(valid_factor_vars)
}


# 전체 데이터 로드와 전처리 함수
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


# 5 fold stratification by fstatus 
create_cv_folds <- function(data, k=5) {
    set.seed(20250521)
    folds <- createFolds(data$fstatus, k=k, list=TRUE)
    return(folds)
}


# RSF 모델 CV 평가 함수
evaluate_rsf_cv <- function(data, folds, variables, params) {
    c_indices <- numeric(length(folds))
    
    formula_str <- paste("Surv(ftime, fstatus) ~", paste(variables, collapse=" + "))
    
    for(i in seq_along(folds)) {
        test_idx <- folds[[i]]
        train_data <- data[-test_idx, ]
        test_data <- data[test_idx, ]
        
        # Select only required columns for RF
        train_data_for_rf <- train_data %>%
            select(ftime, fstatus, all_of(variables))
        test_data_for_rf <- test_data %>%
            select(ftime, fstatus, all_of(variables))
        
        # Adjust mtry if necessary
        current_mtry <- min(params$mtry, length(variables))
        
        # 모델 학습
        suppressWarnings({
            rsf_model <- rfsrc(as.formula(formula_str),
                              data=train_data_for_rf,
                              ntree=params$ntree,
                              nodesize=params$nodesize,
                              mtry=current_mtry,
                              splitrule="logrankCR",
                              na.action="na.omit",
                              importance=FALSE)
        })
        
        # 예측 및 성능 평가 (Wolbers C-index)
        cindex_result <- suppressMessages(pec::cindex(list(rf_model=rsf_model), 
                                                      formula=as.formula(formula_str),
                                                      data=test_data_for_rf, 
                                                      cause=1,
                                                      verbose=FALSE))
        c_indices[i] <- cindex_result$AppCindex$rf_model
        
    }
    mean_c_index = mean(c_indices, na.rm=TRUE)
    return(list(
        mean_c_index = mean_c_index,
        c_indices = c_indices
    ))
}

# Initial baseline performance 계산 함수 (새로 추가)
calculate_initial_performance <- function(data, param_grid, available_vars, k=5) {
    cat("\n=== Calculating Initial Baseline Performance ===\n")
    cat("Using all", length(available_vars), "available variables\n")
    cat("Variables:", paste(available_vars, collapse=", "), "\n")
    
    folds <- create_cv_folds(data, k)
    
    best_mean_c_index <- 0
    best_params <- NULL
    best_c_indices <- NULL
    
    # 각 파라미터 조합으로 초기 성능 계산
    cat("\nTesting parameter combinations for initial baseline:\n")
    for(j in 1:nrow(param_grid)) {
        params <- param_grid[j,]
        cat(sprintf("  Parameters: %s\n", paste(names(params), params, sep="=", collapse=", ")))
        
        result <- evaluate_rsf_cv(data, folds, available_vars, params)
        
        cat(sprintf("    Mean C-index: %.4f (SD: %.4f)\n", 
                   result$mean_c_index, sd(result$c_indices, na.rm=TRUE)))
        
        if(result$mean_c_index > best_mean_c_index) {
            best_mean_c_index <- result$mean_c_index
            best_params <- params
            best_c_indices <- result$c_indices
        }
    }
    
    cat("\n--- Initial Baseline Results ---\n")
    cat(sprintf("Best initial mean C-index: %.4f\n", best_mean_c_index))
    cat(sprintf("Best initial CV SD: %.4f\n", sd(best_c_indices, na.rm=TRUE)))
    cat(sprintf("Best initial parameters: %s\n", paste(names(best_params), best_params, sep="=", collapse=", ")))
    
    return(list(
        mean_c_index = best_mean_c_index,
        c_indices = best_c_indices,
        params = best_params,
        variables = available_vars
    ))
}

# Backward selection with grid search and CV (수정됨)
backward_selection_cv <- function(data, param_grid, k=5) {
    # 유효한 변수들만 선택
    val_factor_vars <- check_n_filter_factor_vars(data, factor_vars)
    all_candidate_vars <- c(cont_vars, val_factor_vars)
    available_vars <- all_candidate_vars[all_candidate_vars %in% colnames(data)]

    cat("length(available_vars) =", length(available_vars), 
        "length(all_candidate_vars) =", length(all_candidate_vars), "\n")
    
    if (length(available_vars) == 0) {
        stop("No valid variables found for analysis.")
    }
    
    cat("Starting analysis with", length(available_vars), "variables:\n")
    cat(paste(available_vars, collapse=", "), "\n")
    
    # 초기 baseline 성능 계산 (새로 추가)
    initial_performance <- calculate_initial_performance(data, param_grid, available_vars, k)
    
    folds <- create_cv_folds(data, k)
    best_results <- list()
    current_vars <- available_vars
    consecutive_decreases <- 0
    previous_best_c_index <- initial_performance$mean_c_index  # 초기값을 baseline으로 설정
    
    overall_best_result <- list(
        mean_c_index = initial_performance$mean_c_index,
        variables = initial_performance$variables,
        params = initial_performance$params,
        c_indices = initial_performance$c_indices,
        iteration = 0  # iteration 0 = initial baseline
    )
    
    # Iteration 시작 전 초기 상태 출력
    cat("\n=== Starting Backward Elimination ===\n")
    cat(sprintf("Initial baseline C-index to beat: %.4f\n", previous_best_c_index))
    
    for(i in 1:length(available_vars)) {
        if (length(current_vars) <= 1) {
            cat("\nStopping: Only one variable remaining\n")
            break
        }
        
        cat(sprintf("\n=== Iteration %d with %d variables ===\n", i, length(current_vars)))
        
        best_mean_c_index <- 0
        best_params <- NULL
        best_vars_to_remove <- NULL
        best_c_indices <- NULL
        
        # 각 파라미터 조합 시도
        for(j in 1:nrow(param_grid)) {
            params <- param_grid[j,]
            cat(sprintf("\nTrying parameters: %s\n", paste(names(params), params, sep="=", collapse=", ")))
            
            # 각 변수를 하나씩 제거해보며 테스트
            cnt = 0
            for(var in current_vars) {
                test_vars <- current_vars[current_vars != var]
                cnt = cnt + 1
                
                if (length(test_vars) == 0) next
                
                result <- evaluate_rsf_cv(data, folds, test_vars, params)
                
                cat(sprintf("  Testing removal of '%s', mean C index: %.4f\n", var, result$mean_c_index))
                
                if(result$mean_c_index > best_mean_c_index) {
                    best_mean_c_index <- result$mean_c_index
                    best_params <- params
                    best_vars_to_remove <- var
                    best_c_indices <- result$c_indices
                }
                
                if(result$mean_c_index > overall_best_result$mean_c_index) {
                    overall_best_result$mean_c_index <- result$mean_c_index
                    overall_best_result$variables <- test_vars
                    overall_best_result$params <- params
                    overall_best_result$c_indices <- result$c_indices
                    overall_best_result$iteration <- i
                }
            }
        }
        
        # C-index 개선 여부 확인
        if(best_mean_c_index <= previous_best_c_index) {
            consecutive_decreases <- consecutive_decreases + 1
            cat(sprintf("\nPerformance did not improve (%.4f <= %.4f). Consecutive decreases: %d\n",
                       best_mean_c_index, previous_best_c_index, consecutive_decreases))
        } else {
            consecutive_decreases <- 0
            cat(sprintf("\nPerformance improved! (%.4f > %.4f)\n", 
                       best_mean_c_index, previous_best_c_index))
        }
        
        # 2번 연속 개선이 없으면 종료
        if(consecutive_decreases >= 2) {
            cat("\nStopping early due to consecutive lack of improvement\n")
            break
        }
        
        # 결과 저장
        best_results[[i]] <- list(
            variables = current_vars,
            removed_var = best_vars_to_remove,
            mean_c_index = best_mean_c_index,
            c_indices = best_c_indices,
            params = best_params
        )
        
        previous_best_c_index <- best_mean_c_index
        current_vars <- current_vars[current_vars != best_vars_to_remove]
        
        # 결과 출력
        cat(sprintf("\n--- Iteration %d Results ---\n", i))
        cat(sprintf("Best mean C-index: %.4f\n", best_mean_c_index))
        cat(sprintf("CV Standard Deviation: %.4f\n", sd(best_c_indices, na.rm=TRUE)))
        cat(sprintf("Removed variable: %s\n", best_vars_to_remove))
        cat(sprintf("Best parameters: %s\n", paste(names(best_params), best_params, sep="=", collapse=", ")))
        cat(sprintf("Remaining variables (%d): %s\n", length(current_vars), paste(current_vars, collapse=", ")))
    }
    
    # 최종 결과 출력
    cat("\n\nFINAL BEST OVERALL PERFORMANCE:\n")
    cat(sprintf("Mean C-index: %.4f\n", overall_best_result$mean_c_index))
    cat(sprintf("CV Standard Deviation: %.4f\n", sd(overall_best_result$c_indices, na.rm=TRUE)))
    cat(sprintf("Found at iteration: %d", overall_best_result$iteration))
    if(overall_best_result$iteration == 0) {
        cat(" (Initial baseline)\n")
    } else {
        cat("\n")
    }
    cat(sprintf("Number of variables: %d\n", length(overall_best_result$variables)))
    cat(sprintf("Variables: %s\n", paste(overall_best_result$variables, collapse=", ")))
    cat("Parameters:\n")
    print(overall_best_result$params)
    
    return(list(
        initial_performance = initial_performance,  # 초기 성능 결과 추가
        results = best_results,
        overall_best = overall_best_result
    ))
}

# 결과 분석 및 저장 함수 (수정됨)
analyze_and_save_results <- function(results, out_base_name) {
    current_time <- format(Sys.time(), "%Y%m%d_%H%M")
    base_filename <- paste(out_base_name, current_time, sep="_")
    
    # 텍스트 결과 파일
    txt_output_file <- paste0(base_filename, ".txt")
    
    sink(txt_output_file, append=FALSE)
    
    cat("=== RSF Unified Cross-Validation Results ===\n\n")
    
    # Initial baseline results (새로 추가)
    initial <- results$initial_performance
    cat("INITIAL BASELINE (ALL VARIABLES):\n")
    cat(sprintf("Mean Wolbers C-index: %.4f (SD = %.4f)\n",
               initial$mean_c_index, sd(initial$c_indices, na.rm=TRUE)))
    cat(sprintf("Number of variables: %d\n", length(initial$variables)))
    cat(sprintf("Variables: %s\n", paste(initial$variables, collapse=", ")))
    cat("Initial parameters:\n")
    for (param_name in names(initial$params)) {
        cat(sprintf("  %s: %s\n", param_name, initial$params[[param_name]]))
    }
    
    # Overall best results
    overall_best <- results$overall_best
    cat("\n\nOVERALL BEST MODEL:\n")
    cat(sprintf("Mean Wolbers C-index: %.4f (SD = %.4f)\n",
               overall_best$mean_c_index, sd(overall_best$c_indices, na.rm=TRUE)))
    cat(sprintf("Selected at iteration: %d", overall_best$iteration))
    if(overall_best$iteration == 0) {
        cat(" (Initial baseline - no improvement found)\n")
    } else {
        cat("\n")
    }
    cat(sprintf("Number of variables: %d\n", length(overall_best$variables)))
    cat(sprintf("Selected variables: %s\n", paste(overall_best$variables, collapse=", ")))
    cat("\nOptimal parameters:\n")
    for (param_name in names(overall_best$params)) {
        cat(sprintf("  %s: %s\n", param_name, overall_best$params[[param_name]]))
    }
    
    # Iteration-by-iteration results
    if (length(results$results) > 0) {
        cat("\n=== ITERATION DETAILS ===\n")
        for (i in 1:length(results$results)) {
            iter_result <- results$results[[i]]
            cat(sprintf("\nIteration %d:\n", i))
            cat(sprintf("  C-index: %.4f (SD: %.4f)\n", 
                       iter_result$mean_c_index, sd(iter_result$c_indices, na.rm=TRUE)))
            cat(sprintf("  Removed variable: %s\n", iter_result$removed_var))
            cat(sprintf("  Remaining variables (%d): %s\n", 
                       length(iter_result$variables) - 1, 
                       paste(setdiff(iter_result$variables, iter_result$removed_var), collapse=", ")))
            cat(sprintf("  Best parameters: %s\n", 
                       paste(names(iter_result$params), iter_result$params, sep="=", collapse=", ")))
        }
    }
    
    sink()
    
    # RDS 파일로 전체 결과 저장
    rds_file <- paste0(base_filename, ".rds")
    saveRDS(results, rds_file)
    
    cat(sprintf("All console output saved to: %s\n", txt_output_file))
    cat(sprintf("Complete results saved to: %s\n", rds_file))
    
    return(list(
        txt_file = txt_output_file,
        rds_file = rds_file
    ))
}

# =============================================================================
# 메인 실행 코드
# =============================================================================

cat("\n=== Starting Unified RSF Cross-Validation Analysis ===\n")

# 1. 데이터 로드 및 전처리
cat("Loading and preprocessing data...\n")
data <- load_data()
cat(sprintf("Loaded data with %d observations and %d variables\n", nrow(data), ncol(data)))

# 2. Backward selection with CV (초기 baseline 계산 포함)
cat("\nStarting backward selection with cross-validation...\n")
results <- backward_selection_cv(data, param_grid, k=5)

# 3. 결과 분석 및 저장
cat("\nAnalyzing and saving results...\n")
output_files <- analyze_and_save_results(results, out_base_name)

cat("\n=== Analysis Complete ===\n")
cat("Output files created:\n")
for (file in output_files) {
    if (!is.null(file)) cat(sprintf("  - %s\n", file))
}
