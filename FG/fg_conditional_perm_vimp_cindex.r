library(survival)
library(riskRegression)
library(prodlim)
library(FNN)
library(dplyr)
#library(pheatmap)


# 결과 해석 함수들 - 모든 변수 포함, 필터링 없음
interpret_importance <- function(cpi_results) {
    
    # 1. 모든 변수의 중요도를 그대로 사용 (필터링 없음)
    cindex_importance <- cpi_results$drop_cindex
    
    # 2. 변수 분류 (유의성 기준)
    classify_variables <- function(cindex_imp, significant) {
        classification <- character(length(cindex_imp))
        names(classification) <- names(cindex_imp)
        
        for (i in seq_along(cindex_imp)) {
            if (significant[i]) {
                if (cindex_imp[i] > 0) {
                    classification[i] <- "Significantly Positive"
                } else {
                    classification[i] <- "Significantly Negative"
                }
            } else {
                classification[i] <- "Not Significant"
            }
        }
        return(classification)
    }
    
    var_classification <- classify_variables(cindex_importance, cpi_results$cindex_significant)
    
    # 3. 결과 정리
    result_df <- data.frame(
        Variable = names(cindex_importance),
        CIndex_Importance = round(cindex_importance, 6),
        P_value = round(cpi_results$cindex_pvalue, 4),
        Classification = var_classification,
        stringsAsFactors = FALSE
    )
    
    # 중요도 순으로 정렬 (내림차순)
    result_df <- result_df[order(result_df$CIndex_Importance, decreasing = TRUE), ]
    
    return(list(
        summary = result_df,
        all_cindex = cindex_importance,
        raw_results = cpi_results
    ))
}


# 시각화 함수 - 모든 변수 포함, 음수값도 표시
plot_importance_barplot <- function(interpretation_results,
                                   filename_prefix = "feature_importance",
                                   width = 12,
                                   height = 9) {
    
    # 모든 변수 사용 (중요도 순으로 이미 정렬됨)
    summary_df <- interpretation_results$summary
    
    if (nrow(summary_df) == 0) {
        cat("No variables to plot.\n")
        return()
    }
    
    # 바 플롯 생성 및 저장
    filename <- paste0(filename_prefix, "_cindex.pdf")
    pdf(filename, width = width, height = height)
    
    # 색상 설정 (양수는 파란색, 음수는 빨간색)
    colors <- ifelse(summary_df$CIndex_Importance > 0, "steelblue", "red")
    
    par(mar = c(8, 4, 4, 2))
    barplot(summary_df$CIndex_Importance,
            names.arg = summary_df$Variable,
            main = "Conditional Permutation Variable Importance",
            ylab = "Wolbers C-index Drop",
            las = 2,  # 수직으로 변수명 표시
            col = colors)
    
    # 0 기준선 추가
    abline(h = 0, lty = 2, col = "gray")
    
    dev.off()
    cat("C-index importance barplot saved as:", filename, "\n")
}


# 추천 변수 선택 함수 - 유의성 기준
recommend_features <- function(interpretation_results,
                              importance_threshold = 0.005,
                              require_significance = TRUE) {
    
    summary_df <- interpretation_results$summary
    
    if (require_significance) {
        # 통계적으로 유의하고 양의 중요도를 가진 변수만 선택
        recommended <- summary_df[summary_df$Classification == "Significantly Positive" &
                                 summary_df$CIndex_Importance >= importance_threshold, ]
    } else {
        # 임계값 이상의 양의 중요도를 가진 모든 변수 선택 (유의성 무관)
        recommended <- summary_df[summary_df$CIndex_Importance >= importance_threshold, ]
    }
    
    # C-index 중요도 기준으로 이미 정렬되어 있음
    
    cat("=== Recommended Features ===\n")
    print(recommended)
    cat("\n=== Feature Selection Summary ===\n")
    cat("Total variables tested:", nrow(summary_df), "\n")
    cat("Significantly positive:", sum(summary_df$Classification == "Significantly Positive"), "\n")
    cat("Significantly negative:", sum(summary_df$Classification == "Significantly Negative"), "\n")
    cat("Not significant:", sum(summary_df$Classification == "Not Significant"), "\n")
    cat("Recommended features:", nrow(recommended), "\n")
    
    return(recommended$Variable)
}


# Wolbers C-index 계산 함수
compute_wolbers_cindex <- function(model, test_data) {

    c_index <- pec::cindex(model, 
			       formula=formula(model$terms), 
			       data=test_data, 
			       cause=1,
			       verbose=FALSE)
	return(c_index$AppCindex$FGR)

}


# kNN 조건부 permutation으로 변수 하나를 교란한 data.frame을 생성
conditional_permute_knn <- function(df, target_col, cond_cols, k = 30, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    n <- nrow(df)
    
    # 수치형/범주형 분리
    is_num <- sapply(df[, cond_cols, drop = FALSE], is.numeric)
    Z_num <- NULL
    
    if (any(is_num)) {
        Z_num <- scale(as.matrix(df[, cond_cols[is_num], drop = FALSE]))
    }
    
    # 범주형 간단 처리: one-hot
    Z_cat <- NULL
    if (any(!is_num)) {
        Z_cat <- model.matrix(~ . - 1, data = df[, cond_cols[!is_num], drop = FALSE])
    }
    Z <- cbind(Z_num, Z_cat)
    
    if (is.null(Z) || ncol(Z) == 0) {
        # 조건집합이 비어있다면 전체에서 단순 무작위 치환(자기 자신 제외)
        idx_choices <- sapply(seq_len(n), function(i) sample(setdiff(seq_len(n), i), 1))
    } else {
        # kNN 이웃 탐색 (자기 자신 포함 k+1 반환 → 자기 자신 제외)
        k_use <- min(k + 1, n)
        nn <- get.knnx(Z, Z, k = k_use)$nn.index
        # 자기 자신 제외한 후보들 중 무작위 선택
        idx_choices <- integer(n)
        for (i in 1:n) {
            cand <- nn[i, ]
            cand <- cand[cand != i]
            if (length(cand) == 0) cand <- i
            idx_choices[i] <- sample(cand, 1)
        }
    }
    
    df_perm <- df
    df_perm[[target_col]] <- df[[target_col]][idx_choices]
    df_perm
}


# 메인 CPI 함수 - Wolbers C-index 버전
cpi_fg_wolbers <- function(fg,
                          data_val,
                          timeVar = "ftime",
                          statusVar = "fstatus",
                          features = NULL,
                          max_time = NULL,
                          B = 50,
                          k = 30,
                          alpha = 0.05,
                          seed = 20250521) {
    set.seed(seed)
    
    # 피처 목록 확정
    if (is.null(features)) {
        features <- setdiff(names(data_val), c(timeVar, statusVar))
    }
    
    
    baseline_cindex <- compute_wolbers_cindex(fg, data_val)
    cat("Baseline C-index:", round(baseline_cindex, 4), "\n")
    
    # 결과를 저장할 벡터 생성
    drop_cindex_vector <- numeric(length(features))
    cindex_pvalue_vector <- numeric(length(features))
    cindex_significant_vector <- logical(length(features))
    
    names(drop_cindex_vector) <- features
    names(cindex_pvalue_vector) <- features
    names(cindex_significant_vector) <- features
    
    for (j in seq_along(features)) {
        varj <- features[j]
        Sj   <- setdiff(features, varj)
        
        # ΔC = C0 - Cperm (클수록 중요)
        drops_c_vector <- numeric(B)
        
        for (b in 1:B) {
            df_perm <- data_val
            # 조건부(kNN)로 varj만 치환
            df_perm[, c(varj, Sj)] <- conditional_permute_knn(df = data_val[, c(varj, Sj), drop = FALSE],
                                                             target_col = varj,
                                                             cond_cols  = Sj,
                                                             k = k,
                                                             seed = as.numeric(Sys.time()))
            
            # calculating the C-index after conditional permutation
            perm_cindex <- compute_wolbers_cindex(fg, df_perm)
            
            drops_c_vector[b] <- baseline_cindex - perm_cindex
        }
        
        # C-index 검정 (two-tailed test로 변경: 양수/음수 모두 의미 있게)
        drop_cindex_vector[varj] <- mean(drops_c_vector)
        
        if (sd(drops_c_vector) < 1e-10 || length(unique(drops_c_vector)) == 1) {
            cindex_pvalue_vector[varj] <- ifelse(abs(mean(drops_c_vector)) > 0, 0.001, 0.999)
            cindex_significant_vector[varj] <- abs(mean(drops_c_vector)) > 0
        } else {
            cindex_test <- t.test(drops_c_vector, mu = 0, alternative = "two.sided")
            cindex_pvalue_vector[varj] <- cindex_test$p.value
            cindex_significant_vector[varj] <- cindex_test$p.value < alpha
        }
        
        # 결과 출력 (유의성 표시)
        cindex_sig_mark <- ifelse(cindex_significant_vector[varj], "***", "")
        
        cat(varj, ": mean(ΔC-index) =", round(drop_cindex_vector[varj], 6), cindex_sig_mark,
            " (p_value=", round(cindex_pvalue_vector[varj], 4), ")\n")
    }
    
    return(list(drop_cindex = drop_cindex_vector,
                cindex_pvalue = cindex_pvalue_vector,
                cindex_significant = cindex_significant_vector,
                baseline_cindex = baseline_cindex,
                alpha = alpha))
}


# 메인 실행 부분
# Loading Model and data
loaded_fg <- readRDS("fitted_fg.rds")

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
               'alb', 'log_tbil', 'cr', 'log_pt_inr', 'na', # lab test
               'age_lt', 'bmi' # demographic
)

# variables to be transformed log scale
vars_to_log <- c('afp', 'pivka', 'tbil', 'pt_inr', 
                 'lr_5_size_ct', 'lr_tr_v_size_whole_ct',
                 'lr_tr_v_size_enhancing_ct', 'lr_3_size_ct',
                 'lr_4_size_ct', 'lr_m_size_ct', 'lr_tr_nv_size_ct',
                 'lr_tr_e_size_ct'
)

# Test_data
cat("Loading test data\n")
test_data_fname <- "../../data/preprocessed_competing_risk_validation_dataset250624.csv"
test_data <- read.csv(test_data_fname)

# Several variables are transformed to log scale
vars_to_log <- intersect(vars_to_log, names(test_data))
for (var_nm in vars_to_log) {
    new_column_nm <- paste0("log_", var_nm)
    test_data[[new_column_nm]] <- log(test_data[[var_nm]]+1)
}
test_data[factor_vars] <- lapply(test_data[factor_vars], factor)
cat("Test data done\n")

# Wolbers C-index 기반 CPI 실행
res <- cpi_fg_wolbers(fg = loaded_fg,
                     data_val = test_data,
                     timeVar = "ftime",
                     statusVar = "fstatus",
                     features = attr(loaded_fg$terms, "term.labels"),
                     max_time = NULL,  # 자동으로 최대 이벤트 시간 사용
                     B = 50,
                     k = 30,
                     alpha = 0.05,
                     seed = 20250521)

# Variable name replacement
var_name_mapping <- c(
  "log_lr_5_size_ct" = "Maximum diameter of LR-5",
  "log_lr_tr_v_size_enhancing_ct" = "Maximum diameter of LR-TR-V",
  "log_lr_m_size_ct" = "Maximum diameter of LR-M",
  "log_lr_tr_nv_size_ct" = "Maximum diameter of LR-TR-NV",
  "log_afp" = "Serum AFP",
  "log_pivka" = "Serum PIVKA-II",
  "log_pt_inr" = "Prothrombin time, INR",
  "sex" = "Sex",
  "ckd_hd" = "Chronic kidney disease"
)

# Replace names for all vectors
current_names <- names(res$drop_cindex)
new_names <- ifelse(current_names %in% names(var_name_mapping),
                   var_name_mapping[current_names],
                   current_names)

names(res$drop_cindex) <- new_names
names(res$cindex_pvalue) <- new_names
names(res$cindex_significant) <- new_names

print(str(res))
saveRDS(res, "cpi_fg_wolbers_results.rds")

# 2. 결과 해석 (필터링 없음)
interpretation <- interpret_importance(res)

# 3. 추천 변수 선택
recommended_vars <- recommend_features(interpretation, 
                                      importance_threshold = 0.0001,
                                      require_significance = TRUE)

# 4. 시각화 (모든 변수 포함)
plot_importance_barplot(interpretation,
                       filename_prefix = "conditional_perm_vimp_wolbers",
                       width = 12,
                       height = 10)

# 5. 결과 요약 출력 (모든 변수 표시)
cat("\n=== 모든 변수의 Wolbers C-index 중요도 (높은 순) ===\n")
cindex_df <- interpretation$summary
print(cindex_df)

cat("\n=== 통계적으로 유의한 양의 중요도 변수 ===\n")
significant_positive <- cindex_df[cindex_df$Classification == "Significantly Positive", ]
if (nrow(significant_positive) > 0) {
    print(significant_positive)
} else {
    cat("No significantly positive variables found.\n")
}

cat("\n=== 통계적으로 유의한 음의 중요도 변수 ===\n")
significant_negative <- cindex_df[cindex_df$Classification == "Significantly Negative", ]
if (nrow(significant_negative) > 0) {
    print(significant_negative)
} else {
    cat("No significantly negative variables found.\n")
}
