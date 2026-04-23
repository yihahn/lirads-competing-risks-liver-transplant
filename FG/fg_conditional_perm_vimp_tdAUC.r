library(survival)
library(riskRegression)
library(prodlim)
library(FNN)
library(dplyr)
library(pheatmap)


# 결과 해석 함수들
interpret_importance <- function(cpi_results, min_effect_size = 0.0001) {

	# 1. 통계적으로 유의한 positive importance만 추출
	get_significant_positive <- function(drops, significant, threshold = min_effect_size) {
		result <- drops
		result[!significant | drops < threshold] <- NA
		return(result)
	}

	significant_auc <- get_significant_positive(cpi_results$drop_auc,
						    cpi_results$auc_significant,
						    min_effect_size)
	significant_brier <- get_significant_positive(cpi_results$drop_brier,
						      cpi_results$brier_significant,
						      min_effect_size)

	# 2. 변수별 종합 중요도 계산 (시간대별 평균)
	auc_overall <- rowMeans(significant_auc, na.rm = TRUE)
	print(auc_overall)
	brier_overall <- rowMeans(significant_brier, na.rm = TRUE)

	# 3. 변수 분류
	classify_variables <- function(auc_imp, brier_imp) {
		classification <- character(length(auc_imp))
		names(classification) <- names(auc_imp)

		for (i in seq_along(auc_imp)) {
			auc_val <- auc_imp[i]
			brier_val <- brier_imp[i]

			if (is.na(auc_val) && is.na(brier_val)) {
				classification[i] <- "Not Important"
			} else if (!is.na(auc_val) && !is.na(brier_val)) {
				classification[i] <- "Consistently Important"
			} else if (!is.na(auc_val)) {
				classification[i] <- "AUC Important Only"
			} else {
				classification[i] <- "Brier Important Only"
			}
		}
		return(classification)
	}

	var_classification <- classify_variables(auc_overall, brier_overall)

	# 4. 결과 정리
	result_df <- data.frame(
		Variable = names(auc_overall),
		AUC_Importance = round(auc_overall, 6),
		Brier_Importance = round(brier_overall, 6),
		Classification = var_classification,
		stringsAsFactors = FALSE
	)

	# NaN을 NA로 변경
	result_df$AUC_Importance[is.nan(result_df$AUC_Importance)] <- NA
	result_df$Brier_Importance[is.nan(result_df$Brier_Importance)] <- NA

	# 중요도 순으로 정렬 (AUC + Brier 평균, NA는 뒤로)
	combined_importance <- rowMeans(cbind(auc_overall, brier_overall), na.rm = TRUE)
	combined_importance[is.nan(combined_importance)] <- -Inf
	result_df <- result_df[order(combined_importance, decreasing = TRUE), ]

	return(list(
		summary = result_df,
		significant_auc = significant_auc,
		significant_brier = significant_brier,
		raw_results = cpi_results
	))
}


plot_importance_heatmap <- function(interpretation_results,
                                  filename_prefix = "feature_importance",
                                  width = 10,
                                  height = 8) {

	# AUC heatmap
	auc_data <- interpretation_results$significant_auc
	auc_data[is.na(auc_data)] <- 0  # NA를 0으로 표시

	# AUC heatmap 저장
	auc_filename <- paste0(filename_prefix, "_auc_fg.pdf")
	pdf(auc_filename, width = width, height = height)
	pheatmap(auc_data,
	         cluster_rows = TRUE, cluster_cols = FALSE,
	         main = "Significant AUC Importance by Time",
	         color = colorRampPalette(c("white", "red"))(50),
	         display_numbers = TRUE,
	         number_format = "%.4f")
	dev.off()
	cat("AUC heatmap saved as:", auc_filename, "\n")

	# Brier heatmap
	brier_data <- interpretation_results$significant_brier
	brier_data[is.na(brier_data)] <- 0

	# Brier heatmap 저장
	brier_filename <- paste0(filename_prefix, "_bs_fg.pdf")
	pdf(brier_filename, width = width, height = height)
	pheatmap(brier_data,
	         cluster_rows = TRUE, cluster_cols = FALSE,
	         main = "Significant Brier Importance by Time",
	         color = colorRampPalette(c("white", "blue"))(50),
	         display_numbers = TRUE,
	         number_format = "%.4f")
	dev.off()
	cat("Brier heatmap saved as:", brier_filename, "\n")
}


# 추천 변수 선택 함수
recommend_features <- function(interpretation_results,
                              importance_threshold = 0.005,
                              require_consistency = TRUE) {

	summary_df <- interpretation_results$summary

	if (require_consistency) {
		# 두 메트릭 모두에서 중요한 변수만 선택
		recommended <- summary_df[summary_df$Classification == "Consistently Important" &
		                         !is.na(summary_df$AUC_Importance) &
		                         !is.na(summary_df$Brier_Importance) &
		                         summary_df$AUC_Importance >= importance_threshold &
		                         summary_df$Brier_Importance >= importance_threshold, ]
	} else {
		# 하나의 메트릭에서라도 중요한 변수 선택
		recommended <- summary_df[(!is.na(summary_df$AUC_Importance) &
		                          summary_df$AUC_Importance >= importance_threshold) |
		                         (!is.na(summary_df$Brier_Importance) &
		                          summary_df$Brier_Importance >= importance_threshold), ]
	}

        # AUC_Importance 기준으로 내림차순 정렬
        recommended <- recommended[order(recommended$AUC_Importance, decreasing = TRUE), ]

	cat("=== Recommended Features ===\n")
	print(recommended)
	cat("\n=== Feature Selection Summary ===\n")
	cat("Total variables tested:", nrow(summary_df), "\n")
	cat("Consistently important:", sum(summary_df$Classification == "Consistently Important"), "\n")
	cat("AUC important only:", sum(summary_df$Classification == "AUC Important Only"), "\n")
	cat("Brier important only:", sum(summary_df$Classification == "Brier Important Only"), "\n")
	cat("Not important:", sum(summary_df$Classification == "Not Important"), "\n")
	cat("Recommended features:", nrow(recommended), "\n")

	return(recommended$Variable)
}

# Score()로 C-index와 Brier 계산
#    - data: validation data
#    - pred: cif of prediction from rsf object
#    - eval_times: evaluation time
compute_metrics <- function(pred, data, eval_times) {

	aucs <- c()
	briers <- c()
	for (i in seq_along(eval_times)) {
		eval_time <- eval_times[i]
		score <- Score(list("model"=pred[, i]),
			       formula = Hist(ftime, fstatus) ~ 1,
			       data = data,
			       times = eval_time,
			       metrics = c("Brier", "AUC"), # AUC는 1-스케일 C-index와 연결됨
			       cause = 1,
			       null.model = FALSE)
		auc <- score$AUC$score %>%
			filter(model=="model") %>%
			select(AUC) %>%
			pull()
		bs <- score$Brier$score %>%
			filter(model=="model") %>%
			select(Brier) %>%
			pull()
		aucs <- c(aucs, auc)
		briers <- c(briers, bs)
	}
	return(list(auc = aucs,
		    brier = briers))
}


# kNN 조건부 permutation으로 변수 하나를 교란한 data.frame을 생성
# - S = 조건집합(기본은 X_{-j})
# - 연속형은 표준화 후 거리 계산
# - 범주형은 간단히 더미화(one-hot) 또는 같은 카테고리 우선 필터 등을 추가 가능
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



# fg      : 학습된 FGR 객체 
# data_val : 검증용 data.frame (time/status/feature 포함)
# timeVar, statusVar : 생존시간/검열변수 컬럼명
# features : CPI를 계산할 설명변수 벡터 (NULL이면 time/status 제외 전부)
# times    : 평가 시점(단일 값 권장; 예: 중간추적기간 또는 임상적 관심 시점)
# B        : 조건부 permutation 반복 수
# k        : kNN 이웃 수
# return   : baseline metric & 변수별 {mean, sd, values} for ΔC, ΔBrier
cpi_fg <- function(fg,
		   data_val,
		   timeVar = "ftime",
		   statusVar = "fstatus",
		   features = NULL,
		   eval_times,
		   B = 50,
		   k = 30,
		   alpha = 0.05,
		   seed = 20250521) {
	set.seed(seed)

	# 피처 목록 확정
	if (is.null(features)) {
		features <- setdiff(names(data_val), c(timeVar, statusVar))
	}

	# prediction after conditional permutation
	sink("/dev/null")
	pred <- predict(fg, newdata=data_val, times=eval_times, cause=1)
	sink()
	# calculating the baseline score before permutation
	base <- compute_metrics(pred, data_val, eval_times)
	M0_a <- base$auc
	M0_b <- base$brier

	# 결과를 저장할 matrix 생성
	drop_auc_matrix <- matrix(NA, nrow = length(features), ncol = length(eval_times))
	drop_brier_matrix <- matrix(NA, nrow = length(features), ncol = length(eval_times))
	auc_pvalue_matrix <- matrix(NA, nrow = length(features), ncol = length(eval_times))
	brier_pvalue_matrix <- matrix(NA, nrow = length(features), ncol = length(eval_times))
	auc_significant_matrix <- matrix(FALSE, nrow = length(features), ncol = length(eval_times))
	brier_significant_matrix <- matrix(FALSE, nrow = length(features), ncol = length(eval_times))

	rownames(drop_auc_matrix) <- features
	rownames(drop_brier_matrix) <- features
	colnames(drop_auc_matrix) <- paste0("t_", eval_times)
	colnames(drop_brier_matrix) <- paste0("t_", eval_times)
	rownames(auc_pvalue_matrix) <- features
	rownames(brier_pvalue_matrix) <- features
	colnames(auc_pvalue_matrix) <- paste0("t_", eval_times)
	colnames(brier_pvalue_matrix) <- paste0("t_", eval_times)
	rownames(auc_significant_matrix) <- features
	rownames(brier_significant_matrix) <- features
	colnames(auc_significant_matrix) <- paste0("t_", eval_times)
	colnames(brier_significant_matrix) <- paste0("t_", eval_times)

	for (j in seq_along(features)) {
		varj <- features[j]
		Sj   <- setdiff(features, varj)

		# ΔC = C0 - Cperm (클수록 중요)
		drops_a_matrix <- matrix(NA, nrow=B, ncol=length(eval_times))
		# ΔB = Bperm - B0 (클수록 중요; Brier는 작을수록 좋음)
		drops_b_matrix <- matrix(NA, nrow=B, ncol=length(eval_times))
		colnames(drops_a_matrix) <- paste0("t_", eval_times)
		colnames(drops_b_matrix) <- paste0("t_", eval_times)

		for (b in 1:B) {
			df_perm <- data_val
			# 조건부(kNN)로 varj만 치환
			df_perm[, c(varj, Sj)] <- conditional_permute_knn(df = data_val[, c(varj, Sj), drop = FALSE],
									  target_col = varj,
									  cond_cols  = Sj,
									  k = k,
									  seed = seed + j * 1000 + b)
			# prediction after conditional permutation
			sink("/dev/null")
			pred <- predict(fg, newdata=df_perm, times=eval_times, cause=1)     
			sink()
			# calculating the scores after conditional permutation
			mb <- compute_metrics(pred, df_perm, eval_times)
			
			drops_a_matrix[b, ] <- M0_a - mb$auc
			drops_b_matrix[b, ] <- mb$brier - M0_b

		}

		for (t_idx in seq_along(eval_times)) {
			# AUC 검정 (one-tailed test: H0: mean <= 0 vs H1: mean > 0)
			auc_drops <- drops_a_matrix[, t_idx]
			drop_auc_matrix[varj, t_idx] <- mean(auc_drops)
			if (sd(auc_drops) < 1e-10 || length(unique(auc_drops))==1) {
				auc_pvalue_matrix[varj, t_idx] <- ifelse(mean(auc_drops) > 0, 0.001, 0.999)
				auc_significant_matrix[varj, t_idx] <- mean(auc_drops) > 0
			} else {
				auc_test <- t.test(auc_drops, mu = 0, alternative = "greater")
				auc_pvalue_matrix[varj, t_idx] <- auc_test$p.value
				auc_significant_matrix[varj, t_idx] <- auc_test$p.value < alpha
			}

			# Brier 검정
			brier_drops <- drops_b_matrix[, t_idx]
			drop_brier_matrix[varj, t_idx] <- mean(brier_drops)
			if (sd(brier_drops) < 1e-10 || length(unique(brier_drops))==1) {
				brier_pvalue_matrix[varj, t_idx] <- ifelse(mean(brier_drops) > 0, 0.001, 0.999)
				brier_significant_matrix[varj, t_idx] <- mean(brier_drops) > 0
			} else {
				brier_test <- t.test(brier_drops, mu = 0, alternative = "greater")
				brier_pvalue_matrix[varj, t_idx] <- brier_test$p.value
				brier_significant_matrix[varj, t_idx] <- brier_test$p.value < alpha
			}

			# 결과 출력 (유의성 표시)
			auc_sig_mark <- ifelse(auc_significant_matrix[varj, t_idx], "***", "")
			brier_sig_mark <- ifelse(brier_significant_matrix[varj, t_idx], "***", "")

			cat(varj, "t_", eval_times[t_idx], ": mean(ΔAUC) =",
			    round(drop_auc_matrix[varj, t_idx], 6), auc_sig_mark,
			    ", mean(ΔBS) =", round(drop_brier_matrix[varj, t_idx], 6), brier_sig_mark,
			    " (p_auc=", round(auc_pvalue_matrix[varj, t_idx], 4),
			    ", p_brier=", round(brier_pvalue_matrix[varj, t_idx], 4), ")\n")
		}


	}
	return(list(drop_auc = drop_auc_matrix,
		    drop_brier = drop_brier_matrix,
		    auc_pvalue = auc_pvalue_matrix,
		    brier_pvalue = brier_pvalue_matrix,
		    auc_significant = auc_significant_matrix,
		    brier_significant = brier_significant_matrix,
		    alpha = alpha
	))
}

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

t_stars <- c(365.24*1, 365.24*3, 365.24*5, 365.24*10)
res <- cpi_fg(fg = loaded_fg,
	      data_val = test_data,
	      timeVar = "ftime",
	      statusVar = "fstatus",
	      features = attr(loaded_fg$terms, "term.labels"),
	      eval_times = t_stars,
	      B = 50,
	      k = 30,
	      alpha = 0.05,
	      seed=20250521)

# Time label replacement
time_labels <- c("1 yr", "3 yrs", "5 yrs", "10 yrs")
# Replace column names for drop_auc
colnames(res$drop_auc) <- time_labels
# Replace column names for auc_pvalue
colnames(res$auc_pvalue) <- time_labels
# Replace column names for auc_significant
colnames(res$auc_significant) <- time_labels

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

# Replace row names for all matrices
current_row_names <- rownames(res$drop_auc)
new_row_names <- ifelse(current_row_names %in% names(var_name_mapping),
                       var_name_mapping[current_row_names],
                       current_row_names)

rownames(res$drop_auc) <- new_row_names
rownames(res$auc_pvalue) <- new_row_names
rownames(res$auc_significant) <- new_row_names

print(str(res))
saveRDS(res, "cpi_fg_results.rds")

# 2. 결과 해석
interpretation <- interpret_importance(res, min_effect_size = 0.0001)

# 3. 추천 변수 선택
recommended_vars <- recommend_features(interpretation, 
				       importance_threshold = 0.0001,
				       require_consistency = FALSE)

# 4. 시각화
plot_importance_heatmap(interpretation,
                       filename_prefix = "conditioned_perm_feat_imp",    # 파일명 접두사
                       width = 12,                         # PDF 너비
                       height = 10)                        # PDF 높이

stop()
# AUC 중요도로 정렬 (내림차순: 클수록 중요)
cat("\n=== AUC 기준 변수 중요도 (높은 순) ===\n")
auc_importance <- sapply(res$importances, function(x) x$auc$mean)
auc_sorted <- sort(auc_importance, decreasing = FALSE)
auc_df <- data.frame(Variable = names(auc_sorted),
		     AUC_Drop_Mean = round(auc_sorted, 4),
		     AUC_Drop_SD = round(sapply(names(auc_sorted), function(v) res$importances[[v]]$auc$sd), 4),
		     stringsAsFactors = FALSE)
print(auc_df)

# Brier Score 중요도로 정렬 (내림차순: 클수록 중요)
cat("\n=== Brier Score 기준 변수 중요도 (높은 순) ===\n")
brier_importance <- sapply(res$importances, function(x) x$brier$mean)
brier_sorted <- sort(brier_importance, decreasing = FALSE)
brier_df <- data.frame(Variable = names(brier_sorted),
		       Brier_Increase_Mean = round(brier_sorted, 4),
		       Brier_Increase_SD = round(sapply(names(brier_sorted), function(v) res$importances[[v]]$brier$sd), 4),  
		       stringsAsFactors = FALSE)
print(brier_df)

# 상위 10개 변수만 출력 (선택사항)
cat("\n=== 상위 10개 중요 변수 ===\n")
cat("AUC 기준:\n")
print(head(auc_df, 10))
cat("\nBrier Score 기준:\n")
print(head(brier_df, 10))

# 두 지표 모두 고려한 종합 순위 (선택사항)
cat("\n=== 종합 순위 (AUC + Brier 평균 순위) ===\n")
auc_ranks <- rank(auc_importance)  # 큰 값이 좋으므로 음수로 변환
brier_ranks <- rank(brier_importance)  # 큰 값이 좋으므로 음수로 변환
combined_ranks <- (auc_ranks + brier_ranks) / 2
combined_sorted <- sort(combined_ranks)
combined_df <- data.frame(Variable = names(combined_sorted),
			  Combined_Rank = round(combined_sorted, 2),
			  AUC_Drop_Mean = round(auc_importance[names(combined_sorted)], 4),
			  Brier_Increase_Mean = round(brier_importance[names(combined_sorted)], 4),
			  stringsAsFactors = FALSE)
print(combined_df)
