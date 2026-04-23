library(caret)
library(dplyr)
library(randomForestSRC)
library(prodlim)
library(survival)
library(riskRegression)
library(timeROC)
library(pracma)
library(jsonlite)

# function for valid factor variable selection 
check_n_filter_factor_vars <- function(data, var_names) {
	valid_factor_vars <- c()
	for (fv in var_names) {
		# 1. 변수가 데이터프레임에 존재하는지 확인
		if (fv %in% colnames(data)) {
			# 2. 변수가 factor 타입인지 확인
			if (is.factor(data[[fv]])) {
				# 3. 고유한 수준의 개수 확인
				num_levels <- length(levels(data[[fv]]))
				if (num_levels < 2) {
					cat(paste0("Notice: '", fv, 
						   "' is removed due to a single value variable.\n"))
				} else {
					valid_factor_vars <- c(valid_factor_vars, fv)
				}
			} else {
				cat(paste0("Notice: '", fv, 
					   "' is removed due to not being a factor.\n"))
			}
		} else {
			cat(paste0("Notice: '", fv, "' does not exist in dataframe.\n"))
		}
	}
	return(valid_factor_vars)
}



# 단일값 변수 확인 및 제거 함수
check_single_value_vars <- function(data, var_names) {
	valid_vars <- c()
	for (var in var_names) {
		if (var %in% colnames(data)) {
			unique_vals <- length(unique(data[[var]][!is.na(data[[var]])]))
			if (unique_vals <= 1) {
				cat(paste0("Notice: '", var, "' is removed due to single value variable.\n"))
			} else {
				valid_vars <- c(valid_vars, var)
			}
		}
	}
	return(valid_vars)
}

# RSF interpolation for evaluation times
interp_predict_risk_rsf <- function(pred, times) {
	time_vec <- pred$time.interest
	interp_cif_matrix <- matrix(NA, nrow=nrow(pred$cif), ncol=length(times))
	for (k in 1:nrow(pred$cif)) {
		individual_cif <- pred$cif[k, , 1]
		for (i in seq_along(times)) {
			# Use approx for interpolation for the current individual
			interpolated_cif_at_time <- approx(x=time_vec, 
							   y=individual_cif, 
							   xout=times[i],
							   rule=2)$y
			interp_cif_matrix[k, i] <- interpolated_cif_at_time
		}
	}

	colnames(interp_cif_matrix) <- paste("t", times, sep="_")
	interp_cif <- as.data.frame(interp_cif_matrix)
	return(interp_cif)
}


set.seed(20250521)

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

#optim_vars_for_ext <- c("log_lr_5_size_ct", "log_lr_tr_v_size_whole_ct", "log_lr_tr_nv_size_ct", 
#			"lr_5_ct", "lr_tr_v_ct", "lr_m_ct", "lr_tr_nv_ct", 
#			"log_afp", "log_pivka", 
#			"cr", "log_pt_inr", 
#			"age_lt", "ckd_hd", "hbv", "ldlt")#, "pre_tx")

optim_vars_for_ext <- c("log_lr_5_size_ct", "log_lr_tr_v_size_whole_ct", "log_lr_tr_nv_size_ct", 
                        "lr_5_ct", "lr_tr_v_ct", "lr_m_ct", "lr_tr_nv_ct", "lr_tr_e_ct",
						"log_afp", "log_pivka",
						"cr", "log_pt_inr", 
                        "age_lt", "ckd_hd", "hbv", "ldlt", "rfa")
optim_ntree <- 500
optim_nodesize <- 10

## Performance test with external validation dataset
eval_times <- c(365.24, 365.24*3, 365.24*5, 365.24*10)
tpr_by_times <- vector("list", length(eval_times))
fpr_by_times <- vector("list", length(eval_times))
names(tpr_by_times) <- paste0("t=", eval_times)
names(fpr_by_times) <- paste0("t=", eval_times)

# Data entry
# Train_data
train_data_fname = "../../data/internal_competing_risk_dataset.csv"
train_data <- read.csv(train_data_fname)

# Several variables are transformed to log scale
for (var_nm in vars_to_log) {
		new_column_nm <- paste0("log_", var_nm)
		train_data[[new_column_nm]] <- log(train_data[[var_nm]]+1)
}
train_data[factor_vars] <- lapply(train_data[factor_vars], factor)
val_factor_vars_train <- check_n_filter_factor_vars(train_data, factor_vars)
val_cont_vars_train <- intersect(cont_vars, names(train_data))
cat("Train data done\n")


# Test_data
test_data_fname = "../../data/preprocessed_competing_risk_validation_dataset250624.csv"
test_data <- read.csv(test_data_fname)
# Several variables are transformed to log scale
vars_to_log <- intersect(vars_to_log, names(test_data))
for (var_nm in vars_to_log) {
		new_column_nm <- paste0("log_", var_nm)
		test_data[[new_column_nm]] <- log(test_data[[var_nm]]+1)
}
test_data[factor_vars] <- lapply(test_data[factor_vars], factor)
val_factor_vars_test <- check_n_filter_factor_vars(test_data, factor_vars)
val_cont_vars_test <- intersect(cont_vars, names(test_data))
cat("Test data done\n")

# intersect of valid factor variables both in train_data and test_data
cmmn_val_factor_vars <- intersect(val_factor_vars_train, val_factor_vars_test)
cmmn_val_cont_vars <- intersect(val_cont_vars_train, val_cont_vars_test)

# combine variables (only valid factor variable included)
all_vars = c(cmmn_val_cont_vars, cmmn_val_factor_vars)
final_vars = intersect(all_vars, unlist(optim_vars_for_ext))

# remove single value variable
final_vars <- check_single_value_vars(train_data, final_vars)
final_vars <- check_single_value_vars(test_data, final_vars)

if (length(final_vars) == 0) {
	stop("No valid variables available")
}
 
# Start to do regression
formula_str <- paste("Surv(ftime, fstatus) ~", 
		     paste(final_vars, collapse=" + "))
print(formula_str)
rsf_formula <- as.formula(formula_str)

# cv unbiased predicted risk for calibration
time_labels <- paste0("t_", eval_times)
create_cv_folds <- function(data, k=3) {
	set.seed(20250521)
	folds <- createFolds(data$fstatus, k=k, list=TRUE)
	return (folds)
}

cv_calibrate_data <- list(predictions=data.frame(),
			  outcomes=data.frame(),
			  eval_times=eval_times,
			  time_labels=time_labels,
			  fold_info=list()
)

folds <- create_cv_folds(train_data)

c_indices <- numeric(length(folds))
auc_per_fold <- list()
brier_per_fold <- list()
ibs_per_fold <- list()
for (i in seq_along(folds)) {
	test_indices <- folds[[i]]
	cv_train_data <- train_data[-test_indices, ]
	cv_test_data <- train_data[test_indices, ]
        # Adjust mtry if necessary
	cv_fit <- rfsrc(rsf_formula,
			data=cv_train_data,
			tree=optim_ntree,
			nodesize=optim_nodesize,
			splitrule = "logrankCR",
			na.action="na.omit")
	cindex_fold <- suppressMessages(pec::cindex(list(rsf=cv_fit),
						    formula=rsf_formula,
						    data=cv_test_data,
						    cause=1,
						    verbose=FALSE))

	c_indices[i] <- cindex_fold$AppCindex$rsf

	score <- Score(list("RSF"=cv_fit),
		       formula=Hist(ftime, fstatus) ~ 1,
		       data=cv_test_data,
		       times=eval_times,
		       metrics=c("AUC", "Brier"),
		       cause=1, # interest event label
		       null.model=FALSE)
	auc_per_fold[[i]] <- score$AUC$score %>%
		filter(model=="RSF") %>%
		select(AUC) %>%
		pull()
	brier_per_fold[[i]] <- score$Brier$score %>%
		filter(model=="RSF") %>%
		select(Brier) %>%
		pull()
		
	# Integrated Brier score by trapz
	total_event_time_ibs <- subset(test_data, fstatus %in% c(1, 2))$ftime
	total_event_time_ibs <- sort(unique(total_event_time_ibs))
	brier_scores <- Score(list("RSF"=cv_fit),
			      formula=Hist(ftime, fstatus) ~ 1,
			      data=cv_test_data,
			      times=total_event_time_ibs,
			      metrics="Brier",
			      cause=1
	)

        brier_scores_dt <- brier_scores$Brier$score
        rsf_brier_scores <- brier_scores_dt[model=="RSF", .(times, Brier)]
        ibs_rsf <- trapz(rsf_brier_scores$times, rsf_brier_scores$Brier)
        ibs_per_fold[[i]] <- ibs_rsf / max(total_event_time_ibs)
		
        cv_pred <- predict(cv_fit,
			   newdata=cv_test_data)
	pred_df <- interp_predict_risk_rsf(cv_pred, eval_times)
	colnames(pred_df) <- time_labels 
	pred_df$row_id <- test_indices
	pred_df$fold <- i

	# extract outcomes for each time points
	outcome_df <- data.frame(row_id=test_indices,
				 fold=i,
				 ftime=cv_test_data$ftime,
				 fstatus=cv_test_data$fstatus
	)
	for (j in seq_along(eval_times)) {
		time_point <- eval_times[j]
		time_label <- time_labels[j]

		# Event occurred before time_point (competing risk event)
		outcome_df[[paste0("event_", time_label)]] <- 
			ifelse(cv_test_data$ftime <= time_point & cv_test_data$fstatus == 1, 1, 0)

		# Censored or survived beyond time_point
		outcome_df[[paste0("censored_", time_label)]] <-
			ifelse(cv_test_data$ftime <= time_point & cv_test_data$fstatus != 1, 1, 0)

		# Available for analysis at this time point (not censored before)
		outcome_df[[paste0("available_", time_label)]] <-
			ifelse(cv_test_data$ftime > time_point |
			       (cv_test_data$ftime <= time_point & cv_test_data$fstatus == 1), 1, 0)
	}

	# Append to main storage
	cv_calibrate_data$predictions <- rbind(cv_calibrate_data$predictions, pred_df)
	cv_calibrate_data$outcomes <- rbind(cv_calibrate_data$outcomes, outcome_df)
	cv_calibrate_data$fold_info[[i]] <- list(fold_number=i,
						   test_indices=test_indices,
						   n_test=length(test_indices)
	)
}

# summary for each fold and each time point
auc_matrix <- do.call(rbind, auc_per_fold)
brier_matrix <- do.call(rbind, brier_per_fold)

year_labels <- round(eval_times / 365.24) 
colnames(auc_matrix) <- paste0("AUC_", year_labels, "yr")
colnames(brier_matrix) <- paste0("Brier_", year_labels, "yr")

cat("\n--- Performance Summary by Time Point ---\n")
#AUC summary
mean_auc_at_times <- colMeans(auc_matrix, na.rm=TRUE)
sd_auc_at_times <- apply(auc_matrix, 2, sd, na.rm=TRUE)
cat("Time-dependent AUROC results:\n")
for (i in seq_along(year_labels)) {
	cat(sprintf("  At %d year(s): Mean = %.3f (SD=%.3f)\n",
		    year_labels[i], mean_auc_at_times[i], sd_auc_at_times[i]))
}

#Brier summary
mean_brier_at_times <- colMeans(brier_matrix, na.rm=TRUE)
sd_brier_at_times <- apply(brier_matrix, 2, sd, na.rm=TRUE)
cat("Time-dependent Brier score results:\n")
for (i in seq_along(year_labels)) {
	cat(sprintf("  At %d year(s): Mean = %.3f (SD=%.3f)\n",
		    year_labels[i], mean_brier_at_times[i], sd_brier_at_times[i]))
}

ibs <- unlist(ibs_per_fold)
cat(sprintf("5-fold mean Integrated Brier score: %.3f (sd=%.3f)\n", mean(ibs), sd(ibs)))
cat(sprintf("5-fold mean Wolbers C-index: %.3f (sd=%.3f)\n", mean(c_indices), sd(c_indices)))


# Sort by row_id for easier access
cv_calibrate_data$predictions <- cv_calibrate_data$predictions[order(cv_calibrate_data$predictions$row_id), ]
cv_calibrate_data$outcomes <- cv_calibrate_data$outcomes[order(cv_calibrate_data$outcomes$row_id), ]

# Save to RDS
#print("saving rsf_cv_data_for_calibration.rds")
#saveRDS(cv_calibrate_data, "rsf_cv_data_for_calibration.rds")
#print(str(cv_calibrate_data))

complete_cases_test <- complete.cases(test_data[c("ftime", "fstatus", final_vars)])
test_data_clean <- test_data[complete_cases_test, ]

test_covariates <- test_data_clean %>%
	select(all_of(final_vars))

# RSF fitting
tryCatch({fit <- rfsrc(rsf_formula,
	data=train_data,
	ntree=optim_ntree,
	nodesize=optim_nodesize,
	splitrule = "logrankCR",
	na.action="na.omit"
		     )
	#saveRDS(fit, "fitted_rsf.rds")

        # Wolbers C-index 계산
	tryCatch({wolbers_cindex <- pec::cindex(list(rsf_model=fit),
		formula=rsf_formula,
		data=test_data_clean, 
		cause=1)
		cindex <- wolbers_cindex$AppCindex$rsf_model
		print(paste("Wolbers C-index =", cindex))
	}, error = function(e) {
		cat("Error in Wolbers C-index calculation:", e$message, "\n")
		cindex <- NA
	})

        # predict test data
        pred <- predict(fit,
			newdata=test_covariates)

	interp_cif <- interp_predict_risk_rsf(pred, eval_times)
	#saveRDS(interp_cif, file="ext_rsf_pred_by_evaltimes.rds")

	# AUROC & Brier from riskRegression
	tryCatch({score <- Score(list("RSF"=fit),
		#formula=rsf_formula,
		formula=Hist(ftime, fstatus) ~ 1,
		data=test_data_clean,
		times=eval_times,
		metrics=c("AUC", "Brier"),
		cause=1, # interest event label
		null.model=FALSE)

		auroc <- score$AUC$score %>%
			filter(model=="RSF") %>%
			select(AUC) %>%
			pull()
		bs <- score$Brier$score %>%
			filter(model=="RSF") %>%
			select(Brier) %>%
			pull()
		}, error = function(e) {
			cat("Error in Score calculation:", e$message, "\n")
			auroc <- rep(NA, length(eval_times))
			bs <- rep(NA, length(eval_times))
	})

	# Integrated Brier score by trapz
	total_event_time_ibs <- subset(test_data, fstatus %in% c(1, 2))$ftime
	brier_scores <- Score(list("RSF"=fit),
			      formula=Hist(ftime, fstatus) ~ 1,
			      data=test_data,
			      times=total_event_time_ibs,
			      metrics="Brier",
			      cause=1
	)

        brier_scores_dt <- brier_scores$Brier$score
        rsf_brier_scores <- brier_scores_dt[model=="RSF", .(times, Brier)]
        ibs_rsf <- trapz(rsf_brier_scores$times, rsf_brier_scores$Brier)
        ibs <- ibs_rsf/max(total_event_time_ibs)

        # timeROC 계산
	
	for (i in seq_along(eval_times)) {
		tryCatch({
			t_pnt <- eval_times[i]
			colname_interp <- paste("t", t_pnt, sep="_")
			current_pred_at_time <- interp_cif[[colname_interp]]
			if (length(current_pred_at_time) == 0) {
				cat(sprintf("At %.2f, no prediction value\n", t_pnt))
			}
			# 길이 확인 및 출력
			cat(sprintf("Time point %s: T length=%d, delta length=%d, marker length=%d\n",
				    t_pnt, length(test_data_clean$ftime),
				    length(test_data_clean$fstatus),
				    length(current_pred_at_time)))
			roc_obj <- timeROC(T=test_data_clean$ftime,
					   delta=test_data_clean$fstatus,
					   marker=current_pred_at_time,
					   cause=1,
					   times=t_pnt,
					   ROC=TRUE
			)
			t_pnt_colname <- paste("t=", t_pnt, sep="")
			current_tpr <- roc_obj$TP[, t_pnt_colname]
			current_fpr <- roc_obj$FP_1[, t_pnt_colname]
			tpr_by_times[[i]] <- c(tpr_by_times[[i]], list(current_tpr))
			fpr_by_times[[i]] <- c(fpr_by_times[[i]], list(current_fpr))
		}, error = function(e) {
			cat("Error in timeROC calculation for time", eval_times[i], ":", e$message, "\n")
               })
	}
}, error = function(e) {
	cat("Error in RSF model fitting :", e$message, "\n")
	## 기본값 설정
	auroc <- rep(NA, length(eval_times))
	bs <- rep(NA, length(eval_times))
	cindex <- NA
})

roc_by_times <- list(tpr=tpr_by_times,
		     fpr=fpr_by_times)
## save tpr and fpr for plotting ROC 
#cat("\n--- Saving tpr and fpr for plotting ROC curve ---\n")
#saveRDS(roc_by_times, file="ext_rsf_roc_by_times.rds")

# summary for each fold and each time point
auc_matrix <- as.matrix(t(auroc))
brier_matrix <- as.matrix(t(bs))

year_labels <- round(eval_times / 365.24) 
colnames(auc_matrix) <- paste0("AUC_", year_labels, "yr")
colnames(brier_matrix) <- paste0("Brier_", year_labels, "yr")

cat("\n--- Performance Summary by Time Point ---\n")

cat("Time-dependent AUROC results:\n")
for (i in seq_along(year_labels)) {
	cat(sprintf("  At %d year(s): AUROC = %.3f\n",
		    year_labels[i], auc_matrix[i]))
}

#Brier summary
cat("Time-dependent Brier score results:\n")
for (i in seq_along(year_labels)) {
	cat(sprintf("  At %d year(s): Brier score = %.3f\n",
		    year_labels[i], brier_matrix[i]))
}

cat(sprintf("\nExternal validation dataset's integrated Brier score: %.3f\n", ibs))
cat(sprintf("\nExternal validation dataset's Wolbers C-index: %.4f\n", cindex))
