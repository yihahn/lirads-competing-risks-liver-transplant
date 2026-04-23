library(caret)
library(dplyr)
library(fastcmprsk)
library(prodlim)
library(survival)
library(riskRegression)
library(timeROC)
library(pracma)


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


set.seed(20250521)

############ MUST BE EDITED ########################################
eval_times <- c(365.24, 365.24*3, 365.24*5, 365.24*10)

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
               'alb', 'log_tbil', 'cr', 'log_pt_inr', 'na', 'log_crp',# lab test
               'age_lt', 'bmi' # demographic
)

# variables to be transformed log scale
vars_to_log <- c('afp', 'pivka', 'tbil', 'pt_inr', 'crp',
		 'lr_5_size_ct', 'lr_tr_v_size_whole_ct',
		 'lr_tr_v_size_enhancing_ct', 'lr_3_size_ct',
		 'lr_4_size_ct', 'lr_m_size_ct', 'lr_tr_nv_size_ct',
		 'lr_tr_e_size_ct'
)

# conservative
#optim_vars_for_ext <- c("log_lr_5_size_ct", "log_lr_tr_v_size_enhancing_ct", "log_lr_m_size_ct", 
#			"log_lr_tr_nv_size_ct", "log_afp", "log_pivka", "log_pt_inr", "ckd_hd", "tace")
# moderate
#optim_vars_for_ext <- c("log_lr_5_size_ct", "log_lr_tr_v_size_enhancing_ct", "log_lr_m_size_ct", 
#			"log_lr_tr_nv_size_ct", "log_afp", "log_pivka", "log_pt_inr", "sex", "ckd_hd", "tace")
# liberal
optim_vars_for_ext <- c("log_lr_5_size_ct", "log_lr_tr_v_size_enhancing_ct", "log_lr_m_size_ct", 
			"log_lr_tr_nv_size_ct", "log_afp", "log_pivka", "log_pt_inr", "sex", "ckd_hd", "ald")

# aic
#optim_vars_for_ext <- c('lr_3_ct', 'lr_m_ct', 'lr_tr_nv_ct', 'log_lr_5_size_ct' , 'log_lr_tr_v_size_enhancing_ct', 
#			'log_lr_3_size_ct', 'log_lr_m_size_ct',  'log_lr_tr_nv_size_ct', 'log_afp', 
#			'log_pivka', 'log_pt_inr', 'sex', 'ckd_hd', 'tace', 'resection')
# bic
#optim_vars_for_ext <- c('lr_3_ct', 'lr_m_ct', 'log_lr_5_size_ct' , 'log_lr_tr_v_size_enhancing_ct', 
#			'log_lr_3_size_ct', 'log_lr_m_size_ct',  'log_lr_tr_nv_size_ct', 'log_afp', 
#			'log_pivka', 'log_pt_inr', 'sex', 'ckd_hd', 'tace', 'resection')


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

# Start to do regression
formula_str <- paste("Hist(ftime, fstatus) ~", 
		     paste(final_vars, collapse=" + "))
cat("\n\n Final formula of Fine-Gray model\n")
print(formula_str)
fgmodel_formula <- as.formula(formula_str)

# cv unbiased predicted risk for calibration
time_labels <- paste0("t_", eval_times)
create_cv_folds <- function(data, k=5) {
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

	cv_fit <- FGR(fgmodel_formula,
		      data=cv_train_data,
		      cause=1
	)

	# Wolbers' c index
	cindex_fold <- suppressMessages(pec::cindex(list(fg=cv_fit),
						    formula=fgmodel_formula,
						    data=cv_test_data,
						    cause=1))
	c_indices[i] <- cindex_fold$AppCindex$fg
	cv_test_covariates <- cv_test_data %>%
		select(all_of(final_vars))

	# AUROC & Brier from riskRegression 
	cv_preds <- predict(cv_fit,
			    newdata=cv_test_covariates,
			    times=eval_times
	)

	score <- Score(list("FG"=cv_preds),
		       formula=Hist(ftime, fstatus) ~ 1,
		       data=cv_test_data,
		       times=eval_times,
		       metrics=c("AUC", "Brier"),
		       cause=1 # interest event label
		       )

	auc_per_fold[[i]] <- score$AUC$score %>%
		filter(model=="FG") %>%
		select(AUC) %>%
		pull()

	brier_per_fold[[i]] <- score$Brier$score %>%
		filter(model=="FG") %>%
		select(Brier) %>%
		pull()

	# Integrated Brier score by trapz
	total_event_time_ibs <- subset(cv_test_data, fstatus %in% c(1, 2))$ftime
	total_event_time_ibs <- sort(unique(total_event_time_ibs))
	pred_ibs <- predict(cv_fit, 
			    newdata=cv_test_covariates, 
			    times=total_event_time_ibs)

	brier_scores <- Score(list("FG"=pred_ibs),
			      formula=Hist(ftime, fstatus) ~ 1,
			      data=cv_test_data,
			      times=total_event_time_ibs,
			      metrics="Brier",
			      cause=1
	)
	brier_scores_dt <- brier_scores$Brier$score
	fg_brier_scores <- brier_scores_dt[model=="FG", .(times, Brier)]
	ibs_fg <- trapz(fg_brier_scores$times, fg_brier_scores$Brier)
	ibs_per_fold[[i]] <- ibs_fg / max(total_event_time_ibs)


	cv_preds <- predict(cv_fit,
			    newdata=cv_test_covariates,
			    times=eval_times
	)
	pred_df <- as.data.frame(cv_preds)
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
#print("saving fg_cv_data_for_calibration.rds")
#saveRDS(cv_calibrate_data, "fg_cv_data_for_calibration.rds")
#print(str(cv_calibrate_data))

test_covariates <- test_data %>%
	select(all_of(final_vars))

fit <- FGR(fgmodel_formula, 
	   data=train_data,
	   cause=1
)

#saveRDS(fit, "fitted_fg.rds")

# Wolbers C-index 
wolbers_cindex_out <- pec::cindex(fit, formula=fgmodel_formula, 
				  data=test_data, cause=1)
wolbers_cindex <- wolbers_cindex_out$AppCindex$FGR

# predict test data
pred <- predict(fit, 
	        newdata=test_covariates, 
	        times=eval_times)


# AUROC & Brier from riskRegression 
score <- Score(list("FG"=pred),
	       formula=Hist(ftime, fstatus) ~ 1,
	       data=test_data,
	       times=eval_times,
	       metrics=c("AUC", "Brier"),
	       cause=1 # interest event label
	       )


auroc <- score$AUC$score %>%
	filter(model=="FG") %>%
	select(AUC) %>%
	pull()

bs <- score$Brier$score %>%
	filter(model=="FG") %>%
	select(Brier) %>%
	pull()

# Integrated Brier score by trapz
total_event_times_ibs <- subset(test_data, fstatus %in% c(1, 2))$ftime
unique_event_times_ibs <- sort(unique(total_event_times_ibs))
pred <- predict(fit, 
	        newdata=test_covariates, 
	        times=unique_event_times_ibs)

brier_scores <- Score(list("FG"=pred),
		      formula=Hist(ftime, fstatus) ~ 1,
		      data=test_data,
		      times=unique_event_times_ibs,
		      metrics="Brier",
		      cause=1
		      )
brier_scores_dt <- brier_scores$Brier$score
fg_brier_scores <- brier_scores_dt[model=="FG", .(times, Brier)]
ibs_fg <- trapz(fg_brier_scores$times, fg_brier_scores$Brier)
ibs <- ibs_fg/max(total_event_times_ibs)
# tpr and fpr to plot ROC curve by timeROC
pred <- predict(fit, 
	        newdata=test_covariates, 
	        times=eval_times)

all_pred_by_time <- matrix(NA, 
			   nrow=nrow(test_data),
			   ncol=length(eval_times),
			   dimnames=list(NULL, paste0("t=", eval_times)))

for (i in seq_along(eval_times)) {
	t_pnt <- eval_times[i]
	current_predictions_for_time <- pred[, which(eval_times==t_pnt)]
	all_pred_by_time[, i] <- current_predictions_for_time
	roc_obj <- timeROC(T=test_data$ftime, 
			   delta=test_data$fstatus,
			   marker=current_predictions_for_time,
			   cause=1,
			   times=t_pnt,
			   ROC=TRUE
			   )
	t_pnt_colname <- paste("t=", t_pnt, sep="")
	current_tpr <- roc_obj$TP[, t_pnt_colname]
	current_fpr <- roc_obj$FP_1[, t_pnt_colname]
	tpr_by_times[[i]] <- c(tpr_by_times[[i]], list(current_tpr))
	fpr_by_times[[i]] <- c(fpr_by_times[[i]], list(current_fpr))
}

# save predictions at eval_time for DeLong's test
#saveRDS(as.data.frame(all_pred_by_time), 
#	file="ext_fg_pred_by_evaltimes.rds")

# t=365.24인 5 fold tpr 모두 접근
roc_by_times <- list(tpr=tpr_by_times,
		     fpr=fpr_by_times
)

## save tpr and fpr for plotting ROC 
#cat("\n--- Saving tpr and fpr for plotting ROC curve ---\n")
#saveRDS(roc_by_times, file="ext_fg_roc_by_times.rds")

#print(paste0("Wolbers' C index", wolbers_cindex))
#print(paste0("AUROC", auroc))
#print(paste0("Brier score", bs))
#print(paste0("integrated Brier score", ibs))

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

cat(sprintf("\nExternal validation dataset's integrated Brier score: %.4f\n", ibs))
cat(sprintf("\nExternal validation dataset's Wolbers C-index: %.4f\n\n", wolbers_cindex))
