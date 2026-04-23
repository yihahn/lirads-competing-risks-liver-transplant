library(dplyr)
library(survival)
library(fastcmprsk)
library(riskRegression)
library(timeROC)
library(pracma)


source("wolbers_cindex_cr_claude.r")

set.seed(20250521)

eval_times <- c(365.24, 365.24*3, 365.24*5, 365.24*10)
tpr_by_times <- vector("list", length(eval_times))
fpr_by_times <- vector("list", length(eval_times))
names(tpr_by_times) <- paste0("t=", eval_times)
names(fpr_by_times) <- paste0("t=", eval_times)

test_data_fname = "../data/preprocessed_competing_risk_validation_dataset250624.csv"
test_data <- read.csv(test_data_fname)
data_cr <- data.frame(ftime=test_data$ftime,
		      fstatus=test_data$fstatus,
		      uts=test_data$uptoseven
)
data_cr$uts <- as.numeric(data_cr$uts)

wolbers_cindex <- wolbers_c_index_competing(times=data_cr$ftime, 
				    events=data_cr$fstatus, 
				    predictions=data_cr$uts,
				    event_of_interest=1
)
cindex <- wolbers_cindex$c_index

score <- Score(list("UTS"=data_cr$uts),
	       formula=Hist(ftime, fstatus) ~ 1,
	       data=data_cr,
	       times=eval_times,
	       metrics=c("AUC", "Brier"),
	       cause=1
)
auroc <- score$AUC$score %>%
	filter(model=="UTS") %>%
	select(AUC) %>%
	pull()

bs <- score$Brier$score %>%
	filter(model=="UTS") %>%
	select(Brier) %>%
	pull()

# Integrated Brier score by trapz
total_event_time_ibs <- subset(data_cr, fstatus %in% c(1, 2))$ftime
brier_scores <- Score(list("UTS"=data_cr$uts),
		      formula=Hist(ftime, fstatus) ~ 1,
		      data=data_cr,
		      times=total_event_time_ibs,
		      metrics="Brier",
		      cause=1
)
brier_scores_dt <- brier_scores$Brier$score
uts_brier_scores <- brier_scores_dt[model=="UTS", .(times, Brier)]
ibs_uts <- trapz(uts_brier_scores$times, uts_brier_scores$Brier)
ibs <- ibs_uts/max(total_event_time_ibs)


for (i in seq_along(eval_times)) {
	t_pnt <- eval_times[i]
	roc_obj <- timeROC(T=data_cr$ftime, 
			   delta=data_cr$fstatus,
			   marker=data_cr$uts,
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
# t=365.24인 5 fold tpr 모두 접근
roc_by_times <- list(tpr=tpr_by_times,
		     fpr=fpr_by_times
)

# save tpr and fpr for plotting ROC 
cat("\n--- Saving tpr and fpr for plotting ROC curve ---\n")
saveRDS(roc_by_times, file="ext_uts_roc_by_times.rds")

# summary for each fold and each time point
auc_matrix <- as.matrix(t(auroc))
brier_matrix <- as.matrix(t(bs))

year_labels <- round(eval_times / 365.24) 
colnames(auc_matrix) <- paste0("AUC_", year_labels, "yr")
colnames(brier_matrix) <- paste0("Brier_", year_labels, "yr")

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

cat(sprintf("External validation dataset's integrated Brier score: %.4f\n", ibs))
cat(sprintf("External validation datasets's Wolbers' C-index: %.4f\n", cindex))
