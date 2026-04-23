library(timeROC)
#library(fastcmprsk)
library(survival)
library(dplyr)
library(riskRegression)
library(prodlim)
library(randomForestSRC)
library(ggplot2)
library(tidyr)


set.seed(20250521)

eval_times <- c(365.24, 365.24*3, 365.24*5, 365.24*10)

loaded_fg <- readRDS("FG/fitted_fg.rds")
loaded_rsf <- readRDS("RSF/fitted_rsf.rds")

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
test_data_fname = "../data/preprocessed_competing_risk_validation_dataset250624.csv"
test_data <- read.csv(test_data_fname)
# Several variables are transformed to log scale
vars_to_log <- intersect(vars_to_log, names(test_data))
for (var_nm in vars_to_log) {
		new_column_nm <- paste0("log_", var_nm)
		test_data[[new_column_nm]] <- log(test_data[[var_nm]]+1)
}
test_data[factor_vars] <- lapply(test_data[factor_vars], factor)
#val_factor_vars_test <- check_n_filter_factor_vars(test_data, factor_vars)
#val_cont_vars_test <- intersect(cont_vars, names(test_data))
cat(sprintf("\n(%d, %d) Test data done\n", 
	    nrow(test_data), ncol(test_data)))
names(test_data)[names(test_data)=="milan"] <- c("Milan")
names(test_data)[names(test_data)=="uptoseven"] <- c("UTS")
names(test_data)[names(test_data)=="mt2.0"] <- c("Mt2")

# Milan
data_milan <- data.frame(ftime=test_data$ftime, 
			 fstatus=test_data$fstatus,
			 score=test_data$Milan
)

# Up-To-Seven
data_uts <- data.frame(ftime=test_data$ftime, 
		       fstatus=test_data$fstatus,
		       score=test_data$UTS
)

# Metroticket 2.0
data_mt2 <- data.frame(ftime=test_data$ftime, 
		       fstatus=test_data$fstatus,
		       score=test_data$Mt2
)

# FGR predict test data
pred.fg <- predict(loaded_fg,
		    newdata=test_data, 
		    times=eval_times)

# RSF predict test data
pred.rsf <- predict(loaded_rsf, 
		   newdata=test_data,
		   times=eval_times)

# RSF는 Event 발생 시점에서의 cumulative incidence function (cif)을 주기 때문에 
# 우리가 원하는 1, 3, 5, 10년에서의 cif를 구하기 위해서는 interpolation 해야 함
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

interp_cif <- interp_predict_risk_rsf(pred.rsf, eval_times)



# Modified DeLong's test from timeROC
roc_objects <- list()
model_names <- c("Milan", "UTS", "Mt2", "FG", "RSF")
for (i in seq_along(eval_times)) {
	t_pnt <- eval_times[i]
	time_key <- paste("t", t_pnt, sep="_")
	roc_objects[[time_key]] <- list()

	# Milan ROC
	roc_objects[[time_key]][["Milan"]] <- timeROC(T=data_milan$ftime,
						      delta=data_milan$fstatus,
						      marker=data_milan$score,
						      cause=1,
						      times=t_pnt,
						      iid=TRUE
	)

	# UTS ROC
	roc_objects[[time_key]][["UTS"]] <- timeROC(T=data_uts$ftime, 
						    delta=data_uts$fstatus,
						    marker=data_uts$score,
						    cause=1,
						    times=t_pnt,
						    iid=TRUE
	)

	# Mt2 ROC
	roc_objects[[time_key]][["Mt2"]] <- timeROC(T=data_mt2$ftime, 
						    delta=data_mt2$fstatus,
						    marker=data_mt2$score,
						    cause=1,
						    times=t_pnt,
						    iid=TRUE
	)

	# Fine Gray ROC
	current_predictions_for_time <- pred.fg[, which(eval_times==t_pnt)]
	roc_objects[[time_key]][["FG"]] <- timeROC(T=test_data$ftime,
						   delta=test_data$fstatus,
						   marker=current_predictions_for_time, 
						   cause=1,
						   times=t_pnt,
						   iid=TRUE
	)

	# RSF ROC
	colname_interp <- paste("t", t_pnt, sep="_")
	current_pred_at_time <- interp_cif[[colname_interp]]
	roc_objects[[time_key]][["RSF"]] <- timeROC(T=test_data$ftime,
						    delta=test_data$fstatus,
						    marker=current_pred_at_time,  
						    cause=1,
						    times=t_pnt,
						    iid=TRUE
	)
}

# Perform pairwise comparison for each time point
for (i in seq_along(eval_times)) {
	t_pnt <- eval_times[i]
	time_key <- paste("t", t_pnt, sep="_")
	time_header <- paste("t", t_pnt, sep="=")
	cat("\n", strrep("=", 50), "\n")
	cat("Time point:", t_pnt, "days\n")
	cat(strrep("=", 50), "\n")

	cat("AUROC values:\n")
	for (model in model_names) {
#		print(roc_objects[[time_key]][[model]])
		auc_val <- roc_objects[[time_key]][[model]]$AUC_1[2]
		cat(sprintf("%s: %.4f\n", model, auc_val))
	}

	# Pairwise comparison (Modified DeLong's test):\n")
	for (j in 1:(length(model_names)-1)) {
		for (k in (j+1):length(model_names)) {
			model1 <- model_names[j]
			model2 <- model_names[k]

			comparison_result <- compare(roc_objects[[time_key]][[model1]],
						     roc_objects[[time_key]][[model2]])

			cat(sprintf("\n%s vs %s:\n", model1, model2))
			cat(sprintf("   p-value: %.4f\n", 
				as.numeric(comparison_result$p_values_AUC_1[time_header])))
#			if (comparison_result%p.value < 0.001) {
#				cat("   Significance: ***\n")
#			} else if (comparison_result$p.value < 0.01) {
#				cat("   Significance: **\n")
#			} else if (comparison_result$p.value <0.05) {
#				cat(" Significance: *\n")
#			} else {
#				cat(" No significance\n")
#			}
		}
	}
}


# Display labels mapping
display_names <- c(Milan="Milan", UTS="UtS", Mt2="Mt2.0", FG="FG", RSF="RSF")
time_labels   <- c("1 year", "3 year", "5 year", "10 year")

# Build long-format p-value table
pval_long_all <- lapply(seq_along(eval_times), function(i) {
  t_pnt      <- eval_times[i]
  time_key   <- paste("t", t_pnt, sep="_")
  time_header <- paste("t", t_pnt, sep="=")
  n <- length(model_names)

  pval_mat <- matrix(NA, n, n, dimnames=list(model_names, model_names))
  for (j in 1:(n-1)) for (k in (j+1):n) {
    pv <- as.numeric(compare(roc_objects[[time_key]][[model_names[j]]],
                             roc_objects[[time_key]][[model_names[k]]])$p_values_AUC_1[time_header])
    pval_mat[j, k] <- pval_mat[k, j] <- pv
  }

  as.data.frame(pval_mat) |>
    mutate(row = model_names) |>
    pivot_longer(-row, names_to="col", values_to="pvalue") |>
    mutate(time = time_labels[i])
})

plot_df <- bind_rows(pval_long_all) |>
  mutate(
    time    = factor(time, levels=time_labels),
    row     = factor(display_names[row], levels=rev(display_names)),
    col     = factor(display_names[col], levels=display_names),
    sig     = case_when(is.na(pvalue) ~ "",
                        pvalue < 0.01 ~ paste0(sprintf("%.2f", pvalue), "**"),
                        pvalue < 0.05 ~ paste0(sprintf("%.2f", pvalue), "*"),
                        TRUE          ~ sprintf("%.2f", pvalue))
  )

ggplot(plot_df, aes(col, row, fill=pvalue)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=sig), size=3.5) +
  scale_fill_distiller(palette="Blues", direction=1, limits=c(0,1),
                       na.value="white", name="p-value") +
  scale_y_discrete(limits=c("Milan","UtS","Mt2.0","FG","RSF")) +
  facet_wrap(~time, nrow=1) +
  labs(tag = "(A)") + 
  theme_minimal(base_size=11) +
  theme(
    axis.title        = element_blank(),
    panel.grid        = element_blank(),
    strip.text        = element_text(face="bold", size=13),
    axis.text.x       = element_text(angle=45, hjust=1, color="black", size="12"),
    axis.text.y       = element_text(color="black", size="12"),
    axis.ticks        = element_line(color="black", linewidth=0.7),
    axis.ticks.length = unit(4, "pt"),
    legend.title      = element_text(size=12),                       # legend 제목
    legend.text       = element_text(size=11),                       # colorbar 숫자
    axis.line         = element_line(color="black", linewidth=0.7),  # 축선 추가
    panel.spacing     = unit(0.8, "lines"),
    plot.tag          = element_text(size=12, face="bold"),
    plot.tag.position = c(0.015, 0.95), # left upper (x = 0, y = 1) 
    aspect.ratio      = 7/4   # coord_fixed 대신 세로로 길게 (5행 4열 비율)
  )

ggsave("Table3-heatmap-auroc.pdf", width=1075/72, height=494/72)
