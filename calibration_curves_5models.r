library(timeROC)
library(fastcmprsk)
library(survival)
library(dplyr)
library(riskRegression)
library(prodlim)
library(randomForestSRC)
library(ggplot2)
library(gridExtra)
library(pec)

set.seed(20250521)

eval_times <- c(365.24, 365.24*3, 365.24*5, 365.24*10)

loaded_fg <- readRDS("FG/fitted_fg.rds")
loaded_rsf <- readRDS("RSF/fitted_rsf.rds")

factor_vars <- c('sex', 'diabetes', 'ckd_hd', # baseline 
		 'hbv', 'hcv', 'image_lc', 'ald', 'masld', 'ldlt',  # cause of liver disease
		 'pre_tx', 'tace', 'rfa', 'resection', 'rtx' # previous HCC treatment
)

#vars_to_log <- c("afp", "pivka", "tbil", "cr", "pt_inr", "crp")
# variables to be transformed log scale
vars_to_log <- c('afp', 'pivka', 'tbil', 'pt_inr', 'crp',
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
cat(sprintf("\n(%d, %d) Test data done\n", 
	    nrow(test_data), ncol(test_data)))
names(test_data)[names(test_data)=="milan"] <- c("Milan")
names(test_data)[names(test_data)=="uptoseven"] <- c("UTS")
names(test_data)[names(test_data)=="mt2.0"] <- c("Mt2")

# FGR predict test data
pred.fg <- predict(loaded_fg,
		    newdata=test_data, 
		    times=eval_times)

# RSF predict test data
pred.rsf <- predict(loaded_rsf, 
		   newdata=test_data)

# RSF interpolation for evaluation times
time_vec <- pred.rsf$time.interest
interp_cif_matrix <- matrix(NA, nrow=nrow(pred.rsf$cif), ncol=length(eval_times))
for (k in 1:nrow(pred.rsf$cif)) {
	individual_cif <- pred.rsf$cif[k, , 1]
	for (i in seq_along(eval_times)) {
		# Use approx for interpolation for the current individual
		interpolated_cif_at_time <- approx(x=time_vec, 
						   y=individual_cif, 
						   xout=eval_times[i],
						   rule=2)$y
		interp_cif_matrix[k, i] <- interpolated_cif_at_time
	}
}
colnames(interp_cif_matrix) <- paste("t", eval_times, sep="_")
interp_cif <- as.data.frame(interp_cif_matrix)

# Function to calculate calibration metrics
calculate_calibration <- function(predicted_risk, time_to_event, event_indicator, 
				  eval_time, n_groups=10) {

	# Test data
	# Test data에서 NA 값 제거
	notnull_idx <- !is.na(predicted_risk) & !is.na(time_to_event) & !is.na(event_indicator)

	if (sum(notnull_idx) == 0) {
		warning("No notnull observation found")
		return(data.frame())
	}

	# notnull sample만 사용
	predicted_risk <- predicted_risk[notnull_idx]
	time_to_event <- time_to_event[notnull_idx]
	event_indicator <- event_indicator[notnull_idx]


	# analysis sample in test data were filtered at eval_time
	analysis_mask <- rep(FALSE, length(time_to_event))
	for (j in seq_along(time_to_event)) {
		if (time_to_event[j] >= eval_time) {
			analysis_mask[j] <- TRUE
		} else {
			if (event_indicator[j]==1) {
				analysis_mask[j] <- TRUE
			} else {
				analysis_mask[j] <- FALSE
			}
		}
	}


	if (sum(analysis_mask) == 0) {
		warning(paste("No valid samples for analysis at eval_time", eval_time))
		return(data.frame())
	}
	# use only analysis_mask=TRUE at eval_time
	predicted_risk <- predicted_risk[analysis_mask]
	time_to_event <- time_to_event[analysis_mask]
	event_indicator <- event_indicator[analysis_mask]

	# re-define event status at eval_time
	event_at_eval_time <- ifelse(time_to_event < eval_time & event_indicator == 1, 1, 0)

	#cat("Test data, len(predicted_risk) =", length(predicted_risk),
	#    "len(event_at_eval_time) =", length(event_at_eval_time), "\n")
	#print(table(event_at_eval_time))

	# Check if predicted_risk is binary (for Milan, UTS, Mt2)
	unique_risks <- unique(predicted_risk)
	unique_risks <- unique_risks[!is.na(unique_risks)]
	
	if (length(unique_risks) <= 2) {
		# For binary predictors, create groups based on unique values
		risk_groups <- as.numeric(as.factor(predicted_risk))
		n_groups <- length(unique_risks)
	} else {
		# For continuous predictors, use quantile-based groups
		quantiles <- quantile(predicted_risk, probs=seq(0, 1, 1/n_groups),
				      na.rm=TRUE)
		# Check if quantiles are unique
		if (length(unique(quantiles)) < length(quantiles)) {
			# If quantiles are not unique, create groups based on unique values
			sorted_unique <- sort(unique_risks)
			if (length(sorted_unique) <= n_groups) {
				# Use all unique values as breaks
				breaks <- c(-Inf, sorted_unique[-1] - 0.001, Inf)
				risk_groups <- cut(predicted_risk, breaks=breaks, 
						   include.lowest=TRUE, labels=FALSE)
				n_groups <- length(sorted_unique)
			} else {
				# Use equally spaced breaks
				breaks <- seq(min(predicted_risk, na.rm=TRUE), 
					      max(predicted_risk, na.rm=TRUE), 
					      length.out=n_groups + 1)
				breaks[1] <- breaks[1] - 0.001  # Ensure all values are included
				risk_groups <- cut(predicted_risk, breaks=breaks, 
						   include.lowest=TRUE, labels=FALSE)
			}
		} else {
			# Normal quantile-based grouping
			risk_groups <- cut(predicted_risk, breaks=quantiles,
					   include.lowest=TRUE, labels=FALSE)
		}
	}

	# Calculate observed and expected rates for each group
	actual_n_groups <- max(risk_groups, na.rm = TRUE)
	calib_data <- data.frame(risk_group=1:actual_n_groups,
				 observed=numeric(actual_n_groups),
				 expected=numeric(actual_n_groups),
				 n_patients=numeric(actual_n_groups)
	)

	for (i in 1:actual_n_groups) {
		group_mask <- risk_groups == i & !is.na(risk_groups)
		if (sum(group_mask) == 0) next

		# Expected risk (mean predicted risk in group)
		calib_data$expected[i] <- mean(predicted_risk[group_mask], na.rm=TRUE)
		calib_data$n_patients[i] <- sum(group_mask)
		calib_data$observed[i] <- mean(event_at_eval_time[group_mask], na.rm=TRUE)

	}
	
	# Remove groups with no patients
	calib_data <- calib_data[calib_data$n_patients > 0, ]
	return(calib_data)
}


create_calibration_plots <- function() {
  model_names    <- c("Milan", "UTS", "Mt2", "FG", "RSF")
  display_names  <- c(Milan="Milan", UTS="UtS", Mt2="Mt2.0", FG="FG", RSF="RSF")
  model_order    <- c("FG", "RSF", "Milan", "UtS", "Mt2.0")
  year_labels    <- c("1 Year", "3 Years", "5 Years", "10 Years")

  # 전체 데이터 구성
  all_data <- data.frame()
  for (i in seq_along(eval_times)) {
    t_pnt <- eval_times[i]
    for (model in model_names) {
      predicted_risk <- switch(model,
        Milan = test_data$Milan,
        UTS   = test_data$UTS,
        Mt2   = test_data$Mt2,
        FG    = pred.fg[, i],
        RSF   = interp_cif[[paste("t", t_pnt, sep="_")]]
      )
      calib_data <- calculate_calibration(
        predicted_risk  = predicted_risk,
        time_to_event   = test_data$ftime,
        event_indicator = test_data$fstatus,
        eval_time       = t_pnt,
        n_groups        = ifelse(model %in% c("Milan","UTS","Mt2"), 2, 5)
      )
      calib_data$model <- display_names[model]
      calib_data$Year  <- year_labels[i]
      all_data <- rbind(all_data, calib_data)
    }
  }

  all_data$model <- factor(all_data$model, levels=model_order)
  all_data$Year  <- factor(all_data$Year,  levels=year_labels)

  ggplot(all_data, aes(x=expected, y=observed, color=model)) +
    geom_line(aes(group=model), alpha=0.6, linewidth=0.6) +
    geom_point(aes(size=n_patients), alpha=0.6) +
    geom_abline(intercept=0, slope=1, linetype="dashed",
                color="gray50", linewidth=0.5) +
    scale_x_continuous(limits=c(-0.05, 1.05), breaks=seq(0, 1, 0.1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.05, 1.05), breaks=seq(0, 1, 0.1), expand=c(0, 0)) +
    scale_size_continuous(range=c(1, 5), name="N") +
    facet_wrap(~Year, nrow=2, ncol=2, scales="fixed") +
    labs(
      title = "Calibration curves of external validation cohort across Years",
      x     = "Predicted Risk",
      y     = "Observed Risk",
      color = "Model"
    ) +
    guides(
      color = guide_legend(
        position       = "bottom",
        title.position = "left", title.hjust = 0.5,
        override.aes   = list(size=3), nrow=1
      ),
      size = guide_legend(
        position       = "inside",
        title.position = "top", title.hjust = 0.5,
        ncol=1
      )
    ) +
    theme_minimal() +
    theme(
      panel.background      = element_blank(),
      panel.grid.major      = element_line(color="gray90", linewidth=0.3),
      panel.grid.minor      = element_blank(),
      panel.border          = element_rect(color="black", fill=NA, linewidth=0.8),
      axis.line             = element_line(color="black", linewidth=0.5),
      axis.ticks            = element_line(color="black", linewidth=0.5),
      axis.ticks.length     = unit(0.1, "cm"),
      axis.title            = element_text(size=12),
      axis.text             = element_text(size=10),
      plot.title            = element_text(hjust=0.5, face="bold", size=16),
      strip.text            = element_text(face="bold", size=13),
      strip.background      = element_rect(fill="gray95", color="black", linewidth=0.5),
      legend.position       = "bottom",
      legend.position.inside = c(0.03, 0.92),   # N legend 왼쪽 상단
      legend.justification.inside = c(0, 1),
      legend.title          = element_text(face="bold", size=11),
      legend.text           = element_text(size=10),
      legend.background     = element_rect(fill="white", color="gray70",
                                           linewidth=0.3, linetype="solid"),
      legend.margin         = margin(3, 5, 3, 5)
    )
}

calibration_plot <- create_calibration_plots()

pdf("calibration_curves_5models.pdf", width=12, height=10)
print(calibration_plot)
dev.off()

tiff("calibration_curves_5models.tiff", width=12, height=10,
     units="in", res=300, bg="white", compression="lzw")
print(calibration_plot)
dev.off()


# Calculate calibration metrics (Brier Score, Calibration slope, C-index)
calculate_calibration_metrics <- function() {
	metrics_df <- data.frame()
	model_names <- c("Milan", "UTS", "Mt2", "FG", "RSF")

	for (i in seq_along(eval_times)) {
		t_pnt <- eval_times[i]

		for (model in model_names) {
			if (model == "Milan") {
				predicted_risk <- test_data$Milan
			} else if (model == "UTS") {
				predicted_risk <- test_data$UTS
			} else if (model == "Mt2") {
				predicted_risk <- test_data$Mt2
			} else if (model == "FG") {
				predicted_risk <- pred.fg[, i]
			} else if (model == "RSF") {
				colname_interp <- paste("t", t_pnt, sep="_")
				predicted_risk <- interp_cif[[colname_interp]]
			}

			# Calculate calibration metrics
			calib_data <- calculate_calibration(predicted_risk=predicted_risk,
							    time_to_event=test_data$ftime,
							    event_indicator=test_data$fstatus,
							    eval_time=t_pnt,
							    n_groups=ifelse(model %in% c("Milan", "UTS", "Mt2"), 2, 5)
			)

			# Calculate calibration slope (slope of observed vs expected)
			if (nrow(calib_data) > 1) {
				calib_slope <- coef(lm(observed ~ expected, data=calib_data))[2]
				calib_intercept <- coef(lm(observed ~ expected, data=calib_data))[1]
				
				# Calculate mean squared error
				mse <- mean((calib_data$observed - calib_data$expected)^2, na.rm=TRUE)
				
				# Calculate mean absolute error
				mae <- mean(abs(calib_data$observed - calib_data$expected), na.rm=TRUE)
				# Expected Calibration Error (ECE)
				total_n_patients <- sum(calib_data$n_patients)
				ece <- sum((calib_data$n_patients/total_n_patients) * abs(calib_data$observed - calib_data$expected))

			} else {
				calib_slope <- NA
				calib_intercept <- NA
				mse <- NA
				mae <- NA
				ece <- NA
			}

			# Add to metrics dataframe
			metrics_df <- rbind(metrics_df, data.frame(Model=model,
								   Time_Point=paste0(round(t_pnt/365.24, 1), " ", ifelse(round(t_pnt/365.24, 1)==1, "year", "years")),
								   Calibration_Slope=calib_slope,
								   Calibration_Intercept=calib_intercept,
								   MSE=mse,
								   MAE=mae,
								   ECE=ece
								   ))
		}
	}

	return(metrics_df)
}

# Calculate and display calibration metrics
calibration_metrics <- calculate_calibration_metrics()
print("Calibration Metrics Summary:")
print(calibration_metrics)

# Save calibration metrics to CSV
write.csv(calibration_metrics, "calibration_metrics_5models.csv", row.names = FALSE)

# Create a summary table for each time point
for (i in seq_along(eval_times)) {
	t_pnt <- eval_times[i]
	time_label <- paste0(round(t_pnt/365.24, 1), " ", 
			     ifelse(round(t_pnt/365.24, 1)==1, "year", "years"))

	cat("\n", strrep("=", 60), "\n")
	cat("Calibration Summary for Time Point:", time_label, "\n")
	cat(strrep("=", 60), "\n")

	time_metrics <- calibration_metrics[calibration_metrics$Time_Point == time_label, ]
	for (j in 1:nrow(time_metrics)) {
		cat(sprintf("%-8s: Slope=%.4f, Intercept=%.4f, MSE=%.4f, MAE=%.4f\n",
			    time_metrics$Model[j],
			    time_metrics$Calibration_Slope[j],
			    time_metrics$Calibration_Intercept[j],
			    time_metrics$MSE[j],
			    time_metrics$MAE[j])
		)
	}
}

cat("\n=== Calibration Analysis Complete ===\n")
cat("Perfect calibration: Slope=1.0, Intercept=0.0\n")
cat("Better calibration: Lower MSE and MAE values\n")
cat("Calibration plots saved as grid arrangement\n")
cat("Calibration metrics saved to 'calibration_metrics_5models.csv'\n")
