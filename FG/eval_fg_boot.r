## eval_fg_boot.R
## Fine-Gray 모델: 이미 fitted된 모델로 bootstrapped test data에서 metrics 계산
## bootstrap_indices.rds (test data 기준) 공유 → 다른 모델과 paired comparison 가능

library(dplyr)
library(survival)
library(riskRegression)
library(timeROC)
library(pracma)
library(pec)

## ── Config ────────────────────────────────────────────────────────────────────
eval_times  <- c(365.24, 365.24*3, 365.24*5, 365.24*10)
year_labels <- round(eval_times / 365.24)

vars_to_log <- c("afp", "pivka", "tbil", "pt_inr",
                 "lr_5_size_ct", "lr_tr_v_size_whole_ct",
                 "lr_tr_v_size_enhancing_ct", "lr_3_size_ct",
                 "lr_4_size_ct", "lr_m_size_ct", "lr_tr_nv_size_ct",
                 "lr_tr_e_size_ct")
factor_vars <- c("sex", "diabetes", "ckd_hd",
                 "hbv", "hcv", "image_lc", "ald", "masld", "ldlt",
                 "pre_tx", "tace", "rfa", "resection", "rtx")

optim_vars_for_ext <- c("log_lr_5_size_ct", "log_lr_tr_v_size_enhancing_ct",
                        "log_lr_m_size_ct", "log_lr_tr_nv_size_ct",
                        "log_afp", "log_pivka", "log_pt_inr", "sex", "ckd_hd")

## ── Load data ─────────────────────────────────────────────────────────────────
test_data_fname <- "../../data/preprocessed_competing_risk_validation_dataset250624.csv"
test_data       <- read.csv(test_data_fname)

vars_to_log_test <- intersect(vars_to_log, names(test_data))
for (var_nm in vars_to_log_test)
  test_data[[paste0("log_", var_nm)]] <- log(test_data[[var_nm]] + 1)
test_data[intersect(factor_vars, names(test_data))] <-
  lapply(test_data[intersect(factor_vars, names(test_data))], factor)

final_vars    <- intersect(optim_vars_for_ext, names(test_data))
test_covs     <- test_data %>% select(all_of(final_vars))

## ── Load fitted model ─────────────────────────────────────────────────────────
fit <- readRDS("fitted_fg.rds")
cat("Loaded: fitted_fg.rds\n")

## ── Load shared bootstrap indices ─────────────────────────────────────────────
boot_indices <- readRDS("../bootstrap_indices.rds")
N_BOOT       <- length(boot_indices)
cat(sprintf("Loaded %d bootstrap indices\n", N_BOOT))

## ── FG formula (needed for pec::cindex) ───────────────────────────────────────
formula_str <- paste("Hist(ftime, fstatus) ~",
                     paste(final_vars, collapse=" + "))
fgmodel_formula <- as.formula(formula_str)

## ── Helper: compute all metrics on a bootstrapped test dataset ────────────────
compute_fg_metrics <- function(d) {

  covs <- d %>% select(all_of(final_vars))

  ## C-index (pec::cindex)
  cindex <- tryCatch({
    pec::cindex(fit,
                formula = fgmodel_formula,
                data    = d,
                cause   = 1)$AppCindex$FGR
  }, error = function(e) NA_real_)

  ## Predict at eval_times
  pred <- tryCatch(
    predict(fit, newdata=covs, times=eval_times),
    error = function(e) NULL
  )
  if (is.null(pred)) return(NULL)

  ## AUROC + Brier
  sc <- tryCatch(
    Score(list("FG" = pred),
          formula = Hist(ftime, fstatus) ~ 1,
          data    = d,
          times   = eval_times,
          metrics = c("AUC", "Brier"),
          cause   = 1),
    error = function(e) NULL
  )
  auroc <- if (!is.null(sc)) sc$AUC$score[model=="FG",   AUC]   else rep(NA_real_, length(eval_times))
  bs    <- if (!is.null(sc)) sc$Brier$score[model=="FG", Brier] else rep(NA_real_, length(eval_times))

  ## IBS via trapz over all event times in bootstrap sample
  ibs <- tryCatch({
    event_times  <- sort(unique(subset(d, fstatus %in% c(1,2))$ftime))
    pred_ibs     <- predict(fit, newdata=covs, times=event_times)
    bs_all       <- Score(list("FG" = pred_ibs),
                          formula = Hist(ftime, fstatus) ~ 1,
                          data    = d,
                          times   = event_times,
                          metrics = "Brier",
                          cause   = 1)
    dt  <- bs_all$Brier$score[model == "FG", .(times, Brier)]
    trapz(dt$times, dt$Brier) / max(event_times)
  }, error = function(e) NA_real_)

  result <- list(cindex = cindex, ibs = ibs)
  for (i in seq_along(eval_times)) {
    result[[paste0("auroc_", year_labels[i], "yr")]] <- auroc[i]
    result[[paste0("brier_", year_labels[i], "yr")]] <- bs[i]
  }
  result
}

## ── Run bootstrap ─────────────────────────────────────────────────────────────
cat("Running bootstrap on test data...\n")
boot_results <- vector("list", N_BOOT)

for (b in seq_len(N_BOOT)) {
  if (b %% 10 == 0) cat(sprintf("  iter %d / %d\n", b, N_BOOT))
  d_boot        <- test_data[boot_indices[[b]], ]
  boot_results[[b]] <- compute_fg_metrics(d_boot)
}

## ── Summarise ─────────────────────────────────────────────────────────────────
boot_results_clean <- Filter(Negate(is.null), boot_results)
summary_df <- do.call(rbind, lapply(boot_results_clean, function(r) {
  as.data.frame(r, stringsAsFactors=FALSE)
}))

cat("\n--- Bootstrap summary (mean ± SD) ---\n")
for (m in names(summary_df)) {
  vals <- summary_df[[m]]
  cat(sprintf("  %-25s %.4f ± %.4f\n", m,
              mean(vals, na.rm=TRUE), sd(vals, na.rm=TRUE)))
}

## ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(list(boot_results = boot_results_clean,
             summary_df   = summary_df,
             eval_times   = eval_times,
             year_labels  = year_labels),
        "fg_boot_results.rds")

cat("\nSaved → fg_boot_results.rds\n")
