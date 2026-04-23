## eval_rsf_boot.R
## RSF 모델: fitted_rsf.rds 로드 후 bootstrapped test data에서 metrics 계산
## (train data 재학습 없음, test data만 bootstrap)

library(dplyr)
library(randomForestSRC)
library(prodlim)
library(survival)
library(riskRegression)
library(timeROC)
library(pracma)
library(pec)

## ── Config ────────────────────────────────────────────────────────────────────
eval_times  <- c(365.24, 365.24*3, 365.24*5, 365.24*10)
year_labels <- round(eval_times / 365.24)

vars_to_log <- c("afp", "pivka", "tbil", "pt_inr", "crp",
                 "lr_5_size_ct", "lr_tr_v_size_whole_ct",
                 "lr_tr_v_size_enhancing_ct", "lr_3_size_ct",
                 "lr_4_size_ct", "lr_m_size_ct", "lr_tr_nv_size_ct",
                 "lr_tr_e_size_ct")

factor_vars <- c("sex", "diabetes", "ckd_hd",
                 "hbv", "hcv", "image_lc", "ald", "masld", "ldlt",
                 "pre_tx", "tace", "rfa", "resection", "rtx",
                 "lr_5_ct", "lr_tr_v_ct", "lr_3_ct", "lr_4_ct",
                 "lr_m_ct", "lr_tr_nv_ct", "lr_tr_e_ct")

#optim_vars_for_ext <- c("log_lr_5_size_ct", "log_lr_tr_v_size_whole_ct",
#                        "log_lr_tr_nv_size_ct",
#                        "lr_5_ct", "lr_tr_v_ct", "lr_m_ct", "lr_tr_nv_ct",
#                        "log_afp", "log_pivka",
#                        "cr", "log_pt_inr",
#                        "age_lt", "ckd_hd", "hbv", "ldlt")

optim_vars_for_ext <- c("log_lr_5_size_ct", "log_lr_tr_v_size_whole_ct", "log_lr_tr_nv_size_ct", 
                        "lr_5_ct", "lr_tr_v_ct", "lr_m_ct", "lr_tr_nv_ct", "lr_tr_e_ct",
						"log_afp", "log_pivka",
						"cr", "log_pt_inr", 
                        "age_lt", "ckd_hd", "hbv", "ldlt", "rfa")


## ── RSF CIF interpolation (same as main script) ───────────────────────────────
interp_predict_risk_rsf <- function(pred, times) {
  time_vec <- pred$time.interest
  interp_cif_matrix <- matrix(NA, nrow=nrow(pred$cif), ncol=length(times))
  for (k in seq_len(nrow(pred$cif))) {
    individual_cif <- pred$cif[k, , 1]
    for (i in seq_along(times)) {
      interp_cif_matrix[k, i] <- approx(x=time_vec, y=individual_cif,
                                         xout=times[i], rule=2)$y
    }
  }
  colnames(interp_cif_matrix) <- paste("t", times, sep="_")
  as.data.frame(interp_cif_matrix)
}

## ── Load test data ────────────────────────────────────────────────────────────
test_data_fname <- "../../data/preprocessed_competing_risk_validation_dataset250624.csv"
test_data       <- read.csv(test_data_fname)

vars_to_log_test <- intersect(vars_to_log, names(test_data))
for (var_nm in vars_to_log_test)
  test_data[[paste0("log_", var_nm)]] <- log(test_data[[var_nm]] + 1)
test_data[intersect(factor_vars, names(test_data))] <-
  lapply(test_data[intersect(factor_vars, names(test_data))], factor)

final_vars <- intersect(optim_vars_for_ext, names(test_data))

# complete cases only (same as main script)
complete_cases <- complete.cases(test_data[c("ftime", "fstatus", final_vars)])
test_data_clean <- test_data[complete_cases, ]
cat(sprintf("Test data: %d rows, %d complete\n", nrow(test_data), nrow(test_data_clean)))

## ── RSF formula (for pec::cindex) ────────────────────────────────────────────
formula_str <- paste("Surv(ftime, fstatus) ~", paste(final_vars, collapse=" + "))
rsf_formula <- as.formula(formula_str)

## ── Load fitted model ─────────────────────────────────────────────────────────
fit <- readRDS("fitted_rsf.rds")
cat("Loaded: fitted_rsf.rds\n")

## ── Load shared bootstrap indices ─────────────────────────────────────────────
boot_indices_all <- readRDS("../bootstrap_indices.rds")
N_BOOT           <- length(boot_indices_all)
cat(sprintf("Loaded %d bootstrap indices\n", N_BOOT))

# bootstrap_indices.rds는 test_data 기준으로 생성됨
# complete_cases 필터 적용: clean data 내부 인덱스로 재매핑
clean_idx <- which(complete_cases)   # original row positions that are clean

## ── Helper: compute all metrics on one bootstrap sample ───────────────────────
compute_rsf_metrics <- function(d) {

  ## C-index (pec::cindex)
  cindex <- tryCatch({
    pec::cindex(list(rsf = fit),
                formula = rsf_formula,
                data    = d,
                cause   = 1,
                verbose = FALSE)$AppCindex$rsf
  }, error = function(e) NA_real_)

  ## Predict CIF at eval_times
  pred <- tryCatch(
    predict(fit, newdata = d),
    error = function(e) NULL
  )
  if (is.null(pred)) return(NULL)

  interp_cif <- interp_predict_risk_rsf(pred, eval_times)

  ## AUROC + Brier
  sc <- tryCatch(
    Score(list("RSF" = fit),
          formula    = Hist(ftime, fstatus) ~ 1,
          data       = d,
          times      = eval_times,
          metrics    = c("AUC", "Brier"),
          cause      = 1,
          null.model = FALSE),
    error = function(e) NULL
  )
  auroc <- if (!is.null(sc)) sc$AUC$score[model=="RSF",   AUC]   else rep(NA_real_, length(eval_times))
  bs    <- if (!is.null(sc)) sc$Brier$score[model=="RSF", Brier] else rep(NA_real_, length(eval_times))

  ## IBS via trapz
  ibs <- tryCatch({
    event_times <- sort(unique(subset(d, fstatus %in% c(1,2))$ftime))
    bs_all <- Score(list("RSF" = fit),
                    formula    = Hist(ftime, fstatus) ~ 1,
                    data       = d,
                    times      = event_times,
                    metrics    = "Brier",
                    cause      = 1,
                    null.model = FALSE)
    dt <- bs_all$Brier$score[model == "RSF", .(times, Brier)]
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
cat("Running bootstrap on test data (no model refitting)...\n")
boot_results <- vector("list", N_BOOT)

for (b in seq_len(N_BOOT)) {
  if (b %% 10 == 0) cat(sprintf("  iter %d / %d\n", b, N_BOOT))

  # boot_indices는 원본 test_data 기준 → clean subset 내 인덱스로 변환
  raw_idx    <- boot_indices_all[[b]]
  # clean_idx 내에서만 선택 (complete cases 유지)
  valid_raw  <- raw_idx[raw_idx %in% clean_idx]
  # clean_idx 내 위치로 변환
  local_idx  <- match(valid_raw, clean_idx)

  if (length(local_idx) < 10) {
    cat(sprintf("  iter %d: too few complete cases (%d), skipping\n", b, length(local_idx)))
    next
  }

  d_boot         <- test_data_clean[local_idx, ]
  boot_results[[b]] <- compute_rsf_metrics(d_boot)
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
#saveRDS(list(boot_results = boot_results_clean,
#             summary_df   = summary_df,
#             eval_times   = eval_times,
#             year_labels  = year_labels),
#        "rsf_boot_results.rds")

#cat("\nSaved → rsf_boot_results.rds\n")
