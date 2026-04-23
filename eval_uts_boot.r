## Step 2: eval_uts_boot.R
## bootstrap_indices.rds를 읽어 UTS.0 모델의 bootstrap 평가 수행

library(dplyr)
library(survival)
library(fastcmprsk)
library(riskRegression)
library(timeROC)
library(pracma)
source("wolbers_cindex_cr_claude.r")

## ── Config ────────────────────────────────────────────────────────────────────
eval_times <- c(365.24, 365.24*3, 365.24*5, 365.24*10)
year_labels <- round(eval_times / 365.24)

test_data_fname <- "../data/preprocessed_competing_risk_validation_dataset250624.csv"
test_data       <- read.csv(test_data_fname)

data_cr <- data.frame(
  ftime   = test_data$ftime,
  fstatus = test_data$fstatus,
  uts     = as.numeric(test_data$uptoseven)
)

## ── Load shared bootstrap indices ─────────────────────────────────────────────
boot_indices <- readRDS("bootstrap_indices.rds")
N_BOOT       <- length(boot_indices)
cat(sprintf("Loaded %d bootstrap indices\n", N_BOOT))

## ── Helper: compute all metrics on a dataset ──────────────────────────────────
compute_metrics <- function(d) {

  ## C-index (Wolbers)
  cindex <- tryCatch({
    wolbers_c_index_competing(
      times              = d$ftime,
      events             = d$fstatus,
      predictions        = d$uts,
      event_of_interest  = 1
    )$c_index
  }, error = function(e) NA_real_)

  ## AUROC + Brier at eval_times
  sc <- tryCatch(
    Score(list("UTS" = d$uts),
          formula = Hist(ftime, fstatus) ~ 1,
          data    = d,
          times   = eval_times,
          metrics = c("AUC", "Brier"),
          cause   = 1),
    error = function(e) NULL
  )

  auroc <- if (!is.null(sc))
    sc$AUC$score[model == "UTS", AUC]   else rep(NA_real_, length(eval_times))
  bs    <- if (!is.null(sc))
    sc$Brier$score[model == "UTS", Brier] else rep(NA_real_, length(eval_times))

  ## IBS via trapz over all event times
  event_times <- subset(d, fstatus %in% c(1, 2))$ftime

  ibs <- tryCatch({
    bs_all <- Score(list("UTS" = d$uts),
                    formula = Hist(ftime, fstatus) ~ 1,
                    data    = d,
                    times   = event_times,
                    metrics = "Brier",
                    cause   = 1)
    dt  <- bs_all$Brier$score[model == "UTS", .(times, Brier)]
    trapz(dt$times, dt$Brier) / max(event_times)
  }, error = function(e) NA_real_)

  ## Return as named list
  result <- list(cindex = cindex, ibs = ibs)
  for (i in seq_along(eval_times)) {
    result[[paste0("auroc_", year_labels[i], "yr")]]  <- auroc[i]
    result[[paste0("brier_", year_labels[i], "yr")]]  <- bs[i]
  }
  result
}

## ── Run bootstrap ─────────────────────────────────────────────────────────────
cat("Running bootstrap...\n")
boot_results <- vector("list", N_BOOT)

for (b in seq_len(N_BOOT)) {
  if (b %% 10 == 0) cat(sprintf("  iter %d / %d\n", b, N_BOOT))
  d_boot         <- data_cr[boot_indices[[b]], ]
  boot_results[[b]] <- compute_metrics(d_boot)
}

## ── Summarise ─────────────────────────────────────────────────────────────────
metrics <- names(boot_results[[1]])
summary_df <- do.call(rbind, lapply(boot_results, function(r) {
  as.data.frame(r, stringsAsFactors = FALSE)
}))

cat("\n--- Bootstrap summary (mean ± SD) ---\n")
for (m in metrics) {
  vals <- summary_df[[m]]
  cat(sprintf("  %-25s %.4f ± %.4f\n", m, mean(vals, na.rm=TRUE), sd(vals, na.rm=TRUE)))
}

## ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(list(boot_results = boot_results,
             summary_df   = summary_df,
             eval_times   = eval_times,
             year_labels  = year_labels),
        "uts_boot_results.rds")

cat("\nSaved → uts_boot_results.rds\n")
