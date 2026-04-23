library(prodlim)
library(riskRegression)

# ── 1. 원본 모델 및 데이터 로드 ──────────────────────
fitted_fg <- readRDS("fitted_fg.rds")

# Train_data
data_fname = "../../data/internal_competing_risk_dataset.csv"
data <- read.csv(data_fname)
print(colnames(data))

data$log_lr_5_size_ct              <- log(data$lr_5_size_ct + 1)
data$log_lr_tr_v_size_enhancing_ct <- log(data$lr_tr_v_size_enhancing_ct + 1)
data$log_lr_m_size_ct              <- log(data$lr_m_size_ct + 1)
data$log_lr_tr_nv_size_ct          <- log(data$lr_tr_nv_size_ct + 1)
data$log_afp                       <- log(data$afp + 1)
data$log_pivka                     <- log(data$pivka + 1)
data$log_pt_inr                    <- log(data$pt_inr + 1)

# factor 변환 (sex2 = reference가 아닌 쪽)
data$sex <- factor(data$sex)
data$ckd_hd <- factor(data$ckd_hd, levels = c(0, 1))
levels(factor(data$sex))

# 원본 sHR (point estimate)
original_coef <- fitted_fg$crrFit$coef
sHR_point     <- exp(original_coef)

# ── formula: terms에서 직접 재구성 ───────────────────
var_names  <- attr(fitted_fg$terms, "term.labels")
fg_formula <- as.formula(
  paste("Hist(ftime, fstatus) ~", paste(var_names, collapse = " + "))
)
print(fg_formula)

# ── 3. Bootstrap (n=100) ─────────────────────────────
set.seed(42)
n_boot    <- 100
boot_coef <- matrix(NA, nrow = n_boot, ncol = length(original_coef),
                    dimnames = list(NULL, names(original_coef)))

for (i in 1:n_boot) {
  cat(sprintf("Bootstrap %d/%d\n", i, n_boot))

  tryCatch({
    idx       <- sample(nrow(data), replace = TRUE)
    boot_data <- data[idx, ]

    fit <- FGR(fg_formula, data = boot_data, cause = 1)

    boot_coef[i, ] <- fit$crrFit$coef  # coef() 사용 안 함

  }, error = function(e) {
    cat(sprintf("Bootstrap %d failed: %s", i, e$message))
  })
}


# ── 4. 95% CI (Percentile) ───────────────────────────
CI_lower <- exp(apply(boot_coef, 2, quantile, 0.025, na.rm = TRUE))
CI_upper <- exp(apply(boot_coef, 2, quantile, 0.975, na.rm = TRUE))

# ── 5. 결과 정리 ─────────────────────────────────────
result <- data.frame(
  Variable = names(sHR_point),
  sHR      = round(sHR_point, 3),
  CI_lower = round(CI_lower, 3),
  CI_upper = round(CI_upper, 3)
)
result$CI_95 <- sprintf("%.3f (%.3f–%.3f)",
                         result$sHR,
                         result$CI_lower,
                         result$CI_upper)

print(result)
write.csv(result, "sHR_bootstrap.csv", row.names = FALSE)
