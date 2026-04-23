## Step 1: generate_bootstrap_indices.R
## 모든 모델이 공유하는 bootstrap 인덱스 생성 (1회만 실행)

set.seed(20250521)
N_BOOT <- 100

test_data_fname <- "../data/preprocessed_competing_risk_validation_dataset250624.csv"
test_data       <- read.csv(test_data_fname)
n               <- nrow(test_data)

boot_indices <- lapply(1:N_BOOT, function(i) sample(n, size=n, replace=TRUE))

saveRDS(boot_indices, "bootstrap_indices.rds")
cat(sprintf("Saved %d bootstrap indices (n=%d each) → bootstrap_indices.rds\n",
            N_BOOT, n))
