## pairwise_bootstrap_comparison.R
## 5개 모델의 bootstrap results로 pairwise p-value 테이블 생성
## Metrics: C-index (overall), Brier score (1/3/5/10yr), IBS (overall)

library(dplyr)

## ── 1. Load all bootstrap results ─────────────────────────────────────────────
model_names <- c("Milan", "UtS", "Mt2.0", "FG", "RSF")
rds_files   <- c("milan_boot_results.rds", "uts_boot_results.rds",
                 "mt2_boot_results.rds",   "FG/fg_boot_results.rds",
                 "RSF/rsf_boot_results.rds")

boot_list <- lapply(rds_files, readRDS)
names(boot_list) <- model_names

# Extract summary_df per model
summ <- lapply(boot_list, function(x) x$summary_df)

cat("Bootstrap iterations per model:\n")
for (nm in model_names) cat(sprintf("  %-8s %d\n", nm, nrow(summ[[nm]])))

## ── 2. Paired bootstrap p-value (two-sided) ───────────────────────────────────
# For metrics where HIGHER = better (C-index, AUROC):
#   diff = A - B;  p = 2 * min(mean(diff>0), mean(diff<0))  [no direction assumed]
# For metrics where LOWER = better (Brier, IBS):
#   same formula — p-value is always two-sided

paired_pval <- function(vec_a, vec_b) {
  n    <- min(length(vec_a), length(vec_b))
  diff <- vec_a[1:n] - vec_b[1:n]
  p    <- 2 * min(mean(diff > 0, na.rm=TRUE),
                  mean(diff < 0, na.rm=TRUE))
  max(p, 1/n)   # floor at 1/N_BOOT
}

fmt_p <- function(p) {
  if (is.na(p))   return("")
  if (p < 0.001)  return("<0.001")
  sprintf("%.3f", p)
}

## ── 3. Build upper-triangular p-value matrix ──────────────────────────────────
make_pval_table <- function(metric_col, model_list=model_names, summ_list=summ) {
  n  <- length(model_list)
  mat <- matrix("", nrow=n, ncol=n,
                dimnames=list(model_list, model_list))
  for (i in seq_len(n-1)) {
    for (j in (i+1):n) {
      a <- summ_list[[model_list[i]]][[metric_col]]
      b <- summ_list[[model_list[j]]][[metric_col]]
      if (is.null(a) || is.null(b)) next
      mat[i, j] <- fmt_p(paired_pval(a, b))
    }
  }
  mat
}

## ── 4. Define metrics to compare ──────────────────────────────────────────────
year_labels <- c(1, 3, 5, 10)

metrics_brier <- paste0("brier_", year_labels, "yr")
metrics_auroc <- paste0("auroc_", year_labels, "yr")

## Check available columns
avail_cols <- names(summ[[1]])
cat("\nAvailable metric columns:\n"); print(avail_cols)

## ── 5. Print & save all tables ────────────────────────────────────────────────
print_upper_tri <- function(mat, title) {
  cat(sprintf("\n══════════════════════════════════════════\n"))
  cat(sprintf(" %s\n", title))
  cat(sprintf("══════════════════════════════════════════\n"))
  n <- nrow(mat)
  # Header
  cat(sprintf("%-10s", "Model"))
  for (j in 2:n) cat(sprintf("%-12s", colnames(mat)[j]))
  cat("\n")
  # Rows
  for (i in 1:(n-1)) {
    cat(sprintf("%-10s", rownames(mat)[i]))
    for (j in 1:n) {
      if (j <= i) cat(sprintf("%-12s", ""))
      else        cat(sprintf("%-12s", mat[i, j]))
    }
    cat("\n")
  }
}

all_tables <- list()

## C-index
if ("cindex" %in% avail_cols) {
  tab <- make_pval_table("cindex")
  print_upper_tri(tab, "C-index (Overall) — pairwise p-value")
  all_tables[["cindex_overall"]] <- tab
}

## IBS
if ("ibs" %in% avail_cols) {
  tab <- make_pval_table("ibs")
  print_upper_tri(tab, "IBS (Overall) — pairwise p-value")
  all_tables[["ibs_overall"]] <- tab
}

## Brier score by year
for (i in seq_along(year_labels)) {
  col <- metrics_brier[i]
  if (col %in% avail_cols) {
    tab <- make_pval_table(col)
    ttl <- sprintf("Brier Score (%d-year) — pairwise p-value", year_labels[i])
    print_upper_tri(tab, ttl)
    all_tables[[sprintf("brier_%dyr", year_labels[i])]] <- tab
  }
}

## AUROC by year (bonus — for cross-check with DeLong's)
for (i in seq_along(year_labels)) {
  col <- metrics_auroc[i]
  if (col %in% avail_cols) {
    tab <- make_pval_table(col)
    ttl <- sprintf("AUROC (%d-year) — pairwise p-value [bootstrap]", year_labels[i])
    print_upper_tri(tab, ttl)
    all_tables[[sprintf("auroc_%dyr", year_labels[i])]] <- tab
  }
}

## ── 6. Save combined result ───────────────────────────────────────────────────
saveRDS(all_tables, "pairwise_pval_tables.rds")
cat("\nSaved → pairwise_pval_tables.rds\n")

## ── 7. Export to CSV (one file per metric) ────────────────────────────────────
for (nm in names(all_tables)) {
  mat <- all_tables[[nm]]
  df  <- as.data.frame(mat)
  df  <- cbind(Model = rownames(mat), df)
  write.csv(df, file=sprintf("pval_%s.csv", nm), row.names=FALSE)
}
cat("CSV files saved for each metric.\n")

## ── 8. Summary: mean ± SD per model per metric ────────────────────────────────
cat("\n══════════════════════════════════════════\n")
cat(" Mean ± SD summary across bootstrap iterations\n")
cat("══════════════════════════════════════════\n")

key_metrics <- c("cindex", "ibs",
                 paste0("auroc_", year_labels, "yr"),
                 paste0("brier_", year_labels, "yr"))
key_metrics <- intersect(key_metrics, avail_cols)

# Header
cat(sprintf("%-12s", "Metric"))
for (nm in model_names) cat(sprintf("%-22s", nm))
cat("\n")

for (m in key_metrics) {
  cat(sprintf("%-12s", m))
  for (nm in model_names) {
    v <- summ[[nm]][[m]]
    if (is.null(v)) { cat(sprintf("%-22s", "N/A")); next }
    cat(sprintf("%-22s", sprintf("%.3f ± %.3f", mean(v,na.rm=T), sd(v,na.rm=T))))
  }
  cat("\n")
}


## ── 9. Heatmap: C-index & IBS ────────────────────────────────────────────────
library(ggplot2)
library(tidyr)

display_order <- c("Milan", "UtS", "Mt2.0", "FG", "RSF")

# String matrix → symmetric numeric matrix 변환
mat_to_long <- function(mat_str, metric_title) {
  # string p-value → numeric
  str2num <- function(x) {
    x[x == ""]       <- NA
    x[x == "<0.001"] <- "0.001"
    as.numeric(x)
  }
  n   <- nrow(mat_str)
  num <- matrix(str2num(mat_str), nrow=n, dimnames=dimnames(mat_str))
  # 대칭 채우기
  for (i in 1:n) for (j in 1:n) {
    if (is.na(num[i,j]) && !is.na(num[j,i])) num[i,j] <- num[j,i]
  }
  # long format
  df <- as.data.frame(num) |>
    mutate(row = rownames(num)) |>
    pivot_longer(-row, names_to="col", values_to="pvalue") |>
    mutate(
      metric = metric_title,
      row    = factor(row, levels=display_order),
      col    = factor(col, levels=display_order),
      sig    = case_when(
        is.na(pvalue)   ~ "",
        pvalue < 0.01   ~ paste0(sprintf("%.2f", pvalue), "**"),
        pvalue < 0.05   ~ paste0(sprintf("%.2f", pvalue), "*"),
        TRUE            ~ sprintf("%.2f", pvalue)
      )
    )
  df
}

draw_heatmap <- function(df_long, filename, fig_width=6, fig_height=6, tagname="(B)") {
  p <- ggplot(df_long, aes(col, row, fill=pvalue)) +
    geom_tile(color="white", linewidth=0.5) +
    geom_text(aes(label=sig), size=4) +
    scale_fill_distiller(palette="Blues", direction=1, limits=c(0,1),
                         na.value="white", name="p-value",
                         breaks=c(0, 0.25, 0.50, 0.75, 1.00)) +
    scale_y_discrete(limits=levels(df_long$row)) +
    coord_fixed() +
    labs(tag = tagname) +
    theme_minimal(base_size=12) +
    theme(
      axis.title        = element_blank(),
      panel.grid        = element_blank(),
      axis.text.x       = element_text(angle=45, hjust=1, color="black", size=12),
      axis.text.y       = element_text(color="black", size=12),
      axis.ticks        = element_line(color="black", linewidth=0.4),
      axis.ticks.length = unit(4, "pt"),
      axis.line         = element_line(color="black", linewidth=0.4),
      plot.tag          = element_text(size=12, face="bold"),
      plot.tag.position = c(0.015, 0.97), # left upper: (x = 0, y = 1) 
      legend.title      = element_text(size=12),
      legend.text       = element_text(size=11)
    )
  ggsave(filename, plot=p, width=fig_width, height=fig_height)
  message("Saved → ", filename)
}

# C-index heatmap
if (!is.null(all_tables[["cindex_overall"]])) {
  df_cindex <- mat_to_long(all_tables[["cindex_overall"]], "C-index")
  draw_heatmap(df_cindex, "table_3_heatmap-Cindex.pdf", 
			   fig_width=506/72, fig_height=401/72, tagname="(B)")
}

# IBS heatmap
if (!is.null(all_tables[["ibs_overall"]])) {
  df_ibs <- mat_to_long(all_tables[["ibs_overall"]], "IBS")
  draw_heatmap(df_ibs, "table_3_heatmap-ibs.pdf", 
			   fig_width=506/72, fig_height=401/72, tagname="(C)")
}
