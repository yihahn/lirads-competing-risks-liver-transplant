library(riskRegression)
library(ggplot2)
library(dplyr)
library(tidyr)

## ── 1. Load & preprocess ──────────────────────────────────────────────────────
train_data_fname <- "../data/internal_competing_risk_dataset.csv"
train_data       <- read.csv(train_data_fname)

vars_to_log <- c("lr_5_size_ct", "lr_tr_v_size_enhancing_ct",
                 "lr_m_size_ct", "lr_tr_nv_size_ct",
                 "afp", "pivka", "pt_inr")
factor_vars <- c("sex", "ckd_hd")

for (var_nm in vars_to_log)
  train_data[[paste0("log_", var_nm)]] <- log(train_data[[var_nm]] + 1)
train_data[factor_vars] <- lapply(train_data[factor_vars], factor)

## ── 2. Load model ─────────────────────────────────────────────────────────────
fit   <- readRDS("FG/fitted_fg.rds")
coefs <- fit$crrFit$coef

## ── 3. Variable meta ──────────────────────────────────────────────────────────
model_vars <- names(coefs)

log_var_map <- c(
  log_lr_5_size_ct              = "lr_5_size_ct",
  log_lr_tr_v_size_enhancing_ct = "lr_tr_v_size_enhancing_ct",
  log_lr_m_size_ct              = "lr_m_size_ct",
  log_lr_tr_nv_size_ct          = "lr_tr_nv_size_ct",
  log_afp                       = "afp",
  log_pivka                     = "pivka",
  log_pt_inr                    = "pt_inr"
)

sci_notation_vars <- c("log_afp", "log_pivka")
size_vars         <- c("log_lr_5_size_ct", "log_lr_tr_v_size_enhancing_ct",
                       "log_lr_m_size_ct",  "log_lr_tr_nv_size_ct")

# !! Manually specified max mm for each size variable
size_max_mm <- c(
  log_lr_5_size_ct              = 80,
  log_lr_tr_v_size_enhancing_ct = 80,
  log_lr_m_size_ct              = 40,
  log_lr_tr_nv_size_ct          = 100
)

factor_coef_map <- list(
  sex2    = c("Male", "Female"),
  ckd_hd1 = c("No",  "Yes")
)

var_labels <- c(
  log_lr_5_size_ct              = "Diameter of LR-5 (mm)",
  log_lr_tr_v_size_enhancing_ct = "Diameter of LR-TR-V (mm)",
  log_lr_m_size_ct              = "Diameter of LR-M (mm)",
  log_lr_tr_nv_size_ct          = "Diameter of LR-TR-NV (mm)",
  log_afp                       = "AFP (ng/mL)",
  log_pivka                     = "PIVKA-II (mAU/mL)",
  log_pt_inr                    = "PT (INR)",
  sex2                          = "Sex",
  ckd_hd1                       = "CKD or on hemodialysis"
)

## ── 4. Ranges & ticks ─────────────────────────────────────────────────────────
var_ranges      <- list()
var_ticks       <- list()
var_tick_labels <- list()
var_tick_parse  <- list()
var_minor_ticks <- list()

for (v in model_vars) {
  if (v %in% names(log_var_map)) {
    orig_x <- train_data[[log_var_map[v]]]
    orig_x <- orig_x[!is.na(orig_x)]

    if (v %in% size_vars) {
      # mm range: always start at 0, end at specified max_mm
      min_mm  <- 0
      max_mm  <- size_max_mm[v]
      # log scale range in cm: log(mm/10 + 1)
      var_ranges[[v]] <- c(log(min_mm/10 + 1), log(max_mm/10 + 1))

      # Major ticks every 20mm, from 0 to max_mm
      major_mm <- seq(0, max_mm, by = 20)
      var_ticks[[v]]       <- log(major_mm/10 + 1)
      var_tick_labels[[v]] <- as.character(major_mm)
      var_tick_parse[[v]]  <- FALSE

      # Minor ticks every 5mm, excluding major
      minor_mm <- seq(0, max_mm, by = 5)
      minor_mm <- minor_mm[!minor_mm %in% major_mm]
      var_minor_ticks[[v]] <- log(minor_mm/10 + 1)

    } else if (v %in% sci_notation_vars) {
      # AFP & PIVKA: log10 ticks
      var_ranges[[v]] <- range(log(orig_x + 1))
      log10_range  <- log10(range(orig_x + 1))
      log10_ticks  <- seq(floor(log10_range[1]), ceiling(log10_range[2]), by = 1)
      log10_ticks  <- log10_ticks[log10_ticks >= log10_range[1] &
                                    log10_ticks <= log10_range[2]]
      orig_tick_vals       <- 10^log10_ticks
      var_ticks[[v]]       <- log(orig_tick_vals + 1)
      var_tick_labels[[v]] <- ifelse(log10_ticks == 0, "1",
                                     paste0("10^", log10_ticks))
      var_tick_parse[[v]]  <- TRUE
      var_minor_ticks[[v]] <- numeric(0)

    } else {
      # PT/INR: plain
      var_ranges[[v]] <- range(log(orig_x + 1))
      orig_ticks <- pretty(orig_x, n = 5)
      orig_ticks <- orig_ticks[orig_ticks >= 0 &
                                 orig_ticks >= min(orig_x) &
                                 orig_ticks <= max(orig_x)]
      var_ticks[[v]]       <- log(orig_ticks + 1)
      var_tick_labels[[v]] <- as.character(orig_ticks)
      var_tick_parse[[v]]  <- FALSE
      var_minor_ticks[[v]] <- numeric(0)
    }

  } else if (v %in% names(factor_coef_map)) {
    var_ranges[[v]]      <- c(0, 1)
    var_ticks[[v]]       <- c(0, 1)
    var_tick_labels[[v]] <- factor_coef_map[[v]]
    var_tick_parse[[v]]  <- FALSE
    var_minor_ticks[[v]] <- numeric(0)
  }
}

# Diagnostics
cat("Variable ranges (log scale):\n")
for (v in model_vars)
  cat(sprintf("  %-42s [%.4f, %.4f]\n", v, var_ranges[[v]][1], var_ranges[[v]][2]))

## ── 5. Points scaling ─────────────────────────────────────────────────────────
MAX_PTS <- 25

lp_min_v <- sapply(model_vars, function(v) min(coefs[v] * var_ranges[[v]]))
lp_max_v <- sapply(model_vars, function(v) max(coefs[v] * var_ranges[[v]]))

max_single_range <- max(lp_max_v - lp_min_v)
scale_f          <- MAX_PTS / max_single_range

total_lp_min  <- sum(lp_min_v)
total_lp_max  <- sum(lp_max_v)
total_pts_max <- (total_lp_max - total_lp_min) * scale_f
cat(sprintf("\nTotal Points max: %.1f\n", total_pts_max))

val_to_pts      <- function(v, val) (coefs[v] * val - lp_min_v[v]) * scale_f
tp_to_phys      <- function(tp) tp / total_pts_max * MAX_PTS
lp_to_total_pts <- function(lp) (lp - total_lp_min) * scale_f
cif_to_lp       <- function(cif, h0) log(-log(1 - cif) / h0)

## ── 6. Baseline CIF ───────────────────────────────────────────────────────────
uftime  <- fit$crrFit$uftime
cum_haz <- cumsum(fit$crrFit$bfitj)

time_points <- c("1-year"=365, "3-year"=1095, "5-year"=1825, "10-year"=3650)
get_H0 <- function(t) {
  idx <- which(uftime <= t)
  if (!length(idx)) { warning(sprintf("t=%g before first event", t)); return(NA) }
  cum_haz[max(idx)]
}
H0 <- sapply(time_points, get_H0)

## ── 7. CIF axes ───────────────────────────────────────────────────────────────
cif_vals <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4,
              0.5,  0.6,  0.7, 0.8, 0.9, 0.95, 0.99)

cif_axis_list <- lapply(names(time_points), function(nm) {
  h0 <- H0[nm]
  if (is.na(h0)) return(NULL)
  result <- lapply(cif_vals, function(cif_val) {
    lp <- tryCatch(cif_to_lp(cif_val, h0), error = function(e) NA)
    if (is.na(lp) || is.infinite(lp)) return(NULL)
    tp <- lp_to_total_pts(lp)
    px <- tp_to_phys(tp)
    if (is.nan(px) || is.infinite(px) || px < -0.1 || px > MAX_PTS + 0.1) return(NULL)
    data.frame(cif=cif_val, px=px, axis_label=nm,
               row.names=NULL, stringsAsFactors=FALSE)
  })
  bind_rows(result)
}) %>% bind_rows()

cif_span <- cif_axis_list %>%
  group_by(axis_label) %>%
  summarise(x_min=min(px), x_max=max(px), .groups="drop")

## ── 8. Build tick_df ──────────────────────────────────────────────────────────
n_vars <- length(model_vars)
y_map  <- setNames(seq(n_vars, 1), model_vars)

tick_list <- lapply(model_vars, function(v) {
  data.frame(
    var      = v,
    label    = unname(var_labels[v]),
    y        = unname(y_map[v]),
    pts      = unname(val_to_pts(v, var_ticks[[v]])),
    tick_lab = unname(var_tick_labels[[v]]),
    do_parse = var_tick_parse[[v]],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
})
tick_df <- bind_rows(tick_list)

minor_list <- lapply(model_vars, function(v) {
  mt <- var_minor_ticks[[v]]
  if (length(mt) == 0) return(NULL)
  data.frame(
    var = v,
    y   = unname(y_map[v]),
    pts = unname(val_to_pts(v, mt)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
})
minor_df <- bind_rows(minor_list)

seg_df <- tick_df %>%
  group_by(var, y, label) %>%
  summarise(xmin=min(pts), xmax=max(pts), .groups="drop")

## ── 9. Total Points axis ──────────────────────────────────────────────────────
tp_major_breaks <- pretty(c(0, total_pts_max), n=6)
tp_major_breaks <- tp_major_breaks[tp_major_breaks >= 0 &
                                     tp_major_breaks <= total_pts_max]
tp_major_phys   <- tp_to_phys(tp_major_breaks)

tp_minor_real <- unlist(lapply(seq_len(length(tp_major_breaks)-1), function(i) {
  seq(tp_major_breaks[i], tp_major_breaks[i+1], length.out=6)[2:5]
}))
tp_minor_phys <- tp_to_phys(tp_minor_real)

## ── 10. Y positions ───────────────────────────────────────────────────────────
y_pts   <- n_vars + 1.5
y_total <- 0
cif_gap <- -1.8
y_cif   <- setNames(cif_gap * seq_along(time_points), names(time_points))

tick_h       <- 0.18
tick_h_minor <- 0.10
txt_off      <- 0.36
x_lab        <- -0.6
pts_breaks   <- seq(0, MAX_PTS, by=5)
pts_minor    <- setdiff(seq(0, MAX_PTS, by=1), pts_breaks)

cif_axis_label_str <- paste0(names(time_points), " HCC-specific mortality (CIF)")

## ── 11. Plot ──────────────────────────────────────────────────────────────────
p <- ggplot() +

  ## Variable lines
  geom_segment(data=seg_df,
               aes(x=xmin,xend=xmax,y=y,yend=y), linewidth=0.5) +
  geom_segment(data=tick_df,
               aes(x=pts,xend=pts,y=y-tick_h,yend=y+tick_h), linewidth=0.4) +

  ## Minor ticks (size vars)
  {if (!is.null(minor_df) && nrow(minor_df) > 0)
    geom_segment(data=minor_df,
                 aes(x=pts,xend=pts,y=y-tick_h_minor,yend=y+tick_h_minor),
                 linewidth=0.3)} +

  ## Tick labels
  geom_text(data=filter(tick_df, !do_parse),
            aes(x=pts,y=y-txt_off,label=tick_lab), size=2.8, vjust=1) +
  geom_text(data=filter(tick_df, do_parse),
            aes(x=pts,y=y-txt_off,label=tick_lab), size=2.8, vjust=1, parse=TRUE) +
  geom_text(data=seg_df,
            aes(x=x_lab,y=y,label=label), hjust=1, size=3.3) +

  ## Points axis
  geom_segment(aes(x=0,xend=MAX_PTS,y=y_pts,yend=y_pts), linewidth=0.5) +
  geom_segment(data=data.frame(x=pts_breaks),
               aes(x=x,xend=x,y=y_pts-tick_h,yend=y_pts+tick_h), linewidth=0.4) +
  geom_segment(data=data.frame(x=pts_minor),
               aes(x=x,xend=x,y=y_pts-tick_h_minor,yend=y_pts+tick_h_minor),
               linewidth=0.3) +
  geom_text(data=data.frame(x=pts_breaks),
            aes(x=x,y=y_pts+txt_off,label=x), size=2.8, vjust=0) +
  geom_text(aes(x=x_lab,y=y_pts,label="Points"),
            hjust=1, size=3.5, fontface="bold") +

  ## Total Points axis
  geom_segment(aes(x=0,xend=MAX_PTS,y=y_total,yend=y_total), linewidth=0.5) +
  geom_segment(data=data.frame(x=tp_major_phys),
               aes(x=x,xend=x,y=y_total-tick_h,yend=y_total+tick_h), linewidth=0.4) +
  geom_segment(data=data.frame(x=tp_minor_phys),
               aes(x=x,xend=x,y=y_total-tick_h_minor,yend=y_total+tick_h_minor),
               linewidth=0.3) +
  geom_text(data=data.frame(x=tp_major_phys,lab=round(tp_major_breaks,0)),
            aes(x=x,y=y_total-txt_off,label=lab), size=2.8, vjust=1) +
  geom_text(aes(x=x_lab,y=y_total,label="Total Points"),
            hjust=1, size=3.5, fontface="bold") +

  ## CIF axes
  geom_segment(
    data=cif_span %>% mutate(yy=as.numeric(y_cif[axis_label])),
    aes(x=x_min,xend=x_max,y=yy,yend=yy), linewidth=0.5) +
  geom_segment(
    data=cif_axis_list %>% mutate(yy=as.numeric(y_cif[axis_label])),
    aes(x=px,xend=px,y=yy-tick_h,yend=yy+tick_h), linewidth=0.4) +
  geom_text(
    data=cif_axis_list %>% mutate(yy=as.numeric(y_cif[axis_label])),
    aes(x=px,y=yy-txt_off,label=as.character(cif)),
    size=2.3, vjust=1, na.rm=TRUE) +
  geom_text(
    data=data.frame(yy=as.numeric(y_cif), lbl=cif_axis_label_str),
    aes(x=x_lab,y=yy,label=lbl), hjust=1, size=3.2) +

  theme_void() +
  theme(
    plot.background  = element_rect(fill="white", color=NA),
    panel.background = element_rect(fill="white", color=NA),
    plot.margin      = margin(15, 20, 20, 160, unit="pt")
  ) +
  coord_cartesian(
    xlim = c(-3.5, MAX_PTS + 0.5),
    ylim = c(min(as.numeric(y_cif)) - 1, y_pts + 1),
    clip = "off"
  )

## ── 12. Save ──────────────────────────────────────────────────────────────────
ggsave("nomogram_FG.pdf",  p, width=18, height=13, bg="white")
ggsave("nomogram_FG.tiff", p, width=18, height=13, dpi=300,
       compression="lzw", device="tiff", bg="white")
cat("Saved: nomogram_FG.pdf / nomogram_FG.tiff\n")
