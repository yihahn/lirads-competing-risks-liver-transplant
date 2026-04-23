library(ggplot2)
library(dplyr)
library(patchwork)

fg_results  <- readRDS("FG/cpi_fg_wolbers_results.rds")
rsf_results <- readRDS("RSF/cpi_rsf_wolbers_results.rds")

prep_data <- function(results) {
  df <- data.frame(
    variable = names(results$drop_cindex),
    vimp     = as.numeric(results$drop_cindex),
    stringsAsFactors = FALSE
  )
  df <- df[order(df$vimp, decreasing = FALSE), ]
  df$variable <- factor(df$variable, levels = df$variable)
  df
}

fg_data  <- prep_data(fg_results)
rsf_data <- prep_data(rsf_results)

base_theme <- theme_minimal(base_size = 9) +
  theme(
    plot.title        = element_text(hjust = 0.5, size = 11, face = "bold"),
    axis.text.y       = element_text(size = 7),
    axis.text.x       = element_text(size = 7),
    axis.title.x      = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(10, 15, 10, 5)
  )

make_plot <- function(df, title) {
  ggplot(df, aes(x = vimp, y = variable,
                 fill = ifelse(vimp >= 0, "#7FB069", "#E63946"))) +
    geom_col(alpha = 0.85, width = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_identity() +
    scale_x_continuous(
      labels = function(x) sprintf("%.3f", x),
      expand = expansion(mult = c(0.05, 0.1))
    ) +
    labs(title = title, x = "Drop in Wolbers C-index", y = "") +
    base_theme
}

p1 <- make_plot(fg_data,  "FG model")
p2 <- make_plot(rsf_data, "RSF model")

# RSF 변수 수가 많으므로 width 1 : 1.5 배분
combined_plot <- p1 + p2 + plot_layout(widths = c(1, 1.5))

ggsave("unified_cpi_vimp_cindex.pdf",  combined_plot,
       width = 10, height = 7.5, device = "pdf")
ggsave("unified_cpi_vimp_cindex.tiff", combined_plot,
       width = 10, height = 7.5, dpi = 300, device = "tiff", compression = "lzw")
