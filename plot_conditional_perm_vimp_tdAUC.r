library(pheatmap)
library(dplyr)
library(RColorBrewer)

# 데이터 로드
cpi_fg_results <- readRDS("FG/cpi_fg_results.rds")
cpi_rsf_results <- readRDS("RSF/cpi_rsf_results.rds")

# 결과 해석 함수 (min_effect_size 조정)
interpret_importance <- function(cpi_results, min_effect_size = 0.0001) {
    # 1. 통계적으로 유의한 positive importance만 추출
    get_significant_positive <- function(drops, significant, threshold = min_effect_size) {
        result <- drops
        result[!significant | drops < threshold] <- NA
        return(result)
    }

    significant_auc <- get_significant_positive(cpi_results$drop_auc,
                                               cpi_results$auc_significant,
                                               min_effect_size)
    significant_brier <- get_significant_positive(cpi_results$drop_brier,
                                                 cpi_results$brier_significant,
                                                 min_effect_size)
    # 2. 변수별 종합 중요도 계산 (시간대별 평균)
    auc_overall <- rowMeans(significant_auc, na.rm = TRUE)
    brier_overall <- rowMeans(significant_brier, na.rm = TRUE)

    # 3. 변수 분류
    classify_variables <- function(auc_imp, brier_imp) {
        classification <- character(length(auc_imp))
        names(classification) <- names(auc_imp)

        for (i in seq_along(auc_imp)) {
            auc_val <- auc_imp[i]
            brier_val <- brier_imp[i]

            if (is.na(auc_val) && is.na(brier_val)) {
                classification[i] <- "Not Important"
            } else if (!is.na(auc_val) && !is.na(brier_val)) {
                classification[i] <- "Consistently Important"
            } else if (!is.na(auc_val)) {
                classification[i] <- "AUC Important Only"
            } else {
                classification[i] <- "Brier Important Only"
            }
        }
        return(classification)
    }

    var_classification <- classify_variables(auc_overall, brier_overall)

    # 4. 결과 정리
    result_df <- data.frame(
        Variable = names(auc_overall),
        AUC_Importance = round(auc_overall, 6),
        Brier_Importance = round(brier_overall, 6),
        Classification = var_classification,
        stringsAsFactors = FALSE
    )

    # NaN을 NA로 변경
    result_df$AUC_Importance[is.nan(result_df$AUC_Importance)] <- NA
    result_df$Brier_Importance[is.nan(result_df$Brier_Importance)] <- NA

    # 중요도 순으로 정렬 (AUC + Brier 평균, NA는 뒤로)
    combined_importance <- rowMeans(cbind(auc_overall, brier_overall), na.rm = TRUE)
    combined_importance[is.nan(combined_importance)] <- -Inf
    result_df <- result_df[order(combined_importance, decreasing = TRUE), ]

    return(list(
        summary = result_df,
        significant_auc = significant_auc,
        significant_brier = significant_brier,
        raw_results = cpi_results,
        auc_overall = auc_overall
    ))
}


# 헬퍼 함수: 양쪽 형식으로 저장하는 함수
save_heatmap_both_formats <- function(unified_data, model_annotation, annotation_colors,
                                    all_display_labels, fg_vars_length,
                                    base_filename, width, height,
                                    title = "Conditional Permutation Variable Importance",
                                    fontsize_number = 10, cellwidth = 60, cellheight = 20,
                                    fontsize = 14, fontsize_row = 11, fontsize_col = 14) {

  # PDF 저장
  pdf_filename <- paste0(tools::file_path_sans_ext(base_filename), ".pdf")
  pdf(pdf_filename, width = width, height = height)
  par(mar = c(1, 1, 2, 1))

  pheatmap(unified_data,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
#           main = title,
           color = colorRampPalette(c("white", "#E6FFE6", "#98FB98", "#66CDAA"))(50), 
#           color = colorRampPalette(c("white", "darkred"))(50),
           display_numbers = TRUE,
           number_format = "%.4f",
           fontsize_number = fontsize_number,
           cellwidth = cellwidth,
           cellheight = cellheight,
           fontsize = fontsize,
           fontsize_row = fontsize_row,
           fontsize_col = fontsize_col,
           annotation_row = model_annotation,
           annotation_colors = annotation_colors,
           annotation_names_row = TRUE,
           labels_row = all_display_labels,
           gaps_row = c(fg_vars_length),
           border_color = "grey60")

  dev.off()

  # TIFF 저장
  tiff_filename <- paste0(tools::file_path_sans_ext(base_filename), ".tiff")
  tiff(tiff_filename, width = width, height = height, units = "in", res = 300, compression = "lzw")
  par(mar = c(1, 1, 2, 1))

  pheatmap(unified_data,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
#           main = title,
           color = colorRampPalette(c("white", "#E6FFE6", "#98FB98", "#66CDAA"))(50),
#           color = colorRampPalette(c("white", "darkred"))(50),
           display_numbers = TRUE,
           number_format = "%.4f",
           fontsize_number = fontsize_number,
           cellwidth = cellwidth,
           cellheight = cellheight,
           fontsize = fontsize,
           fontsize_row = fontsize_row,
           fontsize_col = fontsize_col,
           annotation_row = model_annotation,
           annotation_colors = annotation_colors,
           annotation_names_row = TRUE,
           labels_row = all_display_labels,
           gaps_row = c(fg_vars_length),
           border_color = "grey60")

  dev.off()

  cat("Heatmap saved as PDF:", pdf_filename, "\n")
  cat("Heatmap saved as TIFF:", tiff_filename, "\n")
}


# 기존 함수들을 TIFF와 PDF 모두 저장하도록 수정
plot_unified_heatmap <- function(fg_results, rsf_results,
                               filename = "unified_feature_importance_auc",  # 확장자 제거
                               width = 12, height = 16) {

    # 기존 데이터 준비 코드 동일
    fg_interpretation <- interpret_importance(fg_results, min_effect_size = 0.0001)
    rsf_interpretation <- interpret_importance(rsf_results, min_effect_size = 0.0001)

    fg_auc_data <- fg_interpretation$significant_auc
    rsf_auc_data <- rsf_interpretation$significant_auc

    all_fg_vars <- rownames(fg_auc_data)
    all_rsf_vars <- rownames(rsf_auc_data)
    
    common_vars <- intersect(all_fg_vars, all_rsf_vars)
    fg_only_vars <- setdiff(all_fg_vars, all_rsf_vars)
    rsf_only_vars <- setdiff(all_rsf_vars, all_fg_vars)
    
    cat("Common variables:", length(common_vars), "\n")
    cat("Common variables:", paste(common_vars, collapse=", "), "\n")
    cat("FG-only variables:", length(fg_only_vars), "\n")
    cat("RSF-only variables:", length(rsf_only_vars), "\n")

    fg_all_vars <- all_fg_vars
    fg_all_importance <- rowMeans(fg_auc_data[fg_all_vars, ], na.rm = TRUE)
    fg_vars_ordered <- names(sort(fg_all_importance, decreasing = TRUE))

    rsf_all_vars <- all_rsf_vars
    rsf_all_importance <- rowMeans(rsf_auc_data[rsf_all_vars, ], na.rm = TRUE)
    rsf_vars_ordered <- names(sort(rsf_all_importance, decreasing = TRUE))

    time_points <- colnames(fg_auc_data)
    
    fg_unique_names <- paste0("FG_", seq_along(fg_vars_ordered), "_", make.names(fg_vars_ordered))
    rsf_unique_names <- paste0("RSF_", seq_along(rsf_vars_ordered), "_", make.names(rsf_vars_ordered))
    
    fg_display_labels <- fg_vars_ordered
    rsf_display_labels <- rsf_vars_ordered

    all_unique_names <- c(fg_unique_names, rsf_unique_names)
    all_display_labels <- c(fg_display_labels, rsf_display_labels)
    
    unified_data <- matrix(0, nrow = length(all_unique_names), ncol = length(time_points))
    rownames(unified_data) <- all_unique_names
    colnames(unified_data) <- time_points
    
    # 데이터 채우기
    fg_auc_data[is.na(fg_auc_data)] <- 0
    unified_data[fg_unique_names, ] <- fg_auc_data[fg_vars_ordered, ]
    
    rsf_auc_data[is.na(rsf_auc_data)] <- 0
    unified_data[rsf_unique_names, ] <- rsf_auc_data[rsf_vars_ordered, ]

    model_annotation <- data.frame(
        Model = c(rep("FG", length(fg_vars_ordered)),
                 rep("RSF", length(rsf_vars_ordered))),
        row.names = all_unique_names
    )

    annotation_colors <- list(
        Model = c("FG" = "#E31A1C", "RSF" = "#1F78B4")
    )

    # PDF와 TIFF 모두 저장
    save_heatmap_both_formats(
      unified_data = unified_data,
      model_annotation = model_annotation,
      annotation_colors = annotation_colors,
      all_display_labels = all_display_labels,
      fg_vars_length = length(fg_vars_ordered),
      base_filename = filename,
      width = width,
      height = height,
#      title = "Conditional Permutation Variable Importance\n(Ordered by Individual Model Importance)",
      fontsize_number = 9,
      cellwidth = 50,
      cellheight = 15,
      fontsize = 11,
      fontsize_row = 10,
      fontsize_col = 12
    )

    # 결과 요약 출력
    cat("\n=== Variable Distribution ===\n")
    cat("Total variables in unified heatmap:", nrow(unified_data), "\n")
    cat("- FG variables (ordered by FG importance):", length(fg_vars_ordered), "\n")
    cat("- RSF variables (ordered by RSF importance):", length(rsf_vars_ordered), "\n")
    cat("\nTop 5 FG variables by importance:\n")
    print(head(fg_vars_ordered, 5))
    cat("\nTop 5 RSF variables by importance:\n")
    print(head(rsf_vars_ordered, 5))

    return(list(
        unified_data = unified_data,
        model_annotation = model_annotation,
        fg_interpretation = fg_interpretation,
        rsf_interpretation = rsf_interpretation,
        fg_vars_ordered = fg_vars_ordered,
        rsf_vars_ordered = rsf_vars_ordered,
        common_vars = common_vars,
        fg_only_vars = fg_only_vars,
        rsf_only_vars = rsf_only_vars,
        display_labels = all_display_labels
    ))
}


# compact 버전도 수정
plot_unified_heatmap_compact <- function(fg_results, rsf_results,
                               filename = "unified_feature_importance_auc_compact",
                               width = 10.5, height = 10) {

    # 기존 데이터 준비 코드와 동일
    fg_interpretation <- interpret_importance(fg_results, min_effect_size = 0.0001)
    rsf_interpretation <- interpret_importance(rsf_results, min_effect_size = 0.0001)

    fg_auc_data <- fg_interpretation$significant_auc
    fg_auc_data[is.na(fg_auc_data)] <- 0

    rsf_auc_data <- rsf_interpretation$significant_auc
    rsf_auc_data[is.na(rsf_auc_data)] <- 0

    all_fg_vars <- rownames(fg_auc_data)
    all_rsf_vars <- rownames(rsf_auc_data)
    
    fg_all_importance <- rowMeans(fg_auc_data[all_fg_vars, ], na.rm = TRUE)
    fg_vars_ordered <- names(sort(fg_all_importance, decreasing = TRUE))

    rsf_all_importance <- rowMeans(rsf_auc_data[all_rsf_vars, ], na.rm = TRUE)
    rsf_vars_ordered <- names(sort(rsf_all_importance, decreasing = TRUE))

    time_points <- colnames(fg_auc_data)
    
    fg_unique_names <- paste0("FG_", seq_along(fg_vars_ordered), "_", make.names(fg_vars_ordered))
    rsf_unique_names <- paste0("RSF_", seq_along(rsf_vars_ordered), "_", make.names(rsf_vars_ordered))
    
    fg_display_labels <- fg_vars_ordered
    rsf_display_labels <- rsf_vars_ordered

    all_unique_names <- c(fg_unique_names, rsf_unique_names)
    all_display_labels <- c(fg_display_labels, rsf_display_labels)
    
    unified_data <- matrix(0, nrow = length(all_unique_names), ncol = length(time_points))
    rownames(unified_data) <- all_unique_names
    colnames(unified_data) <- time_points
    
    unified_data[fg_unique_names, ] <- fg_auc_data[fg_vars_ordered, ]
    unified_data[rsf_unique_names, ] <- rsf_auc_data[rsf_vars_ordered, ]

    model_annotation <- data.frame(
        Model = c(rep("FG", length(fg_vars_ordered)),
                 rep("RSF", length(rsf_vars_ordered))),
        row.names = all_unique_names
    )

    annotation_colors <- list(
        Model = c("FG" = "#E31A1C", "RSF" = "#1F78B4")
    )

    # PDF와 TIFF 모두 저장
    save_heatmap_both_formats(
      unified_data = unified_data,
      model_annotation = model_annotation,
      annotation_colors = annotation_colors,
      all_display_labels = all_display_labels,
      fg_vars_length = length(fg_vars_ordered),
      base_filename = filename,
      width = width,
      height = height,
#      title = "Conditional Permutation Variable Importance\n(Ordered by Individual Model Importance)",
      fontsize_number = 10,
      cellwidth = 60,
      cellheight = 20,
      fontsize = 14,
      fontsize_row = 11,
      fontsize_col = 14
    )

    cat("Compact heatmap saved in both formats\n")
    
    return(list(
        unified_data = unified_data,
        model_annotation = model_annotation,
        fg_interpretation = fg_interpretation,
        rsf_interpretation = rsf_interpretation
    ))
}



# minimal 버전도 수정
plot_unified_heatmap_minimal <- function(fg_results, rsf_results,
                               filename = "unified_feature_importance_auc_minimal",
                               width = 11.5, height = 10) {

    # 데이터 준비는 동일...
    fg_interpretation <- interpret_importance(fg_results, min_effect_size = 0.0001)
    rsf_interpretation <- interpret_importance(rsf_results, min_effect_size = 0.0001)

    fg_auc_data <- fg_interpretation$significant_auc
    fg_auc_data[is.na(fg_auc_data)] <- 0

    rsf_auc_data <- rsf_interpretation$significant_auc
    rsf_auc_data[is.na(rsf_auc_data)] <- 0

    all_fg_vars <- rownames(fg_auc_data)
    all_rsf_vars <- rownames(rsf_auc_data)
    
    fg_all_importance <- rowMeans(fg_auc_data[all_fg_vars, ], na.rm = TRUE)
    fg_vars_ordered <- names(sort(fg_all_importance, decreasing = TRUE))

    rsf_all_importance <- rowMeans(rsf_auc_data[all_rsf_vars, ], na.rm = TRUE)
    rsf_vars_ordered <- names(sort(rsf_all_importance, decreasing = TRUE))

    time_points <- colnames(fg_auc_data)
    
    fg_unique_names <- paste0("FG_", seq_along(fg_vars_ordered), "_", make.names(fg_vars_ordered))
    rsf_unique_names <- paste0("RSF_", seq_along(rsf_vars_ordered), "_", make.names(rsf_vars_ordered))
    
    fg_display_labels <- fg_vars_ordered
    rsf_display_labels <- rsf_vars_ordered

    all_unique_names <- c(fg_unique_names, rsf_unique_names)
    all_display_labels <- c(fg_display_labels, rsf_display_labels)
    
    unified_data <- matrix(0, nrow = length(all_unique_names), ncol = length(time_points))
    rownames(unified_data) <- all_unique_names
    colnames(unified_data) <- time_points
    
    unified_data[fg_unique_names, ] <- fg_auc_data[fg_vars_ordered, ]
    unified_data[rsf_unique_names, ] <- rsf_auc_data[rsf_vars_ordered, ]

    model_annotation <- data.frame(
        Model = c(rep("FG", length(fg_vars_ordered)),
                 rep("RSF", length(rsf_vars_ordered))),
        row.names = all_unique_names
    )

    annotation_colors <- list(
        Model = c("FG" = "#E31A1C", "RSF" = "#1F78B4")
    )

    # PDF와 TIFF 모두 저장 (par 설정은 각 디바이스 내에서)
    # PDF 저장
    pdf_filename <- paste0(filename, ".pdf")
    if(require(cairo, quietly = TRUE)) {
        cairo_pdf(pdf_filename, width = width, height = height)
    } else {
        pdf(pdf_filename, width = width, height = height)
    }
    
    old_par <- par(no.readonly = TRUE)
    par(mar = c(0.5, 0.5, 1.5, 0.5), oma = c(0, 0, 0, 0))
    
    pheatmap(unified_data,
             cluster_rows = FALSE, cluster_cols = FALSE,
#             main = "Conditional Permutation Variable Importance",
             color = colorRampPalette(c("white", "#E6FFE6", "#98FB98", "#66CDAA"))(50),
#             color = colorRampPalette(c("white", "darkred"))(50),
             display_numbers = TRUE, number_format = "%.4f",
             fontsize_number = 12, cellwidth = 70, cellheight = 25,
             fontsize = 16, fontsize_row = 13, fontsize_col = 16,
             annotation_row = model_annotation, annotation_colors = annotation_colors,
             annotation_names_row = TRUE, labels_row = all_display_labels,
             gaps_row = c(length(fg_vars_ordered)), border_color = "grey60")
    
    par(old_par)
    dev.off()
    
    # TIFF 저장
    tiff_filename <- paste0(filename, ".tiff")
    tiff(tiff_filename, width = width, height = height, units = "in", res = 300, compression = "lzw")
    
    par(mar = c(0.5, 0.5, 1.5, 0.5), oma = c(0, 0, 0, 0))
    
    pheatmap(unified_data,
             cluster_rows = FALSE, cluster_cols = FALSE,
#             main = "Conditional Permutation Variable Importance",
             color = colorRampPalette(c("white", "#E6FFE6", "#98FB98", "#66CDAA"))(50),
#             color = colorRampPalette(c("white", "darkred"))(50),
             display_numbers = TRUE, number_format = "%.4f",
             fontsize_number = 12, cellwidth = 70, cellheight = 25,
             fontsize = 16, fontsize_row = 13, fontsize_col = 16,
             annotation_row = model_annotation, annotation_colors = annotation_colors,
             annotation_names_row = TRUE, labels_row = all_display_labels,
             gaps_row = c(length(fg_vars_ordered)), border_color = "grey60")

    dev.off()

    cat("Minimal margin heatmap saved as PDF:", pdf_filename, "\n")
    cat("Minimal margin heatmap saved as TIFF:", tiff_filename, "\n")
    
    return(unified_data)
}


# 실행
# 방법 1: 적당히 여백 줄이기
result1 <- plot_unified_heatmap_compact(cpi_fg_results, cpi_rsf_results,
					filename="unified_cpi_vimp_auc_compact",
					width = 10.5, height = 10)

# 방법 2: 최대한 여백 줄이기
result2 <- plot_unified_heatmap_minimal(cpi_fg_results, cpi_rsf_results,
					filename = "unified_cpi_vimp_auc_minimal",
					width = 11.5, height = 10)

# 히트맵 생성 실행
result <- plot_unified_heatmap(cpi_fg_results, cpi_rsf_results,
                              filename = "unified_cpi_vimp_auc_standard",
                              width = 10, height = 7.5)
