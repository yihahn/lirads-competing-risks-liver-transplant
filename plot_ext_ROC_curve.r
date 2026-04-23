# 필요한 패키지 로드
library(ggplot2)
library(dplyr)
library(scales)

# 모델별 파일 이름과 해당 모델의 이름을 매핑합니다.
model_files <- c("FG/ext_fg_roc_by_times.rds"="FG",
		 "RSF/ext_rsf_roc_by_times.rds"="RSF",
		 "ext_milan_roc_by_times.rds"="Milan",
		 "ext_mt2_roc_by_times.rds"="Mt2.0",
		 "ext_uts_roc_by_times.rds"="UtS"
)

# 모든 ROC 데이터를 담을 리스트 초기화
all_roc_data <- list()

# 각 모델 파일에서 데이터 로드 및 전처리
for (file_path in names(model_files)) {
	model_name <- model_files[[file_path]] # 파일 이름에 해당하는 모델 이름 가져오기

	# RDS 파일 로드
	# 파일이 현재 작업 디렉토리에 없으면 전체 경로를 지정해야 합니다.
	current_model_data <- readRDS(file_path)

	# 각 연도별 ROC 데이터 처리
	for (i in seq_along(current_model_data$tpr)) {
		year_label <- names(current_model_data$tpr)[i]

		 # t 값 추출 (예: "t=365.24" -> "365.24")
		t_value <- as.numeric(gsub("t=", "", year_label))

		# 연도를 '년' 단위로 변환 (1년, 3년, 5년, 10년)
		year_in_years <- round(t_value / 365.24)  

		# TPR, FPR 벡터 추출 (List of 1 -> num 벡터)
		tpr_vec <- current_model_data$tpr[[i]][[1]]
		fpr_vec <- current_model_data$fpr[[i]][[1]]
    
		# 데이터 프레임 생성
		roc_df <- data.frame(tpr = tpr_vec,
				     fpr = fpr_vec,
				     model = model_name, # 현재 모델 이름 할당
				     Year = ifelse(year_in_years == 1, "1 Year", paste0(year_in_years, " Years")) # 연도 (예: "1 year", "3 years" 등)
		)  

		# 데이터를 all_roc_data 리스트에 추가 (모델 이름과 연도별로 구분)
		# 리스트의 이름을 사용하여 나중에 쉽게 결합할 수 있도록 합니다.
		all_roc_data[[paste0(model_name, "_", year_in_years, "years")]] <- roc_df
	}
}

# 모든 모델과 연도의 ROC 데이터를 하나의 데이터 프레임으로 결합
final_roc_data <- bind_rows(all_roc_data)

# 'Year' 컬럼을 factor로 변환하여 순서대로 정렬합니다.
# 이렇게 하면 facet_wrap에서 1년, 3년, 5년, 10년 순서로 그래프가 그려집니다.
year_order <- c("1 Year", "3 Years", "5 Years", "10 Years")
final_roc_data$Year <- factor(final_roc_data$Year, levels=year_order)

# 모델 순서 설정 - FG, RSF, Milan, UTS, Mt2 순서로 변경
model_order <- c("FG", "RSF", "Milan", "UtS", "Mt2.0")
n <- length(model_order)
colors <- hue_pal()(n)
names(colors) <- model_order
print(colors)

# hex code를 색상 이름으로 변환
for (nm in names(colors)) {
  cat(sprintf("%s: %s\n", nm, colors[nm]))
}

final_roc_data$model <- factor(final_roc_data$model, levels=model_order)

# 데이터 확인 (디버깅용)
cat("Unique years in data:", unique(final_roc_data$Year), "\n")
cat("Number of rows per year:\n")
print(table(final_roc_data$Year))
cat("Unique models in data:", unique(final_roc_data$model), "\n")

# CSV 파일에서 AUROC 데이터 읽기
auroc_csv <- read.csv("models_auroc_data.csv", stringsAsFactors = FALSE)
print("CSV 데이터:")
print(auroc_csv)

# CSV 데이터를 long format으로 변환
auroc_results <- data.frame()

for(i in 1:nrow(auroc_csv)) {
  model_name <- auroc_csv$model[i]
  
  # 각 연도별 데이터 처리
  years_data <- data.frame(
    model = rep(model_name, 4),
    Year = c("1 Year", "3 Years", "5 Years", "10 Years"),
    AUROC = c(auroc_csv$X1year[i], auroc_csv$X3years[i], 
              auroc_csv$X5years[i], auroc_csv$X10years[i]),
    stringsAsFactors = FALSE
  )
  
  auroc_results <- rbind(auroc_results, years_data)
}

# y 위치 설정 - FG, RSF, Milan, UTS, Mt2 순서에 맞게 조정
auroc_results$y_pos <- ifelse(auroc_results$model == "FG", 0.25,
                       ifelse(auroc_results$model == "RSF", 0.20,
                       ifelse(auroc_results$model == "Milan", 0.15,
                       ifelse(auroc_results$model == "UtS", 0.10, 0.05))))

# 텍스트 라벨 생성
auroc_results$AUROC_text <- sprintf("%.3f", auroc_results$AUROC)
auroc_results$label <- paste0(auroc_results$model, ": ", auroc_results$AUROC_text)

# factor 변환 (원본 데이터와 동일한 레벨 사용)
# CSV의 모델명을 원본 데이터의 모델명과 맞춤
#auroc_results$model <- ifelse(auroc_results$model == "Mt2", "mt2", auroc_results$model)
auroc_results$model <- factor(auroc_results$model, levels = levels(final_roc_data$model))
auroc_results$Year <- factor(auroc_results$Year, levels = levels(final_roc_data$Year))

# AUROC 데이터 확인
print("최종 AUROC 결과:")
print(auroc_results)


# ROC curve 그리기
roc_plot <- ggplot(final_roc_data, aes(x=fpr, y=tpr, color=model)) +
	geom_line(linewidth=0.6, alpha=0.6) + # ROC 곡선 (선 두께 증가, 투명도 추가)
	geom_abline(intercept=0, slope=1, linetype="dashed", color="gray50", linewidth=0.5) + # 기준선 (랜덤 예측 선)
	# AUROC 값 텍스트 추가
	geom_text(data=auroc_results, aes(x=0.87, y=y_pos, label=label, color=model), 
	         size=3, hjust=0, fontface="bold", show.legend=FALSE) +
	labs(title="ROC curves of external validation cohort across Years",
	     x="1 - Specificity (False Positive Rate)",
	     y="Sensitivity (True Positive Rate)",
	     color="Model" # 범례 제목
	     ) +
	scale_x_continuous(limits=c(-0.012, 1.012), breaks=seq(0, 1, by=0.1), expand=c(0, 0)) +
	scale_y_continuous(limits=c(-0.012, 1.012), breaks=seq(0, 1, by=0.1), expand=c(0, 0)) +
	facet_wrap(~ Year, scales="fixed", nrow=2, ncol=2) + # 연도별로 2x2 그림으로 분리
	theme_minimal() + # 깔끔한 테마 적용
	theme(panel.background=element_blank(),
	      panel.grid.major=element_line(color="gray90", linewidth=0.3),
	      panel.grid.minor=element_blank(),
	      panel.border=element_rect(color="black", fill=NA, linewidth=0.8),
	      axis.line=element_line(color="black", linewidth=0.5),
	      axis.ticks=element_line(color="black", linewidth=0.5),
	      axis.ticks.length=unit(0.1, "cm"),
	      plot.title=element_text(hjust=0.5, face="bold", size=16), # 제목 가운데 정렬 및 글꼴 설정
	      legend.position="bottom", # 범례를 아래쪽에 배치
	      legend.title=element_text(face="bold", size=13), # 범례 제목 굵게
	      legend.text=element_text(size=11),
	      strip.text=element_text(face="bold", size=13), # facet 제목 (연도) 굵게 및 크기 설정
	      strip.background=element_rect(fill="gray95", color="black", linewidth=0.5)
	)

tiff("ROC_curves_External_dataset_across_years.tiff",
     width=12, height=10, units="in", res=300, compression="lzw")
print(roc_plot)
dev.off()

pdf("ROC_curves_External_dataset_across_years.pdf",
     width=12, height=10)
print(roc_plot)
dev.off()

