# Short title: LI-RADS-Based Competing Risks in Liver Transplant

> Full title and manuscript details will be disclosed upon journal acceptance.
> **Manuscript in preparation** — Code is provided for reproducibility purposes.


---

## Overview

This repository contains the R code used to develop and validate the **LT-RADS model**, a competing risks-based survival prediction framework for HCC-related mortality following liver transplantation (LT). Pre-transplant clinical and LI-RADS imaging data were collected from LT recipients, and competing risk models were applied to predict HCC-specific mortality under the competing event of non-HCC death.

Two model classes were developed and compared:
- **Fine-Gray (FG)** subdistribution hazard regression
- **Random Survival Forest (RSF)** with competing risks

---

## Requirements

R version 4.4.3

| Package | Version |
|---|---|
| riskRegression | 2025.5.20 |
| randomForestSRC | 3.3.3 |
| survival | 3.8.3 |
| prodlim | 2025.4.28 |
| timeROC | 0.4 |
| pracma | 2.4.4 |
| dplyr | 1.1.4 |
| ggplot2 | 3.5.2 |
| car | 3.1.3 |
| corrplot | 0.95 |

---

## Repository Structure

```
.
├── FG/                                         # Fine-Gray model
│   ├── main_fg.r                               # Model training
│   ├── optim_fg_backward_eliminate_custom_penalty.r  # Variable selection
│   ├── eval_fg_boot.r                          # Bootstrap evaluation
│   ├── fg_conditional_perm_vimp_tdAUC.r        # Variable importance (tdAUC)
│   ├── fg_conditional_perm_vimp_cindex.r       # Variable importance (C-index)
│   ├── compute_sHR.r                           # Subdistribution hazard ratio
│   ├── plus_ald_subanalysis.r                  # ALD subgroup analysis
│   └── plus_ald_subanalysis_boot.r             # ALD subgroup bootstrap
│
├── RSF/                                        # Random Survival Forest model
│   ├── main_rsf.r                              # Model training
│   ├── optim_rsf_backward_eliminate.r          # Variable selection
│   ├── eval_rsf_boot.r                         # Bootstrap evaluation
│   ├── rsf_conditional_perm_vimp_tdAUC.r       # Variable importance (tdAUC)
│   ├── rsf_conditional_perm_vimp_cindex.r      # Variable importance (C-index)
│   └── minus_ald_subanalysis.r                 # ALD subgroup analysis
│
├── generate_bootstrap_indices.r                # Shared bootstrap indices
├── main_{milan,mt2,uts}_external.r             # External validation (criteria-based models)
├── eval_{milan,mt2,uts}_boot.r                 # Bootstrap eval for criteria models
├── calibration_curves_5models.r                # Calibration curves
├── auroc_comparison_figure.r                   # tdAUC comparison figure
├── pairwise_bootstrap_comparison_figure.r      # Pairwise p-value heatmap
├── plot_ext_ROC_curve.r                        # External ROC curves
├── plot_conditional_perm_vimp_tdAUC.r          # Unified VIMP plot (tdAUC)
├── plot_conditional_perm_vimp_cindex.r         # Unified VIMP plot (C-index)
├── create_nomogram.r                           # Nomogram for FG model
├── wolbers_cindex_cr_claude.r                  # Wolbers C-index for competing risks
```

---

## Data Availability

The clinical dataset used in this study contains patient-level information and **cannot be shared publicly** due to patient privacy regulations and institutional data governance policies.

For data access inquiries, please contact the corresponding author.

### Input Variables (n=39)

**Categorical (n=21)**

| Variable | Description |
|---|---|
| sex | Sex |
| diabetes | Diabetes mellitus |
| ckd_hd | Chronic kidney disease |
| hbv | Hepatitis B virus etiology |
| hcv | Hepatitis C virus etiology |
| image_lc | Imaging-based liver cirrhosis |
| ald | Alcoholic liver disease |
| masld | Metabolic dysfunction-associated steatotic liver disease |
| ldlt | Living donor liver transplantation |
| pre_tx | Treatment experience|
| tace | Transarterial chemoembolization |
| rfa | Radiofrequency ablation |
| resection | Resection |
| rtx | Radiation therapy |
| lr_5_ct | Number of LR-5 (CT) |
| lr_tr_v_ct | Number of LR-TR V (CT) |
| lr_3_ct | Number of LR-3 (CT) |
| lr_4_ct | Number of LR-4 (CT) |
| lr_m_ct | Number of LR-M (CT) |
| lr_tr_nv_ct | Number of LR-TR-NV (CT) |
| lr_tr_e_ct | Number of LR-TR-E (CT) |

**Continuous (n=18, log-transformed where indicated)**

| Variable | Description | Transformation |
|---|---|---|
| log_lr_5_size_ct | Maximum diabetes of LR-5 (CT) | log |
| log_lr_tr_v_size_whole_ct | Diameter of largest viable nodule by RECIST v1.1 (CT) | log |
| log_lr_tr_v_size_enhancing_ct | Maximum diabetes of LR-TR-V (CT) | log |
| log_lr_3_size_ct | Maximum diabetes of LR-3 (CT) | log |
| log_lr_4_size_ct | Maximum diabetes of LR-4 (CT) | log |
| log_lr_m_size_ct | Maximum diameter of LR-M (CT) | log |
| log_lr_tr_nv_size_ct | Maximum diameter of LR-TR-NV (CT) | log |
| log_lr_tr_e_size_ct | Maximum diameter of LR-TR-E (CT) | log |
| log_afp | Alpha-fetoprotein | log |
| log_pivka | PIVKA-II | log |
| alb | Serum albumin | — |
| log_tbil | Total bilirubin | log |
| cr | Creatinine | — |
| log_pt_inr | PT-INR | log |
| na | Serum sodium | — |
| log_crp | C-reactive protein | log |
| age_lt | Age at liver transplantation | — |
| bmi | Body mass index | — |

### Outcome Variables

| Variable | Description |
|---|---|
| ftime | Time to event (days from transplantation) |
| fstatus | 1 = HCC-related death, 2 = Non-HCC death (competing), 0 = Censored |

### Benchmark Criteria

- Milan criteria
- Up-to-Seven criteria
- Metroticket 2.0 score

### RSF Hyperparameter Grid

| Parameter | Values searched |
|---|---|
| ntree | 500, 1000 |
| nodesize | 5, 10, 15 |

---

## Usage

### 1. Model Training

```r
# Fine-Gray model
source("FG/main_fg.r")

# RSF model
source("RSF/main_rsf.r")
```

### 2. Variable Selection

```r
# Fine-Gray: penalized backward elimination
source("FG/optim_fg_backward_eliminate_custom_penalty.r")

# RSF: backward elimination
source("RSF/optim_rsf_backward_eliminate.r")
```

### 3. Bootstrap Evaluation (internal & external)

```r
# Generate shared bootstrap indices
source("generate_bootstrap_indices.r")

# Bootstrap evaluation
source("FG/eval_fg_boot.r")
source("RSF/eval_rsf_boot.r")

# Criteria-based models (Milan, UCSF MT2.0, UpToSeven)
source("eval_milan_boot.r")
source("eval_mt2.0_boot.r")
source("eval_uts_boot.r")
```

### 4. External Validation

```r
source("main_milan_external.r")
source("main_mt2_external.r")
source("main_uts_external.r")
```

### 5. Evaluation Metrics

```r
# Wolbers C-index for competing risks
source("wolbers_cindex_cr_claude.r")
```

### 6. ROC & tdAUC Comparison

```r
# External ROC curves across time points
source("plot_ext_ROC_curve.r")

# tdAUC pairwise bootstrap comparison
source("auroc_comparison_figure.r")
source("pairwise_bootstrap_comparison_figure.r")
```

### 7. Calibration Curves

```r
source("calibration_curves_5models.r")
```

### 8. Variable Importance (Conditional Permutation)

```r
source("plot_conditional_perm_vimp_tdAUC.r")
source("plot_conditional_perm_vimp_cindex.r")
```

### 9. Nomogram

```r
source("create_nomogram.r")
```

### 10. Subdistribution Hazard Ratio (sHR) Analysis

```r
source("FG/compute_sHR.r")
```

---

## Results

Figures and aggregate results will be made available upon journal acceptance.

Key outputs generated by this codebase:
- Time-dependent ROC curves (1, 3, 5, 10 years) with DeLong's test
- Calibration curves 
- Conditional permutation variable importance (C-index, tdAUROC)
- Nomogram (Fine-Gray model)
- Pairwise bootstrap comparison heatmaps

---

## Contact

For questions regarding this study or data access, please contact:

**Corresponding Authors**: 

Namkug Kim, PhD
Department of Radiology and Research Institute of Radiology, Asan Medical Center, University of Ulsan College of Medicine, Seoul, Republic of Korea.
Department of Convergence Medicine, Asan Medical Center, University of Ulsan College of Medicine, Seoul, Korea
E-mail: namkugkim@gmail.com
ORCID: [0000-0002-3438-2217](https://orcid.org/0000-0002-3438-2217)

Ju Hyun Shim, MD, PhD
Department of Gastroenterology, Asan Medical Center, University of Ulsan College of Medicine, Seoul, Republic of Korea
E-mail: s5854@amc.seoul.kr
ORCID: [0000-0002-7336-1371](https://orcid.org/0000-0002-7336-1371)

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
