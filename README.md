# Large-scale plasma proteomics identifies molecularly distinct PD endotypes associated with motor and non-motor phenotypes

This repository contains the complete, reproducible computational workflow and the interactive clinical decision portal developed for our study published in *Nature Aging*.

---

## 📁 Repository Structure

- **`scripts/`**: Sequential, reproducible R scripts covering the full analytical workflow:
  - `01_preprocessing_imputation_split_combat.R`: Data preprocessing, 30% completeness filtering, MNAR imputation, and Split-ComBat batch correction (Supplementary Figure S2B).
  - `02_differential_expression_and_meta_analysis.R`: Multivariable limma linear modeling (6 covariates), Stouffer's meta-analysis, and Figure 1 generation (Panels B–G; Tables ST9, ST10, ST11, ST13).
  - `03_diagnostic_classifier_lasso.R`: 16-protein clinical-offset penalized LASSO classifier and Figure 2 generation (Panels A–D; Table ST14a).
  - `04_external_multicohort_replication.R`: Cross-platform replication across PPMI P314, STTT 2026, and UK Biobank cohorts (Supplementary Figure 3; Tables ST19, ST20, ST21).
  - `05_ppmi_diagnostic_validation.R`: Independent cross-platform diagnostic validation of the 16-protein panel across PPMI cohorts (Table ST14b).
  - `06_gtex_tissue_enrichment.R`: GTEx v8 54-tissue transcriptomic deconvolution and clinical partial correlation heatmap (Supplementary Figure S4; Table ST12).
  - `07_consensus_clustering_and_endotype_mapping.R`: 8-covariate residualized consensus clustering, molecular endotyping, and Figure 3 / Supplementary Figure 5 (Tables ST7, ST8, ST15).
  - `08_subtype_mechanisms_and_cellular_deconvolution.R`: Min-T biomarker specificity, PanglaoDB 126 cell-type deconvolution, and Figure 4 / Supplementary Figure 6 (Tables ST16, ST17).
  - `09_subtype_classifier_and_clinical_prediction.R`: Multi-class balanced random forest subtype classifier and 10-fold CV clinical trait predictions (Figure 5; Table ST18).
  - `10_multicohort_endotype_validation.R`: Cross-cohort validation of subtype proportions and phenotypic gradients across 4 independent cohorts (Figure 6A–C; Table ST22).

- **`shiny_app/`**: Standalone source code, trained model artifacts, and test data for the interactive clinical decision portal (`app.R`).

---

## 🌐 Interactive Web Application
Access the live decision support portal online without local installation:
👉 **[Launch Interactive Decision Portal](https://yourlab.shinyapps.io/PD-Proteomics-Predictor/)**

---

## 💻 System Requirements & Dependencies
The analysis pipeline has been tested in R (version >= 4.3.0) on Windows, macOS, and Linux.
Key package dependencies:
- `tidyverse`, `data.table`, `limma`, `sva`, `ConsensusClusterPlus`, `glmnet`, `randomForest`, `pROC`, `ComplexHeatmap`, `clusterProfiler`, `org.Hs.eg.db`, `GSVA`, `ppcor`, `ggraph`, `igraph`, `ggalluvial`.

---

## 📄 License
This project is licensed under the MIT License.# PD-Plasma-Proteomics-Endotypes
Reproducible R code and interactive Shiny portal for plasma proteomic stratification in Parkinson's disease (Nature Aging)
