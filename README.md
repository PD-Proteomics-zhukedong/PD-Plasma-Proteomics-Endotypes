# PD-Plasma-Proteomics-Endotypes

**Reproducible R computational workflows and clinical decision-support portal for plasma proteomic stratification in Parkinson's disease (*Nature Aging*)**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D%204.3.1-blue.svg)](https://www.r-project.org/)
[![Decision Portal](https://img.shields.io/badge/Shiny-Live%20Portal-success.svg)](https://zhukedong.shinyapps.io/PD-Plasma-Proteomics-Portal/)

---

## 📖 Overview

This repository provides the complete, leak-free computational pipeline and interactive decision-support application developed for our study:

> **Large-scale plasma proteomics identifies molecularly distinct PD endotypes associated with motor and non-motor phenotypes.** *Nature Aging*, 2026.

By integrating deep Orbitrap Astral Data-Independent Acquisition mass spectrometry (DIA-MS) across a dual-center East Asian cohort ($N = 1,119$), this project establishes a consensus signature of 823 strictly validated Meta-Differentially Expressed Proteins (Meta-DEPs), defines a parsimonious 16-protein diagnostic model, resolves Parkinson's disease into three reproducible molecular endotypes independent of pharmacological and chronicity confounders, and evaluates translational surrogate classifiers within international external benchmark cohorts ($N = 868$).

---

## 📁 Repository Structure

```text
PD-Plasma-Proteomics-Endotypes/
├── scripts/                      # Sequential, leak-free analytical scripts (01 to 09)
│   ├── 01_preprocessing_imputation_split_combat.R
│   ├── 02_differential_expression_and_meta_analysis.R
│   ├── 03_diagnostic_classifier_lasso.R
│   ├── 04_external_multicohort_replication.R
│   ├── 05_gtex_tissue_enrichment.R
│   ├── 06_consensus_clustering_and_endotype_mapping.R
│   ├── 07_subtype_mechanisms_and_cellular_deconvolution.R
│   ├── 08_subtype_classifier_and_clinical_prediction.R
│   └── 09_multicohort_endotype_validation.R
├── shiny_app/                    # Production code for clinical decision-support portal
│   ├── app.R                     # Bilingual UI and predictive inference engine
│   ├── models/                   # Serialized model artifacts (.rds)
│   └── data/                     # De-identified demonstration cohort data
├── LICENSE                       # MIT Open Source License
└── README.md                     # Pipeline documentation and execution guide
🛠️ Sequential Computational Workflow (scripts/)
The pipeline is organized sequentially. Running scripts 01 through 09 end-to-end reproduces all primary and supplementary findings reported in the manuscript:
Script 01: Preprocessing & Split-ComBat Harmonization
Workflow: 30% completeness filtering, sample-wise down-shifted normal imputation (MNAR), and cohort-independent Split-ComBat empirical Bayes batch harmonization.
Outputs: Normalized quantification matrix (3,946 reliable plasma proteins); Supplementary Figure S2A & S2B.
Script 02: Multivariable DEA & Stouffer Meta-Analysis
Workflow: Two-stage multivariable linear modeling (limma; Basic and 6-covariate Full models), sample-size-weighted Stouffer meta-analysis, background-constrained GO-BP/Reactome enrichment, STRING topological centrality profiling, and covariate-adjusted ssGSEA module burden regressions.
Outputs: Figure 1B–G; Supplementary Tables ST4, ST5, ST9, ST10.
Script 03: 16-Protein Clinical Offset Diagnostic Classifier
Workflow: Unpenalized demographic logit offset modeling, 98% relative efficiency LASSO feature selection on the discovery cohort, independent validation testing, SHAP beeswarm feature attribution, and interactome connectivity mapping.
Outputs: Figure 2A–D; Supplementary Table ST12; trained lasso_diagnostic_model.rds.
Script 04: External Multi-Cohort Concordance Benchmarks
Workflow: Cross-platform effect-size concordance evaluation across PPMI Project 314 (Olink Explore HT), STTT 2026 East Asian cohort (Olink Explore 384 Inflammation), and UK Biobank prospective pre-diagnostic incident cohort (Olink Explore 3072).
Outputs: Supplementary Tables ST6, ST7, ST8.
Script 05: GTEx v8 Transcriptomic Tissue-of-Origin Deconvolution
Workflow: Background-constrained two-sided hypergeometric over-representation tests across 54 GTEx human tissues, characterizing circulating perturbations against central and peripheral organ signatures.
Outputs: Supplementary Figure 3; Supplementary Table ST11.
Script 06: Confounder-Residualized Consensus Clustering & Subtyping
Workflow: 8-covariate multivariable residualization (decoupling LEDD, disease duration, age, sex, BMI, and comorbidities), unsupervised PAM consensus clustering, independent scale-driven clinical clustering, cross-modal Sankey mapping, and dual-model ANCOVA across clinical domains.
Outputs: Figure 3A–I; Supplementary Figure 4A–C; Supplementary Tables ST13, ST14, ST15.
Script 07: Subtype Markers & Cellular Microenvironment Deconvolution
Workflow: 8-covariate pairwise Limma contrasts, bottleneck Min-t exclusivity scoring, topological GO.db directed acyclic graph pruning, PanglaoDB 126 cell lineage deconvolution, and multivariable partial correlation heatmaps with clinical scales.
Outputs: Figure 4A–F; Supplementary Figure 5A–B; Supplementary Tables ST16, ST17, ST18, ST19, ST20.
Script 08: Balanced Random Forest Classifier & Clinical Trait Predictions
Workflow: Multi-compartment STRING PPI network reconstructions (Fig 5A–C), Balanced Random Forest classifier training on the discovery cohort, ROC evaluation in the independent validation cohort, 10-fold cross-validated regression for 7 continuous clinical scales, and model artifact serialization.
Outputs: Figure 5A–F; Supplementary Tables ST19, ST21; production model objects for decision portal integration.
Script 09: External PPMI Cohort Validation & Rank-ANCOVA Gradients
Workflow: Out-of-sample Bayesian prior-adjusted projection onto PPMI Project 314, population proportion stability testing, and unified Conover-Iman non-parametric Rank-ANCOVA evaluating motor (UPDRS-III adjusted for LEDD) and olfactory (UPSIT) phenotypic gradients under planned directional replication hypotheses.
Outputs: Figure 6A, 6B; Supplementary Tables ST22, ST23, ST24.
🌐 Interactive Clinical Decision Support Portal
An interactive web-based implementation of the diagnostic model, molecular endotyping classifier, and multi-domain symptom predictors is deployed for peer review:
🔗 Portal URL: https://zhukedong.shinyapps.io/PD-Plasma-Proteomics-Portal/
🔑 Reviewer access passcode has been provided in the Response to Reviewers document.
Core Modules:
Individual Patient Profiler: Computes real-time 16-protein diagnostic probabilities, predicts molecular endotypes, estimates 7 continuous symptom scores, and provides biology-guided precision intervention pathways.
Cohort Batch Profiling: High-throughput automated classification and 7-dimensional phenotypic radar visualization from user-provided proteomic tables (CSV).
Biomarker Directory: Searchable reference catalog of all 60 unique panel biomarkers mapped to UniProt IDs, gene symbols, and functional annotations.
💻 System Requirements & Dependencies
The pipeline has been developed and validated on R (version >= 4.3.1) across Windows, macOS, and Linux.
Key R Packages:
Data Wrangling & I/O: tidyverse, data.table, readxl, arrow, broom
Differential Expression & Batch Harmonization: limma, sva
Machine Learning & Modeling: glmnet, randomForest, pROC, caret, kernelshap, shapviz, emmeans, ppcor
Clustering & Deconvolution: ConsensusClusterPlus, GSVA, mclust
Network & Functional Enrichment: clusterProfiler, org.Hs.eg.db, ReactomePA, GO.db, igraph, tidygraph, ggraph, ggforce
Data Visualization: ggplot2, ComplexHeatmap, circlize, ggpubr, ggalluvial, UpSetR, patchwork, cowplot
📄 License & Citation
This repository is licensed under the MIT License - see the LICENSE file for terms.
If you utilize this workflow, models, or data in your research, please cite:
Zhu, K. et al. Large-scale plasma proteomics identifies molecularly distinct PD endotypes associated with motor and non-motor phenotypes. Nature Aging (2026).
