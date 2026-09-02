**# Plasma Proteomic Stratification \& Clinical Decision Portal for Parkinson's Disease**



\[!\[R Version](https://img.shields.io/badge/R-%3E%3D4.3.1-blue.svg)](https://www.r-project.org/)

\[!\[License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

\[!\[Journal](https://img.shields.io/badge/Publication-Nature%20Aging%20(2026)-red.svg)](https://www.nature.com/)



An interactive clinical decision-support and translational stratification platform developed using high-throughput Orbitrap Astral Data-Independent Acquisition (DIA) mass spectrometry across a large dual-center East Asian cohort (\*N\* = 1,119).



\---



\## 🌟 Overview \& Key Methodological Features



This web portal translates blood-based molecular endotyping and machine learning models into an operational decision-support tool:

1\. \*\*16-Protein Clinical Offset Diagnostic Model\*\*: Estimates the individualized probability of Parkinson’s disease (PD vs. Healthy Controls) incorporating an unpenalized 6-covariate clinical offset (Age, Sex, BMI, Hepatic disease, Renal impairment, and Cardiovascular/metabolic comorbidities).

2\. \*\*45-Biomarker Molecular Endotype Classifier\*\*: Stratifies sporadic PD patients into three conserved biological endotypes derived from 8-covariate residualized proteomics (stripping medication dosage and disease chronicity):

&#x20;  - \*\*Subtype 1 (Integrin \& Autophagy Axis)\*\*: Characterized by cell adhesion, focal adhesion, and autophagic proteostasis alongside mild motor/non-motor burden.

&#x20;  - \*\*Subtype 2 (Neuroinflammatory \& Humoral Immune Axis)\*\*: Characterized by glial activation (astrocytes/microglia), complement activation, and BBB hyperpermeability, presenting with high motor severity (UPDRS-III) and depressive burden (HAMD).

&#x20;  - \*\*Subtype 3 (Mitochondrial \& Metabolic Axis)\*\*: Characterized by mitochondrial Complex I impairment, ATP synthesis failure, and BCAA/lipid metabolic reprogramming, presenting with pronounced, selective olfactory loss (UPSIT < 10).

3\. \*\*Multi-Domain Continuous Symptom Prediction (10-Fold CV)\*\*: Predicts continuous motor (UPDRS-III), disease staging (H\&Y), depression (HAMD), non-motor total (NMSS), cognitive (MMSE), sleep (PDSS), and olfactory (UPSIT) scores.

4\. \*\*Translational Precision Medicine Roadmap (Figure 6d)\*\*: Dynamically generates subtype-stratified adjunct therapeutic recommendations.

5\. \*\*Universal Identifier \& Automatic Standardization (Reviewer #2 Focus)\*\*:

&#x20;  - Supports \*\*both official HGNC Gene Symbols (e.g., `ITGA6`, `VTN`) and UniProt Accessions (e.g., `P23229`, `P04004`)\*\* with automated recognition.

&#x20;  - Supports input of \*\*Raw MS Log2 Intensities\*\* (automatically standardized against the discovery reference cohort parameters: $Z = \\frac{X - \\mu}{\\sigma}$) or \*\*Standardized Z-scores\*\*.



\---



\## 📁 Repository Structure



```text

PD\_Proteomics\_ShinyApp/

├── app.R                          # Main Shiny application source code

├── README.md                      # Comprehensive user manual and documentation

├── models/

│   ├── diag\_offset\_model.rds      # 16-protein logistic offset diagnostic model

│   ├── subtype\_rf\_classifier.rds  # 45-biomarker balanced random forest subtyping model

│   ├── clinical\_rf\_predictors.rds # 7 continuous clinical symptom regression models

│   └── app\_metadata.rds           # Feature parameters and dual-identifier mapping dictionary

└── data/

&#x20;   └── example\_test\_samples.csv   # Pre-packaged benchmark demonstration cohort



**💻 Local Installation \& Quick Start**

1\. Prerequisites

Ensure you have R (≥ 4.3.0) and RStudio installed.

2\. Install Required Packages

In the R console, execute:

install.packages(c("shiny", "bslib", "tidyverse", "randomForest", "DT", "plotly", "shinyWidgets", "readxl"))

3\. Launch the Application

Clone this repository or download the folder, then run:

library(shiny)

\# Set working directory to the ShinyApp folder

shiny::runApp("path/to/PD\_Proteomics\_ShinyApp")



**📊 User Guide \& Data Input Specifications**

Module 1: Individual Patient Profiler (Default View)

Clinical Input: Enter demographics (Age, Sex, BMI), comorbidities (Hepatic, Renal, Cardiovascular/Metabolic 0/1/2+), disease duration, and levodopa dose (LEDD).

Proteomic Input:

Option A (File Upload / Paste Table): Upload a 1-row CSV or paste two columns (GeneSymbol, Intensity) from Excel. The system automatically extracts all 61 target markers.

Option B (Manual Form Entry): Expand the accordions to enter values across the 16 Diagnostic Markers, Subtype 1 Markers (15), Subtype 2 Markers (15), and Subtype 3 Markers (15).

Output: Real-time diagnostic risk probability (PD vs. HC), predicted molecular endotype(P(S1),P(S2),P(S3)), continuous symptom scores (UPDRS-III, HAMD, UPSIT), and tailored therapeutic roadmap cards.

Module 2: Cohort-Level Batch Profiling

Input Matrix: Upload a standard CSV file where each row represents a patient and columns contain clinical variables (Age, Sex, BMI, LEDD, Duration, hepatic\_disease, kidney.function, CV\_Metabolic\_Cat) and protein markers (Gene Symbols or UniProt IDs).

Demonstration: Click "Load Benchmark Example Cohort" to test 6 multi-subtype reference patients.

Output: Summary table with individualized predictions, interactive endotype distribution pie chart, multi-domain severity radar chart, and downloadable full CSV prediction report.

Module 3: Biomarker Directory

A searchable reference catalog mapping all 61 unique proteins across UniProt Accessions, Gene Symbols, Panel Classifications (Subtype 1, Subtype 2, Subtype 3, Diagnostic Panel Only, or Subtype 1 \& Diagnostic (Dual Role)).



**📄 Citation \& Attribution**

If you utilize this portal, models, or data in your research, please cite:

Large-scale plasma proteomics identifies molecularly distinct PD endotypes associated with motor and non-motor phenotypes. Nature Aging, 2026.



**⚠️ Disclaimer**

This portal is intended solely for scientific research, biological exploration, and translational biomarker evaluation. It does not constitute standalone clinical diagnostic advice or prescription guidelines.

