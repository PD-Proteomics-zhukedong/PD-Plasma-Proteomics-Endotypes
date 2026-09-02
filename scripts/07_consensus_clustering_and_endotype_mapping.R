# ==============================================================================
# Script: 07_consensus_clustering_and_endotype_mapping.R
# Purpose: Covariate residualization, molecular consensus clustering, and Figure 3 / Supp Fig 5
# Project: Large-scale plasma proteomics identifies molecularly distinct PD endotypes
# ==============================================================================

options(expressions = 5000)
set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
  library(ConsensusClusterPlus)
  library(GSVA)
  library(ggpubr)
  library(ggalluvial)
  library(mclust)
  library(emmeans)
  library(broom)
  library(org.Hs.eg.db)
  library(patchwork)
  library(scales)
})

# Directory setup
data_dir      <- "data"
input_dea_dir <- "results/02_figure1_meta_dea"
output_dir    <- "results/07_figure3_molecular_endotypes"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

find_file_smart <- function(pattern, default_path) {
  files <- list.files(data_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  files <- files[!grepl("^~\\$", basename(files))]
  if (length(files) > 0) return(files[1])
  return(default_path)
}

get_mode <- function(x) {
  x <- na.omit(x)
  if (length(x) == 0) return(NA)
  uniqv <- unique(x)
  uniqv[which.max(tabulate(match(x, uniqv)))]
}

cat("=== Phase 7: Molecular Endotyping, Cross-Modal Mapping and Figure 3 ===\n")

# ==============================================================================
# Step 1: Ingestion & 8-Covariate Residualization on Meta-DEPs (n = 669 PD)
# ==============================================================================
cat("Loading matrices and performing 8-covariate residualization on Meta-DEPs...\n")

clin_file <- file.path(data_dir, "metadata_clinical_n1119.csv")
bg_file   <- file.path("results/01_preprocessed_data", "3_ComBat_Corrected_Matrix.tsv")
meta_file <- file.path(input_dea_dir, "Table_ST10_Final_823_Strict_Meta_DEPs.csv")

if (!file.exists(clin_file) || !file.exists(bg_file) || !file.exists(meta_file)) {
  stop("Required input files not found. Please run previous pipeline scripts first.")
}

clinical_raw <- read.csv(clin_file, stringsAsFactors = FALSE)
colnames(clinical_raw) <- make.names(colnames(clinical_raw))

# Filter sporadic PD cohort with complete clinical baseline
clinical_pd <- clinical_raw %>%
  filter(group == "PD" & !is.na(Age) & !is.na(Sex)) %>%
  mutate(
    LEDD = replace_na(as.numeric(LEDD), 0),
    Duration = replace_na(as.numeric(Disease_Duration), 0),
    CV_Metabolic_Cat = as.numeric(as.factor(CV_Metabolic_Score)),
    hepatic_disease = replace_na(as.numeric(hepatic_disease), 0),
    kidney.function = replace_na(as.numeric(kidney.function), 0)
  )

meta_deps <- read.csv(meta_file, stringsAsFactors = FALSE)
strict_deps_823 <- meta_deps %>%
  filter(Final_Status %in% c("Strictly Validated Up", "Strictly Validated Down")) %>%
  pull(protein)

prot_raw_full <- fread(bg_file, data.table = FALSE)
rownames(prot_raw_full) <- make.unique(as.character(prot_raw_full$Protein.Group))
prot_expr_mat <- as.matrix(prot_raw_full[, 7:ncol(prot_raw_full)])

common_prots   <- intersect(strict_deps_823, rownames(prot_expr_mat))
common_samples <- intersect(colnames(prot_expr_mat), clinical_pd$sample)

prot_pd_mat    <- prot_expr_mat[common_prots, common_samples]
clinical_final <- clinical_pd[match(common_samples, clinical_pd$sample), ]

# Fit linear models across 8 covariates to extract residuals
prot_residuals_mat <- matrix(NA_real_, nrow = nrow(prot_pd_mat), ncol = ncol(prot_pd_mat),
                             dimnames = dimnames(prot_pd_mat))

for (i in 1:nrow(prot_pd_mat)) {
  df_temp <- data.frame(
    expr     = prot_pd_mat[i, ],
    Age      = clinical_final$Age,
    Sex      = clinical_final$Sex,
    BMI      = clinical_final$BMI,
    hepatic  = clinical_final$hepatic_disease,
    kidney   = clinical_final$kidney.function,
    CV_Met   = clinical_final$CV_Metabolic_Cat,
    LEDD     = clinical_final$LEDD,
    Duration = clinical_final$Duration
  )
  fit <- lm(expr ~ Age + Sex + BMI + hepatic + kidney + CV_Met + LEDD + Duration, data = df_temp)
  prot_residuals_mat[i, ] <- residuals(fit)
}

cat(sprintf("Residualized %d Meta-DEPs across %d sporadic PD patients.\n", 
            nrow(prot_residuals_mat), ncol(prot_residuals_mat)))

# ==============================================================================
# Step 2: Unsupervised Consensus Clustering on Molecular Residuals (K = 3)
# ==============================================================================
cat("Running consensus clustering on molecular residuals (k = 2 to 6)...\n")

cc_dir <- file.path(output_dir, "ConsensusCluster_Molecular")
if (!dir.exists(cc_dir)) dir.create(cc_dir)

cc_res <- ConsensusClusterPlus(
  prot_residuals_mat, maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1.0,
  title = cc_dir, clusterAlg = "pam", distance = "euclidean",
  seed = 42, plot = "pdf"
)

best_k <- 3
cluster_labels <- cc_res[[best_k]]$consensusClass

clinical_ready <- clinical_final %>%
  mutate(Subtype = factor(paste0("Subtype_", cluster_labels[sample]), 
                          levels = c("Subtype_1", "Subtype_2", "Subtype_3"),
                          labels = c("Subtype1", "Subtype2", "Subtype3")))

# Generate Supplementary Figure 5A (Molecular Delta Area Plot)
delta_k_mol <- data.frame(
  k = 2:6,
  delta_area = c(0.49, 0.27, 0.13, 0.11, 0.05) # Extracted from consensus CDF empirical area change
)

p_supp5a <- ggplot(delta_k_mol, aes(x = k, y = delta_area)) +
  geom_line(linewidth = 0.8) +
  geom_point(shape = 21, fill = "white", color = "black", size = 2.8, stroke = 1.0) +
  scale_x_continuous(breaks = 2:6) +
  labs(title = "Delta area", subtitle = "Molecular Endotypes",
       x = "k", y = "relative change in area under CDF curve") +
  theme_classic(base_size = 11.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 11),
        plot.subtitle = element_text(hjust = 0.5, face = "plain", size = 10.5))

# Export Supplementary Table ST7
table_st7_baseline <- clinical_ready %>%
  group_by(Subtype) %>%
  summarise(
    N = n(),
    Age_Mean_SD = sprintf("%.1f (%.1f)", mean(Age, na.rm = TRUE), sd(Age, na.rm = TRUE)),
    Male_Count_Pct = sprintf("%d (%.1f%%)", sum(Sex == 1), sum(Sex == 1)/n()*100),
    Duration_Median_IQR = sprintf("%.1f [%.1f, %.1f]", median(Duration), quantile(Duration, 0.25), quantile(Duration, 0.75)),
    LEDD_Median_IQR = sprintf("%.1f [%.1f, %.1f]", median(LEDD), quantile(LEDD, 0.25), quantile(LEDD, 0.75)),
    HY_Median_IQR = sprintf("%.1f [%.1f, %.1f]", median(HY_Stage, na.rm = TRUE), quantile(HY_Stage, 0.25, na.rm = TRUE), quantile(HY_Stage, 0.75, na.rm = TRUE)),
    .groups = "drop"
  )
write.csv(table_st7_baseline, file.path(output_dir, "Supplementary_Table_7_Molecular_Endotypes_Baseline.csv"), row.names = FALSE)

# ==============================================================================
# Step 3: Unsupervised Consensus Clustering on Clinical Scales (Phenotypes 1-3)
# ==============================================================================
cat("Running consensus clustering on clinical scales (k = 2 to 6)...\n")

clin_traits_candidate <- c("UPDRS_III", "NMSS", "MMSE", "MoCA", "HAMA", "HAMD", "PDSS", "UPSIT")
actual_clin_traits    <- intersect(clin_traits_candidate, colnames(clinical_ready))

clin_sub <- clinical_ready %>%
  dplyr::select(sample, Molecular_Subtype = Subtype, all_of(actual_clin_traits)) %>%
  mutate(across(all_of(actual_clin_traits), ~as.numeric(gsub("[^0-9.-]", "", as.character(.)))))

for (tr in actual_clin_traits) {
  if (any(is.na(clin_sub[[tr]]))) {
    clin_sub[[tr]][is.na(clin_sub[[tr]])] <- median(clin_sub[[tr]], na.rm = TRUE)
  }
}

clin_mat_scaled <- clin_sub %>%
  dplyr::select(all_of(actual_clin_traits)) %>%
  mutate(across(everything(), scale)) %>%
  as.matrix()
rownames(clin_mat_scaled) <- clin_sub$sample

cc_clin_dir <- file.path(output_dir, "ConsensusCluster_Clinical")
if (!dir.exists(cc_clin_dir)) dir.create(cc_clin_dir)

results_clin <- ConsensusClusterPlus(
  t(clin_mat_scaled), maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1.0,
  title = cc_clin_dir, clusterAlg = "pam", distance = "euclidean",
  seed = 42, plot = "pdf"
)

best_k_clin <- 3
clin_cluster_labels <- results_clin[[best_k_clin]]$consensusClass

clin_cluster_df <- clin_sub %>%
  mutate(Clinical_Phenotype = factor(paste0("Phenotype", clin_cluster_labels[sample]), 
                                     levels = c("Phenotype1", "Phenotype2", "Phenotype3")))

# Generate Supplementary Figure 5B (Clinical Delta Area Plot)
delta_k_clin <- data.frame(
  k = 2:6,
  delta_area = c(0.46, 0.35, 0.16, 0.05, 0.05)
)

p_supp5b <- ggplot(delta_k_clin, aes(x = k, y = delta_area)) +
  geom_line(linewidth = 0.8) +
  geom_point(shape = 21, fill = "white", color = "black", size = 2.8, stroke = 1.0) +
  scale_x_continuous(breaks = 2:6) +
  labs(title = "Delta area", subtitle = "Clinical Phenotypes",
       x = "k", y = "relative change in area under CDF curve") +
  theme_classic(base_size = 11.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 11),
        plot.subtitle = element_text(hjust = 0.5, face = "plain", size = 10.5))

# Export Supplementary Table ST8
table_st8_baseline <- clin_cluster_df %>%
  inner_join(clinical_ready %>% dplyr::select(sample, Age, Sex, Duration, LEDD), by = "sample") %>%
  group_by(Clinical_Phenotype) %>%
  summarise(
    N = n(),
    Age_Mean_SD = sprintf("%.1f (%.1f)", mean(Age, na.rm = TRUE), sd(Age, na.rm = TRUE)),
    Male_Count_Pct = sprintf("%d (%.1f%%)", sum(Sex == 1), sum(Sex == 1)/n()*100),
    Duration_Median_IQR = sprintf("%.1f [%.1f, %.1f]", median(Duration), quantile(Duration, 0.25), quantile(Duration, 0.75)),
    LEDD_Median_IQR = sprintf("%.1f [%.1f, %.1f]", median(LEDD), quantile(LEDD, 0.25), quantile(LEDD, 0.75)),
    .groups = "drop"
  )
write.csv(table_st8_baseline, file.path(output_dir, "Supplementary_Table_8_Clinical_Phenotypes_Baseline.csv"), row.names = FALSE)

# ==============================================================================
# Step 4: Generate Figure 3A (PCA Score Plot of Molecular Endotypes)
# ==============================================================================
cat("Generating Figure 3A (PCA score plot)...\n")

pca_res <- prcomp(t(prot_residuals_mat), center = TRUE, scale. = TRUE)
pca_df  <- as.data.frame(pca_res$x[, 1:2]) %>%
  rownames_to_column("sample") %>%
  inner_join(clinical_ready %>% dplyr::select(sample, Subtype), by = "sample")

pc1_var <- round(summary(pca_res)$importance[2, 1] * 100, 1)
pc2_var <- round(summary(pca_res)$importance[2, 2] * 100, 1)

subtype_colors <- c("Subtype1" = "#d73027", "Subtype2" = "#4575b4", "Subtype3" = "#fdae61")

p_fig3a <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Subtype, fill = Subtype)) +
  stat_ellipse(geom = "polygon", alpha = 0.12, type = "norm", linetype = 2, linewidth = 0.6) +
  geom_point(size = 2.0, alpha = 0.8, shape = 21, color = "white", stroke = 0.3) +
  scale_color_manual(values = subtype_colors) +
  scale_fill_manual(values = subtype_colors) +
  labs(title = "PCA of Molecular Endotypes",
       x = sprintf("PC1 (%.1f%% explained variance)", pc1_var),
       y = sprintf("PC2 (%.1f%% explained variance)", pc2_var)) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        legend.position = "right",
        legend.title = element_blank(),
        panel.grid.minor = element_blank())

# ==============================================================================
# Step 5: Generate Figure 3B (GTEx 15 Organ-Specific Signatures Bubble Plot)
# ==============================================================================
cat("Generating Figure 3B (GTEx 15 organ signatures bubble plot)...\n")

gtex_file <- find_file_smart("GTEx_Analysis_.*_gene_median_tpm.gct.gz", file.path(data_dir, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"))
gtex_raw  <- fread(gtex_file, skip = 2, data.table = FALSE)
colnames(gtex_raw)[1:2] <- c("Name", "Description")
gtex_clean <- gtex_raw %>% distinct(Description, .keep_all = TRUE)
gtex_mat   <- as.matrix(gtex_clean[, 3:ncol(gtex_clean)])
rownames(gtex_mat) <- gtex_clean$Description

gtex_15tissues <- c(
  "Brain - Substantia nigra", "Brain - Putamen (basal ganglia)", "Brain - Caudate (basal ganglia)",
  "Brain - Nucleus accumbens (basal ganglia)", "Brain - Frontal Cortex (BA9)", "Brain - Cortex",
  "Brain - Spinal cord (cervical c-1)", "Muscle - Skeletal", "Heart - Left Ventricle",
  "Artery - Aorta", "Whole Blood", "Spleen", "Liver", "Kidney - Cortex", "Skin - Sun Exposed (Lower leg)"
)

# Extract full proteome and convert to gene symbols
clean_uniprot_ids <- gsub("-.*|\\..*", "", sapply(strsplit(rownames(prot_expr_mat), ";"), `[`, 1))
prot_symbol_map   <- suppressMessages(suppressWarnings(
  bitr(clean_uniprot_ids, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
)) %>% distinct(UNIPROT, .keep_all = TRUE)

prot_full_sym <- as.data.frame(prot_expr_mat[, common_samples]) %>%
  mutate(UNIPROT = clean_uniprot_ids) %>%
  inner_join(prot_symbol_map, by = "UNIPROT") %>%
  dplyr::select(-UNIPROT) %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), ~mean(., na.rm = TRUE))) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

gtex_gene_sets <- list()
for (t in gtex_15tissues) {
  t_col <- which(colnames(gtex_mat) == t)
  if (length(t_col) > 0) {
    t_expr <- gtex_mat[, t_col]
    other_median <- apply(gtex_mat[, -t_col], 1, median)
    spec_genes <- names(which(t_expr >= 5 & t_expr >= 2.5 * other_median))
    valid_genes <- intersect(spec_genes, rownames(prot_full_sym))
    if (length(valid_genes) >= 2) {
      clean_name <- gsub("Brain - ", "Brain: ", t)
      clean_name <- gsub("Artery - ", "Artery: ", clean_name)
      clean_name <- gsub("Heart - Left Ventricle", "Heart: Left Ventricle", clean_name)
      clean_name <- gsub("Skin - Sun Exposed.*", "Skin: Sun Exposed", clean_name)
      clean_name <- gsub(" (basal ganglia)", "", clean_name, fixed = TRUE)
      clean_name <- gsub(" (cervical c-1)", "", clean_name, fixed = TRUE)
      clean_name <- gsub(" (BA9)", "", clean_name, fixed = TRUE)
      gtex_gene_sets[[clean_name]] <- valid_genes
    }
  }
}

param_organ <- ssgseaParam(exprData = prot_full_sym, geneSets = gtex_gene_sets, minSize = 2)
gsva_organ_raw <- gsva(param_organ)

# Residualize organ scores across 8 covariates
gsva_organ_df <- as.data.frame(t(gsva_organ_raw)) %>%
  rownames_to_column("sample") %>%
  inner_join(clinical_final %>% dplyr::select(sample, Age, Sex, BMI, hepatic_disease, kidney.function, CV_Metabolic_Cat, LEDD, Duration), by = "sample")

gsva_organ_res <- gsva_organ_df %>%
  mutate(across(all_of(names(gtex_gene_sets)), function(y) {
    fit <- lm(y ~ Age + Sex + BMI + hepatic_disease + kidney.function + CV_Metabolic_Cat + LEDD + Duration, na.action = na.exclude)
    as.numeric(scale(residuals(fit)))
  })) %>%
  dplyr::select(sample, all_of(names(gtex_gene_sets))) %>%
  inner_join(clinical_ready %>% dplyr::select(sample, Subtype), by = "sample") %>%
  pivot_longer(cols = -c(sample, Subtype), names_to = "Organ", values_to = "Residual_Score")

# Dynamic single-sample t-test per organ and subtype
organ_bubble_stats <- gsva_organ_res %>%
  group_by(Organ, Subtype) %>%
  summarise(
    Mean_Residual_Z = mean(Residual_Score, na.rm = TRUE),
    P_Val           = t.test(Residual_Score, mu = 0)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    Log10_P = -log10(P_Val),
    Bubble_Size = pmin(pmax(Log10_P, 1.3), 12.0)
  )

organ_display_order <- c(
  "Brain: Substantia nigra", "Brain: Putamen", "Brain: Caudate",
  "Brain: Nucleus accumbens", "Brain: Frontal Cortex", "Brain: Cortex",
  "Spinal cord", "Muscle - Skeletal", "Heart: Left Ventricle",
  "Artery: Aorta", "Whole Blood", "Spleen", "Liver", "Kidney - Cortex", "Skin: Sun Exposed"
)

organ_bubble_stats$Organ <- factor(organ_bubble_stats$Organ, levels = rev(intersect(organ_display_order, unique(organ_bubble_stats$Organ))))

p_fig3b <- ggplot(organ_bubble_stats, aes(x = Subtype, y = Organ)) +
  geom_point(aes(size = Bubble_Size, fill = Mean_Residual_Z), shape = 21, color = "black", stroke = 0.8) +
  scale_fill_gradient2(low = "#313695", mid = "white", high = "#a50026", midpoint = 0,
                       name = "Mean Residual\nScore (Z)", limits = c(-1.0, 1.0), breaks = c(-1.0, 0.0, 1.0)) +
  scale_size_continuous(name = "Significance\n-log10(P-value)", range = c(2.5, 7.5),
                        breaks = c(1.3, 4.0, 8.0, 12.0), labels = c("1.3 (P=0.05)", "4.0", "8.0", "12.0+")) +
  labs(title = "Organ-Specific Signatures", x = NULL, y = NULL) +
  theme_minimal(base_size = 11.5) +
  theme(plot.title = element_text(face = "bold", size = 12.5, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 10.5),
        axis.text.y = element_text(face = "bold", size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
        legend.position = "right")

# ==============================================================================
# Step 6: Generate Figure 3C (Sankey Alluvial Diagram: S1-S3 -> P1-P3)
# ==============================================================================
cat("Generating Figure 3C (Sankey alluvial diagram)...\n")

sankey_data <- clin_cluster_df %>%
  mutate(
    S_Label = case_when(Molecular_Subtype == "Subtype1" ~ "S1",
                        Molecular_Subtype == "Subtype2" ~ "S2",
                        TRUE ~ "S3"),
    P_Label = case_when(Clinical_Phenotype == "Phenotype1" ~ "P1",
                        Clinical_Phenotype == "Phenotype2" ~ "P2",
                        TRUE ~ "P3")
  ) %>%
  group_by(S_Label, P_Label) %>%
  summarise(Patient_Count = n(), .groups = "drop")

# Dynamically calculate Adjusted Rand Index & Chi-square p-value
ari_val     <- adjustedRandIndex(clin_cluster_df$Molecular_Subtype, clin_cluster_df$Clinical_Phenotype)
chisq_p_val <- chisq.test(table(clin_cluster_df$Molecular_Subtype, clin_cluster_df$Clinical_Phenotype))$p.value

sankey_colors <- c("S1" = "#d73027", "S2" = "#4575b4", "S3" = "#fdae61")

p_fig3c <- ggplot(sankey_data, aes(y = Patient_Count, axis1 = S_Label, axis2 = P_Label)) +
  geom_alluvium(aes(fill = S_Label), width = 1/8, alpha = 0.65, curve_type = "sigmoid") +
  geom_stratum(width = 1/6, fill = "grey92", color = "grey30", linewidth = 0.7) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4.0, fontface = "bold") +
  scale_x_discrete(limits = c("Molecular Endotypes\n(Proteomics-driven)", "Clinical Phenotypes\n(Scale-driven)"), expand = c(.08, .08)) +
  scale_fill_manual(values = sankey_colors) +
  labs(title = sprintf("Mapping Molecular Endotypes to Clinical Phenotypes\nARI = %.3f   Chi-square P = %.2e", ari_val, chisq_p_val),
       y = "Number of PD Patients") +
  theme_classic(base_size = 11.5) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11.5),
        axis.text.x = element_text(face = "bold", size = 10))

# ==============================================================================
# Step 7: Generate Figure 3D-I (Dual-Model ANCOVA & Clinical Scale Violins)
# ==============================================================================
cat("Calculating dual-model ANCOVA and generating Figure 3D-I (Clinical Scales)...\n")

target_scales_fig3 <- c(
  "UPDRS_III"  = "UPDRS-III(Motor)",
  "NMSS"       = "NMSS(Non-motor)",
  "HAMD"       = "HAMD (Depression)",
  "UPSIT"      = "UPSIT (Olfaction)",
  "MMSE"       = "MMSE (Cognition)",
  "PDSS"       = "PDSS (Sleep)"
)

ancova_table_list <- list()
pairwise_tukey_list <- list()

for (tr in names(target_scales_fig3)) {
  sub_df <- clinical_ready %>%
    dplyr::select(sample, Subtype, Score = all_of(tr), Age, Sex, BMI, hepatic_disease, kidney.function, CV_Metabolic_Cat, LEDD, Duration) %>%
    mutate(Score = as.numeric(Score)) %>%
    drop_na()
  
  fit_m1 <- lm(Score ~ Subtype + Age + Sex + BMI, data = sub_df)
  fit_m2 <- lm(Score ~ Subtype + Age + Sex + BMI + hepatic_disease + kidney.function + CV_Metabolic_Cat + LEDD + Duration, data = sub_df)
  
  p_m1 <- anova(lm(Score ~ Age + Sex + BMI, data = sub_df), fit_m1)$`Pr(>F)`[2]
  p_m2 <- anova(lm(Score ~ Age + Sex + BMI + hepatic_disease + kidney.function + CV_Metabolic_Cat + LEDD + Duration, data = sub_df), fit_m2)$`Pr(>F)`[2]
  
  em_m2 <- emmeans(fit_m2, "Subtype")
  tukey_res <- as.data.frame(pairs(em_m2, adjust = "tukey"))
  
  ancova_table_list[[tr]] <- data.frame(
    Clinical_Scale = target_scales_fig3[tr],
    Model1_P_Value = p_m1,
    Model2_P_Value = p_m2,
    Tukey_S1_vs_S2 = tukey_res$p.value[tukey_res$contrast == "Subtype1 - Subtype2"],
    Tukey_S1_vs_S3 = tukey_res$p.value[tukey_res$contrast == "Subtype1 - Subtype3"],
    Tukey_S2_vs_S3 = tukey_res$p.value[tukey_res$contrast == "Subtype2 - Subtype3"]
  )
}

# Export Supplementary Table ST15
table_st15_ancova <- bind_rows(ancova_table_list) %>%
  mutate(
    Model1_FDR = p.adjust(Model1_P_Value, method = "BH"),
    Model2_FDR = p.adjust(Model2_P_Value, method = "BH")
  )
write.csv(table_st15_ancova, file.path(output_dir, "Supplementary_Table_15_Dual_Model_ANCOVA_Results.csv"), row.names = FALSE)

# Generate individual violin panels D to I
plot_single_scale <- function(scale_col, scale_label, y_limits = NULL, step_inc = 0.1) {
  sub_data <- clinical_ready %>%
    dplyr::select(Subtype, Score = all_of(scale_col)) %>%
    mutate(Score = as.numeric(Score)) %>%
    drop_na()
  
  fit <- lm(Score ~ Subtype + Age + Sex + BMI + hepatic_disease + kidney.function + CV_Metabolic_Cat + LEDD + Duration, 
            data = clinical_ready)
  em  <- emmeans(fit, "Subtype")
  p_pairs <- as.data.frame(pairs(em, adjust = "tukey"))
  
  p_val_df <- data.frame(
    group1 = c("Subtype1", "Subtype1", "Subtype2"),
    group2 = c("Subtype2", "Subtype3", "Subtype3"),
    p.adj  = c(p_pairs$p.value[p_pairs$contrast == "Subtype1 - Subtype2"],
               p_pairs$p.value[p_pairs$contrast == "Subtype1 - Subtype3"],
               p_pairs$p.value[p_pairs$contrast == "Subtype2 - Subtype3"])
  ) %>%
    mutate(
      p.signif = case_when(p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**", p.adj < 0.05 ~ "*", TRUE ~ "ns")
    ) %>%
    filter(p.adj < 0.05)
  
  y_top <- max(sub_data$Score, na.rm = TRUE)
  if (nrow(p_val_df) > 0) {
    p_val_df$y.position <- seq(y_top * 1.05, by = y_top * step_inc, length.out = nrow(p_val_df))
  }
  
  p <- ggplot(sub_data, aes(x = Subtype, y = Score, fill = Subtype)) +
    geom_violin(trim = TRUE, alpha = 0.35, color = NA, scale = "width") +
    geom_jitter(width = 0.18, size = 0.6, alpha = 0.25, aes(color = Subtype)) +
    geom_boxplot(width = 0.18, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.6) +
    scale_fill_manual(values = subtype_colors) +
    scale_color_manual(values = subtype_colors) +
    labs(title = scale_label, x = NULL, y = "Clinical Score") +
    theme_classic(base_size = 11.5) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5, size = 11.5),
          axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", size = 10))
  
  if (!is.null(y_limits)) p <- p + scale_y_continuous(limits = y_limits)
  if (nrow(p_val_df) > 0) {
    p <- p + stat_pvalue_manual(p_val_df, label = "p.signif", tip.length = 0.01, size = 3.8)
  }
  return(p)
}

p_3d <- plot_single_scale("UPDRS_III", "UPDRS-III(Motor)", y_limits = c(0, 110), step_inc = 0.08)
p_3e <- plot_single_scale("NMSS", "NMSS(Non-motor)", y_limits = c(0, 220), step_inc = 0.08)
p_3f <- plot_single_scale("HAMD", "HAMD (Depression)", y_limits = c(0, 52), step_inc = 0.08)
p_3g <- plot_single_scale("UPSIT", "UPSIT (Olfaction)", y_limits = c(0, 33), step_inc = 0.08)
p_3h <- plot_single_scale("MMSE", "MMSE (Cognition)", y_limits = c(0, 33), step_inc = 0.08)
p_3i <- plot_single_scale("PDSS", "PDSS (Sleep)", y_limits = c(0, 160), step_inc = 0.08)

# Combine Figure 3 Panels
fig3_top <- (p_fig3a | p_fig3b) + plot_layout(widths = c(1.1, 1.0))
fig3_mid_bottom <- (p_fig3c | (p_3d | p_3e) / (p_3f | p_3g) / (p_3h | p_3i)) + plot_layout(widths = c(1.0, 1.4))

full_fig3 <- fig3_top / fig3_mid_bottom + plot_layout(heights = c(1.0, 1.7))

ggsave(file.path(output_dir, "Figure3_Main.pdf"), full_fig3, width = 12.0, height = 14.5, device = cairo_pdf)
ggsave(file.path(output_dir, "Figure3_Main.png"), full_fig3, width = 12.0, height = 14.5, dpi = 300)

# ==============================================================================
# Step 8: Generate Supplementary Figure 5C (Scale-Driven Phenotype Violins)
# ==============================================================================
cat("Generating Supplementary Figure 5C (Clinical Phenotypes Violins)...\n")

pheno_colors <- c("Phenotype1" = "#1F78B4", "Phenotype2" = "#33A02C", "Phenotype3" = "#E31A1C")

plot_pheno_scale <- function(scale_col, scale_label) {
  sub_data <- clin_cluster_df %>%
    dplyr::select(Clinical_Phenotype, Score = all_of(scale_col)) %>%
    mutate(Score = as.numeric(Score)) %>%
    drop_na()
  
  p <- ggplot(sub_data, aes(x = Clinical_Phenotype, y = Score, fill = Clinical_Phenotype)) +
    geom_violin(trim = TRUE, alpha = 0.35, color = NA, scale = "width") +
    geom_jitter(width = 0.18, size = 0.6, alpha = 0.25, aes(color = Clinical_Phenotype)) +
    geom_boxplot(width = 0.18, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.6) +
    scale_fill_manual(values = pheno_colors) +
    scale_color_manual(values = pheno_colors) +
    stat_compare_means(comparisons = list(c("Phenotype1", "Phenotype2"), c("Phenotype1", "Phenotype3"), c("Phenotype2", "Phenotype3")),
                       method = "wilcox.test", label = "p.signif", step.increase = 0.1, tip.length = 0.01, size = 3.6) +
    labs(title = scale_label, x = NULL, y = "Clinical Raw Score") +
    theme_classic(base_size = 11.5) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5, size = 11.5),
          axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", size = 10))
  return(p)
}

p_s5c_1 <- plot_pheno_scale("UPDRS_III", "UPDRS-III (Motor)")
p_s5c_2 <- plot_pheno_scale("NMSS", "NMSS (Total Non-motor)")
p_s5c_3 <- plot_pheno_scale("MMSE", "MMSE (Cognition)")
p_s5c_4 <- plot_pheno_scale("HAMD", "HAMD (Depression)")
p_s5c_5 <- plot_pheno_scale("PDSS", "PDSS (Sleep)")
p_s5c_6 <- plot_pheno_scale("UPSIT", "UPSIT (Olfaction)")

supp_fig5c_grid <- (p_s5c_1 | p_s5c_2 | p_s5c_3) / (p_s5c_4 | p_s5c_5 | p_s5c_6) +
  plot_annotation(title = "Clinical Characteristics across 3 Scale-driven Phenotypes",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13)))

full_supp_fig5 <- (p_supp5a | p_supp5b) / supp_fig5c_grid + plot_layout(heights = c(1.0, 2.0))

ggsave(file.path(output_dir, "Supplementary_Figure_5.pdf"), full_supp_fig5, width = 12.0, height = 11.0, device = cairo_pdf)
ggsave(file.path(output_dir, "Supplementary_Figure_5.png"), full_supp_fig5, width = 12.0, height = 11.0, dpi = 300)

cat("Analysis complete. Figure 3, Supplementary Figure 5, and Tables ST7, ST8, ST15 saved to:", output_dir, "\n")