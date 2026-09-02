# ==============================================================================
# Script: 01_preprocessing_imputation_split_combat.R
# Purpose: Data preprocessing, completeness filtering, MNAR imputation, and batch correction
# Project: Large-scale plasma proteomics in Parkinson's disease
# ==============================================================================

options(expressions = 5000)
set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
  library(sva)
  library(patchwork)
})

# Directory setup
data_dir   <- "data"
output_dir <- "results/01_preprocessed_data"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("=== Phase 1: Proteomics Data Preprocessing and Harmonization ===\n")

# ==============================================================================
# Step 1: Load Clinical Metadata and Raw DIA-MS Quantification Matrix
# ==============================================================================
cat("Loading metadata and raw quantification matrix...\n")

clin_raw_file <- file.path(data_dir, "metadata_clinical_n1119.csv")
batch_file    <- file.path(data_dir, "batch_info.csv")
raw_pg_file   <- file.path(data_dir, "raw_dia_protein_matrix.tsv")

if (!file.exists(raw_pg_file)) {
  stop(sprintf("Input matrix file not found: %s", raw_pg_file))
}

clin_raw_data <- read.csv(clin_raw_file, stringsAsFactors = FALSE)
batch_data    <- read.csv(batch_file, stringsAsFactors = FALSE)

# Standardize group identifiers
clin_data_clean <- clin_raw_data %>%
  mutate(group = case_when(
    group %in% c("CTR", "Control", "HC") ~ "HC",
    group == "PD" ~ "PD",
    TRUE ~ as.character(group)
  ))

# Harmonize cohort column
if (!"cohort" %in% colnames(clin_data_clean)) {
  if ("center" %in% colnames(clin_data_clean)) {
    clin_data_clean$cohort <- ifelse(clin_data_clean$center %in% c("Discovery", "Center_1", "Center1"), "Discovery", "Validation")
  } else if ("hospital" %in% colnames(clin_data_clean)) {
    first_site <- unique(clin_data_clean$hospital)[1]
    clin_data_clean$cohort <- ifelse(clin_data_clean$hospital == first_site, "Discovery", "Validation")
  }
}

metadata <- clin_data_clean %>%
  select(sample, group, cohort) %>%
  left_join(batch_data %>% select(sample, board), by = "sample") %>%
  filter(!is.na(board))

raw_matrix   <- fread(raw_pg_file, data.table = FALSE, check.names = FALSE)
protein_anno <- raw_matrix[, 1:6]
expr_data    <- raw_matrix[, 7:ncol(raw_matrix)]

common_samples <- intersect(colnames(expr_data), metadata$sample)
expr_data      <- expr_data[, common_samples]
metadata       <- metadata[match(common_samples, metadata$sample), ]

expr_data_mat <- as.matrix(expr_data)
mode(expr_data_mat) <- "numeric"

# ==============================================================================
# Step 2: Log2 Transformation
# ==============================================================================
cat("Applying log2 transformation...\n")
expr_data_mat[expr_data_mat <= 0] <- NA
expr_log2 <- log2(expr_data_mat)

# ==============================================================================
# Step 3: Feature Completeness Filtering (30% Threshold)
# ==============================================================================
cat("Filtering features based on 30% completeness threshold...\n")
samples_PD <- metadata$sample[metadata$group == "PD"]
samples_HC <- metadata$sample[metadata$group == "HC"]

valid_rate_PD <- rowSums(!is.na(expr_log2[, samples_PD])) / length(samples_PD)
valid_rate_HC <- rowSums(!is.na(expr_log2[, samples_HC])) / length(samples_HC)

# Retain proteins with at least 30% completeness in either group
keep_idx <- which((valid_rate_PD >= 0.3) | (valid_rate_HC >= 0.3))
expr_filtered <- expr_log2[keep_idx, ]
protein_anno_filtered <- protein_anno[keep_idx, ]

filtered_out <- cbind(protein_anno_filtered, expr_filtered)
fwrite(filtered_out, file.path(output_dir, "1_Filtered_Log2_Matrix.tsv"), sep = "\t", quote = FALSE)
cat(sprintf("Retained %d proteins across %d samples after filtering.\n", nrow(expr_filtered), ncol(expr_filtered)))

# ==============================================================================
# Step 4: Down-shifted Normal Imputation (MNAR)
# ==============================================================================
cat("Performing down-shifted normal imputation (MNAR)...\n")
expr_imputed <- expr_filtered
for (i in 1:ncol(expr_imputed)) {
  col_data   <- expr_imputed[, i]
  valid_data <- col_data[!is.na(col_data)]
  n_missing  <- sum(is.na(col_data))
  
  if (n_missing > 0 && length(valid_data) > 0) {
    col_mean <- mean(valid_data)
    col_std  <- sd(valid_data)
    if (is.na(col_std) || col_std == 0) col_std <- 0.1
    
    impute_mean <- col_mean - (1.8 * col_std)
    impute_std  <- 0.3 * col_std
    
    expr_imputed[is.na(col_data), i] <- rnorm(n_missing, mean = impute_mean, sd = impute_std)
  }
}
imputed_out <- cbind(protein_anno_filtered, expr_imputed)
fwrite(imputed_out, file.path(output_dir, "2_Imputed_Matrix.tsv"), sep = "\t", quote = FALSE)
cat("Imputed matrix generated.\n")

# ==============================================================================
# Step 5: Independent Split-ComBat Batch Correction
# ==============================================================================
cat("Running Split-ComBat batch correction per cohort...\n")

expr_combat_all <- expr_imputed 

plot_pca <- function(mat, meta, title) {
  pca_res <- prcomp(t(mat), scale. = TRUE)
  pca_df  <- as.data.frame(pca_res$x[, 1:2])
  pca_df$sample <- rownames(pca_df)
  pca_df  <- left_join(pca_df, meta, by = "sample")
  var_explained <- summary(pca_res)$importance[2, 1:2] * 100
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = as.character(board), shape = group)) +
    geom_point(size = 2.8, alpha = 0.8) +
    theme_bw(base_size = 11) +
    labs(title = title,
         x = sprintf("PC1 (%.1f%%)", var_explained[1]),
         y = sprintf("PC2 (%.1f%%)", var_explained[2]),
         color = "MS Plate/Batch", shape = "Diagnosis") +
    theme(plot.title = element_text(face = "bold", size = 12))
}

pca_plots <- list()
for (cohort_name in c("Discovery", "Validation")) {
  cohort_samples <- metadata$sample[metadata$cohort == cohort_name]
  cohort_meta    <- metadata[metadata$cohort == cohort_name, ]
  cohort_expr    <- expr_imputed[, cohort_samples]
  
  pca_plots[[paste0(cohort_name, "_Before")]] <- plot_pca(cohort_expr, cohort_meta, paste(cohort_name, "Cohort - Raw (Before ComBat)"))
  
  mod <- model.matrix(~ as.factor(group), data = cohort_meta)
  if (length(unique(cohort_meta$board)) > 1) {
    cohort_expr_combat <- ComBat(dat = as.matrix(cohort_expr), 
                                 batch = as.character(cohort_meta$board), 
                                 mod = mod, 
                                 par.prior = TRUE)
  } else {
    cohort_expr_combat <- cohort_expr
  }
  pca_plots[[paste0(cohort_name, "_After")]] <- plot_pca(cohort_expr_combat, cohort_meta, paste(cohort_name, "Cohort - Corrected (After ComBat)"))
  expr_combat_all[, cohort_samples] <- cohort_expr_combat
}

combat_out <- cbind(protein_anno_filtered, expr_combat_all)
fwrite(combat_out, file.path(output_dir, "3_ComBat_Corrected_Matrix.tsv"), sep = "\t", quote = FALSE)

# Export Supplementary Figure S2B
final_pca_grid <- (pca_plots[["Discovery_Before"]] | pca_plots[["Discovery_After"]]) /
  (pca_plots[["Validation_Before"]] | pca_plots[["Validation_After"]])

ggsave(file.path(output_dir, "Supplementary_Figure_S2B_ComBat_PCA_Diagnostics.pdf"), final_pca_grid, width = 14, height = 10, device = cairo_pdf)
ggsave(file.path(output_dir, "Supplementary_Figure_S2B_ComBat_PCA_Diagnostics.png"), final_pca_grid, width = 14, height = 10, dpi = 300)

cat("Preprocessing complete. Processed matrices and diagnostics saved to:", output_dir, "\n")