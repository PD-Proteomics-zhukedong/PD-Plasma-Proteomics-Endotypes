# ==============================================================================
# Script: 10_multicohort_endotype_validation.R
# Purpose: Multi-cohort validation of molecular endotypes across PPMI cohorts (Figure 6A-C)
# Project: Large-scale plasma proteomics identifies molecularly distinct PD endotypes
# ==============================================================================

options(expressions = 5000)
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(arrow)
  library(readxl)
  library(randomForest)
  library(pheatmap)
  library(ggpubr)
  library(patchwork)
  library(scales)
  library(org.Hs.eg.db)
})

# Directory setup
data_dir          <- "data"
input_dea_dir     <- "results/02_figure1_meta_dea"
input_clust_dir   <- "results/07_figure3_molecular_endotypes"
input_mech_dir    <- "results/08_subtype_mechanisms_and_deconvolution"
input_ml_dir      <- "results/09_figure5_machine_learning_translation"
output_dir        <- "results/10_figure6_multicohort_validation"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

find_file_smart <- function(pattern, default_path) {
  files <- list.files(data_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  files <- files[!grepl("^~\\$", basename(files))]
  if (length(files) > 0) return(files[1])
  return(default_path)
}

get_col_safe <- function(df, patterns, type = c("numeric", "character"), default = NA) {
  type <- match.arg(type)
  matched_col <- NA_character_
  for (pat in patterns) {
    hit <- grep(pat, colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) { matched_col <- hit[1]; break }
  }
  if (!is.na(matched_col)) {
    raw_vec <- df[[matched_col]]
    if (type == "numeric") return(as.numeric(gsub("[^0-9.-]", "", as.character(raw_vec))))
    else return(as.character(raw_vec))
  } else {
    if (type == "numeric") return(rep(as.numeric(default), nrow(df)))
    else return(rep(as.character(default), nrow(df)))
  }
}

cat("=== Phase 10: Multi-Cohort Endotype Validation and Figure 6 ===\n")

# ==============================================================================
# Step 1: Ingest In-House Training Matrix (n = 669) and 45-Biomarker Classifier
# ==============================================================================
cat("Loading in-house reference cohort and 45-biomarker panel...\n")

clin_file   <- file.path(data_dir, "metadata_clinical_n1119.csv")
combat_file <- file.path("results/01_preprocessed_data", "3_ComBat_Corrected_Matrix.tsv")
st16_file   <- file.path(input_mech_dir, "Supplementary_Table_16_Subtype_Biomarker_Specificity.csv")

clinical_raw <- read.csv(clin_file, stringsAsFactors = FALSE)
colnames(clinical_raw) <- make.names(colnames(clinical_raw))

clinical_pd <- clinical_raw %>%
  filter(group == "PD" & !is.na(Age) & !is.na(Sex)) %>%
  mutate(
    UPDRS3   = get_col_safe(., c("^UPDRS3|^UPDRS.*3", "NP3TOT"), "numeric", NA),
    UPSIT    = get_col_safe(., c("^UPSIT", "PRTOTSCORE"), "numeric", NA)
  )

prot_raw_full <- fread(combat_file, data.table = FALSE)
rownames(prot_raw_full) <- make.unique(as.character(prot_raw_full$Protein.Group))
prot_expr_mat <- as.matrix(prot_raw_full[, 7:ncol(prot_raw_full)])

common_samples <- intersect(colnames(prot_expr_mat), clinical_pd$sample)
prot_pd_mat    <- prot_expr_mat[, common_samples]
in_house_df    <- clinical_pd[match(common_samples, clinical_pd$sample), ]

clean_uniprots <- gsub("-.*|\\..*", "", rownames(prot_pd_mat))
mapped_symbols <- suppressMessages(suppressWarnings(
  mapIds(org.Hs.eg.db, keys = clean_uniprots, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
))
mapped_symbols[is.na(mapped_symbols)] <- clean_uniprots[is.na(mapped_symbols)]
lookup_global <- setNames(as.character(mapped_symbols), rownames(prot_pd_mat))

st16_df <- read.csv(st16_file, stringsAsFactors = FALSE)
selected_45_panel <- intersect(st16_df$UNIPROT, rownames(prot_pd_mat))

X_tr_45 <- t(prot_pd_mat[selected_45_panel, ])
Y_tr    <- factor(in_house_df$Subtype, levels = c("Subtype1", "Subtype2", "Subtype3"))
priors_vec <- table(Y_tr) / length(Y_tr)

# ==============================================================================
# Step 2: Ingest PPMI Clinical Data and Blinded Subtype Prediction across 3 Projects
# ==============================================================================
cat("Ingesting baseline PPMI clinical data and predicting molecular endotypes...\n")

clin_ppmi_path <- find_file_smart("PPMI_Curated_Data_Cut.*\\.xlsx$", file.path(data_dir, "PPMI_Curated_Data_Cut_Public.xlsx"))
if (!file.exists(clin_ppmi_path)) {
  clin_ppmi_path <- find_file_smart("PPMI_Curated_Data_Cut.*\\.csv$", file.path(data_dir, "PPMI_Curated_Data_Cut_Public.csv"))
}
ppmi_raw <- if (grepl("\\.xlsx$", clin_ppmi_path)) readxl::read_excel(clin_ppmi_path) else read.csv(clin_ppmi_path, stringsAsFactors = FALSE)
colnames(ppmi_raw) <- make.names(colnames(ppmi_raw))

ppmi_baseline_pd <- data.frame(
  PATNO      = get_col_safe(ppmi_raw, c("^PATNO$", "^patno$"), "character", "UNKNOWN"),
  EVENT_ID   = toupper(get_col_safe(ppmi_raw, c("^EVENT_ID$", "^VISIT$"), "character", "BL")),
  Raw_Cohort = get_col_safe(ppmi_raw, c("^COHORT$", "^RESEARCH_GROUP$"), "character", NA),
  UPDRS3     = get_col_safe(ppmi_raw, c("NP3TOT", "UPDRS3", "UPDRS_III"), "numeric", NA),
  UPSIT      = get_col_safe(ppmi_raw, c("UPSIT", "PRTOTSCORE"), "numeric", NA),
  stringsAsFactors = FALSE
) %>%
  filter(EVENT_ID %in% c("BL", "SC", "V00") & (Raw_Cohort %in% c("1", 1) | grepl("^Parkinson|^PD$", Raw_Cohort, ignore.case = TRUE))) %>%
  distinct(PATNO, .keep_all = TRUE)

# Ingest PPMI Proteomics Projects
p314_parquet <- find_file_smart("ppmi_proj314.*\\.parquet$", file.path(data_dir, "ppmi_proj314_plasma_screened_extended_npx.parquet"))
p293_parquet <- find_file_smart("ppmi_proj293.*\\.parquet$", file.path(data_dir, "ppmi_proj293_plasma_screened_extended_npx.parquet"))
p9000_files  <- list.files(data_dir, pattern = "PPMI_Project_9000_Plasma.*\\.csv$", full.names = TRUE, recursive = TRUE)

p314_df  <- if (file.exists(p314_parquet)) read_parquet(p314_parquet) else NULL
p293_df  <- if (file.exists(p293_parquet)) read_parquet(p293_parquet) else NULL
p9000_df <- if (length(p9000_files) > 0) bind_rows(lapply(p9000_files, read.csv, stringsAsFactors = FALSE)) else NULL

predict_ppmi_cohort <- function(p_df, proj_name) {
  if (is.null(p_df) || nrow(p_df) == 0) return(NULL)
  u_col <- intersect(c("UNIPROT", "UniProt", "UniprotID"), colnames(p_df))[1]
  a_col <- intersect(c("ASSAY", "Assay", "GeneSymbol"), colnames(p_df))[1]
  n_col <- intersect(c("NPX", "npx", "PCNormalizedNPX", "value"), colnames(p_df))[1]
  p_col <- intersect(c("PATNO", "patno"), colnames(p_df))[1]
  e_col <- intersect(c("EVENT_ID", "event_id", "VISIT"), colnames(p_df))[1]
  
  p_clean <- p_df %>%
    { if (!is.na(e_col)) filter(., toupper(as.character(.data[[e_col]])) %in% c("BL", "SC", "V00", "BASELINE")) else . } %>%
    mutate(PATNO_clean   = gsub("^PPMI-", "", as.character(.data[[p_col]])),
           UniProt_Clean = gsub("-.*|\\..*", "", as.character(.data[[u_col]])),
           Assay_str     = as.character(.data[[a_col]]),
           NPX_num       = as.numeric(.data[[n_col]])) %>%
    filter(PATNO_clean %in% ppmi_baseline_pd$PATNO & !is.na(NPX_num))
  
  panel_symbols <- lookup_global[selected_45_panel]
  p_matched <- p_clean %>%
    filter(UniProt_Clean %in% selected_45_panel | Assay_str %in% panel_symbols) %>%
    mutate(Feature_ID = ifelse(UniProt_Clean %in% selected_45_panel, UniProt_Clean, names(panel_symbols)[match(Assay_str, panel_symbols)]))
  
  wide_mat <- p_matched %>%
    group_by(PATNO_clean, Feature_ID) %>%
    summarise(val = mean(NPX_num, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Feature_ID, values_from = val)
  
  analytical_df <- ppmi_baseline_pd %>% inner_join(wide_mat, by = c("PATNO" = "PATNO_clean"))
  matched_feats <- intersect(selected_45_panel, colnames(analytical_df))
  if (length(matched_feats) < 4 || nrow(analytical_df) < 15) return(NULL)
  
  X_tr_matched <- scale(X_tr_45[, matched_feats, drop = FALSE])
  min_n <- min(table(Y_tr))
  
  set.seed(42)
  rf_model <- randomForest(x = X_tr_matched, y = Y_tr, ntree = 1000, strata = Y_tr, sampsize = rep(min_n, 3))
  
  X_proj_scaled <- scale(as.matrix(analytical_df[, matched_feats, drop = FALSE]))
  X_proj_scaled[is.na(X_proj_scaled)] <- 0
  
  pred_probs <- predict(rf_model, newdata = X_proj_scaled, type = "prob")
  calibrated_probs <- sweep(pred_probs, 2, priors_vec[colnames(pred_probs)], "/")
  pred_class <- colnames(calibrated_probs)[max.col(calibrated_probs)]
  
  analytical_df$Predicted_Subtype <- factor(pred_class, levels = c("Subtype1", "Subtype2", "Subtype3"))
  analytical_df$Project_Name      <- proj_name
  return(analytical_df)
}

p314_res  <- predict_ppmi_cohort(p314_df, "PPMI P314 (Olink Explore)")
p293_res  <- predict_ppmi_cohort(p293_df, "PPMI P293 (Olink Target)")
p9000_res <- predict_ppmi_cohort(p9000_df, "PPMI P9000 (Targeted MS)")

# ==============================================================================
# Step 3: Generate Figure 6A (Subtype Proportion Stability Stacked Bar Chart)
# ==============================================================================
cat("Generating Figure 6A (Subtype proportion stability stacked bar chart)...\n")

in_house_prop <- data.frame(Cohort = "In-House Cohort", Predicted_Subtype = in_house_df$Subtype)

cohort_prop_list <- list(in_house_prop)
if (!is.null(p314_res))  cohort_prop_list[["P314"]]  <- data.frame(Cohort = "PPMI P314 (Olink Explore)", Predicted_Subtype = p314_res$Predicted_Subtype)
if (!is.null(p293_res))  cohort_prop_list[["P293"]]  <- data.frame(Cohort = "PPMI P293 (Olink Target)", Predicted_Subtype = p293_res$Predicted_Subtype)
if (!is.null(p9000_res)) cohort_prop_list[["P9000"]] <- data.frame(Cohort = "PPMI P9000 (Targeted MS)", Predicted_Subtype = p9000_res$Predicted_Subtype)

stacked_prop_df <- bind_rows(cohort_prop_list) %>%
  group_by(Cohort, Predicted_Subtype) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Cohort) %>%
  mutate(
    Total = sum(N),
    Percentage = (N / Total) * 100,
    Cohort = factor(Cohort, levels = c("In-House Cohort", "PPMI P314 (Olink Explore)", "PPMI P293 (Olink Target)", "PPMI P9000 (Targeted MS)"))
  )

subtype_colors <- c("Subtype1" = "#d73027", "Subtype2" = "#4575b4", "Subtype3" = "#fdae61")

p_6a <- ggplot(stacked_prop_df, aes(x = Cohort, y = Percentage, fill = Predicted_Subtype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.55, color = "black", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%\nn=%d", Percentage, N)),
            position = position_stack(vjust = 0.5), size = 3.3, fontface = "bold", color = "white") +
  scale_fill_manual(values = subtype_colors, name = "Endotype") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  labs(title = "Stability of Subtype Proportions Across Independent Cohorts",
       x = NULL, y = "Percentage of PD Patients") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12.5, hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 9.5, color = "black", angle = 18, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    plot.margin = ggplot2::margin(15, 15, 15, 20)
  )

# ==============================================================================
# Step 4: Generate Figure 6B (Cross-Platform Clinical Severity Fingerprint Heatmap)
# ==============================================================================
cat("Generating Figure 6B (Cross-platform clinical severity fingerprint heatmap)...\n")

compute_cohort_z <- function(df, c_name) {
  df %>%
    dplyr::select(Predicted_Subtype, UPDRS3, UPSIT) %>%
    pivot_longer(cols = c("UPDRS3", "UPSIT"), names_to = "Trait", values_to = "Score") %>%
    drop_na(Score) %>%
    group_by(Trait) %>%
    mutate(Score_Z = scale(Score)[, 1]) %>%
    group_by(Predicted_Subtype, Trait) %>%
    summarise(Mean_Z = mean(Score_Z, na.rm = TRUE), .groups = "drop") %>%
    mutate(Cohort = c_name, Column_ID = paste(c_name, Predicted_Subtype, sep = " - "))
}

in_house_formatted <- in_house_df %>%
  mutate(Predicted_Subtype = Subtype) %>%
  dplyr::select(Predicted_Subtype, UPDRS3, UPSIT)

heatmap_data_list <- list()
heatmap_data_list[["InHouse"]] <- compute_cohort_z(in_house_formatted, "In-House Cohort")
if (!is.null(p314_res))  heatmap_data_list[["P314"]]  <- compute_cohort_z(p314_res, "PPMI P314 (Olink Explore)")
if (!is.null(p293_res))  heatmap_data_list[["P293"]]  <- compute_cohort_z(p293_res, "PPMI P293 (Olink Target)")
if (!is.null(p9000_res)) heatmap_data_list[["P9000"]] <- compute_cohort_z(p9000_res, "PPMI P9000 (Targeted MS)")

heatmap_combined <- bind_rows(heatmap_data_list) %>%
  mutate(Trait_Label = ifelse(Trait == "UPDRS3", "Motor Severity\n(UPDRS-III)", "Olfactory Function\n(UPSIT)"))

heatmap_matrix <- heatmap_combined %>%
  dplyr::select(Trait_Label, Column_ID, Mean_Z) %>%
  pivot_wider(names_from = Column_ID, values_from = Mean_Z) %>%
  column_to_rownames("Trait_Label") %>%
  as.matrix()

col_annotation <- data.frame(
  Endotype = gsub(".* - ", "", colnames(heatmap_matrix)),
  Cohort   = gsub(" - Subtype.*", "", colnames(heatmap_matrix)),
  row.names = colnames(heatmap_matrix)
)

annotation_colors <- list(
  Endotype = c("Subtype1" = "#d73027", "Subtype2" = "#4575b4", "Subtype3" = "#fdae61"),
  Cohort   = c("In-House Cohort" = "#2b5c8f", "PPMI P314 (Olink Explore)" = "#81b29a", 
               "PPMI P293 (Olink Target)" = "#7ab6d6", "PPMI P9000 (Targeted MS)" = "#e07a5f")
)

pdf(file.path(output_dir, "Figure6B_Clinical_Severity_Fingerprint_Heatmap.pdf"), width = 11.5, height = 4.8)
pheatmap(
  heatmap_matrix,
  annotation_col    = col_annotation,
  annotation_colors = annotation_colors,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  gaps_col          = c(3, 6, 9),
  breaks            = seq(-0.40, 0.40, length.out = 100),
  color             = colorRampPalette(c("#313695", "#FAF8F8", "#A50026"))(100),
  main              = "Cross-Platform Reproducibility of Baseline Clinical Severity Fingerprints",
  fontsize_row      = 10.5,
  fontsize_col      = 9.0,
  angle_col         = "45",
  border_color      = "grey80",
  name              = "Clinical Z-score"
)
dev.off()

# ==============================================================================
# Step 5: Generate Figure 6C (Phenotypic Severity Gradients PointRange Plot)
# ==============================================================================
cat("Generating Figure 6C (Phenotypic severity gradients pointrange plot)...\n")

format_gradient_df <- function(df, cohort_lbl) {
  if (is.null(df)) return(NULL)
  df %>%
    mutate(Cohort_Label = cohort_lbl) %>%
    dplyr::select(Cohort_Label, Predicted_Subtype, UPDRS3, UPSIT)
}

all_cohorts_gradient <- bind_rows(
  format_gradient_df(in_house_formatted, "In-House Cohort(N=669)"),
  format_gradient_df(p314_res, "PPMI P314 (Olink Explore, N=756)"),
  format_gradient_df(p293_res, "PPMI P293 (Olink Target, N=181)"),
  format_gradient_df(p9000_res, "PPMI P9000 (Targeted MS, N=72)")
) %>%
  mutate(
    Cohort_Label = factor(Cohort_Label, levels = c(
      "In-House Cohort(N=669)", "PPMI P314 (Olink Explore, N=756)", 
      "PPMI P293 (Olink Target, N=181)", "PPMI P9000 (Targeted MS, N=72)"
    ))
  )

gradient_summary <- all_cohorts_gradient %>%
  pivot_longer(cols = c("UPDRS3", "UPSIT"), names_to = "Trait", values_to = "Score") %>%
  drop_na(Score) %>%
  group_by(Cohort_Label, Trait) %>%
  mutate(Score_Z = scale(Score)[, 1]) %>%
  group_by(Cohort_Label, Trait, Predicted_Subtype) %>%
  summarise(
    Mean_Z = mean(Score_Z, na.rm = TRUE),
    SE_Z   = sd(Score_Z, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    Trait_Name = ifelse(Trait == "UPDRS3", "UPDRS-III", "UPSIT")
  )

p_6c <- ggplot(gradient_summary, aes(x = Predicted_Subtype, y = Mean_Z, color = Predicted_Subtype, group = Predicted_Subtype)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_errorbar(aes(ymin = Mean_Z - 1.96 * SE_Z, ymax = Mean_Z + 1.96 * SE_Z), width = 0.20, linewidth = 0.85) +
  geom_point(aes(fill = Predicted_Subtype), shape = 21, size = 3.5, stroke = 1.0, color = "black") +
  facet_grid(Trait_Name ~ Cohort_Label, scales = "free_y") +
  scale_color_manual(values = subtype_colors) +
  scale_fill_manual(values = subtype_colors) +
  labs(
    title = "Cross-Platform Validation of Subtype Phenotypic Gradients across 4 Cohorts",
    x = "Predicted Molecular Endotype",
    y = "Standardized Phenotype Deviation (Z-score)"
  ) +
  theme_bw(base_size = 11.5) +
  theme(
    plot.title    = element_text(face = "bold", size = 12.5, hjust = 0.5),
    strip.background = element_rect(fill = "grey93", color = "black", linewidth = 0.6),
    strip.text.x  = element_text(face = "bold", size = 9.0),
    strip.text.y  = element_text(face = "bold", size = 10.0),
    axis.text.x   = element_text(face = "bold", size = 9.0, color = "black"),
    axis.text.y   = element_text(size = 9.0, color = "black"),
    legend.position = "none",
    plot.margin   = ggplot2::margin(12, 20, 12, 12)
  )

# ==============================================================================
# Step 6: Export Supplementary Table 22 (ST22: Cross-Cohort Gradients)
# ==============================================================================
cat("Calculating cross-cohort statistical parameters for Supplementary Table 22...\n")

calc_gradient_stats <- function(df) {
  stats_rows <- list()
  for (c_lbl in levels(df$Cohort_Label)) {
    for (tr in c("UPDRS3", "UPSIT")) {
      sub <- df %>% filter(Cohort_Label == c_lbl) %>% drop_na(all_of(tr))
      if (nrow(sub) >= 15) {
        kw_p <- kruskal.test(as.formula(paste(tr, "~ Predicted_Subtype")), data = sub)$p.value
        m1 <- sub[[tr]][sub$Predicted_Subtype == "Subtype1"]
        m2 <- sub[[tr]][sub$Predicted_Subtype == "Subtype2"]
        m3 <- sub[[tr]][sub$Predicted_Subtype == "Subtype3"]
        
        p12 <- if (length(m1) >= 3 && length(m2) >= 3) wilcox.test(m1, m2)$p.value else NA
        p23 <- if (length(m2) >= 3 && length(m3) >= 3) wilcox.test(m2, m3)$p.value else NA
        p13 <- if (length(m1) >= 3 && length(m3) >= 3) wilcox.test(m1, m3)$p.value else NA
        
        fmt <- function(x) sprintf("%.2f ± %.2f (%.1f [%.1f–%.1f])", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE), median(x, na.rm=TRUE), quantile(x, 0.25, na.rm=TRUE), quantile(x, 0.75, na.rm=TRUE))
        
        stats_rows[[paste(c_lbl, tr)]] <- data.frame(
          Cohort = c_lbl,
          Phenotypic_Trait = ifelse(tr == "UPDRS3", "Motor Impairment (UPDRS-III)", "Olfactory Deficit (UPSIT)"),
          Evaluated_N = nrow(sub),
          Subtype1_Stats = fmt(m1),
          Subtype2_Stats = fmt(m2),
          Subtype3_Stats = fmt(m3),
          Global_Kruskal_Wallis_P = kw_p,
          Pairwise_P_S1_vs_S2 = p12,
          Pairwise_P_S2_vs_S3 = p23,
          Pairwise_P_S1_vs_S3 = p13
        )
      }
    }
  }
  bind_rows(stats_rows)
}

table_st22_export <- calc_gradient_stats(all_cohorts_gradient)
write.csv(table_st22_export, file.path(output_dir, "Supplementary_Table_22_Cross_Cohort_Subtype_Gradients.csv"), row.names = FALSE)

# ==============================================================================
# Step 7: Export Figure 6 Vector Assets
# ==============================================================================
cat("Exporting Figure 6 panels...\n")

ggsave(file.path(output_dir, "Figure6A_Subtype_Proportions.pdf"), p_6a, width = 7.5, height = 5.8, device = cairo_pdf)
ggsave(file.path(output_dir, "Figure6C_Phenotypic_Gradients.pdf"), p_6c, width = 11.5, height = 6.0, device = cairo_pdf)

cat("Analysis complete. Figure 6 panels and Supplementary Table 22 saved to:", output_dir, "\n")