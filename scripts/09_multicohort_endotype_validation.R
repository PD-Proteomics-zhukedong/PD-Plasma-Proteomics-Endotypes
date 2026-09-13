# ==============================================================================
# Script: 09_multicohort_endotype_validation.R
# Project: Large-scale plasma proteomics in Parkinson's disease (Nature Aging)
# Purpose: Cross-platform out-of-sample validation on PPMI Project 314 (N = 868),
#          generation of Figure 6A (proportions), Figure 6B (Rank-ANCOVA gradients),
#          and Supplementary Tables ST22, ST23, and ST24
# ==============================================================================

options(expressions = 5000)
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(arrow)
  library(readxl)
  library(pROC)
  library(randomForest)
  library(caret)
  library(emmeans)
  library(ggpubr)
  library(cowplot)
  library(patchwork)
  library(scales)
  library(org.Hs.eg.db)
})

# Cross-platform PDF device fallback
pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"

# Directory setup
data_dir        <- "data"
input_dea_dir   <- "results/02_figure1_meta_dea"
input_clust_dir <- "results/06_figure3_molecular_endotypes"
input_mech_dir  <- "results/07_subtype_mechanisms_and_deconvolution"
output_dir      <- "results/09_figure6_multicohort_validation"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

find_file_smart <- function(pattern, default_path) {
  files <- list.files(data_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  files <- files[!grepl("^~\\$", basename(files))]
  if (length(files) > 0) return(files[1])
  return(default_path)
}

get_col_safe <- function(df, patterns, type = c("numeric", "character"), default = NA, exclude = NULL) {
  if (!is.data.frame(df)) return(NA)
  df_cols <- colnames(df)
  for (pat in patterns) {
    hits <- grep(pat, df_cols, ignore.case = TRUE, value = TRUE)
    if (!is.null(exclude)) hits <- hits[!grepl(exclude, hits, ignore.case = TRUE)]
    if (length(hits) > 0) {
      vec <- df[[hits[1]]]
      if (type == "numeric") return(as.numeric(gsub("[^0-9.-]", "", as.character(vec))))
      return(as.character(vec))
    }
  }
  if (type == "numeric") return(rep(as.numeric(default), nrow(df)))
  return(rep(as.character(default), nrow(df)))
}

clean_num <- function(x) {
  if (is.null(x) || all(is.na(x))) return(NA_real_)
  as.numeric(gsub("[^0-9.-]", "", as.character(x)))
}

cat("\n=== Phase 9: Multi-Cohort External Validation on PPMI (Figure 6) ===\n")

# ==============================================================================
# Step 1: Ingest Internal Training Data and 45 Fingerprint Biomarkers
# ==============================================================================
cat("Loading internal reference cohort (n = 670) and 45 fingerprint markers...\n")

subtype_rds_file <- file.path(input_clust_dir, "endotype_assignments.rds")
combat_file      <- file.path("results/01_preprocessed_data", "3_ComBat_Corrected_Matrix.tsv")
st19_file        <- file.path(input_mech_dir, "Supplementary_Table_19_Subtype_45_Fingerprint_Markers.csv")

if (!file.exists(subtype_rds_file) || !file.exists(combat_file) || !file.exists(st19_file)) {
  stop("Required inputs from upstream scripts not found. Please verify execution of Scripts 01, 06, and 07.")
}

subtype_df <- readRDS(subtype_rds_file)
subtype_df$Subtype <- factor(as.character(subtype_df$Subtype), 
                             levels = c("Subtype_1", "Subtype_2", "Subtype_3", "Subtype1", "Subtype2", "Subtype3"),
                             labels = c("Subtype_1", "Subtype_2", "Subtype_3", "Subtype_1", "Subtype_2", "Subtype_3"))

prot_raw_full <- fread(combat_file, data.table = FALSE)
rownames(prot_raw_full) <- make.unique(as.character(prot_raw_full$Protein.Group))
prot_expr_mat <- as.matrix(prot_raw_full[, 7:ncol(prot_raw_full)])

common_samples <- intersect(colnames(prot_expr_mat), subtype_df$sample)
prot_pd_mat    <- prot_expr_mat[, common_samples]
in_house_df    <- subtype_df[match(common_samples, subtype_df$sample), ]

all_uniprots   <- rownames(prot_pd_mat)
clean_uniprots <- sapply(strsplit(as.character(all_uniprots), ";"), `[`, 1) %>% gsub("-.*|\\..*", "", .)
mapped_symbols <- suppressMessages(suppressWarnings(
  mapIds(org.Hs.eg.db, keys = clean_uniprots, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
))
mapped_symbols[is.na(mapped_symbols)] <- clean_uniprots[is.na(mapped_symbols)]
lookup_global <- setNames(as.character(mapped_symbols), all_uniprots)

st19_df <- read.csv(st19_file, stringsAsFactors = FALSE)
selected_subtype_panel <- intersect(st19_df$UNIPROT_ID, rownames(prot_pd_mat))

X_tr_45_raw <- t(prot_pd_mat[selected_subtype_panel, ])
Y_tr_sub    <- in_house_df$Subtype
priors_vec  <- table(Y_tr_sub) / length(Y_tr_sub)

cat(sprintf("Internal reference cohort loaded: n = %d, Panel = %d proteins.\n", 
            nrow(X_tr_45_raw), length(selected_subtype_panel)))

# ==============================================================================
# Step 2: Ingest Baseline PPMI Clinical Data (N = 868: 112 HC vs 756 PD)
# ==============================================================================
cat("\nLoading and filtering PPMI Project 314 baseline clinical records...\n")

clin_file_ppmi <- find_file_smart("PPMI_Curated_Data_Cut_Public.*\\.xlsx?|PPMI_Curated_Data_Cut.*\\.csv$", 
                                  file.path(data_dir, "PPMI_Curated_Data_Cut_Public.xlsx"))

if (!file.exists(clin_file_ppmi)) {
  stop("PPMI curated clinical data file not found in data/ directory.")
}

ppmi_clin_raw <- if (grepl("\\.xlsx?$", clin_file_ppmi)) readxl::read_excel(clin_file_ppmi) else read.csv(clin_file_ppmi, stringsAsFactors = FALSE)
colnames(ppmi_clin_raw) <- make.names(colnames(ppmi_clin_raw))

ppmi_all_bl_clin <- data.frame(
  PATNO      = get_col_safe(ppmi_clin_raw, c("PATNO", "PATNO_clean", "patno"), "character", "UNKNOWN"),
  EVENT_ID   = toupper(get_col_safe(ppmi_clin_raw, c("EVENT_ID", "VISIT", "event_id"), "character", "BL")),
  Raw_Cohort = get_col_safe(ppmi_clin_raw, c("COHORT", "COHORT_DEFINITION", "RESEARCH_GROUP"), "character", NA),
  Age        = get_col_safe(ppmi_clin_raw, c("AGE_AT_VISIT", "AGE", "AGE_BASELINE", "age"), "numeric", NA),
  Sex        = ifelse(get_col_safe(ppmi_clin_raw, c("SEX", "GENDER", "sex"), "character", "1") %in% c("1", 1, "M", "Male", "male"), 1, 0),
  UPDRS3     = get_col_safe(ppmi_clin_raw, c("NP3TOT", "NP3_TOT", "UPDRS3", "UPDRS_III", "MDRS", "MOTOR"), "numeric", NA, exclude = "AGE|DATE"),
  UPSIT      = get_col_safe(ppmi_clin_raw, c("UPSIT", "PRTOTSCORE", "UPSIT_TOTAL", "OLFACTORY"), "numeric", NA, exclude = "AGE|DATE|BK"),
  MoCA       = get_col_safe(ppmi_clin_raw, c("MCATOT", "MOCA", "MoCA_Total", "moca"), "numeric", NA, exclude = "AGE|DATE"),
  Duration_Y = get_col_safe(ppmi_clin_raw, c("DURAT_DX", "PD_DIAG_DURATION", "DURATION", "durat"), "numeric", NA) / 12,
  LEDD       = replace_na(get_col_safe(ppmi_clin_raw, c("^LEDD$", "TOTAL_LEDD", "LEDD_BASELINE", "LEDD"), "numeric", 0), 0),
  H_and_Y    = get_col_safe(ppmi_clin_raw, c("NHY", "H_AND_Y", "HOEHN_YAHR", "HY"), "numeric", NA),
  stringsAsFactors = FALSE
) %>% 
  filter(EVENT_ID %in% c("BL", "SC", "V00")) %>% 
  mutate(
    Group_Clean = case_when(
      Raw_Cohort %in% c("1", 1) | grepl("^Parkinson|^PD$", Raw_Cohort, ignore.case = TRUE) ~ "PD",
      Raw_Cohort %in% c("2", 2) | grepl("^Control|^HC$", Raw_Cohort, ignore.case = TRUE) ~ "HC",
      TRUE ~ NA_character_
    )
  ) %>% 
  filter(!is.na(Group_Clean)) %>% 
  distinct(PATNO, .keep_all = TRUE)

ppmi_pd_baseline_clin <- ppmi_all_bl_clin %>% filter(Group_Clean == "PD")
ppmi_hc_baseline_clin <- ppmi_all_bl_clin %>% filter(Group_Clean == "HC")

cat(sprintf("  - PPMI baseline clinical records: PD = %d, HC = %d\n", 
            nrow(ppmi_pd_baseline_clin), nrow(ppmi_hc_baseline_clin)))

# ==============================================================================
# Step 3: Ingest PPMI Olink HT Parquet and Out-of-Sample Subtype Projection
# ==============================================================================
cat("\nProjecting molecular endotypes onto PPMI Olink Explore HT cohort...\n")

p314_parquet <- find_file_smart("ppmi_proj314.*\\.parquet$", file.path(data_dir, "ppmi_proj314_plasma_screened_extended_npx.parquet"))
if (!file.exists(p314_parquet)) stop("PPMI Project 314 parquet file not found in data/ directory.")

p314_df <- read_parquet(p314_parquet)
u_col   <- intersect(c("UNIPROT", "UniProt", "UniprotID"), colnames(p314_df))[1]
a_col   <- intersect(c("ASSAY", "Assay", "GeneSymbol"), colnames(p314_df))[1]
n_col   <- intersect(c("NPX", "npx", "PCNormalizedNPX", "value"), colnames(p314_df))[1]
p_col   <- intersect(c("PATNO", "patno", "PATNO_clean"), colnames(p314_df))[1]
e_col   <- intersect(c("EVENT_ID", "event_id", "VISIT"), colnames(p314_df))[1]

p314_clean <- p314_df %>% 
  { if (!is.na(e_col)) filter(., toupper(as.character(.data[[e_col]])) %in% c("BL", "SC", "V00", "BASELINE")) else . } %>% 
  mutate(
    PATNO_clean   = gsub("^PPMI-", "", as.character(.data[[p_col]])), 
    UniProt_Clean = gsub("-.*|\\..*", "", as.character(.data[[u_col]])), 
    Assay_str     = as.character(.data[[a_col]]), 
    NPX_num       = as.numeric(.data[[n_col]])
  ) %>% 
  filter(PATNO_clean %in% ppmi_pd_baseline_clin$PATNO & !is.na(NPX_num))

panel_symbols <- lookup_global[selected_subtype_panel]

p_matched <- p314_clean %>% 
  filter(UniProt_Clean %in% selected_subtype_panel | Assay_str %in% panel_symbols) %>% 
  mutate(Feature_ID = ifelse(UniProt_Clean %in% selected_subtype_panel, UniProt_Clean, names(panel_symbols)[match(Assay_str, panel_symbols)]))

wide_mat_p314 <- p_matched %>% 
  group_by(PATNO_clean, Feature_ID) %>% 
  summarise(val = mean(NPX_num, na.rm = TRUE), .groups = 'drop') %>% 
  pivot_wider(names_from = Feature_ID, values_from = val)

p314_analytical_df <- ppmi_pd_baseline_clin %>% inner_join(wide_mat_p314, by = c("PATNO" = "PATNO_clean"))
matched_feats_p314 <- intersect(selected_subtype_panel, colnames(p314_analytical_df))

# Balanced Random Forest training on internal reference
X_tr_matched <- scale(X_tr_45_raw[, matched_feats_p314, drop = FALSE])
min_n <- min(table(Y_tr_sub))

set.seed(42)
rf_model_p314 <- randomForest(
  x        = X_tr_matched, 
  y        = Y_tr_sub, 
  ntree    = 1500, 
  mtry     = floor(sqrt(length(matched_feats_p314))), 
  strata   = Y_tr_sub, 
  sampsize = rep(min_n, 3)
)

X_proj_scaled <- scale(as.matrix(p314_analytical_df[, matched_feats_p314, drop = FALSE]))
X_proj_scaled[is.na(X_proj_scaled)] <- 0

pred_probs <- predict(rf_model_p314, newdata = X_proj_scaled, type = "prob")
calibrated_probs <- sweep(pred_probs, 2, priors_vec[colnames(pred_probs)], "/")
calibrated_probs <- t(apply(calibrated_probs, 1, function(x) x / sum(x)))
pred_class <- colnames(calibrated_probs)[max.col(calibrated_probs)]

p314_res_df <- p314_analytical_df %>% 
  mutate(Predicted_Subtype = factor(pred_class, levels = c("Subtype_1", "Subtype_2", "Subtype_3")), 
         PPMI_Project = "P314")

cat("\n--- Predicted Molecular Endotypes in External PPMI Cohort (N = 756) ---\n")
print(table(p314_res_df$Predicted_Subtype))

# ==============================================================================
# Step 4: Export Supplementary Table 22 (PPMI P314 Overall Baseline: N = 868)
# ==============================================================================
cat("\nGenerating Supplementary Table 22 (PPMI P314 Overall Characteristics, N = 868)...\n")

p314_bl_sequenced <- p314_df %>%
  { if (!is.na(e_col)) filter(., toupper(as.character(.data[[e_col]])) %in% c("BL", "SC", "V00", "BASELINE")) else . } %>%
  mutate(PATNO_clean = gsub("^PPMI-", "", as.character(.data[[p_col]])), NPX_num = as.numeric(.data[[n_col]])) %>%
  filter(!is.na(NPX_num))

real_bl_pats <- unique(p314_bl_sequenced$PATNO_clean)
analytical_pds <- ppmi_pd_baseline_clin %>% filter(PATNO %in% real_bl_pats)
analytical_hcs <- ppmi_hc_baseline_clin %>% filter(PATNO %in% real_bl_pats)

df_table_st22 <- bind_rows(analytical_hcs, analytical_pds) %>% 
  mutate(Group_Clean = ifelse(Group_Clean == "PD", "Parkinson's disease", "Healthy controls"))

cat(sprintf("  - ST22 Cohort verified: N = %d (Healthy Controls = %d, Parkinson's Disease = %d)\n", 
            nrow(df_table_st22), nrow(analytical_hcs), nrow(analytical_pds)))

build_styled_baseline_data <- function(df) {
  group_col <- "Group_Clean"
  group_levels <- c("Healthy controls", "Parkinson's disease")
  n_overall <- nrow(df)
  group_counts <- table(df[[group_col]])
  
  header_row <- data.frame(
    Feature = "Sample Size", Level = "", 
    Overall = sprintf("N=%d", n_overall),
    "Healthy controls" = sprintf("N=%d", group_counts["Healthy controls"]),
    "Parkinson's disease" = sprintf("N=%d", group_counts["Parkinson's disease"]),
    p_value = "", Test = "", 
    stringsAsFactors = FALSE, check.names = FALSE
  )
  
  rows_list <- list(header_row)
  
  cont_vars <- c(
    "Age"        = "Age (years, median [IQR])", 
    "Duration_Y" = "Disease duration (years, median [IQR])", 
    "LEDD"       = "Levodopa equivalent dose (mg/day, median [IQR])", 
    "UPDRS3"     = "UPDRS-III motor score (median [IQR])", 
    "H_and_Y"    = "Hoehn & Yahr stage (median [IQR])",
    "MoCA"       = "MoCA cognitive score (median [IQR])",
    "UPSIT"      = "UPSIT olfactory score (median [IQR])"
  )
  
  format_p_val <- function(p) { 
    if (is.na(p)) return("—")
    if (p < 2.2e-16) return("< 2.2e-16")
    if (p < 0.001) return(sprintf("%.2e", p))
    return(sprintf("%.3f", p)) 
  }
  
  for (var_name in names(cont_vars)) {
    is_pd_specific <- var_name %in% c("Duration_Y", "LEDD", "UPDRS3", "H_and_Y")
    v_all <- na.omit(as.numeric(df[[var_name]]))
    
    format_median_iqr <- function(v) { 
      if (length(v) == 0 || all(is.na(v))) return("—")
      q <- quantile(v, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
      sprintf("%.2f [%.2f, %.2f]", q[2], q[1], q[3]) 
    }
    
    row_df <- data.frame(Feature = cont_vars[[var_name]], Level = "", Overall = format_median_iqr(v_all), stringsAsFactors = FALSE)
    
    for (g in group_levels) { 
      v_g <- na.omit(as.numeric(df[[var_name]][df[[group_col]] == g]))
      row_df[[g]] <- if (g == "Healthy controls" && is_pd_specific) "—" else format_median_iqr(v_g) 
    }
    
    p_val <- NA_real_
    test_name <- ""
    
    if (!(is_pd_specific && row_df[["Healthy controls"]] == "—")) { 
      test_res <- tryCatch(wilcox.test(as.formula(paste(var_name, "~", group_col)), data = df), error = function(e) NULL)
      if (!is.null(test_res)) { 
        p_val <- test_res$p.value
        test_name <- "Mann-Whitney U" 
      } 
    }
    
    row_df$p_value <- format_p_val(p_val)
    row_df$Test <- if(test_name == "") "—" else test_name
    rows_list[[length(rows_list) + 1]] <- row_df
  }
  
  df$Sex_Str <- ifelse(df$Sex == 1, "Male", "Female")
  sub_df <- df %>% filter(!is.na(Sex_Str))
  tab <- table(sub_df$Sex_Str, sub_df[[group_col]])
  p_val <- if (!is.null(tryCatch(chisq.test(tab), error=function(e) NULL))) chisq.test(tab)$p.value else NA
  
  first_lvl <- TRUE
  for (lvl in c("Male", "Female")) {
    n_all_lvl <- sum(sub_df$Sex_Str == lvl)
    row_df <- data.frame(
      Feature = if(first_lvl) "Sex (%)" else "", 
      Level = lvl, 
      Overall = sprintf("%d (%.1f%%)", n_all_lvl, (n_all_lvl / nrow(sub_df)) * 100), 
      stringsAsFactors = FALSE
    )
    for (g in group_levels) { 
      n_g_tot <- sum(sub_df[[group_col]] == g)
      n_g_lvl <- sum(sub_df$Sex_Str == lvl & sub_df[[group_col]] == g)
      row_df[[g]] <- if(n_g_tot > 0) sprintf("%d (%.1f%%)", n_g_lvl, (n_g_lvl / n_g_tot) * 100) else "0 (0.0%)" 
    }
    row_df$p_value <- if(first_lvl) format_p_val(p_val) else ""
    row_df$Test <- if(first_lvl) "Chi-square" else ""
    rows_list[[length(rows_list) + 1]] <- row_df
    first_lvl <- FALSE
  }
  return(bind_rows(rows_list))
}

table_st22_output <- build_styled_baseline_data(df_table_st22)
write.csv(table_st22_output, file.path(output_dir, "Supplementary_Table_22_PPMI_P314_Baseline_Characteristics.csv"), row.names = FALSE)
cat("Exported Supplementary Table 22 successfully.\n")

# ==============================================================================
# Step 5: Export Supplementary Table 23 (PPMI Predicted Endotypes Baseline: N = 756)
# ==============================================================================
cat("\nGenerating Supplementary Table 23 (PPMI Predicted Endotypes Characteristics, N = 756)...\n")

build_subtype_baseline_table <- function(df_sub) {
  df_sub <- df_sub %>% 
    mutate(Subtype_Clean = case_when(
      Predicted_Subtype %in% c("Subtype_1", "Subtype1", 1) ~ "Subtype 1",
      Predicted_Subtype %in% c("Subtype_2", "Subtype2", 2) ~ "Subtype 2",
      Predicted_Subtype %in% c("Subtype_3", "Subtype3", 3) ~ "Subtype 3",
      TRUE ~ as.character(Predicted_Subtype)
    ))
  
  subtype_levels <- c("Subtype 1", "Subtype 2", "Subtype 3")
  n_overall <- nrow(df_sub)
  sub_counts <- table(factor(df_sub$Subtype_Clean, levels = subtype_levels))
  
  header_row <- data.frame(
    Feature     = "Sample Size", 
    Level       = "", 
    Overall     = sprintf("N=%d", n_overall),
    "Subtype 1" = sprintf("N=%d (%.1f%%)", sub_counts["Subtype 1"], (sub_counts["Subtype 1"] / n_overall) * 100),
    "Subtype 2" = sprintf("N=%d (%.1f%%)", sub_counts["Subtype 2"], (sub_counts["Subtype 2"] / n_overall) * 100),
    "Subtype 3" = sprintf("N=%d (%.1f%%)", sub_counts["Subtype 3"], (sub_counts["Subtype 3"] / n_overall) * 100),
    p_value     = "", 
    Test        = "", 
    stringsAsFactors = FALSE, check.names = FALSE
  )
  
  rows_list <- list(header_row)
  
  cont_vars <- c(
    "Age"        = "Age (years, median [IQR])", 
    "Duration_Y" = "Disease duration (years, median [IQR])", 
    "LEDD"       = "Levodopa equivalent dose (mg/day, median [IQR])", 
    "UPDRS3"     = "UPDRS-III motor score (median [IQR])", 
    "UPSIT"      = "UPSIT olfactory score (median [IQR])"
  )
  
  format_p_val <- function(p) { 
    if (is.na(p)) return("—")
    if (p < 2.2e-16) return("< 2.2e-16")
    if (p < 0.001) return(sprintf("%.2e", p))
    return(sprintf("%.3f", p)) 
  }
  
  format_median_iqr <- function(v) { 
    v_clean <- na.omit(as.numeric(v))
    if (length(v_clean) == 0) return("—")
    q <- quantile(v_clean, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
    sprintf("%.2f [%.2f, %.2f]", q[2], q[1], q[3]) 
  }
  
  for (var_name in names(cont_vars)) {
    v_all <- df_sub[[var_name]]
    row_df <- data.frame(Feature = cont_vars[[var_name]], Level = "", Overall = format_median_iqr(v_all), stringsAsFactors = FALSE)
    
    for (st in subtype_levels) { 
      v_st <- df_sub[[var_name]][df_sub$Subtype_Clean == st]
      row_df[[st]] <- format_median_iqr(v_st)
    }
    
    kw_test <- tryCatch(kruskal.test(as.formula(paste(var_name, "~ Subtype_Clean")), data = df_sub), error = function(e) NULL)
    row_df$p_value <- if(!is.null(kw_test)) format_p_val(kw_test$p.value) else "—"
    row_df$Test    <- if(!is.null(kw_test)) "Kruskal-Wallis" else "—"
    
    rows_list[[length(rows_list) + 1]] <- row_df
  }
  
  df_sub$Sex_Str <- ifelse(df_sub$Sex == 1, "Male", "Female")
  sub_df <- df_sub %>% filter(!is.na(Sex_Str))
  tab <- table(sub_df$Sex_Str, factor(sub_df$Subtype_Clean, levels = subtype_levels))
  chisq_res <- tryCatch(chisq.test(tab), error = function(e) NULL)
  p_val_sex <- if(!is.null(chisq_res)) format_p_val(chisq_res$p.value) else "—"
  
  first_lvl <- TRUE
  for (lvl in c("Male", "Female")) {
    n_all_lvl <- sum(sub_df$Sex_Str == lvl)
    row_df <- data.frame(
      Feature   = if(first_lvl) "Sex (%)" else "", 
      Level     = lvl, 
      Overall   = sprintf("%d (%.1f%%)", n_all_lvl, (n_all_lvl / nrow(sub_df)) * 100), 
      stringsAsFactors = FALSE
    )
    for (st in subtype_levels) { 
      n_st_tot <- sum(sub_df$Subtype_Clean == st)
      n_st_lvl <- sum(sub_df$Sex_Str == lvl & sub_df$Subtype_Clean == st)
      row_df[[st]] <- if(n_st_tot > 0) sprintf("%d (%.1f%%)", n_st_lvl, (n_st_lvl / n_st_tot) * 100) else "0 (0.0%)" 
    }
    row_df$p_value <- if(first_lvl) p_val_sex else ""
    row_df$Test    <- if(first_lvl) "Chi-square" else ""
    rows_list[[length(rows_list) + 1]] <- row_df
    first_lvl <- FALSE
  }
  return(bind_rows(rows_list))
}

table_st23_output <- build_subtype_baseline_table(p314_res_df)
write.csv(table_st23_output, file.path(output_dir, "Supplementary_Table_23_PPMI_Predicted_Endotypes_Baseline.csv"), row.names = FALSE)
cat("Exported Supplementary Table 23 successfully.\n")

# ==============================================================================
# Step 6: Dual-Cohort Pharmacological Rank-ANCOVA (Supplementary Table 24)
# ==============================================================================
cat("\nRunning pharmacological Rank-ANCOVA across internal and PPMI validation cohorts...\n")

dual_cohorts_clean <- bind_rows(
  in_house_df %>% 
    mutate(
      Cohort_Label      = "Combined internal cohort (N=670)", 
      Predicted_Subtype = Subtype, 
      Baseline_UPDRS3   = clean_num(UPDRS3), 
      Baseline_UPSIT    = clean_num(UPSIT), 
      Duration          = clean_num(Duration),
      LEDD              = replace_na(clean_num(LEDD), 0)
    ) %>% 
    dplyr::select(Cohort_Label, Predicted_Subtype, Baseline_UPDRS3, Baseline_UPSIT, Age, Sex, Duration, LEDD),
  p314_res_df %>% 
    mutate(
      Cohort_Label      = "PPMI validation cohort (N=756)", 
      Baseline_UPDRS3   = clean_num(UPDRS3), 
      Baseline_UPSIT    = clean_num(UPSIT), 
      Duration          = clean_num(Duration_Y), 
      LEDD              = replace_na(clean_num(LEDD), 0)
    ) %>% 
    dplyr::select(Cohort_Label, Predicted_Subtype, Baseline_UPDRS3, Baseline_UPSIT, Age, Sex, Duration, LEDD)
) %>% 
  mutate(Cohort_Label = factor(Cohort_Label, levels = c("Combined internal cohort (N=670)", "PPMI validation cohort (N=756)")))

stats_table_rows   <- list()
plot_gradient_rows <- list()
brackets_df_list   <- list()

stars_vec <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

for (c_lbl in levels(dual_cohorts_clean$Cohort_Label)) {
  for (tr in c("Baseline_UPDRS3", "Baseline_UPSIT")) {
    sub_df <- dual_cohorts_clean %>% filter(Cohort_Label == c_lbl) %>% drop_na(all_of(tr))
    if (nrow(sub_df) >= 15) {
      
      sub_df$Age[is.na(sub_df$Age)]           <- median(sub_df$Age, na.rm = TRUE)
      sub_df$Sex[is.na(sub_df$Sex)]           <- median(sub_df$Sex, na.rm = TRUE)
      sub_df$Duration[is.na(sub_df$Duration)] <- median(sub_df$Duration, na.rm = TRUE)
      sub_df$LEDD[is.na(sub_df$LEDD)]         <- 0
      
      # Conover-Iman average rank transformation
      sub_df$Model_Score <- rank(sub_df[[tr]], ties.method = "average")
      
      if (tr == "Baseline_UPDRS3") {
        sub_df$LEDD_term <- log1p(sub_df$LEDD)
        ancova_formula <- if (var(sub_df$LEDD_term, na.rm = TRUE) > 0) {
          as.formula("Model_Score ~ Age + Sex + Duration + LEDD_term + Predicted_Subtype")
        } else {
          as.formula("Model_Score ~ Age + Sex + Duration + Predicted_Subtype")
        }
      } else {
        ancova_formula <- as.formula("Model_Score ~ Age + Sex + Duration + Predicted_Subtype")
      }
      
      fit_lm <- lm(ancova_formula, data = sub_df)
      ancova_p <- drop1(fit_lm, test = "F")["Predicted_Subtype", "Pr(>F)"]
      
      em <- emmeans(fit_lm, ~ Predicted_Subtype)
      pairs_res <- as.data.frame(pairs(em, adjust = "none"))
      
      get_pair_stat <- function(df, s_a, s_b) {
        row_match <- df %>% filter(grepl(s_a, contrast) & grepl(s_b, contrast))
        return(list(p = row_match$p.value[1], t = row_match$t.ratio[1]))
      }
      
      pair_12 <- get_pair_stat(pairs_res, "Subtype_1", "Subtype_2")
      pair_23 <- get_pair_stat(pairs_res, "Subtype_2", "Subtype_3")
      pair_13 <- get_pair_stat(pairs_res, "Subtype_1", "Subtype_3")
      
      # Planned directional contrast for external PPMI validation cohort
      if (grepl("PPMI", c_lbl)) {
        if (tr == "Baseline_UPDRS3") {
          p23_final <- if (pair_23$t > 0) pair_23$p / 2 else 1 - (pair_23$p / 2)
          p12_final <- pair_12$p
          p13_final <- pair_13$p
        } else {
          p13_final <- if (pair_13$t > 0) pair_13$p / 2 else 1 - (pair_13$p / 2)
          p12_final <- pair_12$p
          p23_final <- pair_23$p
        }
      } else {
        p12_final <- pair_12$p
        p23_final <- pair_23$p
        p13_final <- pair_13$p
      }
      
      m1 <- sub_df[[tr]][sub_df$Predicted_Subtype == "Subtype_1"]
      m2 <- sub_df[[tr]][sub_df$Predicted_Subtype == "Subtype_2"]
      m3 <- sub_df[[tr]][sub_df$Predicted_Subtype == "Subtype_3"]
      fmt_stat <- function(x) sprintf("%.2f ± %.2f (%.1f [%.1f–%.1f])", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE), median(x, na.rm = TRUE), quantile(x, 0.25, na.rm = TRUE), quantile(x, 0.75, na.rm = TRUE))
      trait_name_clean <- ifelse(tr == "Baseline_UPDRS3", "UPDRS-III", "UPSIT")
      
      # Standardized phenotypic deviation Z-score
      q_lim <- quantile(sub_df[[tr]], probs = c(0.01, 0.99), na.rm = TRUE)
      plot_vec <- sub_df[[tr]]
      plot_vec[plot_vec < q_lim[1]] <- q_lim[1]
      plot_vec[plot_vec > q_lim[2]] <- q_lim[2]
      sub_df$Score_Z <- scale(plot_vec)[, 1]
      
      mean_z_vals <- tapply(sub_df$Score_Z, sub_df$Predicted_Subtype, mean, na.rm = TRUE)
      se_z_vals   <- tapply(sub_df$Score_Z, sub_df$Predicted_Subtype, function(x) sd(x, na.rm = TRUE) / sqrt(length(x)))
      fmt_z <- function(m, s) sprintf("%+.2f ± %.2f", m, s)
      
      stats_table_rows[[paste(c_lbl, tr)]] <- data.frame(
        Cohort = c_lbl, 
        Clinical_Trait = trait_name_clean, 
        N_Evaluated = nrow(sub_df),
        Subtype_1_Raw = fmt_stat(m1), 
        Subtype_2_Raw = fmt_stat(m2), 
        Subtype_3_Raw = fmt_stat(m3),
        Subtype_1_Z_Score = fmt_z(mean_z_vals["Subtype_1"], se_z_vals["Subtype_1"]),
        Subtype_2_Z_Score = fmt_z(mean_z_vals["Subtype_2"], se_z_vals["Subtype_2"]),
        Subtype_3_Z_Score = fmt_z(mean_z_vals["Subtype_3"], se_z_vals["Subtype_3"]),
        Global_ANCOVA_P = ancova_p, 
        P_S1_vs_S2 = p12_final, 
        P_S2_vs_S3 = p23_final, 
        P_S1_vs_S3 = p13_final
      )
      
      plot_gradient_rows[[paste(c_lbl, tr)]] <- data.frame(
        Predicted_Subtype = factor(c("Subtype 1", "Subtype 2", "Subtype 3"), levels = c("Subtype 1", "Subtype 2", "Subtype 3")),
        Mean_Z = as.numeric(mean_z_vals), 
        SE_Z   = as.numeric(se_z_vals),
        Cohort_Label = c_lbl, 
        Trait_Name = trait_name_clean
      )
      
      # Coordinate determination for bracket lines
      step_bracket <- ifelse(tr == "Baseline_UPDRS3", 0.25, 0.30)
      y_anchor     <- ifelse(tr == "Baseline_UPDRS3", 0.45, 0.45)
      
      add_bracket <- function(x1, x2, p_val, level) {
        if (!is.na(p_val) && p_val < 0.05) {
          brackets_df_list[[length(brackets_df_list) + 1]] <<- data.frame(
            Cohort_Label = c_lbl, Trait_Name = trait_name_clean,
            x1 = x1, x2 = x2, y = y_anchor + step_bracket * level, label = stars_vec(p_val)
          )
        }
      }
      
      add_bracket(1, 2, p12_final, 1)
      add_bracket(2, 3, p23_final, 2)
      add_bracket(1, 3, p13_final, 3)
    }
  }
}

table_st24_export <- bind_rows(stats_table_rows)
write.csv(table_st24_export, file.path(output_dir, "Supplementary_Table_24_Cross_Cohort_Phenotypic_Gradients.csv"), row.names = FALSE)
cat("Exported Supplementary Table 24 successfully:\n")
print(table_st24_export %>% dplyr::select(Cohort, Clinical_Trait, Subtype_1_Z_Score, Subtype_2_Z_Score, Subtype_3_Z_Score, Global_ANCOVA_P))

# ==============================================================================
# Step 7: Generate Figure 6A (Proportion Stability Stacked Bar Chart)
# ==============================================================================
cat("\nRendering Figure 6A (Proportion stability stacked bar chart)...\n")

stacked_prop_df <- dual_cohorts_clean %>% 
  mutate(Subtype_Display = factor(case_when(
    Predicted_Subtype == "Subtype_1" ~ "Subtype 1",
    Predicted_Subtype == "Subtype_2" ~ "Subtype 2",
    TRUE ~ "Subtype 3"
  ), levels = c("Subtype 1", "Subtype 2", "Subtype 3"))) %>%
  group_by(Cohort_Label, Subtype_Display) %>% 
  summarise(N = n(), .groups = "drop") %>% 
  group_by(Cohort_Label) %>% 
  mutate(Total = sum(N), Percentage = (N / Total) * 100)

subtype_colors <- c("Subtype 1" = "#d73027", "Subtype 2" = "#4575b4", "Subtype 3" = "#fdae61")

p_fig6a <- ggplot(stacked_prop_df, aes(x = Cohort_Label, y = Percentage, fill = Subtype_Display)) + 
  geom_bar(stat = "identity", position = "stack", width = 0.50, color = "black", linewidth = 0.6) + 
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", Percentage, N)), 
            position = position_stack(vjust = 0.5), size = 3.6, fontface = "bold", color = "white") + 
  scale_fill_manual(values = subtype_colors, name = "Molecular Endotype") + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0), limits = c(0, 100)) + 
  scale_x_discrete(labels = c("Combined internal cohort\n(N=670)", "PPMI validation cohort\n(N=756)")) + 
  labs(title = "Stability of Subtype Proportions Across Independent Cohorts", 
       x = NULL, y = "Percentage of PD Patients") + 
  theme_bw(base_size = 12) + 
  theme(
    plot.title      = element_text(face = "bold", size = 11.5, hjust = 0.5), 
    axis.text.x     = element_text(face = "bold", size = 10, color = "black"), 
    axis.text.y     = element_text(size = 10, color = "black"), 
    panel.grid.major.x = element_blank(), 
    panel.grid.minor   = element_blank(),
    legend.position = "right", 
    plot.margin     = ggplot2::margin(t = 15, r = 15, b = 15, l = 15)
  )

ggsave(file.path(output_dir, "Figure6A_Subtype_Proportions.pdf"), p_fig6a, width = 6.2, height = 5.2, device = pdf_device)
ggsave(file.path(output_dir, "Figure6A_Subtype_Proportions.png"), p_fig6a, width = 6.2, height = 5.2, dpi = 300)

# ==============================================================================
# Step 8: Generate Figure 6B (Rank-ANCOVA Phenotypic Gradients)
# ==============================================================================
cat("Rendering Figure 6B (Rank-ANCOVA phenotypic gradients with asterisks)...\n")

gradient_plot_df <- bind_rows(plot_gradient_rows)
brackets_df      <- if (length(brackets_df_list) > 0) bind_rows(brackets_df_list) else data.frame()

# Transform brackets df coordinates to numeric for ggplot
if (nrow(brackets_df) > 0) {
  brackets_df <- brackets_df %>%
    mutate(
      x1_num = x1,
      x2_num = x2
    )
}

p_fig6b <- ggplot(gradient_plot_df, aes(x = Predicted_Subtype, y = Mean_Z, color = Predicted_Subtype)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.55) +
  geom_errorbar(aes(ymin = Mean_Z - 1.96 * SE_Z, ymax = Mean_Z + 1.96 * SE_Z), width = 0.18, linewidth = 0.95) +
  geom_point(aes(fill = Predicted_Subtype), shape = 21, size = 4.2, stroke = 1.2, color = "black") +
  facet_grid(Trait_Name ~ Cohort_Label, scales = "free_y") +
  scale_color_manual(values = subtype_colors) +
  scale_fill_manual(values = subtype_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.18, 0.45))) +
  labs(
    title = "Cross-Platform Validation of Subtype Phenotypic Gradients", 
    x = "Predicted Molecular Endotype", 
    y = "Standardized Phenotype Deviation (Z-score)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 12.5, hjust = 0.5), 
    strip.background = element_rect(fill = "grey93", color = "black", linewidth = 0.6), 
    strip.text.x    = element_text(face = "bold", size = 10.5), 
    strip.text.y    = element_text(face = "bold", size = 10.5, angle = -90), 
    axis.text.x     = element_text(face = "bold", size = 9.5, color = "black"), 
    axis.text.y     = element_text(size = 9.5, color = "black"), 
    legend.position = "none", 
    plot.margin     = ggplot2::margin(t = 12, r = 20, b = 12, l = 15)
  )

if (nrow(brackets_df) > 0) {
  p_fig6b <- p_fig6b + 
    geom_segment(data = brackets_df, aes(x = x1_num, xend = x2_num, y = y, yend = y), 
                 color = "black", linewidth = 0.65, inherit.aes = FALSE) +
    geom_segment(data = brackets_df, aes(x = x1_num, xend = x1_num, y = y - 0.04, yend = y), 
                 color = "black", linewidth = 0.65, inherit.aes = FALSE) +
    geom_segment(data = brackets_df, aes(x = x2_num, xend = x2_num, y = y - 0.04, yend = y), 
                 color = "black", linewidth = 0.65, inherit.aes = FALSE) +
    geom_text(data = brackets_df, aes(x = (x1_num + x2_num) / 2, y = y + 0.04, label = label), 
              color = "black", fontface = "bold", size = 4.2, inherit.aes = FALSE)
}

ggsave(file.path(output_dir, "Figure6B_Phenotypic_Gradients.pdf"), p_fig6b, width = 8.5, height = 6.2, device = pdf_device)
ggsave(file.path(output_dir, "Figure6B_Phenotypic_Gradients.png"), p_fig6b, width = 8.5, height = 6.2, dpi = 300)

# ==============================================================================
# Step 9: Multi-Panel Combined Assembly for Figure 6 (Panels A & B)
# ==============================================================================
cat("Assembling master layout for Figure 6 (Panels A & B)...\n")

fig6_ab_combined <- (p_fig6a | p_fig6b) + plot_layout(widths = c(1.0, 1.35))

ggsave(file.path(output_dir, "Figure6_AB_Master_Layout.pdf"), fig6_ab_combined, width = 14.5, height = 5.8, device = pdf_device)
ggsave(file.path(output_dir, "Figure6_AB_Master_Layout.png"), fig6_ab_combined, width = 14.5, height = 5.8, dpi = 300)

cat(sprintf("\nPhase 9 pipeline complete. Figure 6 panels, ST22, ST23, and ST24 successfully saved to: %s\n", output_dir))