# ==============================================================================
# Script: 05_ppmi_diagnostic_validation.R
# Purpose: Cross-platform diagnostic validation across independent PPMI cohorts (ST14b)
# Project: Large-scale plasma proteomics identifies molecularly distinct PD endotypes
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
  library(ggpubr)
  library(org.Hs.eg.db)
})

# Directory setup
data_dir       <- "data"
lasso_dir      <- "results/03_figure2_diagnostic_model"
output_dir     <- "results/05_ppmi_diagnostic_validation"
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

cat("\n=== Phase 4: PPMI Cross-Platform Diagnostic Validation ===\n")

# ==============================================================================
# Step 1: Load Diagnostic Panel Biomarkers (Table ST14a)
# ==============================================================================
cat("Loading diagnostic panel biomarkers...\n")

s4_file <- file.path(lasso_dir, "Table_ST14a_Diagnostic_Panel_OR_Parameters.csv")
if (!file.exists(s4_file)) {
  s4_candidates <- list.files(lasso_dir, pattern = "Table_.*Protein_Panel_OR.*\\.csv$", full.names = TRUE)
  s4_file <- if (length(s4_candidates) > 0) s4_candidates[1] else ""
}

if (!file.exists(s4_file)) stop("Required diagnostic panel file not found. Please run Script 03 first.")

s4_raw <- read.csv(s4_file, stringsAsFactors = FALSE)
clin_covars_pattern <- "Intercept|age|sex|BMI|hepatic|kidney|CV_Met"
s4_clean <- s4_raw %>%
  filter(!grepl(clin_covars_pattern, Symbol, ignore.case = TRUE) & 
           !grepl(clin_covars_pattern, UNIPROT_ID, ignore.case = TRUE))

panel_uniprots <- unique(gsub("`", "", as.character(s4_clean$UNIPROT_ID)))
panel_symbols  <- unique(as.character(s4_clean$Symbol))

cat(sprintf("Loaded %d panel biomarkers for external cross-platform matching.\n", length(panel_uniprots)))

# ==============================================================================
# Step 2: Ingest PPMI Baseline Clinical Metadata
# ==============================================================================
cat("Parsing baseline PPMI clinical records...\n")

clin_file_path <- find_file_smart("PPMI_Curated_Data_Cut.*\\.xlsx$", file.path(data_dir, "PPMI_Curated_Data_Cut_Public.xlsx"))
if (!file.exists(clin_file_path)) {
  clin_file_path <- find_file_smart("PPMI_Curated_Data_Cut.*\\.csv$", file.path(data_dir, "PPMI_Curated_Data_Cut_Public.csv"))
}

ppmi_clin_raw <- if (grepl("\\.xlsx$", clin_file_path)) readxl::read_excel(clin_file_path) else read.csv(clin_file_path, stringsAsFactors = FALSE)
colnames(ppmi_clin_raw) <- make.names(colnames(ppmi_clin_raw))

raw_patno  <- get_col_safe(ppmi_clin_raw, c("^PATNO$", "^patno$", "PATNO_clean"), "character", "UNKNOWN")
raw_event  <- get_col_safe(ppmi_clin_raw, c("^EVENT_ID$", "^event_id$", "^VISIT$"), "character", "BL")
raw_cohort <- get_col_safe(ppmi_clin_raw, c("^COHORT$", "^COHORT_DEFINITION$", "^RESEARCH_GROUP$"), "character", NA)
raw_age    <- get_col_safe(ppmi_clin_raw, c("AGE_AT_VISIT", "^AGE$", "AGE_BASELINE", "BL_AGE"), "numeric", 65)
raw_sex    <- get_col_safe(ppmi_clin_raw, c("^SEX$", "^GENDER$"), "character", "1")

# Restrict to baseline events (BL) and confirmed PD vs HC
ppmi_pure_baseline_clin <- data.frame(
  PATNO      = raw_patno,
  EVENT_ID   = toupper(raw_event),
  Raw_Cohort = raw_cohort,
  Age        = ifelse(is.na(raw_age), 65, raw_age),
  Sex        = ifelse(raw_sex %in% c("1", 1, "M", "Male", "male"), 1, 0),
  stringsAsFactors = FALSE
) %>%
  filter(EVENT_ID %in% c("BL", "SC", "V00", "BASELINE")) %>%
  mutate(
    Group_Clean = case_when(
      Raw_Cohort %in% c("1", 1) | grepl("^Parkinson|^PD$", Raw_Cohort, ignore.case = TRUE) ~ "PD",
      Raw_Cohort %in% c("2", 2) | grepl("^Control|^HC$", Raw_Cohort, ignore.case = TRUE) ~ "HC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group_Clean)) %>%
  distinct(PATNO, .keep_all = TRUE) %>%
  mutate(group = factor(Group_Clean, levels = c("HC", "PD")))

cat(sprintf("Identified %d baseline PPMI clinical participants (HC: %d, PD: %d).\n",
            nrow(ppmi_pure_baseline_clin), 
            sum(ppmi_pure_baseline_clin$group == "HC", na.rm = TRUE), 
            sum(ppmi_pure_baseline_clin$group == "PD", na.rm = TRUE)))

# ==============================================================================
# Step 3: Load Independent PPMI Projects
# ==============================================================================
cat("Ingesting proteomics data across independent PPMI cohorts...\n")

p314_path   <- find_file_smart("ppmi_proj314.*\\.parquet$", file.path(data_dir, "ppmi_proj314_plasma_screened_extended_npx.parquet"))
p293_path   <- find_file_smart("ppmi_proj293.*\\.parquet$", file.path(data_dir, "ppmi_proj293_plasma_screened_extended_npx.parquet"))
p9000_files <- list.files(data_dir, pattern = "PPMI_Project_9000_Plasma.*\\.csv$", full.names = TRUE, recursive = TRUE)

p314_df  <- if(file.exists(p314_path)) read_parquet(p314_path) else NULL
p293_df  <- if(file.exists(p293_path)) read_parquet(p293_path) else NULL
p9000_df <- if(length(p9000_files) > 0) bind_rows(lapply(p9000_files, read.csv, stringsAsFactors = FALSE)) else NULL

# Build gene mapping lookup
mapped_symbols_db <- suppressMessages(suppressWarnings(
  mapIds(org.Hs.eg.db, keys = panel_uniprots, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
))
lookup_panel <- setNames(as.character(mapped_symbols_db), panel_uniprots)
lookup_panel[is.na(lookup_panel)] <- names(lookup_panel)[is.na(lookup_panel)]

ppmi_projects <- list(
  "Project 314 (Olink Flagship)" = p314_df,
  "Project 293 (Olink Pilot)"    = p293_df,
  "Project 9000 (Targeted MS)"   = p9000_df
)

# ==============================================================================
# Step 4: Cross-Platform Diagnostic Performance Evaluation Engine
# ==============================================================================
cat("Evaluating incremental diagnostic utility against clinical baselines...\n")

evaluate_diagnostic_performance <- function(proj_df, proj_name) {
  if (is.null(proj_df) || nrow(proj_df) == 0) return(NULL)
  
  u_col <- intersect(c("UNIPROT", "UniProt", "UniprotID"), colnames(proj_df))[1]
  a_col <- intersect(c("ASSAY", "Assay", "GeneSymbol"), colnames(proj_df))[1]
  n_col <- intersect(c("NPX", "npx", "PCNormalizedNPX", "value"), colnames(proj_df))[1]
  p_col <- intersect(c("PATNO", "patno", "PATNO_clean"), colnames(proj_df))[1]
  e_col <- intersect(c("EVENT_ID", "event_id", "VISIT"), colnames(proj_df))[1]
  
  # Filter strictly for baseline plasma draws
  p_clean <- proj_df %>%
    { if (!is.na(e_col)) filter(., toupper(as.character(.data[[e_col]])) %in% c("BL", "SC", "V00", "BASELINE")) else . } %>%
    mutate(
      PATNO_clean   = gsub("^PPMI-", "", as.character(.data[[p_col]])),
      UniProt_Clean = gsub("-.*|\\..*", "", as.character(.data[[u_col]])),
      Assay_str     = as.character(.data[[a_col]]),
      NPX_num       = as.numeric(.data[[n_col]])
    ) %>%
    filter(PATNO_clean %in% ppmi_pure_baseline_clin$PATNO & !is.na(NPX_num))
  
  # Match panel biomarkers
  p_matched <- p_clean %>%
    filter(UniProt_Clean %in% panel_uniprots | Assay_str %in% panel_symbols) %>%
    mutate(
      Feature_ID = ifelse(UniProt_Clean %in% panel_uniprots, UniProt_Clean, names(lookup_panel)[match(Assay_str, lookup_panel)])
    )
  
  wide_mat <- p_matched %>%
    group_by(PATNO_clean, Feature_ID) %>%
    summarise(val = mean(NPX_num, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Feature_ID, values_from = val)
  
  analytical_df <- ppmi_pure_baseline_clin %>% inner_join(wide_mat, by = c("PATNO" = "PATNO_clean"))
  
  n_pd <- sum(analytical_df$group == "PD", na.rm = TRUE)
  n_hc <- sum(analytical_df$group == "HC", na.rm = TRUE)
  matched_feats <- intersect(panel_uniprots, colnames(analytical_df))
  
  if (n_hc < 5 || n_pd < 5 || length(matched_feats) == 0) return(NULL)
  
  for (f in matched_feats) {
    analytical_df[[f]] <- as.numeric(scale(analytical_df[[f]]))
    analytical_df[[f]][is.na(analytical_df[[f]])] <- 0
  }
  
  # Model 1: Clinical baseline (Age + Sex)
  fit_clin <- glm(group ~ Age + Sex, data = analytical_df, family = binomial)
  prob_clin <- predict(fit_clin, type = "response")
  roc_clin <- roc(response = analytical_df$group, predictor = prob_clin, levels = c("HC", "PD"), direction = "<", quiet = TRUE, ci = TRUE)
  
  # Model 2: Available panel proteins only
  form_prot <- as.formula(paste("group ~", paste(paste0("`", matched_feats, "`"), collapse = " + ")))
  fit_prot <- glm(form_prot, data = analytical_df, family = binomial)
  prob_prot <- predict(fit_prot, type = "response")
  roc_prot <- roc(response = analytical_df$group, predictor = prob_prot, levels = c("HC", "PD"), direction = "<", quiet = TRUE, ci = TRUE)
  
  # Model 3: Combined model (Clinical + Available Proteins)
  form_comb <- as.formula(paste("group ~ Age + Sex +", paste(paste0("`", matched_feats, "`"), collapse = " + ")))
  fit_comb <- glm(form_comb, data = analytical_df, family = binomial)
  prob_comb <- predict(fit_comb, type = "response")
  roc_comb <- roc(response = analytical_df$group, predictor = prob_comb, levels = c("HC", "PD"), direction = "<", quiet = TRUE, ci = TRUE)
  
  delong_res <- roc.test(roc_comb, roc_clin, method = "delong")
  
  return(list(
    Project         = proj_name,
    N_Total         = nrow(analytical_df),
    N_PD            = n_pd,
    N_HC            = n_hc,
    Matched_Count   = length(matched_feats),
    Matched_Genes   = paste(lookup_panel[matched_feats], collapse = ", "),
    ROC_Clin        = roc_clin,
    ROC_Prot        = roc_prot,
    ROC_Comb        = roc_comb,
    DeLong_P        = delong_res$p.value
  ))
}

p314_eval <- evaluate_diagnostic_performance(ppmi_projects[["Project 314 (Olink Flagship)"]], "Project 314 (Olink Flagship)")
p293_eval <- evaluate_diagnostic_performance(ppmi_projects[["Project 293 (Olink Pilot)"]], "Project 293 (Olink Pilot)")
p9000_eval <- evaluate_diagnostic_performance(ppmi_projects[["Project 9000 (Targeted MS)"]], "Project 9000 (Targeted MS)")

diag_eval_list <- list(p314_eval, p293_eval, p9000_eval) %>% compact()

# ==============================================================================
# Step 5: Export Supplementary Table 14b and Multi-Project ROC Plots
# ==============================================================================
cat("Generating Supplementary Table 14b and multi-project ROC curves...\n")

table_s14b_summary <- bind_rows(lapply(diag_eval_list, function(x) {
  data.frame(
    Project = x$Project,
    Platform = ifelse(grepl("Olink", x$Project), "Olink PEA", "Targeted MS"),
    Total_Patients = x$N_Total,
    PD_Cases = x$N_PD,
    HC_Controls = x$N_HC,
    Matched_Proteins = sprintf("%d / %d", x$Matched_Count, length(panel_uniprots)),
    Matched_Gene_List = x$Matched_Genes,
    Clinical_Baseline_AUC = sprintf("%.3f (%.3f–%.3f)", as.numeric(auc(x$ROC_Clin)), x$ROC_Clin$ci[1], x$ROC_Clin$ci[3]),
    Proteins_Only_AUC     = sprintf("%.3f (%.3f–%.3f)", as.numeric(auc(x$ROC_Prot)), x$ROC_Prot$ci[1], x$ROC_Prot$ci[3]),
    Combined_Model_AUC    = sprintf("%.3f (%.3f–%.3f)", as.numeric(auc(x$ROC_Comb)), x$ROC_Comb$ci[1], x$ROC_Comb$ci[3]),
    DeLong_P_Value        = sprintf("%.4e", x$DeLong_P)
  )
}))

table_s14b_path <- file.path(output_dir, "Supplementary_Table_14b_PPMI_Diagnostic_Validation.csv")
write.csv(table_s14b_summary, table_s14b_path, row.names = FALSE)

# Generate multi-panel ROC plot
pdf(file.path(output_dir, "Supplementary_Figure_PPMI_Diagnostic_ROC.pdf"), width = 12.5, height = 4.5)
par(mfrow = c(1, length(diag_eval_list)), mar = c(4.5, 4.5, 3.5, 1.5))

for (ev in diag_eval_list) {
  plot(ev$ROC_Comb, col = "#d73027", lwd = 3.2, 
       main = sprintf("%s (N=%d)", ev$Project, ev$N_Total),
       cex.main = 1.1, font.main = 2,
       xlab = "1 - Specificity (FPR)", ylab = "Sensitivity (TPR)")
  plot(ev$ROC_Prot, col = "#377eb8", lwd = 2.2, lty = 2, add = TRUE)
  plot(ev$ROC_Clin, col = "grey60", lwd = 2.0, lty = 3, add = TRUE)
  
  legend("bottomright", legend = c(
    sprintf("Combined: AUC %.3f", as.numeric(auc(ev$ROC_Comb))),
    sprintf("Proteins Only (%d DEPs): AUC %.3f", ev$Matched_Count, as.numeric(auc(ev$ROC_Prot))),
    sprintf("Clinical Baseline: AUC %.3f", as.numeric(auc(ev$ROC_Clin)))
  ), col = c("#d73027", "#377eb8", "grey60"), lwd = c(3.2, 2.2, 2.0), lty = c(1, 2, 3), bty = "n", cex = 0.8)
}
dev.off()

cat("Analysis complete. Supplementary Table 14b and ROC curves saved to:", output_dir, "\n")