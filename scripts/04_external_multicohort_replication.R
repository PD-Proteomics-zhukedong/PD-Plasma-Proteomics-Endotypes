# ==============================================================================
# Script: 04_external_multicohort_replication.R
# Purpose: Cross-platform effect size concordance and external cohort validation
#          across PPMI (P314), STTT 2026, and UK Biobank (Outputs: ST6, ST7, ST8)
# Project: Large-scale plasma proteomics in Parkinson's disease (Nature Aging)
# ==============================================================================

options(expressions = 5000)
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(arrow)
  library(readxl)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# Directory setup
data_dir      <- "data"
input_dea_dir <- "results/02_figure1_meta_dea"
output_dir    <- "results/04_external_multicohort_replication"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

find_file_smart <- function(pattern, default_path) {
  files <- list.files(data_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  files <- files[!grepl("^~\\$", basename(files))]
  if (length(files) > 0) return(files[1])
  return(default_path)
}

extract_numeric_safe <- function(x) {
  x_str <- trimws(as.character(x))
  num_str <- str_extract(x_str, "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?")
  suppressWarnings(as.numeric(num_str))
}

cat("\n=== Phase 4: External Multi-Cohort Replication (ST6, ST7, ST8) ===\n")

# ==============================================================================
# Step 1: Load Consensus Meta-DEPs (823 Proteins)
# ==============================================================================
cat("Loading consensus Meta-DEPs...\n")

meta_file_st5 <- file.path(input_dea_dir, "Supplementary_Table_5_Strict_MetaDEPs.csv")
meta_file_leg1 <- file.path(input_dea_dir, "Table_ST10_Final_823_Strict_Meta_DEPs.csv")
meta_file_leg2 <- file.path(input_dea_dir, "Table_3_Final_Strict_DEPs_StoufferMeta.csv")

if (file.exists(meta_file_st5)) {
  meta_file <- meta_file_st5
} else if (file.exists(meta_file_leg1)) {
  meta_file <- meta_file_leg1
} else if (file.exists(meta_file_leg2)) {
  meta_file <- meta_file_leg2
} else {
  stop("Required Meta-DEPs table not found. Please run Script 02 first.")
}

meta_df <- read.csv(meta_file, stringsAsFactors = FALSE)

# Harmonize column names depending on source table
if ("Consensus_Regulation" %in% colnames(meta_df)) {
  meta_df$Final_Status <- meta_df$Consensus_Regulation
}
if (!"meta_est" %in% colnames(meta_df) && "Meta_log2FC" %in% colnames(meta_df)) {
  meta_df$meta_est <- meta_df$Meta_log2FC
}
if (!"protein" %in% colnames(meta_df) && "Protein_Group" %in% colnames(meta_df)) {
  meta_df$protein <- meta_df$Protein_Group
}

meta_df <- meta_df %>%
  mutate(
    clean_uniprot = gsub("-.*|\\..*", "", sapply(strsplit(as.character(protein), ";"), `[`, 1)),
    Is_Meta_DEP   = grepl("Strictly Validated", Final_Status),
    My_Direction  = ifelse(meta_est > 0, "Up", "Down")
  )

if (!"Symbol" %in% colnames(meta_df) || any(is.na(meta_df$Symbol))) {
  mapped_symbols <- suppressMessages(suppressWarnings(
    mapIds(org.Hs.eg.db, keys = meta_df$clean_uniprot, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
  ))
  meta_df$Symbol <- as.character(mapped_symbols)
  meta_df$Symbol <- ifelse(is.na(meta_df$Symbol) | meta_df$Symbol == "", meta_df$clean_uniprot, meta_df$Symbol)
}

meta_df$Clean_Symbol <- toupper(gsub("[-_.]", "", meta_df$Symbol))
meta_deps_only <- meta_df %>% filter(Is_Meta_DEP)
cat(sprintf("Loaded %d consensus Meta-DEPs for external replication benchmarking.\n", nrow(meta_deps_only)))

# ==============================================================================
# Step 2: PPMI Project 314 Cross-Platform Replication (Supplementary Table 6)
# ==============================================================================
cat("\nAnalyzing PPMI Project 314 (Olink Explore HT, N = 868)...\n")

p314_parquet <- find_file_smart("ppmi_proj314.*\\.parquet$", file.path(data_dir, "ppmi_proj314_plasma_screened_extended_npx.parquet"))
clin_file    <- find_file_smart("PPMI_Curated_Data_Cut.*\\.xlsx$", file.path(data_dir, "PPMI_Curated_Data_Cut_Public.xlsx"))

if (!file.exists(p314_parquet) || !file.exists(clin_file)) {
  stop("PPMI Project 314 files not found in data directory.")
}

ppmi_clin <- readxl::read_excel(clin_file)
colnames(ppmi_clin) <- make.names(colnames(ppmi_clin))

col_age <- grep("^age", colnames(ppmi_clin), ignore.case = TRUE, value = TRUE)[1]
col_sex <- grep("^sex|^gender", colnames(ppmi_clin), ignore.case = TRUE, value = TRUE)[1]

clin_clean <- ppmi_clin %>%
  mutate(
    PATNO = as.character(PATNO),
    EVENT_ID = toupper(as.character(EVENT_ID)),
    COHORT = as.character(COHORT),
    Age = as.numeric(.data[[col_age]]),
    Sex = ifelse(.data[[col_sex]] %in% c("1", 1, "M", "Male", "male"), 1, 0)
  ) %>%
  filter(EVENT_ID %in% c("BL", "SC", "V00") & COHORT %in% c("1", "2")) %>%
  mutate(Group = ifelse(COHORT == "1", "PD", "HC")) %>%
  distinct(PATNO, .keep_all = TRUE)

p314_raw <- read_parquet(p314_parquet)
p314_clean <- p314_raw %>%
  mutate(
    PATNO = gsub("^PPMI-", "", as.character(PATNO)),
    UniProt = gsub("-.*|\\..*", "", as.character(UniProt)),
    NPX = as.numeric(NPX)
  ) %>%
  filter(PATNO %in% clin_clean$PATNO & !is.na(NPX))

common_prots <- intersect(unique(p314_clean$UniProt), meta_deps_only$clean_uniprot)

p314_results <- map_df(common_prots, function(u) {
  sub_data <- p314_clean %>% filter(UniProt == u) %>% inner_join(clin_clean, by = "PATNO")
  if (nrow(sub_data) < 50) return(NULL)
  fit <- lm(NPX ~ Group + Age + Sex, data = sub_data)
  tidy_fit <- broom::tidy(fit) %>% filter(term == "GroupPD")
  if (nrow(tidy_fit) == 0) return(NULL)
  data.frame(
    clean_uniprot = u,
    PPMI_Log2FC   = tidy_fit$estimate,
    PPMI_P_Val    = tidy_fit$p.value
  )
})

ppmi_matched <- meta_deps_only %>%
  inner_join(p314_results, by = "clean_uniprot") %>%
  mutate(
    Replicated_P05 = PPMI_P_Val < 0.05,
    Direction_Match = case_when(
      (PPMI_Log2FC > 0 & meta_est > 0) ~ "Consistent Up",
      (PPMI_Log2FC < 0 & meta_est < 0) ~ "Consistent Down",
      TRUE ~ "Opposite Direction"
    )
  ) %>%
  dplyr::select(
    Protein_Group = protein,
    UNIPROT_ID = clean_uniprot,
    Gene_Symbol = Symbol,
    Discovery_Meta_log2FC = meta_est,
    Discovery_Regulation = Final_Status,
    PPMI_P314_Log2FC = PPMI_Log2FC,
    PPMI_P314_Pval = PPMI_P_Val,
    Replicated_Nominal_P05 = Replicated_P05,
    Directional_Concordance = Direction_Match
  )

# Export Supplementary Table 6
write.csv(ppmi_matched, file.path(output_dir, "Supplementary_Table_6_PPMI_P314_Concordance.csv"), row.names = FALSE)

ppmi_n_matched <- nrow(ppmi_matched)
ppmi_n_sig     <- sum(ppmi_matched$Replicated_Nominal_P05)
ppmi_n_concord <- sum(ppmi_matched$Replicated_Nominal_P05 & ppmi_matched$Directional_Concordance != "Opposite Direction")
cat(sprintf("  - ST6 Exported: %d matched proteins, %d reaching P < 0.05 (%d concordant).\n", 
            ppmi_n_matched, ppmi_n_sig, ppmi_n_concord))

# ==============================================================================
# Step 3: STTT 2026 East Asian Cohort Validation (Supplementary Table 7)
# ==============================================================================
cat("\nAnalyzing STTT 2026 East Asian cohort (N = 226)...\n")

npx_file  <- find_file_smart("NPX data.*\\.xlsx$", file.path(data_dir, "NPX data.xlsx"))
prot_file <- find_file_smart("Protein list.*\\.xlsx$", file.path(data_dir, "Protein list.xlsx"))

if (!file.exists(npx_file) || !file.exists(prot_file)) {
  stop("STTT 2026 dataset files not found in data directory.")
}

prot_anno <- read_excel(prot_file)
col_sym <- grep("gene|symbol|assay", colnames(prot_anno), ignore.case = TRUE, value = TRUE)[1]

sttt_anno <- prot_anno %>%
  mutate(Clean_Symbol = toupper(gsub("[-_.]", "", .data[[col_sym]]))) %>%
  distinct(Clean_Symbol, .keep_all = TRUE)

npx_raw <- read_excel(npx_file)
protein_ids <- as.character(npx_raw[[1]])
sample_ids  <- colnames(npx_raw)[2:ncol(npx_raw)]
npx_mat     <- as.matrix(npx_raw[, 2:ncol(npx_raw)])
rownames(npx_mat) <- protein_ids

sample_groups <- ifelse(grepl("Ctrl|Control|HC|^C", sample_ids, ignore.case = TRUE), "Ctrl", "PD")

sttt_dea_list <- list()
for (i in 1:nrow(npx_mat)) {
  p_name <- rownames(npx_mat)[i]
  p_vals <- as.numeric(npx_mat[i, ])
  v_ctrl <- p_vals[sample_groups == "Ctrl"]
  v_pd   <- p_vals[sample_groups == "PD"]
  if (sum(!is.na(v_ctrl)) >= 5 && sum(!is.na(v_pd)) >= 5) {
    t_res <- t.test(v_pd, v_ctrl)
    sttt_dea_list[[p_name]] <- data.frame(
      Assay       = p_name,
      STTT_Log2FC = mean(v_pd, na.rm = TRUE) - mean(v_ctrl, na.rm = TRUE),
      STTT_P_Val  = t_res$p.value
    )
  }
}

sttt_dea_df <- bind_rows(sttt_dea_list) %>%
  mutate(Clean_Symbol = toupper(gsub("[-_.]", "", Assay)))

sttt_matched <- meta_deps_only %>%
  inner_join(sttt_dea_df, by = "Clean_Symbol") %>%
  mutate(
    Direction_Match = case_when(
      (STTT_Log2FC > 0 & meta_est > 0) ~ "Consistent Up",
      (STTT_Log2FC < 0 & meta_est < 0) ~ "Consistent Down",
      TRUE ~ "Opposite Direction"
    )
  ) %>%
  dplyr::select(
    Protein_Group = protein,
    UNIPROT_ID = clean_uniprot,
    Gene_Symbol = Symbol,
    Discovery_Meta_log2FC = meta_est,
    Discovery_Regulation = Final_Status,
    STTT_2026_Log2FC = STTT_Log2FC,
    STTT_2026_Pval = STTT_P_Val,
    Directional_Concordance = Direction_Match
  )

# Export Supplementary Table 7
write.csv(sttt_matched, file.path(output_dir, "Supplementary_Table_7_STTT_2026_Validation.csv"), row.names = FALSE)

sttt_n_matched <- nrow(sttt_matched)
sttt_n_concord <- sum(sttt_matched$Directional_Concordance != "Opposite Direction")
cat(sprintf("  - ST7 Exported: %d overlapping inflammation proteins (%d directionally concordant).\n", 
            sttt_n_matched, sttt_n_concord))

# ==============================================================================
# Step 4: UK Biobank Prospective Pre-diagnostic Risk (Supplementary Table 8)
# ==============================================================================
cat("\nAnalyzing UK Biobank prospective incident PD cohort (N > 49,000)...\n")

ukb_excel <- find_file_smart("43587_2025_818_MOESM3_ESM.*\\.xlsx$", file.path(data_dir, "43587_2025_818_MOESM3_ESM.xlsx"))

if (!file.exists(ukb_excel)) {
  stop("UK Biobank dataset file not found in data directory.")
}

raw_ukb_st4 <- read_excel(ukb_excel, sheet = "ST4", col_names = FALSE)
start_row <- which(apply(raw_ukb_st4[, 1:min(3, ncol(raw_ukb_st4))], 1, function(r) any(grepl("^[A-Za-z0-9-]{2,12}$", trimws(r)))))[1]
if (is.na(start_row)) start_row <- 4

ukb_st4_data <- raw_ukb_st4[start_row:nrow(raw_ukb_st4), ] %>%
  filter(!is.na(`...1`) & trimws(`...1`) != "") %>%
  mutate(
    Protein      = trimws(as.character(`...1`)),
    Clean_Symbol = toupper(gsub("[-_.]", "", Protein)),
    UKB_HR       = extract_numeric_safe(`...2`),
    UKB_Log2HR   = log2(UKB_HR),
    UKB_P_Val    = extract_numeric_safe(`...3`)
  ) %>%
  filter(!is.na(UKB_Log2HR) & is.finite(UKB_Log2HR))

ukb_matched <- meta_deps_only %>%
  inner_join(ukb_st4_data, by = "Clean_Symbol") %>%
  mutate(
    Is_Significant = !is.na(UKB_P_Val) & UKB_P_Val < 0.05,
    Direction_Match = case_when(
      (UKB_Log2HR > 0 & meta_est > 0) ~ "Consistent Up",
      (UKB_Log2HR < 0 & meta_est < 0) ~ "Consistent Down",
      TRUE ~ "Opposite Direction"
    )
  ) %>%
  dplyr::select(
    Protein_Group = protein,
    UNIPROT_ID = clean_uniprot,
    Gene_Symbol = Symbol,
    Discovery_Meta_log2FC = meta_est,
    Discovery_Regulation = Final_Status,
    UKB_Hazard_Ratio = UKB_HR,
    UKB_Log2_HR = UKB_Log2HR,
    UKB_Cox_Pval = UKB_P_Val,
    Significant_Incident_Risk = Is_Significant,
    Directional_Concordance = Direction_Match
  )

# Export Supplementary Table 8
write.csv(ukb_matched, file.path(output_dir, "Supplementary_Table_8_UK_Biobank_Prospective_Risk.csv"), row.names = FALSE)

ukb_n_matched <- nrow(ukb_matched)
ukb_n_sig     <- sum(ukb_matched$Significant_Incident_Risk)
ukb_n_concord <- sum(ukb_matched$Significant_Incident_Risk & ukb_matched$Directional_Concordance != "Opposite Direction")
cat(sprintf("  - ST8 Exported: %d overlapping proteins, %d significant incident risk predictors (%d concordant, %.1f%%).\n", 
            ukb_n_matched, ukb_n_sig, ukb_n_concord, (ukb_n_concord / ukb_n_sig) * 100))

cat("\nPhase 4 replication analysis complete. Supplementary Tables ST6, ST7, and ST8 successfully generated in:", output_dir, "\n")