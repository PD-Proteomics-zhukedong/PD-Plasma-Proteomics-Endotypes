# ==============================================================================
# Script: 04_external_multicohort_replication.R
# Purpose: Cross-platform replication of Meta-DEPs across PPMI, STTT, and UK Biobank
# Project: Large-scale plasma proteomics identifies molecularly distinct PD endotypes
# ==============================================================================

options(expressions = 5000)
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(arrow)
  library(readxl)
  library(ggrepel)
  library(patchwork)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# Directory setup
data_dir      <- "data"
input_dea_dir <- "results/02_figure1_meta_dea"
output_dir    <- "results/04_external_validation_supp_fig_s3"
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

cat("=== Phase 4: External Multi-Cohort Replication and Supplementary Figure 3 ===\n")

# ==============================================================================
# Step 1: Load Consensus Meta-DEPs (823 Proteins)
# ==============================================================================
cat("Loading consensus Meta-DEPs...\n")

meta_file <- file.path(input_dea_dir, "Table_ST10_Final_823_Strict_Meta_DEPs.csv")
if (!file.exists(meta_file)) {
  meta_file <- file.path(input_dea_dir, "Table_3_Final_Strict_DEPs_StoufferMeta.csv")
}
if (!file.exists(meta_file)) stop("Required Meta-DEPs table not found. Please run Script 02 first.")

meta_df <- read.csv(meta_file, stringsAsFactors = FALSE) %>%
  mutate(
    clean_uniprot = gsub("-.*|\\..*", "", sapply(strsplit(as.character(protein), ";"), `[`, 1)),
    Is_Meta_DEP   = grepl("Strictly Validated", Final_Status),
    My_Direction  = ifelse(meta_est > 0, "Up", "Down")
  )

mapped_symbols <- suppressMessages(suppressWarnings(
  mapIds(org.Hs.eg.db, keys = meta_df$clean_uniprot, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
))
meta_df$Symbol <- as.character(mapped_symbols)
meta_df$Symbol <- ifelse(is.na(meta_df$Symbol) | meta_df$Symbol == "", meta_df$clean_uniprot, meta_df$Symbol)
meta_df$Clean_Symbol <- toupper(gsub("[-_.]", "", meta_df$Symbol))

meta_deps_only <- meta_df %>% filter(Is_Meta_DEP)

# ==============================================================================
# Step 2: PPMI Project 314 Cross-Platform Replication (Panel A & Table ST21)
# ==============================================================================
cat("Analyzing PPMI Project 314 (Olink Explore HT, N = 868)...\n")

p314_parquet <- find_file_smart("ppmi_proj314.*\\.parquet$", file.path(data_dir, "ppmi_proj314_plasma_screened_extended_npx.parquet"))
clin_file    <- find_file_smart("PPMI_Curated_Data_Cut.*\\.xlsx$", file.path(data_dir, "PPMI_Curated_Data_Cut_Public.xlsx"))

if (!file.exists(p314_parquet) || !file.exists(clin_file)) {
  stop("PPMI Project 314 files not found in data directory.")
}

ppmi_clin <- readxl::read_excel(clin_file)
colnames(ppmi_clin) <- make.names(colnames(ppmi_clin))

# Filter baseline visits
clin_clean <- ppmi_clin %>%
  mutate(
    PATNO = as.character(PATNO),
    EVENT_ID = toupper(as.character(EVENT_ID)),
    COHORT = as.character(COHORT),
    Age = as.numeric(AGE_AT_VISIT),
    Sex = ifelse(SEX %in% c("1", 1, "M", "Male", "male"), 1, 0)
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

# Calculate per-protein multivariable regression in PPMI (Group + Age + Sex)
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
  )

write.csv(ppmi_matched, file.path(output_dir, "Supplementary_Table_21_PPMI_P314_Concordance.csv"), row.names = FALSE)

# Dynamically compute statistical counts for annotation
ppmi_n_matched <- nrow(ppmi_matched)
ppmi_n_sig     <- sum(ppmi_matched$Replicated_P05)
ppmi_n_concord <- sum(ppmi_matched$Replicated_P05 & ppmi_matched$Direction_Match != "Opposite Direction")

ppmi_annotation_label <- sprintf(
  "PPMI Project 314 (N = 868, Olink Explore HT)\nMatched Meta-DEPs: %d\nReplicated in PPMI (P < 0.05): %d\nDirectional Concordance: %d",
  ppmi_n_matched, ppmi_n_sig, ppmi_n_concord
)

# Subset nominally significant proteins for visualization
plot_ppmi_df <- ppmi_matched %>% filter(Replicated_P05)

# Generate Supplementary Figure 3A
p_supp3a <- ggplot(plot_ppmi_df, aes(x = PPMI_Log2FC, y = meta_est)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  geom_point(aes(fill = Direction_Match), shape = 21, size = 4.0, color = "black", stroke = 0.8) +
  scale_fill_manual(values = c("Consistent Up" = "#d73027", "Consistent Down" = "#2b5c8f", "Opposite Direction" = "grey70"),
                    name = NULL) +
  geom_text_repel(data = filter(plot_ppmi_df, Direction_Match != "Opposite Direction"),
                  aes(label = Symbol), size = 3.6, fontface = "bold.italic", box.padding = 0.35, max.overlaps = 30) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
           label = ppmi_annotation_label, size = 3.6, fontface = "plain") +
  labs(title = "Cross-Platform Validation with PPMI Project 314",
       x = expression(bold(paste("PPMI P314 Effect Size (Log"[2], " FC)"))),
       y = expression(bold(paste("Discovery Meta Effect Size (Log"[2], " FC)")))) +
  scale_x_continuous(limits = c(-0.35, 0.45)) +
  scale_y_continuous(limits = c(-1.25, 1.0)) +
  theme_bw(base_size = 11.5) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

# ==============================================================================
# Step 3: STTT 2026 East Asian Cohort Validation (Panel B & Table ST20)
# ==============================================================================
cat("Analyzing STTT 2026 East Asian cohort (N = 226)...\n")

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
  )

write.csv(sttt_matched, file.path(output_dir, "Supplementary_Table_20_STTT_2026_Validation.csv"), row.names = FALSE)

# Dynamically compute statistical counts for annotation
sttt_n_matched <- nrow(sttt_matched)
sttt_n_concord <- sum(sttt_matched$Direction_Match != "Opposite Direction")

sttt_annotation_label <- sprintf(
  "STTT 2026 (East Asian Olink, N = 226)\nMatched Meta-DEPs: %d\nDirectional Concordance: %d",
  sttt_n_matched, sttt_n_concord
)

plot_sttt_df <- sttt_matched

# Generate Supplementary Figure 3B
p_supp3b <- ggplot(plot_sttt_df, aes(x = STTT_Log2FC, y = meta_est)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  geom_point(aes(fill = Direction_Match), shape = 21, size = 4.0, color = "black", stroke = 0.8) +
  scale_fill_manual(values = c("Consistent Up" = "#d73027", "Consistent Down" = "#2b5c8f", "Opposite Direction" = "grey70"),
                    name = NULL) +
  geom_text_repel(data = filter(plot_sttt_df, Direction_Match != "Opposite Direction"),
                  aes(label = Symbol), size = 3.6, fontface = "bold.italic", box.padding = 0.35, max.overlaps = 30) +
  annotate("text", x = -Inf, y = 0.15, hjust = -0.05, vjust = 0,
           label = sttt_annotation_label, size = 3.6, fontface = "plain") +
  labs(title = "Cross-Platform Validation with STTT 2026",
       x = expression(bold(paste("STTT 2026 Effect Size (Log"[2], " FC)"))),
       y = expression(bold(paste("Discovery Meta Effect Size (Log"[2], " FC)")))) +
  scale_x_continuous(limits = c(-0.7, 0.5)) +
  scale_y_continuous(limits = c(-1.25, 0.85)) +
  theme_bw(base_size = 11.5) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

# ==============================================================================
# Step 4: UK Biobank Prospective Pre-diagnostic Risk (Panel C & Table ST19)
# ==============================================================================
cat("Analyzing UK Biobank prospective incident PD cohort (N > 49,000)...\n")

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
  )

write.csv(ukb_matched, file.path(output_dir, "Supplementary_Table_19_UK_Biobank_Prospective_Risk.csv"), row.names = FALSE)

# Dynamically compute statistical counts for annotation
ukb_n_matched <- nrow(ukb_matched)
ukb_n_sig     <- sum(ukb_matched$Is_Significant)
ukb_n_concord <- sum(ukb_matched$Is_Significant & ukb_matched$Direction_Match != "Opposite Direction")

ukb_annotation_label <- sprintf(
  "UK Biobank (Pre-diagnostic PD)\nMatched Meta-DEPs: %d\nSignificant Prospective Predictors (P < 0.05): %d\nDirectional Concordance: %d",
  ukb_n_matched, ukb_n_sig, ukb_n_concord
)

plot_ukb_df <- ukb_matched %>% filter(Is_Significant)

# Generate Supplementary Figure 3C
p_supp3c <- ggplot(plot_ukb_df, aes(x = UKB_Log2HR, y = meta_est)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  geom_point(aes(fill = Direction_Match), shape = 21, size = 4.0, color = "black", stroke = 0.8) +
  scale_fill_manual(values = c("Consistent Up" = "#d73027", "Consistent Down" = "#2b5c8f", "Opposite Direction" = "grey70"),
                    name = NULL) +
  geom_text_repel(data = filter(plot_ukb_df, Direction_Match != "Opposite Direction"),
                  aes(label = Symbol), size = 3.2, fontface = "bold.italic", box.padding = 0.35, max.overlaps = 50) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
           label = ukb_annotation_label, size = 3.6, fontface = "plain") +
  labs(title = "Prospective Pre-diagnostic Risk Concordance (UK Biobank)",
       x = expression(bold(paste("UK Biobank Pre-diagnostic Risk (Log"[2], " Hazard Ratio)"))),
       y = expression(bold(paste("Discovery Meta Effect Size (Log"[2], " FC)")))) +
  scale_x_continuous(limits = c(-1.8, 0.8)) +
  scale_y_continuous(limits = c(-1.3, 1.0)) +
  theme_bw(base_size = 11.5) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

# ==============================================================================
# Step 5: Export Supplementary Figure 3 Multi-Panel Layout
# ==============================================================================
cat("Composing and exporting Supplementary Figure 3 (Panels A, B, C)...\n")

supp_fig_3 <- (p_supp3a | p_supp3b) / (p_supp3c | plot_spacer()) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "Supplementary_Figure_3.pdf"), supp_fig_3, width = 12.5, height = 11.0, device = cairo_pdf)
ggsave(file.path(output_dir, "Supplementary_Figure_3.png"), supp_fig_3, width = 12.5, height = 11.0, dpi = 300)

cat("Analysis complete. Supplementary Figure 3 and Tables ST19, ST20, ST21 saved to:", output_dir, "\n")