# ==============================================================================
# Script: 06_gtex_tissue_enrichment.R
# Purpose: GTEx v8 multi-tissue enrichment analysis and Supplementary Figure S4
# Project: Large-scale plasma proteomics identifies molecularly distinct PD endotypes
# ==============================================================================

options(expressions = 5000)
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(pheatmap)
  library(ppcor)
  library(scales)
  library(grid)
  library(patchwork)
})

# Directory setup
data_dir       <- "data"
input_dea_dir  <- "results/02_figure1_meta_dea"
output_dir     <- "results/06_gtex_tissue_enrichment"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("=== Phase 6: GTEx Tissue-of-Origin Profiling and Supplementary Figure S4 ===\n")

# ==============================================================================
# Step 1: Load Input Data and Clinical Metadata
# ==============================================================================
cat("Loading input matrices and metadata...\n")

meta_deps_file <- file.path(input_dea_dir, "Table_ST10_Final_823_Strict_Meta_DEPs.csv")
combat_file    <- file.path("results/01_preprocessed_data", "3_ComBat_Corrected_Matrix.tsv")
clinical_file  <- file.path(data_dir, "metadata_clinical_n1119.csv")

if (!file.exists(meta_deps_file) || !file.exists(combat_file)) {
  stop("Required DEA or ComBat files not found. Please run scripts 01 and 02 first.")
}

meta_final <- read.csv(meta_deps_file, stringsAsFactors = FALSE)
gold_deps_df <- meta_final %>% 
  filter(Final_Status %in% c("Strictly Validated Up", "Strictly Validated Down")) %>%
  mutate(Protein_ID = gsub("-.*|\\..*", "", protein),
         Direction = ifelse(Final_Status == "Strictly Validated Up", "UP", "DOWN"))

# Map protein UniProt IDs to Gene Symbols
dep_symbols <- suppressMessages(suppressWarnings(
  bitr(gold_deps_df$Protein_ID, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
))
gold_deps_df <- gold_deps_df %>% inner_join(dep_symbols, by = c("Protein_ID" = "UNIPROT"))

# Load quantified proteins as the tested background universe (3,940 unique genes)
bg_data <- fread(combat_file, data.table = FALSE)
bg_ids  <- unique(gsub("-.*|\\..*", "", sapply(strsplit(as.character(bg_data[[1]]), ";"), `[`, 1)))
bg_symbols <- suppressMessages(suppressWarnings(
  bitr(bg_ids, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
))$SYMBOL %>% unique()

# Filter clinical metadata for sporadic PD cohort (n = 669)
clinical_clean <- read.csv(clinical_file, stringsAsFactors = FALSE) %>%
  mutate(group = ifelse(group %in% c("CTR", "Control", "HC"), "HC", "PD"))

clinical_pd <- clinical_clean %>% 
  filter(group == "PD" & !is.na(Age) & !is.na(Sex)) %>%
  mutate(
    LEDD = replace_na(as.numeric(LEDD), 0),
    Duration = replace_na(as.numeric(Disease_Duration), 0),
    CV_Metabolic_Cat = as.numeric(as.factor(CV_Metabolic_Score)),
    hepatic_disease = replace_na(as.numeric(hepatic_disease), 0),
    kidney.function = replace_na(as.numeric(kidney.function), 0)
  )

cat(sprintf("Loaded %d mapped Meta-DEPs, %d background universe genes, and %d PD participants.\n", 
            nrow(gold_deps_df), length(bg_symbols), nrow(clinical_pd)))

# ==============================================================================
# Step 2: Load GTEx Expression Matrix and Define Tissue Signatures
# ==============================================================================
cat("Loading GTEx v8 median TPM matrix and extracting tissue-specific signatures...\n")

gtex_local_file <- file.path(data_dir, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")

if (!file.exists(gtex_local_file)) {
  cat("Downloading GTEx reference matrix...\n")
  gtex_url <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
  download.file(gtex_url, destfile = gtex_local_file, mode = "wb")
}

gtex_raw <- fread(gtex_local_file, skip = 2, data.table = FALSE)
colnames(gtex_raw)[1:2] <- c("Name", "Description")

gtex_bg <- gtex_raw %>% 
  filter(Description %in% bg_symbols) %>%
  distinct(Description, .keep_all = TRUE)

rownames(gtex_bg) <- gtex_bg$Description
gtex_mat <- as.matrix(gtex_bg[, 3:ncol(gtex_bg)])
mode(gtex_mat) <- "numeric"

# Define tissue-specific signatures (TPM >= 5 and >= 2.5x median of other tissues)
tissue_specific_list <- list()
for (t in colnames(gtex_mat)) {
  t_expr <- gtex_mat[, t]
  other_expr <- gtex_mat[, colnames(gtex_mat) != t]
  other_median <- apply(other_expr, 1, median)
  
  is_specific <- (t_expr >= 5) & (t_expr >= 2.5 * other_median)
  specific_genes <- names(which(is_specific))
  if (length(specific_genes) >= 3) {
    tissue_specific_list[[t]] <- specific_genes
  }
}

# ==============================================================================
# Step 3: Multi-Tissue Hypergeometric Enrichment Analysis (Table ST12)
# ==============================================================================
cat("Performing hypergeometric enrichment across 54 tissues...\n")

up_symbols   <- gold_deps_df %>% filter(Direction == "UP") %>% pull(SYMBOL) %>% unique()
down_symbols <- gold_deps_df %>% filter(Direction == "DOWN") %>% pull(SYMBOL) %>% unique()
N_bg_total   <- length(bg_symbols)

calc_gtex_enrichment <- function(query_genes, direction_label) {
  N_fg <- length(query_genes)
  res_list <- lapply(names(tissue_specific_list), function(t) {
    t_genes <- tissue_specific_list[[t]]
    k <- length(intersect(query_genes, t_genes))
    K <- length(intersect(bg_symbols, t_genes))
    if (k == 0 || K == 0) return(NULL)
    
    cont_mat <- matrix(c(k, K - k, N_fg - k, N_bg_total - K - (N_fg - k)), nrow = 2)
    if (any(cont_mat < 0)) return(NULL)
    
    p_val <- fisher.test(cont_mat, alternative = "greater")$p.value
    fe    <- (k / N_fg) / (K / N_bg_total)
    data.frame(GTEx_Tissue = t, Count = k, Bg_Count = K, Fold_Enrichment = fe, P_Value = p_val, Direction = direction_label)
  })
  bind_rows(res_list) %>% mutate(FDR = p.adjust(P_Value, method = "BH")) %>% arrange(P_Value)
}

enrich_gtex_up   <- calc_gtex_enrichment(up_symbols, "UP")
enrich_gtex_down <- calc_gtex_enrichment(down_symbols, "DOWN")
enrich_gtex_all  <- bind_rows(enrich_gtex_up, enrich_gtex_down)

clean_gtex_name <- function(t_str) {
  t_str <- gsub("Brain - ", "Brain: ", t_str)
  t_str <- gsub("Skin - ", "Skin: ", t_str)
  t_str <- gsub("Heart - ", "Heart: ", t_str)
  t_str <- gsub("Artery - ", "Artery: ", t_str)
  t_str <- gsub("_", " ", t_str)
  return(t_str)
}
enrich_gtex_all$Clean_Tissue <- clean_gtex_name(enrich_gtex_all$GTEx_Tissue)

st12_file <- file.path(output_dir, "Supplementary_Table_12_GTEx_v8_54_Tissue_Enrichment.csv")
write.csv(enrich_gtex_all, st12_file, row.names = FALSE)
cat("Supplementary Table 12 exported.\n")

# ==============================================================================
# Step 4: Generate Supplementary Figure S4A (Core Tissues Lollipop Plot)
# ==============================================================================
cat("Generating Supplementary Figure S4A...\n")

pd_core_tissues <- c(
  "Brain: Substantia nigra",
  "Brain: Putamen (basal ganglia)",
  "Brain: Caudate (basal ganglia)",
  "Brain: Nucleus accumbens (basal ganglia)",
  "Brain: Frontal Cortex (BA9)",
  "Brain: Cortex",
  "Brain: Spinal cord (cervical c-1)",
  "Muscle - Skeletal",
  "Heart: Left Ventricle",
  "Artery: Aorta",
  "Whole Blood",
  "Spleen",
  "Liver",
  "Kidney - Cortex",
  "Skin: Sun Exposed (Lower leg)"
)

plot_main_df <- enrich_gtex_all %>%
  filter(Clean_Tissue %in% pd_core_tissues) %>%
  group_by(Clean_Tissue) %>%
  slice_min(P_Value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    logP = -log10(P_Value),
    plot_val = ifelse(Direction == "DOWN", -logP, logP),
    plot_val = pmax(pmin(plot_val, 4.5), -4.5),
    sig_star = case_when(P_Value < 0.001 ~ "***", P_Value < 0.01 ~ "**", P_Value < 0.05 ~ "*", TRUE ~ "")
  )

order_main <- plot_main_df %>% arrange(plot_val) %>% pull(Clean_Tissue)
plot_main_df$Clean_Tissue <- factor(plot_main_df$Clean_Tissue, levels = order_main)
n_main <- nlevels(plot_main_df$Clean_Tissue)

p_supp_s4a <- ggplot(plot_main_df, aes(x = Clean_Tissue, y = plot_val, color = Direction)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  geom_hline(yintercept = c(-1.301, 1.301), color = "grey65", linetype = "dashed", linewidth = 0.5) +
  geom_segment(aes(x = Clean_Tissue, xend = Clean_Tissue, y = 0, yend = plot_val), linewidth = 1.15, show.legend = FALSE) +
  geom_point(aes(size = Count, fill = Direction), shape = 21, stroke = 1.15, color = "white") +
  geom_text(aes(label = sig_star, y = plot_val + sign(plot_val) * 0.35), size = 5.0, color = "black", vjust = 0.75) +
  scale_color_manual(values = c("UP" = "#d73027", "DOWN" = "#313695"), labels = c("DOWN" = "DOWN", "UP" = "UP")) +
  scale_fill_manual(values = c("UP" = "#d73027", "DOWN" = "#313695"), labels = c("DOWN" = "DOWN", "UP" = "UP")) +
  scale_size_continuous(range = c(4.0, 9.0), name = "Number of DEPs", breaks = c(5, 15, 30, 60)) +
  coord_flip(clip = "off") +
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4), labels = c("4", "2", "0", "2", "4"), limits = c(-5.5, 5.5)) +
  annotate("text", x = n_main + 0.85, y = 1.6, label = "Upregulated in PD →", color = "#d73027", fontface = "bold", size = 3.6, hjust = 0) +
  annotate("text", x = n_main + 0.85, y = -1.6, label = "← Downregulated in PD", color = "#313695", fontface = "bold", size = 3.6, hjust = 1) +
  labs(
    title = "Target Tissue of Origin Profiling of Circulating DEPs (GTEx v8 Benchmark)",
    x = NULL, y = expression(bold(paste("Statistical Significance (-Log"[10], " P-value)")))) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12.5, margin = ggplot2::margin(b = 10)),
    axis.text.y = element_text(face = "bold", color = "black", size = 10.5),
    axis.text.x = element_text(color = "black", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linetype = "dotted"),
    legend.position = "right",
    plot.margin = ggplot2::margin(15, 15, 15, 15)
  )

# ==============================================================================
# Step 5: Partial Correlation Analysis and Supplementary Figure S4B Heatmap
# ==============================================================================
cat("Calculating partial correlations and generating Supplementary Figure S4B...\n")

raw_ids <- make.unique(gsub("-.*|\\..*", "", sapply(strsplit(as.character(bg_data[[1]]), ";"), `[`, 1)))
mapped_syms <- suppressMessages(suppressWarnings(
  mapIds(org.Hs.eg.db, keys = raw_ids, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
))
mapped_syms[is.na(mapped_syms)] <- raw_ids[is.na(mapped_syms)]

prot_expr_df <- as.data.frame(bg_data[, 7:ncol(bg_data)])
prot_expr_df$SYMBOL <- as.character(mapped_syms)

prot_mat_sym <- prot_expr_df %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), ~mean(., na.rm = TRUE))) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

common_samples <- intersect(colnames(prot_mat_sym), clinical_pd$sample)
prot_final_sym <- prot_mat_sym[, common_samples]
clinical_final <- clinical_pd[match(common_samples, clinical_pd$sample), ]
prot_z_mat     <- t(scale(t(prot_final_sym)))

gtex_tissues_dict <- list(
  "Brain: Substantia nigra"        = "Brain - Substantia nigra",
  "Brain: Putamen (basal ganglia)" = "Brain - Putamen (basal ganglia)",
  "Brain: Caudate (basal ganglia)" = "Brain - Caudate (basal ganglia)",
  "Brain: Cortex"                  = "Brain - Cortex",
  "Brain: Frontal Cortex (BA9)"    = "Brain - Frontal Cortex (BA9)",
  "Muscle – Skeletal"              = "Muscle - Skeletal",
  "Heart: Left Ventricle"          = "Heart - Left Ventricle",
  "Liver"                          = "Liver",
  "Whole Blood"                    = "Whole Blood",
  "Kidney – Cortex"                = "Kidney - Cortex"
)

tissue_protein_df_list <- list()
for (label in names(gtex_tissues_dict)) {
  t_name <- gtex_tissues_dict[[label]]
  t_col  <- which(colnames(gtex_mat) == t_name)
  if (length(t_col) > 0) {
    t_expr <- gtex_mat[, t_col]
    other_med <- apply(gtex_mat[, -t_col], 1, median)
    spec_genes <- names(which(t_expr >= 5 & t_expr >= 2.0 * other_med))
    matched_deps <- gold_deps_df %>% 
      filter(SYMBOL %in% spec_genes & SYMBOL %in% rownames(prot_z_mat)) %>%
      mutate(Tissue_Origin = label) %>%
      dplyr::select(SYMBOL, Direction, Tissue_Origin)
    if (nrow(matched_deps) > 0) {
      tissue_protein_df_list[[label]] <- matched_deps
    }
  }
}
tissue_protein_all <- bind_rows(tissue_protein_df_list) %>% distinct(SYMBOL, .keep_all = TRUE)

traits <- c("UPDRS_III", "NMSS", "MMSE", "HAMA", "HAMD", "PDSS", "UPSIT")
actual_traits <- intersect(traits, colnames(clinical_final))

trait_labels_display <- c(
  "UPDRS_III" = "Motor (UPDRS-III)", 
  "NMSS"      = "Total Non-motor (NMSS)", 
  "MMSE"      = "Cognition (MMSE)", 
  "HAMA"      = "Anxiety (HAMA)", 
  "HAMD"      = "Depression (HAMD)", 
  "PDSS"      = "Sleep (PDSS)", 
  "UPSIT"     = "Olfaction (UPSIT)"
)

r_all_mat <- matrix(NA_real_, nrow = nrow(tissue_protein_all), ncol = length(actual_traits),
                    dimnames = list(tissue_protein_all$SYMBOL, actual_traits))
p_all_mat <- matrix(NA_real_, nrow = nrow(tissue_protein_all), ncol = length(actual_traits),
                    dimnames = list(tissue_protein_all$SYMBOL, actual_traits))

for (gene in tissue_protein_all$SYMBOL) {
  for (tr in actual_traits) {
    temp_df <- data.frame(
      score     = as.numeric(prot_z_mat[gene, ]),
      trait_val = as.numeric(clinical_final[[tr]]),
      Age       = clinical_final$Age,
      Sex       = clinical_final$Sex,
      BMI       = clinical_final$BMI,
      hepatic   = clinical_final$hepatic_disease,
      kidney    = clinical_final$kidney.function,
      CV_Met    = clinical_final$CV_Metabolic_Cat,
      LEDD      = clinical_final$LEDD,
      Duration  = clinical_final$Duration
    ) %>% drop_na()
    
    if (nrow(temp_df) >= 30) {
      covars <- temp_df[, c("Age", "Sex", "BMI", "hepatic", "kidney", "CV_Met", "LEDD", "Duration")]
      pcor_res <- tryCatch({
        pcor.test(temp_df$score, temp_df$trait_val, covars, method = "spearman")
      }, error = function(e) NULL)
      
      if (!is.null(pcor_res)) {
        r_all_mat[gene, tr] <- pcor_res$estimate
        p_all_mat[gene, tr] <- pcor_res$p.value
      }
    }
  }
}

colnames(r_all_mat) <- trait_labels_display[colnames(r_all_mat)]
colnames(p_all_mat) <- trait_labels_display[colnames(p_all_mat)]

# Select top 2 representative proteins per tissue
tissue_protein_ranked <- tissue_protein_all %>%
  rowwise() %>%
  mutate(
    Min_Clinical_P = min(p_all_mat[SYMBOL, ], na.rm = TRUE),
    Max_Clinical_R = max(abs(r_all_mat[SYMBOL, ]), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  group_by(Tissue_Origin) %>%
  arrange(Min_Clinical_P, desc(Max_Clinical_R)) %>%
  slice_head(n = 2) %>%
  ungroup()

r_sub_mat <- r_all_mat[tissue_protein_ranked$SYMBOL, ]
p_sub_mat <- p_all_mat[tissue_protein_ranked$SYMBOL, ]

sig_stars_sub <- matrix(
  ifelse(is.na(p_sub_mat), "",
         ifelse(p_sub_mat < 0.001, "***",
                ifelse(p_sub_mat < 0.01, "**",
                       ifelse(p_sub_mat < 0.05, "*", "")))),
  nrow = nrow(p_sub_mat), dimnames = dimnames(p_sub_mat)
)

row_ann_sub <- data.frame(
  PD_Regulation = ifelse(tissue_protein_ranked$Direction == "UP", "Upregulated (Leakage)", "Downregulated (Depletion)"),
  Tissue_Origin = tissue_protein_ranked$Tissue_Origin,
  row.names     = tissue_protein_ranked$SYMBOL
)

ann_colors_sub <- list(
  PD_Regulation = c("Upregulated (Leakage)" = "#d73027", "Downregulated (Depletion)" = "#313695"),
  Tissue_Origin = setNames(
    scales::hue_pal()(length(unique(tissue_protein_ranked$Tissue_Origin))),
    unique(tissue_protein_ranked$Tissue_Origin)
  )
)

limit_val <- 0.40

pdf_s4b_path <- file.path(output_dir, "Supplementary_Figure_S4B_Heatmap.pdf")
pdf(pdf_s4b_path, width = 11.0, height = 7.5)
pheatmap(
  r_sub_mat,
  annotation_row    = row_ann_sub,
  annotation_colors = ann_colors_sub,
  display_numbers   = sig_stars_sub,
  number_color      = "black",
  fontsize_number   = 11.0,
  color             = colorRampPalette(c("#313695", "#FFFFFF", "#A50026"))(100),
  breaks            = seq(-limit_val, limit_val, length.out = 100),
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  main              = "Representative Organ-Specific Biomarkers Correlated with Clinical Phenotypes\n(8-Covariate Adjusted Partial Correlation: LEDD, Duration, Age, Sex, BMI, Hepatic, Renal, CV)",
  angle_col         = 45,
  border_color      = "grey88",
  fontsize_row      = 10.5,
  fontsize_col      = 10.5
)
dev.off()

ggsave(file.path(output_dir, "Supplementary_Figure_S4A_Lollipop.pdf"), p_supp_s4a, width = 9.5, height = 6.5, device = cairo_pdf)

cat("Analysis complete. Supplementary Figure S4 and Table ST12 saved to:", output_dir, "\n")