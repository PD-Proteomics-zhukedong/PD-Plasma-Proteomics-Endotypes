# ==============================================================================
# Script: 05_gtex_tissue_enrichment.R
# Purpose: GTEx v8 multi-tissue hypergeometric enrichment analysis,
#          export of Supplementary Table 11, and generation of Supplementary Figure 3
# Project: Large-scale plasma proteomics in Parkinson's disease (Nature Aging)
# ==============================================================================

options(expressions = 5000)
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(scales)
})

# Cross-platform PDF device fallback
pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"

# Directory setup
data_dir       <- "data"
input_dea_dir  <- "results/02_figure1_meta_dea"
output_dir     <- "results/05_gtex_tissue_enrichment"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

find_file_smart <- function(pattern, default_path) {
  files <- list.files(data_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  files <- files[!grepl("^~\\$", basename(files))]
  if (length(files) > 0) return(files[1])
  return(default_path)
}

cat("=== Phase 5: GTEx Tissue-of-Origin Profiling (ST11 & Supp Fig 3) ===\n")

# ==============================================================================
# Step 1: Load Input Data and Quantified Plasma Universe
# ==============================================================================
cat("Loading consensus Meta-DEPs and plasma proteome background universe...\n")

meta_st5_file <- file.path(input_dea_dir, "Supplementary_Table_5_Strict_MetaDEPs.csv")
meta_leg_file <- file.path(input_dea_dir, "Table_ST10_Final_823_Strict_Meta_DEPs.csv")

if (file.exists(meta_st5_file)) {
  meta_deps_file <- meta_st5_file
} else if (file.exists(meta_leg_file)) {
  meta_deps_file <- meta_leg_file
} else {
  stop("Required Meta-DEPs table not found. Please run Script 02 first.")
}

combat_file <- file.path("results/01_preprocessed_data", "3_ComBat_Corrected_Matrix.tsv")
if (!file.exists(combat_file)) {
  stop("ComBat matrix not found: results/01_preprocessed_data/3_ComBat_Corrected_Matrix.tsv")
}

meta_final <- read.csv(meta_deps_file, stringsAsFactors = FALSE)

# Harmonize column names
if ("Consensus_Regulation" %in% colnames(meta_final)) meta_final$Final_Status <- meta_final$Consensus_Regulation
if (!"protein" %in% colnames(meta_final) && "Protein_Group" %in% colnames(meta_final)) meta_final$protein <- meta_final$Protein_Group

gold_deps_df <- meta_final %>% 
  filter(grepl("Strictly Validated", Final_Status)) %>%
  mutate(
    Protein_ID = gsub("-.*|\\..*", "", sapply(strsplit(as.character(protein), ";"), `[`, 1)),
    Direction  = ifelse(grepl("Up", Final_Status), "UP", "DOWN")
  )

# Map protein UniProt IDs to Gene Symbols
dep_symbols <- suppressMessages(suppressWarnings(
  bitr(gold_deps_df$Protein_ID, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
))
gold_deps_df <- gold_deps_df %>% inner_join(dep_symbols, by = c("Protein_ID" = "UNIPROT"))

# Define reliably quantified plasma proteome background universe (3,946 proteins)
bg_data <- fread(combat_file, data.table = FALSE)
bg_ids  <- unique(gsub("-.*|\\..*", "", sapply(strsplit(as.character(bg_data[[1]]), ";"), `[`, 1)))
bg_symbols <- suppressMessages(suppressWarnings(
  bitr(bg_ids, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
))$SYMBOL %>% unique()

cat(sprintf("Loaded %d mapped Meta-DEPs constrained against %d tested plasma background universe genes.\n", 
            nrow(gold_deps_df), length(bg_symbols)))

# ==============================================================================
# Step 2: Load GTEx Expression Matrix and Define Tissue Signatures
# ==============================================================================
cat("\nLoading GTEx v8 median TPM matrix and defining tissue-specific signatures...\n")

gtex_local_file <- find_file_smart("GTEx_Analysis_2017-06-05_v8.*gene_median_tpm.*", 
                                   file.path(data_dir, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"))

if (!file.exists(gtex_local_file)) {
  cat("Local GTEx file not found. Downloading reference matrix...\n")
  gtex_url <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
  download.file(gtex_url, destfile = file.path(data_dir, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"), mode = "wb")
  gtex_local_file <- file.path(data_dir, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")
}

gtex_raw <- fread(gtex_local_file, skip = 2, data.table = FALSE)
colnames(gtex_raw)[1:2] <- c("Name", "Description")

gtex_bg <- gtex_raw %>% 
  filter(Description %in% bg_symbols) %>%
  distinct(Description, .keep_all = TRUE)

rownames(gtex_bg) <- gtex_bg$Description
gtex_mat <- as.matrix(gtex_bg[, 3:ncol(gtex_bg)])
mode(gtex_mat) <- "numeric"

# Define tissue-specific signatures: median TPM >= 5 and >= 2.5x median of all other 53 tissues
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
cat(sprintf("Extracted tissue-specific gene signatures across %d GTEx human tissues.\n", length(tissue_specific_list)))

# ==============================================================================
# Step 3: Multi-Tissue Hypergeometric Enrichment Analysis (Supplementary Table 11)
# ==============================================================================
cat("\nPerforming background-constrained hypergeometric over-representation tests across 54 tissues...\n")

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
    
    # Two-sided Fisher's exact test (as documented in Methods & Legends)
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

# Export Supplementary Table 11 (Strict ST11 alignment)
st11_file <- file.path(output_dir, "Supplementary_Table_11_GTEx_v8_54_Tissue_Enrichment.csv")
write.csv(enrich_gtex_all, st11_file, row.names = FALSE)
cat("Supplementary Table 11 successfully exported to:", st11_file, "\n")

# ==============================================================================
# Step 4: Generate Supplementary Figure 3 (Target Tissue-of-Origin Lollipop Plot)
# ==============================================================================
cat("\nGenerating Supplementary Figure 3 (Exact Publication Format)...\n")

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

p_supp_fig3 <- ggplot(plot_main_df, aes(x = Clean_Tissue, y = plot_val, color = Direction)) +
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
    title = "Target Tissue-of-Origin Deconvolution of Circulating PD Meta-DEPs",
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

# Save Supplementary Figure 3 in PDF and 300 dpi PNG
ggsave(file.path(output_dir, "Supplementary_Figure_3_GTEx_Tissue_Enrichment.pdf"), 
       p_supp_fig3, width = 9.5, height = 6.5, device = pdf_device)
ggsave(file.path(output_dir, "Supplementary_Figure_3_GTEx_Tissue_Enrichment.png"), 
       p_supp_fig3, width = 9.5, height = 6.5, dpi = 300)

cat("Phase 5 analysis complete. Supplementary Figure 3 and Table ST11 successfully saved to:", output_dir, "\n")