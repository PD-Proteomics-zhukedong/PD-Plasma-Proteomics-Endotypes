# ==============================================================================
# Script: 07_subtype_mechanisms_and_cellular_deconvolution.R
# Project: Large-scale plasma proteomics in Parkinson's disease (Nature Aging)
# Purpose: Subtype-specific pairwise Limma contrasts, Min-T bottleneck scoring,
#          PanglaoDB cellular deconvolution, 8-covariate partial correlations,
#          Figure 4, Supplementary Figure 5, and Tables ST16–ST20 generation
# ==============================================================================

options(expressions = 5000)
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(limma)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(ggpubr)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(GO.db)
  library(stringr)
  library(readxl)
  library(broom)
  library(UpSetR)
  library(ppcor)
  library(GSVA)
  library(patchwork)
  library(scales)
})

# Cross-platform PDF device fallback
pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"

# Directory setup
data_dir        <- "data"
input_dea_dir   <- "results/02_figure1_meta_dea"
input_clust_dir <- "results/06_figure3_molecular_endotypes"
output_dir      <- "results/07_subtype_mechanisms_and_deconvolution"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

find_file_smart <- function(pattern, default_path) {
  files <- list.files(data_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  files <- files[!grepl("^~\\$", basename(files))]
  if (length(files) > 0) return(files[1])
  return(default_path)
}

clean_num <- function(x) {
  if (is.null(x) || all(is.na(x))) return(NA_real_)
  as.numeric(gsub("[^0-9.-]", "", as.character(x)))
}

BASE_COVARIATES <- c("Age", "Sex", "BMI", "hepatic_disease", "kidney.function", "CV_Metabolic_Cat", "LEDD", "Duration")

cat("\n=== Phase 7: Subtype Mechanisms, Deconvolution, and Figure 4 Pipeline ===\n")

# ==============================================================================
# Step 1: Load Proteomics Matrices, Meta-DEPs, and Subtype Assignments
# ==============================================================================
cat("Loading preprocessed matrices, consensus Meta-DEPs, and subtype assignments...\n")

combat_file <- file.path("results/01_preprocessed_data", "3_ComBat_Corrected_Matrix.tsv")
meta_st5_file <- file.path(input_dea_dir, "Supplementary_Table_5_Strict_MetaDEPs.csv")
meta_leg_file <- file.path(input_dea_dir, "Table_3_Final_Strict_DEPs_StoufferMeta.csv")
meta_file     <- if (file.exists(meta_st5_file)) meta_st5_file else meta_leg_file

subtype_rds_file <- file.path(input_clust_dir, "endotype_assignments.rds")

if (!file.exists(combat_file) || !file.exists(meta_file)) {
  stop("Required inputs from Scripts 01 or 02 not found. Please verify upstream pipeline runs.")
}

# Ingest molecular subtype metadata generated in Script 06
if (file.exists(subtype_rds_file)) {
  subtype_df <- readRDS(subtype_rds_file)
} else {
  stop("Subtype assignments RDS not found in results/06_figure3_molecular_endotypes/.")
}

# Harmonize Subtype factor levels
if ("Subtype" %in% colnames(subtype_df)) {
  subtype_df$Subtype <- factor(as.character(subtype_df$Subtype), 
                               levels = c("Subtype_1", "Subtype_2", "Subtype_3", "Subtype1", "Subtype2", "Subtype3"),
                               labels = c("Subtype_1", "Subtype_2", "Subtype_3", "Subtype_1", "Subtype_2", "Subtype_3"))
}

# Ingest consensus Meta-DEPs (823 strict proteins)
meta_deps <- read.csv(meta_file, stringsAsFactors = FALSE)
col_stat  <- if ("Consensus_Regulation" %in% colnames(meta_deps)) "Consensus_Regulation" else "Final_Status"
col_prot  <- if ("Protein_Group" %in% colnames(meta_deps)) "Protein_Group" else "protein"

gold_823_prots <- meta_deps %>%
  filter(grepl("Strictly Validated", .data[[col_stat]])) %>%
  pull(.data[[col_prot]]) %>% unique()

# Ingest ComBat-corrected quantification matrix
prot_raw_full <- fread(combat_file, data.table = FALSE)
rownames(prot_raw_full) <- make.unique(as.character(prot_raw_full$Protein.Group))
prot_expr_mat <- as.matrix(prot_raw_full[, 7:ncol(prot_raw_full)])

common_samples <- intersect(colnames(prot_expr_mat), subtype_df$sample)
prot_mat   <- prot_expr_mat[, common_samples]
sub_clin   <- subtype_df[match(common_samples, subtype_df$sample), ]

# Global protein-to-gene mapping with semicolon splitting
all_uniprots <- rownames(prot_mat)
clean_uniprots <- sapply(strsplit(as.character(all_uniprots), ";"), `[`, 1) %>% 
  gsub("-.*|\\..*", "", .)

mapped_symbols <- suppressMessages(suppressWarnings(
  mapIds(org.Hs.eg.db, keys = clean_uniprots, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
))
mapped_symbols[is.na(mapped_symbols)] <- clean_uniprots[is.na(mapped_symbols)]
lookup_global <- setNames(as.character(mapped_symbols), all_uniprots)

# Structural and hematological noise filtering
blacklist_pattern <- "^(KRT|RPL|RPS|HBA|HBB|HBD|IGH|IGK|IGL|ALB)"
is_clean_protein  <- !grepl(blacklist_pattern, lookup_global[all_uniprots])

prot_mat <- prot_mat[is_clean_protein, ]
gold_823_clean_prots <- intersect(gold_823_prots[!grepl(blacklist_pattern, lookup_global[gold_823_prots])], rownames(prot_mat))

cat(sprintf("Loaded %d clean proteins across %d sporadic PD patients (823 Meta-DEPs retained: %d).\n", 
            nrow(prot_mat), ncol(prot_mat), length(gold_823_clean_prots)))

# ==============================================================================
# Step 2: 8-Covariate Pairwise Limma & Min-T Specificity (ST16, ST17, ST18, ST19)
# ==============================================================================
cat("\nExecuting 8-covariate pairwise Limma contrasts and Min-T bottleneck scoring...\n")

valid_covs <- c()
for (cv in BASE_COVARIATES) {
  if (cv %in% colnames(sub_clin) && length(unique(as.character(sub_clin[[cv]]))) > 1) {
    valid_covs <- c(valid_covs, cv)
  }
}

design_formula <- as.formula(paste("~ 0 + Subtype +", paste(valid_covs, collapse = " + ")))
design <- model.matrix(design_formula, data = sub_clin)
colnames(design)[1:3] <- c("Subtype_1", "Subtype_2", "Subtype_3")

fit <- lmFit(prot_mat, design)

contrast_matrix_pairwise <- makeContrasts(
  S1_vs_S2 = Subtype_1 - Subtype_2,
  S1_vs_S3 = Subtype_1 - Subtype_3,
  S2_vs_S3 = Subtype_2 - Subtype_3,
  levels = design
)
fit_pw <- eBayes(contrasts.fit(fit, contrast_matrix_pairwise))

res_12 <- topTable(fit_pw, coef = "S1_vs_S2", number = Inf, sort.by = "none") %>% mutate(UNIPROT = rownames(.))
res_13 <- topTable(fit_pw, coef = "S1_vs_S3", number = Inf, sort.by = "none") %>% mutate(UNIPROT = rownames(.))
res_23 <- topTable(fit_pw, coef = "S2_vs_S3", number = Inf, sort.by = "none") %>% mutate(UNIPROT = rownames(.))

# Calculate subtype-exclusive markers based on minimum t-statistic (Min-T)
s1_combined_all <- res_12 %>%
  inner_join(res_13, by = "UNIPROT", suffix = c("_12", "_13")) %>%
  mutate(Symbol = lookup_global[UNIPROT], Min_T_Score = pmin(t_12, t_13)) %>%
  filter(!grepl(blacklist_pattern, Symbol)) %>%
  filter(logFC_12 > 0 & logFC_13 > 0 & (P.Value_12 < 0.05 | P.Value_13 < 0.05) & Min_T_Score > 0.5) %>%
  arrange(desc(Min_T_Score))

s2_combined_all <- res_12 %>%
  inner_join(res_23, by = "UNIPROT", suffix = c("_12", "_23")) %>%
  mutate(Symbol = lookup_global[UNIPROT], Min_T_Score = pmin(-t_12, t_23)) %>%
  filter(!grepl(blacklist_pattern, Symbol)) %>%
  filter(logFC_12 < 0 & logFC_23 > 0 & (P.Value_12 < 0.05 | P.Value_23 < 0.05) & Min_T_Score > 0.5) %>%
  arrange(desc(Min_T_Score))

s3_combined_all <- res_13 %>%
  inner_join(res_23, by = "UNIPROT", suffix = c("_13", "_23")) %>%
  mutate(Symbol = lookup_global[UNIPROT], Min_T_Score = pmin(-t_13, -t_23)) %>%
  filter(!grepl(blacklist_pattern, Symbol)) %>%
  filter(logFC_13 < 0 & logFC_23 < 0 & (P.Value_13 < 0.05 | P.Value_23 < 0.05) & Min_T_Score > 0.5) %>%
  arrange(desc(Min_T_Score))

# Export Supplementary Tables ST16, ST17, ST18
write.csv(s1_combined_all, file.path(output_dir, "Supplementary_Table_16_Subtype1_Markers.csv"), row.names = FALSE)
write.csv(s2_combined_all, file.path(output_dir, "Supplementary_Table_17_Subtype2_Markers.csv"), row.names = FALSE)
write.csv(s3_combined_all, file.path(output_dir, "Supplementary_Table_18_Subtype3_Markers.csv"), row.names = FALSE)
cat("Exported Supplementary Tables ST16, ST17, and ST18 successfully.\n")

# Intersect with consensus Meta-DEPs
s1_intersect_df <- s1_combined_all %>% filter(UNIPROT %in% gold_823_clean_prots)
s2_intersect_df <- s2_combined_all %>% filter(UNIPROT %in% gold_823_clean_prots)
s3_intersect_df <- s3_combined_all %>% filter(UNIPROT %in% gold_823_clean_prots)

# Select top 15 non-redundant biomarkers per subtype for Supplementary Table 19 (45-biomarker panel)
top15_s1 <- head(s1_intersect_df %>% arrange(desc(Min_T_Score)) %>% pull(UNIPROT), 15)
top15_s2 <- head(s2_intersect_df %>% arrange(desc(Min_T_Score)) %>% pull(UNIPROT), 15)
top15_s3 <- head(s3_intersect_df %>% arrange(desc(Min_T_Score)) %>% pull(UNIPROT), 15)

table_st19_export <- bind_rows(
  s1_intersect_df %>% filter(UNIPROT %in% top15_s1) %>% mutate(Subtype_Assignment = "Subtype 1 Drivers", Rank_In_Subtype = 1:n()),
  s2_intersect_df %>% filter(UNIPROT %in% top15_s2) %>% mutate(Subtype_Assignment = "Subtype 2 Drivers", Rank_In_Subtype = 1:n()),
  s3_intersect_df %>% filter(UNIPROT %in% top15_s3) %>% mutate(Subtype_Assignment = "Subtype 3 Drivers", Rank_In_Subtype = 1:n())
) %>%
  dplyr::select(
    UNIPROT_ID = UNIPROT,
    Gene_Symbol = Symbol,
    Subtype_Assignment,
    Rank_In_Subtype,
    Min_T_Score,
    starts_with("logFC"),
    starts_with("P.Value")
  )

write.csv(table_st19_export, file.path(output_dir, "Supplementary_Table_19_Subtype_45_Fingerprint_Markers.csv"), row.names = FALSE)
cat("Exported Supplementary Table 19 (45 Fingerprint Biomarkers).\n")

# ==============================================================================
# Step 3: Figure 4A (UpSet Diagram) & Figure 4B (Top 15 Fingerprints Heatmap)
# ==============================================================================
cat("\nGenerating Figure 4A (UpSet diagram) and Figure 4B (fingerprints heatmap)...\n")

universal_list_all <- list(
  "Universal PD (823)" = gold_823_clean_prots, 
  "S1 Drivers"         = s1_combined_all$UNIPROT, 
  "S2 Drivers"         = s2_combined_all$UNIPROT, 
  "S3 Drivers"         = s3_combined_all$UNIPROT
)

pdf(file.path(output_dir, "Figure4A_UpSet_Subtypes_vs_GlobalPD.pdf"), width = 8.5, height = 6.0)
upset(fromList(universal_list_all), 
      sets = c("Universal PD (823)", "S3 Drivers", "S2 Drivers", "S1 Drivers"),
      keep.order = TRUE, order.by = "freq", nsets = 4,
      mainbar.y.label = "Number of Intersecting Proteins", 
      sets.x.label = "Total Proteins per Category",
      text.scale = c(1.3, 1.3, 1.1, 1.1, 1.2, 1.2), point.size = 3.5, line.size = 1, mb.ratio = c(0.7, 0.3),
      queries = list(
        list(query = intersects, params = list("Universal PD (823)", "S1 Drivers"), color = "#E64B35", active = TRUE),
        list(query = intersects, params = list("Universal PD (823)", "S2 Drivers"), color = "#4DBBD5", active = TRUE),
        list(query = intersects, params = list("Universal PD (823)", "S3 Drivers"), color = "#00A087", active = TRUE)
      ))
dev.off()

# Extract top 5 per subtype for Figure 4B main text heatmap (15 markers)
n_top_main <- 5
top5_s1 <- head(top15_s1, n_top_main)
top5_s2 <- head(top15_s2, n_top_main)
top5_s3 <- head(top15_s3, n_top_main)
top_prots_main <- c(top5_s1, top5_s2, top5_s3)

raw_mat_main <- prot_mat[top_prots_main, ]
mat_scaled_main <- t(scale(t(raw_mat_main)))
mat_capped_main <- pmax(pmin(mat_scaled_main, 2.5), -3.0)
rownames(mat_capped_main) <- lookup_global[top_prots_main]

smooth_reflective <- function(x, w = 5) {
  n <- length(x)
  if (n <= w) return(x)
  kernel <- c(0.08, 0.24, 0.36, 0.24, 0.08)
  padded <- c(x[2], x[1], x, x[n], x[n-1])
  smoothed <- stats::filter(padded, filter = kernel, sides = 2)
  return(as.numeric(smoothed[3:(n + 2)]))
}

ordered_cols_main <- c()
mat_smoothed_list_main <- list()
sub_gene_dict_main <- list(Subtype_1 = top5_s1, Subtype_2 = top5_s2, Subtype_3 = top5_s3)

for (st in c("Subtype_1", "Subtype_2", "Subtype_3")) {
  st_samps <- sub_clin$sample[sub_clin$Subtype == st]
  st_targets <- lookup_global[sub_gene_dict_main[[st]]]
  st_mat <- mat_capped_main[st_targets, st_samps, drop = FALSE]
  
  pca <- prcomp(t(st_mat), center = TRUE, scale. = FALSE)
  if (cor(pca$x[, 1], colMeans(st_mat)) < 0) pca$x[, 1] <- -pca$x[, 1]
  st_ordered <- st_samps[order(pca$x[, 1], decreasing = TRUE)]
  ordered_cols_main <- c(ordered_cols_main, st_ordered)
  
  st_full_mat <- mat_capped_main[, st_ordered]
  if (ncol(st_full_mat) > 6) {
    st_full_smooth <- t(apply(st_full_mat, 1, smooth_reflective))
    colnames(st_full_smooth) <- st_ordered
    mat_smoothed_list_main[[st]] <- st_full_smooth
  } else {
    mat_smoothed_list_main[[st]] <- st_full_mat
  }
}

heatmap_mat_final_main <- do.call(cbind, mat_smoothed_list_main)

col_fun_heatmap <- colorRamp2(
  c(-3.0, -1.8, -0.8, 0, 0.5, 1.3, 2.5),
  c("#3A4866", "#5D6E94", "#92A1C2", "#FAF8F8", "#ECC3CC", "#D77C91", "#A83853")
)

column_split <- factor(sub_clin$Subtype[match(ordered_cols_main, sub_clin$sample)], 
                       levels = c("Subtype_1", "Subtype_2", "Subtype_3"))
row_split <- factor(rep(c("Subtype 1 Drivers", "Subtype 2 Drivers", "Subtype 3 Drivers"), each = n_top_main),
                    levels = c("Subtype 1 Drivers", "Subtype 2 Drivers", "Subtype 3 Drivers"))

top_ann <- HeatmapAnnotation(
  Subtype = column_split,
  col = list(Subtype = c("Subtype_1" = "#d73027", "Subtype_2" = "#3182bd", "Subtype_3" = "#fdae61")),
  show_annotation_name = TRUE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  simple_anno_size = unit(4, "mm"),
  border = FALSE
)

pdf(file.path(output_dir, "Figure4B_Subtype_Exclusive_Fingerprints_Top15_MainText.pdf"), 
    width = 10.5, height = 4.8)

ht_main <- Heatmap(
  heatmap_mat_final_main,
  name = "Z-score",
  col = col_fun_heatmap,
  column_split = column_split,
  row_split = row_split,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  show_row_dend = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 11, fontface = "italic"), 
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  top_annotation = top_ann,
  column_gap = unit(1.2, "mm"),
  row_gap = unit(1.5, "mm"),
  border = TRUE,
  heatmap_legend_param = list(
    title = "Z-score",
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8),
    at = c(-2, 0, 2),
    labels = c("-2 (Down)", " 0", "+2 (Up)"),
    legend_height = unit(3.5, "cm"),
    grid_width = unit(4, "mm")
  )
)
draw(ht_main, padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()

# ==============================================================================
# Step 4: Subtype Pathway Enrichment with GO.db Topological Pruning (Figure 4C)
# ==============================================================================
cat("\nRunning background-constrained GO-BP enrichment with topological filtering...\n")

bg_raw_ids <- unique(gsub("-.*|\\..*", "", sapply(strsplit(as.character(prot_raw_full[[1]]), ";"), `[`, 1)))
universe_entrez <- suppressMessages(
  bitr(bg_raw_ids, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
)$ENTREZID %>% unique()

n_go_limit <- 120
s1_entrez <- na.omit(suppressMessages(bitr(gsub("-.*|\\..*", "", head(s1_intersect_df$UNIPROT, n_go_limit)), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db))$ENTREZID) %>% unique()
s2_entrez <- na.omit(suppressMessages(bitr(gsub("-.*|\\..*", "", head(s2_intersect_df$UNIPROT, n_go_limit)), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db))$ENTREZID) %>% unique()
s3_entrez <- na.omit(suppressMessages(bitr(gsub("-.*|\\..*", "", head(s3_intersect_df$UNIPROT, n_go_limit)), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db))$ENTREZID) %>% unique()

gene_clusters_subtypes <- list("S1" = s1_entrez, "S2" = s2_entrez, "S3" = s3_entrez)

comp_bp_subtypes <- compareCluster(
  geneCluster   = gene_clusters_subtypes, 
  fun           = "enrichGO", 
  OrgDb         = org.Hs.eg.db, 
  ont           = "BP", 
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  minGSSize     = 10,
  maxGSSize     = 500,
  universe      = universe_entrez,
  readable      = TRUE
)

df_go_raw <- as.data.frame(comp_bp_subtypes)

# Filter developmental ancestor terms (GO:0048856 / GO:0040007)
bp_ancestors <- as.list(GOBPANCESTOR)
is_developmental_node <- function(go_id) {
  anc <- bp_ancestors[[go_id]]
  if (is.null(anc)) return(FALSE)
  return("GO:0048856" %in% anc || "GO:0040007" %in% anc)
}

df_go_filtered <- df_go_raw %>% filter(!sapply(ID, is_developmental_node))

# Prune overlapping terms per cluster (Overlap >= 0.60 | Jaccard >= 0.35)
prune_pure_mathematical <- function(df_cluster, max_terms = 6, overlap_cut = 0.60, jaccard_cut = 0.35) {
  df_sig <- df_cluster %>% filter(pvalue < 0.05) %>% arrange(pvalue)
  if (nrow(df_sig) <= 1) return(df_sig)
  selected_rows <- c(1)
  for (i in 2:nrow(df_sig)) {
    if (length(selected_rows) >= max_terms) break
    current_genes <- unlist(strsplit(df_sig$geneID[i], "/"))
    is_redundant  <- FALSE
    for (sel_idx in selected_rows) {
      existing_genes <- unlist(strsplit(df_sig$geneID[sel_idx], "/"))
      inter_len <- length(intersect(current_genes, existing_genes))
      union_len <- length(union(current_genes, existing_genes))
      min_len   <- min(length(current_genes), length(existing_genes))
      if ((inter_len / min_len >= overlap_cut) || (inter_len / union_len >= jaccard_cut)) {
        is_redundant <- TRUE
        break
      }
    }
    # Bug fixed: using is_redundant (consistent variable naming)
    if (!is_redundant) selected_rows <- c(selected_rows, i)
  }
  return(df_sig[selected_rows, ])
}

df_pruned_list <- list()
for (cl in c("S1", "S2", "S3")) {
  sub_df <- df_go_filtered %>% filter(Cluster == cl)
  if (nrow(sub_df) > 0) {
    df_pruned_list[[cl]] <- prune_pure_mathematical(sub_df, max_terms = 6, overlap_cut = 0.60, jaccard_cut = 0.35)
  }
}

top_terms_df <- do.call(rbind, df_pruned_list) %>%
  separate(GeneRatio, into = c("n", "N"), sep = "/", remove = FALSE) %>%
  mutate(
    GeneRatio_num = as.numeric(n) / as.numeric(N),
    Log_P = -log10(pvalue),
    Description_Clean = str_wrap(str_to_title(Description), width = 36)
  )

write.csv(top_terms_df, file.path(output_dir, "Supplementary_Table_Subtype_Mechanisms_Pruned.csv"), row.names = FALSE)

all_selected_terms <- unique(top_terms_df$Description_Clean)
plot_data_final <- expand.grid(
  Description_Clean = all_selected_terms,
  Cluster = c("S1", "S2", "S3"),
  stringsAsFactors = FALSE
) %>%
  left_join(top_terms_df %>% dplyr::select(Description_Clean, Cluster, Log_P, GeneRatio_num), 
            by = c("Description_Clean", "Cluster")) %>%
  mutate(
    Log_P = replace_na(Log_P, 0),
    GeneRatio_num = replace_na(GeneRatio_num, 0),
    Cluster = factor(Cluster, levels = c("S1", "S2", "S3"))
  )

order_y <- top_terms_df %>% arrange(Cluster, desc(Log_P)) %>% pull(Description_Clean) %>% unique() %>% rev()
plot_data_final$Description_Clean <- factor(plot_data_final$Description_Clean, levels = order_y)

p_fig4c <- ggplot(plot_data_final, aes(x = Cluster, y = Description_Clean)) +
  geom_point(data = filter(plot_data_final, Log_P > 0), 
             aes(size = GeneRatio_num, fill = Log_P), shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradientn(colors = c("#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026"),
                       name = expression(bold(-log[10](P-value))), limits = c(1.3, 5.0), oob = scales::squish) +
  scale_size_continuous(range = c(3.5, 8.5), name = "Gene Ratio") +
  labs(title = "Core Mechanistic Divergence across Molecular Endotypes", x = NULL, y = NULL) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13.5),
    axis.text.x = element_text(face = "bold", size = 12.5),
    axis.text.y = element_text(size = 10, lineheight = 0.85),
    panel.grid.major.y = element_line(color = "grey88", linetype = "dashed"),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(output_dir, "Figure4C_Subtype_Mechanisms_Dotplot.pdf"), p_fig4c, width = 9.5, height = 8.0, device = pdf_device)

# ==============================================================================
# Step 5: Cellular Deconvolution & Figure 4D/E (Supplementary Table 20)
# ==============================================================================
cat("\nPerforming PanglaoDB cellular microenvironment deconvolution...\n")

# Increase download timeout for large reference files
options(timeout = max(600, getOption("timeout")))

# Smart candidates search: project data/ directory, root directory, or common download folders
panglao_candidates <- c(
  file.path(data_dir, "PanglaoDB_markers_27_Mar_2020.tsv.gz"),
  file.path(data_dir, "PanglaoDB_markers_27_Mar_2020.tsv"),
  "PanglaoDB_markers_27_Mar_2020.tsv.gz"
)

panglao_path <- NA_character_
for (cand in panglao_candidates) {
  if (file.exists(cand)) {
    panglao_path <- cand
    break
  }
}

# Fuzzy search in data directory if not matched directly
if (is.na(panglao_path)) {
  fuzzy_files <- list.files(data_dir, pattern = "PanglaoDB.*\\.tsv(\\.gz)?$", full.names = TRUE, recursive = TRUE)
  if (length(fuzzy_files) > 0) panglao_path <- fuzzy_files[1]
}

# Automated download fallback for external reviewers if completely absent
if (is.na(panglao_path) || !file.exists(panglao_path)) {
  dest_file <- file.path(data_dir, "PanglaoDB_markers_27_Mar_2020.tsv.gz")
  cat("PanglaoDB markers file not found locally. Attempting automatic download...\n")
  panglao_url <- "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"
  
  dl_status <- tryCatch({
    download.file(panglao_url, destfile = dest_file, mode = "wb", method = "libcurl")
  }, error = function(e) 1)
  
  if (dl_status == 0 && file.exists(dest_file) && file.info(dest_file)$size > 100000) {
    panglao_path <- dest_file
    cat("PanglaoDB markers file downloaded successfully.\n")
  } else {
    stop("\n[Data Notice]: PanglaoDB reference file not found.\n",
         "Please place 'PanglaoDB_markers_27_Mar_2020.tsv.gz' into the 'data/' directory to proceed.")
  }
}

cat(sprintf("Using PanglaoDB database file: %s\n", panglao_path))
panglao_db <- fread(panglao_path, sep = "\t", data.table = FALSE) %>% filter(str_detect(species, "Hs"))

# Convert proteomic matrix to Gene Symbols
gene_map_full <- suppressMessages(suppressWarnings(
  bitr(clean_uniprots, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
)) %>% distinct(UNIPROT, .keep_all = TRUE)

prot_symbol_mat <- prot_raw_full %>% 
  mutate(UNIPROT = sapply(strsplit(as.character(prot_raw_full[[1]]), ";"), `[`, 1) %>% gsub("-.*|\\..*", "", .)) %>% 
  inner_join(gene_map_full, by = "UNIPROT") %>%
  dplyr::select(-1, -UNIPROT) %>% 
  group_by(SYMBOL) %>% 
  summarise(across(where(is.numeric), ~mean(., na.rm = TRUE))) %>% 
  column_to_rownames("SYMBOL") %>% 
  as.matrix()

prot_symbol_mat <- prot_symbol_mat[, common_samples]

cell_signatures <- list()
for (cell in unique(panglao_db$`cell type`)) {
  genes <- panglao_db %>% filter(`cell type` == cell) %>% pull(`official gene symbol`) %>% unique()
  valid_g <- intersect(genes, rownames(prot_symbol_mat))
  if (length(valid_g) >= 3) cell_signatures[[cell]] <- valid_g
}

param_cells <- ssgseaParam(exprData = prot_symbol_mat, geneSets = cell_signatures, minSize = 3)
gsva_res_cells <- gsva(param_cells)

score_df_cells <- as.data.frame(t(gsva_res_cells), check.names = FALSE) %>% 
  rownames_to_column("sample") %>% 
  inner_join(sub_clin %>% dplyr::select(sample, Subtype), by = "sample")

# Export Supplementary Table 20 (All 126 lineages)
table_st20_stats <- map_df(names(cell_signatures), function(ct) {
  s_vec <- score_df_cells[[ct]]
  sub_vec <- score_df_cells$Subtype
  kw_res <- kruskal.test(s_vec ~ sub_vec)
  
  data.frame(
    Cell_Type = ct,
    Gene_Count = length(cell_signatures[[ct]]),
    Subtype_1_Mean_Score = round(mean(s_vec[sub_vec == "Subtype_1"], na.rm = TRUE), 4),
    Subtype_2_Mean_Score = round(mean(s_vec[sub_vec == "Subtype_2"], na.rm = TRUE), 4),
    Subtype_3_Mean_Score = round(mean(s_vec[sub_vec == "Subtype_3"], na.rm = TRUE), 4),
    Kruskal_Wallis_Chi2 = round(as.numeric(kw_res$statistic), 3),
    Kruskal_Wallis_P = kw_res$p.value
  )
}) %>%
  mutate(Global_FDR = p.adjust(Kruskal_Wallis_P, method = "BH")) %>%
  arrange(Kruskal_Wallis_P)

write.csv(table_st20_stats, file.path(output_dir, "Supplementary_Table_20_Cellular_Deconvolution_126_Lineages.csv"), row.names = FALSE)
cat("Exported Supplementary Table 20 (PanglaoDB 126 Lineages Full Stats).\n")

# Main Figure 4D/E & Supp Fig S5A cell definitions
main_cns_cells    <- c("Neurons", "Astrocytes", "Microglia", "Oligodendrocytes")
main_periph_cells <- c("Endothelial cells", "Monocytes", "Macrophages", "T cells")
main_all_cells    <- c(main_cns_cells, main_periph_cells)

supp_cns_cells    <- c("Purkinje neurons", "Ependymal cells", "Neural stem/precursor cells", "Schwann cells")
supp_periph_cells <- c("Pericytes", "Neutrophils", "Platelets", "B cells")
supp_all_cells    <- c(supp_cns_cells, supp_periph_cells)

subtype_colors <- c("Subtype_1" = "#d73027", "Subtype_2" = "#4575b4", "Subtype_3" = "#fdae61")
my_comparisons <- list(c("Subtype_1", "Subtype_2"), c("Subtype_2", "Subtype_3"), c("Subtype_1", "Subtype_3"))

draw_cell_violins <- function(cell_vec, title_text) {
  valid_cells <- intersect(cell_vec, colnames(score_df_cells))
  df_plot <- score_df_cells %>%
    dplyr::select(sample, Subtype, all_of(valid_cells)) %>%
    pivot_longer(cols = all_of(valid_cells), names_to = "Cell_Type", values_to = "Score")
  df_plot$Cell_Type <- factor(df_plot$Cell_Type, levels = cell_vec)
  
  ggplot(df_plot, aes(x = Subtype, y = Score, fill = Subtype)) +
    geom_violin(trim = FALSE, alpha = 0.45, color = NA, scale = "width", width = 0.85) +
    geom_boxplot(width = 0.16, fill = "white", outlier.shape = NA, color = "grey25", linewidth = 0.55) +
    facet_wrap(~ Cell_Type, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = subtype_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", 
                       step.increase = 0.08, tip.length = 0.01, size = 3.6) +
    labs(title = title_text, x = NULL, y = "Deconvoluted Score (ssGSEA)") +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 11.5, hjust = 0.5),
      strip.background = element_rect(fill = "grey93", color = "black", linewidth = 0.6),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(face = "bold", size = 9.5, angle = 20, hjust = 1)
    )
}

p_fig4d <- draw_cell_violins(main_cns_cells, "CNS Cellular Microenvironment (Figure 4D)")
p_fig4e <- draw_cell_violins(main_periph_cells, "Peripheral Immune & Vascular Lineages (Figure 4E)")

# ==============================================================================
# Step 6: Multivariable Partial Correlation Heatmap (Figure 4F)
# ==============================================================================
cat("\nComputing 8-covariate partial correlation for Figure 4F...\n")

clinical_traits <- c("UPDRS3", "NMSS_Total", "MMSE", "HAMA", "HAMD", "PDSS", "UPSIT")
actual_traits   <- intersect(clinical_traits, colnames(sub_clin))

trait_labels_final <- c(
  "UPDRS3"     = "Motor (UPDRS-III)", 
  "NMSS_Total" = "Non-motor (NMSS)", 
  "MMSE"       = "Cognition (MMSE)", 
  "HAMA"       = "Anxiety (HAMA)", 
  "HAMD"       = "Depression (HAMD)", 
  "PDSS"       = "Sleep (PDSS)", 
  "UPSIT"      = "Olfaction (UPSIT)"
)

calculate_8cov_pcor <- function(feature_mat, feat_names, clin_df, traits_vector) {
  r_mat <- matrix(NA_real_, nrow = length(feat_names), ncol = length(traits_vector),
                  dimnames = list(feat_names, traits_vector))
  p_mat <- matrix(NA_real_, nrow = length(feat_names), ncol = length(traits_vector),
                  dimnames = list(feat_names, traits_vector))
  
  df_feat <- as.data.frame(t(feature_mat[feat_names, , drop = FALSE]), check.names = FALSE) %>% rownames_to_column("sample")
  merged_df <- df_feat %>% inner_join(clin_df, by = "sample")
  
  for (feat in feat_names) {
    for (tr in traits_vector) {
      clean_data <- data.frame(
        score     = as.numeric(merged_df[[feat]]), 
        trait_val = as.numeric(merged_df[[tr]]),
        Age       = as.numeric(merged_df$Age), 
        Sex       = as.numeric(merged_df$Sex), 
        BMI       = as.numeric(merged_df$BMI),
        hepatic   = as.numeric(merged_df$hepatic_disease), 
        kidney    = as.numeric(merged_df$kidney.function),
        CV_Met    = as.numeric(merged_df$CV_Metabolic_Cat), 
        LEDD      = as.numeric(merged_df$LEDD), 
        Duration  = as.numeric(merged_df$Duration)
      ) %>% drop_na()
      
      if (nrow(clean_data) >= 30) {
        covars <- clean_data[, c("Age", "Sex", "BMI", "hepatic", "kidney", "CV_Met", "LEDD", "Duration")]
        pcor_res <- tryCatch({
          pcor.test(clean_data$score, clean_data$trait_val, covars, method = "spearman")
        }, error = function(e) NULL)
        
        if (!is.null(pcor_res)) {
          r_mat[feat, tr] <- pcor_res$estimate
          p_mat[feat, tr] <- pcor_res$p.value
        }
      }
    }
  }
  colnames(r_mat) <- trait_labels_final[colnames(r_mat)]
  colnames(p_mat) <- trait_labels_final[colnames(p_mat)]
  return(list(r = r_mat, p = p_mat))
}

pcor_cells_res <- calculate_8cov_pcor(gsva_res_cells, main_all_cells, sub_clin, actual_traits)
r_cells <- pcor_cells_res$r
p_cells <- pcor_cells_res$p

col_fun_pcor <- colorRamp2(c(-0.35, 0, 0.35), c("#313695", "#FAF8F8", "#a50026"))
cell_split_tags <- factor(ifelse(main_all_cells %in% main_cns_cells, "CNS Lineages", "Peripheral Lineages"),
                          levels = c("CNS Lineages", "Peripheral Lineages"))

pdf(file.path(output_dir, "Figure4F_Cellular_Clinical_Pcor_Heatmap.pdf"), width = 9.2, height = 6.2)
ht_cells <- Heatmap(
  r_cells,
  name = "Partial r",
  col = col_fun_pcor,
  row_split = cell_split_tags,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "plain"),
  column_names_gp = gpar(fontsize = 10.5, fontface = "bold"),
  column_names_rot = 45,
  row_title_gp = gpar(fontsize = 10.5, fontface = "bold"),
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(p_cells[i, j]) && p_cells[i, j] < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 13, fontface = "bold", col = "black"))
    }
  },
  column_title = "Figure 4F: Cellular Microenvironment vs Clinical Traits (8-Covariate Adjusted)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)
draw(ht_cells, padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()

# ==============================================================================
# Step 7: Assemble Main Figure 4 Layout
# ==============================================================================
cat("\nAssembling Main Figure 4...\n")

fig4_row_deconv <- (p_fig4d | p_fig4e) + plot_layout(widths = c(1.0, 1.0))

ggsave(file.path(output_dir, "Figure4DE_Deconvolution_Row.pdf"), fig4_row_deconv, width = 14.0, height = 4.2, device = pdf_device)
ggsave(file.path(output_dir, "Figure4DE_Deconvolution_Row.png"), fig4_row_deconv, width = 14.0, height = 4.2, dpi = 300)

# ==============================================================================
# Step 8: Generate Supplementary Figure 5 (Panels S5A & S5B)
# ==============================================================================
cat("Generating Supplementary Figure 5 (S5A Extended Lineages & S5B Pathway Pcor)...\n")

# Supp Fig 5A: Extended cell lineages violins
p_supp5a_top <- draw_cell_violins(supp_cns_cells, "Extended CNS & Neuroglial Support Lineages")
p_supp5a_bot <- draw_cell_violins(supp_periph_cells, "Extended Peripheral Vascular & Immune Lineages")
supp_fig5a_grid <- p_supp5a_top / p_supp5a_bot

# Supp Fig 5B: 17 Core Subtype Pathways partial correlation heatmap
mech_gene_sets  <- list()
valid_term_rows <- list()

for (i in 1:nrow(top_terms_df)) {
  term_name <- str_to_title(top_terms_df$Description[i])
  raw_genes <- unlist(strsplit(top_terms_df$geneID[i], "/"))
  valid_genes <- intersect(raw_genes, rownames(prot_symbol_mat))
  if (length(valid_genes) >= 3) {
    mech_gene_sets[[term_name]]  <- valid_genes
    valid_term_rows[[term_name]] <- top_terms_df$Cluster[i]
  }
}

param_mech <- ssgseaParam(exprData = prot_symbol_mat, geneSets = mech_gene_sets, minSize = 3)
gsva_res_mech <- gsva(param_mech)
actual_mechs  <- rownames(gsva_res_mech)

pcor_mech_res <- calculate_8cov_pcor(gsva_res_mech, actual_mechs, sub_clin, actual_traits)
r_mech <- pcor_mech_res$r
p_mech <- pcor_mech_res$p

pathway_subtype_tags <- factor(paste0("Subtype ", gsub("S", "", unlist(valid_term_rows[actual_mechs])), " Signatures"),
                               levels = c("Subtype 1 Signatures", "Subtype 2 Signatures", "Subtype 3 Signatures"))

pdf(file.path(output_dir, "Supplementary_Figure_5B_Pathway_Pcor_Heatmap.pdf"), width = 10.2, height = 7.8)
ht_supp5b <- Heatmap(
  r_mech,
  name = "Partial r",
  col = col_fun_pcor,
  row_split = pathway_subtype_tags,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 9.5, fontface = "plain"),
  column_names_gp = gpar(fontsize = 10.5, fontface = "bold"),
  column_names_rot = 45,
  row_title_gp = gpar(fontsize = 10.5, fontface = "bold"),
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(p_mech[i, j]) && p_mech[i, j] < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 13, fontface = "bold", col = "black"))
    }
  },
  column_title = "Supplementary Figure 5B: Subtype Pathways vs Clinical Traits (8-Covariate Adjusted)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)
draw(ht_supp5b, padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()

# Save combined Supplementary Figure 5A
ggsave(file.path(output_dir, "Supplementary_Figure_5A_Extended_Cell_Lineages.pdf"), 
       supp_fig5a_grid, width = 12.0, height = 7.2, device = pdf_device)
ggsave(file.path(output_dir, "Supplementary_Figure_5A_Extended_Cell_Lineages.png"), 
       supp_fig5a_grid, width = 12.0, height = 7.2, dpi = 300)

# Save intermediate objects for downstream classifier training (Script 08)
saveRDS(top15_s1, file.path(output_dir, "top15_s1_markers.rds"))
saveRDS(top15_s2, file.path(output_dir, "top15_s2_markers.rds"))
saveRDS(top15_s3, file.path(output_dir, "top15_s3_markers.rds"))
saveRDS(table_st19_export, file.path(output_dir, "locked_45_fingerprint_markers.rds"))
saveRDS(gsva_res_cells, file.path(output_dir, "gsva_res_cells.rds"))

cat("\nPhase 7 pipeline complete. Figure 4, Supp Fig 5, and ST16–ST20 successfully exported to:", output_dir, "\n")