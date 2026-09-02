# ==============================================================================
# Script: 08_subtype_mechanisms_and_cellular_deconvolution.R
# Purpose: Subtype biomarker specificity, pathway enrichment, cellular deconvolution,
#          and generation of Figure 4, Supplementary Figure 6, and Tables ST16 & ST17
# Project: Large-scale plasma proteomics identifies molecularly distinct PD endotypes
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

# Directory setup
data_dir        <- "data"
input_dea_dir   <- "results/02_figure1_meta_dea"
input_clust_dir <- "results/07_figure3_molecular_endotypes"
output_dir      <- "results/08_subtype_mechanisms_and_deconvolution"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

find_file_smart <- function(pattern, default_path) {
  files <- list.files(data_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  files <- files[!grepl("^~\\$", basename(files))]
  if (length(files) > 0) return(files[1])
  return(default_path)
}

BASE_COVARIATES <- c("Age", "Sex", "BMI", "hepatic_disease", "kidney.function", "CV_Metabolic_Cat", "LEDD", "Duration")

cat("=== Phase 8: Subtype Mechanism Profiling and Cellular Deconvolution ===\n")

# ==============================================================================
# Step 1: Ingest Preprocessed Matrices and Subtype Metadata
# ==============================================================================
cat("Loading preprocessed expression matrices and molecular subtype assignments...\n")

clin_file    <- file.path(data_dir, "metadata_clinical_n1119.csv")
combat_file  <- file.path("results/01_preprocessed_data", "3_ComBat_Corrected_Matrix.tsv")
meta_file    <- file.path(input_dea_dir, "Table_ST10_Final_823_Strict_Meta_DEPs.csv")
subtype_file <- file.path(input_clust_dir, "Supplementary_Table_7_Molecular_Endotypes_Baseline.csv")

if (!file.exists(clin_file) || !file.exists(combat_file) || !file.exists(meta_file)) {
  stop("Required input files not found. Please verify previous pipeline scripts have run.")
}

clinical_raw <- read.csv(clin_file, stringsAsFactors = FALSE)
colnames(clinical_raw) <- make.names(colnames(clinical_raw))

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
gold_823_prots <- meta_deps %>%
  filter(Final_Status %in% c("Strictly Validated Up", "Strictly Validated Down")) %>%
  pull(protein) %>% unique()

prot_raw_full <- fread(combat_file, data.table = FALSE)
rownames(prot_raw_full) <- make.unique(as.character(prot_raw_full$Protein.Group))
prot_expr_mat <- as.matrix(prot_raw_full[, 7:ncol(prot_raw_full)])

common_samples <- intersect(colnames(prot_expr_mat), clinical_pd$sample)
prot_mat   <- prot_expr_mat[, common_samples]
subtype_df <- clinical_pd[match(common_samples, clinical_pd$sample), ]

# Retrieve Subtype labels from Phase 7 output
if ("Subtype" %in% colnames(subtype_df)) {
  subtype_df$Subtype <- factor(subtype_df$Subtype, levels = c("Subtype1", "Subtype2", "Subtype3"))
} else {
  subtype_map_file <- file.path(input_clust_dir, "ConsensusCluster_Molecular")
  # Fallback to direct clinical table if present
}

# Gene mapping and blacklist filtering
clean_uniprots <- gsub("-.*|\\..*", "", rownames(prot_mat))
mapped_symbols <- suppressMessages(suppressWarnings(
  mapIds(org.Hs.eg.db, keys = clean_uniprots, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
))
mapped_symbols[is.na(mapped_symbols)] <- clean_uniprots[is.na(mapped_symbols)]
lookup_global <- setNames(as.character(mapped_symbols), rownames(prot_mat))

blacklist_pattern <- "^(KRT|RPL|RPS|HBA|HBB|HBD|IGH|IGK|IGL|ALB)"
is_clean_protein  <- !grepl(blacklist_pattern, lookup_global[rownames(prot_mat)])
prot_mat <- prot_mat[is_clean_protein, ]

gold_823_clean_prots <- intersect(gold_823_prots[!grepl(blacklist_pattern, lookup_global[gold_823_prots])], rownames(prot_mat))

# ==============================================================================
# Step 2: Subtype-Specific Pairwise Limma & Min-T Metric (Table ST16)
# ==============================================================================
cat("Executing 8-covariate pairwise Limma models and calculating Min-T specificity scores...\n")

sub_clin <- subtype_df %>% drop_na(all_of(BASE_COVARIATES))
sub_samps <- intersect(colnames(prot_mat), sub_clin$sample)
sub_clin <- sub_clin[match(sub_samps, sub_clin$sample), ]
sub_expr <- prot_mat[, sub_samps]

design_formula <- as.formula(paste("~ 0 + Subtype +", paste(BASE_COVARIATES, collapse = " + ")))
design <- model.matrix(design_formula, data = sub_clin)
colnames(design)[1:3] <- c("Subtype1", "Subtype2", "Subtype3")

fit <- lmFit(sub_expr, design)

contrast_matrix <- makeContrasts(
  S1_vs_S2 = Subtype1 - Subtype2,
  S1_vs_S3 = Subtype1 - Subtype3,
  S2_vs_S3 = Subtype2 - Subtype3,
  levels = design
)
fit_pw <- eBayes(contrasts.fit(fit, contrast_matrix))

res_12 <- topTable(fit_pw, coef = "S1_vs_S2", number = Inf, sort.by = "none") %>% mutate(UNIPROT = rownames(.))
res_13 <- topTable(fit_pw, coef = "S1_vs_S3", number = Inf, sort.by = "none") %>% mutate(UNIPROT = rownames(.))
res_23 <- topTable(fit_pw, coef = "S2_vs_S3", number = Inf, sort.by = "none") %>% mutate(UNIPROT = rownames(.))

# Bottleneck Min-T specificity metrics
s1_markers <- res_12 %>%
  inner_join(res_13, by = "UNIPROT", suffix = c("_12", "_13")) %>%
  mutate(Symbol = lookup_global[UNIPROT], Min_T_Score = pmin(t_12, t_13)) %>%
  filter(!grepl(blacklist_pattern, Symbol) & logFC_12 > 0 & logFC_13 > 0 & Min_T_Score > 0.5) %>%
  arrange(desc(Min_T_Score))

s2_markers <- res_12 %>%
  inner_join(res_23, by = "UNIPROT", suffix = c("_12", "_23")) %>%
  mutate(Symbol = lookup_global[UNIPROT], Min_T_Score = pmin(-t_12, t_23)) %>%
  filter(!grepl(blacklist_pattern, Symbol) & logFC_12 < 0 & logFC_23 > 0 & Min_T_Score > 0.5) %>%
  arrange(desc(Min_T_Score))

s3_markers <- res_13 %>%
  inner_join(res_23, by = "UNIPROT", suffix = c("_13", "_23")) %>%
  mutate(Symbol = lookup_global[UNIPROT], Min_T_Score = pmin(-t_13, -t_23)) %>%
  filter(!grepl(blacklist_pattern, Symbol) & logFC_13 < 0 & logFC_23 < 0 & Min_T_Score > 0.5) %>%
  arrange(desc(Min_T_Score))

table_st16_export <- bind_rows(
  s1_markers %>% mutate(Specific_Endotype = "Subtype 1"),
  s2_markers %>% mutate(Specific_Endotype = "Subtype 2"),
  s3_markers %>% mutate(Specific_Endotype = "Subtype 3")
)
write.csv(table_st16_export, file.path(output_dir, "Supplementary_Table_16_Subtype_Biomarker_Specificity.csv"), row.names = FALSE)

# ==============================================================================
# Step 3: Generate Figure 4A (UpSet Diagram) & Figure 4B (Top 15 Biomarkers Heatmap)
# ==============================================================================
cat("Generating Figure 4A (UpSet diagram) and Figure 4B (fingerprint heatmap)...\n")

upset_input_list <- list(
  "Universal PD (823)" = gold_823_clean_prots,
  "S1 Drivers"         = s1_markers$UNIPROT,
  "S2 Drivers"         = s2_markers$UNIPROT,
  "S3 Drivers"         = s3_markers$UNIPROT
)

pdf(file.path(output_dir, "Figure4A_UpSet_Subtypes_vs_GlobalPD.pdf"), width = 8.5, height = 5.8)
upset(fromList(upset_input_list),
      sets = c("Universal PD (823)", "S3 Drivers", "S2 Drivers", "S1 Drivers"),
      keep.order = TRUE, order.by = "freq", nsets = 4,
      mainbar.y.label = "Number of Intersecting Proteins",
      sets.x.label = "Total Proteins per Category",
      text.scale = c(1.2, 1.2, 1.0, 1.0, 1.1, 1.1), point.size = 3.2, line.size = 0.9,
      queries = list(
        list(query = intersects, params = list("Universal PD (823)", "S1 Drivers"), color = "#E64B35", active = TRUE),
        list(query = intersects, params = list("Universal PD (823)", "S2 Drivers"), color = "#4DBBD5", active = TRUE),
        list(query = intersects, params = list("Universal PD (823)", "S3 Drivers"), color = "#00A087", active = TRUE)
      ))
dev.off()

# Extract top 5 markers per subtype for main text Figure 4B
s1_intersect <- s1_markers %>% filter(UNIPROT %in% gold_823_clean_prots)
s2_intersect <- s2_markers %>% filter(UNIPROT %in% gold_823_clean_prots)
s3_intersect <- s3_markers %>% filter(UNIPROT %in% gold_823_clean_prots)

top5_s1 <- head(s1_intersect$UNIPROT, 5)
top5_s2 <- head(s2_intersect$UNIPROT, 5)
top5_s3 <- head(s3_intersect$UNIPROT, 5)
top15_prots <- c(top5_s1, top5_s2, top5_s3)

raw_mat_top15 <- prot_mat[top15_prots, sub_samps]
mat_scaled_top15 <- t(scale(t(raw_mat_top15)))
mat_capped_top15 <- pmax(pmin(mat_scaled_top15, 2.5), -2.5)
rownames(mat_capped_top15) <- lookup_global[top15_prots]

col_fun_heatmap <- colorRamp2(
  c(-2.5, -1.2, 0, 1.2, 2.5),
  c("#3A4866", "#92A1C2", "#FAF8F8", "#D77C91", "#A83853")
)

col_split <- factor(sub_clin$Subtype, levels = c("Subtype1", "Subtype2", "Subtype3"))
row_split <- factor(rep(c("Subtype 1 Drivers", "Subtype 2 Drivers", "Subtype 3 Drivers"), each = 5),
                    levels = c("Subtype 1 Drivers", "Subtype 2 Drivers", "Subtype 3 Drivers"))

top_ann <- HeatmapAnnotation(
  Subtype = col_split,
  col = list(Subtype = c("Subtype1" = "#d73027", "Subtype2" = "#3182bd", "Subtype3" = "#fdae61")),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  simple_anno_size = unit(3.5, "mm"),
  border = FALSE
)

pdf(file.path(output_dir, "Figure4B_Subtype_Fingerprints_Top15.pdf"), width = 10.0, height = 4.8)
ht_top15 <- Heatmap(
  mat_capped_top15,
  name = "Z-score",
  col = col_fun_heatmap,
  column_split = col_split,
  row_split = row_split,
  cluster_columns = TRUE,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 10.5, fontface = "italic"),
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  top_annotation = top_ann,
  column_gap = unit(1.2, "mm"),
  row_gap = unit(1.5, "mm"),
  border = TRUE
)
draw(ht_top15, padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()

# ==============================================================================
# Step 4: Subtype Pathway Enrichment & Pruning (Figure 4C)
# ==============================================================================
cat("Running GO-BP enrichment with topological filtering for Figure 4C...\n")

bg_raw_ids <- unique(gsub("-.*|\\..*", "", rownames(prot_raw_full)))
universe_entrez <- suppressMessages(
  bitr(bg_raw_ids, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
)$ENTREZID %>% unique()

s1_entrez <- na.omit(suppressMessages(bitr(gsub("-.*|\\..*", "", s1_intersect$UNIPROT), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db))$ENTREZID) %>% unique()
s2_entrez <- na.omit(suppressMessages(bitr(gsub("-.*|\\..*", "", s2_intersect$UNIPROT), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db))$ENTREZID) %>% unique()
s3_entrez <- na.omit(suppressMessages(bitr(gsub("-.*|\\..*", "", s3_intersect$UNIPROT), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db))$ENTREZID) %>% unique()

gene_clusters <- list("S1" = head(s1_entrez, 120), "S2" = head(s2_entrez, 120), "S3" = head(s3_entrez, 120))

comp_bp <- compareCluster(
  geneCluster   = gene_clusters,
  fun           = "enrichGO",
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pvalueCutoff  = 1.0,
  minGSSize     = 10,
  maxGSSize     = 500,
  universe      = universe_entrez,
  readable      = TRUE
)

df_go_raw <- as.data.frame(comp_bp)

# Filter anatomical structure development ancestor nodes (GO:0048856 / GO:0040007)
bp_ancestors <- as.list(GOBPANCESTOR)
is_dev_node <- function(go_id) {
  anc <- bp_ancestors[[go_id]]
  if (is.null(anc)) return(FALSE)
  return("GO:0048856" %in% anc || "GO:0040007" %in% anc)
}

df_go_filtered <- df_go_raw %>% filter(!sapply(ID, is_dev_node))

# Prune overlapping terms per cluster (Overlap >= 0.60 | Jaccard >= 0.35)
prune_terms <- function(df_sub, max_terms = 6, overlap_cut = 0.60, jaccard_cut = 0.35) {
  df_sig <- df_sub %>% filter(pvalue < 0.05) %>% arrange(pvalue)
  if (nrow(df_sig) <= 1) return(df_sig)
  selected <- c(1)
  for (i in 2:nrow(df_sig)) {
    if (length(selected) >= max_terms) break
    genes_curr <- unlist(strsplit(df_sig$geneID[i], "/"))
    redundant <- FALSE
    for (s in selected) {
      genes_prev <- unlist(strsplit(df_sig$geneID[s], "/"))
      inter_len  <- length(intersect(genes_curr, genes_prev))
      union_len  <- length(union(genes_curr, genes_prev))
      min_len    <- min(length(genes_curr), length(genes_prev))
      if ((inter_len / min_len >= overlap_cut) || (inter_len / union_len >= jaccard_cut)) {
        redundant <- TRUE
        break
      }
    }
    if (!redundant) selected <- c(selected, i)
  }
  return(df_sig[selected, ])
}

pruned_list <- list()
for (cl in c("S1", "S2", "S3")) {
  sub_cl <- df_go_filtered %>% filter(Cluster == cl)
  if (nrow(sub_cl) > 0) pruned_list[[cl]] <- prune_terms(sub_cl, max_terms = 6)
}

top_pathways_df <- do.call(rbind, pruned_list) %>%
  separate(GeneRatio, into = c("n", "N"), sep = "/", remove = FALSE) %>%
  mutate(
    GeneRatio_num = as.numeric(n) / as.numeric(N),
    Log10_P = -log10(pvalue),
    Description_Clean = str_wrap(str_to_title(Description), width = 36)
  )

# Build plotting grid for Figure 4C
plot_data_4c <- expand.grid(
  Description_Clean = unique(top_pathways_df$Description_Clean),
  Cluster = c("S1", "S2", "S3"),
  stringsAsFactors = FALSE
) %>%
  left_join(top_pathways_df %>% dplyr::select(Description_Clean, Cluster, Log10_P, GeneRatio_num), by = c("Description_Clean", "Cluster")) %>%
  mutate(
    Log10_P = replace_na(Log10_P, 0),
    GeneRatio_num = replace_na(GeneRatio_num, 0),
    Cluster = factor(Cluster, levels = c("S1", "S2", "S3"))
  )

y_order_4c <- top_pathways_df %>% arrange(Cluster, desc(Log10_P)) %>% pull(Description_Clean) %>% unique() %>% rev()
plot_data_4c$Description_Clean <- factor(plot_data_4c$Description_Clean, levels = y_order_4c)

p_fig4c <- ggplot(plot_data_4c, aes(x = Cluster, y = Description_Clean)) +
  geom_point(data = filter(plot_data_4c, Log10_P > 0),
             aes(size = GeneRatio_num, fill = Log10_P), shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradientn(colors = c("#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026"),
                       name = expression(bold(-log[10](P-value))), limits = c(1.3, 5.0), oob = scales::squish) +
  scale_size_continuous(range = c(3.0, 7.5), name = "Gene Ratio") +
  labs(title = "Core Mechanistic Divergence across Molecular Endotypes", x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12.5),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 9.5, lineheight = 0.85),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey88", linetype = "dashed"))

ggsave(file.path(output_dir, "Figure4C_Subtype_Mechanisms_Dotplot.pdf"), p_fig4c, width = 8.8, height = 7.5, device = cairo_pdf)

# ==============================================================================
# Step 5: Cellular Microenvironment Deconvolution (Figure 4D, 4E & Table ST17)
# ==============================================================================
cat("Performing PanglaoDB cellular microenvironment deconvolution...\n")

panglao_file <- find_file_smart("PanglaoDB_markers.*\\.tsv(\\.gz)?$", file.path(data_dir, "PanglaoDB_markers_27_Mar_2020.tsv.gz"))
if (!file.exists(panglao_file)) stop("PanglaoDB reference file not found in data directory.")

panglao_db <- fread(panglao_file, sep = "\t", data.table = FALSE) %>% filter(str_detect(species, "Hs"))

# Convert proteome matrix to Gene Symbols
prot_symbol_df <- as.data.frame(prot_mat) %>%
  mutate(SYMBOL = lookup_global[rownames(prot_mat)]) %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), ~mean(., na.rm = TRUE))) %>%
  column_to_rownames("SYMBOL")
prot_symbol_mat <- as.matrix(prot_symbol_df)

cell_types_all <- unique(panglao_db$`cell type`)
cell_signatures <- list()
for (ct in cell_types_all) {
  g_set <- panglao_db %>% filter(`cell type` == ct) %>% pull(`official gene symbol`) %>% unique()
  valid_g <- intersect(g_set, rownames(prot_symbol_mat))
  if (length(valid_g) >= 3) cell_signatures[[ct]] <- valid_g
}

param_cells <- ssgseaParam(exprData = prot_symbol_mat, geneSets = cell_signatures, minSize = 3)
gsva_cells  <- gsva(param_cells)
score_cells_df <- as.data.frame(t(gsva_cells), check.names = FALSE) %>%
  rownames_to_column("sample") %>%
  inner_join(sub_clin %>% dplyr::select(sample, Subtype), by = "sample")

# Export Supplementary Table ST17 across all 126 lineages
table_st17_stats <- map_df(names(cell_signatures), function(ct) {
  s_vec <- score_cells_df[[ct]]
  sub_vec <- score_cells_df$Subtype
  kw_res <- kruskal.test(s_vec ~ sub_vec)
  data.frame(
    Cell_Type = ct,
    Gene_Count = length(cell_signatures[[ct]]),
    Subtype1_Mean = mean(s_vec[sub_vec == "Subtype1"], na.rm = TRUE),
    Subtype2_Mean = mean(s_vec[sub_vec == "Subtype2"], na.rm = TRUE),
    Subtype3_Mean = mean(s_vec[sub_vec == "Subtype3"], na.rm = TRUE),
    Kruskal_Wallis_P = kw_res$p.value
  )
}) %>%
  mutate(Global_FDR = p.adjust(Kruskal_Wallis_P, method = "BH")) %>%
  arrange(Kruskal_Wallis_P)
write.csv(table_st17_stats, file.path(output_dir, "Supplementary_Table_17_Cellular_Deconvolution_Full.csv"), row.names = FALSE)

# Generate Figure 4D (CNS Lineages) & Figure 4E (Peripheral Lineages)
subtype_palette <- c("Subtype1" = "#d73027", "Subtype2" = "#4575b4", "Subtype3" = "#fdae61")
pairwise_comparisons <- list(c("Subtype1", "Subtype2"), c("Subtype2", "Subtype3"), c("Subtype1", "Subtype3"))

plot_cell_violins <- function(cell_vec, title_text) {
  valid_c <- intersect(cell_vec, colnames(score_cells_df))
  df_plot <- score_cells_df %>%
    dplyr::select(sample, Subtype, all_of(valid_c)) %>%
    pivot_longer(cols = all_of(valid_c), names_to = "Cell_Type", values_to = "Score")
  df_plot$Cell_Type <- factor(df_plot$Cell_Type, levels = cell_vec)
  
  ggplot(df_plot, aes(x = Subtype, y = Score, fill = Subtype)) +
    geom_violin(trim = FALSE, alpha = 0.45, color = NA, scale = "width", width = 0.85) +
    geom_boxplot(width = 0.16, fill = "white", outlier.shape = NA, color = "grey25", linewidth = 0.55) +
    facet_wrap(~ Cell_Type, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = subtype_palette) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
    stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.signif",
                       step.increase = 0.08, tip.length = 0.01, size = 3.6) +
    labs(title = title_text, x = NULL, y = "Deconvoluted Score (ssGSEA)") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 11.5, hjust = 0.5),
          strip.background = element_rect(fill = "grey93", color = "black", linewidth = 0.6),
          strip.text = element_text(face = "bold", size = 9.5),
          axis.text.x = element_text(face = "bold", size = 9.5, angle = 20, hjust = 1))
}

cns_cells_main    <- c("Neurons", "Astrocytes", "Microglia", "Oligodendrocytes")
periph_cells_main <- c("Endothelial cells", "Monocytes", "Macrophages", "T cells")

p_fig4d <- plot_cell_violins(cns_cells_main, "CNS Cellular Microenvironment")
p_fig4e <- plot_cell_violins(periph_cells_main, "Peripheral Immune & Vascular Lineages")

# ==============================================================================
# Step 6: Multivariable Partial Correlation Heatmap (Figure 4F)
# ==============================================================================
cat("Computing 8-covariate partial correlation for Figure 4F...\n")

clinical_traits <- c("UPDRS_III", "NMSS", "MMSE", "HAMA", "HAMD", "PDSS", "UPSIT")
actual_traits   <- intersect(clinical_traits, colnames(sub_clin))

trait_labels_display <- c(
  "UPDRS_III" = "Motor (UPDRS-III)",
  "NMSS"      = "Non-motor (NMSS)",
  "MMSE"      = "Cognition (MMSE)",
  "HAMA"      = "Anxiety (HAMA)",
  "HAMD"      = "Depression (HAMD)",
  "PDSS"      = "Sleep (PDSS)",
  "UPSIT"     = "Olfaction (UPSIT)"
)

core_heatmap_cells <- c(cns_cells_main, periph_cells_main)

calc_pcor_matrix <- function(mat_features, feat_names, clin_data, trait_names) {
  r_mat <- matrix(NA_real_, nrow = length(feat_names), ncol = length(trait_names),
                  dimnames = list(feat_names, trait_names))
  p_mat <- matrix(NA_real_, nrow = length(feat_names), ncol = length(trait_names),
                  dimnames = list(feat_names, trait_names))
  
  df_feat <- as.data.frame(t(mat_features[feat_names, , drop = FALSE])) %>% rownames_to_column("sample")
  df_merged <- df_feat %>% inner_join(clin_data, by = "sample")
  
  for (f in feat_names) {
    for (tr in trait_names) {
      sub_pcor <- df_merged %>%
        dplyr::select(feat_val = all_of(f), trait_val = all_of(tr), all_of(BASE_COVARIATES)) %>%
        drop_na()
      if (nrow(sub_pcor) >= 30) {
        cov_sub <- sub_pcor[, BASE_COVARIATES]
        res <- tryCatch({ pcor.test(sub_pcor$feat_val, sub_pcor$trait_val, cov_sub, method = "spearman") }, error = function(e) NULL)
        if (!is.null(res)) {
          r_mat[f, tr] <- res$estimate
          p_mat[f, tr] <- res$p.value
        }
      }
    }
  }
  colnames(r_mat) <- trait_labels_display[colnames(r_mat)]
  colnames(p_mat) <- trait_labels_display[colnames(p_mat)]
  return(list(r = r_mat, p = p_mat))
}

pcor_cells_res <- calc_pcor_matrix(gsva_cells, core_heatmap_cells, sub_clin, actual_traits)
r_cells_mat <- pcor_cells_res$r
p_cells_mat <- pcor_cells_res$p

cell_split_tags <- factor(ifelse(core_heatmap_cells %in% cns_cells_main, "CNS Lineages", "Peripheral Lineages"),
                          levels = c("CNS Lineages", "Peripheral Lineages"))

col_pcor <- colorRamp2(c(-0.35, 0, 0.35), c("#313695", "#FAF8F8", "#A50026"))

pdf(file.path(output_dir, "Figure4F_Cellular_Clinical_Pcor_Heatmap.pdf"), width = 9.2, height = 6.2)
ht_cells_pcor <- Heatmap(
  r_cells_mat,
  name = "Partial r",
  col = col_pcor,
  row_split = cell_split_tags,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "plain"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45,
  row_title_gp = gpar(fontsize = 10.5, fontface = "bold"),
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(p_cells_mat[i, j]) && p_cells_mat[i, j] < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 12, fontface = "bold"))
    }
  }
)
draw(ht_cells_pcor, padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()

# ==============================================================================
# Step 7: Generate Supplementary Figure 6 (Panels A & B)
# ==============================================================================
cat("Generating Supplementary Figure 6 (Extended lineages & 17 Pathway Heatmap)...\n")

# Supp Fig 6A: 8 Extended cell lineages
supp_cns_cells    <- c("Purkinje neurons", "Ependymal cells", "Neural stem/precursor cells", "Schwann cells")
supp_periph_cells <- c("Pericytes", "Neutrophils", "Platelets", "B cells")
supp_all_cells    <- c(supp_cns_cells, supp_periph_cells)

p_supp6a_top <- plot_cell_violins(supp_cns_cells, "Extended CNS & Neuroglial Support Lineages")
p_supp6a_bot <- plot_cell_violins(supp_periph_cells, "Extended Peripheral Vascular & Immune Lineages")
supp_fig6a <- p_supp6a_top / p_supp6a_bot

# Supp Fig 6B: 17 Core subtype-specific biological pathways partial correlation heatmap
pathway_gene_sets <- list()
pathway_subtype_origin <- c()

for (i in 1:nrow(top_pathways_df)) {
  p_name <- str_to_title(top_pathways_df$Description[i])
  raw_g  <- unlist(strsplit(top_pathways_df$geneID[i], "/"))
  val_g  <- intersect(raw_g, rownames(prot_symbol_mat))
  if (length(val_g) >= 3) {
    pathway_gene_sets[[p_name]] <- val_g
    pathway_subtype_origin[p_name] <- top_pathways_df$Cluster[i]
  }
}

param_paths <- ssgseaParam(exprData = prot_symbol_mat, geneSets = pathway_gene_sets, minSize = 3)
gsva_paths  <- gsva(param_paths)
actual_pathways <- rownames(gsva_paths)

pcor_paths_res <- calc_pcor_matrix(gsva_paths, actual_pathways, sub_clin, actual_traits)
r_paths_mat <- pcor_paths_res$r
p_paths_mat <- pcor_paths_res$p

path_subtype_tags <- factor(paste0("Subtype ", gsub("S", "", pathway_subtype_origin[actual_pathways]), " Signatures"),
                            levels = c("Subtype 1 Signatures", "Subtype 2 Signatures", "Subtype 3 Signatures"))

pdf(file.path(output_dir, "Supplementary_Figure_6B_Pathway_Pcor_Heatmap.pdf"), width = 9.8, height = 7.5)
ht_supp6b <- Heatmap(
  r_paths_mat,
  name = "Partial r",
  col = col_pcor,
  row_split = path_subtype_tags,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 9.5, fontface = "plain"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45,
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(p_paths_mat[i, j]) && p_paths_mat[i, j] < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 12, fontface = "bold"))
    }
  }
)
draw(ht_supp6b, padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()

# Export combined Supplementary Figure 6A
ggsave(file.path(output_dir, "Supplementary_Figure_6A_Extended_Lineages.pdf"), supp_fig6a, width = 12.0, height = 7.0, device = cairo_pdf)
ggsave(file.path(output_dir, "Supplementary_Figure_6A_Extended_Lineages.png"), supp_fig6a, width = 12.0, height = 7.0, dpi = 300)

cat("Analysis complete. Figure 4, Supplementary Figure 6, and Tables ST16, ST17 saved to:", output_dir, "\n")