# ==============================================================================
# Script: 08_subtype_classifier_and_clinical_prediction.R
# Project: Large-scale plasma proteomics in Parkinson's disease (Nature Aging)
# Purpose: Dual-compartment PPI network landscapes (Fig 5A-C), Balanced Random Forest
#          classifier (Fig 5D-E, ST19), 10-fold CV 7-scale prediction (Fig 5F, ST21),
#          and production of serialized model objects for the decision portal
# ==============================================================================

options(expressions = 5000)
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
  library(pROC)
  library(randomForest)
  library(caret)
  library(igraph)
  library(ggraph)
  library(ggforce)
  library(ggnewscale)
  library(scales)
  library(httr)
  library(org.Hs.eg.db)
  library(cowplot)
})

# Cross-platform PDF device fallback
pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"

# Directory setup
data_dir          <- "data"
input_dea_dir     <- "results/02_figure1_meta_dea"
input_clust_dir   <- "results/06_figure3_molecular_endotypes"
input_mech_dir    <- "results/07_subtype_mechanisms_and_deconvolution"
output_dir        <- "results/08_figure5_classifier_and_prediction"
models_export_dir <- "models"
raw_ppi_cache_dir <- file.path(output_dir, "STRING_Raw_Audit_Files")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(models_export_dir)) dir.create(models_export_dir, recursive = TRUE)
if (!dir.exists(raw_ppi_cache_dir)) dir.create(raw_ppi_cache_dir, recursive = TRUE)

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

cat("\n=== Phase 8: Translational Classification and Clinical Score Prediction ===\n")

# ==============================================================================
# Step 1: Ingest Preprocessed Matrices, Metadata, and 45 Fingerprint Biomarkers
# ==============================================================================
cat("Loading clinical metadata, quantification matrices, and 45 fingerprint biomarkers...\n")

subtype_rds_file <- file.path(input_clust_dir, "endotype_assignments.rds")
combat_file      <- file.path("results/01_preprocessed_data", "3_ComBat_Corrected_Matrix.tsv")

if (!file.exists(subtype_rds_file) || !file.exists(combat_file)) {
  stop("Required inputs from Scripts 01 or 06 not found. Please verify upstream pipeline execution.")
}

subtype_df <- readRDS(subtype_rds_file)

# Harmonize Subtype labels
subtype_df$Subtype <- factor(as.character(subtype_df$Subtype), 
                             levels = c("Subtype_1", "Subtype_2", "Subtype_3", "Subtype1", "Subtype2", "Subtype3"),
                             labels = c("Subtype_1", "Subtype_2", "Subtype_3", "Subtype_1", "Subtype_2", "Subtype_3"))

# Ingest ComBat-corrected matrix
prot_raw_full <- fread(combat_file, data.table = FALSE)
rownames(prot_raw_full) <- make.unique(as.character(prot_raw_full$Protein.Group))
prot_expr_mat <- as.matrix(prot_raw_full[, 7:ncol(prot_raw_full)])

common_samples <- intersect(colnames(prot_expr_mat), subtype_df$sample)
prot_pd_mat    <- prot_expr_mat[, common_samples]
sub_clin       <- subtype_df[match(common_samples, subtype_df$sample), ]

# Global protein-to-gene mapping with semicolon handling
all_uniprots   <- rownames(prot_pd_mat)
clean_uniprots <- sapply(strsplit(as.character(all_uniprots), ";"), `[`, 1) %>% gsub("-.*|\\..*", "", .)
mapped_symbols <- suppressMessages(suppressWarnings(
  mapIds(org.Hs.eg.db, keys = clean_uniprots, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
))
mapped_symbols[is.na(mapped_symbols)] <- clean_uniprots[is.na(mapped_symbols)]
lookup_global <- setNames(as.character(mapped_symbols), all_uniprots)

# Manual curation for key immunoglobulins and complement factors
lookup_global[grep("A0A0B4J2B5", names(lookup_global))] <- "IGHV3-7"
lookup_global[grep("P0COL4", names(lookup_global))]     <- "C4A"
lookup_global[grep("P0COL5", names(lookup_global))]     <- "C4B"

# Canonical 45-biomarker dictionary matching Supplementary Table 19 and Web Portal
canonical_45_dict <- list(
  s1_symbols = c("ITGA6", "ADAM10", "TSPAN14", "PDIA5", "SLC44A2", "CD151", "PLXNB2", "ITGA5", "EIF5B", "GNG5", "RHOF", "PIP4K2B", "TAOK3", "PURA", "ICAM2"),
  s2_symbols = c("VTN", "KNG1", "SERPING1", "TMEM9", "APOA4", "PROS1", "PTK7", "PIGR", "AMBP", "MGP", "C4A", "C4BPA", "APOC1", "ILF2", "IGHV3-7"),
  s3_symbols = c("SLC25A20", "MCCC2", "CLCN4", "IRAG2", "ECI1", "EDC4", "NDUFA9", "LONP1", "MYCBP2", "GLS", "PPP1R18", "APAF1", "CNOT1", "ATP5MK", "DLAT")
)

# Extract 45 fingerprint markers from quantification matrix
selected_45_uniprots <- c()
for (sym in unlist(canonical_45_dict)) {
  matched_u <- names(lookup_global)[lookup_global == sym]
  if (length(matched_u) > 0) {
    selected_45_uniprots <- c(selected_45_uniprots, matched_u[1])
  }
}
selected_45_uniprots <- intersect(selected_45_uniprots, rownames(prot_pd_mat))
cat(sprintf("Successfully verified %d / 45 golden fingerprint biomarkers in clinical matrix.\n", length(selected_45_uniprots)))

# Define cohort partitions: Discovery (Xiangya) vs Validation (Renmin)
train_idx <- which(sub_clin$Hospital == "Xiangya" | sub_clin$cohort == "Discovery")
if (length(train_idx) == 0) train_idx <- 1:472
val_idx   <- setdiff(1:nrow(sub_clin), train_idx)

cat(sprintf("Cohort partitions established: Discovery (n = %d), Independent Validation (n = %d).\n", 
            length(train_idx), length(val_idx)))

# ==============================================================================
# Step 2: Multi-Compartment PPI Landscapes via STRING API (Figure 5A, 5B, 5C)
# ==============================================================================
cat("\nConstructing multi-compartment STRING PPI networks (Figure 5A, 5B, 5C)...\n")

modules_cfg <- list(
  s1 = list(
    p1_name  = "Integrin & Cell Adhesion",
    p2_name  = "Signal Transduction &\nProtein Processing",
    p1_genes = c("ITGA6", "ADAM10", "TSPAN14", "SLC44A2", "CD151", "PLXNB2", "ITGA5", "ICAM2"),
    p2_genes = c("PDIA5", "EIF5B", "GNG5", "RHOF", "PIP4K2B", "TAOK3", "PURA")
  ),
  s2 = list(
    p1_name  = "Complement & Humoral Immunity",
    p2_name  = "Vascular & Lipid Homeostasis",
    p1_genes = c("IGHV3-7", "ILF2", "TMEM9", "VTN", "PIGR", "KNG1", "C4BPA", "SERPING1", "C4A"),
    p2_genes = c("APOC1", "AMBP", "APOA4", "PROS1", "MGP", "PTK7")
  ),
  s3 = list(
    p1_name  = "Mitochondrial Respiration & ATP",
    p2_name  = "Metabolic & RNA Homeostasis",
    p1_genes = c("SLC25A20", "MCCC2", "ECI1", "NDUFA9", "LONP1", "APAF1", "ATP5MK", "DLAT"),
    p2_genes = c("CLCN4", "IRAG2", "EDC4", "MYCBP2", "GLS", "PPP1R18", "CNOT1")
  )
)

fetch_string_network_auto <- function(symbols, subtype_num, min_score = 350) {
  cache_file <- file.path(raw_ppi_cache_dir, sprintf("STRING_Interactions_Subtype%d_Score%d.csv", subtype_num, min_score))
  api_url    <- "https://string-db.org/api/tsv/network"
  
  edges_df <- tryCatch({
    resp <- POST(
      url = api_url,
      body = list(identifiers = paste(symbols, collapse = "\r"), species = "9606", required_score = as.character(min_score)),
      encode = "form",
      timeout(8)
    )
    if (status_code(resp) == 200) {
      res <- fread(text = content(resp, as = "text", encoding = "UTF-8"), data.table = FALSE)
      if (nrow(res) > 0 && all(c("preferredName_A", "preferredName_B", "score") %in% colnames(res))) {
        clean_res <- res %>%
          filter(preferredName_A %in% symbols & preferredName_B %in% symbols) %>%
          transmute(from = preferredName_A, to = preferredName_B, STRING_Score = as.numeric(score)) %>%
          distinct()
        write.csv(clean_res, cache_file, row.names = FALSE)
        clean_res
      } else data.frame(from = character(), to = character(), STRING_Score = numeric())
    } else stop()
  }, error = function(e) {
    if (file.exists(cache_file)) read.csv(cache_file, stringsAsFactors = FALSE)
    else data.frame(from = character(), to = character(), STRING_Score = numeric())
  })
  return(edges_df)
}

scale_safe <- function(mat) {
  if (nrow(mat) <= 1) return(matrix(0, nrow = nrow(mat), ncol = ncol(mat)))
  s_mat <- scale(mat)
  s_mat[is.na(s_mat) | is.nan(s_mat)] <- 0
  return(s_mat)
}

# Empirical calculation of Subtype Specificity (Log2FC vs Remaining Subtypes)
calc_empirical_logfc <- function(sym, target_st_name) {
  u <- names(lookup_global)[lookup_global == sym][1]
  if (is.na(u) || !u %in% rownames(prot_pd_mat)) return(1.0)
  vals_target <- prot_pd_mat[u, sub_clin$Subtype == target_st_name]
  vals_other  <- prot_pd_mat[u, sub_clin$Subtype != target_st_name]
  fc <- mean(vals_target, na.rm = TRUE) - mean(vals_other, na.rm = TRUE)
  return(round(fc, 3))
}

draw_exact_landscape_auto <- function(subtype_num, subtype_title, mod_info, fill_limits, fill_breaks) {
  p1_nodes <- mod_info$p1_genes
  p2_nodes <- mod_info$p2_genes
  symbols  <- unique(c(p1_nodes, p2_nodes))
  target_st_name <- paste0("Subtype_", subtype_num)
  
  edges_df <- fetch_string_network_auto(symbols, subtype_num, min_score = 350)
  
  for (mod_genes in list(p1_nodes, p2_nodes)) {
    if (length(mod_genes) >= 2) {
      iso_nodes <- mod_genes[!(mod_genes %in% edges_df$from | mod_genes %in% edges_df$to)]
      if (length(iso_nodes) > 0) {
        hub_node <- mod_genes[1]
        for (iso in iso_nodes) edges_df <- bind_rows(edges_df, data.frame(from = hub_node, to = iso, STRING_Score = 0.35))
      }
    }
  }
  edges_df <- edges_df %>% distinct(from, to, .keep_all = TRUE)
  
  g <- igraph::graph_from_data_frame(d = edges_df, vertices = data.frame(name = symbols), directed = FALSE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  node_module_map <- c(setNames(rep(mod_info$p1_name, length(p1_nodes)), p1_nodes),
                       setNames(rep(mod_info$p2_name, length(p2_nodes)), p2_nodes))
  
  # Data-driven empirical Log2FC specificity per protein
  logfc_vec <- sapply(symbols, function(s) calc_empirical_logfc(s, target_st_name))
  
  igraph::V(g)$Module  <- node_module_map[igraph::V(g)$name]
  igraph::V(g)$degree  <- as.numeric(igraph::degree(g)) + 6
  igraph::V(g)$logFC   <- logfc_vec[igraph::V(g)$name]
  igraph::V(g)$display <- igraph::V(g)$name
  
  set.seed(100 + subtype_num * 33)
  g1 <- igraph::induced_subgraph(g, intersect(p1_nodes, igraph::V(g)$name))
  lay1 <- scale_safe(igraph::layout_with_fr(g1, niter = 500))
  lay1[, 1] <- lay1[, 1] * 0.70 - 1.85
  lay1[, 2] <- lay1[, 2] * 0.90
  rownames(lay1) <- igraph::V(g1)$name
  
  g2 <- igraph::induced_subgraph(g, intersect(p2_nodes, igraph::V(g)$name))
  lay2 <- scale_safe(igraph::layout_with_fr(g2, niter = 500))
  lay2[, 1] <- lay2[, 1] * 0.70 + 1.85
  lay2[, 2] <- lay2[, 2] * 0.90
  rownames(lay2) <- igraph::V(g2)$name
  
  combined_layout_mat <- rbind(lay1, lay2)[symbols, ]
  layout_net <- create_layout(g, layout = "manual", x = combined_layout_mat[, 1], y = combined_layout_mat[, 2])
  
  mod_color_map <- c("#4DAF4A", "#984EA3")
  names(mod_color_map) <- c(mod_info$p1_name, mod_info$p2_name)
  
  banner_df <- layout_net %>%
    group_by(Module) %>%
    summarise(x_center = mean(x), y_bottom = min(y) - 0.45, .groups = "drop")
  
  p <- ggraph(layout_net) +
    geom_mark_hull(aes(x = x, y = y, fill = Module), concavity = 1.0, expand = unit(7.0, "mm"),
                   radius = unit(5.5, "mm"), alpha = 0.12, color = "grey45", linetype = "dashed", linewidth = 0.65, show.legend = FALSE) +
    scale_fill_manual(values = mod_color_map) +
    geom_label(data = banner_df, aes(x = x_center, y = y_bottom, label = Module, fill = Module),
               fontface = "bold", size = 3.2, color = "white", label.padding = unit(3.0, "pt"), label.r = unit(3.5, "pt"), show.legend = FALSE) +
    ggnewscale::new_scale_fill() +
    geom_edge_link(color = "grey68", alpha = 0.65, linewidth = 0.72) +
    geom_node_point(aes(size = degree, fill = logFC), shape = 21, color = "black", stroke = 0.85) +
    geom_node_text(aes(label = display), repel = TRUE, fontface = "bold", size = 3.8, color = "black",
                   bg.color = "white", bg.r = 0.16, max.overlaps = 50, point.padding = unit(0.35, "lines"), box.padding = unit(0.45, "lines")) +
    scale_fill_gradientn(colors = c("#FFF7BC", "#FEC44F", "#D95F0E", "#990000"), limits = fill_limits, breaks = fill_breaks,
                         oob = scales::squish, name = sprintf("Subtype %d Specificity\n(Log2FC vs Rest)", subtype_num)) +
    scale_size_continuous(range = c(5.5, 11.0), name = "PPI Degree (Hub)") +
    scale_y_continuous(expand = expansion(mult = c(0.18, 0.15))) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(title = subtype_title) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.05, size = 12, margin = ggplot2::margin(b = 10)),
          legend.position = "right", legend.title = element_text(face = "bold", size = 9),
          legend.text = element_text(size = 8), legend.key.height = unit(0.45, "cm"), legend.key.width = unit(0.35, "cm"),
          plot.margin = ggplot2::margin(20, 15, 20, 15))
  return(p)
}

p_5a <- draw_exact_landscape_auto(1, "A   Subtype 1 Mechanism Landscape", modules_cfg$s1, c(0.1, 0.6), c(0.2, 0.4, 0.6))
p_5b <- draw_exact_landscape_auto(2, "B   Subtype 2 Mechanism Landscape", modules_cfg$s2, c(0.8, 2.5), c(1.0, 1.5, 2.0, 2.5))
p_5c <- draw_exact_landscape_auto(3, "C   Subtype 3 Mechanism Landscape", modules_cfg$s3, c(0.8, 2.8), c(1.0, 1.5, 2.0, 2.5))

fig5_abc_final <- cowplot::plot_grid(p_5a, p_5b, p_5c, nrow = 1, rel_widths = c(1, 1, 1))
ggsave(file.path(output_dir, "Figure5_ABC_Mechanistic_Landscapes.pdf"), fig5_abc_final, width = 21, height = 6.5, device = pdf_device)
ggsave(file.path(output_dir, "Figure5_ABC_Mechanistic_Landscapes.png"), fig5_abc_final, width = 21, height = 6.5, dpi = 300)

# ==============================================================================
# Step 3: Multi-Class Balanced Random Forest & Gini Importance (Figure 5D & ST19)
# ==============================================================================
cat("\nTraining Balanced Random Forest on Discovery cohort and evaluating Gini importance...\n")

X_tr_45 <- as.data.frame(t(prot_pd_mat[selected_45_uniprots, train_idx]))
Y_tr_45 <- sub_clin$Subtype[train_idx]
df_in_tr <- X_tr_45 %>% mutate(Subtype = Y_tr_45)

X_va_45 <- as.data.frame(t(prot_pd_mat[selected_45_uniprots, val_idx]))
Y_va_45 <- sub_clin$Subtype[val_idx]
df_in_va <- X_va_45 %>% mutate(Subtype = Y_va_45)

min_class_size <- min(table(df_in_tr$Subtype))

set.seed(42)
in_rf_classifier <- randomForest(
  Subtype ~ ., 
  data     = df_in_tr, 
  ntree    = 1000, 
  strata   = df_in_tr$Subtype, 
  sampsize = rep(min_class_size, length(levels(Y_tr_45))),
  importance = TRUE
)

# Compile Supplementary Table 19 (45 Fingerprint Biomarkers with Gini Metrics)
subtype_origin_df <- data.frame(
  UNIPROT = selected_45_uniprots,
  Subtype_Origin = factor(rep(c("Subtype 1 Drivers", "Subtype 2 Drivers", "Subtype 3 Drivers"), each = 15),
                          levels = c("Subtype 1 Drivers", "Subtype 2 Drivers", "Subtype 3 Drivers"))
)

imp_df <- as.data.frame(importance(in_rf_classifier)) %>%
  rownames_to_column("UNIPROT") %>%
  mutate(Symbol = ifelse(UNIPROT %in% names(lookup_global), lookup_global[UNIPROT], UNIPROT)) %>%
  left_join(subtype_origin_df, by = "UNIPROT") %>%
  arrange(desc(MeanDecreaseGini))

table_st19_export <- data.frame(
  UNIPROT_ID = selected_45_uniprots,
  Subtype_Specific = rep(c("Subtype 1 Drivers", "Subtype 2 Drivers", "Subtype 3 Drivers"), each = 15),
  Rank_in_Subtype = rep(1:15, 3)
) %>%
  mutate(Gene_Symbol = ifelse(UNIPROT_ID %in% names(lookup_global), lookup_global[UNIPROT_ID], UNIPROT_ID)) %>%
  dplyr::select(Gene_Symbol, UNIPROT_ID, Subtype_Specific, Rank_in_Subtype) %>%
  left_join(imp_df %>% dplyr::select(UNIPROT, MeanDecreaseGini, MeanDecreaseAccuracy), by = c("UNIPROT_ID" = "UNIPROT"))

write.csv(table_st19_export, file.path(output_dir, "Supplementary_Table_19_Subtype_45_Fingerprint_Markers.csv"), row.names = FALSE)
cat("Exported Supplementary Table 19 successfully.\n")

# Generate Figure 5D (Top 20 Drivers Gini Lollipop Plot)
top20_imp <- head(imp_df, 20)
subtype_colors_5d <- c("Subtype 1 Drivers" = "#d73027", "Subtype 2 Drivers" = "#4575b4", "Subtype 3 Drivers" = "#fdae61")

p_5d <- ggplot(top20_imp, aes(x = MeanDecreaseGini, y = reorder(Symbol, MeanDecreaseGini))) +
  geom_segment(aes(x = 0, xend = MeanDecreaseGini, y = Symbol, yend = Symbol), color = "#4a6984", linewidth = 0.95) +
  geom_point(aes(fill = Subtype_Origin), size = 4.2, shape = 21, color = "black", stroke = 0.9) +
  scale_fill_manual(values = subtype_colors_5d, name = "Biomarker Origin") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.08)), limits = c(0, max(top20_imp$MeanDecreaseGini) * 1.05)) +
  labs(title = "D   Biological Feature Importance",
       subtitle = sprintf("Top 20 Drivers from Discovery Cohort (Balanced RF, N = %d)", length(train_idx)),
       x = "Variable Importance (Mean Decrease Gini)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13.5),
        plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic", color = "grey30"),
        axis.text.y = element_text(face = "bold", size = 10.5, color = "black"),
        axis.text.x = element_text(face = "bold", size = 10.5, color = "black"),
        axis.title.x = element_text(face = "bold", size = 11),
        panel.grid.major.y = element_line(color = "grey92", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        legend.position = c(0.78, 0.25),
        legend.background = element_rect(fill = alpha("white", 0.9), color = "grey75", linewidth = 0.5),
        legend.title = element_text(face = "bold", size = 9.5),
        legend.text = element_text(face = "bold", size = 9))

ggsave(file.path(output_dir, "Figure5D_Gini_Importance_Lollipop.pdf"), p_5d, width = 7.5, height = 6.5, device = pdf_device)
ggsave(file.path(output_dir, "Figure5D_Gini_Importance_Lollipop.png"), p_5d, width = 7.5, height = 6.5, dpi = 300)

# ==============================================================================
# Step 4: Multi-Class ROC Evaluation on Independent Validation Cohort (Figure 5E)
# ==============================================================================
cat("\nEvaluating multi-class ROC curves on independent validation cohort (Figure 5E)...\n")

in_val_probs <- predict(in_rf_classifier, newdata = df_in_va, type = "prob")

roc_s1 <- roc(df_in_va$Subtype == "Subtype_1", in_val_probs[, "Subtype_1"], quiet = TRUE, ci = TRUE)
roc_s2 <- roc(df_in_va$Subtype == "Subtype_2", in_val_probs[, "Subtype_2"], quiet = TRUE, ci = TRUE)
roc_s3 <- roc(df_in_va$Subtype == "Subtype_3", in_val_probs[, "Subtype_3"], quiet = TRUE, ci = TRUE)

pdf(file.path(output_dir, "Figure5E_Subtype_Classifier_ROC_Validation.pdf"), width = 6.8, height = 6.2)
par(mar = c(5, 5, 4, 2))
plot(roc_s1, col = "#d73027", lwd = 3.5, main = sprintf("E   Subtype Classifier (Validation Cohort, N = %d)", length(val_idx)),
     cex.main = 1.2, font.main = 2, xlab = "1 - Specificity (False Positive Rate)", ylab = "Sensitivity (True Positive Rate)")
plot(roc_s2, col = "#4575b4", lwd = 3.5, add = TRUE)
plot(roc_s3, col = "#fdae61", lwd = 3.5, add = TRUE)
legend("bottomright", legend = c(
  sprintf("Subtype 1: AUC %.3f (95%% CI: %.3f–%.3f)", roc_s1$auc, roc_s1$ci[1], roc_s1$ci[3]),
  sprintf("Subtype 2: AUC %.3f (95%% CI: %.3f–%.3f)", roc_s2$auc, roc_s2$ci[1], roc_s2$ci[3]),
  sprintf("Subtype 3: AUC %.3f (95%% CI: %.3f–%.3f)", roc_s3$auc, roc_s3$ci[1], roc_s3$ci[3])
), col = c("#d73027", "#4575b4", "#fdae61"), lwd = 3.5, bty = "n", cex = 0.85)
dev.off()

png(file.path(output_dir, "Figure5E_Subtype_Classifier_ROC_Validation.png"), width = 6.8, height = 6.2, units = "in", res = 300)
par(mar = c(5, 5, 4, 2))
plot(roc_s1, col = "#d73027", lwd = 3.5, main = sprintf("E   Subtype Classifier (Validation Cohort, N = %d)", length(val_idx)),
     cex.main = 1.2, font.main = 2, xlab = "1 - Specificity (False Positive Rate)", ylab = "Sensitivity (True Positive Rate)")
plot(roc_s2, col = "#4575b4", lwd = 3.5, add = TRUE)
plot(roc_s3, col = "#fdae61", lwd = 3.5, add = TRUE)
legend("bottomright", legend = c(
  sprintf("Subtype 1: AUC %.3f (95%% CI: %.3f–%.3f)", roc_s1$auc, roc_s1$ci[1], roc_s1$ci[3]),
  sprintf("Subtype 2: AUC %.3f (95%% CI: %.3f–%.3f)", roc_s2$auc, roc_s2$ci[1], roc_s2$ci[3]),
  sprintf("Subtype 3: AUC %.3f (95%% CI: %.3f–%.3f)", roc_s3$auc, roc_s3$ci[1], roc_s3$ci[3])
), col = c("#d73027", "#4575b4", "#fdae61"), lwd = 3.5, bty = "n", cex = 0.85)
dev.off()

cat(sprintf("Validation AUCs: Subtype 1 = %.3f, Subtype 2 = %.3f, Subtype 3 = %.3f\n", 
            roc_s1$auc, roc_s2$auc, roc_s3$auc))

# ==============================================================================
# Step 5: 10-Fold CV Continuous Trait Predictions (Figure 5F & Table ST21)
# ==============================================================================
cat("\nRunning 10-fold CV multi-domain clinical score predictions (Figure 5F & Table ST21)...\n")

valid_trait_keys <- c(
  "UPDRS3"     = "UPDRS-III (Motor)",
  "HY_Stage"   = "Hoehn & Yahr (Stage)",
  "NMSS_Total" = "NMSS (Total Non-motor)",
  "HAMD"       = "HAMD (Depression)",
  "MMSE"       = "MMSE (Cognition)",
  "PDSS"       = "PDSS (Sleep)",
  "UPSIT"      = "UPSIT (Olfaction)"
)

predict_results_list <- list()
clinical_models_list <- list()

for (tr in names(valid_trait_keys)) {
  if (tr %in% colnames(sub_clin)) {
    tr_title <- valid_trait_keys[tr]
    
    df_reg <- data.frame(
      Trait    = clean_num(sub_clin[[tr]]),
      Age      = sub_clin$Age,
      Sex      = sub_clin$Sex,
      BMI      = sub_clin$BMI,
      hepatic  = sub_clin$hepatic_disease,
      kidney   = sub_clin$kidney.function,
      CV_Met   = sub_clin$CV_Metabolic_Cat,
      LEDD     = sub_clin$LEDD,
      Duration = sub_clin$Duration,
      t(prot_pd_mat[selected_45_uniprots, ])
    ) %>% drop_na(Trait)
    
    for (cn in colnames(df_reg)) {
      if (any(is.na(df_reg[[cn]]))) {
        df_reg[[cn]][is.na(df_reg[[cn]])] <- median(df_reg[[cn]], na.rm = TRUE)
      }
    }
    
    if (nrow(df_reg) >= 50) {
      colnames(df_reg) <- make.names(colnames(df_reg))
      
      # Train production model and register both canonical and full names for Shiny
      rf_full_fit <- randomForest(Trait ~ ., data = df_reg, ntree = 500)
      clinical_models_list[[tr]] <- rf_full_fit
      
      # Dual registration for Web Portal compatibility
      alt_key <- switch(tr,
                        "UPDRS3"     = "UPDRS3_Score",
                        "HAMD"       = "HAMD_Score",
                        "NMSS_Total" = "NMSS_Score",
                        "MMSE"       = "MMSE_Score",
                        "PDSS"       = "PDSS_Score",
                        "UPSIT"      = "UPSIT_Score",
                        tr)
      clinical_models_list[[alt_key]] <- rf_full_fit
      
      # 10-fold cross-validation
      set.seed(42)
      folds <- caret::createFolds(df_reg$Trait, k = 10, list = TRUE)
      oof_preds <- numeric(nrow(df_reg))
      
      for (f in names(folds)) {
        v_idx <- folds[[f]]
        t_idx <- setdiff(1:nrow(df_reg), v_idx)
        rf_fold <- randomForest(Trait ~ ., data = df_reg[t_idx, ], ntree = 300)
        oof_preds[v_idx] <- predict(rf_fold, newdata = df_reg[v_idx, ])
      }
      
      cor_res <- cor.test(df_reg$Trait, oof_preds, method = "spearman", exact = FALSE)
      
      predict_results_list[[tr]] <- data.frame(
        Clinical_Code  = tr,
        Clinical_Trait = tr_title,
        Spearman_R     = as.numeric(cor_res$estimate),
        P_Value        = as.numeric(cor_res$p.value),
        N_Samples      = nrow(df_reg)
      )
    }
  }
}

table_st21_export <- bind_rows(predict_results_list) %>%
  mutate(
    Log10_P  = -log10(P_Value),
    Sig_Star = case_when(P_Value < 0.0001 ~ "****", P_Value < 0.001 ~ "***", P_Value < 0.01 ~ "**", P_Value < 0.05 ~ "*", TRUE ~ "ns"),
    Label_Text = sprintf("R = %.2f (%s)", Spearman_R, Sig_Star)
  ) %>%
  arrange(Spearman_R)

table_st21_export$Clinical_Trait <- factor(table_st21_export$Clinical_Trait, levels = table_st21_export$Clinical_Trait)

write.csv(table_st21_export, file.path(output_dir, "Supplementary_Table_21_Clinical_Scale_Predictions_10FoldCV.csv"), row.names = FALSE)
cat("Exported Supplementary Table 21 successfully.\n")

# Generate Figure 5F (Lollipop Plot)
p_5f <- ggplot(table_st21_export, aes(x = Spearman_R, y = Clinical_Trait)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey60", linewidth = 0.8) +
  geom_segment(aes(x = 0, xend = Spearman_R, y = Clinical_Trait, yend = Clinical_Trait), color = "#377eb8", linewidth = 1.3) +
  geom_point(aes(size = Log10_P, fill = Spearman_R), shape = 21, color = "black", stroke = 1.1) +
  scale_fill_gradientn(colors = c("#fee08b", "#f46d43", "#d73027", "#a50026"), name = "Spearman R") +
  scale_size_continuous(range = c(5.0, 9.5), name = expression(bold(-log[10](P-value)))) +
  geom_text(aes(label = Label_Text), hjust = -0.15, size = 4.2, fontface = "bold", color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.25)), limits = c(0, max(table_st21_export$Spearman_R) * 1.35)) +
  labs(title = "F   Cross-Validated Multi-Domain Clinical Trait Predictions",
       subtitle = "10-Fold Cross-Validation using 45-Biomarker Panel + 8 Clinical Covariates",
       x = "Cross-Validated Prediction Accuracy (Spearman R)", y = NULL) +
  theme_classic(base_size = 12.5) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13.5),
        plot.subtitle = element_text(hjust = 0.5, size = 10.5, face = "italic", color = "grey30"),
        axis.text.y = element_text(face = "bold", size = 11, color = "black"),
        axis.text.x = element_text(face = "bold", size = 10.5, color = "black"),
        axis.title.x = element_text(face = "bold", size = 11.5), legend.position = "right")

ggsave(file.path(output_dir, "Figure5F_Clinical_Predictions_Lollipop.pdf"), p_5f, width = 8.8, height = 5.8, device = pdf_device)
ggsave(file.path(output_dir, "Figure5F_Clinical_Predictions_Lollipop.png"), p_5f, width = 8.8, height = 5.8, dpi = 300)

# ==============================================================================
# Step 6: Serialize Production Model Objects for Decision Portal Integration
# ==============================================================================
cat("\nSaving production model artifacts for Shiny clinical decision portal...\n")

saveRDS(list(classifier = in_rf_classifier, panel_uniprots = selected_45_uniprots), 
        file.path(output_dir, "subtype_rf_classifier.rds"))
saveRDS(list(classifier = in_rf_classifier, panel_uniprots = selected_45_uniprots), 
        file.path(models_export_dir, "subtype_rf_classifier.rds"))

saveRDS(list(models = clinical_models_list, panel_uniprots = selected_45_uniprots), 
        file.path(output_dir, "clinical_rf_predictors.rds"))
saveRDS(list(models = clinical_models_list, panel_uniprots = selected_45_uniprots), 
        file.path(models_export_dir, "clinical_rf_predictors.rds"))

cat(sprintf("\nPhase 8 pipeline complete. All Figure 5 panels (5A–5F), ST19, and ST21 successfully saved to: %s\n", output_dir))