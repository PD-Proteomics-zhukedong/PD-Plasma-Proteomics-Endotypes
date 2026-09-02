# ==============================================================================
# Script: 09_subtype_classifier_and_clinical_prediction.R
# Purpose: Network topology, subtype classification, and 10-fold CV clinical trait predictions
# Project: Large-scale plasma proteomics identifies molecularly distinct PD endotypes
# ==============================================================================

options(expressions = 5000)
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(pROC)
  library(randomForest)
  library(caret)
  library(ggpubr)
  library(org.Hs.eg.db)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggforce)
  library(ggnewscale)
  library(patchwork)
  library(scales)
  library(stringr)
})

# Directory setup
data_dir          <- "data"
input_dea_dir     <- "results/02_figure1_meta_dea"
input_clust_dir   <- "results/07_figure3_molecular_endotypes"
input_mech_dir    <- "results/08_subtype_mechanisms_and_deconvolution"
output_dir        <- "results/09_figure5_machine_learning_translation"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

BASE_COVARIATES <- c("Age", "Sex", "BMI", "hepatic_disease", "kidney.function", "CV_Metabolic_Cat", "LEDD", "Duration")

cat("=== Phase 9: Network Topology, Subtype Classifier, and Clinical Trait Prediction ===\n")

# ==============================================================================
# Step 1: Ingest Preprocessed Matrices, Clinical Data, and 45 Fingerprint Biomarkers
# ==============================================================================
cat("Loading preprocessed matrices, clinical metadata, and 45 fingerprint biomarkers...\n")

clin_file    <- file.path(data_dir, "metadata_clinical_n1119.csv")
combat_file  <- file.path("results/01_preprocessed_data", "3_ComBat_Corrected_Matrix.tsv")
st16_file    <- file.path(input_mech_dir, "Supplementary_Table_16_Subtype_Biomarker_Specificity.csv")

if (!file.exists(clin_file) || !file.exists(combat_file) || !file.exists(st16_file)) {
  stop("Required input files not found. Please verify previous pipeline scripts have run.")
}

clinical_raw <- read.csv(clin_file, stringsAsFactors = FALSE)
colnames(clinical_raw) <- make.names(colnames(clinical_raw))

# Filter sporadic PD cohort with complete baseline records
clinical_pd <- clinical_raw %>%
  filter(group == "PD" & !is.na(Age) & !is.na(Sex)) %>%
  mutate(
    LEDD = replace_na(as.numeric(LEDD), 0),
    Duration = replace_na(as.numeric(Disease_Duration), 0),
    CV_Metabolic_Cat = as.numeric(as.factor(CV_Metabolic_Score)),
    hepatic_disease = replace_na(as.numeric(hepatic_disease), 0),
    kidney.function = replace_na(as.numeric(kidney.function), 0)
  )

prot_raw_full <- fread(combat_file, data.table = FALSE)
rownames(prot_raw_full) <- make.unique(as.character(prot_raw_full$Protein.Group))
prot_expr_mat <- as.matrix(prot_raw_full[, 7:ncol(prot_raw_full)])

common_samples <- intersect(colnames(prot_expr_mat), clinical_pd$sample)
prot_pd_mat    <- prot_expr_mat[, common_samples]
clinical_final <- clinical_pd[match(common_samples, clinical_pd$sample), ]

# Retrieve Subtype labels
if (!"Subtype" %in% colnames(clinical_final)) {
  st7_file <- file.path(input_clust_dir, "Supplementary_Table_7_Molecular_Endotypes_Baseline.csv")
  # Fallback to direct clinical table if present
}

# Gene mapping lookup
clean_uniprots <- gsub("-.*|\\..*", "", rownames(prot_pd_mat))
mapped_symbols <- suppressMessages(suppressWarnings(
  mapIds(org.Hs.eg.db, keys = clean_uniprots, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
))
mapped_symbols[is.na(mapped_symbols)] <- clean_uniprots[is.na(mapped_symbols)]
lookup_global  <- setNames(as.character(mapped_symbols), rownames(prot_pd_mat))
lookup_reverse <- setNames(names(lookup_global), lookup_global)

# Ingest the 45 subtype fingerprint markers from Table ST16 (Top 15 per subtype)
st16_df <- read.csv(st16_file, stringsAsFactors = FALSE)

selected_45_df <- st16_df %>%
  group_by(Specific_Endotype) %>%
  arrange(desc(Min_T_Score)) %>%
  slice_head(n = 15) %>%
  ungroup()

selected_45_uniprots <- intersect(selected_45_df$UNIPROT, rownames(prot_pd_mat))
cat(sprintf("Locked %d fingerprint biomarkers across Subtypes 1, 2, and 3.\n", length(selected_45_uniprots)))

# Partition independent cohorts: Discovery (n = 471) vs Validation (n = 198)
train_idx <- which(clinical_final$cohort == "Discovery")
val_idx   <- which(clinical_final$cohort == "Validation")

cat(sprintf("Cohort partitions: Discovery (n = %d), Validation (n = %d).\n", length(train_idx), length(val_idx)))

# ==============================================================================
# Step 2: Dual-Pillar Subtype PPI Networks (Figure 5A, 5B, 5C)
# ==============================================================================
cat("Constructing dual-pillar PPI networks for Subtypes 1, 2, and 3 (Figure 5A-C)...\n")

fetch_string_edges <- function(symbols) {
  genes_param <- paste(symbols, collapse = "%0d")
  string_url <- paste0("https://string-db.org/api/tsv/network?identifiers=", 
                       genes_param, "&species=9606&required_score=150")
  edges_df <- tryCatch({
    res <- read.delim(string_url, stringsAsFactors = FALSE)
    if (nrow(res) > 0 && "preferredName_A" %in% colnames(res)) {
      res %>% 
        dplyr::select(from = preferredName_A, to = preferredName_B, score = score) %>%
        filter(from %in% symbols & to %in% symbols & from != to) %>%
        distinct()
    } else NULL
  }, error = function(e) NULL)
  if (is.null(edges_df)) edges_df <- data.frame(from=character(), to=character(), score=numeric())
  return(edges_df)
}

draw_subtype_ppi_panel <- function(subtype_num, p1_name, p2_name, p1_genes, p2_genes) {
  sub_label <- sprintf("Subtype %d", subtype_num)
  sub_markers <- selected_45_df %>% filter(Specific_Endotype == sub_label)
  symbols <- sub_markers$Symbol
  
  mod_map <- setNames(rep(p2_name, length(symbols)), symbols)
  mod_map[intersect(symbols, p1_genes)] <- p1_name
  
  edges_df <- fetch_string_edges(symbols)
  
  g <- graph_from_data_frame(d = edges_df, vertices = data.frame(name = symbols), directed = FALSE)
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  V(g)$Module       <- mod_map[V(g)$name]
  V(g)$degree       <- as.numeric(degree(g)) + 7
  V(g)$logFC        <- pmax(sub_markers$Min_T_Score[match(V(g)$name, sub_markers$Symbol)], 0.2)
  V(g)$display_name <- ifelse(V(g)$name == "A0A0B4J2B5", "IGHV3-7", V(g)$name)
  
  set.seed(42 + subtype_num * 10)
  g1 <- induced_subgraph(g, intersect(p1_genes, V(g)$name))
  lay1 <- layout_with_fr(g1, niter = 300)
  lay1 <- scale(lay1)
  lay1[, 1] <- lay1[, 1] * 0.7 - 1.8
  rownames(lay1) <- V(g1)$name
  
  g2 <- induced_subgraph(g, intersect(p2_genes, V(g)$name))
  lay2 <- layout_with_fr(g2, niter = 300)
  lay2 <- scale(lay2)
  lay2[, 1] <- lay2[, 1] * 0.7 + 1.8
  rownames(lay2) <- V(g2)$name
  
  combined_lay <- rbind(lay1, lay2)[symbols, ]
  layout_net   <- create_layout(g, layout = "manual", x = combined_lay[,1], y = combined_lay[,2])
  
  banner_df <- layout_net %>%
    group_by(Module) %>%
    summarise(x_center = mean(x), y_top = max(y) + 0.35, .groups = "drop")
  
  mod_colors <- c("#4DAF4A", "#984EA3")
  mod_color_map <- setNames(mod_colors, c(p1_name, p2_name))
  
  p <- ggraph(layout_net) +
    geom_mark_hull(aes(x = x, y = y, fill = Module), concavity = 1.0, expand = unit(6.0, "mm"),
                   radius = unit(5.0, "mm"), alpha = 0.12, color = "grey45", linetype = "dashed", linewidth = 0.65, show.legend = FALSE) +
    scale_fill_manual(values = mod_color_map) +
    geom_label(data = banner_df, aes(x = x_center, y = y_top, label = Module), fontface = "bold", size = 3.5,
               fill = "white", color = "black", label.padding = unit(2.5, "pt"), label.r = unit(3.0, "pt")) +
    new_scale_fill() +
    geom_edge_link(color = "grey65", alpha = 0.50, linewidth = 0.70) +
    geom_node_point(aes(size = degree, fill = logFC), shape = 21, color = "black", stroke = 0.85) +
    geom_node_text(aes(label = display_name), repel = TRUE, fontface = "bold", size = 3.8, color = "black",
                   bg.color = "white", bg.r = 0.15) +
    scale_fill_gradientn(colors = c("#FFF7BC", "#FEC44F", "#D95F0E", "#990000"),
                         name = sprintf("Subtype %d Specificity\n(Log2FC vs Rest)", subtype_num)) +
    scale_size_continuous(range = c(6.0, 11.5), name = "PPI Degree (Hub)") +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.25))) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = 11) +
    theme(legend.position = "right",
          plot.margin = ggplot2::margin(20, 15, 15, 15))
  return(p)
}

# Subtype 1 PPI
p_5a <- draw_subtype_ppi_panel(
  subtype_num = 1,
  p1_name     = "Integrin & Cell Adhesion",
  p2_name     = "Signaling & Proteostasis",
  p1_genes    = c("ITGA6", "ITGA5", "CD151", "ICAM2", "ADAM10", "TSPAN14", "PLXNB2", "SLC44A2"),
  p2_genes    = c("RHOF", "TAOK3", "EIF5B", "PDIA5", "PURA", "GNG5", "CSNK2A2")
)

# Subtype 2 PPI
p_5b <- draw_subtype_ppi_panel(
  subtype_num = 2,
  p1_name     = "Complement & Humoral Immunity",
  p2_name     = "Vascular & Lipid Homeostasis",
  p1_genes    = c("VTN", "C4A", "C4BPA", "SERPING1", "PIGR", "A0A0B4J2B5", "ILF2", "IGHV3-7"),
  p2_genes    = c("PROS1", "MGP", "APOC1", "AMBP", "KNG1", "APOA4", "PTK7", "TMEM9")
)

# Subtype 3 PPI
p_5c <- draw_subtype_ppi_panel(
  subtype_num = 3,
  p1_name     = "Mitochondrial Respiration & ATP",
  p2_name     = "Metabolic & RNA Homeostasis",
  p1_genes    = c("NDUFA9", "ATP5MK", "STOML2", "LONP1", "SLC25A20", "ECI1", "MCCC2"),
  p2_genes    = c("EDC4", "PPP1R18", "CNOT1", "CLCN4", "APAF1", "IRAG2", "MYCBP2", "GLS")
)

# ==============================================================================
# Step 3: Multi-Class Balanced Random Forest & Gini Importance (Figure 5D)
# ==============================================================================
cat("Training Balanced Random Forest on Discovery cohort and calculating Gini importance (Figure 5D)...\n")

X_train <- t(prot_pd_mat[selected_45_uniprots, train_idx])
Y_train <- factor(clinical_final$Subtype[train_idx])

X_val   <- t(prot_pd_mat[selected_45_uniprots, val_idx])
Y_val   <- factor(clinical_final$Subtype[val_idx])

df_train <- as.data.frame(X_train) %>% mutate(Subtype = Y_train)
df_val   <- as.data.frame(X_val) %>% mutate(Subtype = Y_val)

min_size <- min(table(Y_train))

set.seed(42)
rf_classifier <- randomForest(
  Subtype ~ ., 
  data     = df_train, 
  ntree    = 1000, 
  strata   = df_train$Subtype, 
  sampsize = rep(min_size, length(levels(Y_train))),
  importance = TRUE
)

# Extract variable importance (Mean Decrease Gini)
imp_df <- as.data.frame(importance(rf_classifier)) %>%
  rownames_to_column("UNIPROT") %>%
  mutate(Symbol = ifelse(UNIPROT == "A0A0B4J2B5", "IGHV3-7", lookup_global[UNIPROT])) %>%
  arrange(desc(MeanDecreaseGini))

top20_imp <- head(imp_df, 20)

p_5d <- ggplot(top20_imp, aes(x = MeanDecreaseGini, y = reorder(Symbol, MeanDecreaseGini))) +
  geom_segment(aes(x = 0, xend = MeanDecreaseGini, y = Symbol, yend = Symbol), color = "#377eb8", linewidth = 1.1) +
  geom_point(size = 3.8, color = "#d73027", fill = "#f4a582", shape = 21, stroke = 1.1) +
  labs(title = "Biological Feature Importance",
       x = "Variable Importance (Mean Decrease Gini)", y = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_bw(base_size = 11.5) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.text.y = element_text(face = "bold", size = 9.5, color = "black"),
        panel.grid.minor = element_blank())

# ==============================================================================
# Step 4: Multi-Class ROC Evaluation on Validation Cohort (Figure 5E)
# ==============================================================================
cat("Evaluating multi-class ROC curves on independent Validation cohort (Figure 5E)...\n")

val_probs <- predict(rf_classifier, newdata = df_val, type = "prob")

roc_s1 <- roc(df_val$Subtype == "Subtype1", val_probs[, "Subtype1"], quiet = TRUE, ci = TRUE)
roc_s2 <- roc(df_val$Subtype == "Subtype2", val_probs[, "Subtype2"], quiet = TRUE, ci = TRUE)
roc_s3 <- roc(df_val$Subtype == "Subtype3", val_probs[, "Subtype3"], quiet = TRUE, ci = TRUE)

pdf(file.path(output_dir, "Figure5E_Subtype_Classifier_ROC.pdf"), width = 6.5, height = 6.2)
par(mar = c(4.5, 4.5, 2.5, 2.0))

plot(roc_s1, col = "#d73027", lwd = 3.2, main = "Subtype Classifier (Validation Cohort)",
     cex.main = 1.15, font.main = 2, xlab = "1 - Specificity (False Positive Rate)", ylab = "Sensitivity (True Positive Rate)")
plot(roc_s2, col = "#4575b4", lwd = 3.2, add = TRUE)
plot(roc_s3, col = "#fdae61", lwd = 3.2, add = TRUE)
abline(a = 1, b = -1, col = "grey65", lwd = 1.0)

legend("bottomright", legend = c(
  sprintf("Subtype 1: AUC %.3f (95%% CI: %.3f–%.3f)", roc_s1$auc, roc_s1$ci[1], roc_s1$ci[3]),
  sprintf("Subtype 2: AUC %.3f (95%% CI: %.3f–%.3f)", roc_s2$auc, roc_s2$ci[1], roc_s2$ci[3]),
  sprintf("Subtype 3: AUC %.3f (95%% CI: %.3f–%.3f)", roc_s3$auc, roc_s3$ci[1], roc_s3$ci[3])
), col = c("#d73027", "#4575b4", "#fdae61"), lwd = 3.2, bty = "n", cex = 0.82)
dev.off()

# ==============================================================================
# Step 5: 10-Fold CV Continuous Clinical Trait Predictions (Figure 5F & Table ST18)
# ==============================================================================
cat("Running 10-fold CV regression predictions for 7 continuous clinical scales (Figure 5F & Table ST18)...\n")

clinical_traits_target <- c(
  "UPDRS_III"  = "UPDRS-III (Motor)",
  "HY_Stage"   = "Hoehn & Yahr (Stage)",
  "HAMD"       = "HAMD (Depression)",
  "NMSS"       = "NMSS (Total Non-motor)",
  "MMSE"       = "MMSE (Cognition)",
  "PDSS"       = "PDSS (Sleep)",
  "UPSIT"      = "UPSIT (Olfaction)"
)

cv_prediction_list <- list()

for (tr in names(clinical_traits_target)) {
  if (tr %in% colnames(clinical_final)) {
    tr_vec <- as.numeric(clinical_final[[tr]])
    
    df_reg <- data.frame(
      Trait = tr_vec,
      clinical_final %>% dplyr::select(all_of(BASE_COVARIATES)),
      t(prot_pd_mat[selected_45_uniprots, ])
    ) %>% drop_na(Trait)
    
    if (nrow(df_reg) >= 50) {
      colnames(df_reg) <- make.names(colnames(df_reg))
      set.seed(42)
      folds <- caret::createFolds(df_reg$Trait, k = 10, list = TRUE)
      oof_predictions <- numeric(nrow(df_reg))
      
      for (f in names(folds)) {
        v_idx <- folds[[f]]
        t_idx <- setdiff(1:nrow(df_reg), v_idx)
        rf_reg <- randomForest(Trait ~ ., data = df_reg[t_idx, ], ntree = 300)
        oof_predictions[v_idx] <- predict(rf_reg, newdata = df_reg[v_idx, ])
      }
      
      cor_res <- cor.test(df_reg$Trait, oof_predictions, method = "spearman", exact = FALSE)
      
      cv_prediction_list[[tr]] <- data.frame(
        Clinical_Trait = clinical_traits_target[tr],
        Spearman_R     = as.numeric(cor_res$estimate),
        P_Value        = as.numeric(cor_res$p.value),
        N_Samples      = nrow(df_reg)
      )
    }
  }
}

table_st18_export <- bind_rows(cv_prediction_list) %>%
  mutate(
    Log10_P  = -log10(P_Value),
    Sig_Star = case_when(
      P_Value < 0.0001 ~ "****",
      P_Value < 0.001  ~ "***",
      P_Value < 0.01   ~ "**",
      P_Value < 0.05   ~ "*",
      TRUE ~ "ns"
    ),
    Label_Text = sprintf("R = %.2f (%s)", Spearman_R, Sig_Star)
  ) %>%
  arrange(Spearman_R)

table_st18_export$Clinical_Trait <- factor(table_st18_export$Clinical_Trait, levels = table_st18_export$Clinical_Trait)

write.csv(table_st18_export, file.path(output_dir, "Supplementary_Table_18_Clinical_Scale_Predictions_10FoldCV.csv"), row.names = FALSE)

# Generate Figure 5F
p_5f <- ggplot(table_st18_export, aes(x = Spearman_R, y = Clinical_Trait)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey60", linewidth = 0.8) +
  geom_segment(aes(x = 0, xend = Spearman_R, y = Clinical_Trait, yend = Clinical_Trait), color = "#377eb8", linewidth = 1.3) +
  geom_point(aes(size = Log10_P, fill = Spearman_R), shape = 21, color = "black", stroke = 1.1) +
  scale_fill_gradientn(colors = c("#fee08b", "#f46d43", "#d73027", "#a50026"), name = "Spearman R") +
  scale_size_continuous(range = c(4.5, 9.5), name = expression(bold(-log[10](P-value))), breaks = c(10, 20, 30, 40)) +
  geom_text(aes(label = Label_Text), hjust = -0.15, size = 3.8, fontface = "bold", color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.28)), limits = c(0, 0.85)) +
  labs(title = "Cross-Validated Multi-Domain Clinical Trait Predictions",
       x = "Cross-Validated Prediction Accuracy (Spearman R)", y = NULL) +
  theme_classic(base_size = 11.5) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        axis.text.x = element_text(color = "black", size = 10),
        legend.position = "right")

# ==============================================================================
# Step 6: Multi-Panel Master Layout for Figure 5
# ==============================================================================
cat("Assembling master layout for Figure 5...\n")

fig5_top <- (p_5a | p_5b)
fig5_mid <- (p_5c | p_5d)

# Export individual panels and combined layouts
ggsave(file.path(output_dir, "Figure5A_Subtype1_PPI.pdf"), p_5a, width = 7.5, height = 7.5, device = cairo_pdf)
ggsave(file.path(output_dir, "Figure5B_Subtype2_PPI.pdf"), p_5b, width = 7.5, height = 7.5, device = cairo_pdf)
ggsave(file.path(output_dir, "Figure5C_Subtype3_PPI.pdf"), p_5c, width = 7.5, height = 7.5, device = cairo_pdf)
ggsave(file.path(output_dir, "Figure5D_Gini_Importance.pdf"), p_5d, width = 6.8, height = 6.2, device = cairo_pdf)
ggsave(file.path(output_dir, "Figure5F_Clinical_Predictions_Lollipop.pdf"), p_5f, width = 8.5, height = 5.5, device = cairo_pdf)

# Save final trained subtype classifier model object for interactive Shiny App
saveRDS(rf_classifier, file.path(output_dir, "balanced_rf_subtype_classifier.rds"))

cat("Analysis complete. Figure 5 panels, Supplementary Table 18, and serialized model saved to:", output_dir, "\n")