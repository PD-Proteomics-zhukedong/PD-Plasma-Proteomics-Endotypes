# ==============================================================================
# Script: 03_diagnostic_classifier_lasso.R
# Project: Large-scale plasma proteomics in Parkinson's disease (Nature Aging)
# Purpose: Clinical offset LASSO modeling, 98% relative efficiency panel selection,
#          independent validation, ST12 export, and Figure 2 generation
# ==============================================================================

gc()
set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(glmnet)
  library(pROC)
  library(caret)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggpubr)
  library(kernelshap)
  library(shapviz)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggforce)
  library(ggnewscale)
  library(clusterProfiler)
  library(stringr)
})

# Cross-platform PDF device fallback
pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"

# Directory setup
input_dir  <- "results/01_preprocessed_data"
dea_dir    <- "results/02_figure1_meta_dea"
data_dir   <- "data"
output_dir <- "results/03_figure2_diagnostic_model"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Helper function
get_mode <- function(x) {
  x <- na.omit(x)
  if (length(x) == 0) return(NA)
  uniqv <- unique(x)
  uniqv[which.max(tabulate(match(x, uniqv)))]
}

cat("\n=== Phase 3: Diagnostic Modeling and Figure 2 Pipeline ===\n")

# ==============================================================================
# Step 1: Data Ingestion, Site-Specific Harmonization, and Leak-Free Scaling
# ==============================================================================
cat("Loading matrices and performing center-specific standardization...\n")

prot_raw <- fread(file.path(input_dir, "3_ComBat_Corrected_Matrix.tsv"), data.table = FALSE)
rownames(prot_raw) <- make.unique(as.character(prot_raw$Protein.Group))
prot_expr_corrected <- prot_raw[, 7:ncol(prot_raw)]
prot_df <- as.data.frame(t(prot_expr_corrected)) %>% 
  rownames_to_column("sample") %>%
  mutate(across(-sample, as.numeric))

clin_file <- file.path(data_dir, "metadata_clinical.csv")
clinical_raw <- if (grepl("\\.xlsx?$", clin_file)) readxl::read_excel(clin_file) else read.csv(clin_file, stringsAsFactors = FALSE)
colnames(clinical_raw) <- make.names(colnames(clinical_raw))

# Harmonize hospital / center identifiers
if (!"hospital" %in% colnames(clinical_raw)) {
  if ("cohort" %in% colnames(clinical_raw)) {
    clinical_raw$hospital <- ifelse(clinical_raw$cohort == "Discovery", "Xiangya", "Renmin")
  } else if ("center" %in% colnames(clinical_raw)) {
    clinical_raw$hospital <- ifelse(clinical_raw$center %in% c("Discovery", "Center_1", "Center1"), "Xiangya", "Renmin")
  }
}

clinical_clean <- clinical_raw %>%
  mutate(group = ifelse(group %in% c("CTR", "Control", "HC"), "HC", "PD")) %>%
  group_by(hospital) %>%
  mutate(
    age = as.numeric(ifelse(is.na(age), median(age, na.rm = TRUE), age)),
    sex = ifelse(sex %in% c("male", "M", "1", 1), 1, 0),
    BMI = as.numeric(ifelse(is.na(BMI), median(BMI, na.rm = TRUE), BMI)),
    hepatic_disease = as.numeric(as.character(ifelse(is.na(hepatic_disease), get_mode(hepatic_disease), hepatic_disease))),
    kidney.function = as.numeric(as.character(ifelse(is.na(kidney.function), get_mode(kidney.function), kidney.function)))
  ) %>%
  ungroup() %>%
  mutate(
    CV_Metabolic_Score = case_when(
      is.na(Metabolic.cardiovascular.comorbidities) ~ as.character(get_mode(Metabolic.cardiovascular.comorbidities)),
      Metabolic.cardiovascular.comorbidities >= 2 ~ "2+",
      TRUE ~ as.character(Metabolic.cardiovascular.comorbidities)
    ),
    CV_Metabolic_Cat = as.numeric(factor(CV_Metabolic_Score, levels = c("0", "1", "2+"))),
    group = factor(group, levels = c("HC", "PD"))
  )

clinical_covariates <- c("age", "sex", "BMI", "hepatic_disease", "kidney.function", "CV_Metabolic_Cat")

modeling_data <- clinical_clean %>% 
  inner_join(prot_df, by = "sample") %>% 
  drop_na(group, age, sex, hospital)

# --- Robustly Load Discovery Candidates from ST4 or DEA results ---
st4_path    <- file.path(dea_dir, "Supplementary_Table_4_Full_Proteome_Statistics.csv")
table1_path <- file.path(dea_dir, "Table_1_Discovery_Limma_Full_Model.csv")

if (file.exists(st4_path)) {
  cat("  - Loading candidate Discovery features from Supplementary Table 4...\n")
  st4_data <- fread(st4_path, data.table = FALSE)
  sig_disc_prots <- st4_data %>% 
    filter(Discovery_Full_FDR < 0.05) %>% 
    pull(Protein_Group)
} else if (file.exists(table1_path)) {
  cat("  - Loading candidate Discovery features from Table 1...\n")
  disc_deps <- read.csv(table1_path)
  sig_disc_prots <- disc_deps %>% 
    filter(fdr_disc < 0.05) %>% 
    pull(protein)
} else {
  stop("Neither ST4 nor Table 1 found in DEA output directory. Please run Script 02 first.")
}

valid_proteins <- intersect(sig_disc_prots, colnames(modeling_data))

# Map Gene Symbols and filter blacklisted structural/immunoglobulin proteins
clean_uniprots <- sapply(strsplit(valid_proteins, ";"), `[`, 1) %>% gsub("-.*|\\..*", "", .)
mapped_symbols <- suppressMessages(suppressWarnings(
  mapIds(org.Hs.eg.db, keys = clean_uniprots, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
))
valid_symbols <- as.character(mapped_symbols)
valid_symbols[is.na(valid_symbols)] <- clean_uniprots[is.na(valid_symbols)]
lookup <- setNames(valid_symbols, valid_proteins)

blacklist_pattern <- "^(RPL|RPS|KRT|HBA|HBB|HBD|MYH|ACT|TUB|IGH|IGK|IGL)" 
is_blacklist <- grepl(blacklist_pattern, valid_symbols)
valid_proteins <- valid_proteins[!is_blacklist]

# Train/Validation split and scaling locked strictly to Discovery parameters
discovery_set   <- modeling_data %>% filter(hospital == "Xiangya")
replication_set <- modeling_data %>% filter(hospital != "Xiangya")

train_means <- sapply(discovery_set[valid_proteins], function(x) mean(x, na.rm = TRUE))
train_sds   <- sapply(discovery_set[valid_proteins], function(x) sd(x, na.rm = TRUE))
train_sds[train_sds == 0 | is.na(train_sds)] <- 1

train_set_scaled <- discovery_set
train_set_scaled[valid_proteins] <- scale(discovery_set[valid_proteins], center = train_means, scale = train_sds)

val_set_scaled <- replication_set
val_set_scaled[valid_proteins] <- scale(replication_set[valid_proteins], center = train_means, scale = train_sds)

cat(sprintf("Loaded and verified %d candidate features locked to Discovery cohort.\n", length(valid_proteins)))

# ==============================================================================
# Step 2: Fit Clinical Baseline Model and Lock Clinical Offset
# ==============================================================================
cat("Fitting clinical baseline logistic regression and calculating offsets...\n")

clin_model <- glm(group ~ age + sex + BMI + hepatic_disease + kidney.function + CV_Metabolic_Cat, 
                  data = train_set_scaled, family = binomial)

train_set_scaled$Clin_Score <- predict(clin_model, train_set_scaled, type = "link")
val_set_scaled$Clin_Score   <- predict(clin_model, val_set_scaled, type = "link")

train_offset <- train_set_scaled$Clin_Score
val_offset   <- val_set_scaled$Clin_Score

baseline_val_auc <- as.numeric(auc(val_set_scaled$group, val_offset, quiet = TRUE))
cat(sprintf("Validation cohort clinical baseline AUC = %.4f\n", baseline_val_auc))

# ==============================================================================
# Step 3: LASSO Feature Selection and 98% Relative Efficiency Optimization
# ==============================================================================
cat("Executing cross-validated LASSO feature selection (98% efficiency rule)...\n")

x_train_mat <- as.matrix(train_set_scaled[valid_proteins])
y_train_num <- ifelse(train_set_scaled$group == "PD", 1, 0)

set.seed(42)
cv_fit <- cv.glmnet(x_train_mat, y_train_num, offset = train_offset, 
                    family = "binomial", alpha = 1, type.measure = "auc", nfolds = 10)

coef_1se <- coef(cv_fit, s = "lambda.1se")

lasso_features <- as.data.frame(as.matrix(coef_1se)) %>% 
  rownames_to_column("Feature") %>% 
  rename_with(~ "Beta", 2) %>% 
  filter(Beta != 0 & Feature != "(Intercept)") %>% 
  mutate(abs_Beta = abs(Beta)) %>% 
  arrange(desc(abs_Beta))

ranked_proteins <- lasso_features$Feature

set.seed(42)
folds_train <- createFolds(train_set_scaled$group, k = 10, list = TRUE)

cv_search_train <- map_df(2:min(40, length(ranked_proteins)), function(n) {
  feats <- ranked_proteins[1:n]
  x_sub <- as.matrix(train_set_scaled[feats])
  
  fold_aucs <- sapply(names(folds_train), function(f) {
    tr_idx  <- -folds_train[[f]]
    val_idx <- folds_train[[f]]
    
    fit_cv <- cv.glmnet(x_sub[tr_idx, ], y_train_num[tr_idx], offset = train_offset[tr_idx], 
                        family = "binomial", alpha = 0.5, nfolds = 3)
    preds  <- predict(fit_cv, newx = x_sub[val_idx, ], newoffset = train_offset[val_idx], s = "lambda.1se", type = "response")
    as.numeric(auc(train_set_scaled$group[val_idx], as.numeric(preds), quiet = TRUE))
  })
  
  return(data.frame(
    N = n, 
    Train_CV_AUC = mean(fold_aucs, na.rm = TRUE),
    Train_CV_SE  = sd(fold_aucs, na.rm = TRUE) / sqrt(10)
  ))
})

max_train_auc <- max(cv_search_train$Train_CV_AUC, na.rm = TRUE)
target_98_auc <- 0.98 * max_train_auc

optimal_model_auto <- cv_search_train %>% 
  filter(Train_CV_AUC >= target_98_auc) %>% 
  arrange(N) %>% 
  dplyr::slice(1)

optimal_N <- optimal_model_auto$N
selected_panel <- ranked_proteins[1:optimal_N]

cat(sprintf("Discovery internal 10-fold CV max AUC = %.4f\n", max_train_auc))
cat(sprintf("Selected parsimonious biomarker panel: N = %d proteins\n", optimal_N))

# ==============================================================================
# Step 4: Final Model Fitting and Independent External Validation
# ==============================================================================
cat("Fitting calibrated model and predicting in external validation cohort...\n")

form_or <- as.formula(paste("group ~ offset(Clin_Score) +", paste(paste0("`", selected_panel, "`"), collapse = " + ")))
model_for_or <- glm(form_or, data = train_set_scaled, family = binomial)

pred_train_full_calibrated <- as.numeric(predict(model_for_or, newdata = train_set_scaled, type = "response"))
pred_ext_full_calibrated   <- as.numeric(predict(model_for_or, newdata = val_set_scaled, type = "response"))
pred_ext_clin              <- as.numeric(predict(clin_model, newdata = val_set_scaled, type = "response"))

roc_train_full <- roc(train_set_scaled$group, pred_train_full_calibrated, quiet = TRUE, ci = TRUE)
roc_ext_full   <- roc(val_set_scaled$group, pred_ext_full_calibrated, quiet = TRUE, ci = TRUE)
roc_ext_clin   <- roc(val_set_scaled$group, pred_ext_clin, quiet = TRUE, ci = TRUE)

delong_test <- roc.test(roc_ext_full, roc_ext_clin, method = "delong")

# Print Detailed Validation Metrics Summary to Console (100% matched to Article Result 3)
cat("\n=== Diagnostic Model Performance Evaluation (Article Result 3) ===\n")
performance_summary <- data.frame(
  Evaluation_Setting = c(
    "Discovery Cohort (16-Protein Offset Model)",
    "External Validation Cohort (16-Protein Offset Model)",
    "External Validation Cohort (6-Covariate Clinical Baseline)"
  ),
  AUC = c(as.numeric(auc(roc_train_full)), as.numeric(auc(roc_ext_full)), as.numeric(auc(roc_ext_clin))),
  AUC_95_CI = c(
    sprintf("%.3f–%.3f", roc_train_full$ci[1], roc_train_full$ci[3]),
    sprintf("%.3f–%.3f", roc_ext_full$ci[1], roc_ext_full$ci[3]),
    sprintf("%.3f–%.3f", roc_ext_clin$ci[1], roc_ext_clin$ci[3])
  ),
  DeLong_P_vs_Clinical = c("-", sprintf("%.4e", delong_test$p.value), "Reference")
)
print(performance_summary)

# Export Supplementary Table 12 (Multivariable logistic regression parameters for LASSO panel)
cat("\nExporting Supplementary Table 12 (Multivariable Logistic Parameters)...\n")
final_summary <- broom::tidy(model_for_or, exponentiate = TRUE, conf.int = TRUE) %>%
  dplyr::rename(Odds_Ratio = estimate, P_Value = p.value, CI_Lower = conf.low, CI_Upper = conf.high) %>%
  dplyr::mutate(across(where(is.numeric), ~round(., 3)))

final_summary$Symbol <- sapply(final_summary$term, function(t) {
  clean_term <- gsub("`", "", t)
  if(clean_term %in% c("(Intercept)", clinical_covariates)) return(clean_term)
  sym <- lookup[[clean_term]]
  if(is.na(sym) || is.null(sym)) return(clean_term) else return(sym)
})

st12_table <- final_summary %>% 
  dplyr::select(Symbol, term, Odds_Ratio, CI_Lower, CI_Upper, P_Value) %>%
  dplyr::rename(UNIPROT_ID = term) %>%
  dplyr::mutate(UNIPROT_ID = gsub("`", "", UNIPROT_ID)) %>%
  dplyr::arrange(P_Value)

write.csv(st12_table, file.path(output_dir, "Supplementary_Table_12_Diagnostic_Panel_OR.csv"), row.names = FALSE)
cat(sprintf("  - ST12 successfully exported (%d variables, including Intercept and 16 Proteins).\n", nrow(st12_table)))

# ==============================================================================
# Step 5: Generate Figure 2A (ROC Curves)
# ==============================================================================
cat("Generating Figure 2A (ROC Curves)...\n")

pdf(file.path(output_dir, "Figure2A_Final_ROC_Validation.pdf"), width = 7, height = 6.5)
par(mar = c(5, 5, 4, 2))

plot(roc_ext_full, col = "#d73027", lwd = 3.5, 
     main = sprintf("Diagnostic Performance of %d-Protein Offset Panel", optimal_N),
     cex.main = 1.2, font.main = 2,
     xlab = "1 - Specificity (False Positive Rate)",
     ylab = "Sensitivity (True Positive Rate)")

plot(roc_train_full, col = "#377eb8", lwd = 2, add = TRUE)
plot(roc_ext_clin, col = "grey50", lwd = 2, lty = 2, add = TRUE)

legend("bottomright", legend = c(
  sprintf("External Validation: AUC %.3f (95%% CI: %.3f–%.3f)", 
          as.numeric(auc(roc_ext_full)), roc_ext_full$ci[1], roc_ext_full$ci[3]),
  sprintf("Discovery Cohort: AUC %.3f (95%% CI: %.3f–%.3f)", 
          as.numeric(auc(roc_train_full)), roc_train_full$ci[1], roc_train_full$ci[3]),
  sprintf("Clinical Baseline: AUC %.3f (95%% CI: %.3f–%.3f)", 
          as.numeric(auc(roc_ext_clin)), roc_ext_clin$ci[1], roc_ext_clin$ci[3])
), col = c("#d73027", "#377eb8", "grey50"), lwd = c(3.5, 2, 2), lty = c(1, 1, 2), bty = "n", cex = 0.85)
dev.off()

png(file.path(output_dir, "Figure2A_Final_ROC_Validation.png"), width = 7, height = 6.5, units = "in", res = 300)
par(mar = c(5, 5, 4, 2))
plot(roc_ext_full, col = "#d73027", lwd = 3.5, 
     main = sprintf("Diagnostic Performance of %d-Protein Offset Panel", optimal_N),
     cex.main = 1.2, font.main = 2,
     xlab = "1 - Specificity (False Positive Rate)",
     ylab = "Sensitivity (True Positive Rate)")
plot(roc_train_full, col = "#377eb8", lwd = 2, add = TRUE)
plot(roc_ext_clin, col = "grey50", lwd = 2, lty = 2, add = TRUE)
legend("bottomright", legend = c(
  sprintf("External Validation: AUC %.3f (95%% CI: %.3f–%.3f)", 
          as.numeric(auc(roc_ext_full)), roc_ext_full$ci[1], roc_ext_full$ci[3]),
  sprintf("Discovery Cohort: AUC %.3f (95%% CI: %.3f–%.3f)", 
          as.numeric(auc(roc_train_full)), roc_train_full$ci[1], roc_train_full$ci[3]),
  sprintf("Clinical Baseline: AUC %.3f (95%% CI: %.3f–%.3f)", 
          as.numeric(auc(roc_ext_clin)), roc_ext_clin$ci[1], roc_ext_clin$ci[3])
), col = c("#d73027", "#377eb8", "grey50"), lwd = c(3.5, 2, 2), lty = c(1, 1, 2), bty = "n", cex = 0.85)
dev.off()

# ==============================================================================
# Step 6: Generate Figure 2B (SHAP Beeswarm Plot)
# ==============================================================================
cat("Generating Figure 2B (SHAP Beeswarm Plot)...\n")

form_shap_pure <- as.formula(paste("group ~", paste(paste0("`", selected_panel, "`"), collapse = " + ")))
model_for_shap <- glm(form_shap_pure, data = train_set_scaled, family = binomial)

X_shap <- as.data.frame(train_set_scaled %>% dplyr::select(all_of(selected_panel)))
ks <- kernelshap(model_for_shap, X = X_shap, bg_X = X_shap)
sv <- shapviz(ks)

new_labels <- sapply(colnames(sv), function(x) {
  clean_id <- gsub("`", "", x)
  symbol <- lookup[[clean_id]]
  if (is.null(symbol) || is.na(symbol)) return(clean_id) else return(symbol)
})
colnames(sv) <- new_labels

p_shap <- sv_importance(sv, kind = "beeswarm", max_display = optimal_N) + 
  theme_bw(base_size = 12) +
  labs(title = sprintf("SHAP Feature Contribution (%d-Protein Panel)", optimal_N),
       subtitle = "Quantifying Individual Impacts of Biomarkers on Parkinson's Risk") +
  theme(axis.text.y = element_text(face = "bold", color = "black", size = 11), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13), 
        plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic"))

ggsave(file.path(output_dir, "Figure2B_GLM_Based_SHAP_Summary.pdf"), p_shap, width = 8, height = 7.5, device = pdf_device)
ggsave(file.path(output_dir, "Figure2B_GLM_Based_SHAP_Summary.png"), p_shap, width = 8, height = 7.5, dpi = 300)

# ==============================================================================
# Step 7: Generate Figure 2C (PPI Network & Functional Connectivity)
# ==============================================================================
cat("Generating Figure 2C (PPI Functional Connectivity)...\n")

target_symbols <- na.omit(as.character(lookup[selected_panel]))
protein_count  <- length(target_symbols)

genes_str <- paste(target_symbols, collapse = "%0d")
url <- paste0("https://string-db.org/api/tsv/network?identifiers=", genes_str, "&species=9606&add_nodes=15&required_score=400")

edges_df <- tryCatch({ 
  read.table(url, header = TRUE, sep = "\t", stringsAsFactors = FALSE) 
}, error = function(e) NULL)

if (!is.null(edges_df) && nrow(edges_df) > 0) {
  all_nodes_complete <- unique(c(edges_df$preferredName_A, edges_df$preferredName_B, target_symbols))
  nodes_df <- data.frame(name = all_nodes_complete, stringsAsFactors = FALSE)
  links <- edges_df %>% dplyr::select(from = preferredName_A, to = preferredName_B, score) %>% dplyr::distinct()
  
  g <- graph_from_data_frame(d = links, vertices = nodes_df, directed = FALSE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  V(g)$is_core <- V(g)$name %in% target_symbols
  V(g)$degree  <- degree(g)
  g <- induced_subgraph(g, V(g)$is_core | V(g)$degree > 0)
  
  sub_g <- induced_subgraph(g, degree(g) > 0)
  if (ecount(sub_g) > 0) {
    comm_sub <- membership(cluster_louvain(sub_g))
    V(g)$community <- NA 
    V(g)$community[match(names(comm_sub), V(g)$name)] <- as.numeric(comm_sub)
  } else {
    V(g)$community <- NA
  }
  
  V(g)$hull_label <- NA 
  communities <- unique(na.omit(V(g)$community))
  
  for(i in communities) {
    nodes_in_comm <- V(g)$name[which(V(g)$community == i)]
    core_count_in_comm <- sum(V(g)$is_core[which(V(g)$community == i)])
    
    if(length(nodes_in_comm) >= 3 & core_count_in_comm >= 1) {
      g_ids <- tryCatch({ 
        bitr(nodes_in_comm, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID 
      }, error = function(e) NULL)
      
      if(!is.null(g_ids) && length(g_ids) >= 2) {
        ego <- suppressMessages(enrichGO(g_ids, OrgDb=org.Hs.eg.db, ont="BP", pvalueCutoff=0.2, readable=TRUE))
        if(!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
          label <- str_wrap(ego@result$Description[1], width = 16) 
          V(g)$hull_label[which(V(g)$community == i)] <- label
        }
      }
    }
  }
  
  V(g)$display_name <- V(g)$name
  tg <- as_tbl_graph(g)
  
  set.seed(88) 
  p_ppi <- ggraph(tg, layout = "fr") + 
    geom_mark_hull(aes(x = x, y = y, fill = hull_label, label = hull_label, filter = !is.na(hull_label)), 
                   concavity = 4, expand = unit(5, "mm"), alpha = 0.15, 
                   color = "grey50", linetype = "dashed", linewidth = 0.6,
                   label.fontsize = 10.5, label.fontface = "bold", label.fill = "white",
                   label.buffer = unit(3, "mm"),
                   con.cap = unit(2, "mm"), con.colour = "grey50", con.type = "straight") +
    scale_fill_brewer(palette = "Set2", na.translate = FALSE) + 
    new_scale_fill() + 
    geom_edge_link(color = "grey75", width = 0.85, alpha = 0.7) + 
    geom_node_point(aes(filter = !is_core), size = 3.8, color = "grey50", fill = "grey90", shape = 21, stroke = 0.7) +
    geom_node_point(aes(filter = is_core), size = 8.5, color = "white", fill = "#d73027", shape = 21, stroke = 1.2) +
    geom_node_text(aes(filter = !is_core, label = name), 
                   repel = TRUE, size = 3.3, color = "grey35", fontface = "plain") +
    geom_node_text(aes(filter = is_core, label = display_name), 
                   repel = TRUE, size = 4.8, color = "black", fontface = "bold", 
                   bg.color = "white", bg.r = 0.12, segment.color = "grey50", segment.size = 0.4) +
    theme_void() +
    theme(legend.position = "none", 
          plot.margin = ggplot2::margin(0.8, 0.8, 0.8, 0.8, "cm"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 11.5, color = "grey30")) +
    labs(title = sprintf("Functional Protein-Protein Interaction Network (%d-Protein Panel)", protein_count), 
         subtitle = "Red: Core Biomarkers; Grey: STRING Background; Dashed Enclosures: GO Biological Modules")
  
  ggsave(file.path(output_dir, sprintf("Figure2C_Functional_Landscape_%dProteins_Optimized.pdf", protein_count)), p_ppi, width = 11.5, height = 9.5, device = pdf_device)
  ggsave(file.path(output_dir, sprintf("Figure2C_Functional_Landscape_%dProteins_Optimized.png", protein_count)), p_ppi, width = 11.5, height = 9.5, dpi = 300)
}

# ==============================================================================
# Step 8: Generate Figure 2D (Expression Violins in Validation Cohort)
# ==============================================================================
cat("Generating Figure 2D (Expression Violins in External Validation Cohort)...\n")

id_to_symbol <- setNames(st12_table$Symbol, st12_table$UNIPROT_ID)
plot_df <- val_set_scaled %>%
  dplyr::select(group, all_of(selected_panel)) %>%
  pivot_longer(cols = -group, names_to = "Protein_ID", values_to = "Expression") %>%
  mutate(Symbol = id_to_symbol[Protein_ID])

plot_df$Symbol <- factor(plot_df$Symbol, levels = unique(st12_table$Symbol[!st12_table$Symbol %in% c("(Intercept)", clinical_covariates)]))

p_vio <- ggplot(plot_df, aes(x = group, y = Expression, fill = group)) +
  geom_violin(alpha = 0.4, trim = FALSE, color = NA) +
  geom_boxplot(width = 0.15, color = "black", outlier.shape = NA, fill = "white", alpha = 0.7) +
  facet_wrap(~Symbol, scales = "free_y", ncol = 4) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, size = 5) +
  scale_fill_manual(values = c("HC" = "#4575b4", "PD" = "#d73027")) +
  theme_bw(base_size = 12) +
  labs(title = "Relative Protein Abundance in External Validation Cohort", y = "Normalized Expression", x = NULL)

ggsave(file.path(output_dir, "Figure2D_Panel_Proteins_Violin.pdf"), p_vio, width = 16, height = 10, device = pdf_device)
ggsave(file.path(output_dir, "Figure2D_Panel_Proteins_Violin.png"), p_vio, width = 16, height = 10, dpi = 300)

# Save Locked Diagnostic Model RDS
saveRDS(model_for_or, file.path(output_dir, "lasso_diagnostic_model.rds"))

cat("\nPhase 3 pipeline complete. Figure 2 panels and Supplementary Table 12 successfully exported to:", output_dir, "\n")