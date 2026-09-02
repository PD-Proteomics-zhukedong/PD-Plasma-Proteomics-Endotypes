# ==============================================================================
# Script: 03_diagnostic_classifier_lasso.R
# Purpose: LASSO regression modeling with clinical offset and Figure 2 generation
# Project: Large-scale plasma proteomics in Parkinson's disease
# ==============================================================================

options(expressions = 5000)
set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(glmnet)
  library(pROC)
  library(caret)
  library(org.Hs.eg.db)
  library(ggpubr)
  library(kernelshap)
  library(shapviz)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggforce)
  library(patchwork)
  library(scales)
  library(stringr)
  library(e1071)
  library(randomForest)
})

# Directory setup
input_dir  <- "results/01_preprocessed_data"
dea_dir    <- "results/02_figure1_meta_dea"
data_dir   <- "data"
output_dir <- "results/03_figure2_diagnostic_model"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

BASE_COVARIATES <- c("age", "sex", "BMI", "hepatic_disease", "kidney.function", "CV_Metabolic_Cat")

# Helper functions
get_mode <- function(x) {
  x <- na.omit(x)
  if(length(x) == 0) return(NA)
  uniqv <- unique(x)
  uniqv[which.max(tabulate(match(x, uniqv)))]
}

calc_nri_idi <- function(p_old, p_new, y_true) {
  y <- as.numeric(y_true) - 1
  p_new_case <- p_new[y == 1]; p_old_case <- p_old[y == 1]
  p_new_ctrl <- p_new[y == 0]; p_old_ctrl <- p_old[y == 0]
  
  nri_case <- (sum(p_new_case > p_old_case) - sum(p_new_case < p_old_case)) / length(p_new_case)
  nri_ctrl <- (sum(p_new_ctrl < p_old_ctrl) - sum(p_new_ctrl > p_old_ctrl)) / length(p_new_ctrl)
  nri <- nri_case + nri_ctrl
  se_nri <- sqrt((1 - nri_case^2)/length(p_new_case) + (1 - nri_ctrl^2)/length(p_new_ctrl))
  p_nri <- 2 * (1 - pnorm(abs(nri / se_nri)))
  
  idi_case <- mean(p_new_case) - mean(p_old_case)
  idi_ctrl <- mean(p_new_ctrl) - mean(p_old_ctrl)
  idi <- idi_case - idi_ctrl
  se_idi <- sqrt((var(p_new_case - p_old_case)/length(p_new_case)) + (var(p_new_ctrl - p_old_ctrl)/length(p_new_ctrl)))
  p_idi <- 2 * (1 - pnorm(abs(idi / se_idi)))
  
  return(data.frame(
    Metric = c("Continuous NRI", "IDI"), 
    Estimate = c(nri, idi), 
    P_Value = c(p_nri, p_idi)
  ))
}

cat("\n=== Phase 3: Diagnostic Modeling and Figure 2 Generation ===\n")

# ==============================================================================
# Step 1: Data Ingestion and Standardization
# ==============================================================================
cat("Loading matrices and standardizing features...\n")

combat_raw <- fread(file.path(input_dir, "3_ComBat_Corrected_Matrix.tsv"), data.table = FALSE)
rownames(combat_raw) <- make.unique(as.character(combat_raw$Protein.Group))
prot_df <- as.data.frame(t(combat_raw[, 7:ncol(combat_raw)])) %>% 
  rownames_to_column("sample") %>% mutate(across(-sample, as.numeric))

clinical_raw <- read.csv(file.path(data_dir, "metadata_clinical_n1119.csv"), stringsAsFactors = FALSE)
colnames(clinical_raw) <- make.names(colnames(clinical_raw))

if (!"cohort" %in% colnames(clinical_raw)) {
  if ("center" %in% colnames(clinical_raw)) {
    clinical_raw$cohort <- ifelse(clinical_raw$center %in% c("Discovery", "Center_1", "Center1"), "Discovery", "Validation")
  } else if ("hospital" %in% colnames(clinical_raw)) {
    first_site <- unique(clinical_raw$hospital)[1]
    clinical_raw$cohort <- ifelse(clinical_raw$hospital == first_site, "Discovery", "Validation")
  }
}

clinical_clean <- clinical_raw %>%
  mutate(group = factor(ifelse(group %in% c("CTR", "Control", "HC"), "HC", "PD"), levels = c("HC", "PD")),
         sex = as.factor(sex),
         CV_Metabolic_Cat = as.numeric(as.factor(Metabolic.cardiovascular.comorbidities)))

modeling_data <- clinical_clean %>% inner_join(prot_df, by = "sample") %>% drop_na(group, age, sex, cohort)

# Load Discovery DEA results
disc_deps_file <- file.path(dea_dir, "Table_1_Discovery_Limma_Full_Model.csv")
if (!file.exists(disc_deps_file)) {
  disc_deps_file <- file.path(dea_dir, "Table_ST10_Final_823_Strict_Meta_DEPs.csv")
}
disc_deps <- read.csv(disc_deps_file, stringsAsFactors = FALSE)
sig_disc_prots <- disc_deps %>% filter(fdr_disc < 0.05) %>% pull(protein)
valid_proteins <- intersect(sig_disc_prots, colnames(modeling_data))

# Map gene symbols and filter immunoglobulin/keratin blacklists
clean_uniprots <- sapply(strsplit(valid_proteins, ";"), `[`, 1) %>% gsub("-.*|\\..*", "", .)
mapped_symbols <- suppressMessages(suppressWarnings(
  mapIds(org.Hs.eg.db, keys = clean_uniprots, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
))
valid_symbols <- as.character(mapped_symbols)
valid_symbols[is.na(valid_symbols)] <- clean_uniprots[is.na(valid_symbols)]
lookup <- setNames(valid_symbols, valid_proteins)

blacklist_pattern <- "^(RPL|RPS|KRT|HBA|HBB|HBD|MYH|ACT|TUB|IGH|IGK|IGL)" 
valid_proteins <- valid_proteins[!grepl(blacklist_pattern, valid_symbols)]

discovery_set   <- modeling_data %>% filter(cohort == "Discovery")
replication_set <- modeling_data %>% filter(cohort == "Validation")

train_means <- sapply(discovery_set[valid_proteins], function(x) mean(x, na.rm = TRUE))
train_sds   <- sapply(discovery_set[valid_proteins], function(x) sd(x, na.rm = TRUE))
train_sds[train_sds == 0 | is.na(train_sds)] <- 1

train_set_scaled <- discovery_set
train_set_scaled[valid_proteins] <- scale(discovery_set[valid_proteins], center = train_means, scale = train_sds)

val_set_scaled <- replication_set
val_set_scaled[valid_proteins] <- scale(replication_set[valid_proteins], center = train_means, scale = train_sds)

# ==============================================================================
# Step 2: Fit Clinical Baseline Model
# ==============================================================================
cat("Fitting clinical baseline logistic model...\n")

clin_model <- glm(group ~ age + sex + BMI + hepatic_disease + kidney.function + CV_Metabolic_Cat, 
                  data = train_set_scaled, family = binomial)

train_set_scaled$Clin_Score <- predict(clin_model, train_set_scaled, type = "link")
val_set_scaled$Clin_Score   <- predict(clin_model, val_set_scaled, type = "link")

train_offset <- train_set_scaled$Clin_Score
val_offset   <- val_set_scaled$Clin_Score

# ==============================================================================
# Step 3: Cross-Validation and Biomarker Panel Selection
# ==============================================================================
cat("Executing cross-validation and feature selection...\n")

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

# Evaluate panel size via 10-fold cross-validation
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
  slice(1)

optimal_N <- optimal_model_auto$N
selected_panel <- ranked_proteins[1:optimal_N]
selected_symbols <- unname(lookup[selected_panel])

# Fit final logistic model with clinical offset
form_or <- as.formula(paste("group ~ offset(Clin_Score) +", paste(paste0("`", selected_panel, "`"), collapse = " + ")))
final_model <- glm(form_or, data = train_set_scaled, family = binomial)

pred_train_full <- as.numeric(predict(final_model, newdata = train_set_scaled, type = "response"))
pred_ext_full   <- as.numeric(predict(final_model, newdata = val_set_scaled, type = "response"))
pred_ext_clin   <- as.numeric(predict(clin_model, newdata = val_set_scaled, type = "response"))

roc_train_full <- roc(train_set_scaled$group, pred_train_full, quiet = TRUE, ci = TRUE)
roc_ext_full   <- roc(val_set_scaled$group, pred_ext_full, quiet = TRUE, ci = TRUE)
roc_ext_clin   <- roc(val_set_scaled$group, pred_ext_clin, quiet = TRUE, ci = TRUE)

delong_test <- roc.test(roc_ext_full, roc_ext_clin, method = "delong")
nri_idi_res <- calc_nri_idi(p_old = pred_ext_clin, p_new = pred_ext_full, y_true = val_set_scaled$group)

# Benchmarks with Random Forest and SVM
set.seed(42)
ml_cols <- c("group", BASE_COVARIATES, selected_panel)
rf_data_train <- train_set_scaled %>% dplyr::select(all_of(ml_cols))
rf_data_val   <- val_set_scaled   %>% dplyr::select(all_of(ml_cols))
colnames(rf_data_train) <- make.names(colnames(rf_data_train))
colnames(rf_data_val)   <- make.names(colnames(rf_data_val))

rf_model <- randomForest(group ~ ., data = rf_data_train, ntree = 500)
pred_ext_rf <- as.numeric(predict(rf_model, newdata = rf_data_val, type = "prob")[, "PD"])
roc_ext_rf  <- roc(val_set_scaled$group, pred_ext_rf, quiet = TRUE)

svm_model <- svm(group ~ ., data = rf_data_train, probability = TRUE)
svm_pred_obj <- predict(svm_model, newdata = rf_data_val, probability = TRUE)
svm_probs <- attr(svm_pred_obj, "probabilities")
pred_ext_svm_prob <- as.numeric(svm_probs[, "PD"])
roc_ext_svm <- roc(val_set_scaled$group, pred_ext_svm_prob, quiet = TRUE)

table2_df <- data.frame(
  Model = c("Clinical Baseline", sprintf("Calibrated Offset Model (Clin + %d DEPs)", optimal_N), "Random Forest (Benchmark)", "SVM (Benchmark)"),
  AUC = c(as.numeric(auc(roc_ext_clin)), as.numeric(auc(roc_ext_full)), as.numeric(auc(roc_ext_rf)), as.numeric(auc(roc_ext_svm))),
  AUC_95_CI = c(
    sprintf("%.3f–%.3f", roc_ext_clin$ci[1], roc_ext_clin$ci[3]),
    sprintf("%.3f–%.3f", roc_ext_full$ci[1], roc_ext_full$ci[3]),
    "-", "-"
  ),
  DeLong_P = c("-", sprintf("%.4e", delong_test$p.value), "-", "-"),
  NRI = c("-", sprintf("+%.3f (P=%.4e)", nri_idi_res$Estimate[1], nri_idi_res$P_Value[1]), "-", "-"),
  IDI = c("-", sprintf("+%.3f (P=%.4e)", nri_idi_res$Estimate[2], nri_idi_res$P_Value[2]), "-", "-")
)
write.csv(table2_df, file.path(output_dir, "Table_2_External_Validation_Metrics.csv"), row.names = FALSE)

final_summary <- broom::tidy(final_model, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(Odds_Ratio = estimate, P_Value = p.value, CI_Lower = conf.low, CI_Upper = conf.high) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

final_summary$Symbol <- sapply(final_summary$term, function(t) {
  clean_term <- gsub("`", "", t)
  if(clean_term %in% c("(Intercept)", BASE_COVARIATES)) return(clean_term)
  sym <- lookup[[clean_term]]
  if(is.na(sym) || is.null(sym)) return(clean_term) else return(sym)
})

final_table <- final_summary %>% 
  dplyr::select(Symbol, term, Odds_Ratio, CI_Lower, CI_Upper, P_Value) %>%
  rename(UNIPROT_ID = term) %>%
  mutate(UNIPROT_ID = gsub("`", "", UNIPROT_ID)) %>%
  arrange(P_Value)
write.csv(final_table, file.path(output_dir, "Table_ST14a_Diagnostic_Panel_OR_Parameters.csv"), row.names = FALSE)

# ==============================================================================
# Step 4: Generate Figure 2A (ROC Curves)
# ==============================================================================
cat("Generating Figure 2A (ROC curves)...\n")

pdf(file.path(output_dir, "Figure2A_ROC_Validation.pdf"), width = 6.8, height = 6.5)
par(mar = c(4.5, 4.5, 2.5, 2.0))

plot(roc_ext_full, col = "#e41a1c", lwd = 2.8, legacy.axes = FALSE,
     xlab = "Specificity", ylab = "Sensitivity",
     cex.lab = 1.15, font.lab = 2)
plot(roc_train_full, col = "#377eb8", lwd = 2.0, add = TRUE)
plot(roc_ext_clin, col = "grey55", lwd = 2.0, lty = 2, add = TRUE)
abline(a = 1, b = -1, col = "grey75", lwd = 1.2)

legend("bottomright", legend = c(
  sprintf("Validation Cohort: %.3f (%.3f–%.3f)", roc_ext_full$auc, roc_ext_full$ci[1], roc_ext_full$ci[3]),
  sprintf("Discovery Cohort: %.3f (%.3f–%.3f)", roc_train_full$auc, roc_train_full$ci[1], roc_train_full$ci[3]),
  sprintf("Clinical Baseline: %.3f (%.3f–%.3f)", roc_ext_clin$auc, roc_ext_clin$ci[1], roc_ext_clin$ci[3])
), col = c("#e41a1c", "#377eb8", "grey55"), lwd = c(2.8, 2.0, 2.0), lty = c(1, 1, 2), bty = "n", cex = 0.95)
dev.off()

# ==============================================================================
# Step 5: Generate Figure 2B (SHAP Beeswarm Plot)
# ==============================================================================
cat("Generating Figure 2B (SHAP beeswarm plot)...\n")

form_shap_pure <- as.formula(paste("group ~", paste(paste0("`", selected_panel, "`"), collapse = " + ")))
model_for_shap <- glm(form_shap_pure, data = train_set_scaled, family = binomial)

X_shap <- as.data.frame(train_set_scaled %>% dplyr::select(all_of(selected_panel)))
ks <- kernelshap(model_for_shap, X = X_shap, bg_X = X_shap)
sv <- shapviz(ks)
colnames(sv) <- unname(lookup[colnames(sv)])

pdf(file.path(output_dir, "Figure2B_SHAP_Summary.pdf"), width = 7.5, height = 7.2)
p_shap <- sv_importance(sv, kind = "beeswarm", max_display = optimal_N) + 
  theme_bw(base_size = 12) +
  labs(title = "Global Feature Importance", x = "SHAP value") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.text.y = element_text(face = "bold", color = "black", size = 11),
    axis.title = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )
print(p_shap)
dev.off()

# ==============================================================================
# Step 6: Generate Figure 2C (PPI Network) & Figure 2D (Validation Violins)
# ==============================================================================
cat("Generating Figure 2C (PPI connectivity) and Figure 2D (expression violins)...\n")

# PPI Network
genes_str_panel <- paste(selected_symbols, collapse = "%0d")
url_panel <- paste0("https://string-db.org/api/tsv/network?identifiers=", genes_str_panel, "&species=9606&add_nodes=15&required_score=400")

edges_df_panel <- tryCatch({ 
  read.table(url_panel, header = TRUE, sep = "\t", stringsAsFactors = FALSE) 
}, error = function(e) NULL)

if (!is.null(edges_df_panel) && nrow(edges_df_panel) > 0) {
  all_network_nodes <- unique(c(edges_df_panel$preferredName_A, edges_df_panel$preferredName_B, selected_symbols))
  nodes_df <- data.frame(name = all_network_nodes, stringsAsFactors = FALSE)
  links_df <- edges_df_panel %>% dplyr::select(from = preferredName_A, to = preferredName_B, score) %>% distinct()
  
  g_panel <- graph_from_data_frame(d = links_df, vertices = nodes_df, directed = FALSE)
  g_panel <- igraph::simplify(g_panel, remove.multiple = TRUE, remove.loops = TRUE)
  
  V(g_panel)$is_core <- V(g_panel)$name %in% selected_symbols
  V(g_panel)$degree  <- degree(g_panel)
  
  sub_g <- induced_subgraph(g_panel, degree(g_panel) > 0)
  if (ecount(sub_g) > 0) {
    comm_sub <- membership(cluster_louvain(sub_g))
    V(g_panel)$community <- NA 
    V(g_panel)$community[match(names(comm_sub), V(g_panel)$name)] <- as.numeric(comm_sub)
  } else {
    V(g_panel)$community <- NA
  }
  
  V(g_panel)$hull_label <- NA_character_
  communities <- unique(na.omit(V(g_panel)$community))
  
  for(i in communities) {
    nodes_in_comm <- V(g_panel)$name[which(V(g_panel)$community == i)]
    core_count_in_comm <- sum(V(g_panel)$is_core[which(V(g_panel)$community == i)])
    
    if(length(nodes_in_comm) >= 3 && core_count_in_comm >= 1) {
      g_ids <- tryCatch({ 
        bitr(nodes_in_comm, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID 
      }, error = function(e) NULL)
      
      if(!is.null(g_ids) && length(g_ids) >= 2) {
        ego <- suppressMessages(enrichGO(g_ids, OrgDb=org.Hs.eg.db, ont="BP", pvalueCutoff=0.5, readable=TRUE))
        if(!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
          V(g_panel)$hull_label[which(V(g_panel)$community == i)] <- stringr::str_wrap(stringr::str_to_title(ego@result$Description[1]), width = 16)
        }
      }
    }
  }
  
  tg_panel <- as_tbl_graph(g_panel)
  
  pdf(file.path(output_dir, "Figure2C_Functional_Connectivity_Diagnostic_Panel.pdf"), width = 11.5, height = 8.5)
  set.seed(88)
  p_fig2c <- ggraph(tg_panel, layout = "fr") +
    geom_mark_hull(aes(x = x, y = y, fill = hull_label, label = hull_label, filter = !is.na(hull_label)),
                   concavity = 3.5, expand = unit(4.5, "mm"), alpha = 0.14,
                   color = "grey55", linetype = "dashed", linewidth = 0.6,
                   label.fontsize = 10.5, label.fontface = "bold", label.fill = "white",
                   label.buffer = unit(2.5, "mm"), con.colour = "grey50", con.type = "straight") +
    scale_fill_brewer(palette = "Set2", na.translate = FALSE) +
    new_scale_fill() +
    geom_edge_link(color = "grey75", width = 0.8, alpha = 0.7) +
    geom_node_point(aes(filter = !is_core), size = 3.8, color = "grey50", fill = "grey90", shape = 21, stroke = 0.8) +
    geom_node_point(aes(filter = is_core), size = 8.5, color = "white", fill = "#d73027", shape = 21, stroke = 1.2) +
    geom_node_text(aes(filter = !is_core, label = name), repel = TRUE, size = 3.4, color = "grey35") +
    geom_node_text(aes(filter = is_core, label = name), repel = TRUE, size = 4.8, color = "black", fontface = "bold",
                   bg.color = "white", bg.r = 0.12) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15, margin = ggplot2::margin(b = 10)),
      plot.margin = ggplot2::margin(0.8, 0.8, 0.8, 0.8, "cm")
    ) +
    labs(title = "Functional Connectivity of the Diagnostic Biomarker Panel",
         subtitle = "Red: Core Diagnostic Biomarkers; Grey: First-shell Interactors; Dashed: Biological Modules")
  print(p_fig2c)
  dev.off()
}

# Expression violins for top panel biomarkers in validation cohort
top_8_features <- selected_panel[1:min(8, optimal_N)]
top_8_symbols  <- unname(lookup[top_8_features])

plot_val_df <- val_set_scaled %>%
  dplyr::select(group, all_of(top_8_features)) %>%
  pivot_longer(cols = -group, names_to = "Protein_ID", values_to = "Relative_Abundance") %>%
  mutate(Symbol = factor(lookup[Protein_ID], levels = top_8_symbols))

p_fig2d <- ggplot(plot_val_df, aes(x = group, y = Relative_Abundance, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.45, color = NA, scale = "width", width = 0.85) +
  geom_boxplot(width = 0.16, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.65) +
  scale_fill_manual(values = c("HC" = "#9ecae1", "PD" = "#fc9272"), guide = "none") +
  facet_wrap(~ Symbol, scales = "free_y", nrow = 2) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5, size = 4.8, fontface = "bold", color = "black") +
  labs(title = "Expression of Core Proteins (Validation Cohort)",
       y = "Relative Protein Abundance (Z-score)", x = NULL) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13.5, margin = ggplot2::margin(b = 10)),
    strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.7),
    strip.text = element_text(face = "bold", size = 11, color = "black"),
    axis.text.x = element_text(face = "bold", size = 11, color = "black"),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(face = "bold", size = 11.5, margin = ggplot2::margin(r = 8)),
    panel.grid.major.y = element_line(linetype = "dotted", color = "grey85", linewidth = 0.35)
  )

ggsave(file.path(output_dir, "Figure2D_Core_Proteins_Violin.pdf"), p_fig2d, width = 12.0, height = 6.5, device = cairo_pdf)
ggsave(file.path(output_dir, "Figure2D_Core_Proteins_Violin.png"), p_fig2d, width = 12.0, height = 6.5, dpi = 300)

saveRDS(final_model, file.path(output_dir, "lasso_diagnostic_model.rds"))

cat("Analysis complete. Figure 2 panels and model parameters saved to:", output_dir, "\n")