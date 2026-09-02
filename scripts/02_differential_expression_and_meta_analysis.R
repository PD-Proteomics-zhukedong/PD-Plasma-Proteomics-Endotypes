# ==============================================================================
# Script: 02_differential_expression_and_meta_analysis.R
# Purpose: Limma DEA (6 Covariates), Stouffer Meta-Analysis, and Figure 1 (B-G)
# Project: Large-scale plasma proteomics in Parkinson's disease
# ==============================================================================

options(expressions = 5000)
set.seed(42)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(limma)
  library(patchwork)
  library(ggrepel)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ReactomePA)
  library(ggnewscale)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggforce)
  library(GSVA)
  library(ggpubr)
  library(broom)
  library(scales)
  library(stringr)
})

# Directory setup
input_dir  <- "results/01_preprocessed_data"
data_dir   <- "data"
output_dir <- "results/02_figure1_meta_dea"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

FDR_CUTOFF   <- 0.05
P_CUTOFF     <- 0.05
FC_THRESHOLD <- log2(1.5)
BASE_COVARIATES <- c("age", "sex", "BMI", "hepatic_disease", "kidney.function", "CV_Metabolic_Cat")

# Helper functions
get_mode <- function(x) {
  x <- na.omit(x)
  if(length(x) == 0) return(NA)
  uniqv <- unique(x)
  uniqv[which.max(tabulate(match(x, uniqv)))]
}

format_ontology_label <- function(text_vec, max_width = 34) {
  if (length(text_vec) == 0) return(character(0))
  text_vec %>%
    as.character() %>%
    stringr::str_to_title() %>%
    stringr::str_wrap(width = max_width)
}

cat("\n=== Phase 2: Differential Expression and Figure 1 Pipeline ===\n")

# ==============================================================================
# Step 1: Clinical Metadata Harmonization
# ==============================================================================
cat("Harmonizing clinical covariates per cohort...\n")

clin_raw_file <- file.path(data_dir, "metadata_clinical_n1119.csv")
clinical_raw  <- read.csv(clin_raw_file, stringsAsFactors = FALSE)
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
  mutate(group = ifelse(group %in% c("CTR", "Control", "HC"), "HC", "PD")) %>%
  group_by(cohort) %>%
  mutate(
    age = as.numeric(ifelse(is.na(age), median(age, na.rm = TRUE), age)),
    sex = ifelse(sex %in% c("male", "M", "1", 1), 1, 0),
    BMI = as.numeric(ifelse(is.na(BMI), median(BMI, na.rm = TRUE), BMI)),
    hepatic_disease = as.numeric(as.character(ifelse(is.na(hepatic_disease), get_mode(hepatic_disease), hepatic_disease))),
    kidney.function = as.numeric(as.character(ifelse(is.na(kidney.function), get_mode(kidney.function), kidney.function))),
    CV_Metabolic_Score = case_when(
      is.na(Metabolic.cardiovascular.comorbidities) ~ as.character(get_mode(Metabolic.cardiovascular.comorbidities)),
      Metabolic.cardiovascular.comorbidities >= 2 ~ "2+",
      TRUE ~ as.character(Metabolic.cardiovascular.comorbidities)
    )
  ) %>%
  ungroup() %>%
  mutate(
    group = factor(group, levels = c("HC", "PD")), 
    sex = as.factor(sex),
    hepatic_disease = as.factor(hepatic_disease),
    kidney.function = as.factor(kidney.function),
    CV_Metabolic_Cat = factor(CV_Metabolic_Score, levels = c("0", "1", "2+")),
    cohort = as.factor(cohort)
  )

# ==============================================================================
# Step 2: Limma Modeling and Meta-Analysis
# ==============================================================================
cat("Running linear modeling with empirical Bayes (limma) and meta-analysis...\n")

combat_raw <- fread(file.path(input_dir, "3_ComBat_Corrected_Matrix.tsv"), data.table = FALSE)
rownames(combat_raw) <- make.unique(as.character(combat_raw$Protein.Group))
expr_combat_all <- as.matrix(combat_raw[, 7:ncol(combat_raw)])

train_clin <- clinical_clean %>% filter(cohort == "Discovery" & sample %in% colnames(expr_combat_all))
val_clin   <- clinical_clean %>% filter(cohort == "Validation" & sample %in% colnames(expr_combat_all))

run_limma_analysis <- function(expr_mat, clin_df, base_covariates) {
  common_samps <- intersect(colnames(expr_mat), clin_df$sample)
  expr_sub <- expr_mat[, common_samps]
  clin_sub <- clin_df[match(common_samps, clin_df$sample), ]
  
  valid_covs <- c()
  for (cov in base_covariates) {
    if (cov %in% colnames(clin_sub)) {
      vals <- na.omit(clin_sub[[cov]])
      if (length(unique(vals)) > 1) valid_covs <- c(valid_covs, cov)
    }
  }
  form_str <- paste0("~ group", if(length(valid_covs) > 0) paste0(" + ", paste(valid_covs, collapse = " + ")) else "")
  design <- model.matrix(as.formula(form_str), data = clin_sub)
  
  fit <- lmFit(expr_sub, design)
  fit <- eBayes(fit)
  res_tab <- topTable(fit, coef = "groupPD", number = Inf, sort.by = "none")
  res_tab$protein <- rownames(res_tab)
  
  res_tab %>% dplyr::select(protein, estimate = logFC, p.value = P.Value, fdr = adj.P.Val, t_stat = t)
}

limma_disc <- run_limma_analysis(expr_combat_all, train_clin, BASE_COVARIATES) %>%
  rename(est_disc = estimate, p_disc = p.value, fdr_disc = fdr, t_disc = t_stat) %>%
  mutate(Status_Disc = case_when(
    fdr_disc < FDR_CUTOFF & est_disc > FC_THRESHOLD ~ "Up in Discovery",
    fdr_disc < FDR_CUTOFF & est_disc < -FC_THRESHOLD ~ "Down in Discovery",
    TRUE ~ "NS"
  ))
write.csv(limma_disc, file.path(output_dir, "Table_1_Discovery_Limma_Full_Model.csv"), row.names = FALSE)

limma_rep <- run_limma_analysis(expr_combat_all, val_clin, BASE_COVARIATES) %>%
  rename(est_rep = estimate, p_rep = p.value, fdr_rep = fdr, t_rep = t_stat) %>%
  mutate(Status_Rep = case_when(
    p_rep < P_CUTOFF & est_rep > 0 ~ "Up in Validation",
    p_rep < P_CUTOFF & est_rep < 0 ~ "Down in Validation",
    TRUE ~ "NS"
  ))
write.csv(limma_rep, file.path(output_dir, "Table_2_Validation_Limma_Full_Model.csv"), row.names = FALSE)

N_disc <- nrow(train_clin)
N_rep  <- nrow(val_clin)

meta_final <- limma_disc %>%
  inner_join(limma_rep, by = "protein") %>%
  mutate(
    p_disc_safe = pmin(pmax(p_disc, 1e-300), 1 - 1e-16),
    p_rep_safe  = pmin(pmax(p_rep, 1e-300), 1 - 1e-16),
    z_disc = sign(est_disc) * qnorm(p_disc_safe / 2, lower.tail = FALSE),
    z_rep  = sign(est_rep) * qnorm(p_rep_safe / 2, lower.tail = FALSE),
    z_meta = (sqrt(N_disc) * z_disc + sqrt(N_rep) * z_rep) / sqrt(N_disc + N_rep),
    meta_p = pmax(2 * pnorm(abs(z_meta), lower.tail = FALSE), 1e-300),
    meta_fdr = p.adjust(meta_p, method = "BH"),
    meta_est = (est_disc + est_rep) / 2
  ) %>%
  mutate(Final_Status = case_when(
    meta_fdr < FDR_CUTOFF & est_disc > FC_THRESHOLD & est_rep > 0 ~ "Strictly Validated Up",
    meta_fdr < FDR_CUTOFF & est_disc < -FC_THRESHOLD & est_rep < 0 ~ "Strictly Validated Down",
    TRUE ~ "Filtered Out"
  ))
write.csv(meta_final, file.path(output_dir, "Table_3_Final_Strict_DEPs_StoufferMeta.csv"), row.names = FALSE)

n_meta_up   <- sum(meta_final$Final_Status == "Strictly Validated Up")
n_meta_down <- sum(meta_final$Final_Status == "Strictly Validated Down")
cat(sprintf("Identified %d Meta-DEPs (Up: %d, Down: %d).\n", n_meta_up + n_meta_down, n_meta_up, n_meta_down))

# ==============================================================================
# Step 3: Generate Figure 1B-D (Volcano Plots)
# ==============================================================================
cat("Generating Figure 1B-D volcano plots...\n")

clean_uniprots <- sapply(strsplit(meta_final$protein, ";"), `[`, 1) %>% gsub("-.*|\\..*", "", .)
mapped_symbols <- suppressMessages(suppressWarnings(
  mapIds(org.Hs.eg.db, keys = clean_uniprots, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
))
mapped_symbols[is.na(mapped_symbols)] <- clean_uniprots[is.na(mapped_symbols)]
lookup <- setNames(as.character(mapped_symbols), meta_final$protein)
meta_final$Symbol <- unname(lookup[meta_final$protein])

# Identify top representative biomarkers for labeling
blacklist_regex <- "^(RPL|RPS|KRT|HBA|HBB|HBD|MYH|ACT|TUB|IGH|IGK|IGL|ALB)"
valid_meta_candidates <- meta_final %>% filter(!grepl(blacklist_regex, Symbol))

top_meta_up <- valid_meta_candidates %>%
  filter(Final_Status == "Strictly Validated Up" & p_disc < 1e-4 & p_rep < 0.01) %>%
  arrange(desc(abs(meta_est)), meta_fdr) %>%
  slice_head(n = 5) %>%
  pull(protein)

top_meta_down <- valid_meta_candidates %>%
  filter(Final_Status == "Strictly Validated Down" & p_disc < 1e-4 & p_rep < 0.01) %>%
  arrange(desc(abs(meta_est)), meta_fdr) %>%
  slice_head(n = 5) %>%
  pull(protein)

landmark_proteins <- c(top_meta_up, top_meta_down)
landmark_symbols  <- unname(lookup[landmark_proteins])
cat(sprintf("Landmark biomarkers for annotation:\n  Up: %s\n  Down: %s\n", 
            paste(landmark_symbols[1:5], collapse = ", "), 
            paste(landmark_symbols[6:10], collapse = ", ")))

limma_disc_plot <- limma_disc %>%
  mutate(log10_fdr = -log10(pmin(pmax(fdr_disc, 1e-85), 1)), Symbol = unname(lookup[protein]))

limma_rep_plot <- limma_rep %>%
  mutate(log10_p = -log10(pmin(pmax(p_rep, 1e-35), 1)), Symbol = unname(lookup[protein]))

disc_p_safe <- pmin(pmax(limma_disc$p_disc, 1e-300), 1 - 1e-16)
rep_p_safe  <- pmin(pmax(limma_rep$p_rep, 1e-300), 1 - 1e-16)
z1 <- sign(limma_disc$est_disc) * qnorm(disc_p_safe / 2, lower.tail = FALSE)
z2 <- sign(limma_rep$est_rep)   * qnorm(rep_p_safe / 2, lower.tail = FALSE)
z_meta_vec <- (sqrt(N_disc) * z1 + sqrt(N_rep) * z2) / sqrt(N_disc + N_rep)

log_p_raw <- pnorm(abs(z_meta_vec), lower.tail = FALSE, log.p = TRUE) + log(2)
meta_log10_p <- - (log_p_raw / log(10))
meta_log10_fdr_smooth <- meta_log10_p - log10(length(z_meta_vec) / rank(log_p_raw))

meta_final_plot <- meta_final
meta_final_plot$meta_log10_fdr <- pmax(meta_log10_fdr_smooth, 0)

n_disc_up   <- sum(limma_disc$Status_Disc == "Up in Discovery")
n_disc_down <- sum(limma_disc$Status_Disc == "Down in Discovery")
n_rep_up    <- sum(limma_rep$Status_Rep == "Up in Validation")
n_rep_down  <- sum(limma_rep$Status_Rep == "Down in Validation")

plot_volcano_publication <- function(df, x_col, y_col, status_col, title_text, y_title_expr, target_prots, 
                                     n_down, n_up, fc_cut = NULL, p_cut_val = 0.05, y_limit = 70) {
  df$X_val <- df[[x_col]]
  df$Y_val <- df[[y_col]]
  df$Status <- df[[status_col]]
  df$is_highlight <- df$protein %in% target_prots
  
  ggplot(df, aes(x = X_val, y = Y_val)) +
    geom_hline(yintercept = -log10(p_cut_val), linetype = "dashed", color = "grey75", linewidth = 0.5) +
    {if(!is.null(fc_cut)) geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "grey75", linewidth = 0.5)} +
    geom_point(data = filter(df, !grepl("Up|Down", Status)), color = "#E0E0E0", size = 1.3, alpha = 0.45) +
    geom_point(data = filter(df, grepl("Down", Status)), color = "#2B5C8F", size = 2.0, alpha = 0.7) +
    geom_point(data = filter(df, grepl("Up", Status)), color = "#D73027", size = 2.0, alpha = 0.7) +
    geom_point(data = filter(df, is_highlight), shape = 21, color = "black", fill = "gold", size = 3.6, stroke = 1.3) +
    geom_text_repel(data = filter(df, is_highlight), aes(label = Symbol),
                    size = 3.8, fontface = "bold.italic", box.padding = 0.5, point.padding = 0.3, 
                    max.overlaps = 50, segment.color = "grey30", segment.size = 0.35, bg.color = "white", bg.r = 0.1) +
    annotate("label", x = -Inf, y = Inf, label = sprintf("Down: %d", n_down), 
             hjust = -0.15, vjust = 1.35, color = "#2B5C8F", fill = alpha("white", 0.88), fontface = "bold", size = 4.2) +
    annotate("label", x = Inf, y = Inf, label = sprintf("Up: %d", n_up), 
             hjust = 1.15, vjust = 1.35, color = "#D73027", fill = alpha("white", 0.88), fontface = "bold", size = 4.2) +
    theme_classic(base_size = 13) +
    labs(title = title_text, x = expression(bold(log[2]("Fold Change"))), y = y_title_expr) +
    coord_cartesian(ylim = c(0, y_limit)) +
    theme(
      plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(color = "black", size = 11),
      panel.grid.major = element_line(color = "grey95", linetype = "dashed")
    )
}

p1 <- plot_volcano_publication(limma_disc_plot, "est_disc", "log10_fdr", "Status_Disc", 
                               "Discovery Cohort", expression(bold(-log[10]("FDR"))), 
                               landmark_proteins, n_disc_down, n_disc_up, FC_THRESHOLD, FDR_CUTOFF, 70)

p2 <- plot_volcano_publication(limma_rep_plot, "est_rep", "log10_p", "Status_Rep", 
                               "Validation Cohort", expression(bold(-log[10]("P-value"))), 
                               landmark_proteins, n_rep_down, n_rep_up, NULL, P_CUTOFF, 25)

p3 <- plot_volcano_publication(meta_final_plot, "meta_est", "meta_log10_fdr", "Final_Status", 
                               "Stouffer Meta-Analysis", expression(bold(-log[10]("Meta-FDR"))), 
                               landmark_proteins, n_meta_down, n_meta_up, FC_THRESHOLD, FDR_CUTOFF, 45)

final_volcano <- p1 | p2 | p3
ggsave(file.path(output_dir, "Figure1BCD_Volcano_Panels.pdf"), final_volcano, width = 16, height = 5.3, device = cairo_pdf)
ggsave(file.path(output_dir, "Figure1BCD_Volcano_Panels.png"), final_volcano, width = 16, height = 5.3, dpi = 300)

# ==============================================================================
# Step 4: Generate Figure 1E (GO-BP Functional Landscape)
# ==============================================================================
cat("Performing GO-BP enrichment and generating Figure 1E...\n")

tested_uniprots <- unique(sapply(strsplit(rownames(expr_combat_all), ";"), `[`, 1)) %>% gsub("-.*|\\..*", "", .)
bg_map <- suppressMessages(suppressWarnings(
  bitr(tested_uniprots, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
))
universe_entrez <- unique(bg_map$ENTREZID)

gold_deps <- meta_final %>% filter(Final_Status %in% c("Strictly Validated Up", "Strictly Validated Down"))
gold_deps$Clean_ID <- sapply(strsplit(gold_deps$protein, ";"), `[`, 1) %>% gsub("-.*|\\..*", "", .)
gene_map <- suppressMessages(suppressWarnings(
  bitr(gold_deps$Clean_ID, fromType="UNIPROT", toType="ENTREZID", OrgDb=org.Hs.eg.db)
))
gold_mapped <- gold_deps %>% inner_join(gene_map, by = c("Clean_ID" = "UNIPROT"))

gene_clusters <- list(
  "Up-regulated"   = unique(gold_mapped$ENTREZID[gold_mapped$Final_Status == "Strictly Validated Up"]),
  "Down-regulated" = unique(gold_mapped$ENTREZID[gold_mapped$Final_Status == "Strictly Validated Down"])
)

comp_go <- compareCluster(geneCluster = gene_clusters, fun = "enrichGO", OrgDb = org.Hs.eg.db, ont = "BP", 
                          pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500, 
                          universe = universe_entrez, readable = TRUE)
df_go <- as.data.frame(comp_go)
write.csv(df_go, file.path(output_dir, "Table_ST11_GO_BP_MetaDEPs_Full.csv"), row.names = FALSE)

# Functional ontology themes
go_theme_dict <- list(
  "Mitochondria &\nEnergy" = c("atp synthesis", "oxidative phosphorylation", "respiration", "electron transport", "mitochondri", "energy"),
  "Neurotrophic &\nSynaptic Signaling" = c("golgi vesicle", "post-golgi", "vesicle organization", "synap", "neuro", "axon", "transmembrane transport", "retromer", "trafficking"),
  "Autophagy &\nProteostasis" = c("proteasomal", "ubiquitin", "deubiquitin", "starvation", "autophag", "lysosom", "catabolic process", "proteoly"),
  "Immunity" = c("antigen processing", "antigen presentation", "cytotoxicity", "immune", "complement", "t cell", "leukocyte", "cytokine"),
  "Platelet &\nCoagulation" = c("fibrinolysis", "platelet", "coagulat", "hemostasis", "blood coagulation", "clot")
)

theme_colors_nature <- c(
  "Mitochondria &\nEnergy"             = "#2ca02c", 
  "Neurotrophic &\nSynaptic Signaling" = "#1f78b4", 
  "Autophagy &\nProteostasis"          = "#6a3d9a", 
  "Immunity"                           = "#e31a1c", 
  "Platelet &\nCoagulation"            = "#ff7f00"
)

processed_all <- df_go %>%
  mutate(logP = -log10(pvalue)) %>%
  rename_with(~"Cluster_Raw", contains("Cluster")) %>%
  mutate(Theme = "Others")

for(thm in names(go_theme_dict)) {
  keys <- go_theme_dict[[thm]]
  match_idx <- str_detect(tolower(processed_all$Description), paste(keys, collapse="|"))
  processed_all$Theme[match_idx] <- thm
}

# Select top pathways per functional theme
selected_df <- processed_all %>%
  filter(Theme != "Others" & pvalue < 0.05) %>%
  group_by(Theme) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  mutate(Display_Name = format_ontology_label(Description, max_width = 34)) %>%
  distinct(Display_Name, .keep_all = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup()

plot_data_final <- expand.grid(
  Description = unique(selected_df$Description),
  Cluster_Raw = c("Up-regulated", "Down-regulated"),
  stringsAsFactors = FALSE
) %>%
  left_join(selected_df %>% dplyr::select(Description, Theme) %>% distinct(), by = "Description") %>%
  left_join(processed_all %>% dplyr::select(Description, Cluster_Raw, logP, GeneRatio, pvalue), by = c("Description", "Cluster_Raw")) %>%
  mutate(
    logP = replace_na(logP, 0),
    GeneRatio = replace_na(GeneRatio, "0/1"),
    Display_Name = format_ontology_label(Description, max_width = 34),
    Cluster_Plot = factor(Cluster_Raw, levels = c("Up-regulated", "Down-regulated"), labels = c("UP", "DOWN"))
  ) %>%
  separate(GeneRatio, into = c("n", "N"), sep = "/") %>%
  mutate(
    GeneRatio_num = as.numeric(n) / as.numeric(N),
    Theme = factor(Theme, levels = names(theme_colors_nature))
  ) %>%
  filter(!is.na(Theme))

order_df <- plot_data_final %>% 
  group_by(Display_Name) %>% 
  summarise(Theme = dplyr::first(Theme), MaxP = max(logP, na.rm = TRUE)) %>% 
  arrange(Theme, MaxP)

plot_data_final$Display_Name <- factor(plot_data_final$Display_Name, levels = unique(order_df$Display_Name))

p_fig1e <- ggplot(plot_data_final, aes(x = Cluster_Plot, y = Display_Name)) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "grey88", linewidth = 0.75) +
  geom_tile(aes(x = 2.28, y = Display_Name, fill = Theme), width = 0.05, show.legend = FALSE) +
  scale_fill_manual(values = theme_colors_nature, guide = "none") +
  new_scale_fill() +
  geom_point(aes(size = GeneRatio_num, fill = logP, alpha = logP > 0), shape = 21, color = "black", stroke = 0.45) +
  scale_alpha_manual(values = c(0.08, 1), guide = "none") +
  scale_fill_gradientn(colors = c("#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026"), 
                       name = expression(bold(-log[10](FDR)))) +
  scale_size_continuous(range = c(2.8, 6.8), name = "Gene Ratio", breaks = c(0, 0.05, 0.10, 0.15)) +
  facet_grid(Theme ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
  scale_x_discrete(position = "top", expand = expansion(mult = c(0.2, 0.2))) +
  labs(title = "Functional Landscape of PD Plasma", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.spacing.y = unit(0.35, "lines"),
    strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 9.0, color = "grey15"),
    strip.background = element_blank(),
    axis.text.y = element_text(color = "black", size = 9.0, hjust = 1, lineheight = 0.9),
    axis.text.x = element_text(face = "bold", size = 11, color = "black"),
    panel.grid.major.y = element_line(linetype = "dotted", color = "grey90", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    plot.margin = ggplot2::margin(10, 15, 10, 10),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
  )

ggsave(file.path(output_dir, "Figure1E_Functional_Landscape.pdf"), p_fig1e, width = 7.8, height = 7.5, device = cairo_pdf)

# ==============================================================================
# Step 5: Generate Figure 1F (Louvain PPI Hub Network)
# ==============================================================================
cat("Constructing PPI network and generating Figure 1F...\n")

gold_deps_df <- meta_final %>% 
  filter(Final_Status %in% c("Strictly Validated Up", "Strictly Validated Down")) %>%
  mutate(Clean_UNIPROT = sapply(strsplit(protein, ";"), `[`, 1) %>% gsub("-.*|\\..*", "", .))

gene_map_ppi <- suppressMessages(suppressWarnings(
  bitr(gold_deps_df$Clean_UNIPROT, fromType="UNIPROT", toType="SYMBOL", OrgDb=org.Hs.eg.db)
)) %>% distinct(UNIPROT, .keep_all = TRUE)

genes_str <- paste(unique(gene_map_ppi$SYMBOL)[1:min(400, length(unique(gene_map_ppi$SYMBOL)))], collapse = "%0d")
url <- paste0("https://string-db.org/api/tsv/network?identifiers=", genes_str, "&species=9606&required_score=700")

edges_df <- tryCatch({ read.table(url, header = TRUE, sep = "\t", stringsAsFactors = FALSE) }, error = function(e) NULL)

if (!is.null(edges_df) && nrow(edges_df) > 0) {
  sub_links <- edges_df %>% dplyr::select(from = preferredName_A, to = preferredName_B, score = score) %>% distinct()
  g_initial <- simplify(graph_from_data_frame(sub_links, directed = FALSE))
  
  degree_scores      <- degree(g_initial)
  betweenness_scores <- betweenness(g_initial, directed = FALSE, normalized = TRUE)
  
  mcc_scores <- setNames(numeric(vcount(g_initial)), V(g_initial)$name)
  clqs <- cliques(g_initial, min = 2, max = 5)
  for (clq in clqs) {
    size <- length(clq)
    mcc_scores[names(clq)] <- mcc_scores[names(clq)] + factorial(size - 1)
  }
  
  hub_comparison_df <- data.frame(
    Symbol = V(g_initial)$name,
    Degree = as.numeric(degree_scores),
    Betweenness = as.numeric(betweenness_scores),
    MCC = as.numeric(mcc_scores)
  ) %>%
    mutate(
      Rank_Deg = min_rank(desc(Degree)),
      Rank_Bet = min_rank(desc(Betweenness)),
      Rank_MCC = min_rank(desc(MCC)),
      Composite_Rank_Score = (Rank_Deg + Rank_Bet + Rank_MCC) / 3
    ) %>%
    left_join(gene_map_ppi, by = c("Symbol" = "SYMBOL")) %>%
    left_join(gold_deps_df, by = c("UNIPROT" = "Clean_UNIPROT")) %>%
    arrange(Composite_Rank_Score) %>%
    dplyr::select(Symbol, UNIPROT, Composite_Rank_Score, Degree, MCC, Betweenness, Final_Status, meta_est)
  
  write.csv(hub_comparison_df, file.path(output_dir, "Table_4_Full_PPI_Algorithm_Comparison.csv"), row.names = FALSE)
  
  top_30_candidates <- hub_comparison_df %>%
    filter(!str_detect(Symbol, blacklist_regex)) %>%
    arrange(Composite_Rank_Score) %>%
    slice_head(n = 30)
  
  sub_links_30 <- sub_links %>% filter(from %in% top_30_candidates$Symbol & to %in% top_30_candidates$Symbol)
  g_sub <- graph_from_data_frame(sub_links_30, directed = FALSE, vertices = top_30_candidates)
  g_sub <- induced_subgraph(g_sub, degree(g_sub) > 0)
  
  set.seed(42)
  comm_res <- cluster_louvain(g_sub)
  V(g_sub)$community_id <- as.character(membership(comm_res))
  
  # Community functional annotation via GO enrichment
  community_labels <- c()
  for (cid in unique(V(g_sub)$community_id)) {
    genes_in_comm <- V(g_sub)$name[V(g_sub)$community_id == cid]
    if (length(genes_in_comm) >= 2) {
      g_entrez <- suppressMessages(suppressWarnings(bitr(genes_in_comm, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID))
      ego <- tryCatch({ suppressMessages(enrichGO(g_entrez, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.8)) }, error = function(e) NULL)
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        community_labels[cid] <- format_ontology_label(as.data.frame(ego)$Description[1], max_width = 18)
      } else {
        community_labels[cid] <- paste("Cluster", cid)
      }
    } else {
      community_labels[cid] <- paste("Cluster", cid)
    }
  }
  V(g_sub)$Module_Annotation <- community_labels[V(g_sub)$community_id]
  
  set.seed(88)
  layout_raw <- layout_with_fr(g_sub, niter = 3500)
  layout_df <- as.data.frame(layout_raw) %>% rename(x = 1, y = 2) %>% mutate(community = V(g_sub)$community_id)
  cluster_centers <- layout_df %>% group_by(community) %>% summarise(cx = mean(x), cy = mean(y))
  
  layout_compact <- layout_df %>%
    left_join(cluster_centers, by = "community") %>%
    mutate(new_x = cx * 0.46 + (x - cx), new_y = cy * 0.46 + (y - cy)) %>%
    dplyr::select(new_x, new_y) %>% as.matrix()
  
  p_fig1f <- ggraph(as_tbl_graph(g_sub), layout = layout_compact) +
    geom_mark_hull(aes(x, y, group = Module_Annotation, fill = Module_Annotation, label = Module_Annotation),
                   color = "grey55", linetype = "dashed", linewidth = 0.55,
                   concavity = 1.0, expand = unit(3.0, "mm"), radius = unit(2.8, "mm"),
                   alpha = 0.16, show.legend = FALSE, 
                   label.fontsize = 9.5, label.fontface = "bold", label.fill = "white",
                   label.buffer = unit(2.0, "mm"), con.type = "straight", con.colour = "grey50") +
    scale_fill_brewer(palette = "Set2") + 
    new_scale_fill() + 
    geom_edge_link(color = "grey70", width = 0.75, alpha = 0.65) +
    geom_node_point(aes(fill = meta_est, size = Degree), shape = 21, stroke = 1.15, color = "white") +
    scale_fill_gradient2(low = "#313695", mid = "white", high = "#a50026", midpoint = 0, name = "Effect Size\n(Log2FC)") +
    scale_size_continuous(range = c(5.2, 10.5), name = "PPI Degree", breaks = c(10, 15, 20, 25, 30, 35)) +
    geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 3.8, bg.color = "white", bg.r = 0.12, color = "black") +
    theme_graph() +
    theme(legend.position = "right", plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, "cm"))
  
  ggsave(file.path(output_dir, "Figure1F_PPI_Network.pdf"), p_fig1f, width = 8.5, height = 7.5, device = cairo_pdf)
}

# ==============================================================================
# Step 6: Generate Figure 1G (ssGSEA Module Burden Scores)
# ==============================================================================
cat("Calculating ssGSEA module burden scores and generating Figure 1G...\n")

theme_specs <- list(
  "Mitochondria & Energy"       = list(pat = "mitochondr|atp synthesis|electron transport|respiration", dir = "Up-regulated"),
  "Neurotrophic & Synaptic"     = list(pat = "synap|neuro|axon|dendrit|golgi transport|snare|transmembrane transport", dir = "Up-regulated"),
  "Protein Modification"        = list(pat = "proteasom|ubiquitin|deubiquitin|autophag|lysosom|protein modification", dir = "Up-regulated"),
  "Immune Activation"           = list(pat = "immune|cytotoxicity|antigen presentation|complement", dir = "Up-regulated"),
  "Platelet & Coagulation"      = list(pat = "platelet|coagulat|hemostasis|fibrinoly", dir = "Down-regulated")
)

clean_expr_symbols <- combat_raw[, 7:ncol(combat_raw)]
clean_expr_symbols$SYMBOL <- unname(lookup[rownames(combat_raw)])
expr_symbol_mat <- clean_expr_symbols %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), ~mean(., na.rm = TRUE))) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

dep_module_genes <- list()
for(thm in names(theme_specs)) {
  spec <- theme_specs[[thm]]
  matched_rows <- df_go %>% filter(str_detect(tolower(Description), spec$pat) & Cluster == spec$dir)
  if(nrow(matched_rows) == 0) matched_rows <- df_go %>% filter(str_detect(tolower(Description), spec$pat))
  if (nrow(matched_rows) > 0) {
    genes_in_theme <- unique(unlist(strsplit(matched_rows$geneID, "/")))
    dep_module_genes[[thm]] <- intersect(genes_in_theme, rownames(expr_symbol_mat))
  }
}

param <- ssgseaParam(exprData = expr_symbol_mat, geneSets = dep_module_genes, minSize = 3)
burden_res <- gsva(param)

score_df <- as.data.frame(t(burden_res)) %>% 
  rownames_to_column("sample") %>%
  inner_join(clinical_clean, by = "sample") %>%
  pivot_longer(cols = all_of(names(dep_module_genes)), names_to = "Module", values_to = "Burden_Score")

score_df$group  <- factor(score_df$group, levels = c("HC", "PD"))
score_df$Module <- factor(score_df$Module, levels = names(theme_specs))

# Fit 6-covariate adjusted linear regression models
cov_adjusted_results <- map_df(names(dep_module_genes), function(mod) {
  sub_df <- score_df %>% filter(Module == mod)
  fit_adj <- lm(Burden_Score ~ group + age + sex + BMI + hepatic_disease + kidney.function + CV_Metabolic_Cat, data = sub_df)
  tidy_res <- broom::tidy(fit_adj) %>% filter(term == "groupPD")
  wilcox_p <- wilcox.test(Burden_Score ~ group, data = sub_df)$p.value
  data.frame(Module = mod, Raw_Wilcox_P = wilcox_p, Adjusted_LM_Beta = tidy_res$estimate, Adjusted_LM_P = tidy_res$p.value)
}) %>% mutate(Adjusted_LM_FDR = p.adjust(Adjusted_LM_P, method = "BH"))

write.csv(cov_adjusted_results, file.path(output_dir, "Table_6_ssGSEA_Covariate_Adjusted_Results.csv"), row.names = FALSE)

p_fig1g <- ggplot(score_df, aes(x = group, y = Burden_Score, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.65, color = NA, scale = "width", width = 0.85) +
  geom_boxplot(width = 0.16, fill = "white", color = "black", outlier.shape = NA, linewidth = 0.65) +
  scale_fill_manual(values = c("HC" = "#9ecae1", "PD" = "#fc9272"), guide = "none") +
  facet_wrap(~ Module, scales = "free_y", nrow = 1) + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.6, size = 4.8, fontface = "bold", color = "#b2182b") +
  labs(y = "Module Burden Score (ssGSEA)", x = NULL) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.7),
    strip.text = element_text(face = "bold", size = 8.8, color = "black"),
    axis.text.x = element_text(face = "bold", size = 11.5, color = "black"),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(face = "bold", size = 12, margin = ggplot2::margin(r = 10)),
    axis.line = element_line(color = "black", linewidth = 0.6),
    plot.margin = ggplot2::margin(10, 15, 10, 10),
    panel.spacing = unit(0.6, "lines")
  )

ggsave(file.path(output_dir, "Figure1G_Pathway_Activity.pdf"), p_fig1g, width = 14, height = 4.2, device = cairo_pdf)
ggsave(file.path(output_dir, "Figure1G_Pathway_Activity.png"), p_fig1g, width = 14, height = 4.2, dpi = 300)

cat("Phase 2 analysis complete. Figure 1 panels and tables saved to:", output_dir, "\n")