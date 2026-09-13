# ==============================================================================
# Plasma Proteomic Stratification & Clinical Decision Portal for Parkinson's Disease
# Production-Grade Public Application (Nature Aging 2026) - Complete 7-Scale Edition
# ==============================================================================

library(shiny)
library(bslib)
library(tidyverse)
library(randomForest)
library(DT)
library(plotly)
library(shinyWidgets)
library(readxl)

# --- 1. Load Pretrained Models & Global Metadata ---
models_dir <- "models"
data_dir   <- "data"

diag_data     <- readRDS(file.path(models_dir, "diag_offset_model.rds"))
subtype_data  <- readRDS(file.path(models_dir, "subtype_rf_classifier.rds"))
clinical_data <- readRDS(file.path(models_dir, "clinical_rf_predictors.rds"))
metadata      <- readRDS(file.path(models_dir, "app_metadata.rds"))

# 统一 UniProt 与 Symbol 映射目录
all_syms <- metadata$all_target_symbols
all_uids <- sapply(all_syms, function(s) {
  if (s %in% names(metadata$symbol_to_uniprot)) metadata$symbol_to_uniprot[[s]]
  else if (s %in% names(diag_data$scale_params$means)) s
  else "Unknown"
})

s1_official_symbols <- c("ITGA6", "ADAM10", "TSPAN14", "PDIA5", "SLC44A2", "CD151", "PLXNB2", "ITGA5", "EIF5B", "GNG5", "RHOF", "PIP4K2B", "TAOK3", "PURA", "ICAM2")
s2_official_symbols <- c("VTN", "KNG1", "SERPING1", "TMEM9", "APOA4", "PROS1", "PTK7", "PIGR", "AMBP", "MGP", "C4A", "C4BPA", "APOC1", "ILF2", "IGHV3-7")
s3_official_symbols <- c("SLC25A20", "MCCC2", "CLCN4", "IRAG2", "ECI1", "EDC4", "NDUFA9", "LONP1", "MYCBP2", "GLS", "PPP1R18", "APAF1", "CNOT1", "ATP5MK", "DLAT")

catalog_df_full <- data.frame(
  Gene_Symbol = all_syms,
  UniProt_Accession = unname(all_uids),
  Panel_Classification = case_when(
    all_syms == "ITGA6" ~ "Subtype 1 & Diagnostic (Dual Role)",
    all_syms %in% s1_official_symbols ~ "Subtype 1 Drivers",
    all_syms %in% s2_official_symbols ~ "Subtype 2 Drivers",
    all_syms %in% s3_official_symbols ~ "Subtype 3 Drivers",
    TRUE ~ "Diagnostic Panel (16-Protein)"
  ),
  Diagnostic_16_Panel = ifelse(all_syms %in% diag_data$panel_symbols, "Yes", "No"),
  Subtyping_45_Panel  = ifelse(all_syms %in% subtype_data$panel_symbols | all_syms %in% c("IGHV3-7", s1_official_symbols, s2_official_symbols, s3_official_symbols), "Yes", "No"),
  stringsAsFactors = FALSE
)

# --- 2. Bilingual Translation Dictionary ---
dict <- list(
  en = list(
    title = "PD Proteomic Decision Portal",
    tab_single = "Individual Patient Profiler",
    tab_batch = "Cohort Batch Profiling",
    tab_catalog = "Biomarker Directory",
    tab_about = "About & Methodology",
    
    # Tab Single (Individual)
    demo_header = "1. Demographics & Clinical Covariates",
    prot_header = "2. Plasma Proteomics Data Provision",
    mode_switch = "Choose How to Provide Proteomic Data:",
    mode_file = "Option A: Upload Full File or Paste Table",
    mode_form = "Option B: Manual Entry (Key Biomarkers)",
    upload_single_label = "Upload Patient Proteomic File (.csv, .xlsx, .tsv):",
    upload_single_help = "Supports full DIA proteomics matrices. Target markers will be automatically extracted.",
    paste_single_label = "Or Directly Paste Data from Excel (2 Columns: Protein, Value):",
    scale_switch = "Measurement Scale of Manual Inputs:",
    scale_raw = "Raw MS Log2 Intensity (Auto-Standardized)",
    scale_z = "Standardized Relative Level (Z-score)",
    btn_predict = "Run Comprehensive Patient Prediction",
    diag_card_title = "Diagnostic Risk Assessment (16-Protein Offset Model)",
    endotype_card_title = "Predicted Molecular Endotype (45-Biomarker Classifier)",
    scores_card_title = "Predicted Multi-Domain Clinical Trait Scores (Figure 5F, 10-Fold CV)",
    therapy_card_title = "Tailored Precision Medicine Roadmap (Figure 6D)",
    
    # Tab Batch (Cohort)
    upload_header = "Cohort Batch Data Input",
    upload_help = "Upload normalized plasma proteomic matrix (CSV). Gene Symbols and UniProt IDs are automatically recognized.",
    upload_btn = "Choose Cohort CSV File",
    demo_btn = "Load Benchmark Example Cohort",
    covar_help = "Default baseline values will be applied if clinical variables are missing.",
    download_btn = "Download Prediction Report (CSV)",
    results_header = "1. Predicted Diagnostics, Endotypes & Multi-Domain Scores (7 Scales)",
    pie_header = "2. Molecular Endotype Distribution",
    radar_header = "3. Multi-Domain Phenotypic Severity Radar (7 Dimensions)",
    
    disclaimer = "Disclaimer: This software is designed strictly for scientific research and translational exploration. It does not constitute standalone clinical diagnostic advice."
  ),
  zh = list(
    title = "帕金森病血浆蛋白质组学智能辅助决策系统",
    tab_single = "单患者个体精准决策",
    tab_batch = "队列级批量分型预测",
    tab_catalog = "标志物全景目录",
    tab_about = "关于平台与方法学",
    
    # Tab Single (Individual)
    demo_header = "1. 患者人口学与临床基线特征",
    prot_header = "2. 血浆蛋白质组学数据提供",
    mode_switch = "请选择蛋白质组数据提供方式：",
    mode_file = "方式一：上传完整组学文件 / 粘贴表格",
    mode_form = "方式二：手动表单录入 (核心标志物)",
    upload_single_label = "上传单患者完整组学文件 (.csv, .xlsx, .tsv)：",
    upload_single_help = "支持上传完整 DIA 质谱矩阵。系统会自动检索并提取目标标志物。",
    paste_single_label = "或直接从 Excel 复制两列数据粘贴于此 (蛋白名, 表达量)：",
    scale_switch = "手动录入的蛋白量纲类型：",
    scale_raw = "质谱原始 Log2 丰度 (系统自动执行基准归一化)",
    scale_z = "相对标准化分值 (Z-score)",
    btn_predict = "运行个体化精准分层分析",
    diag_card_title = "疾病诊断风险评估 (16 蛋白 Offset 模型)",
    endotype_card_title = "预测分子内型归属 (45 标志物随机森林分类器)",
    scores_card_title = "预测多维度临床症状评分 (完整 7 大量表，Figure 5F)",
    therapy_card_title = "匹配的精准医疗干预路线图 (Figure 6D)",
    
    # Tab Batch (Cohort)
    upload_header = "队列数据上传与配置",
    upload_help = "上传标准化血浆蛋白质组表达矩阵（CSV 格式）。系统自动识别 Gene Symbol 和 UniProt ID。",
    upload_btn = "选择队列 CSV 文件",
    demo_btn = "加载基准测试队列样例",
    covar_help = "未提供的临床变量将自动采用人群基线均值填补。",
    download_btn = "导出完整预测报告 (CSV)",
    results_header = "1. 预测诊断、分子内型与 7 大临床量表总览",
    pie_header = "2. 分子内型人群构成比例",
    radar_header = "3. 多维度临床表型指纹雷达图 (7 维度)",
    
    disclaimer = "免责声明：本系统仅供科学研究与学术探讨使用，不作为独立临床诊断或处方依据。"
  )
)

# --- 3. 内部主界面布局 ---
main_app_ui <- page_navbar(
  id = "main_navbar",
  title = textOutput("ui_app_title"),
  theme = bs_theme(version = 5, bootswatch = "cosmo", primary = "#1F78B4"),
  
  nav_spacer(),
  nav_item(
    div(style = "padding-top: 6px; margin-right: 15px;",
        radioGroupButtons(
          inputId = "app_lang", label = NULL,
          choices = c("English" = "en", "中文" = "zh"),
          selected = "en", size = "sm", status = "primary"
        )
    )
  ),
  
  # TAB 1: 单患者精准决策
  nav_panel(
    title = textOutput("ui_tab_single"),
    icon = icon("user-md"),
    layout_sidebar(
      sidebar = sidebar(
        width = 390,
        h6(textOutput("ui_demo_header"), class = "fw-bold text-primary"),
        fluidRow(
          column(6, numericInput("p_age", "Age / 年龄", value = 65, min = 30, max = 95)),
          column(6, selectInput("p_sex", "Sex / 性别", choices = c("Male / 男" = 1, "Female / 女" = 0), selected = 1))
        ),
        fluidRow(
          column(6, numericInput("p_bmi", "BMI (kg/m²)", value = 23.5, min = 15, max = 40)),
          column(6, numericInput("p_duration", "Duration / 病程(年)", value = 2.0, min = 0, max = 20))
        ),
        numericInput("p_ledd", "Levodopa Dose (LEDD, mg/day) / 用药", value = 0, min = 0, max = 2000),
        fluidRow(
          column(4, selectInput("p_hepatic", "Hepatic / 肝病", choices = c("No"=0, "Yes"=1), selected=0)),
          column(4, selectInput("p_kidney", "Renal / 肾病", choices = c("Normal"=0, "Abnormal"=1), selected=0)),
          column(4, selectInput("p_cv", "CV / 共病", choices = c("0"=0, "1"=1, "2+"=2), selected=0))
        ),
        hr(),
        h6(textOutput("ui_prot_header"), class = "fw-bold text-primary"),
        radioButtons("single_input_mode", textOutput("ui_mode_switch"),
                     choices = c("file_mode", "form_mode"),
                     choiceNames = list(textOutput("ui_mode_file"), textOutput("ui_mode_form")),
                     selected = "file_mode"),
        conditionalPanel(
          condition = "input.single_input_mode == 'file_mode'",
          fileInput("single_file_input", label = textOutput("ui_upload_single_label"), accept = c(".csv", ".tsv", ".xlsx")),
          p(textOutput("ui_upload_single_help"), class = "text-muted small mb-2"),
          textAreaInput("single_paste_area", label = textOutput("ui_paste_single_label"), 
                        placeholder = "GeneSymbol\tValue\nITGA6\t16.4\nVTN\t18.2\nNDUFA9\t14.8", rows = 4)
        ),
        conditionalPanel(
          condition = "input.single_input_mode == 'form_mode'",
          radioButtons("input_scale_type", textOutput("ui_scale_switch"),
                       choices = c("scale_z", "scale_raw"),
                       choiceNames = list(textOutput("ui_scale_z"), textOutput("ui_scale_raw")),
                       selected = "scale_z"),
          accordion(
            id = "accordion_full_biomarkers",
            accordion_panel(
              "16 Diagnostic Panel Markers",
              fluidRow(column(6, numericInput("in_CSRP1", "CSRP1", value = 0.0)), column(6, numericInput("in_ARL8B", "ARL8B", value = 0.0))),
              fluidRow(column(6, numericInput("in_RAP1A", "RAP1A", value = 0.0)), column(6, numericInput("in_LECT2", "LECT2", value = 0.0))),
              fluidRow(column(6, numericInput("in_SNRPD3", "SNRPD3", value = 0.0)), column(6, numericInput("in_STAG2", "STAG2", value = 0.0))),
              fluidRow(column(6, numericInput("in_BCAM", "BCAM", value = 0.0)), column(6, numericInput("in_VPS26C", "VPS26C", value = 0.0))),
              fluidRow(column(6, numericInput("in_GABARAPL2", "GABARAPL2", value = 0.0)), column(6, numericInput("in_SLFN14", "SLFN14", value = 0.0))),
              fluidRow(column(6, numericInput("in_CSNK2B", "CSNK2B", value = 0.0)), column(6, numericInput("in_PLEKHF1", "PLEKHF1", value = 0.0))),
              fluidRow(column(6, numericInput("in_SPTLC1", "SPTLC1", value = 0.0)), column(6, numericInput("in_CYBB", "CYBB", value = 0.0))),
              fluidRow(column(6, numericInput("in_DHX29", "DHX29", value = 0.0)))
            ),
            accordion_panel(
              "Subtype 1 Markers (15 Prots)",
              fluidRow(column(6, numericInput("in_ITGA6", "ITGA6", value = 0.0)), column(6, numericInput("in_ADAM10", "ADAM10", value = 0.0))),
              fluidRow(column(6, numericInput("in_TSPAN14", "TSPAN14", value = 0.0)), column(6, numericInput("in_PDIA5", "PDIA5", value = 0.0))),
              fluidRow(column(6, numericInput("in_SLC44A2", "SLC44A2", value = 0.0)), column(6, numericInput("in_CD151", "CD151", value = 0.0))),
              fluidRow(column(6, numericInput("in_PLXNB2", "PLXNB2", value = 0.0)), column(6, numericInput("in_EIF5B", "EIF5B", value = 0.0))),
              fluidRow(column(6, numericInput("in_ITGA5", "ITGA5", value = 0.0)), column(6, numericInput("in_GNG5", "GNG5", value = 0.0))),
              fluidRow(column(6, numericInput("in_RHOF", "RHOF", value = 0.0)), column(6, numericInput("in_PIP4K2B", "PIP4K2B", value = 0.0))),
              fluidRow(column(6, numericInput("in_TAOK3", "TAOK3", value = 0.0)), column(6, numericInput("in_PURA", "PURA", value = 0.0))),
              fluidRow(column(6, numericInput("in_ICAM2", "ICAM2", value = 0.0)))
            ),
            accordion_panel(
              "Subtype 2 Markers (15 Prots)",
              fluidRow(column(6, numericInput("in_VTN", "VTN", value = 0.0)), column(6, numericInput("in_KNG1", "KNG1", value = 0.0))),
              fluidRow(column(6, numericInput("in_SERPING1", "SERPING1", value = 0.0)), column(6, numericInput("in_TMEM9", "TMEM9", value = 0.0))),
              fluidRow(column(6, numericInput("in_APOA4", "APOA4", value = 0.0)), column(6, numericInput("in_PROS1", "PROS1", value = 0.0))),
              fluidRow(column(6, numericInput("in_PTK7", "PTK7", value = 0.0)), column(6, numericInput("in_PIGR", "PIGR", value = 0.0))),
              fluidRow(column(6, numericInput("in_AMBP", "AMBP", value = 0.0)), column(6, numericInput("in_MGP", "MGP", value = 0.0))),
              fluidRow(column(6, numericInput("in_C4A", "C4A", value = 0.0)), column(6, numericInput("in_C4BPA", "C4BPA", value = 0.0))),
              fluidRow(column(6, numericInput("in_APOC1", "APOC1", value = 0.0)), column(6, numericInput("in_ILF2", "ILF2", value = 0.0))),
              fluidRow(column(6, numericInput("in_IGHV37", "IGHV3-7", value = 0.0)))
            ),
            accordion_panel(
              "Subtype 3 Markers (15 Prots)",
              fluidRow(column(6, numericInput("in_SLC25A20", "SLC25A20", value = 0.0)), column(6, numericInput("in_MCCC2", "MCCC2", value = 0.0))),
              fluidRow(column(6, numericInput("in_CLCN4", "CLCN4", value = 0.0)), column(6, numericInput("in_IRAG2", "IRAG2", value = 0.0))),
              fluidRow(column(6, numericInput("in_ECI1", "ECI1", value = 0.0)), column(6, numericInput("in_EDC4", "EDC4", value = 0.0))),
              fluidRow(column(6, numericInput("in_NDUFA9", "NDUFA9", value = 0.0)), column(6, numericInput("in_LONP1", "LONP1", value = 0.0))),
              fluidRow(column(6, numericInput("in_MYCBP2", "MYCBP2", value = 0.0)), column(6, numericInput("in_GLS", "GLS", value = 0.0))),
              fluidRow(column(6, numericInput("in_PPP1R18", "PPP1R18", value = 0.0)), column(6, numericInput("in_APAF1", "APAF1", value = 0.0))),
              fluidRow(column(6, numericInput("in_CNOT1", "CNOT1", value = 0.0)), column(6, numericInput("in_ATP5MK", "ATP5MK", value = 0.0))),
              fluidRow(column(6, numericInput("in_DLAT", "DLAT", value = 0.0)))
            )
          )
        ),
        actionButton("btn_predict_single", label = textOutput("ui_btn_predict"), icon = icon("calculator"), class = "btn-primary w-100 mt-3")
      ),
      layout_columns(
        card(
          card_header(class = "bg-primary text-white", textOutput("ui_diag_card_title")),
          uiOutput("card_single_diag"),
          hr(),
          h6(textOutput("ui_endotype_card_title"), class = "fw-bold"),
          uiOutput("card_single_endotype"),
          hr(),
          h6(textOutput("ui_scores_card_title"), class = "fw-bold text-primary"),
          uiOutput("card_single_scales") # 💡 7 大量表卡片渲染区
        ),
        card(
          card_header(class = "bg-success text-white", textOutput("ui_therapy_card_title")),
          uiOutput("card_single_therapy")
        )
      )
    )
  ),
  
  # TAB 2: 批量队列预测
  nav_panel(
    title = textOutput("ui_tab_batch"),
    icon = icon("users"),
    layout_sidebar(
      sidebar = sidebar(
        width = 330,
        h5(textOutput("ui_upload_header"), class = "fw-bold text-primary"),
        p(textOutput("ui_upload_help"), class = "text-muted small"),
        fileInput("file_input", label = textOutput("ui_upload_btn"), accept = c(".csv")),
        actionButton("btn_load_demo", label = textOutput("ui_demo_btn"), icon = icon("flask"), class = "btn-outline-primary w-100 mb-2"),
        hr(),
        p(textOutput("ui_covar_help"), class = "text-muted small"),
        hr(),
        downloadButton("btn_download_csv", label = textOutput("ui_download_btn"), class = "btn-success w-100")
      ),
      card(
        card_header(class = "bg-primary text-white", textOutput("ui_results_header")),
        DTOutput("table_batch_results")
      ),
      layout_columns(
        card(
          card_header(textOutput("ui_pie_header")),
          plotlyOutput("plot_endotype_pie", height = "320px")
        ),
        card(
          card_header(textOutput("ui_radar_header")),
          plotlyOutput("plot_severity_radar", height = "320px")
        )
      )
    )
  ),
  
  # TAB 3: 标志物目录
  nav_panel(
    title = textOutput("ui_tab_catalog"),
    icon = icon("list-check"),
    card(
      card_header("Complete Biomarker Panel Directory & Reference Mapping (N = 60 Unique Proteins)"),
      DTOutput("table_marker_catalog")
    )
  ),
  
  # TAB 4: 关于
  nav_panel(
    title = textOutput("ui_tab_about"),
    icon = icon("circle-info"),
    card(
      card_header("About this Platform & Clinical Methodology"),
      uiOutput("about_content")
    )
  )
)

# --- 4. 外层包含密码门禁与【显示/隐藏密码】切换的 UI ---
ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "cosmo", primary = "#1F78B4"),
  
  # 门禁弹窗
  conditionalPanel(
    condition = "output.authenticated == false",
    div(style = "max-width: 440px; margin: 120px auto; padding: 35px; background: #ffffff; border-radius: 12px; box-shadow: 0 8px 25px rgba(0,0,0,0.08); border: 1px solid #e3e6f0;",
        h4("Peer-Review Access Portal", class = "text-primary fw-bold text-center mb-2"),
        h6("Nature Aging (2026) Dedicated Reviewer Channel", class = "text-muted text-center small mb-4"),
        p("This clinical decision-support portal is currently under active peer-review. Please enter the access passcode provided in the Response to Reviewers letter.", class = "text-muted small text-center mb-3"),
        passwordInput("passcode", "Enter Reviewer Passcode / 访问密码:", placeholder = "Passcode"),
        div(style = "margin-top: -10px; margin-bottom: 18px;",
            tags$input(type = "checkbox", id = "toggle_pw", 
                       onclick = "var x = document.getElementById('passcode'); x.type = (x.type === 'password') ? 'text' : 'password';"),
            tags$label(`for` = "toggle_pw", " 👁️ Show Passcode / 显示明文密码", style = "font-size: 0.85rem; color: #495057; cursor: pointer; user-select: none;")
        ),
        actionButton("btn_login", "Enter System / 解锁访问", class = "btn-primary w-100", icon = icon("key"))
    )
  ),
  
  conditionalPanel(
    condition = "output.authenticated == true",
    main_app_ui
  )
)

# --- 5. Server 计算引擎 ---
server <- function(input, output, session) {
  
  auth          <- reactiveVal(FALSE)
  dataset_input <- reactiveVal(NULL)
  
  observeEvent(input$btn_login, {
    if (trimws(input$passcode) == "NatureAging2026") {
      auth(TRUE)
    } else {
      showNotification("Incorrect Passcode. Please check the Response Letter.", type = "error")
    }
  })
  
  output$authenticated <- reactive({ auth() })
  outputOptions(output, "authenticated", suspendWhenHidden = FALSE)
  
  L <- reactive({
    lang <- input$app_lang
    if (is.null(lang) || !lang %in% c("en", "zh")) lang <- "en"
    dict[[lang]]
  })
  
  output$ui_app_title       <- renderText({ L()$title })
  output$ui_tab_single      <- renderText({ L()$tab_single })
  output$ui_tab_batch       <- renderText({ L()$tab_batch })
  output$ui_tab_catalog     <- renderText({ L()$tab_catalog })
  output$ui_tab_about       <- renderText({ L()$tab_about })
  
  output$ui_demo_header     <- renderText({ L()$demo_header })
  output$ui_prot_header     <- renderText({ L()$prot_header })
  output$ui_mode_switch     <- renderText({ L()$mode_switch })
  output$ui_mode_file       <- renderText({ L()$mode_file })
  output$ui_mode_form       <- renderText({ L()$mode_form })
  output$ui_upload_single_label <- renderText({ L()$upload_single_label })
  output$ui_upload_single_help  <- renderText({ L()$upload_single_help })
  output$ui_paste_single_label  <- renderText({ L()$paste_single_label })
  
  output$ui_scale_switch    <- renderText({ L()$scale_switch })
  output$ui_scale_raw       <- renderText({ L()$scale_raw })
  output$ui_scale_z         <- renderText({ L()$scale_z })
  output$ui_btn_predict     <- renderText({ L()$btn_predict })
  output$ui_diag_card_title <- renderText({ L()$diag_card_title })
  output$ui_endotype_card_title <- renderText({ L()$endotype_card_title })
  output$ui_scores_card_title   <- renderText({ L()$scores_card_title })
  output$ui_therapy_card_title  <- renderText({ L()$therapy_card_title })
  
  output$ui_upload_header   <- renderText({ L()$upload_header })
  output$ui_upload_help     <- renderText({ L()$upload_help })
  output$ui_upload_btn      <- renderText({ L()$upload_btn })
  output$ui_demo_btn        <- renderText({ L()$demo_btn })
  output$ui_covar_help      <- renderText({ L()$covar_help })
  output$ui_download_btn    <- renderText({ L()$download_btn })
  output$ui_results_header  <- renderText({ L()$results_header })
  output$ui_pie_header      <- renderText({ L()$pie_header })
  output$ui_radar_header    <- renderText({ L()$radar_header })
  
  # --- 1. 单患者全量预测推理 (集成完整 7 大量表) ---
  single_prediction <- reactive({
    req(input$btn_predict_single)
    
    age  <- input$p_age
    sex  <- as.numeric(input$p_sex)
    bmi  <- input$p_bmi
    dur  <- input$p_duration
    ledd <- input$p_ledd
    hep  <- as.numeric(input$p_hepatic)
    kid  <- as.numeric(input$p_kidney)
    cv   <- as.numeric(input$p_cv)
    
    prot_values <- list()
    
    if (input$single_input_mode == "file_mode") {
      parsed_df <- NULL
      if (!is.null(input$single_file_input)) {
        fpath <- input$single_file_input$datapath
        if (grepl("\\.xlsx$", fpath)) parsed_df <- read_excel(fpath)
        else parsed_df <- read.csv(fpath, stringsAsFactors = FALSE)
      } else if (!is.null(input$single_paste_area) && trimws(input$single_paste_area) != "") {
        txt <- trimws(input$single_paste_area)
        parsed_df <- tryCatch({
          read.table(text = txt, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
        }, error = function(e) {
          read.csv(text = txt, stringsAsFactors = FALSE)
        })
      }
      
      if (!is.null(parsed_df) && nrow(parsed_df) > 0) {
        if (ncol(parsed_df) >= 2 && nrow(parsed_df) > 5) {
          id_vec  <- as.character(parsed_df[[1]])
          val_vec <- as.numeric(parsed_df[[2]])
          for (k in seq_along(id_vec)) {
            prot_name <- id_vec[k]
            u <- if (prot_name %in% names(metadata$symbol_to_uniprot)) metadata$symbol_to_uniprot[[prot_name]] else prot_name
            prot_values[[u]] <- val_vec[k]
          }
        } else {
          for (cn in colnames(parsed_df)) {
            u <- if (cn %in% names(metadata$symbol_to_uniprot)) metadata$symbol_to_uniprot[[cn]] else cn
            prot_values[[u]] <- as.numeric(parsed_df[[cn]][1])
          }
        }
      }
    } else {
      is_raw <- input$input_scale_type == "scale_raw"
      read_prot <- function(sym, input_val) {
        val <- as.numeric(input_val)
        u <- metadata$symbol_to_uniprot[sym]
        if (is_raw && !is.na(u) && u %in% names(diag_data$scale_params$means)) {
          m <- diag_data$scale_params$means[u]
          s <- diag_data$scale_params$sds[u]
          if (!is.na(m) && !is.na(s) && s > 0) val <- (val - m) / s
        }
        if (!is.na(u)) prot_values[[u]] <<- val
      }
      
      # 16 Diagnostic Panel
      read_prot("CSRP1", input$in_CSRP1); read_prot("ARL8B", input$in_ARL8B)
      read_prot("RAP1A", input$in_RAP1A); read_prot("LECT2", input$in_LECT2)
      read_prot("SNRPD3", input$in_SNRPD3); read_prot("STAG2", input$in_STAG2)
      read_prot("BCAM", input$in_BCAM); read_prot("VPS26C", input$in_VPS26C)
      read_prot("GABARAPL2", input$in_GABARAPL2); read_prot("SLFN14", input$in_SLFN14)
      read_prot("CSNK2B", input$in_CSNK2B); read_prot("PLEKHF1", input$in_PLEKHF1)
      read_prot("SPTLC1", input$in_SPTLC1); read_prot("CYBB", input$in_CYBB); read_prot("DHX29", input$in_DHX29)
      
      # Subtype 1 (15 Proteins, including PIP4K2B)
      read_prot("ITGA6", input$in_ITGA6); read_prot("ADAM10", input$in_ADAM10); read_prot("TSPAN14", input$in_TSPAN14)
      read_prot("PDIA5", input$in_PDIA5); read_prot("SLC44A2", input$in_SLC44A2); read_prot("CD151", input$in_CD151)
      read_prot("PLXNB2", input$in_PLXNB2); read_prot("EIF5B", input$in_EIF5B); read_prot("ITGA5", input$in_ITGA5)
      read_prot("GNG5", input$in_GNG5); read_prot("RHOF", input$in_RHOF); read_prot("PIP4K2B", input$in_PIP4K2B)
      read_prot("TAOK3", input$in_TAOK3); read_prot("PURA", input$in_PURA); read_prot("ICAM2", input$in_ICAM2)
      
      # Subtype 2 (15 Proteins, including IGHV3-7)
      read_prot("VTN", input$in_VTN); read_prot("KNG1", input$in_KNG1); read_prot("SERPING1", input$in_SERPING1)
      read_prot("TMEM9", input$in_TMEM9); read_prot("APOA4", input$in_APOA4); read_prot("PROS1", input$in_PROS1)
      read_prot("PTK7", input$in_PTK7); read_prot("PIGR", input$in_PIGR); read_prot("AMBP", input$in_AMBP)
      read_prot("MGP", input$in_MGP); read_prot("C4A", input$in_C4A); read_prot("C4BPA", input$in_C4BPA)
      read_prot("APOC1", input$in_APOC1); read_prot("ILF2", input$in_ILF2); read_prot("A0A0B4J2B5", input$in_IGHV37)
      
      # Subtype 3 (15 Proteins, including DLAT)
      read_prot("SLC25A20", input$in_SLC25A20); read_prot("MCCC2", input$in_MCCC2); read_prot("CLCN4", input$in_CLCN4)
      read_prot("IRAG2", input$in_IRAG2); read_prot("ECI1", input$in_ECI1); read_prot("EDC4", input$in_EDC4)
      read_prot("NDUFA9", input$in_NDUFA9); read_prot("LONP1", input$in_LONP1); read_prot("MYCBP2", input$in_MYCBP2)
      read_prot("GLS", input$in_GLS); read_prot("PPP1R18", input$in_PPP1R18); read_prot("APAF1", input$in_APAF1)
      read_prot("CNOT1", input$in_CNOT1); read_prot("ATP5MK", input$in_ATP5MK); read_prot("DLAT", input$in_DLAT)
    }
    
    # 诊断模型计算
    clin_score <- predict(diag_data$clin_model, newdata = data.frame(
      age = age, sex = sex, BMI = bmi, hepatic_disease = hep, kidney.function = kid, CV_Metabolic_Cat = cv
    ), type = "link")
    
    diag_input <- data.frame(matrix(0, nrow=1, ncol=length(diag_data$panel_uniprots), dimnames=list(NULL, diag_data$panel_uniprots)))
    for (u in diag_data$panel_uniprots) {
      if (u %in% names(prot_values)) {
        val <- prot_values[[u]]
        if (input$single_input_mode == "file_mode" && !is.na(val) && val > 5) {
          m <- diag_data$scale_params$means[u]
          s <- diag_data$scale_params$sds[u]
          if (!is.na(m) && !is.na(s) && s > 0) val <- (val - m) / s
        }
        diag_input[[u]] <- val
      }
    }
    
    pd_prob <- as.numeric(predict(diag_data$diag_model, newdata = cbind(Clin_Score = clin_score, diag_input), type = "response"))
    
    # 亚型模型计算
    sub_input <- data.frame(matrix(0, nrow=1, ncol=length(subtype_data$panel_uniprots), dimnames=list(NULL, subtype_data$panel_uniprots)))
    for (u in subtype_data$panel_uniprots) {
      if (u %in% names(prot_values)) {
        sub_input[[u]] <- prot_values[[u]]
      }
    }
    
    sub_probs <- predict(subtype_data$classifier, newdata = sub_input, type = "prob")[1, ]
    assigned_st <- names(sub_probs)[which.max(sub_probs)]
    
    # 💡 7 大连续多维量表预测
    rf_single_feat <- data.frame(
      Age = age, Sex = sex, BMI = bmi, hepatic = hep, kidney = kid, CV_Met = cv, LEDD = ledd, Duration = dur, sub_input
    )
    colnames(rf_single_feat) <- make.names(colnames(rf_single_feat))
    
    pred_upd3  <- predict(clinical_data$models$UPDRS3_Score, newdata = rf_single_feat)
    pred_hy    <- predict(clinical_data$models$HY_Stage, newdata = rf_single_feat)
    pred_hamd  <- predict(clinical_data$models$HAMD_Score, newdata = rf_single_feat)
    pred_nmss  <- predict(clinical_data$models$NMSS_Score, newdata = rf_single_feat)
    pred_mmse  <- predict(clinical_data$models$MMSE_Score, newdata = rf_single_feat)
    pred_pdss  <- predict(clinical_data$models$PDSS_Score, newdata = rf_single_feat)
    pred_upsit <- predict(clinical_data$models$UPSIT_Score, newdata = rf_single_feat)
    
    list(pd_prob = pd_prob, sub_probs = sub_probs, assigned_st = assigned_st,
         upd3 = pred_upd3, hy = pred_hy, hamd = pred_hamd, nmss = pred_nmss,
         mmse = pred_mmse, pdss = pred_pdss, upsit = pred_upsit)
  })
  
  output$card_single_diag <- renderUI({
    res <- tryCatch(single_prediction(), error = function(e) NULL)
    is_zh <- input$app_lang == "zh"
    if (is.null(res)) {
      return(p(if(is_zh) "请配置左侧临床特征并提供蛋白质数据，点击下方按钮运行综合分析。" else "Configure parameters on the left and click 'Run Comprehensive Patient Prediction'.", class="text-muted"))
    }
    prob_val <- res$pd_prob * 100
    risk_class <- if(prob_val >= 50) "alert alert-danger" else "alert alert-success"
    call_text <- if(prob_val >= 50) {
      if(is_zh) sprintf("高风险：支持帕金森病诊断 (患病概率: %.1f%%)", prob_val) 
      else sprintf("High Risk: Consistent with Parkinson's Disease (Probability: %.1f%%)", prob_val)
    } else {
      if(is_zh) sprintf("低风险：倾向健康对照基线 (患病概率: %.1f%%)", prob_val) 
      else sprintf("Low Risk: Consistent with Control Baseline (Probability: %.1f%%)", prob_val)
    }
    div(class = risk_class, style="font-size:1.05rem; font-weight:bold;", call_text)
  })
  
  output$card_single_endotype <- renderUI({
    res <- tryCatch(single_prediction(), error = function(e) NULL)
    is_zh <- input$app_lang == "zh"
    if (is.null(res)) return(NULL)
    
    st <- res$assigned_st
    st_title <- case_when(
      st == "Subtype_1" ~ if(is_zh) "预测归属：Subtype 1 (整合素与自噬稳态轴)" else "Predicted Endotype: Subtype 1 (Integrin & Autophagy Axis)",
      st == "Subtype_2" ~ if(is_zh) "预测归属：Subtype 2 (神经炎症与补体免疫轴)" else "Predicted Endotype: Subtype 2 (Neuroinflammatory & Immune Axis)",
      TRUE ~ if(is_zh) "预测归属：Subtype 3 (线粒体代谢脆弱轴)" else "Predicted Endotype: Subtype 3 (Mitochondrial & Metabolic Axis)"
    )
    st_color <- case_when(st == "Subtype_1" ~ "alert-danger", st == "Subtype_2" ~ "alert-primary", TRUE ~ "alert-warning")
    
    tagList(
      div(class = paste("alert", st_color), style="font-size:1.05rem; font-weight:bold;", st_title),
      p(sprintf("Posterior Probabilities: Subtype 1 (%.1f%%) | Subtype 2 (%.1f%%) | Subtype 3 (%.1f%%)", 
                res$sub_probs["Subtype_1"]*100, res$sub_probs["Subtype_2"]*100, res$sub_probs["Subtype_3"]*100),
        class = "text-muted small")
    )
  })
  
  # 💡 7 大量表整齐卡片布局
  output$card_single_scales <- renderUI({
    res <- tryCatch(single_prediction(), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    is_zh <- input$app_lang == "zh"
    
    tagList(
      fluidRow(
        column(3, div(class="card p-2 text-center bg-light", 
                      div(class="text-muted small", if(is_zh) "运动 (UPDRS-III)" else "Motor (UPDRS-III)"),
                      h5(round(res$upd3, 1), class="fw-bold text-primary"))),
        column(3, div(class="card p-2 text-center bg-light", 
                      div(class="text-muted small", if(is_zh) "疾病分期 (H&Y)" else "Stage (H&Y)"),
                      h5(round(res$hy, 1), class="fw-bold text-primary"))),
        column(3, div(class="card p-2 text-center bg-light", 
                      div(class="text-muted small", if(is_zh) "抑郁 (HAMD)" else "Depression (HAMD)"),
                      h5(round(res$hamd, 1), class="fw-bold text-danger"))),
        column(3, div(class="card p-2 text-center bg-light", 
                      div(class="text-muted small", if(is_zh) "总非运动 (NMSS)" else "Non-motor (NMSS)"),
                      h5(round(res$nmss, 1), class="fw-bold text-danger")))
      ),
      div(style="height:8px;"),
      fluidRow(
        column(4, div(class="card p-2 text-center bg-light", 
                      div(class="text-muted small", if(is_zh) "认知功能 (MMSE)" else "Cognition (MMSE)"),
                      h5(round(res$mmse, 1), class="fw-bold text-success"))),
        column(4, div(class="card p-2 text-center bg-light", 
                      div(class="text-muted small", if(is_zh) "睡眠质量 (PDSS)" else "Sleep (PDSS)"),
                      h5(round(res$pdss, 1), class="fw-bold text-info"))),
        column(4, div(class="card p-2 text-center bg-light", 
                      div(class="text-muted small", if(is_zh) "嗅觉功能 (UPSIT)" else "Olfaction (UPSIT)"),
                      h5(round(res$upsit, 1), class="fw-bold text-warning")))
      )
    )
  })
  
  output$card_single_therapy <- renderUI({
    res <- tryCatch(single_prediction(), error = function(e) NULL)
    is_zh <- input$app_lang == "zh"
    if (is.null(res)) {
      return(p(if(is_zh) "分析完成后将自动生成匹配的精准医疗干预路线图。" else "Stratified therapeutic roadmap will appear here after profiling.", class="text-muted"))
    }
    
    st <- res$assigned_st
    if (st == "Subtype_1") {
      tagList(
        h6(if(is_zh) "Subtype 1 推荐精准干预路线图：" else "Subtype 1 Stratified Clinical Roadmap:", class="fw-bold text-success"),
        tags$ul(
          tags$li(strong(if(is_zh) "标准护理：" else "Standard Care: "), if(is_zh) "常规多巴胺替代治疗（多巴胺能反应良好）。" else "Standard dopaminergic replacement therapy (favorable response)."),
          tags$li(strong(if(is_zh) "自噬稳态增强：" else "Autophagy Enhancement: "), if(is_zh) "优先考虑自噬-溶酶体功能激活剂与蛋白稳态维持疗法。" else "Prioritize autophagy-lysosomal enhancers and proteostasis stabilizers."),
          tags$li(strong(if(is_zh) "基质黏附保护：" else "Structural Integrity: "), if(is_zh) "维持整合素-局域黏附受体信号稳态与突触囊泡转运。" else "Support integrin-mediated focal adhesion signaling and synaptic vesicle trafficking.")
        )
      )
    } else if (st == "Subtype_2") {
      tagList(
        h6(if(is_zh) "Subtype 2 推荐精准干预路线图：" else "Subtype 2 Stratified Clinical Roadmap:", class="fw-bold text-primary"),
        tags$ul(
          tags$li(strong(if(is_zh) "标准护理：" else "Standard Care: "), if(is_zh) "优化多巴胺能药物，密切监测严重运动并发症与运动波动。" else "Optimize dopaminergic regimen and monitor closely for motor fluctuations."),
          tags$li(strong(if(is_zh) "靶向免疫与补体抑制：" else "Targeted Immunomodulation: "), if(is_zh) "优先考虑抗神经炎症干预、补体级联反应抑制剂（如 C4/C1-INH 靶向）。" else "Prioritize anti-neuroinflammatory strategies and complement cascade inhibitors."),
          tags$li(strong(if(is_zh) "血脑屏障保护：" else "Vascular/BBB Protection: "), if(is_zh) "内皮功能稳定剂与降低血脑屏障高通透性疗法。" else "Endothelial stabilizing agents to mitigate blood-brain barrier hyperpermeability."),
          tags$li(strong(if(is_zh) "精神心理预警：" else "Psychiatric Monitoring: "), if(is_zh) "早期开展抑郁（HAMD）筛查与预防性干预。" else "Early affective screening and prophylactic intervention for depressive burden.")
        )
      )
    } else {
      tagList(
        h6(if(is_zh) "Subtype 3 推荐精准干预路线图：" else "Subtype 3 Stratified Clinical Roadmap:", class="fw-bold text-warning"),
        tags$ul(
          tags$li(strong(if(is_zh) "标准护理：" else "Standard Care: "), if(is_zh) "标准多巴胺能治疗。" else "Standard dopaminergic therapy."),
          tags$li(strong(if(is_zh) "线粒体能量复苏：" else "Mitochondrial Rescue: "), if(is_zh) "补充线粒体复合物 I 辅因子（如 CoQ10、NAD+ 前体、核黄素）以改善 ATP 合成受损。" else "Mitochondrial Complex I cofactors (e.g., CoQ10, NAD+ precursors) to support bioenergetics."),
          tags$li(strong(if(is_zh) "抗氧化与代谢重塑：" else "Metabolic Reprogramming: "), if(is_zh) "线粒体靶向抗氧化剂、支链氨基酸（BCAA）与脂代谢代偿营养支持。" else "Mitochondrial-targeted antioxidants and lipid/BCAA metabolic support."),
          tags$li(strong(if(is_zh) "感觉障碍管理：" else "Sensory Care: "), if(is_zh) "重度嗅觉丧失（UPSIT）的生活安全预警与嗅觉康复训练。" else "Comprehensive olfactory rehabilitation and environmental safety guidance.")
        )
      )
    }
  })
  
  # --- 2. 批量队列计算 (包含 7 大量表预测) ---
  observeEvent(input$btn_load_demo, {
    demo_file <- file.path(data_dir, "example_test_samples.csv")
    if (file.exists(demo_file)) {
      dataset_input(read.csv(demo_file, stringsAsFactors = FALSE))
      showNotification(
        if(input$app_lang=="zh") "已成功载入基准测试队列样本！" else "Example benchmark cohort loaded successfully!", 
        type = "message"
      )
    }
  })
  
  observeEvent(input$file_input, {
    req(input$file_input)
    df <- read.csv(input$file_input$datapath, stringsAsFactors = FALSE)
    dataset_input(df)
  })
  
  batch_predictions <- reactive({
    req(dataset_input())
    raw_df <- dataset_input()
    
    col_names <- colnames(raw_df)
    for (i in seq_along(col_names)) {
      cn <- col_names[i]
      if (cn %in% names(metadata$symbol_to_uniprot)) {
        colnames(raw_df)[i] <- metadata$symbol_to_uniprot[[cn]]
      }
    }
    
    n_samps <- nrow(raw_df)
    age_vec  <- if ("Age" %in% colnames(raw_df)) as.numeric(raw_df$Age) else rep(65, n_samps)
    sex_vec  <- if ("Sex" %in% colnames(raw_df)) as.numeric(raw_df$Sex) else rep(1, n_samps)
    bmi_vec  <- if ("BMI" %in% colnames(raw_df)) as.numeric(raw_df$BMI) else rep(23.5, n_samps)
    hep_vec  <- if ("hepatic_disease" %in% colnames(raw_df)) as.numeric(raw_df$hepatic_disease) else rep(0, n_samps)
    kid_vec  <- if ("kidney.function" %in% colnames(raw_df)) as.numeric(raw_df$kidney.function) else rep(0, n_samps)
    cv_vec   <- if ("CV_Metabolic_Cat" %in% colnames(raw_df)) as.numeric(raw_df$CV_Metabolic_Cat) else rep(0, n_samps)
    ledd_vec <- if ("LEDD" %in% colnames(raw_df)) as.numeric(raw_df$LEDD) else rep(0, n_samps)
    dur_vec  <- if ("Duration" %in% colnames(raw_df)) as.numeric(raw_df$Duration) else rep(2.0, n_samps)
    
    diag_uniprots <- diag_data$panel_uniprots
    diag_mat <- data.frame(matrix(0, nrow = n_samps, ncol = length(diag_uniprots), dimnames = list(NULL, diag_uniprots)))
    for (u in diag_uniprots) {
      if (u %in% colnames(raw_df)) {
        raw_v <- as.numeric(raw_df[[u]])
        m <- diag_data$scale_params$means[u]
        s <- diag_data$scale_params$sds[u]
        diag_mat[[u]] <- if(!is.na(m) && !is.na(s) && s > 0) (raw_v - m) / s else raw_v
      }
    }
    
    clin_score_vec <- predict(diag_data$clin_model, newdata = data.frame(
      age = age_vec, sex = sex_vec, BMI = bmi_vec, hepatic_disease = hep_vec, kidney.function = kid_vec, CV_Metabolic_Cat = cv_vec
    ), type = "link")
    
    diag_probs <- predict(diag_data$diag_model, newdata = cbind(Clin_Score = clin_score_vec, diag_mat), type = "response")
    
    sub_uniprots <- subtype_data$panel_uniprots
    sub_mat <- data.frame(matrix(0, nrow = n_samps, ncol = length(sub_uniprots), dimnames = list(NULL, sub_uniprots)))
    for (u in sub_uniprots) {
      if (u %in% colnames(raw_df)) {
        sub_mat[[u]] <- as.numeric(raw_df[[u]])
      }
    }
    
    sub_probs <- predict(subtype_data$classifier, newdata = sub_mat, type = "prob")
    sub_pred  <- colnames(sub_probs)[max.col(sub_probs)]
    
    rf_feat_df <- data.frame(
      Age = age_vec, Sex = sex_vec, BMI = bmi_vec, hepatic = hep_vec, kidney = kid_vec, CV_Met = cv_vec, LEDD = ledd_vec, Duration = dur_vec, sub_mat
    )
    colnames(rf_feat_df) <- make.names(colnames(rf_feat_df))
    
    # 💡 7 大多维量表批量预测
    pred_upd3  <- predict(clinical_data$models$UPDRS3_Score, newdata = rf_feat_df)
    pred_hy    <- predict(clinical_data$models$HY_Stage, newdata = rf_feat_df)
    pred_hamd  <- predict(clinical_data$models$HAMD_Score, newdata = rf_feat_df)
    pred_nmss  <- predict(clinical_data$models$NMSS_Score, newdata = rf_feat_df)
    pred_mmse  <- predict(clinical_data$models$MMSE_Score, newdata = rf_feat_df)
    pred_pdss  <- predict(clinical_data$models$PDSS_Score, newdata = rf_feat_df)
    pred_upsit <- predict(clinical_data$models$UPSIT_Score, newdata = rf_feat_df)
    
    data.frame(
      Sample_ID = if("sample" %in% colnames(raw_df)) raw_df$sample else paste0("Sample_", 1:n_samps),
      PD_Risk_Probability = sprintf("%.1f%%", diag_probs * 100),
      Diagnostic_Call = ifelse(diag_probs >= 0.5, "PD (High Risk)", "HC (Low Risk)"),
      Assigned_Endotype = sub_pred,
      Prob_Subtype1 = sprintf("%.1f%%", sub_probs[, "Subtype_1"] * 100),
      Prob_Subtype2 = sprintf("%.1f%%", sub_probs[, "Subtype_2"] * 100),
      Prob_Subtype3 = sprintf("%.1f%%", sub_probs[, "Subtype_3"] * 100),
      Pred_Motor_UPDRS3   = round(pred_upd3, 1),
      Pred_Stage_HY       = round(pred_hy, 1),
      Pred_Depression_HAMD= round(pred_hamd, 1),
      Pred_NonMotor_NMSS  = round(pred_nmss, 1),
      Pred_Cognition_MMSE = round(pred_mmse, 1),
      Pred_Sleep_PDSS     = round(pred_pdss, 1),
      Pred_Olfaction_UPSIT= round(pred_upsit, 1)
    )
  })
  
  output$table_batch_results <- renderDT({
    req(batch_predictions())
    datatable(batch_predictions(), options = list(pageLength = 6, scrollX = TRUE), rownames = FALSE) %>%
      formatStyle('Diagnostic_Call', color = styleEqual(c("PD (High Risk)", "HC (Low Risk)"), c("#d73027", "#1a9850")), fontWeight = 'bold') %>%
      formatStyle('Assigned_Endotype', backgroundColor = styleEqual(c("Subtype_1", "Subtype_2", "Subtype_3"), c("#ffcccc", "#cce5ff", "#fff2cc")), fontWeight = 'bold')
  })
  
  output$plot_endotype_pie <- renderPlotly({
    req(batch_predictions())
    df <- batch_predictions()
    pie_data <- as.data.frame(table(df$Assigned_Endotype))
    colnames(pie_data) <- c("Endotype", "Count")
    
    plot_ly(pie_data, labels = ~Endotype, values = ~Count, type = 'pie',
            marker = list(colors = c("#d73027", "#4575b4", "#fdae61")),
            textinfo = "label+percent", hole = 0.4) %>%
      layout(showlegend = TRUE, margin = list(t=10, b=10, l=10, r=10))
  })
  
  output$plot_severity_radar <- renderPlotly({
    req(batch_predictions())
    df <- batch_predictions()
    
    m_upd3  <- mean(df$Pred_Motor_UPDRS3, na.rm = TRUE)
    m_hy    <- mean(df$Pred_Stage_HY, na.rm = TRUE) * 10
    m_hamd  <- mean(df$Pred_Depression_HAMD, na.rm = TRUE) * 2
    m_nmss  <- mean(df$Pred_NonMotor_NMSS, na.rm = TRUE) / 2
    m_mmse  <- (30 - mean(df$Pred_Cognition_MMSE, na.rm = TRUE)) * 5
    m_pdss  <- (150 - mean(df$Pred_Sleep_PDSS, na.rm = TRUE)) / 2
    m_upsit <- (40 - mean(df$Pred_Olfaction_UPSIT, na.rm = TRUE))
    
    r_vals <- c(m_upd3, m_hy, m_hamd, m_nmss, m_mmse, m_pdss, m_upsit, m_upd3)
    t_labels <- c('Motor (UPDRS-III)', 'Stage (H&Y)', 'Depression (HAMD)', 'Non-Motor (NMSS)', 
                  'Cognitive Decline', 'Sleep Deficit', 'Olfactory Loss', 'Motor (UPDRS-III)')
    
    plot_ly(type = 'scatterpolar', fill = 'toself',
            r = r_vals, theta = t_labels,
            line = list(color = '#1F78B4')) %>%
      layout(polar = list(radialaxis = list(visible = TRUE, range = c(0, max(50, max(r_vals)*1.15)))),
             margin = list(t=25, b=25, l=25, r=25))
  })
  
  output$btn_download_csv <- downloadHandler(
    filename = function() { paste0("PD_Proteomics_Stratification_Report_7Scales_", Sys.Date(), ".csv") },
    content = function(file) { write.csv(batch_predictions(), file, row.names = FALSE) }
  )
  
  output$table_marker_catalog <- renderDT({
    datatable(catalog_df_full, options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE) %>%
      formatStyle('Panel_Classification', 
                  color = styleEqual(
                    c("Subtype 1 Drivers", "Subtype 2 Drivers", "Subtype 3 Drivers", "Diagnostic Panel (16-Protein)", "Subtype 1 & Diagnostic (Dual Role)"),
                    c("#d73027", "#1F78B4", "#d95f02", "#756bb1", "#6a3d9a")
                  ), fontWeight = 'bold')
  })
  
  output$about_content <- renderUI({
    is_zh <- input$app_lang == "zh"
    if (is_zh) {
      tagList(
        p("本平台基于高通量 Orbitrap Astral DIA 质谱技术对大型双中心东亚队列（N = 1,119）进行系统研究所构建。"),
        p("整合了校正 6 协变量的 16 蛋白临床 Offset 诊断模型，以及校正 8 协变量的 45 标志物分子内型分类器与 7 大临床量表预测引擎（Figure 5F）。"),
        p(strong("参考文献："), "Large-scale plasma proteomics identifies molecularly distinct PD endotypes associated with motor and non-motor phenotypes. Nature Aging, 2026."),
        p(class="text-muted small", L()$disclaimer)
      )
    } else {
      tagList(
        p("This interactive clinical decision-support portal was established using high-throughput Orbitrap Astral DIA mass spectrometry across a dual-center East Asian cohort (N = 1,119)."),
        p("It integrates a 16-protein clinical-offset diagnostic model alongside an 8-covariate residualized 45-biomarker molecular endotyping classifier and 7 continuous clinical trait predictors (Figure 5F)."),
        p(strong("Citation: "), "Large-scale plasma proteomics identifies molecularly distinct PD endotypes associated with motor and non-motor phenotypes. Nature Aging, 2026."),
        p(class="text-muted small", L()$disclaimer)
      )
    }
  })
}

# 启动应用
shinyApp(ui = ui, server = server)