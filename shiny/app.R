#
# Synergy Analysis Shiny App
# Single-file app for transcriptomics synergy analysis
#
# Launch: shiny::runApp("synergy_framework/shiny/")
#

library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)
library(tidyr)
library(readr)

# Increase file upload limit (default 5MB, NHF files are ~65MB each)
options(shiny.maxRequestSize = 500 * 1024^2)

# Robust source paths — resolve framework root relative to this app file
app_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
  error = function(e) getwd()
)
framework_root <- dirname(app_dir)
source(file.path(framework_root, "R", "synergy_io.R"))
source(file.path(framework_root, "R", "synergy_core.R"))
source(file.path(framework_root, "R", "synergy_plot.R"))
source(file.path(framework_root, "R", "synergy_report.R"))

# ── Color palette ─────────────────────────────────────────────────────────────
COLORS <- list(
  primary   = "#5A6EA8",
  secondary = "#6B8EAD",
  up        = "#D64545",
  down      = "#4A7FB5",
  synergy   = "#E8923F",
  nt        = "#6B8EAD",
  a         = "#78A55A",
  b         = "#D64545",
  c         = "#8A62A8",
  bg        = "#F7F8FA",
  card_bg   = "#FFFFFF"
)

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- page_sidebar(
  title = tags$span(
    tags$img(src = "data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 24 24' fill='none' stroke='%235A6EA8' stroke-width='2'%3E%3Ccircle cx='12' cy='12' r='3'/%3E%3Ccircle cx='12' cy='5' r='2'/%3E%3Ccircle cx='12' cy='19' r='2'/%3E%3Ccircle cx='5' cy='8' r='2'/%3E%3Ccircle cx='19' cy='8' r='2'/%3E%3Ccircle cx='5' cy='16' r='2'/%3E%3Ccircle cx='19' cy='16' r='2'/%3E%3Cline x1='12' y1='8' x2='12' y2='9'/%3E%3Cline x1='7.05' y1='9.05' x2='8.05' y2='9.95'/%3E%3Cline x1='16.95' y1='9.95' x2='15.95' y2='9.05'/%3E%3Cline x1='6.05' y1='15.05' x2='7.15' y2='14.15'/%3E%3Cline x1='17.95' y1='14.15' x2='16.85' y2='15.05'/%3E%3C/svg%3E",
             height = "24px", style = "margin-right:8px; vertical-align:middle;"),
    "Synergy Analysis"
  ),
  theme = bs_theme(
    version      = 5,
    primary      = COLORS$primary,
    base_font    = bslib::font_google("Inter", wght = c(300, 400, 500, 600)),
    heading_font = bslib::font_google("Inter", wght = c(500, 600, 700)),
    font_scale   = 0.92,
    bg           = COLORS$bg,
    fg           = "#2C3E50",
    "border-color" = "#E5E8EC"
  ),
  tags$head(tags$style(HTML(sprintf("
    /* Typography */
    body { font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif; }
    .navbar-brand { font-weight: 700; letter-spacing: -0.01em; color: %s !important; }

    /* Cards */
    .card {
      border: 1px solid #E5E8EC;
      border-radius: 12px;
      box-shadow: 0 1px 4px rgba(0,0,0,0.04);
      transition: box-shadow 0.2s;
    }
    .card:hover { box-shadow: 0 2px 8px rgba(0,0,0,0.06); }
    .card-header {
      background: #FAFBFC;
      border-bottom: 1px solid #EEF0F2;
      font-weight: 600;
      font-size: 0.9rem;
      color: #3A4A5C;
      padding: 0.75rem 1rem;
    }

    /* Form controls */
    .form-label { font-weight: 500; font-size: 0.8rem; color: #556677; margin-bottom: 0.15rem; }
    .form-control, .form-select { border-color: #D5DAE0; border-radius: 6px; font-size: 0.85rem; }
    .form-control:focus, .form-select:focus { border-color: %s; box-shadow: 0 0 0 0.15rem rgba(90,110,168,0.15); }

    /* Sidebar */
    .sidebar { background: #FAFBFC; border-right: 1px solid #E5E8EC; }
    .sidebar-title { font-weight: 700; font-size: 0.85rem; color: #3A4A5C; text-transform: uppercase; letter-spacing: 0.04em; }

    /* Value boxes */
    .bslib-value-box { border-radius: 10px; }

    /* Buttons */
    .btn { border-radius: 6px; font-weight: 500; font-size: 0.85rem; transition: all 0.15s; }
    .btn-primary { background: %s; border-color: %s; }
    .btn-primary:hover { background: #4A5E98; border-color: #42568C; }

    /* Tables */
    .dataTables_wrapper { font-size: 0.82rem; }
    table.dataTable thead th { font-weight: 600; color: #3A4A5C; background: #FAFBFC; }

    /* Radio buttons */
    .radio-inline, .checkbox-inline { margin-right: 1rem; font-weight: 500; }

    /* Section divider */
    .section-divider { border-top: 1px solid #EEF0F2; margin: 1rem 0; }

    /* Status badge */
    .status-ready { color: #3A8C5C; font-weight: 600; }
    .status-waiting { color: #8899AA; }

    /* Input mode tabs */
    .input-mode-active { font-weight: 600; }
  ", COLORS$primary, COLORS$primary, COLORS$primary, COLORS$primary)))
  ),

  sidebar = sidebar(
    width = 290,
    padding = "0.6rem",
    tags$div(class = "sidebar-title", "Parameters"),

    accordion(
      open = c("params", "run"),
      multiple = TRUE,

      accordion_panel(
        "Synergy Filter",
        value = "params",
        icon = icon("sliders", lib = "font-awesome"),

        radioButtons("mode", "Synergy mode",
                     choices = c("Strict (4 criteria)" = "strict",
                                 "Relaxed (C vs NT only)" = "relaxed"),
                     selected = "strict"),
        tags$div(style = "font-size: 0.72rem; color: #6B7A8C; margin-top:-0.25rem; margin-bottom: 0.6rem;",
                 "Strict: all 3 comparisons sig + magnitude.",
                 tags$br(),
                 "Relaxed: only C vs NT sig + magnitude."),
        selectInput("pval_col", "Significance column",
                    choices = c("Adjusted P (Qvalue)" = "Qvalue",
                                "Raw P-value"          = "Pvalue")),
        sliderInput("p_cutoff", "P-value threshold",
                    min = 0.0001, max = 0.2, value = 0.05, step = 0.001),
        sliderInput("fc_cutoff", "|log2FC| minimum",
                    min = 0, max = 3, value = 0, step = 0.1)
      ),

      accordion_panel(
        "Group Labels",
        value = "labels",
        icon = icon("tags", lib = "font-awesome"),
        tags$div(
          class = "row g-1",
          tags$div(class = "col-6", textInput("label_nt", "NT", "NT", width = "100%")),
          tags$div(class = "col-6", textInput("label_c",  "A+B", "A+B", width = "100%")),
          tags$div(class = "col-6", textInput("label_a",  "A",   "A",   width = "100%")),
          tags$div(class = "col-6", textInput("label_b",  "B",   "B",   width = "100%"))
        )
      )
    ),

    tags$div(style = "height: 0.4rem;"),
    actionButton("run_btn", "Run Analysis",
                 class = "btn-primary w-100",
                 icon = icon("play-circle", lib = "font-awesome")),

    conditionalPanel(
      "output.analysis_done == true",
      tags$div(style = "height: 0.6rem;"),
      tags$div(class = "sidebar-title", "Export"),
      downloadButton("download_excel", "Excel  (incl. FPKM)",
                     class = "btn-outline-primary w-100 mb-2"),
      downloadButton("download_csv",   "Full CSV",
                     class = "btn-outline-secondary w-100 mb-2"),
      downloadButton("download_report", "HTML Report",
                     class = "btn-outline-secondary w-100")
    )
  ),

  navset_card_tab(
    # ── Tab 1: Data ────────────────────────────────────────────────────────
    nav_panel(
      "Data Input",
      icon = icon("folder-open", lib = "font-awesome"),

      layout_columns(
        col_widths = c(7, 5),
        card(
          card_header("Comparison Files"),
          tags$div(
            style = "padding: 0.5rem; background: #F0F4FF; border-radius: 8px; margin-bottom: 1rem;",
            tags$span(icon("info-circle", lib = "font-awesome"), style = "color: #5A6EA8;"),
            " Upload 5 pairwise edgeR/DESeq2 result files, or provide local file paths."
          ),
          radioButtons("input_mode", NULL,
                       choices = c("Upload from browser" = "upload",
                                   "Read from local disk" = "local"),
                       inline = TRUE),
          # ── Upload mode ──────────────────────────────────────────────
          conditionalPanel(
            "input.input_mode == 'upload'",
            fileInput("file_c_vs_nt", "C vs NT",  accept = c(".xls", ".tsv", ".txt", ".csv")),
            fileInput("file_a_vs_nt", "A vs NT",  accept = c(".xls", ".tsv", ".txt", ".csv")),
            fileInput("file_b_vs_nt", "B vs NT",  accept = c(".xls", ".tsv", ".txt", ".csv")),
            fileInput("file_c_vs_a",  "C vs A",   accept = c(".xls", ".tsv", ".txt", ".csv")),
            fileInput("file_c_vs_b",  "C vs B",   accept = c(".xls", ".tsv", ".txt", ".csv"))
          ),
          # ── Local path mode ──────────────────────────────────────────
          conditionalPanel(
            "input.input_mode == 'local'",
            textInput("path_c_vs_nt", "C vs NT",  placeholder = "/full/path/to/C_vs_NT.xls"),
            textInput("path_a_vs_nt", "A vs NT",  placeholder = "/full/path/to/A_vs_NT.xls"),
            textInput("path_b_vs_nt", "B vs NT",  placeholder = "/full/path/to/B_vs_NT.xls"),
            textInput("path_c_vs_a",  "C vs A",   placeholder = "/full/path/to/C_vs_A.xls"),
            textInput("path_c_vs_b",  "C vs B",   placeholder = "/full/path/to/C_vs_B.xls")
          )
        ),
        card(
          card_header("Data Status"),
          uiOutput("data_status_ui"),
          hr(),
          selectInput("preview_file", "Quick preview",
                      choices = c("C vs NT", "A vs NT", "B vs NT", "C vs A", "C vs B")),
          DTOutput("preview_table", height = "320px")
        )
      ),
      card(
        card_header("Column Name Mapping"),
        tags$p(class = "text-muted", style = "font-size:0.82rem;",
               "If your file uses different column names, adjust them here. ",
               "Changes apply when you re-load data or re-run analysis."),
        tags$div(
          class = "row g-2 mb-2",
          tags$div(class = "col-3", textInput("col_gene_id",   "Gene ID",     "gene_id",   width = "100%")),
          tags$div(class = "col-3", textInput("col_gene_name", "Gene name",   "gene_name", width = "100%")),
          tags$div(class = "col-3", textInput("col_log2fc",    "log2FC",       "log2FC",    width = "100%")),
          tags$div(class = "col-3", textInput("col_pvalue",    "P-value",      "Pvalue",    width = "100%"))
        ),
        tags$div(
          class = "row g-2",
          tags$div(class = "col-3", textInput("col_qvalue",    "Q-value",      "Qvalue",    width = "100%")),
          tags$div(class = "col-3", textInput("col_updown",    "Direction",    "updown",    width = "100%")),
          tags$div(class = "col-3", textInput("col_up_value",  "UP value",     "UP",        width = "100%")),
          tags$div(class = "col-3", textInput("col_down_value","DOWN value",   "DOWN",      width = "100%"))
        )
      )
    ),

    # ── Tab 2: Results ────────────────────────────────────────────────────
    nav_panel(
      "Results",
      icon = icon("table", lib = "font-awesome"),

      layout_columns(
        col_widths = c(4, 4, 4),
        value_box(
          title = "Genes tested",
          value = textOutput("total_genes"),
          showcase = icon("dna", lib = "font-awesome"),
          theme = "primary",
          style = "min-height: 110px;"
        ),
        value_box(
          title = "Synergy UP",
          value = textOutput("n_synergy_up"),
          showcase = icon("arrow-up", lib = "font-awesome"),
          theme = "danger",
          style = "min-height: 110px;"
        ),
        value_box(
          title = "Synergy DOWN",
          value = textOutput("n_synergy_down"),
          showcase = icon("arrow-down", lib = "font-awesome"),
          theme = "secondary",
          style = "min-height: 110px;"
        )
      ),
      card(
        card_header(
          tags$span("Synergistic Genes"),
          tags$span(style = "float:right;",
                    radioButtons("synergy_direction", NULL,
                                 choices = c("UP", "DOWN"), inline = TRUE))
        ),
        DTOutput("synergy_table", height = "480px")
      )
    ),

    # ── Tab 3: Visualization ──────────────────────────────────────────────
    nav_panel(
      "Visualization",
      icon = icon("chart-bar", lib = "font-awesome"),

      layout_columns(
        col_widths = c(7, 5),
        # Left column: Volcano on top, Effector Contribution below
        tags$div(
          card(
            card_header("Volcano Plot  (C vs NT, synergy genes highlighted)"),
            plotlyOutput("volcano_plot", height = "380px")
          ),
          tags$div(style = "height: 0.75rem;"),
          card(
            card_header(
              tags$span("Effect Contribution"),
              tags$span(style = "float:right;",
                        radioButtons("contrib_dir", NULL,
                                     choices = c("UP", "DOWN"), inline = TRUE))
            ),
            plotOutput("contrib_plot", height = "380px")
          )
        ),
        # Right column: Gene Inspector
        card(
          fillable = FALSE,
          full_screen = TRUE,
          card_header("Gene Inspector"),
          tags$div(style = "padding: 0.25rem 0.5rem;",
            selectizeInput("detail_gene",
                           "Select a synergistic gene:",
                           choices = NULL, width = "100%")
          ),
          tags$div(
            style = "padding: 0 0.5rem;",
            plotOutput("gene_detail_plot", height = "520px")
          ),
          tags$div(style = "padding: 0 0.5rem 0.5rem 0.5rem;",
            tableOutput("gene_detail_table")
          )
        )
      )
    ),

    # ── Tab 4: Enrichment Analysis ────────────────────────────────────────
    nav_panel(
      "Enrichment",
      icon = icon("diagram-project", lib = "font-awesome"),

      card(
        card_header("ORA Settings"),
        tags$div(
          class = "row g-2",
          tags$div(class = "col-2",
                   radioButtons("ora_dir", "Gene set",
                                choices = c("Synergy UP" = "up",
                                            "Synergy DOWN" = "down"),
                                inline = TRUE)),
          tags$div(class = "col-2",
                   radioButtons("ora_db", "Database",
                                choices = c("GO BP" = "gobp", "KEGG" = "kegg"),
                                inline = TRUE)),
          tags$div(class = "col-2",
                   selectInput("ora_keytype", "Gene ID",
                               choices = c("SYMBOL", "ENSEMBL"),
                               selected = "SYMBOL")),
          tags$div(class = "col-2",
                   selectInput("ora_filter_col", "Filter by",
                               choices = c("Raw p-value" = "pvalue",
                                           "Adjusted (p.adjust)" = "p.adjust"),
                               selected = "pvalue")),
          tags$div(class = "col-2",
                   numericInput("ora_pcut", "Cutoff",
                                value = 0.05, min = 0.001, max = 1, step = 0.01)),
          tags$div(class = "col-2",
                   tags$br(),
                   actionButton("ora_run", "Run ORA",
                                class = "btn-primary w-100",
                                icon = icon("play", lib = "font-awesome")))
        )
      ),
      layout_columns(
        col_widths = c(7, 5),
        card(
          card_header("Top Enriched Terms"),
          plotOutput("ora_dotplot", height = "620px")
        ),
        card(
          card_header(
            tags$span("Enrichment Table"),
            tags$span(style = "float:right;",
                      downloadButton("download_ora", ".xlsx",
                                     class = "btn-outline-primary btn-sm",
                                     icon = icon("download", lib = "font-awesome")))
          ),
          DTOutput("ora_table", height = "620px")
        )
      )
    ),

    # ── Tab 5: Help ───────────────────────────────────────────────────────
    nav_panel(
      "Help",
      icon = icon("circle-question", lib = "font-awesome"),

      card(
        card_header("Synergy Criteria"),
        tags$div(style = "padding: 0.5rem;",
          tags$h6("A gene is called", tags$strong("synergistic UP"), "when ALL 4 criteria are met:"),
          tags$ol(
            tags$li(tags$strong("C vs NT significant:"), " Q < p_cutoff  and  direction = UP"),
            tags$li(tags$strong("C vs A significant:"),  " Q < p_cutoff  and  direction = UP"),
            tags$li(tags$strong("C vs B significant:"),  " Q < p_cutoff  and  direction = UP"),
            tags$li(tags$strong("Increase additivity:"),
                    " Increase(C vs NT) > Increase(A vs NT) + Increase(B vs NT),",
                    " where Increase = FC − 1, FC = 2",
                    tags$sup("log2FC"))
          ),
          tags$h6("For", tags$strong("synergistic DOWN"), ", replace UP with DOWN, and Increase with Decrease (Decrease = 1 − FC)."),
          tags$hr(),
          tags$h6("Required Input Files (5 pairwise comparisons):"),
          tags$ul(
            tags$li(tags$code("C vs NT"), " — Combination treatment vs Control"),
            tags$li(tags$code("A vs NT"), " — Treatment A alone vs Control"),
            tags$li(tags$code("B vs NT"), " — Treatment B alone vs Control"),
            tags$li(tags$code("C vs A"),  " — Combination vs A alone"),
            tags$li(tags$code("C vs B"),  " — Combination vs B alone")
          ),
          tags$hr(),
          tags$h6("Expected File Format:"),
          tags$p("Tab-delimited text file with columns:", tags$code("gene_id"), ",",
                 tags$code("gene_name"), ",", tags$code("log2FC"), ",",
                 tags$code("Pvalue"), ",", tags$code("Qvalue"), ",",
                 tags$code("updown"), ". Compatible with edgeR and DESeq2 output.")
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # ── Data loading ────────────────────────────────────────────────────────────

  # Reactive file path list
  file_paths <- reactive({
    if (input$input_mode == "upload") {
      paths <- list(
        c_vs_nt = input$file_c_vs_nt$datapath,
        a_vs_nt = input$file_a_vs_nt$datapath,
        b_vs_nt = input$file_b_vs_nt$datapath,
        c_vs_a  = input$file_c_vs_a$datapath,
        c_vs_b  = input$file_c_vs_b$datapath
      )
    } else {
      paths <- list(
        c_vs_nt = input$path_c_vs_nt,
        a_vs_nt = input$path_a_vs_nt,
        b_vs_nt = input$path_b_vs_nt,
        c_vs_a  = input$path_c_vs_a,
        c_vs_b  = input$path_c_vs_b
      )
    }
    paths
  })

  # File validation summary
  file_status <- reactive({
    paths <- file_paths()
    vals <- sapply(paths, function(p) {
      if (is.null(p) || nchar(p) == 0) return("empty")
      if (input$input_mode == "local" && !file.exists(p)) return("missing")
      return("ready")
    })
    names(vals) <- c("C vs NT", "A vs NT", "B vs NT", "C vs A", "C vs B")
    vals
  })

  # Read data reactively
  uploaded_data <- reactive({
    paths <- file_paths()
    if (any(sapply(paths, function(x) is.null(x) || nchar(x) == 0))) return(NULL)

    # In upload mode, datapath is set once file is uploaded
    # In local mode, check file existence
    if (input$input_mode == "local") {
      missing <- !sapply(paths, file.exists)
      if (any(missing)) {
        msg <- paste("File(s) not found:", paste(names(paths)[missing], collapse = ", "))
        showNotification(msg, type = "error", duration = 6)
        return(NULL)
      }
    }

    tryCatch({
      loaded <- lapply(paths, function(p) {
        raw <- read.delim(p, header = TRUE, stringsAsFactors = FALSE,
                          check.names = FALSE, quote = "", fill = TRUE)

        # Apply column name mapping: rename user columns -> standard names
        col_map <- c(
          setNames(input$col_gene_id,   "gene_id"),
          setNames(input$col_gene_name, "gene_name"),
          setNames(input$col_log2fc,    "log2FC"),
          setNames(input$col_pvalue,    "Pvalue"),
          setNames(input$col_qvalue,    "Qvalue"),
          setNames(input$col_updown,    "updown")
        )
        # Only rename columns that exist and are different from target
        for (target in names(col_map)) {
          src <- col_map[target]
          if (src != target && src %in% colnames(raw)) {
            colnames(raw)[colnames(raw) == src] <- target
          }
        }

        # Coerce types
        raw$log2FC <- as.numeric(raw$log2FC)
        raw$Pvalue  <- as.numeric(raw$Pvalue)
        raw$Qvalue  <- as.numeric(raw$Qvalue)
        raw$updown  <- as.character(raw$updown)

        # Check required columns
        required <- c("gene_id", "log2FC", "Pvalue", "Qvalue", "updown")
        missing_cols <- setdiff(required, colnames(raw))
        if (length(missing_cols) > 0) {
          stop("Missing columns: ", paste(missing_cols, collapse = ", "))
        }

        # Fill missing gene_name
        if (!"gene_name" %in% colnames(raw)) {
          raw$gene_name <- raw$gene_id
        } else {
          raw$gene_name <- as.character(raw$gene_name)
          raw$gene_name[is.na(raw$gene_name) | raw$gene_name == "" | raw$gene_name == "-"] <-
            raw$gene_id[is.na(raw$gene_name) | raw$gene_name == "" | raw$gene_name == "-"]
        }

        raw
      })
      names(loaded) <- names(paths)
      showNotification("Data loaded successfully", type = "message", duration = 2)
      loaded
    }, error = function(e) {
      showNotification(paste("Error reading files:", e$message),
                       type = "error", duration = 10)
      return(NULL)
    })
  })

  # Data status UI
  output$data_status_ui <- renderUI({
    data <- uploaded_data()
    status <- file_status()
    n_ready <- sum(status == "ready")
    n_total <- length(status)
    labels <- c("C vs NT", "A vs NT", "B vs NT", "C vs A", "C vs B")

    if (!is.null(data)) {
      n_genes <- nrow(data$c_vs_nt)
      tags$div(
        tags$p(class = "status-ready",
               icon("check-circle", lib = "font-awesome"),
               sprintf(" %d files loaded  ·  %s genes detected", n_ready,
                       format(n_genes, big.mark = ","))),
        tags$small(class = "text-muted",
                   paste("Ready for analysis. Adjust parameters in sidebar and click 'Run Analysis'."))
      )
    } else if (n_ready == n_total) {
      tags$p(class = "text-muted", icon("hourglass-half", lib = "font-awesome"),
             " Loading files...")
    } else {
      tags$div(
        lapply(seq_along(status), function(i) {
          s <- status[i]
          icon_name <- switch(s, ready = "check", empty = "circle", missing = "times-circle")
          icon_color <- switch(s, ready = "#3A8C5C", empty = "#CCCCCC", missing = "#D64545")
          tags$p(style = sprintf("margin: 0.15rem 0; color: %s;", icon_color),
                 icon(icon_name, lib = "font-awesome"), " ", labels[i])
        })
      )
    }
  })

  # Preview table
  output$preview_table <- renderDT({
    data <- uploaded_data()
    if (is.null(data)) return(NULL)

    idx <- switch(input$preview_file,
      "C vs NT" = "c_vs_nt", "A vs NT" = "a_vs_nt", "B vs NT" = "b_vs_nt",
      "C vs A"  = "c_vs_a",  "C vs B"  = "c_vs_b"
    )
    df <- data[[idx]]
    display_cols <- intersect(
      c("gene_id", "gene_name", "log2FC", "Pvalue", "Qvalue", "updown"),
      colnames(df)
    )
    datatable(df[, display_cols, drop = FALSE],
              options = list(pageLength = 8, scrollX = TRUE, dom = "tip"),
              rownames = FALSE, class = "compact") |>
      formatRound("log2FC", 4) |>
      formatSignif(c("Pvalue", "Qvalue"), 3)
  })

  # ── Synergy analysis ────────────────────────────────────────────────────────

  analysis_done <- reactiveVal(FALSE)

  synergy_res <- eventReactive(input$run_btn, {
    data <- uploaded_data()
    if (is.null(data)) {
      showNotification("Please provide all 5 comparison files first.", type = "warning")
      analysis_done(FALSE)
      return(NULL)
    }

    labels <- c(nt = input$label_nt, a = input$label_a,
                b = input$label_b, c = input$label_c)

    withProgress(message = "Running synergy analysis...", value = 0.5, {
      tryCatch({
        res <- calculate_synergy(
          results_list = data,
          p_cutoff     = input$p_cutoff,
          fc_cutoff    = input$fc_cutoff,
          use_qvalue   = (input$pval_col == "Qvalue"),
          mode         = input$mode,
          labels       = labels
        )
        analysis_done(TRUE)
        showNotification(
          sprintf("Found %d UP + %d DOWN synergistic genes",
                  res$summary$n_synergy_up, res$summary$n_synergy_down),
          type = "message", duration = 5
        )
        res
      }, error = function(e) {
        analysis_done(FALSE)
        showNotification(paste("Analysis failed:", e$message),
                         type = "error", duration = 10)
        return(NULL)
      })
    })
  })

  output$analysis_done <- reactive({ analysis_done() })
  outputOptions(output, "analysis_done", suspendWhenHidden = FALSE)

  # ── Value boxes ─────────────────────────────────────────────────────────────

  output$total_genes <- renderText({
    res <- synergy_res()
    if (is.null(res)) "—" else format(res$summary$n_total_genes, big.mark = ",")
  })
  output$n_synergy_up <- renderText({
    res <- synergy_res()
    if (is.null(res)) "—" else res$summary$n_synergy_up
  })
  output$n_synergy_down <- renderText({
    res <- synergy_res()
    if (is.null(res)) "—" else res$summary$n_synergy_down
  })

  # ── Synergy table ───────────────────────────────────────────────────────────

  output$synergy_table <- renderDT({
    res <- synergy_res()
    if (is.null(res)) return(NULL)

    df <- if (input$synergy_direction == "UP") res$synergy_up else res$synergy_down
    if (nrow(df) == 0) {
      return(datatable(
        data.frame(Message = sprintf("No synergistic %s genes found at this threshold.",
                                      input$synergy_direction)),
        options = list(dom = "t"), rownames = FALSE
      ))
    }

    labels <- res$params$labels
    pval_col <- if (res$params$use_qvalue) "Qvalue" else "Pvalue"
    show_cols <- c("gene_name",
                   paste0("log2FC_c_vs_nt"), paste0("log2FC_a_vs_nt"), paste0("log2FC_b_vs_nt"),
                   paste0("log2FC_c_vs_a"), paste0("log2FC_c_vs_b"),
                   paste0(pval_col, "_c_vs_nt"),
                   paste0(pval_col, "_c_vs_a"),
                   paste0(pval_col, "_c_vs_b"))

    if (input$synergy_direction == "UP") {
      show_cols <- c(show_cols, "Increase_c_vs_nt", "Increase_sum_ab")
    } else {
      show_cols <- c(show_cols, "Decrease_c_vs_nt", "Decrease_sum_ab")
    }

    existing <- intersect(show_cols, colnames(df))
    log2fc_cols <- grep("^log2FC", existing, value = TRUE)
    pq_cols     <- grep("^(Qvalue|Pvalue)", existing, value = TRUE)
    effect_cols <- intersect(existing,
                             c("Increase_c_vs_nt", "Increase_sum_ab",
                               "Decrease_c_vs_nt", "Decrease_sum_ab"))

    dt <- datatable(df[, existing, drop = FALSE],
                    options = list(pageLength = 15, scrollX = TRUE, dom = "ltip"),
                    rownames = FALSE, selection = "single", class = "compact")
    if (length(log2fc_cols) > 0) dt <- formatRound(dt, log2fc_cols, 3)
    if (length(pq_cols)     > 0) dt <- formatSignif(dt, pq_cols, 3)
    if (length(effect_cols) > 0) dt <- formatRound(dt, effect_cols, 4)
    dt
  })

  # ── Volcano plot ────────────────────────────────────────────────────────────

  # Pre-compute volcano data once per analysis result (heavy work cached here)
  volcano_data <- reactive({
    res <- synergy_res()
    if (is.null(res)) return(NULL)

    merged <- res$merged_all
    pval_col <- if (res$params$use_qvalue) "Qvalue" else "Pvalue"
    pval_name <- if (res$params$use_qvalue) "Q-value" else "P-value"
    pval_c <- paste0(pval_col, "_c_vs_nt")

    category <- rep("Other", nrow(merged))
    category[merged$gene_id %in% res$synergy_up$gene_id]   <- "Synergy UP"
    category[merged$gene_id %in% res$synergy_down$gene_id] <- "Synergy DOWN"

    log2fc <- merged$log2FC_c_vs_nt
    log2fc[is.infinite(log2fc)] <- NA
    nlog10 <- -log10(merged[[pval_c]])

    hover <- sprintf("%s\nlog2FC: %.3f\n%s: %.2e",
                     merged$gene_name, merged$log2FC_c_vs_nt,
                     pval_name, merged[[pval_c]])

    list(
      df = data.frame(
        log2fc = log2fc,
        nlog10 = nlog10,
        category = factor(category,
                          levels = c("Synergy UP", "Synergy DOWN", "Other")),
        hover = hover,
        stringsAsFactors = FALSE
      ),
      pval_name = pval_name,
      p_cutoff  = res$params$p_cutoff
    )
  })

  output$volcano_plot <- renderPlotly({
    vd <- volcano_data()
    if (is.null(vd)) return(NULL)

    df <- vd$df
    col_map <- c("Synergy UP"   = COLORS$up,
                 "Synergy DOWN" = COLORS$down,
                 "Other"        = "#C8CFD8")

    # Use scattergl (WebGL) — orders of magnitude faster than SVG for 20K+ points
    plot_ly(data = df, x = ~log2fc, y = ~nlog10,
            color = ~category, colors = col_map,
            type = "scattergl", mode = "markers",
            marker = list(size = 7, opacity = 0.65,
                          line = list(width = 0.4, color = "white")),
            text = ~hover, hoverinfo = "text",
            hovertemplate = "%{text}<extra></extra>") |>
      layout(
        xaxis = list(title = "log2FC  (C vs NT)", zeroline = TRUE,
                     zerolinecolor = "#DDD", gridcolor = "#F0F0F0"),
        yaxis = list(title = sprintf("-log10(%s)", vd$pval_name),
                     gridcolor = "#F0F0F0"),
        shapes = list(list(
          type = "line",
          x0 = min(df$log2fc, na.rm = TRUE),
          x1 = max(df$log2fc, na.rm = TRUE),
          y0 = -log10(vd$p_cutoff), y1 = -log10(vd$p_cutoff),
          line = list(dash = "dash", color = "#999", width = 1)
        )),
        legend = list(orientation = "h", y = 1.1, x = 0.5, xanchor = "center",
                      font = list(size = 11)),
        margin = list(t = 30),
        font = list(family = "Inter, sans-serif")
      ) |>
      config(displaylogo = FALSE,
             modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d"))
  }) |>
    bindCache(synergy_res())

  # ── Contribution chart ──────────────────────────────────────────────────────

  output$contrib_plot <- renderPlot({
    res <- synergy_res()
    if (is.null(res)) return(NULL)
    dir <- tolower(input$contrib_dir)
    plot_synergy_contrib(res, n = 20, direction = dir)
  }) |>
    bindCache(synergy_res(), input$contrib_dir)

  # ── Gene detail inspector ───────────────────────────────────────────────────

  observe({
    res <- synergy_res()
    if (is.null(res)) return()
    genes <- unique(c(res$synergy_up$gene_name, res$synergy_down$gene_name))
    updateSelectizeInput(session, "detail_gene", choices = genes, server = TRUE)
  })

  # Visual: bar chart of log2FC across the three vs-NT comparisons,
  # with a reference line for Increase(A)+Increase(B) — the additivity bar
  output$gene_detail_plot <- renderPlot({
    res <- synergy_res(); gene <- input$detail_gene
    if (is.null(res) || is.null(gene) || gene == "") return(NULL)

    row <- res$merged_all[res$merged_all$gene_name == gene, ]
    if (nrow(row) == 0) return(NULL)
    row <- row[1, ]

    lbl <- res$params$labels
    is_up <- gene %in% res$synergy_up$gene_name

    if (is_up) {
      metric_label <- "Increase (FC − 1)"
      vals <- c(row$Increase_a_vs_nt, row$Increase_b_vs_nt, row$Increase_c_vs_nt)
      sum_ab <- row$Increase_sum_ab
    } else {
      metric_label <- "Decrease (1 − FC)"
      vals <- c(row$Decrease_a_vs_nt, row$Decrease_b_vs_nt, row$Decrease_c_vs_nt)
      sum_ab <- row$Decrease_sum_ab
    }

    comp_names <- c(paste0(lbl["a"], " vs ", lbl["nt"]),
                    paste0(lbl["b"], " vs ", lbl["nt"]),
                    paste0(lbl["c"], " vs ", lbl["nt"]))
    bar_df <- data.frame(
      comp = factor(comp_names, levels = comp_names),
      val  = vals,
      fill = c("A", "B", "C"),
      stringsAsFactors = FALSE
    )

    fill_cols <- c(A = "#76B041", B = "#E45756", C = "#9467BD")

    ggplot2::ggplot(bar_df, ggplot2::aes(x = comp, y = val, fill = fill)) +
      ggplot2::geom_col(width = 0.6, show.legend = FALSE) +
      ggplot2::geom_hline(yintercept = sum_ab, linetype = "dashed",
                          color = "#444", linewidth = 0.7) +
      ggplot2::annotate("text", x = 3, y = sum_ab,
                       label = sprintf("Sum (A+B) = %.2f", sum_ab),
                       hjust = 1.05, vjust = -0.6, size = 5, color = "#444") +
      ggplot2::scale_fill_manual(values = fill_cols) +
      ggplot2::labs(x = NULL, y = metric_label,
                    title = sprintf("%s — %s",
                                    gene,
                                    if (is_up) "Synergy UP" else "Synergy DOWN")) +
      ggplot2::theme_minimal(base_size = 15) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "italic", size = 17),
        axis.title.y = ggplot2::element_text(size = 14),
        axis.text.x = ggplot2::element_text(size = 13),
        axis.text.y = ggplot2::element_text(size = 12),
        panel.grid.major.x = ggplot2::element_blank()
      )
  }) |>
    bindCache(synergy_res(), input$detail_gene)

  output$gene_detail_table <- renderTable({
    res <- synergy_res(); gene <- input$detail_gene
    if (is.null(res) || is.null(gene) || gene == "") return(NULL)

    row <- res$merged_all[res$merged_all$gene_name == gene, ]
    if (nrow(row) == 0) return(data.frame(Note = "Gene not found"))
    row <- row[1, ]
    lbl <- res$params$labels
    pval_col   <- if (res$params$use_qvalue) "Qvalue" else "Pvalue"
    pval_label <- if (res$params$use_qvalue) "Q" else "P"

    is_up <- gene %in% res$synergy_up$gene_name

    if (is_up) {
      effect_label <- "Increase"
      eff_c  <- row$Increase_c_vs_nt
      eff_a  <- row$Increase_a_vs_nt
      eff_b  <- row$Increase_b_vs_nt
      eff_ab <- row$Increase_sum_ab
    } else {
      effect_label <- "Decrease"
      eff_c  <- row$Decrease_c_vs_nt
      eff_a  <- row$Decrease_a_vs_nt
      eff_b  <- row$Decrease_b_vs_nt
      eff_ab <- row$Decrease_sum_ab
    }
    synergy_effect <- eff_c - eff_ab

    data.frame(
      Metric = c(
        sprintf("%s (%s vs %s)", effect_label, lbl["c"], lbl["nt"]),
        sprintf("%s (%s vs %s)", effect_label, lbl["a"], lbl["nt"]),
        sprintf("%s (%s vs %s)", effect_label, lbl["b"], lbl["nt"]),
        sprintf("Synergy Effect  [%s(%s) − %s(%s)+%s(%s)]",
                effect_label, lbl["c"], effect_label, lbl["a"], effect_label, lbl["b"]),
        sprintf("%s (%s vs %s)", pval_label, lbl["c"], lbl["nt"]),
        sprintf("%s (%s vs %s)", pval_label, lbl["c"], lbl["a"]),
        sprintf("%s (%s vs %s)", pval_label, lbl["c"], lbl["b"])
      ),
      Value = c(
        sprintf("%.3f", eff_c),
        sprintf("%.3f", eff_a),
        sprintf("%.3f", eff_b),
        sprintf("%.3f", synergy_effect),
        format(row[[paste0(pval_col, "_c_vs_nt")]], digits = 3),
        format(row[[paste0(pval_col, "_c_vs_a")]],  digits = 3),
        format(row[[paste0(pval_col, "_c_vs_b")]],  digits = 3)
      )
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s", width = "100%")

  # ── Enrichment Analysis (ORA) ───────────────────────────────────────────────

  ora_result <- eventReactive(input$ora_run, {
    res <- synergy_res()
    if (is.null(res)) {
      showNotification("Run synergy analysis first.", type = "warning")
      return(NULL)
    }

    gene_df <- if (input$ora_dir == "up") res$synergy_up else res$synergy_down
    if (nrow(gene_df) == 0) {
      showNotification(sprintf("No synergy %s genes to enrich.", toupper(input$ora_dir)),
                       type = "warning")
      return(NULL)
    }

    keytype <- input$ora_keytype
    gene_ids <- if (keytype == "SYMBOL") gene_df$gene_name else gene_df$gene_id

    # Strip ENSEMBL version suffix if present
    if (keytype == "ENSEMBL") gene_ids <- sub("\\..*$", "", gene_ids)

    gene_ids <- unique(gene_ids[!is.na(gene_ids) & nzchar(gene_ids)])
    if (length(gene_ids) < 3) {
      showNotification("Too few genes for enrichment.", type = "warning")
      return(NULL)
    }

    withProgress(message = "Running ORA...", value = 0.3, {
      tryCatch({
        # Convert to ENTREZ if needed (for KEGG) or pass through for GO
        entrez <- clusterProfiler::bitr(
          gene_ids, fromType = keytype, toType = "ENTREZID",
          OrgDb = org.Hs.eg.db::org.Hs.eg.db
        )$ENTREZID
        entrez <- unique(entrez)

        incProgress(0.4)

        # Run with no internal cutoff; we apply user-chosen filter post-hoc
        ego <- if (input$ora_db == "gobp") {
          clusterProfiler::enrichGO(
            gene = entrez, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
            keyType = "ENTREZID", ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 1, qvalueCutoff = 1,
            readable = TRUE
          )
        } else {
          clusterProfiler::enrichKEGG(
            gene = entrez, organism = "hsa",
            keyType = "kegg",
            pAdjustMethod = "BH",
            pvalueCutoff = 1, qvalueCutoff = 1
          )
        }

        incProgress(0.3)

        if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
          showNotification("Enrichment returned no terms.", type = "warning")
          return(NULL)
        }

        # Apply user-chosen filter (pvalue or p.adjust)
        filter_col <- input$ora_filter_col
        keep <- ego@result[[filter_col]] <= input$ora_pcut
        keep[is.na(keep)] <- FALSE
        ego@result <- ego@result[keep, , drop = FALSE]

        if (nrow(ego@result) == 0) {
          showNotification(sprintf("No terms with %s ≤ %.3g",
                                    filter_col, input$ora_pcut),
                            type = "warning")
          return(NULL)
        }

        # Order by chosen filter column (ascending)
        ego@result <- ego@result[order(ego@result[[filter_col]]), , drop = FALSE]

        # For KEGG, map ENTREZ → SYMBOL in geneID column for readability
        if (input$ora_db == "kegg") {
          ego <- DOSE::setReadable(ego, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                    keyType = "ENTREZID")
        }
        ego
      }, error = function(e) {
        showNotification(paste("ORA failed:", e$message), type = "error", duration = 10)
        NULL
      })
    })
  })

  output$ora_dotplot <- renderPlot({
    ego <- ora_result()
    if (is.null(ego)) {
      plot.new()
      text(0.5, 0.5,
           "Click 'Run ORA' to compute enrichment.\nResults will appear here.",
           cex = 1.3, col = "#888")
      return()
    }
    enrichplot::dotplot(ego, showCategory = 15,
                        label_format = 60) +
      ggplot2::scale_y_discrete(
        labels = function(x) stringr::str_wrap(x, width = 55)
      ) +
      # Narrow the dot panel: term labels get most of the horizontal room
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0.05, 0.10))
      ) +
      ggplot2::theme(
        plot.margin  = ggplot2::margin(8, 12, 8, 12),
        axis.text.y  = ggplot2::element_text(size = 15, lineheight = 0.95),
        axis.text.x  = ggplot2::element_text(size = 13),
        axis.title.x = ggplot2::element_text(size = 15),
        legend.text  = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 13),
        plot.title   = ggplot2::element_text(size = 16)
      )
  }) |>
    bindCache(input$ora_run, input$ora_dir, input$ora_db,
              input$ora_keytype, input$ora_filter_col, input$ora_pcut)

  output$ora_table <- renderDT({
    ego <- ora_result()
    if (is.null(ego)) return(NULL)
    df <- as.data.frame(ego)
    keep <- intersect(c("ID", "Description", "GeneRatio", "BgRatio",
                        "pvalue", "p.adjust", "qvalue", "Count", "geneID"),
                      colnames(df))
    df <- df[, keep]
    datatable(df,
              options = list(pageLength = 15, scrollX = TRUE, dom = "ltip"),
              rownames = FALSE, class = "compact") |>
      formatSignif(intersect(c("pvalue", "p.adjust", "qvalue"), keep), 3)
  })

  output$download_ora <- downloadHandler(
    filename = function() {
      sprintf("synergy_%s_%s_%s.xlsx",
              input$ora_dir, input$ora_db, Sys.Date())
    },
    content = function(file) {
      ego <- ora_result()
      if (is.null(ego)) return()
      df <- as.data.frame(ego)
      if (!requireNamespace("openxlsx", quietly = TRUE)) {
        utils::write.csv(df, file, row.names = FALSE)
      } else {
        openxlsx::write.xlsx(df, file)
      }
    }
  )

  # ── Download handlers ───────────────────────────────────────────────────────

  output$download_excel <- downloadHandler(
    filename = function() paste0("synergy_results_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      res <- synergy_res()
      if (is.null(res)) return()
      export_synergy_excel(res, file)
    }
  )

  output$download_csv <- downloadHandler(
    filename = function() paste0("synergy_merged_", Sys.Date(), ".csv"),
    content  = function(file) {
      res <- synergy_res()
      if (is.null(res)) return()
      write.csv(res$merged_all, file, row.names = FALSE)
    }
  )

  output$download_report <- downloadHandler(
    filename = function() paste0("synergy_report_", Sys.Date(), ".html"),
    content  = function(file) {
      res <- synergy_res()
      if (is.null(res)) return()
      render_synergy_report(res, file)
    }
  )
}

# ── Run ───────────────────────────────────────────────────────────────────────

shinyApp(ui = ui, server = server)
