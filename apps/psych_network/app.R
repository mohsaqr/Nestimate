## Psychometric Networks — minimal Shiny front-end for Nestimate
## Standalone app. Not part of the package build.
## Run with: shiny::runApp("apps/psych_network")

library(shiny)
library(Nestimate)
library(cograph)

## Built-in demo: SRL constructs (5,420 obs of 5 Likert-derived scales).
## Resolve relative to the app file so it works under both runApp() and shiny-server.
.app_dir   <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
.demo_path <- file.path(.app_dir, "data", "srl_demo.csv")
if (!file.exists(.demo_path)) .demo_path <- file.path("data", "srl_demo.csv")
srl_demo <- read.csv(.demo_path)

## Drop columns that are constant or all-NA — those break cor/glasso/ising.
.clean_numeric <- function(df) {
  keep <- vapply(df, function(x) {
    if (!is.numeric(x)) return(FALSE)
    if (all(is.na(x)))  return(FALSE)
    sd_val <- stats::sd(x, na.rm = TRUE)
    is.finite(sd_val) && sd_val > 0
  }, logical(1))
  df[, keep, drop = FALSE]
}

## Long labels for the variable picker. Column names stay short for plot nodes.
srl_labels <- c(
  "Intrinsic Value (IV)"          = "IV",
  "Self-Efficacy (SE)"            = "SE",
  "Test Anxiety (TA)"             = "TA",
  "Cognitive Strategy Use (CSU)"  = "CSU",
  "Self-Regulation (SR)"          = "SR"
)

guide_md <- "
### What this app does

It estimates a **psychometric network** from a numeric data table — one column per
variable, one row per observation — and draws it with `cograph::splot()`. Three
estimators are exposed, all standard in the network-psychometrics literature:

- **`cor`** — zero-order Pearson correlations. Useful as a baseline; every pair has an
  edge, so this is the *unconditional* association structure.
- **`pcor`** — partial correlations. Each edge is the association between two
  variables *after controlling for every other variable in the table*. This is what
  most psychometricians mean by 'a network'.
- **`glasso`** — graphical lasso with EBIC tuning (Epskamp & Fried, 2018). Like
  `pcor` but with regularization that shrinks small partials to exactly zero,
  giving a sparser, more replicable network. The `gamma` slider controls the EBIC
  hyperparameter (0 = BIC, 0.5 = standard psychometric default, higher = sparser).
- **`ising`** — L1-penalized logistic regressions for *binary* data only (van
  Borkulo et al., 2014). Provided for symptom-level networks. Continuous columns
  are auto-binarized at the median when this method is chosen.

### How to use

1. Pick a data source: the built-in **SRL** dataset (5,420 observations of five
   self-regulated-learning constructs — Intrinsic Value (IV), Self-Efficacy (SE),
   Test Anxiety (TA), Cognitive Strategy Use (CSU), Self-Regulation (SR)) or
   upload your own CSV. The app keeps only the numeric columns.
2. Tick the **variables** to include in the network. All are selected by default;
   untick any you want to leave out without modifying the source data.
3. Choose an estimator. Start with `glasso` if you don't have a strong prior.
4. Set a threshold to hide weak edges in the plot (does not refit; purely visual).
5. Click **Estimate**. The view jumps to **Summary**; switch to **Plot** to see
   the network and **Edges** for the sorted edge list.

### How to read the plot

- Node size and color follow `splot()` defaults. Edge thickness is the absolute
  edge weight; edge color encodes sign (positive vs negative partial / correlation).
- Two strongly connected variables share a thick edge. A node with many connections
  is *central* — see `Nestimate::centrality()` for numeric centrality measures
  (not exposed in this UI).

### Layouts

- **Oval / Ellipse / Circle** — geometric layouts that place nodes on a closed
  curve. Useful when you want a clean, deterministic plot and centrality is
  carried by edges, not coordinates.
- **Spring / Fruchterman-Reingold / Gephi** — force-directed layouts. Strongly
  connected nodes pull together, weakly connected nodes drift to the periphery.
  Best for revealing clusters. These are stochastic — change the **seed** to
  re-roll a layout you don't like.
- **Grid / Star / Random** — diagnostic layouts mostly useful for inspection or
  presentation, not inference.

### What this app does **not** do

Bootstrap stability, case-drop CS-coefficients, network comparison tests, and
multilevel decomposition all live in Nestimate but are out of scope for this
single-file demo. Call `bootstrap_network()`, `centrality_stability()`, or
`nct()` from the R console for those.
"

footer <- tags$footer(
  style = paste(
    "text-align: center;",
    "padding: 14px 12px;",
    "margin-top: 28px;",
    "border-top: 1px solid #e5e5e5;",
    "font-size: 0.85em;",
    "color: #666;"
  ),
  HTML(paste0(
    "Built with <strong>Nestimate</strong> and <strong>cograph</strong> &middot; ",
    "Developed by <a href='https://saqr.me' target='_blank' rel='noopener'>Mohammed Saqr</a> &middot; ",
    "<a href='https://saqr.me/Nestimate/' target='_blank' rel='noopener'>Docs</a> &middot; ",
    "<a href='https://github.com/mohsaqr/Nestimate' target='_blank' rel='noopener'>GitHub</a>"
  ))
)

ui <- fluidPage(
  titlePanel("Psychometric Networks"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      radioButtons(
        "source", "Data source",
        choices = c("Built-in SRL (IV, SE, TA, CSU, SR)" = "builtin",
                    "Upload CSV" = "upload"),
        selected = "builtin"
      ),
      conditionalPanel(
        "input.source == 'upload'",
        fileInput("file", "CSV file", accept = ".csv"),
        checkboxInput("header", "Header row", value = TRUE)
      ),
      uiOutput("var_selector"),
      selectInput(
        "method", "Estimator",
        choices = c("Pearson correlation" = "cor",
                    "Partial correlation" = "pcor",
                    "EBIC-glasso" = "glasso",
                    "Ising (binary)" = "ising"),
        selected = "glasso"
      ),
      conditionalPanel(
        "input.method == 'glasso'",
        sliderInput("gamma", "EBIC gamma", min = 0, max = 1, value = 0.5, step = 0.05)
      ),
      sliderInput("threshold", "Plot threshold (|edge|)",
                  min = 0, max = 1, value = 0, step = 0.01),
      selectInput(
        "layout", "Plot layout",
        choices = c(
          "Oval (psych default)"     = "oval",
          "Circle"                   = "circle",
          "Ellipse"                  = "ellipse",
          "Spring"                   = "spring",
          "Fruchterman-Reingold"     = "fr",
          "Gephi (ForceAtlas-style)" = "gephi",
          "Grid"                     = "grid",
          "Star"                     = "star",
          "Random"                   = "random"
        ),
        selected = "oval"
      ),
      conditionalPanel(
        "['spring','fr','gephi','random'].indexOf(input.layout) >= 0",
        numericInput("seed", "Layout seed", value = 42, min = 0, step = 1)
      ),
      tags$hr(),
      actionButton("go", "Estimate", class = "btn-primary"),
      tags$hr(),
      tags$strong("Plot appearance"),
      sliderInput("node_size",   "Node size",        min = 4,   max = 25,  value = 10,  step = 1),
      sliderInput("label_size",  "Node label size",  min = 0.5, max = 3,   value = 1.5, step = 0.1),
      checkboxInput("show_edge_labels", "Edge weight", value = TRUE),
      conditionalPanel(
        "input.show_edge_labels == true",
        sliderInput("edge_label_size", "Edge label size",
                    min = 1, max = 6, value = 3, step = 0.5),
        sliderInput("weight_digits", "Edge label digits",
                    min = 1, max = 4, value = 2, step = 1)
      ),
      hr(),
      helpText("Numeric columns only. Non-numeric columns are dropped.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", selected = "guide",
        tabPanel("How to use", value = "guide",
                 div(style = "max-width: 780px; padding: 12px;",
                     HTML(markdown::markdownToHTML(text = guide_md, fragment.only = TRUE)))),
        tabPanel("Summary", value = "summary",
                 verbatimTextOutput("net_summary")),
        tabPanel("Plot", value = "plot",
                 plotOutput("net_plot", height = "640px")),
        tabPanel("Edges", value = "edges",
                 tableOutput("edge_table"))
      )
    )
  ),
  footer
)

server <- function(input, output, session) {

  # All numeric columns of the active dataset — universe for the variable picker.
  # Drops zero-variance / all-NA columns so downstream estimators don't blow up.
  raw_data_all <- reactive({
    df <- if (input$source == "builtin") {
      srl_demo
    } else {
      req(input$file)
      tryCatch(
        read.csv(input$file$datapath, header = input$header),
        error = function(e) {
          validate(need(FALSE, paste0(
            "Could not read CSV: ", conditionMessage(e),
            ". Make sure the file is comma-separated with a single header row."
          )))
        }
      )
    }
    df <- .clean_numeric(df)
    if (ncol(df) < 2L) {
      validate(need(FALSE,
        "Need at least two numeric, non-constant columns. ",
        "Constant or all-NA columns were dropped."
      ))
    }
    df
  })

  # Variable picker: long labels for the SRL built-in, plain column names for uploads.
  output$var_selector <- renderUI({
    df <- raw_data_all()
    cols <- names(df)
    choices <- if (input$source == "builtin") {
      keep <- srl_labels[srl_labels %in% cols]
      c(keep, setNames(setdiff(cols, keep), setdiff(cols, keep)))
    } else {
      setNames(cols, cols)
    }
    checkboxGroupInput(
      "vars", "Variables",
      choices  = choices,
      selected = unname(choices)
    )
  })

  # Selected subset — feeds the fit.
  raw_data <- reactive({
    df <- raw_data_all()
    sel <- intersect(input$vars, names(df))
    if (length(sel) < 2L) {
      validate("Pick at least two variables to estimate a network.")
    }
    df[, sel, drop = FALSE]
  })

  # Fit on Estimate click; isolate inputs so slider tweaks don't refit.
  fit <- eventReactive(input$go, {
    df <- raw_data()
    method <- input$method
    if (method == "ising") {
      df <- as.data.frame(lapply(df, function(x) as.integer(x > stats::median(x, na.rm = TRUE))))
      # Median-binarization can produce all-zero or all-one columns when ties dominate.
      keep <- vapply(df, function(x) length(unique(x)) > 1L, logical(1))
      if (sum(keep) < 2L) {
        validate(need(FALSE,
          "Ising needs at least two columns with both 0 and 1 after median-binarization. ",
          "Pick continuous variables with enough variation, or switch estimator."
        ))
      }
      df <- df[, keep, drop = FALSE]
    }
    # Drop rows with any NA — cor/glasso/ising all error on NaNs.
    if (anyNA(df)) df <- df[stats::complete.cases(df), , drop = FALSE]
    if (nrow(df) < 3L) {
      validate(need(FALSE,
        "Need at least three complete rows after dropping missing values."
      ))
    }
    args <- list(data = df, method = method)
    if (method == "glasso") args$gamma <- input$gamma
    withProgress(message = "Estimating network", value = 0.4, {
      tryCatch(
        do.call(Nestimate::build_network, args),
        error = function(e) {
          validate(need(FALSE, paste0(
            "Estimation failed: ", conditionMessage(e),
            ". Try a different estimator or different variables."
          )))
        }
      )
    })
  })

  # Jump to the Summary tab the moment the user clicks Build.
  observeEvent(input$go, {
    updateTabsetPanel(session, "tabs", selected = "summary")
  })

  # Apply visual threshold without refitting: clone netobject, zero small edges.
  fit_thresh <- reactive({
    net <- fit()
    thr <- input$threshold
    if (thr > 0) {
      W <- net$weights
      W[abs(W) < thr] <- 0
      net$weights <- W
    }
    net
  })

  output$net_plot <- renderPlot({
    validate(need(input$go > 0L,
      "Click Estimate in the sidebar to compute the network."))
    net    <- fit_thresh()
    layout <- input$layout
    seed   <- if (is.null(input$seed) || is.na(input$seed)) 42L else as.integer(input$seed)

    plot_args <- list(
      node_size  = input$node_size,
      label_size = input$label_size,
      seed       = seed
    )
    if (isTRUE(input$show_edge_labels)) {
      plot_args$edge_labels     <- TRUE
      plot_args$edge_label_size <- input$edge_label_size
      plot_args$weight_digits   <- input$weight_digits
    }

    # igraph's Fruchterman-Reingold rejects signed weights (partial correlations
    # are signed). Precompute FR coords on |weights| and pass via layout="custom".
    if (identical(layout, "fr") && any(net$weights < 0, na.rm = TRUE)) {
      net_abs <- net
      net_abs$weights <- abs(net_abs$weights)
      if (!is.null(net_abs$edges) && "weight" %in% names(net_abs$edges)) {
        net_abs$edges$weight <- abs(net_abs$edges$weight)
      }
      pos    <- cograph:::compute_layout_for_cograph(net_abs, layout = "fr", seed = seed)
      coords <- as.matrix(pos$nodes[, c("x", "y")])
      do.call(cograph::splot, c(list(net, layout = "custom", coords = coords), plot_args))
    } else {
      do.call(cograph::splot, c(list(net, layout = layout), plot_args))
    }
  })

  output$net_summary <- renderPrint({
    validate(need(input$go > 0L,
      "Click Estimate in the sidebar to compute the network."))
    print(fit_thresh())
  })

  output$edge_table <- renderTable({
    validate(need(input$go > 0L,
      "Click Estimate in the sidebar to compute the network."))
    net <- fit_thresh()
    W <- net$weights
    nm <- rownames(W)
    pairs <- which(W != 0 & upper.tri(W), arr.ind = TRUE)
    if (!nrow(pairs)) return(data.frame(message = "No edges above threshold."))
    out <- data.frame(
      from   = nm[pairs[, 1]],
      to     = nm[pairs[, 2]],
      weight = W[pairs],
      stringsAsFactors = FALSE
    )
    out[order(-abs(out$weight)), ]
  }, digits = 4)
}

shinyApp(ui, server)
