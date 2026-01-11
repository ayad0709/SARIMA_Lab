# ui.R — SARIMA Scientific Writing Lab (v2)
# FIX v2: Manual SARIMA validation forecast is forced to align with the test set horizon (h = test length).

library(shiny)
library(shinythemes)
library(shinyjs)



ui <- fluidPage(
  theme = shinytheme("flatly"),
  useShinyjs(),
  tags$head(tags$title("SARIMA Scientific Writing Lab (v2)")),
  
  tags$style(HTML("
                    pre {
                      white-space: pre-wrap;
                      word-wrap: break-word;
                      overflow-wrap: break-word;
                    }
                  ")),
  
  titlePanel("SARIMA Modeling Lab (v2)"),

  sidebarLayout(
    sidebarPanel(
      width = 2,
      h4("A) Data"),
      fileInput("fileData", "Upload data (CSV/XLSX)", accept = c(".csv", ".xls", ".xlsx")),
      uiOutput("dateColUI"),
      uiOutput("valueColUI"),
      helpText("Use a Date column and one numeric column."),

      hr(),
      h4("B) Time structure"),
      selectInput(
        "frequency", "Seasonal frequency (s)",
        choices = c(
          "Monthly (12)" = 12,
          "Quarterly (4)" = 4,
          "Weekly (52)" = 52,
          "Daily (365)" = 365,
          "Daily with weekly seasonality (7)" = 7,
          "Custom" = "other"
        ),
        selected = 12
      ),
      conditionalPanel("input.frequency == 'other'",
        numericInput("customFrequency", "Custom s", value = 12, min = 1, step = 1)
      ),
      checkboxInput("align_regular", "Align to regular calendar grid (recommended)", value = TRUE),
      selectInput(
        "missing_policy", "Missing values handling",
        choices = c(
          "Seasonal interpolation (na.interp)" = "seasonal",
          "Linear interpolation" = "linear",
          "Carry forward/backward (LOCF)" = "locf",
          "Drop missing rows (may break spacing)" = "drop"
        ),
        selected = "seasonal"
      ),

      hr(),
      h4("C) Transform"),
      radioButtons(
        "transform", "Transformation",
        choices = c("None" = "none", "Log (ln y)" = "log", "Box-Cox" = "boxcox"),
        selected = "none"
      ),
      conditionalPanel("input.transform == 'boxcox'",
        numericInput("lambda", "λ (Box-Cox), NA = auto", value = NA, step = 0.1)
      ),

      hr(),
      h4("D) Train split"),
      sliderInput("train_prop", "Training size", min = 0.10, max = 1.00, value = 0.80, step = 0.01),
      helpText("When training = 100%: no test set; forecasts extend into the future."),

      hr(),
      h4("E) Diagnostics defaults"),
      numericInput("diag_lag", "Residual test lag (Ljung-Box)", value = 24, min = 1, step = 1),
      checkboxInput("show_teaching_notes", "Show teaching notes", value = TRUE),

      hr(),
      actionButton("refresh_all", "Refresh workflow", icon = icon("sync"))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        
        tabPanel("0 — Roadmap",
          tabsetPanel(
            tabPanel(
              title = "Roadmap",
              uiOutput("roadmap_ui")
            ),
            
            tabPanel(
              title = "Roadmap (Detailed)",
              uiOutput("roadmap_Detailed_ui")
            )
          )),

        # tabPanel("0 — Roadmap", uiOutput("roadmap_ui")),

        tabPanel("1 — Data & descriptives",
          uiOutput("step1_notes"),
          tabsetPanel(

            tabPanel("Preview", tableOutput("data_preview")),
            tabPanel("Descriptives", tableOutput("basic_stats"), plotOutput("hist_plot", height = 300)),
            tabPanel("Missing & outliers", verbatimTextOutput("missing_text"), tableOutput("outlier_table")),
            tabPanel(
              "APA paragraph",
              div(
                style = "width: 100%; height: 1900px; overflow: auto;",
                verbatimTextOutput("apa_data_paragraph")
              )
            ),
          )
        ),

        tabPanel("2 — Exploration",
          uiOutput("step2_notes"),
          tabsetPanel(
            tabPanel("Time series", plotOutput("plot_series", height = 420)),
            tabPanel("Seasonality", plotOutput("season_plot", height = 340), plotOutput("subseries_plot", height = 320)),
            tabPanel("Autocorrelation",
              fluidRow(
                column(6, plotOutput("acf_plot", height = 300)),
                column(6, plotOutput("pacf_plot", height = 300))
              )
            ),
            

            tabPanel("APA paragraph", verbatimTextOutput("apa_explore_paragraph"))
          )
        ),
        
        
        
        
        # --- ADD: New transformation + diagnostics panel (d, D, log) ---
        
        
        tabPanel("2 — Exploration", br(),
                 sidebarLayout(
                   # ======= OUTER SIDEBAR (keeps only data/transformation + stationarity inputs) =======
                   sidebarPanel(
                     width = 2,
                     
                     # Transformation controls
                     numericInput("d_n",  label = "d (non-seasonal differencing):", min = 0, value = 0),
                     numericInput("DS_n", label = "D (seasonal differencing):",     min = 0, value = 0),
                     checkboxInput("check_box", HTML("<b>log(S(t))</b>"), value = FALSE),
                     
                     # hr(),
                     tags$strong("Stationarity Test (ADF/KPSS)"),
                     
                     selectInput("adfTypeSt2", "ADF model type:",
                                 choices = c("none", "drift", "trend"), selected = "drift"),
                     
                     selectInput("alternd2St", "ADF alternative (tseries):",
                                 choices = c("stationary", "explosive", "regression"),
                                 selected = "stationary"),
                     
                     numericInput("LagOrderADFd2St", "Lag order (k):", min = 0, step = 1, value = 10),
                     selectInput("alphaSt2", "Significance level (α):",
                                 choices = c("0.01", "0.05", "0.1"), selected = "0.05")
                     
                     # NOTE: all plot-related controls moved into the Plot (*) tab below
                   ),
                   
                   # ======= MAIN PANEL with tabset =======
                   mainPanel(
                     width = 10,
                     tabsetPanel(
                       id = "transform_tabs",
                       
                       # First tab (unchanged): diagnostics display
                       tabPanel("d?D?log?(St) (*)", uiOutput("d_D_Log_ts_Choice_UI")),
                       
                       # ======= PLOT TAB WITH ITS OWN SIDEBAR =======
                       tabPanel(
                         "Plot (*)",
                         sidebarLayout(
                           sidebarPanel(
                             width = 3,
                             
                             # tags$strong("Plot type"),
                             selectInput(
                               "plot_type_choice", "Plot type:",
                               choices = c(
                                 "Line",
                                 "Points",
                                 "Line + Points",
                                 "Smoothed (LOESS)",
                                 "Moving average",
                                 "Cumulative sum",
                                 "Seasonal plot",
                                 "Seasonal subseries",
                                 "Polar seasonal",
                                 "Seasonal boxplot",
                                 "Classical decomposition (additive)",
                                 "Classical decomposition (multiplicative)",
                                 "STL decomposition",
                                 "Histogram",
                                 "Density",
                                 "QQ plot",
                                 "Lag-1 scatter",
                                 "Lag plot (1..m)",
                                 "ACF",
                                 "PACF",
                                 "ACF+PACF",
                                 "Time + ACF+PACF",
                                 "Periodogram"
                               ),
                               selected = "Line"
                             ),
                             
                             # hr(),
                             
                             # tags$strong("Plot size (optional)"),
                             textInput("plot_width",  "Width (e.g., 800 or 100%)",  value = "800"),
                             textInput("plot_height", "Height (e.g., 500 or 100%)", value = "500"),
                            
                             # tags$strong("Theme"),
                             selectInput(
                               "plot_theme", "Theme:",
                               choices = c("Gray", "Minimal", "Classic", "Light", "Dark", "BW", "Void"),
                               selected = "Minimal"
                             ),
                             
                             # hr(),

                             numericInput("ma_k",  "Moving average window k:", min = 2, step = 1, value = 5),
                             numericInput("lag_m", "Lag plot: number of lags (m):", min = 1, step = 1, value = 12),
                             

                             
                             hr(),
                             uiOutput("ts_color_ui"),
                             numericInput("tickSize",
                                          label = HTML("<span style='color:red;'>Tick size:</span>"),
                                          min = 6, value = 12),
                             # color picker produced server-side
                             
                           ),
                           mainPanel(
                             width = 9,
                             uiOutput("tsPlot_Choice_UI")
                           )
                         )
                       ),
                       
                       # ACF+PACF tab
                       tabPanel("ACF+PACF (*)", uiOutput("difference2ACFPACF_UI")),
                       
                       # Stationarity report
                       tabPanel("Stationarity [ADF + KPSS]",
                                tags$head(tags$style(HTML("
                                  #teststationarited3St{
                                    height: 745px !important;
                                    max-height: 800px;
                                    width: 100% !important;
                                    white-space: pre;
                                    overflow-y: auto;
                                    border: 2px solid #cccccc;
                                    font-size: 13px;
                                  }
                                "))),
                                verbatimTextOutput("teststationarited3St")
                       ),
                       
                       # Checklist tab
                       tabPanel("CHECKLIST & STEPS",
                                tags$head(tags$style(HTML("
                                  #CHECKLIST{
                                    height: 745px !important;
                                    max-height: 800px;
                                    width: 100% !important;
                                    white-space: pre;
                                    overflow-y: auto;
                                    border: 2px solid #cccccc;
                                    font-size: 13px;
                                  }
                                "))),
                                verbatimTextOutput("CHECKLIST")
                       ),
                       
                       # Auto-ARIMA tab
                       # tabPanel("Auto-ARIMA", verbatimTextOutput("ARIMA_d_D_log"))
                     )
                   )
                 )
        ),
        


        tabPanel("3 — Decomposition",
          uiOutput("step3_notes"),
          tabsetPanel(
            tabPanel("Classical additive", plotOutput("decomp_add", height = 420)),
            tabPanel("Classical multiplicative", plotOutput("decomp_mul", height = 420)),
            tabPanel("STL (robust)", plotOutput("stl_plot", height = 520)),
            tabPanel("APA paragraph", verbatimTextOutput("apa_decomp_paragraph"))
          )
        ),

        tabPanel("4 — Stationarity & differencing",
          uiOutput("step4_notes"),
          sidebarLayout(
            sidebarPanel(
              width = 2,
              h5("Stationarity tests"),
              selectInput("adf_type", "ADF type (ur.df)", choices = c("none", "drift", "trend"), selected = "trend"),
              numericInput("adf_lags", "ADF lags (k)", value = 10, min = 0, step = 1),
              selectInput("kpss_type", "KPSS null", choices = c("mu", "tau"), selected = "mu"),
              actionButton("run_tests", "Run tests", icon = icon("vial")),
              hr(),
              h5("Differencing preview"),
              numericInput("d_preview", "Non-seasonal difference (d)", value = 0, min = 0, step = 1),
              numericInput("D_preview", "Seasonal difference (D)", value = 0, min = 0, step = 1),
              # actionButton("preview_diff", "Preview differenced  series", icon = icon("chart-line")),
              actionButton(
                "preview_diff",
                label = HTML("Preview<br>differenced<br>series"),
                icon = icon("chart-line")
              )
            ),
            mainPanel(
              width = 8,
              tabsetPanel(
                tabPanel("Results", verbatimTextOutput("stationarity_results")),
                tabPanel("Interpretation", verbatimTextOutput("stationarity_interpretation")),
                tabPanel("Suggested differencing", verbatimTextOutput("diff_suggestion")),
                tabPanel("Differenced plot", plotOutput("diff_plot", height = 360)),
                tabPanel("APA paragraph", verbatimTextOutput("apa_stationarity_paragraph"))
              )
            )
          )
        ),

        tabPanel("5 — Auto-ARIMA (baseline)",
          uiOutput("step5_notes"),
          sidebarLayout(
            sidebarPanel(
              width = 3,
              checkboxInput("auto_seasonal", "Allow seasonal terms", value = TRUE),
              checkboxInput("auto_stepwise", "Stepwise search", value = TRUE),
              checkboxInput("auto_approx", "Approximation", value = FALSE),
              checkboxInput("auto_allow_mean", "Allow mean/drift", value = TRUE),
              numericInput("auto_max_order", "max.order", value = 10, min = 3, step = 1),
              numericInput("auto_h", "Future horizon h (used when no test set)", value = NA, min = 1, step = 1),
              actionButton("fit_auto", "Fit Auto-ARIMA", icon = icon("magic"))
            ),
            mainPanel(
              width = 8,
              tabsetPanel(
                tabPanel("Model specification", verbatimTextOutput("auto_model_spec"), tableOutput("auto_coef_table")),
                tabPanel("Diagnostics",
                  fluidRow(
                    column(6, plotOutput("auto_resid_ts", height = 240)),
                    column(6, plotOutput("auto_resid_acf", height = 240))
                  ),
                  fluidRow(
                    column(6, plotOutput("auto_resid_hist", height = 240)),
                    column(6, plotOutput("auto_resid_qq", height = 240))
                  ),
                  verbatimTextOutput("auto_diag_tests")
                ),
                tabPanel("Forecast & accuracy",
                  verbatimTextOutput("auto_horizon_note"),
                  plotOutput("auto_forecast_plot", height = 420),
                  tableOutput("auto_accuracy_table"),
                  tableOutput("auto_forecast_table")
                ),
                tabPanel("APA paragraph", verbatimTextOutput("apa_auto_paragraph"))
              )
            )
          )
        ),

        tabPanel("6 — Manual SARIMA (theory-driven)",
          uiOutput("step6_notes"),
          sidebarLayout(
            sidebarPanel(
              width = 2,
                            actionButton(
                "fit_manual",
                label = HTML("__Fit__"),
                # label = HTML("Fit<br>Manual<br>SARIMA"),
                icon = icon("cogs")
              ),
              h5("Non-seasonal (p,d,q)"),
              numericInput("p", "p", value = 1, min = 0, step = 1),
              numericInput("d", "d", value = 0, min = 0, step = 1),
              numericInput("q", "q", value = 1, min = 0, step = 1),
              hr(),
              h5("Seasonal (P,D,Q,s)"),
              numericInput("P", "P", value = 0, min = 0, step = 1),
              numericInput("D", "D", value = 0, min = 0, step = 1),
              numericInput("Q", "Q", value = 0, min = 0, step = 1),
              numericInput("s", "Season length s (NA = sidebar frequency)", value = NA, min = 1, step = 1),
              checkboxInput("manual_drift", "Include drift/mean", value = FALSE),
              hr(),
              numericInput("manual_h", "Future horizon h (used when no test set)", value = NA, min = 1, step = 1),
             
              # actionButton("fit_manual", "Fit Manual SARIMA", icon = icon("cogs"))
              

            ),
            mainPanel(
              width = 8,
              tabsetPanel(
                tabPanel("Train/Test split",
                  verbatimTextOutput("manual_split_text"),
                  plotOutput("manual_split_plot", height = 320)
                ),
                tabPanel("Model specification", verbatimTextOutput("manual_model_spec"), tableOutput("manual_coef_table")),
                
                # --- MOD: equation panel ---
                tabPanel(
                  "Model equation",
                  withMathJax(
                    uiOutput("manual_model_equation")
                  )
                ),
                
                tabPanel("Diagnostics",
                  fluidRow(
                    column(6, plotOutput("manual_resid_ts", height = 240)),
                    column(6, plotOutput("manual_resid_acf", height = 240))
                  ),
                  fluidRow(
                    column(6, plotOutput("manual_resid_hist", height = 240)),
                    column(6, plotOutput("manual_resid_qq", height = 240))
                  ),
                  verbatimTextOutput("manual_diag_tests")
                ),
                tabPanel("Forecast & accuracy",
                  verbatimTextOutput("manual_horizon_note"),
                  plotOutput("manual_forecast_plot", height = 420),
                  tableOutput("manual_accuracy_table"),
                  tableOutput("manual_forecast_table")
                ),
                tabPanel("APA paragraph", verbatimTextOutput("apa_manual_paragraph"))
              )
            )
          )
        ),

        tabPanel("7 — Comparison & Paper builder",
          uiOutput("step7_notes"),
          tabsetPanel(
            tabPanel("Model comparison", tableOutput("comparison_table"), verbatimTextOutput("comparison_interpretation")),
            tabPanel("Methods (APA draft)", verbatimTextOutput("apa_methods_draft")),
            tabPanel("Results (APA draft)", verbatimTextOutput("apa_results_draft")),
            tabPanel("Checklist", uiOutput("paper_checklist_ui"))
          )
        )
      )
    )
  )
)
