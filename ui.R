# ui.R — SARIMA Scientific Writing Lab (v2) 
# FIX v2: Manual SARIMA validation forecast is forced to align with the test set horizon (h = test length).

library(shiny)
library(shinythemes)
library(shinyjs)
library(DiagrammeR)

ui <- fluidPage(
  theme = shinytheme("flatly"),
  useShinyjs(),
  
  #  Load MathJax once
  withMathJax(),
  
  #  Everything that belongs in <head> goes here (merged)
  tags$head(
    tags$title("SARIMA Scientific Writing Lab (v2)"),
    
    tags$style(HTML("
    pre {
      white-space: pre-wrap;
      word-wrap: break-word;
      overflow-wrap: break-word;
    }

    /* Reduce spacing between stacked wellPanels */
    .well { 
      margin-bottom: 2px; 
      padding: 4px;
    }

    /* Reduce spacing between inputs inside a wellPanel */
    .form-group {
      margin-bottom: 8px;
    }

    /* Reduce hr spacing */
    hr {
      margin-top: 8px;
      margin-bottom: 8px;
    }
    
    
    /* “left border Color + background Color” callout wrapper around your text */
    .callout {
      padding: 10px;
      border-left: 4px solid;
      border-radius: 6px;
      margin-top: 8px;
      background: #f8fbff;
    }
    .callout-teal   { border-left-color:#14B8A6; background:#F0FDFA; }
    .callout-green  { border-left-color:#22C55E; background:#F0FDF4; }
    .callout-amber  { border-left-color:#F59E0B; background:#FFFBEB; }
    .callout-rose   { border-left-color:#F43F5E; background:#FFF1F2; }
    .callout-purple { border-left-color:#8B5CF6; background:#F5F3FF; }
    .callout-slate  { border-left-color:#64748B; background:#F8FAFC; }
    .callout-blue   { border-left-color:#3B82F6; background:#F8FBFF; }
    .callout-cyan    { border-left-color:#06B6D4; background:#ECFEFF; }
    .callout-sky     { border-left-color:#0EA5E9; background:#F0F9FF; }
    .callout-indigo  { border-left-color:#4F46E5; background:#EEF2FF; }
    .callout-violet  { border-left-color:#7C3AED; background:#F5F3FF; }
    .callout-fuchsia { border-left-color:#D946EF; background:#FDF4FF; }
    .callout-lime    { border-left-color:#84CC16; background:#F7FEE7; }
    .callout-orange  { border-left-color:#F97316; background:#FFF7ED; }
    .callout-red     { border-left-color:#EF4444; background:#FEF2F2; }
    .callout-stone   { border-left-color:#78716C; background:#FAFAF9; }
    .callout-zinc    { border-left-color:#71717A; background:#FAFAFA; }
  ")),
    
    tags$script(HTML("
    Shiny.addCustomMessageHandler('mathjax-typeset', function(id){
      // MathJax v2 (Shiny default)
      if (window.MathJax && MathJax.Hub && MathJax.Hub.Queue) {
        MathJax.Hub.Queue(['Typeset', MathJax.Hub, id]);
      }
      // MathJax v3 (future-proof)
      if (window.MathJax && MathJax.typesetPromise) {
        var el = document.getElementById(id);
        if (el) MathJax.typesetPromise([el]);
      }
    });
  ")),
    
  ),
  
  
  
  tags$style(HTML("
            /* Only style THIS details' own summary (not nested ones) */
            .defs-details > summary {
              list-style: none;
              cursor: pointer;
              font-weight: 600;
              display: block;
              text-align: left;
              padding-left: 1.6em;
              position: relative;
            }
            
            .defs-details > summary::-webkit-details-marker { display: none; }
            
            /* closed */
            .defs-details > summary::before {
              content: '▶';
              position: absolute;
              left: 0;
              top: 50%;
              transform: translateY(-50%);
              line-height: 1;
            }
            
            /* open (IMPORTANT: direct child selector) */
            .defs-details[open] > summary::before {
              content: '▼';
            }
            
            /* optional spacing for content */
            .defs-details > *:not(summary) {
              margin-left: 0.6em;
              margin-top: 6px;
            }
            
            
            /* ============= */
            /* Buttons Color */
            /* ============= */

            /* ===== BLUE ===== */
            .blue-btn-1{background-color:#0d6efd;color:white;border:none}.blue-btn-1:hover{background-color:#0b5ed7;color:white}
            .blue-btn-2{background-color:#2C7BE5;color:white;border:none}.blue-btn-2:hover{background-color:#1A5BB8;color:white}
            .blue-btn-3{background-color:#1976D2;color:white;border:none}.blue-btn-3:hover{background-color:#115293;color:white}
            .blue-btn-4{background-color:#2563EB;color:white;border:none}.blue-btn-4:hover{background-color:#1D4ED8;color:white}
            .blue-btn-5{background-color:#0EA5E9;color:white;border:none}.blue-btn-5:hover{background-color:#0284C7;color:white}
            
            /* ===== GREEN ===== */
            .green-btn-1{background-color:#198754;color:white;border:none}.green-btn-1:hover{background-color:#157347;color:white}
            .green-btn-2{background-color:#2E7D32;color:white;border:none}.green-btn-2:hover{background-color:#1B5E20;color:white}
            .green-btn-3{background-color:#22C55E;color:white;border:none}.green-btn-3:hover{background-color:#16A34A;color:white}
            .green-btn-4{background-color:#16A34A;color:white;border:none}.green-btn-4:hover{background-color:#15803D;color:white}
            .green-btn-5{background-color:#059669;color:white;border:none}.green-btn-5:hover{background-color:#047857;color:white}
            
            /* ===== ORANGE ===== */
            .orange-btn-1{background-color:#fd7e14;color:white;border:none}.orange-btn-1:hover{background-color:#e96b03;color:white}
            .orange-btn-2{background-color:#F59E0B;color:white;border:none}.orange-btn-2:hover{background-color:#D97706;color:white}
            .orange-btn-3{background-color:#FB923C;color:white;border:none}.orange-btn-3:hover{background-color:#F97316;color:white}
            .orange-btn-4{background-color:#EA580C;color:white;border:none}.orange-btn-4:hover{background-color:#C2410C;color:white}
            .orange-btn-5{background-color:#FF9800;color:white;border:none}.orange-btn-5:hover{background-color:#F57C00;color:white}
            
            /* ===== RED ===== */
            .red-btn-1{background-color:#dc3545;color:white;border:none}.red-btn-1:hover{background-color:#bb2d3b;color:white}
            .red-btn-2{background-color:#B91C1C;color:white;border:none}.red-btn-2:hover{background-color:#991B1B;color:white}
            .red-btn-3{background-color:#EF4444;color:white;border:none}.red-btn-3:hover{background-color:#DC2626;color:white}
            .red-btn-4{background-color:#F87171;color:white;border:none}.red-btn-4:hover{background-color:#EF4444;color:white}
            .red-btn-5{background-color:#E53935;color:white;border:none}.red-btn-5:hover{background-color:#C62828;color:white}
            
            /* ===== GRAY ===== */
            .gray-btn-1{background-color:#6c757d;color:white;border:none}.gray-btn-1:hover{background-color:#5c636a;color:white}
            .gray-btn-2{background-color:#4B5563;color:white;border:none}.gray-btn-2:hover{background-color:#374151;color:white}
            .gray-btn-3{background-color:#64748B;color:white;border:none}.gray-btn-3:hover{background-color:#475569;color:white}
            .gray-btn-4{background-color:#334155;color:white;border:none}.gray-btn-4:hover{background-color:#1E293B;color:white}
            .gray-btn-5{background-color:#9CA3AF;color:white;border:none}.gray-btn-5:hover{background-color:#6B7280;color:white}
            
            /* ===== PURPLE ===== */
            .purple-btn-1{background-color:#6f42c1;color:white;border:none}.purple-btn-1:hover{background-color:#59339d;color:white}
            .purple-btn-2{background-color:#7C3AED;color:white;border:none}.purple-btn-2:hover{background-color:#5B21B6;color:white}
            .purple-btn-3{background-color:#8B5CF6;color:white;border:none}.purple-btn-3:hover{background-color:#7C3AED;color:white}
            .purple-btn-4{background-color:#4F46E5;color:white;border:none}.purple-btn-4:hover{background-color:#4338CA;color:white}
            .purple-btn-5{background-color:#A855F7;color:white;border:none}.purple-btn-5:hover{background-color:#9333EA;color:white}
            

            
            
        ")),
  
  
  
  
  
  
  #   🟡 SARIMA Modeling Lab Main Panel 🟡

 
  titlePanel("SARIMA Modeling Lab (v2)"),
  
  sidebarLayout(
    sidebarPanel(
      width = 2,
      
      # h4("A) Data"),
      
      fileInput("fileData", "Upload data (CSV/XLSX)", accept = c(".csv", ".xls", ".xlsx")),
      uiOutput("dateColUI"),
      uiOutput("valueColUI"),
      
      # h4("B) Time structure"),
      
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
      conditionalPanel(
        "input.frequency == 'other'",
        numericInput("customFrequency", "Custom s", value = 12, min = 1, step = 1)
      ),
      checkboxInput("align_regular", "Align to regular calendar grid (recommended)", value = TRUE),
      selectInput(
        "missing_policy", "Missing values handling",
        choices = c(
          "Seasonal interpolation" = "seasonal",
          "Linear interpolation" = "linear",
          "Carry forward/backward" = "locf",
          "Drop missing rows" = "drop"
        ),
        selected = "seasonal"
      ),
      
      # h4("C) Transform"),
      
      radioButtons(
        "transform", "Transformation",
        choices = c("None" = "none", "Log (ln y)" = "log", "Box-Cox" = "boxcox"),
        selected = "none"
      ),
      conditionalPanel(
        "input.transform == 'boxcox'",
        numericInput("lambda", "λ (Box-Cox), NA = auto", value = NA, step = 0.1)
      ),
      
      # h4("D) Train split"),
      hr(),
      h4("Train split"),
      
      sliderInput("train_prop", "Training size", min = 0.10, max = 1.00, value = 1, step = 0.01),
      
      # h4("E) Diagnostics defaults"),
      
      numericInput("diag_lag", "Residual test lag (Ljung-Box)", value = 12, min = 1, step = 1),
      checkboxInput("show_teaching_notes", "Show teaching notes", value = TRUE),
      
      hr(),
      
      # tags$div(
      #   style = "text-align: center;",
      #   actionButton(
      #     "refresh_all",
      #     label = HTML("Refresh workflow"),
      #     icon  = icon("sync"),
      #     width = "90%"
      #   )
      # ),
      
      tags$div(
        style = "text-align: center;",
        actionButton(
          "refresh_all",
          label = HTML("Refresh workflow"),
          icon  = icon("sync"),
          width = "90%",
          # class = "gray-btn-1"
        )
      ),

      
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        
        
        #--------------
        #  0️ Roadmap
        #--------------
        
        tabPanel(
          "Roadmap",
          tabsetPanel(
            tabPanel(title = "Roadmap", uiOutput("roadmap_ui")),
            # tabPanel(title = "Roadmap (Detailed Ang)", uiOutput("roadmap_Detailed_Ang_ui")),
            # tabPanel(title = "Roadmap (Detailed Fr)", uiOutput("roadmap_Detailed_Fr_ui")),
            # tabPanel(title = "Roadmap (Detailed Fr)", uiOutput("roadmap_Detailed_Fr_ui2")),
            # tabPanel(title = "Roadmap (Detailed Fr)", uiOutput("roadmap_Detailed_Fr_ui3")),
            tabPanel(title = "Roadmap Fr", uiOutput("roadmap_Detailed_Fr_ui4")),
            # tabPanel(title = "Roadmap Fr", uiOutput("roadmap_Detailed_Fr_ui4_2")),
            # tabPanel(title = "Roadmap (Detailed)", uiOutput("roadmap_Detailed_Fr_ui5")),
            # tabPanel(title = "Roadmap (Detailed)", uiOutput("roadmap_Detailed_Fr_ui6")),
            tabPanel(title = "Roadmap (Detailed)", uiOutput("roadmap_Detailed_Fr_ui7")),
            # tabPanel(title = "Roadmap (Detailed)", uiOutput("roadmap_Detailed_Fr_ui8")),
            # tabPanel(title = "Roadmap (Detailed)", uiOutput("roadmap_Detailed_Fr_ui9")),
            
            # tabPanel(title = "Roadmap Ar", uiOutput("roadmap_Detailed_Ar_ui3")),
            
            tabPanel(title = "Environment",uiOutput("package_status")),
            
            
          )
        ),
        

        
        # ---------------------------------------------------
        # -2️️️ #️⃣- Advanced Exploration tab (unchanged layout) ---
        # ---------------------------------------------------
        
        tabPanel(
          "—> Exploration",
          br(),
          sidebarLayout(
            sidebarPanel(
              width = 2,
              
              numericInput("d_n", label = "d (non-seasonal differencing):", min = 0, value = 0),
              numericInput("DS_n", label = "D (seasonal differencing):", min = 0, value = 0),
              checkboxInput("check_box", HTML("<b>log(S(t))</b>"), value = FALSE),
              
              tags$strong("Stationarity Test (ADF/KPSS)"),
              
              selectInput(
                "adfTypeSt2", "ADF model type:",
                choices = c("none", "drift", "trend"),
                selected = "drift"
              ),
              
              selectInput(
                "alternd2St", "ADF alternative (tseries):",
                choices = c("stationary", "explosive", "regression"),
                selected = "stationary"
              ),
              
              numericInput("LagOrderADFd2St", "Lag order (k):", min = 0, step = 1, value = 10),
              selectInput("alphaSt2", "Significance level (α):",
                          choices = c("0.01", "0.05", "0.1"), selected = "0.05"),
              
              checkboxInput(
                "use_train_explore",
                HTML("<b>Use training split (train only)</b>"),
                value = FALSE
              ),
              helpText("Uncheck to use the full series on this tab.")
            ),
            
            mainPanel(
              width = 10,
              tabsetPanel(
                id = "transform_tabs",
                
                tabPanel("d?D?log?(St) (*)", uiOutput("d_D_Log_ts_Choice_UI")),
                
                tabPanel(
                  "Plot (*)",
                  sidebarLayout(
                    sidebarPanel(
                      width = 2,
                      
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
                      
                      # ✅ NEW: conditional numeric input (shows only for ACF/PACF family)
                      # conditionalPanel(
                      #   condition = "input.plot_type_choice == 'ACF' || input.plot_type_choice == 'PACF' || input.plot_type_choice == 'ACF+PACF' || input.plot_type_choice == 'Time + ACF+PACF'",
                      #   numericInput("St_Lag", "Lag (ACF/PACF)", value = 40, min = 1, step = 1)
                      # ),
                      
                      
                      textInput("plot_width", "Width:", value = "800"),
                      textInput("plot_height", "Height:", value = "500"),
                      
                      selectInput(
                        "plot_theme", "Theme:",
                        choices = c("Gray", "Minimal", "Classic", "Light", "Dark", "BW", "Void"),
                        selected = "Minimal"
                      ),
                      
                      numericInput("St_Lag", "Lag# (ACF/PACF)", value = 40, min = 1, step = 1),
                      numericInput("ma_k", "Moving average window k:", min = 2, step = 1, value = 5),
                      numericInput("lag_m", "Lag plot: number of lags (m):", min = 1, step = 1, value = 12),
                      
                      uiOutput("ts_color_ui")
                    ),
                    mainPanel(
                      width = 9,
                      uiOutput("tsPlot_Choice_UI")
                    )
                  )
                ),
                
                tabPanel("ACF+PACF (*)", uiOutput("difference2ACFPACF_UI")),
                
                tabPanel(
                  "Stationarity [ADF + KPSS]",
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
                
                tabPanel(
                  "CHECKLIST & STEPS",
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
                )
              )
            )
          )
        ),
        
        
        
        #---------------------
        #  1️ 1—Descriptives
        #---------------------
        
        tabPanel(
          "1—Descriptives",
          tabsetPanel(
            tabPanel("What to do", uiOutput("step1_notes")),
            
            # ✅ NEW: Full dataset tab (before Preview)
            tabPanel("Full data", DT::dataTableOutput("full_data_table")),
            
            tabPanel("Preview", tableOutput("data_preview")),
            tabPanel("Descriptives", tableOutput("basic_stats"), plotOutput("hist_plot", height = 300)),
            tabPanel("Missing & outliers", verbatimTextOutput("missing_text"), tableOutput("outlier_table")),
            tabPanel(
              "APA paragraph",
              div(
                style = "width: 100%; height: 1900px; overflow: auto;",
                verbatimTextOutput("apa_data_paragraph")
              )
            )
          )
        ),
        
        
        #------------
        #  2️ —S(t)——
        #------------
        
        tabPanel(
          "2—S(t)——",
          tabsetPanel(
            tabPanel("What to do", uiOutput("step2_notes")),
            
            tabPanel("Time series", plotOutput("plot_series", height = 420)),
            

            #----------------------
            # ---- -S(t) plots——
            #----------------------
            
            tabPanel(
              "S(t) plots",
              sidebarLayout(
                sidebarPanel(
                  width = 3,
                  
                  # Sticky + scrollable sidebar so plots remain visible
                  div(
                    style = "position: sticky; top: 70px; max-height: calc(100vh - 120px); overflow-y: auto; padding-right: 6px;",
                    
                    # =============== Scope ===============
                    wellPanel(
                      tags$b("Scope"),
                      uiOutput("stp_scope_ui"),
                      uiOutput("stp_scope_warning")
                    ),
                    
                    # =============== Plot type ===============
                    wellPanel(
                      tags$b("Plot type"),
                      selectInput(
                        "stp_plot_type",
                        label = NULL,
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
                      
                      # extra params
                      conditionalPanel(
                        "input.stp_plot_type == 'ACF' || input.stp_plot_type == 'PACF' || input.stp_plot_type == 'ACF+PACF' || input.stp_plot_type == 'Time + ACF+PACF'",
                        sliderInput("stp_corr_lag", "Max lag for ACF/PACF", min = 5, max = 200, value = 60, step = 1)
                      ),
                      
                      conditionalPanel(
                        "input.stp_plot_type == 'Moving average'",
                        numericInput("stp_ma_k", "MA window (k)", value = 5, min = 1, step = 1),
                        checkboxInput("stp_ma_show_raw", "Show raw series too", value = TRUE)
                      ),
                      conditionalPanel(
                        "input.stp_plot_type == 'Smoothed (LOESS)'",
                        sliderInput("stp_loess_span", "LOESS span", min = 0.1, max = 1, value = 0.4, step = 0.05)
                      ),
                      conditionalPanel(
                        "input.stp_plot_type == 'Lag plot (1..m)'",
                        numericInput("stp_lag_m", "Number of lags (m)", value = 12, min = 1, step = 1)
                      ),
                      conditionalPanel(
                        "input.stp_plot_type == 'Histogram'",
                        numericInput("stp_hist_bins", "Bins", value = 30, min = 5, step = 1)
                      ),
                      conditionalPanel(
                        "input.stp_plot_type == 'Density'",
                        sliderInput("stp_bw_adj", "Bandwidth adjust", min = 0.2, max = 3, value = 1, step = 0.1)
                      ),
                      conditionalPanel(
                        "input.stp_plot_type == 'Periodogram'",
                        sliderInput("stp_spec_taper", "Taper", min = 0, max = 0.5, value = 0.1, step = 0.05)
                      ),
                      
                      # global transform note (read-only, uses main sidebar inputs)
                      tags$hr(),
                      uiOutput("stp_transform_note")
                    ),
                    
                    # =============== Style (conditional UI) ===============
                    wellPanel(
                      uiOutput("stp_style_ui"),
                      sliderInput("stp_alpha", "Alpha", min = 0, max = 1, value = 1, step = 0.05)
                    ),
                    
                    # =============== Labels size ===============
                    sliderInput("stp_base_size", "Base font size", min = 8, max = 40, value = 12, step = 1),
                    
                    # =============== Output size ===============
                    wellPanel(
                      numericInput("stp_plot_width_px", "Width (px)", value = 980, min = 400, step = 10),
                      numericInput("stp_plot_height_px", "Height (px)", value = 520, min = 250, step = 10)
                    ),
                    
                    
                    # =============== X-axis (Date) ===============
                    wellPanel(
                      tags$b("X-axis (Date)"),
                      checkboxInput("stp_date_smart", "Smart date axis", value = TRUE),
                      
                      # number of ticks (approx.)
                      sliderInput("stp_x_ticks", "Approx. number of x ticks", min = 3, max = 100, value = 8, step = 1),
                      

                      selectInput(
                        "stp_date_format",
                        "Date label format",
                        choices = c(
                          "Auto" = "auto",
                          "Year" = "%Y",
                          "Year-Month" = "%Y-%m",
                          "Month/Year" = "%m/%Y",
                          "Month-Year" = "%b %Y",
                          "Day/Month/Year" = "%d/%m/%Y",
                          "Day-Month-Year" = "%d %b %Y",
                          "ISO" = "%Y-%m-%d",
                          "Custom (type below)" = "custom"
                        ),
                        selected = "auto"
                      ),
                      
                      selectInput(
                        "stp_date_lang",
                        "Date language",
                        choices = c(
                          "Auto / English" = "en",
                          "French" = "fr",
                          "Arabic" = "ar"
                        ),
                        selected = "en"
                      ),
                      
                      
                      conditionalPanel(
                        "input.stp_date_format == 'custom'",
                        textInput("stp_date_format_custom", "Custom format (strftime)", value = "%Y-%m")
                      ),
                      
                      checkboxInput("stp_x_rotate", "Rotate x labels", value = TRUE),
                      sliderInput("stp_x_angle", "Rotation angle", min = -90, max = 90, value = 0, step = 5)
                    ),
                    
                    # =============== Labels ===============
                    wellPanel(
                      textInput("stp_title", "Title", value = ""),
                      textInput("stp_subtitle", "Subtitle", value = ""),
                      textInput("stp_xlab", "X label", value = ""),
                      textInput("stp_ylab", "Y label", value = "")
                    ),
                    
                    # =============== Palette + Theme ===============
                    wellPanel(
                      selectInput(
                        "stp_theme",
                        "Theme",
                        choices = c("Minimal", "Classic", "Light", "Dark", "BW", "Gray", "Void"),
                        selected = "BW"
                      )
                    ),
                    
                    # =============== Palette ===============
                    selectInput(
                      "stp_palette",
                      "Palette",
                      choices = c(
                        "Paired (brewer)",
                        "Viridis"
                      ),
                      selected = "Paired (brewer)"
                    ),
                    checkboxInput("stp_palette_rev", "Reverse palette", value = FALSE),
                    
                    hr(),
                    hr(),
                    tags$b("------------"),
                    hr(),
                  )
                ),
                
                mainPanel(
                  width = 9,
                  uiOutput("stp_plot_ui")
                )
              )
            ),
            
            tabPanel("Seasonality",
                     plotOutput("season_plot", height = 340),
                     plotOutput("subseries_plot", height = 320)
            ),
            tabPanel(
              "Autocorrelation",
              fluidRow(
                column(6, plotOutput("acf_plot", height = 300)),
                column(6, plotOutput("pacf_plot", height = 300))
              )
            ),
            tabPanel("APA paragraph", verbatimTextOutput("apa_explore_paragraph"))
          )
        ),
        
        
        
        
        # ---------------------------------------------------
        #  3️ 3—Decomposition
        # ---------------------------------------------------
        
        tabPanel(
          "3—Decomposition",
          tabsetPanel(
            tabPanel("What to do", uiOutput("step3_notes")),
            tabPanel("Classical additive", plotOutput("decomp_add", height = 420)),
            tabPanel("Classical multiplicative", plotOutput("decomp_mul", height = 420)),
            tabPanel("STL (robust)", plotOutput("stl_plot", height = 520)),
            tabPanel("APA paragraph", verbatimTextOutput("apa_decomp_paragraph"))
          )
        ),
        
        
        # ---------------------------------------------------
        # ️4️ 4—Stationarity
        # ---------------------------------------------------
        
        tabPanel(
          "4—Stationarity",
          sidebarLayout(
            sidebarPanel(
              width = 2,
              h5("Stationarity tests"),
              selectInput("adf_type", "ADF type (ur.df)", choices = c("none", "drift", "trend"), selected = "trend"),
              numericInput("adf_lags", "ADF lags (k)", value = 10, min = 0, step = 1),
              selectInput("kpss_type", "KPSS null", choices = c("mu", "tau"), selected = "mu"),
              
              # actionButton("run_tests", "Run tests", icon = icon("vial")),
              
              tags$div(
                style = "text-align: center;",
                actionButton(
                  "run_tests",
                  label = HTML("Run tests"),
                  icon  = icon("vial"),
                  width = "90%",
                  # class = "blue-btn-1"
                )
              ),
              
              
              
              hr(),
              h5("Differencing preview"),
              numericInput("d_preview", "Non-seasonal difference (d)", value = 0, min = 0, step = 1),
              numericInput("D_preview", "Seasonal difference (D)", value = 0, min = 0, step = 1),
              
              # actionButton(
              #   "preview_diff",
              #   label = HTML("Preview<br>differenced<br>series"),
              #   icon = icon("chart-line")
              # ),
              
              
              tags$div(
                style = "text-align: center;",
                actionButton(
                  "preview_diff",
                  label = HTML("Preview<br>differenced<br>series"),
                  icon  = icon("chart-line"),
                  width = "90%",
                  # class = "blue-btn-1"
                )
              ),
              
            ),
            mainPanel(
              width = 8,
              tabsetPanel(
                tabPanel("What to do", uiOutput("step4_notes")),
                
                # tabPanel("Results", verbatimTextOutput("stationarity_results")),
                tabPanel(
                  "Results",
                  tags$head(tags$style(HTML("
                    #stationarity_results{
                      height: 745px !important;
                      max-height: 800px;
                      width: 120% !important;
                      white-space: pre;
                      overflow-y: auto;
                      border: 2px solid #cccccc;
                      font-size: 13px;
                    }
                  "))),
                  verbatimTextOutput("stationarity_results")
                ),
                
                # tabPanel("Interpretation", verbatimTextOutput("stationarity_interpretation")),
                 tabPanel(
                  "Interpretation",
                  tags$head(tags$style(HTML("
                    #stationarity_interpretation{
                      height: 745px !important;
                      max-height: 800px;
                      width: 120% !important;
                      white-space: pre;
                      overflow-y: auto;
                      border: 2px solid #cccccc;
                      font-size: 13px;
                    }
                  "))),
                  verbatimTextOutput("stationarity_interpretation")
                ),
                
                tabPanel("Suggested diff.", verbatimTextOutput("diff_suggestion")),
                tabPanel("Differenced plot", plotOutput("diff_plot", height = 360)),
                tabPanel("APA paragraph", verbatimTextOutput("apa_stationarity_paragraph"))
              )
            )
          )
        ),
        
        
        # ---------------------------------------------------
        #  ️5️ 5—Auto-ARIMA
        # ---------------------------------------------------
        
        tabPanel(
          "5—Auto-ARIMA",
          sidebarLayout(
            sidebarPanel(
              width = 2,
              checkboxInput("auto_seasonal", "Allow seasonal terms", value = TRUE),
              checkboxInput("auto_stepwise", "Stepwise search", value = TRUE),
              checkboxInput("auto_approx", "Approximation", value = FALSE),
              checkboxInput("auto_allow_mean", "Allow mean/drift", value = TRUE),
              numericInput("auto_max_order", "max.order", value = 10, min = 3, step = 1),
              numericInput("auto_h", "Future horizon h (used when no test set)", value = NA, min = 1, step = 1),
              # actionButton("fit_auto", "Fit Auto-ARIMA", icon = icon("magic")),
              
              
              # actionButton(
              #   "fit_auto",
              #   label = HTML("&nbsp;&nbsp;Fit<br>Auto-ARIMA"),
              #   # label = HTML("Fit<br>Auto-ARIMA"),
              #   icon = icon("magic")
              # ),
              
              tags$div(
                style = "text-align: center;",
                actionButton(
                  "fit_auto",
                  label = HTML("Fit<br>Auto-ARIMA"),
                  icon  = icon("magic"),
                  width = "90%",
                  # class = "blue-btn-1"
                )
              ),
              
              
              
              
              
            ),
            mainPanel(
              width = 8,
              tabsetPanel(
                tabPanel("What to do", uiOutput("step5_notes")),
                
                # tabPanel(
                #   "Model specification",
                #   verbatimTextOutput("auto_model_spec"),
                #   tableOutput("auto_coef_table")
                # ),
                
                tabPanel(
                  "Model equation",
                  withMathJax(
                    tags$div(
                      style = "padding: 8px; line-height: 1.4;",
                      uiOutput("auto_model_equation")
                    )
                  )
                ),
                
                tabPanel(
                  "Diagnostics",
                  fluidRow(
                    column(6, plotOutput("auto_resid_ts",   height = 240)),
                    column(6, plotOutput("auto_resid_acf",  height = 240))
                  ),
                  fluidRow(
                    column(6, plotOutput("auto_resid_hist", height = 240)),
                    column(6, plotOutput("auto_resid_qq",   height = 240))
                  ),
                  ## NEW: Ljung–Box p-values by lag (Auto-ARIMA)
                  fluidRow(
                    column(12, plotOutput("auto_resid_lb_pvals", height = 260))
                  )
                ),
                
                tabPanel(
                  "Residual tests",
                  tags$head(tags$style(HTML("
                    #auto_diag_tests{
                      height: 745px !important;
                      max-height: 800px;
                      width: 100% !important;
                      white-space: pre;
                      overflow-y: auto;
                      border: 2px solid #cccccc;
                      font-size: 13px;
                      padding: 6px;
                      background-color: #fafafa;
                    }
                  "))),
                  verbatimTextOutput("auto_diag_tests")
                ),
                
                tabPanel(
                  "Forecast",
                  verbatimTextOutput("auto_horizon_note"),
                  plotOutput("auto_forecast_plot", height = 420),
                  tableOutput("auto_accuracy_table"),
                  tableOutput("auto_forecast_table")
                ),
                
                tabPanel("APA", verbatimTextOutput("apa_auto_paragraph")),
                
                tabPanel(
                  "Conclusion",
                  tags$head(tags$style(HTML("
                      #auto_conclusion_box{
                        height: 760px !important;
                        max-height: 900px;
                        width: 100% !important;
                        overflow-y: auto;
                        border: 2px solid #cccccc;
                        font-size: 13px;
                        padding: 12px;
                        background-color: #ffffff;
                        line-height: 1.5;
                      }
                    "))),
                  tags$div(
                    id = "auto_conclusion_box",
                    uiOutput("auto_conclusion_full")
                  )
                ),
                
              )
            )
          )
        ),
        
        
        # ---------------------------------------------------
        #  6️ 6—Manual
        # ---------------------------------------------------
        
        tabPanel(
          "6—Manual",
          sidebarLayout(
            sidebarPanel(
              width = 2,
              
              tags$div(
                style = "text-align: center;",
                actionButton(
                  "fit_manual",
                  label = HTML("__Fit__"),
                  icon  = icon("cogs"),
                  width = "90%"
                )
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
              numericInput("manual_h", "Future horizon h (used when no test set)", value = NA, min = 1, step = 1)
            ),
            mainPanel(
              width = 8,
              tabsetPanel(
                # tabPanel("What to do", uiOutput("step6_notes")),
                
                tabPanel(
                  "What to do",
                  tags$div(
                    style = "
                              width: 120%;
                            ",
                    uiOutput("step6_notes")
                  )
                ),
                
                
                tabPanel(
                  "Infos",
                  verbatimTextOutput("manual_split_text"),
                  plotOutput("manual_split_plot", height = 320)
                ),
                
                # tabPanel(
                #   "Model specification",
                #   verbatimTextOutput("manual_model_spec"),
                #   tableOutput("manual_coef_table")
                # ),
                
                tabPanel(
                  "Model equation",
                  withMathJax(uiOutput("manual_model_equation"))
                ),

                tabPanel(
                  "Diagnostics",
                  
                  # Allow page horizontal scroll when container width exceeds screen
                  tags$head(tags$style(HTML("body { overflow-x: auto; }"))),
                  
                  # The whole Diagnostics layout will be created by server so width can be reactive
                  uiOutput("diag_container_ui")
                ),
                
 
                tabPanel(
                  "Residual tests",
                  tags$head(tags$style(HTML("
                    #manual_diag_tests{
                      height: 750px !important;
                      max-height: 800px;
                      width: 120% !important;
                      white-space: pre;
                      overflow-y: auto;
                      border: 2px solid #cccccc;
                      font-size: 14px;
                    }
                  "))),
                  verbatimTextOutput("manual_diag_tests")
                ),
                
                tabPanel(
                  "Forecast",
                  verbatimTextOutput("manual_horizon_note"),
                  plotOutput("manual_forecast_plot", height = 420),
                  tableOutput("manual_accuracy_table"),
                  tableOutput("manual_forecast_table")
                ),
                
                # tabPanel("APA paragraph", verbatimTextOutput("apa_manual_paragraph")),
                
                tabPanel(
                  "Conclusion",
                  tags$head(tags$style(HTML("
                      #manual_conclusion_box{
                        height: 760px !important;
                        max-height: 900px;
                        width: 120% !important;
                        overflow-y: auto;
                        border: 2px solid #cccccc;
                        font-size: 13px;
                        padding: 12px;
                        background-color: #ffffff;
                        line-height: 1.5;
                      }
                    "))),
                  tags$div(
                    id = "manual_conclusion_box",
                    uiOutput("manual_conclusion_full")
                  )
                ),
                
                
              )
            )
          )
        ),
        
        
        # ---------------------------------------------------
        #  7️ 7—Compare
        # ---------------------------------------------------
        
        tabPanel(
          "7—Compare",
          tabsetPanel(
            tabPanel("What to do", uiOutput("step7_notes")),
            tabPanel("Model comparison", tableOutput("comparison_table"), verbatimTextOutput("comparison_interpretation")),
            tabPanel("Methods (APA draft)", verbatimTextOutput("apa_methods_draft")),
            tabPanel("Results (APA draft)", verbatimTextOutput("apa_results_draft")),
            tabPanel("Checklist", uiOutput("paper_checklist_ui"))
          )
        ),
        
        

      )
    )
  )
)


