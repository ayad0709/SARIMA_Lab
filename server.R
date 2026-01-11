# server.R — SARIMA Scientific Writing Lab (v2)
# FIX v2: Manual SARIMA validation forecast aligns with test horizon (h = test length).

library(shiny)
library(ggplot2)
library(forecast)
library(lubridate)
library(zoo)
library(tseries)
library(urca)


has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)






# ============================================================
# --- MOD: Robust date parsing (R Date / POSIX / Excel serial / text) ---
# ============================================================
parse_dates <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  
  if (is.numeric(x)) {
    # MOD: heuristic to distinguish R Date numeric vs Excel serial
    # - R Date numeric typically around [-60000, 60000] (days since 1970-01-01)
    # - Excel serial is typically large positive (days since 1899-12-30)
    med <- suppressWarnings(stats::median(x, na.rm = TRUE))
    if (is.finite(med) && med > -60000 && med < 60000) {
      return(as.Date(x, origin = "1970-01-01"))
    }
    return(as.Date(x, origin = "1899-12-30"))
  }
  
  x_chr <- as.character(x)
  
  d <- suppressWarnings(as.Date(zoo::as.yearmon(x_chr)))
  if (all(is.na(d))) {
    d <- suppressWarnings(lubridate::parse_date_time(
      x_chr,
      orders = c(
        "ymd", "dmy", "mdy", "Ymd", "Y-m-d", "d-m-Y", "m/d/Y",
        "Y", "ym", "my", "bY", "Y-b", "any"
      )
    ))
    d <- as.Date(d)
  }
  
  d
}




# ---------------- APA helpers ----------------

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < .001")
  s <- format(p, digits = 3, nsmall = 3)
  s <- sub("^0\\.", ".", s)
  paste0("p = ", s)
}

fmt_num <- function(x, digits = 2) {
  if (is.na(x)) return("NA")
  format(round(x, digits), nsmall = digits, trim = TRUE)
}

fmt_pct <- function(x, digits = 1) {
  if (is.na(x)) return("NA")
  paste0(fmt_num(100 * x, digits), "%")
}

# ---------------- Dates & time grid ----------------

parse_dates <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  x_chr <- as.character(x)

  d <- suppressWarnings(as.Date(zoo::as.yearmon(x_chr)))
  if (all(is.na(d))) {
    d <- suppressWarnings(lubridate::parse_date_time(
      x_chr,
      orders = c("ymd", "dmy", "mdy", "Ymd", "Y-m-d", "d-m-Y", "m/d/Y", "Y", "ym", "my", "bY", "Y-b", "any")
    ))
    d <- as.Date(d)
  }
  d
}

freq_value <- function(input) {
  if (identical(input$frequency, "other")) as.numeric(input$customFrequency) else as.numeric(input$frequency)
}

freq_to_by <- function(freq) {
  switch(
    as.character(freq),
    "12" = "month",
    "4" = "quarter",
    "52" = "week",
    "365" = "day",
    "7" = "day",
    NULL
  )
}

make_regular_grid <- function(dates, by) {
  if (is.null(by)) return(sort(unique(dates)))
  seq.Date(from = min(dates), to = max(dates), by = by)
}

extend_grid <- function(last_x, h, by) {
  if (inherits(last_x, "Date") && !is.null(by)) {
    seq.Date(from = last_x, by = by, length.out = h + 1)[-1]
  } else {
    (as.numeric(last_x) + seq_len(h))
  }
}

# ---------------- Missing values & transforms ----------------

fill_missing <- function(y, policy, freq) {
  if (policy == "drop") return(y)
  if (policy == "locf") {
    y <- zoo::na.locf(y, na.rm = FALSE)
    y <- zoo::na.locf(y, fromLast = TRUE, na.rm = FALSE)
    return(y)
  }
  if (policy == "linear") return(zoo::na.approx(y, na.rm = FALSE))
  if (policy == "seasonal") return(as.numeric(forecast::na.interp(ts(y, frequency = freq))))
  y
}

apply_transform <- function(y, transform, lambda) {
  if (transform == "none") return(y)
  if (any(y <= 0, na.rm = TRUE)) stop("Chosen transformation requires strictly positive values.")
  if (transform == "log") return(log(y))
  if (transform == "boxcox") {
    lam <- if (!is.na(lambda)) lambda else forecast::BoxCox.lambda(y, lower = 0)
    return(forecast::BoxCox(y, lam))
  }
  y
}

# ---------------- Descriptives & diagnostics ----------------

basic_stats_df <- function(y) {
  y <- as.numeric(y)
  y_ok <- y[is.finite(y)]
  n <- length(y)
  n_ok <- length(y_ok)
  miss <- n - n_ok
  if (n_ok == 0) {
    return(data.frame(Metric = c("N", "Missing"), Value = c(n, miss), stringsAsFactors = FALSE))
  }
  mu <- mean(y_ok)
  sdv <- sd(y_ok)
  med <- median(y_ok)
  mn <- min(y_ok); mx <- max(y_ok)
  cv <- if (mu != 0) sdv / mu else NA_real_
  m3 <- mean((y_ok - mu)^3)
  m4 <- mean((y_ok - mu)^4)
  skew <- if (sdv > 0) m3 / (sdv^3) else NA_real_
  kurt <- if (sdv > 0) m4 / (sdv^4) else NA_real_
  data.frame(
    Metric = c("N", "Valid", "Missing", "Mean", "Median", "SD", "Min", "Max", "CV", "Skewness", "Kurtosis"),
    Value = c(n, n_ok, miss, mu, med, sdv, mn, mx, cv, skew, kurt),
    stringsAsFactors = FALSE
  )
}

z_outliers <- function(y, z_cut = 3) {
  y <- as.numeric(y)
  ok <- is.finite(y)
  mu <- mean(y[ok])
  sdv <- sd(y[ok])
  if (!is.finite(sdv) || sdv == 0) return(integer(0))
  which(abs((y - mu) / sdv) >= z_cut)
}

coef_table <- function(fit) {
  s <- tryCatch(summary(fit), error = function(e) NULL)
  if (is.null(s) || is.null(s$coef)) return(data.frame())
  cf <- as.data.frame(s$coef)
  cf$term <- rownames(cf)
  rownames(cf) <- NULL
  names(cf) <- sub("Pr\\(>\\|t\\|\\)", "p_value", names(cf))
  names(cf) <- sub("Pr\\(>\\|z\\|\\)", "p_value", names(cf))
  cf
}

diag_tests_text <- function(resid, lag, fitdf = 0) {
  r <- as.numeric(na.omit(resid))
  out <- c("Residual diagnostics")

  lb <- tryCatch(Box.test(r, lag = lag, type = "Ljung-Box", fitdf = fitdf), error = function(e) NULL)
  if (!is.null(lb)) out <- c(out, paste0("- Ljung-Box: Q(", lag, ") = ", fmt_num(lb$statistic, 2), ", ", fmt_p(lb$p.value)))

  jb <- tryCatch(tseries::jarque.bera.test(r), error = function(e) NULL)
  if (!is.null(jb)) out <- c(out, paste0("- Jarque–Bera: JB = ", fmt_num(jb$statistic, 2), ", ", fmt_p(jb$p.value)))

  sh <- tryCatch(stats::shapiro.test(if (length(r) > 5000) sample(r, 5000) else r), error = function(e) NULL)
  if (!is.null(sh)) out <- c(out, paste0("- Shapiro–Wilk: W = ", fmt_num(sh$statistic, 3), ", ", fmt_p(sh$p.value)))

  if (has_pkg("nortest")) {
    ad <- tryCatch(nortest::ad.test(r), error = function(e) NULL)
    if (!is.null(ad)) out <- c(out, paste0("- Anderson–Darling: A = ", fmt_num(ad$statistic, 2), ", ", fmt_p(ad$p.value)))
  }

  if (has_pkg("FinTS")) {
    arch <- tryCatch(FinTS::ArchTest(r, lags = min(12, max(1, floor(lag / 2)))), error = function(e) NULL)
    if (!is.null(arch)) out <- c(out, paste0("- ARCH LM: TR^2 = ", fmt_num(arch$statistic, 2), ", ", fmt_p(arch$p.value)))
  }

  out <- c(out, "", "Interpretation (typical): Ljung-Box p > .05 suggests residuals are approximately white noise.")
  paste(out, collapse = "\n")
}

accuracy_df <- function(actual, forecast_mean) {
  a <- as.numeric(actual)
  f <- as.numeric(forecast_mean)
  n <- min(length(a), length(f))
  a <- a[seq_len(n)]
  f <- f[seq_len(n)]
  e <- a - f

  rmse <- sqrt(mean(e^2, na.rm = TRUE))
  mae <- mean(abs(e), na.rm = TRUE)
  mape <- mean(abs(e / a), na.rm = TRUE)
  smape <- mean(2 * abs(e) / (abs(a) + abs(f)), na.rm = TRUE)

  data.frame(
    Metric = c("RMSE", "MAE", "MAPE", "sMAPE"),
    Value = c(rmse, mae, mape, smape),
    stringsAsFactors = FALSE
  )
}

forecast_table <- function(fc) {
  out <- data.frame(step = seq_along(fc$mean), mean = as.numeric(fc$mean))
  if (!is.null(fc$lower)) {
    out$lo80 <- as.numeric(fc$lower[, 1])
    out$hi80 <- as.numeric(fc$upper[, 1])
    if (ncol(fc$lower) >= 2) {
      out$lo95 <- as.numeric(fc$lower[, 2])
      out$hi95 <- as.numeric(fc$upper[, 2])
    }
  }
  out
}

# ---------------- Plot helpers (Date-safe) ----------------

plot_series_df <- function(df, train_n) {
  df$set <- ifelse(seq_len(nrow(df)) <= train_n, "Train", "Test/Future")
  df
}

plot_forecast_df <- function(obs_df, train_n, fc, by) {
  n_obs <- nrow(obs_df)
  test_n <- n_obs - train_n
  h <- length(fc$mean)

  # FIX v2:
  # If a test set exists, align the forecast x values with the test period.
  # Only extend beyond the observed sample if h > test length.
  if (test_n > 0) {
    x_test <- obs_df$x[(train_n + 1):n_obs]
    n_align <- min(h, length(x_test))
    x_fc <- x_test[seq_len(n_align)]
    if (h > n_align) {
      x_extra <- extend_grid(obs_df$x[n_obs], h - n_align, by)
      x_fc <- c(x_fc, x_extra)
    }
  } else {
    x_fc <- extend_grid(obs_df$x[n_obs], h, by)
  }

  fc_tab <- forecast_table(fc)
  fc_tab$x <- x_fc
  fc_tab
}

gg_forecast_plot <- function(obs_df, train_n, fc_df, title) {
  p <- ggplot() +
    geom_line(data = obs_df[seq_len(train_n), ], aes(x = x, y = y, color = "Train")) +
    theme_minimal() +
    labs(title = title, x = "Time", y = "Value", color = NULL) +
    theme(legend.position = "bottom")

  if (nrow(obs_df) > train_n) {
    p <- p + geom_line(data = obs_df[(train_n + 1):nrow(obs_df), ], aes(x = x, y = y, color = "Test"))
  }

  if ("lo95" %in% names(fc_df)) p <- p + geom_ribbon(data = fc_df, aes(x = x, ymin = lo95, ymax = hi95), alpha = 0.15)
  if ("lo80" %in% names(fc_df)) p <- p + geom_ribbon(data = fc_df, aes(x = x, ymin = lo80, ymax = hi80), alpha = 0.25)

  p + geom_line(data = fc_df, aes(x = x, y = mean, color = "Forecast"), linewidth = 1)
}

# ---------------- Shiny server ----------------

server <- function(input, output, session) {

  # ---- Roadmap & teaching notes ----

  output$roadmap_ui <- renderUI({
    tags$div(
      style = "background:#f7f7f7;padding:12px;border-radius:8px;",
      tags$h4("Roadmap (what students do, what they write)"),
      tags$ol(
        tags$li(tags$b("Describe the data"), ": sample size, missing values, descriptive statistics."),
        tags$li(tags$b("Explore visually"), ": trend/seasonality/outliers; report observations."),
        tags$li(tags$b("Decompose"), ": justify additive vs multiplicative; use STL when robust needed."),
        tags$li(tags$b("Check stationarity"), ": ADF/KPSS/PP; justify differencing (d and D)."),
        tags$li(tags$b("Fit a baseline model"), ": Auto-ARIMA to obtain a strong starting SARIMA."),
        tags$li(tags$b("Fit a theory-driven model"), ": Manual SARIMA using ACF/PACF + tests."),
        tags$li(tags$b("Diagnose & compare"), ": residual tests + forecast accuracy; choose final model."),
        tags$li(tags$b("Write your paper"), ": use APA paragraphs in each step; assemble Methods/Results.")
      )
    )
  })
  

  note_box <- function(items) {
    if (!isTRUE(input$show_teaching_notes)) return(NULL)
    tags$div(
      style = "background:#eef5ff;padding:10px;border-radius:6px;margin-bottom:10px;",
      tags$b("What to do:"),
      tags$ul(lapply(items, tags$li))
    )
  }

  output$step1_notes <- renderUI({ note_box(list(
    "Check the data preview; confirm correct date/value columns.",
    "Describe missingness and how you handled it.",
    "Report mean/SD and distribution shape (skewness/kurtosis)."
  ))})

  output$step2_notes <- renderUI({ note_box(list(
    "Inspect trend and seasonality in the time series plot.",
    "Use seasonal/subseries plots to justify seasonality.",
    "Use ACF/PACF to propose candidate SARIMA orders."
  ))})

  output$step3_notes <- renderUI({ note_box(list(
    "Compare additive vs multiplicative decomposition.",
    "Use STL (robust) if you suspect outliers or changing seasonality.",
    "Interpret trend/seasonal/remainder components."
  ))})

  output$step4_notes <- renderUI({ note_box(list(
    "Run stationarity tests (ADF/KPSS/PP).",
    "Use suggested differencing as hints, not rules.",
    "Preview differenced series to confirm stationarity visually."
  ))})

  output$step5_notes <- renderUI({ note_box(list(
    "Fit Auto-ARIMA as a baseline.",
    "Report selected (p,d,q)(P,D,Q)[s] and key diagnostics.",
    "Evaluate forecast accuracy on the test set (if available)."
  ))})

  output$step6_notes <- renderUI({ note_box(list(
    "Select orders guided by ACF/PACF and differencing evidence.",
    "Check residual diagnostics: white-noise residuals are expected.",
    "Compare with Auto-ARIMA to justify your final choice.",
    "IMPORTANT: With a test set, the validation forecast is forced to h = test length (overlays the test period)."
  ))})

  output$step7_notes <- renderUI({ note_box(list(
    "Compare models using AICc/BIC and test-set accuracy (if available).",
    "Use Methods/Results drafts as a starting point; edit for your dataset.",
    "Use the Checklist to avoid common reporting omissions."
  ))})

  output$paper_checklist_ui <- renderUI({
    tags$div(
      tags$h4("Paper checklist (SARIMA reporting essentials)"),
      tags$ul(
        tags$li("Data: variable, unit, frequency, time range, N."),
        tags$li("Missing values: amount and handling approach."),
        tags$li("Exploration: plots + seasonality indicators + ACF/PACF rationale."),
        tags$li("Stationarity: ADF/KPSS/PP results and differencing decisions."),
        tags$li("Model: (p,d,q)(P,D,Q)[s], software, estimation method."),
        tags$li("Diagnostics: Ljung-Box (+ normality/ARCH if reported)."),
        tags$li("Forecast evaluation: horizon, metrics, interpretation.")
      )
    )
  })

  # ---- Data ingest ----

  raw_data <- reactive({
    req(input$fileData)
    ext <- tolower(tools::file_ext(input$fileData$name))
    if (ext == "csv") {
      read.csv(input$fileData$datapath, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (ext %in% c("xls", "xlsx")) {
      validate(need(has_pkg("readxl"), "Install package 'readxl' to read Excel files."))
      readxl::read_excel(input$fileData$datapath)
    } else {
      validate("Unsupported file type. Use CSV or XLSX.")
    }
  })

  output$dateColUI <- renderUI({
    req(raw_data())
    cols <- names(raw_data())
    selectInput("dateCol", "Date column", choices = cols, selected = cols[1])
  })

  output$valueColUI <- renderUI({
    req(raw_data())
    cols <- names(raw_data())
    selectInput("valueCol", "Value column", choices = cols, selected = cols[min(2, length(cols))])
  })

  
  # ============================================================
  # --- MOD: Use robust date parsing inside prepared() ---
  # ============================================================
  prepared <- reactive({
    req(raw_data(), input$dateCol, input$valueCol)
    
    f <- freq_value(input)
    by <- freq_to_by(f)
    
    df <- raw_data()
    
    # MOD: robust conversion to Date
    d <- parse_dates(df[[input$dateCol]])
    
    y <- suppressWarnings(as.numeric(df[[input$valueCol]]))
    
    keep <- !is.na(d)
    df2 <- data.frame(date = as.Date(d[keep]), y_raw = y[keep])
    df2 <- df2[order(df2$date), , drop = FALSE]
    
    if (isTRUE(input$align_regular) && !is.null(by)) {
      grid <- make_regular_grid(df2$date, by = by)
      df2 <- merge(data.frame(date = grid), df2, by = "date", all.x = TRUE, sort = TRUE)
    }
    
    df2$y_filled <- fill_missing(df2$y_raw, input$missing_policy, f)
    
    df2$y_trans <- tryCatch(
      apply_transform(df2$y_filled, input$transform, input$lambda),
      error = function(e) { validate(e$message); df2$y_filled }
    )
    
    # MOD: If we have a Date-based grid, use Date on x
    if (!is.null(by)) {
      df2$x <- df2$date
      x_label <- "Date"
    } else {
      df2$x <- seq_len(nrow(df2))
      x_label <- "Index"
    }
    
    list(df = df2, freq = f, by = by, x_label = x_label)
  })
  
  

  # ============================================================
  # --- MOD: Fix Train split = 100% (no test set) bug in ts_train_test() ---
  # Reason: ts(numeric(0)) is invalid in R ("ts must have at least one observation")
  # Fix: store ts_test = NULL when there is no test set, and guard all downstream uses.
  # ============================================================
  ts_train_test <- reactive({
    p <- prepared()
    df <- p$df
    
    ok <- is.finite(df$y_trans)
    dfm <- df[ok, , drop = FALSE]
    validate(need(nrow(dfm) >= 10, "Not enough valid observations after cleaning."))
    
    train_n <- max(2, floor(nrow(dfm) * as.numeric(input$train_prop)))
    train_n <- min(train_n, nrow(dfm)) # MOD: hard cap so train_n never exceeds n
    
    y_tr <- dfm$y_trans[seq_len(train_n)]
    y_te <- if (train_n < nrow(dfm)) dfm$y_trans[(train_n + 1):nrow(dfm)] else numeric(0)
    
    list(
      dfm = dfm,
      train_n = train_n,
      test_n = length(y_te),
      ts_train = ts(y_tr, start = 1, frequency = p$freq),
      
      # --- MOD: when test is empty, keep it NULL (do NOT create ts(numeric(0))) ---
      ts_test = if (length(y_te) > 0) {
        ts(y_te, start = train_n + 1, frequency = p$freq)
      } else {
        NULL
      }
    )
  })
  
  # ============================================================
  # --- MOD: Helper to safely reconstruct the full series (train + test) ---
  # Used in plots requiring the full observed sample.
  # ============================================================
  full_ts <- function(s) {
    if (!is.null(s$ts_test) && length(s$ts_test) > 0) {
      ts(
        c(as.numeric(s$ts_train), as.numeric(s$ts_test)),
        start = 1,
        frequency = frequency(s$ts_train)
      )
    } else {
      s$ts_train
    }
  }
  
  # ============================================================
  # --- MOD: Replace all "ts(c(train, test))" calls by full_ts(s) ---
  # This prevents errors when train_prop = 1 (test set absent).
  # ============================================================
  
  output$season_plot <- renderPlot({
    req(ts_train_test())
    s <- ts_train_test()
    x <- full_ts(s) # MOD
    validate(need(frequency(x) >= 2, "Seasonal plots need frequency >= 2."))
    forecast::seasonplot(x, s = frequency(x))
  })
  
  output$subseries_plot <- renderPlot({
    req(ts_train_test())
    s <- ts_train_test()
    x <- full_ts(s) # MOD
    validate(need(frequency(x) >= 2, "Subseries plot needs frequency >= 2."))
    forecast::ggsubseriesplot(x) + theme_minimal()
  })
  
  output$decomp_add <- renderPlot({
    req(ts_train_test())
    s <- ts_train_test()
    x <- full_ts(s) # MOD
    validate(need(frequency(x) >= 2, "Decomposition needs frequency >= 2."))
    validate(need(length(x) >= 2 * frequency(x), "Need at least 2 seasonal cycles."))
    plot(decompose(x, type = "additive"))
  })
  
  output$decomp_mult <- renderPlot({
    req(ts_train_test())
    s <- ts_train_test()
    x <- full_ts(s) # MOD
    validate(need(frequency(x) >= 2, "Decomposition needs frequency >= 2."))
    validate(need(length(x) >= 2 * frequency(x), "Need at least 2 seasonal cycles."))
    validate(need(all(x > 0, na.rm = TRUE), "Multiplicative decomposition requires strictly positive values."))
    plot(decompose(x, type = "multiplicative"))
  })
  
  output$diff_suggestion <- renderPrint({
    req(ts_train_test())
    s <- ts_train_test()
    x <- full_ts(s) # MOD
    d_rec <- tryCatch(forecast::ndiffs(x), error = function(e) NA_integer_)
    D_rec <- tryCatch(forecast::nsdiffs(x), error = function(e) NA_integer_)
    cat("Suggested differencing (heuristics):\n")
    cat("- ndiffs (d):", d_rec, "\n")
    cat("- nsdiffs (D):", D_rec, "\n")
  })
  
  
  
  # ---- Step 1 outputs ----
  
  
  # ============================================================
  # --- MOD: display Date columns as readable strings in the preview table ---
  # ============================================================
  output$data_preview <- renderTable({
    req(prepared())
    df <- head(prepared()$df, 12)
    
    # MOD: ensure Date columns print nicely (Shiny sometimes prints Date as numeric)
    for (nm in intersect(c("date", "x"), names(df))) {
      if (inherits(df[[nm]], "Date")) df[[nm]] <- format(df[[nm]], "%Y-%m-%d")
    }
    
    df
  }, rownames = FALSE)
  
  output$basic_stats <- renderTable({
    req(prepared())
    basic_stats_df(prepared()$df$y_filled)
  }, rownames = FALSE)
  
  output$hist_plot <- renderPlot({
    req(prepared())
    df <- prepared()$df
    
    ggplot(df, aes(x = y_filled)) +
      geom_histogram(bins = 30) +
      theme_minimal() +
      labs(title = "Distribution (filled values)", x = "Value", y = "Count")
  })
  

  # output$data_preview <- renderTable({ req(prepared()); head(prepared()$df, 12) }, rownames = FALSE)
  
  
  output$basic_stats <- renderTable({ req(prepared()); basic_stats_df(prepared()$df$y_filled) }, rownames = FALSE)

  output$hist_plot <- renderPlot({
    req(prepared())
    df <- prepared()$df
    ggplot(df, aes(x = y_filled)) +
      geom_histogram(bins = 30) +
      theme_minimal() +
      labs(title = "Distribution (filled values)", x = "Value", y = "Count")
  })

  output$missing_text <- renderPrint({
    req(prepared())
    p <- prepared()
    df <- p$df
    n <- nrow(df)
    miss_raw <- sum(is.na(df$y_raw))
    cat("Missing values\n")
    cat("- N rows:", n, "\n")
    cat("- Missing in raw y:", miss_raw, " (", fmt_pct(miss_raw / n), ")\n", sep = "")
    cat("- Handling method:", input$missing_policy, "\n")
    cat("- Date range:", format(min(df$date)), "to", format(max(df$date)), "\n")
  })

  output$outlier_table <- renderTable({
    req(prepared())
    df <- prepared()$df
    idx <- z_outliers(df$y_filled, z_cut = 3)
    if (length(idx) == 0) return(data.frame(message = "No |z| ≥ 3 outliers detected (filled values)."))
    data.frame(index = idx, date = df$date[idx], value = df$y_filled[idx], stringsAsFactors = FALSE)
  }, rownames = FALSE)

  output$apa_data_paragraph <- renderPrint({
    req(prepared())
    p <- prepared()
    df <- p$df
    n <- nrow(df)
    miss_raw <- sum(is.na(df$y_raw))
    miss_pct <- miss_raw / n
    trans <- switch(input$transform, none = "no transformation", log = "a log transformation", boxcox = "a Box–Cox transformation")
    cat(
      "APA-ready paragraph (edit variable names/unit as needed):\n\n",
      "The dataset comprised ", n, " observations collected from ", format(min(df$date)), " to ", format(max(df$date)), ". ",
      "Missing values were present in ", miss_raw, " observations (", fmt_pct(miss_pct), "). ",
      "Missingness was handled using the ", input$missing_policy, " approach. ",
      "Prior to modeling, ", trans, " was applied to the series when appropriate.\n",
      sep = ""
    )
  })

  # ---- Step 2 outputs ----

  output$plot_series <- renderPlot({
    req(prepared(), ts_train_test())
    p <- prepared()
    s <- ts_train_test()
    df <- s$dfm
    df$set <- ifelse(seq_len(nrow(df)) <= s$train_n, "Train", "Test/Future")
    ggplot(df, aes(x = x, y = y_trans, color = set)) +
      geom_line(linewidth = 0.9) +
      theme_minimal() +
      labs(title = "Time series (transformed)", x = p$x_label, y = "Value", color = NULL) +
      theme(legend.position = "bottom")
  })

  output$season_plot <- renderPlot({
    req(ts_train_test())
    s <- ts_train_test()
    x <- ts(c(as.numeric(s$ts_train), as.numeric(s$ts_test)), start = 1, frequency = frequency(s$ts_train))
    validate(need(frequency(x) >= 2, "Seasonal plots need frequency >= 2."))
    forecast::seasonplot(x, s = frequency(x))
  })

  output$subseries_plot <- renderPlot({
    req(ts_train_test())
    s <- ts_train_test()
    x <- ts(c(as.numeric(s$ts_train), as.numeric(s$ts_test)), start = 1, frequency = frequency(s$ts_train))
    validate(need(frequency(x) >= 2, "Subseries plot needs frequency >= 2."))
    forecast::ggsubseriesplot(x) + theme_minimal()
  })

  output$acf_plot <- renderPlot({ req(ts_train_test()); x <- ts_train_test()$ts_train; plot(acf(x, lag.max = min(60, length(x) - 1)), main = "ACF (training)") })
  output$pacf_plot <- renderPlot({ req(ts_train_test()); x <- ts_train_test()$ts_train; plot(pacf(x, lag.max = min(60, length(x) - 1)), main = "PACF (training)") })

  output$apa_explore_paragraph <- renderPrint({
    req(prepared())
    p <- prepared()
    cat(
      "APA-ready paragraph (edit based on what you observed):\n\n",
      "Visual inspection of the time series suggested the presence of trend and/or seasonal patterns. ",
      "Seasonal and subseries plots were examined to evaluate recurring periodic behavior (s = ", p$freq, "). ",
      "Autocorrelation (ACF) and partial autocorrelation (PACF) plots were inspected to inform candidate ARIMA and seasonal ARIMA orders.\n",
      sep = ""
    )
  })

  # ---- Step 3 outputs ----

  output$decomp_add <- renderPlot({
    req(ts_train_test())
    s <- ts_train_test()
    x <- ts(c(as.numeric(s$ts_train), as.numeric(s$ts_test)), start = 1, frequency = frequency(s$ts_train))
    validate(need(frequency(x) >= 2, "Decomposition needs frequency >= 2."))
    validate(need(length(x) >= 2 * frequency(x), "Need at least 2 seasonal cycles."))
    plot(decompose(x, type = "additive"))
  })

  output$decomp_mul <- renderPlot({
    req(ts_train_test())
    s <- ts_train_test()
    x <- ts(c(as.numeric(s$ts_train), as.numeric(s$ts_test)), start = 1, frequency = frequency(s$ts_train))
    validate(need(frequency(x) >= 2, "Decomposition needs frequency >= 2."))
    validate(need(length(x) >= 2 * frequency(x), "Need at least 2 seasonal cycles."))
    validate(need(all(as.numeric(x) > 0), "Multiplicative decomposition requires strictly positive values."))
    plot(decompose(x, type = "multiplicative"))
  })

  output$stl_plot <- renderPlot({
    req(ts_train_test())
    s <- ts_train_test()
    x <- ts(c(as.numeric(s$ts_train), as.numeric(s$ts_test)), start = 1, frequency = frequency(s$ts_train))
    validate(need(frequency(x) >= 2, "STL needs frequency >= 2."))
    validate(need(length(x) >= 2 * frequency(x), "Need at least 2 seasonal cycles."))
    fit <- stl(x, s.window = "periodic", robust = TRUE)
    plot(fit, main = "STL decomposition (robust)")
  })

  output$apa_decomp_paragraph <- renderPrint({
    cat(
      "APA-ready paragraph (edit based on which decomposition you used):\n\n",
      "Decomposition methods were used to separate the series into trend, seasonal, and remainder components. ",
      "Classical decomposition was evaluated under additive and multiplicative assumptions, and robust STL decomposition ",
      "was also examined to reduce the influence of potential outliers. ",
      "The decomposition results provided evidence regarding the stability of seasonality and the magnitude of irregular variation.\n",
      sep = ""
    )
  })

  # ---- Step 4 outputs ----

  output$diff_suggestion <- renderPrint({
    req(ts_train_test())
    s <- ts_train_test()
    x <- ts(c(as.numeric(s$ts_train), as.numeric(s$ts_test)), start = 1, frequency = frequency(s$ts_train))
    d_rec <- tryCatch(forecast::ndiffs(x), error = function(e) NA_integer_)
    D_rec <- tryCatch(forecast::nsdiffs(x), error = function(e) NA_integer_)
    cat("Suggested differencing (heuristics):\n")
    cat("- ndiffs (d):", d_rec, "\n")
    cat("- nsdiffs (D):", D_rec, "\n")
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  stationarity <- eventReactive(input$run_tests, {
    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
    to_int <- function(x, d = 0L) { v <- suppressWarnings(as.integer(x)); if (!is.finite(v)) d else v }
    
    s <- ts_train_test()
    x <- as.numeric(stats::na.omit(s$ts_train))
    
    k <- max(0L, to_int(input$adf_lags, 10L))
    
    # single source of truth for the deterministic type
    type_used <- tolower(as.character(input$adf_type %||% input$adfTypeSt2 %||% "drift"))
    type_used <- if (type_used %in% c("none","drift","trend")) type_used else "drift"
    
    alpha_used <- suppressWarnings(as.numeric(input$alphaSt2 %||% 0.05))
    if (!is.finite(alpha_used)) alpha_used <- 0.05
    
    out <- list(type_used = type_used, alpha_used = alpha_used)
    
    # run tests (keep errors if any)
    out$adf <- tryCatch(tseries::adf.test(x, k = k, alternative = input$alternd2St %||% "stationary"),
                        error = function(e) { out$adf_error <<- e$message; NULL })
    out$kpss <- tryCatch({
      null_kpss <- if (identical(type_used, "trend")) "Trend" else "Level"
      tseries::kpss.test(x, null = null_kpss)
    }, error = function(e) { out$kpss_error <<- e$message; NULL })
    out$pp <- tryCatch(tseries::pp.test(x),
                       error = function(e) { out$pp_error <<- e$message; NULL })
    out$ur <- tryCatch(urca::ur.df(x, type = type_used, lags = k),
                       error = function(e) { out$ur_error <<- e$message; NULL })
    
    out
  })
  
  alpha_to_col <- function(a) {
    if (!is.finite(a)) return("5pct")
    if (abs(a - 0.01) < 1e-8) return("1pct")
    if (abs(a - 0.05) < 1e-8) return("5pct")
    if (abs(a - 0.10) < 1e-8 || abs(a - 0.1) < 1e-8) return("10pct")
    "5pct"
  }
  
  extract_tau_urdf <- function(ur_obj, alpha_used = 0.05, type_used = NULL) {
    out <- list(tau_obs = NA_real_, tau_crit = NA_real_, tau_row = NA_character_, alpha_col = NA_character_)
    if (is.null(ur_obj)) return(out)
    
    # tau observed: first entry whose name starts with "tau" (fallback: first element)
    ts_vec <- tryCatch(ur_obj@teststat, error = function(e) NULL)
    if (!is.null(ts_vec)) {
      nms <- names(ts_vec)
      if (!is.null(nms)) {
        idx <- grep("^tau", nms)
        out$tau_obs <- suppressWarnings(as.numeric(ts_vec[ if (length(idx)) idx[1] else 1 ]))
      } else {
        out$tau_obs <- suppressWarnings(as.numeric(ts_vec[1]))
      }
    }
    
    cv <- tryCatch(ur_obj@cval, error = function(e) NULL)
    if (is.null(cv)) return(out)
    
    # Work with matrix cval (usual case)
    if (is.matrix(cv)) {
      rnames <- rownames(cv); cnames <- colnames(cv)
      
      # preferred tau row from declared type
      tau_pref <- switch(tolower(type_used %||% ""), "none" = "tau1", "drift" = "tau2", "trend" = "tau3", NA_character_)
      tau_rows <- grep("^tau", rnames, value = TRUE)
      
      row_pick <- if (!is.na(tau_pref) && tau_pref %in% rnames) tau_pref else if (length(tau_rows)) tau_rows[1] else NA_character_
      out$tau_row <- row_pick
      
      col_pick <- alpha_to_col(alpha_used)
      if (!is.null(cnames) && !(col_pick %in% cnames)) {
        # choose the nearest available alpha column
        pct <- suppressWarnings(as.numeric(gsub("pct", "", cnames)))/100
        if (any(is.finite(pct))) col_pick <- cnames[ which.min(abs(pct - alpha_used)) ] else col_pick <- cnames[1]
      }
      out$alpha_col <- col_pick
      
      if (!is.na(row_pick) && !is.null(col_pick) && row_pick %in% rnames && col_pick %in% cnames) {
        out$tau_crit <- suppressWarnings(as.numeric(cv[row_pick, col_pick]))
      }
      return(out)
    }
    
    # Named vector fallback (rare)
    col_pick <- alpha_to_col(alpha_used)
    out$alpha_col <- col_pick
    if (!is.null(names(cv)) && col_pick %in% names(cv)) {
      out$tau_crit <- suppressWarnings(as.numeric(cv[[col_pick]]))
    }
    out
  }
  
  
  
  # stationarity <- eventReactive(input$run_tests, {
  #   # --- helpers ---
  #   `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
  #   to_int <- function(x, d = 0L) {
  #     v <- suppressWarnings(as.integer(x)); if (!is.finite(v)) d else v
  #   }
  #   # Map many possible UI values to tseries KPSS 'null' values
  #   norm_kpss_null <- function(val, adf_type = NULL) {
  #     v <- tolower(as.character(val))
  #     if (v %in% c("level", "mu", "l", "lev"))  return("Level")
  #     if (v %in% c("trend", "tau", "t"))        return("Trend")
  #     # fall back: if ADF type has trend, use Trend; else Level
  #     if (!is.null(adf_type) && tolower(adf_type) == "trend") return("Trend")
  #     "Level"
  #   }
  #   
  #   # --- get training series ---
  #   s <- ts_train_test()
  #   x_raw <- s$ts_train
  #   
  #   # coerce to numeric, drop NA/Inf
  #   x <- suppressWarnings(as.numeric(x_raw))
  #   x <- x[is.finite(x)]
  #   
  #   # metadata & guards
  #   N  <- length(x)
  #   sdv <- if (N > 1) stats::sd(x) else NA_real_
  #   k  <- max(0L, to_int(input$adf_lags, 10L))
  #   
  #   adf_type_in  <- tolower(as.character(input$adf_type %||% "drift"))
  #   adf_type     <- if (adf_type_in %in% c("none","drift","trend")) adf_type_in else "drift"
  #   kpss_null    <- norm_kpss_null(input$kpss_type, adf_type)
  #   
  #   out <- list(meta = list(N = N, sd = sdv, k = k, kpss_null = kpss_null, adf_type = adf_type))
  #   
  #   # if series unusable, bail early with reasons
  #   if (N < 8 || !is.finite(sdv) || sdv == 0) {
  #     msg <- sprintf("Invalid input for tests (N=%d, sd=%s).", N, as.character(sdv))
  #     out$adf_error  <- msg
  #     out$kpss_error <- msg
  #     out$pp_error   <- msg
  #     out$ur_error   <- msg
  #     return(out)
  #   }
  #   
  #   have_tseries <- requireNamespace("tseries", quietly = TRUE)
  #   have_urca    <- requireNamespace("urca",    quietly = TRUE)
  #   
  #   # --- tseries tests ---
  #   if (have_tseries) {
  #     out$adf <- tryCatch(
  #       tseries::adf.test(x, k = k, alternative = (input$alternd2St %||% "stationary")),
  #       error = function(e) { out$adf_error <<- e$message; NULL }
  #     )
  #     out$kpss <- tryCatch(
  #       tseries::kpss.test(x, null = kpss_null),                # "Level" or "Trend"
  #       error = function(e) { out$kpss_error <<- e$message; NULL }
  #     )
  #     out$pp <- tryCatch(
  #       tseries::pp.test(x),
  #       error = function(e) { out$pp_error <<- e$message; NULL }
  #     )
  #   } else {
  #     out$adf_error  <- "Package 'tseries' not installed."
  #     out$kpss_error <- "Package 'tseries' not installed."
  #     out$pp_error   <- "Package 'tseries' not installed."
  #   }
  #   
  #   # --- urca (ADF with tau/critical) ---
  #   if (have_urca) {
  #     out$ur <- tryCatch(
  #       urca::ur.df(x, type = adf_type, lags = k),
  #       error = function(e) { out$ur_error <<- e$message; NULL }
  #     )
  #   } else {
  #     out$ur_error <- "Package 'urca' not installed."
  #   }
  #   
  #   out
  # })
  
  
  

  # stationarity <- eventReactive(input$run_tests, {
  #   s <- ts_train_test()
  #   x <- s$ts_train
  #   k <- as.numeric(input$adf_lags)
  #   list(
  #     adf = tryCatch(tseries::adf.test(x, k = k), error = function(e) NULL),
  #     kpss = tryCatch(tseries::kpss.test(x, null = input$kpss_type), error = function(e) NULL),
  #     pp = tryCatch(tseries::pp.test(x), error = function(e) NULL),
  #     ur = tryCatch(urca::ur.df(x, type = input$adf_type, lags = k), error = function(e) NULL)
  #   )
  # })
  
  
  # stationarity <- reactive({
  #   x_raw <- tryCatch(myData_Choice(), error = function(e) NULL)
  #   out <- list()
  #   
  #   # early exits if data missing
  #   if (is.null(x_raw)) return(out)
  #   
  #   # clean & basic guards
  #   x <- as.numeric(stats::na.omit(x_raw))
  #   if (length(x) < 8 || !is.finite(stats::sd(x)) || stats::sd(x) == 0) {
  #     out$kpss_error <- sprintf("Insufficient/invalid data for KPSS (N=%d, sd=%s).",
  #                               length(x), as.character(stats::sd(x)))
  #   }
  #   
  #   # map UI -> KPSS null (tseries uses 'Level' or 'Trend')
  #   kpss_null <- if ((input$adfTypeSt2 %||% "drift") == "trend") "Trend" else "Level"
  #   
  #   # ADF (tseries)
  #   if (requireNamespace("tseries", quietly = TRUE)) {
  #     out$adf <- tryCatch(
  #       tseries::adf.test(x, alternative = (input$alternd2St %||% "stationary"),
  #                         k = as.integer(input$LagOrderADFd2St %||% 10)),
  #       error = function(e) { out$adf_error <<- e$message; NULL }
  #     )
  #     # KPSS (only attempt if basic guards pass)
  #     if (is.null(out$kpss_error)) {
  #       out$kpss <- tryCatch(
  #         tseries::kpss.test(x, null = kpss_null),
  #         error = function(e) { out$kpss_error <<- e$message; NULL }
  #       )
  #     }
  #     # PP
  #     out$pp <- tryCatch(
  #       tseries::pp.test(x),
  #       error = function(e) { out$pp_error <<- e$message; NULL }
  #     )
  #   } else {
  #     out$adf_error  <- "Package 'tseries' not installed."
  #     out$kpss_error <- "Package 'tseries' not installed."
  #     out$pp_error   <- "Package 'tseries' not installed."
  #   }
  #   
  #   # ADF (urca)
  #   if (requireNamespace("urca", quietly = TRUE)) {
  #     out$ur <- tryCatch(
  #       urca::ur.df(x, type = (input$adfTypeSt2 %||% "drift"),
  #                   lags = as.integer(input$LagOrderADFd2St %||% 10)),
  #       error = function(e) { out$ur_error <<- e$message; NULL }
  #     )
  #   } else {
  #     out$ur_error <- "Package 'urca' not installed."
  #   }
  #   
  #   out
  # })
  
  
  
  
  
  output$stationarity_results <- renderPrint({
    req(stationarity())
    st <- stationarity()
    
    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
    to_num <- function(x, d = NA_real_) {
      y <- suppressWarnings(as.numeric(x))
      if (length(y) == 0 || all(is.na(y)) || !is.finite(y[1])) d else y[1]
    }
    fmt_p   <- if (exists("fmt_p", inherits = TRUE)) get("fmt_p") else function(p) if (!is.finite(p)) "NA" else format.pval(p, digits = 4, eps = .Machine$double.eps)
    fmt_num <- if (exists("fmt_num", inherits = TRUE)) get("fmt_num") else function(x, d = 4) if (!is.finite(x)) "NA" else format(round(x, d), nsmall = d)
    
    adf_type_ui  <- input$adfTypeSt2 %||% "trend"
    k_lags       <- to_num(input$LagOrderADFd2St, 10)
    alpha        <- to_num(input$alphaSt2, 0.05)
    kpss_null_ui <- if (identical(adf_type_ui, "trend")) "tau" else "mu"
    
    tau_row   <- switch(adf_type_ui, "none"="tau1", "drift"="tau2", "trend"="tau3", "tau3")
    alpha_col <- switch(as.character(alpha), "0.01"="1pct", "0.05"="5pct", "0.1"="10pct", "0.10"="10pct", "5pct")
    
    adf_p     <- if (!is.null(st$adf))  to_num(st$adf$p.value)    else NA_real_
    adf_stat  <- if (!is.null(st$adf))  to_num(st$adf$statistic)  else NA_real_
    kpss_p    <- if (!is.null(st$kpss)) to_num(st$kpss$p.value)   else NA_real_
    kpss_stat <- if (!is.null(st$kpss)) to_num(st$kpss$statistic) else NA_real_
    pp_p      <- if (!is.null(st$pp))   to_num(st$pp$p.value)     else NA_real_
    pp_stat   <- if (!is.null(st$pp))   to_num(st$pp$statistic)   else NA_real_
    
    tau_obs  <- NA_real_
    tau_crit <- NA_real_
    if (!is.null(st$ur)) {
      tau_obs <- suppressWarnings(as.numeric(st$ur@teststat[tau_row]))
      if (!is.finite(tau_obs)) tau_obs <- suppressWarnings(as.numeric(st$ur@teststat[1]))
      cv <- tryCatch(st$ur@cval, error = function(e) NULL)
      if (!is.null(cv) && is.matrix(cv) && !is.null(rownames(cv)) && !is.null(colnames(cv)) &&
          tau_row %in% rownames(cv) && alpha_col %in% colnames(cv)) {
        tau_crit <- suppressWarnings(as.numeric(cv[tau_row, alpha_col]))
      } else if (!is.null(cv) && !is.matrix(cv) && !is.null(names(cv)) && alpha_col %in% names(cv)) {
        tau_crit <- suppressWarnings(as.numeric(cv[[alpha_col]]))
      }
    }
    
    adf_dec_ts  <- if (is.finite(adf_p))  (adf_p  < alpha) else NA
    adf_dec_ur  <- if (is.finite(tau_obs) && is.finite(tau_crit)) (tau_obs < tau_crit) else NA
    kpss_reject <- if (is.finite(kpss_p)) (kpss_p < alpha) else NA
    pp_reject   <- if (is.finite(pp_p))   (pp_p   < alpha) else NA
    
    cat("==========================================================================\n")
    cat("           ACADEMIC REPORT: STATIONARITY TESTS (TRAINING SERIES)          \n")
    cat("==========================================================================\n")
    cat(sprintf(" Significance Level (α): %.2f\n", alpha))
    cat("--------------------------------------------------------------------------\n")
    
    # 1) ADF (tseries)
    cat(" 1) AUGMENTED DICKEY–FULLER — tseries::adf.test\n")
    if (is.null(st$adf)) {
      cat("    Not available.\n")
    } else {
      cat(" USE CASE:\n")
      cat("  • Parametric unit-root test with lagged differences (need for differencing d).\n")
      cat(" HYPOTHESES:\n")
      cat("  • H0: Unit root (non-stationary)\n")
      cat("  • Ha: Stationary\n")
      cat(sprintf(" DECISION RULE (α=%.2f): Reject H0 if P-value < α.\n", alpha))
      cat(" ASSUMPTIONS & CAVEATS:\n")
      cat("  • Proper lag length k; low power near unit root.\n\n")
      cat(" STATISTICS:\n")
      cat(sprintf("  • DF statistic : %s\n", fmt_num(adf_stat, 4)))
      cat(sprintf("  • P-value      : %s\n\n", fmt_p(adf_p)))
      cat(" DECISION:\n")
      if (isTRUE(adf_dec_ts)) {
        cat("  → P-value < α ⇒ REJECT H0 ⇒ Stationarity suggested.\n")
      } else if (identical(adf_dec_ts, FALSE)) {
        cat("  → P-value ≥ α ⇒ FAIL TO REJECT H0 ⇒ Possible unit root.\n")
      } else cat("  → Inconclusive (p-value NA/Inf).\n")
      cat(" SUGGESTIONS:\n")
      cat("  • If ADF fails but KPSS rejects, prefer d=1 and re-test.\n")
      cat("--------------------------------------------------------------------------\n")
    }
    
    # 2) ADF (ur.df)
    cat(" 2) ADF — urca::ur.df (tau vs. critical)\n")
    if (is.null(st$ur)) {
      cat("    Not available.\n")
    } else {
      cat(" USE CASE:\n")
      cat("  • ADF with explicit deterministic component and chosen lag order.\n")
      cat(" HYPOTHESES:\n")
      cat("  • H0: Unit root (non-stationary)\n")
      cat("  • Ha: Stationary\n")
      cat(sprintf(" DECISION RULE (α=%.2f, %s): Reject H0 if Tau-Obs < Tau-Crit.\n", alpha, alpha_col))
      cat(" SPECIFICATION:\n")
      cat(sprintf("  • Model type : %s   • Lags (k): %s\n", adf_type_ui, as.character(k_lags)))
      cat(" ASSUMPTIONS & CAVEATS:\n")
      cat("  • Wrong type (none/drift/trend) or too-small k can bias inference.\n\n")
      cat(" STATISTICS:\n")
      cat(sprintf("  • Tau-Obs    : %s\n", fmt_num(tau_obs, 4)))
      cat(sprintf("  • Tau-Crit   : %s (%s)\n\n", fmt_num(tau_crit, 4), alpha_col))
      cat(" DECISION:\n")
      if (isTRUE(adf_dec_ur)) {
        cat("  → Tau-Obs < Tau-Crit ⇒ REJECT H0 ⇒ Stationarity supported.\n")
      } else if (identical(adf_dec_ur, FALSE)) {
        cat("  → Tau-Obs ≥ Tau-Crit ⇒ FAIL TO REJECT H0 ⇒ Possible unit root.\n")
      } else cat("  → Inconclusive (tau/critical NA).\n")
      cat(" SUGGESTIONS:\n")
      cat("  • If residuals autocorrelate (Ljung–Box), increase k or difference the series.\n")
      cat("--------------------------------------------------------------------------\n")
    }
    
    # 3) KPSS
    cat(" 3) KPSS — tseries::kpss.test\n")
    if (is.null(st$kpss)) {
      cat("    Not available.\n")
    } else {
      cat(" USE CASE:\n")
      cat("  • Tests null of stationarity (level or trend); complements ADF/PP.\n")
      cat(" HYPOTHESES:\n")
      cat(sprintf("  • H0: Stationary around a %s\n", ifelse(kpss_null_ui == "tau", "trend", "level")))
      cat("  • Ha: Non-stationary (unit root)\n")
      cat(sprintf(" DECISION RULE (α=%.2f): Reject H0 if P-value < α.\n", alpha))
      cat(" ASSUMPTIONS & CAVEATS:\n")
      cat("  • Sensitive to unremoved trend/seasonality; small-sample distortions.\n\n")
      cat(" STATISTICS:\n")
      cat(sprintf("  • KPSS η     : %s\n", fmt_num(kpss_stat, 4)))
      cat(sprintf("  • P-value    : %s\n\n", fmt_p(kpss_p)))
      cat(" DECISION:\n")
      if (isTRUE(kpss_reject)) {
        cat("  → P-value < α ⇒ REJECT H0 ⇒ Evidence of NON-STATIONARITY.\n")
      } else if (identical(kpss_reject, FALSE)) {
        cat("  → P-value ≥ α ⇒ FAIL TO REJECT H0 ⇒ Stationarity is plausible.\n")
      } else cat("  → Inconclusive (p-value NA/Inf).\n")
      cat(" SUGGESTIONS:\n")
      cat("  • If KPSS rejects and ADF/PP don’t, set d=1 (and D=1 if seasonal), then re-test.\n")
      cat("--------------------------------------------------------------------------\n")
    }
    
    # 4) PP
    cat(" 4) PHILLIPS–PERRON — tseries::pp.test\n")
    if (is.null(st$pp)) {
      cat("    Not available.\n")
    } else {
      cat(" USE CASE:\n")
      cat("  • Unit-root test with nonparametric correction (robust to serial correlation/heteroskedasticity).\n")
      cat(" HYPOTHESES:\n")
      cat("  • H0: Unit root (non-stationary)\n")
      cat("  • Ha: Stationary\n")
      cat(sprintf(" DECISION RULE (α=%.2f): Reject H0 if P-value < α.\n", alpha))
      cat(" ASSUMPTIONS & CAVEATS:\n")
      cat("  • Low power near unit root; interpret with KPSS.\n\n")
      cat(" STATISTICS:\n")
      cat(sprintf("  • PP statistic : %s\n", fmt_num(pp_stat, 4)))
      cat(sprintf("  • P-value      : %s\n\n", fmt_p(pp_p)))
      cat(" DECISION:\n")
      if (isTRUE(pp_reject)) {
        cat("  → P-value < α ⇒ REJECT H0 ⇒ Stationarity suggested.\n")
      } else if (identical(pp_reject, FALSE)) {
        cat("  → P-value ≥ α ⇒ FAIL TO REJECT H0 ⇒ Possible unit root.\n")
      } else cat("  → Inconclusive (p-value NA/Inf).\n")
      cat(" SUGGESTIONS:\n")
      cat("  • If PP & ADF fail but KPSS rejects, difference (d=1; D=1 if seasonal) and re-test.\n")
      cat("--------------------------------------------------------------------------\n")
    }
    
    # Synthesis
    cat(" SYNTHESIS & PRACTICAL RECOMMENDATIONS\n")
    if ((isTRUE(adf_dec_ts) || isTRUE(adf_dec_ur) || isTRUE(pp_reject)) && !isTRUE(kpss_reject)) {
      cat(" • Convergent STATIONARITY: use d=0; if seasonal, model with SARIMA terms rather than differencing more.\n")
    } else if ((identical(adf_dec_ts, FALSE) || identical(adf_dec_ur, FALSE) || identical(pp_reject, FALSE)) && isTRUE(kpss_reject)) {
      cat(" • Convergent UNIT ROOT: set d=1 (and D=1 if seasonal), then re-test before identification.\n")
    } else {
      cat(" • Mixed results: align deterministic spec (none/drift/trend), tune ADF lag k (check Ljung–Box),\n")
      cat("   treat seasonality first (D=1 if needed), consider log/Box–Cox for variance, and check for breaks.\n")
    }
    cat("\n==========================================================================\n\n")
  })
  
  
  
  # output$stationarity_results <- renderPrint({
  #   req(stationarity())
  #   st <- stationarity()
  #   
  #   # ---- helpers ----
  #   `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
  #   to_num  <- function(x, d = NA_real_) { y <- suppressWarnings(as.numeric(x)); ifelse(is.finite(y), y, d) }
  #   fmt_p   <- if (exists("fmt_p", inherits = TRUE)) get("fmt_p") else function(p) {
  #     if (!is.finite(p)) return("NA"); format.pval(p, digits = 4, eps = .Machine$double.eps)
  #   }
  #   fmt_num <- if (exists("fmt_num", inherits = TRUE)) get("fmt_num") else function(x, d = 4) {
  #     if (!is.finite(x)) return("NA"); format(round(x, d), nsmall = d)
  #   }
  #   
  #   # UI-driven settings (with safe defaults)
  #   adf_type_ui  <- input$adfTypeSt2 %||% "trend"  # ur.df model type
  #   k_lags       <- to_num(input$LagOrderADFd2St, 10)
  #   alpha        <- to_num(input$alphaSt2, 0.05)
  #   kpss_null_ui <- if (identical(adf_type_ui, "trend")) "tau" else "mu"
  #   
  #   # Map for ur.df critical values
  #   tau_row <- switch(adf_type_ui, "none" = "tau1", "drift" = "tau2", "trend" = "tau3", "tau3")
  #   alpha_col <- switch(as.character(alpha), "0.01"="1pct","0.05"="5pct","0.1"="10pct","0.10"="10pct","5pct")
  #   
  #   # Extract ur.df tau & critical
  #   tau_obs  <- NA_real_; tau_crit <- NA_real_
  #   if (!is.null(st$ur)) {
  #     tau_obs <- suppressWarnings(as.numeric(st$ur@teststat[tau_row]))
  #     if (!is.finite(tau_obs)) tau_obs <- suppressWarnings(as.numeric(st$ur@teststat[1]))
  #     cv <- tryCatch(st$ur@cval, error = function(e) NULL)
  #     if (!is.null(cv)) {
  #       if (!is.null(dim(cv)) && tau_row %in% rownames(cv) && alpha_col %in% colnames(cv)) {
  #         tau_crit <- suppressWarnings(as.numeric(cv[tau_row, alpha_col]))
  #       } else if (is.null(dim(cv)) && !is.null(names(cv)) && alpha_col %in% names(cv)) {
  #         tau_crit <- suppressWarnings(as.numeric(cv[[alpha_col]]))
  #       }
  #     }
  #   }
  #   
  #   # Gather main numbers
  #   adf_p    <- if (!is.null(st$adf))  to_num(st$adf$p.value)     else NA_real_
  #   adf_stat <- if (!is.null(st$adf))  to_num(st$adf$statistic)   else NA_real_
  #   kpss_p   <- if (!is.null(st$kpss)) to_num(st$kpss$p.value)    else NA_real_
  #   kpss_stat<- if (!is.null(st$kpss)) to_num(st$kpss$statistic)  else NA_real_
  #   pp_p     <- if (!is.null(st$pp))   to_num(st$pp$p.value)      else NA_real_
  #   pp_stat  <- if (!is.null(st$pp))   to_num(st$pp$statistic)    else NA_real_
  #   
  #   # Decisions
  #   adf_dec_ts  <- if (is.finite(adf_p))  (adf_p  < alpha) else NA                # reject unit root?
  #   adf_dec_ur  <- if (is.finite(tau_obs) && is.finite(tau_crit)) (tau_obs < tau_crit) else NA
  #   kpss_reject <- if (is.finite(kpss_p)) (kpss_p < alpha) else NA               # reject stationarity?
  #   pp_reject   <- if (is.finite(pp_p))   (pp_p   < alpha) else NA               # reject unit root?
  #   
  #   # ---- HEADER ----
  #   cat("==========================================================================\n")
  #   cat("                ACADEMIC REPORT: STATIONARITY TEST RESULTS                \n")
  #   cat("==========================================================================\n")
  #   cat(" Significance Level (α): ", format(alpha, nsmall = 2), "\n", sep = "")
  #   cat("--------------------------------------------------------------------------\n")
  #   
  #   # ---- 1) ADF (tseries::adf.test) ----
  #   cat(" 1. AUGMENTED DICKEY–FULLER — tseries::adf.test\n")
  #   if (is.null(st$adf)) {
  #     cat("    Not available.\n")
  #   } else {
  #     cat(" DECISION RULE:\n")
  #     cat(" • Use case: parametric unit-root test (with lagged differences) to assess if differencing is needed.\n")
  #     cat(" • H0: Unit root (non-stationary)\n")
  #     cat(" • Ha: Stationary (mean-reverting)\n")
  #     cat(" • α :", format(alpha, nsmall = 2), "\n")
  #     cat(" • Reject H0 if P-value < α.\n\n")
  #     
  #     cat(" STATISTICS:\n")
  #     cat("    - DF statistic       : ", fmt_num(adf_stat, 4), "\n", sep = "")
  #     cat("    - P-Value            : ", fmt_p(adf_p), "\n\n", sep = "")
  #     cat(" DECISION:\n")
  #     if (isTRUE(adf_dec_ts)) {
  #       cat("  P-value < α ⇒ REJECT H0: Evidence suggests STATIONARITY.\n")
  #     } else if (identical(adf_dec_ts, FALSE)) {
  #       cat("  P-value ≥ α ⇒ FAIL TO REJECT H0: Evidence suggests NON-STATIONARITY (unit root).\n")
  #     } else {
  #       cat("  ADF p-value is NA/Inf; decision is inconclusive.\n")
  #     }
  #   }
  #   cat("--------------------------------------------------------------------------\n")
  #   
  #   # ---- 2) ADF (urca::ur.df) tau vs critical ----
  #   cat(" 2. ADF — urca::ur.df (tau vs. critical)\n")
  #   if (is.null(st$ur)) {
  #     cat("    Not available.\n")
  #   } else {
  #     cat(" DECISION RULE:\n")
  #     cat(" • Use case: unit-root test with explicit deterministic component (type = none/drift/trend) and chosen k.\n")
  #     cat(" • H0: Unit root (non-stationary)\n")
  #     cat(" • Ha: Stationary\n")
  #     cat(" • α :", format(alpha, nsmall = 2), " (", alpha_col, ")\n", sep = "")
  #     cat(" • Reject H0 if Tau-Observed < Tau-Critical.\n\n")
  #     
  #     cat(" SPECIFICATION:\n")
  #     cat("    - Model type         : ", adf_type_ui, "\n", sep = "")
  #     cat("    - Lags (k)           : ", k_lags, "\n", sep = "")
  #     cat(" STATISTICS:\n")
  #     cat("    - Tau Observed       : ", fmt_num(tau_obs, 4), "\n", sep = "")
  #     cat("    - Tau Critical       : ", fmt_num(tau_crit, 4), "\n\n", sep = "")
  #     cat(" DECISION:\n")
  #     if (isTRUE(adf_dec_ur)) {
  #       cat("  Tau-Obs < Tau-Crit ⇒ REJECT H0: Stationarity supported.\n")
  #     } else if (identical(adf_dec_ur, FALSE)) {
  #       cat("  Tau-Obs ≥ Tau-Crit ⇒ FAIL TO REJECT H0: Possible unit root (non-stationary).\n")
  #     } else {
  #       cat("  Tau or critical value missing; cannot form a decision.\n")
  #     }
  #   }
  #   cat("--------------------------------------------------------------------------\n")
  #   
  #   # ---- 3) KPSS ----
  #   cat(" 3. KPSS — tseries::kpss.test\n")
  #   if (is.null(st$kpss)) {
  #     cat("    Not available.\n")
  #   } else {
  #     cat(" DECISION RULE:\n")
  #     cat(" • Use case: tests the *null* of stationarity (level/trend), complementing ADF/PP.\n")
  #     cat(" • H0: Stationary around a ", ifelse(kpss_null_ui == "tau", "trend", "level"), "\n", sep = "")
  #     cat(" • Ha: Non-stationary (unit root)\n")
  #     cat(" • α :", format(alpha, nsmall = 2), "\n")
  #     cat(" • Reject H0 if P-value < α.\n\n")
  #     
  #     cat(" STATISTICS:\n")
  #     cat("    - KPSS statistic (η) : ", fmt_num(kpss_stat, 4), "\n", sep = "")
  #     cat("    - P-Value            : ", fmt_p(kpss_p), "\n\n", sep = "")
  #     cat(" DECISION:\n")
  #     if (isTRUE(kpss_reject)) {
  #       cat("  P-value < α ⇒ REJECT H0: Evidence of NON-STATIONARITY.\n")
  #     } else if (identical(kpss_reject, FALSE)) {
  #       cat("  P-value ≥ α ⇒ FAIL TO REJECT H0: Stationarity is plausible.\n")
  #     } else {
  #       cat("  KPSS p-value missing; decision is inconclusive.\n")
  #     }
  #   }
  #   cat("--------------------------------------------------------------------------\n")
  #   
  #   # ---- 4) Phillips–Perron ----
  #   cat(" 4. PHILLIPS–PERRON — tseries::pp.test\n")
  #   if (is.null(st$pp)) {
  #     cat("    Not available.\n")
  #   } else {
  #     cat(" DECISION RULE:\n")
  #     cat(" • Use case: unit-root test with nonparametric correction for serial correlation/heteroskedasticity.\n")
  #     cat(" • H0: Unit root (non-stationary)\n")
  #     cat(" • Ha: Stationary\n")
  #     cat(" • α :", format(alpha, nsmall = 2), "\n")
  #     cat(" • Reject H0 if P-value < α.\n\n")
  #     
  #     cat(" STATISTICS:\n")
  #     cat("    - PP statistic       : ", fmt_num(pp_stat, 4), "\n", sep = "")
  #     cat("    - P-Value            : ", fmt_p(pp_p), "\n\n", sep = "")
  #     cat(" DECISION:\n")
  #     if (isTRUE(pp_reject)) {
  #       cat("  P-value < α ⇒ REJECT H0: Stationarity suggested.\n")
  #     } else if (identical(pp_reject, FALSE)) {
  #       cat("  P-value ≥ α ⇒ FAIL TO REJECT H0: Possible unit root (non-stationary).\n")
  #     } else {
  #       cat("  PP p-value missing; decision is inconclusive.\n")
  #     }
  #   }
  #   cat("--------------------------------------------------------------------------\n")
  #   
  #   # ---- Synthesis (optional short) ----
  #   cat(" SYNTHESIS:\n")
  #   cat(" • Combine test evidence with diagnostic plots to select differencing orders d and D.\n")
  #   cat("==========================================================================\n\n")
  # })
  
  
  

  # output$stationarity_results <- renderPrint({
  #   req(stationarity())
  #   st <- stationarity()
  #   cat("Stationarity tests (training series)\n\n")
  #   cat("ADF (tseries::adf.test)\n"); if (is.null(st$adf)) cat("  Not available.\n") else print(st$adf)
  #   cat("\nKPSS (tseries::kpss.test)\n"); if (is.null(st$kpss)) cat("  Not available.\n") else print(st$kpss)
  #   cat("\nPhillips–Perron (tseries::pp.test)\n"); if (is.null(st$pp)) cat("  Not available.\n") else print(st$pp)
  #   cat("\nADF (urca::ur.df) summary\n"); if (is.null(st$ur)) cat("  Not available.\n") else print(summary(st$ur))
  # })
  
  
  
  
  output$stationarity_interpretation <- renderPrint({
    req(stationarity())
    st <- stationarity()
    
    # ---- helpers ----
    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
    to_num  <- function(x, d = NA_real_) { y <- suppressWarnings(as.numeric(x)); ifelse(is.finite(y), y, d) }
    # prefer your global helpers if they exist; otherwise fall back
    fmt_p   <- if (exists("fmt_p", inherits = TRUE)) get("fmt_p") else function(p) {
      if (!is.finite(p)) return("NA"); format.pval(p, digits = 4, eps = .Machine$double.eps)
    }
    fmt_num <- if (exists("fmt_num", inherits = TRUE)) get("fmt_num") else function(x, d = 4) {
      if (!is.finite(x)) return("NA"); format(round(x, d), nsmall = d)
    }
    
    # UI-driven settings (with safe defaults)
    adf_type_ui  <- input$adfTypeSt2 %||% "trend"  # ur.df model type
    k_lags       <- to_num(input$LagOrderADFd2St, 10)
    kpss_null_ui <- if (adf_type_ui == "trend") "tau" else "mu"  # typical pairing
    alpha        <- to_num(input$alphaSt2, 0.05)
    
    # Map for ur.df critical values
    tau_row <- switch(adf_type_ui, "none" = "tau1", "drift" = "tau2", "trend" = "tau3", "tau3")
    alpha_col <- switch(as.character(alpha), "0.01"="1pct","0.05"="5pct","0.1"="10pct","0.10"="10pct","5pct")
    
    # Try to extract ur.df tau stat & crit
    tau_obs  <- NA_real_; tau_crit <- NA_real_
    if (!is.null(st$ur)) {
      tau_obs <- suppressWarnings(as.numeric(st$ur@teststat[tau_row]))
      if (!is.finite(tau_obs)) tau_obs <- suppressWarnings(as.numeric(st$ur@teststat[1]))
      cv <- tryCatch(st$ur@cval, error = function(e) NULL)
      if (!is.null(cv)) {
        if (!is.null(dim(cv)) && tau_row %in% rownames(cv) && alpha_col %in% colnames(cv)) {
          tau_crit <- suppressWarnings(as.numeric(cv[tau_row, alpha_col]))
        } else if (is.null(dim(cv)) && !is.null(names(cv)) && alpha_col %in% names(cv)) {
          tau_crit <- suppressWarnings(as.numeric(cv[[alpha_col]]))
        }
      }
    }
    
    # Gather main numbers
    adf_p    <- if (!is.null(st$adf))  to_num(st$adf$p.value)     else NA_real_
    adf_stat <- if (!is.null(st$adf))  to_num(st$adf$statistic)   else NA_real_
    kpss_p   <- if (!is.null(st$kpss)) to_num(st$kpss$p.value)    else NA_real_
    kpss_stat<- if (!is.null(st$kpss)) to_num(st$kpss$statistic)  else NA_real_
    pp_p     <- if (!is.null(st$pp))   to_num(st$pp$p.value)      else NA_real_
    pp_stat  <- if (!is.null(st$pp))   to_num(st$pp$statistic)    else NA_real_
    
    # Decisions
    adf_dec_ts  <- if (is.finite(adf_p))  (adf_p  < alpha) else NA
    adf_dec_ur  <- if (is.finite(tau_obs) && is.finite(tau_crit)) (tau_obs < tau_crit) else NA
    kpss_reject <- if (is.finite(kpss_p)) (kpss_p < alpha) else NA      # reject stationarity
    pp_reject   <- if (is.finite(pp_p))   (pp_p   < alpha) else NA      # reject unit root
    
    # ---- HEADER ----
    cat("==========================================================================\n")
    cat("                 ACADEMIC REPORT: STATIONARITY ANALYSIS                   \n")
    cat("==========================================================================\n")
    cat(" Significance Level (α): ", format(alpha, nsmall = 2), "\n", sep = "")
    cat("--------------------------------------------------------------------------\n")
    
    # ---- 1) ADF (tseries::adf.test) ----
    cat(" 1. AUGMENTED DICKEY–FULLER — tseries::adf.test\n")
    if (is.null(st$adf)) {
      cat("    Not available.\n")
    } else {
      cat(" DECISION RULE:\n")
      cat(" • Use case: detects a unit root via a parametric regression with lagged differences;\n")
      cat("   useful to confirm if differencing (d) is needed to achieve stationarity.\n")
      cat(" • H0: the series has a unit root (non-stationary)\n")
      cat(" • Ha: the series is stationary (mean-reverting)\n")
      cat(" • α :", format(alpha, nsmall = 2), "\n")
      cat(" • Reject H0 if P-value < α.\n\n")
      
      cat(" STATISTICS:\n")
      cat("    - Test Statistic (DF) :", fmt_num(adf_stat, 4), "\n")
      cat("    - P-Value             :", fmt_p(adf_p), "\n\n")
      
      cat(" DECISION:\n")
      if (isTRUE(adf_dec_ts)) {
        cat("  P-value < α ⇒ REJECT H0: Evidence suggests STATIONARITY.\n")
      } else if (identical(adf_dec_ts, FALSE)) {
        cat("  P-value ≥ α ⇒ FAIL TO REJECT H0: Evidence suggests NON-STATIONARITY (unit root).\n")
      } else {
        cat("  ADF p-value is NA/Inf; decision is inconclusive.\n")
      }
    }
    cat("--------------------------------------------------------------------------\n")
    
    # ---- 2) ADF (urca::ur.df) tau vs critical ----
    cat(" 2. ADF — urca::ur.df (tau vs. critical)\n")
    if (is.null(st$ur)) {
      cat("    Not available.\n")
    } else {
      cat(" DECISION RULE:\n")
      cat(" • Use case: complementary ADF formulation with explicit deterministic component\n")
      cat("   (type = none/drift/trend) and chosen lag order k; compares τ to critical values.\n")
      cat(" • H0: the series has a unit root (non-stationary)\n")
      cat(" • Ha: the series is stationary\n")
      cat(" • α :", format(alpha, nsmall = 2), " (", alpha_col, ")\n", sep = "")
      cat(" • Reject H0 if Tau-Observed < Tau-Critical.\n\n")
      
      cat(" SPECIFICATION:\n")
      cat("    - Model type         : ", adf_type_ui, "\n", sep = "")
      cat("    - Lags (k)           : ", k_lags, "\n", sep = "")
      cat(" STATISTICS:\n")
      cat("    - Tau Observed       : ", fmt_num(tau_obs, 4), "\n", sep = "")
      cat("    - Tau Critical       : ", fmt_num(tau_crit, 4), "\n\n", sep = "")
      
      cat(" DECISION:\n")
      if (isTRUE(adf_dec_ur)) {
        cat("  Tau-Obs < Tau-Crit ⇒ REJECT H0: Stationarity supported.\n")
      } else if (identical(adf_dec_ur, FALSE)) {
        cat("  Tau-Obs ≥ Tau-Crit ⇒ FAIL TO REJECT H0: Possible unit root (non-stationary).\n")
      } else {
        cat("  Tau or critical value missing; cannot form a decision.\n")
      }
    }
    cat("--------------------------------------------------------------------------\n")
    
    # ---- 3) KPSS ----
    cat(" 3. KPSS — tseries::kpss.test\n")
    if (is.null(st$kpss)) {
      cat("    Not available.\n")
    } else {
      cat(" DECISION RULE:\n")
      cat(" • Use case: tests the *null* of stationarity (around a level or trend);\n")
      cat("   complements ADF/PP by flagging residual non-stationarity when ADF/PP are inconclusive.\n")
      cat(" • H0: the series is stationary around a ", ifelse(kpss_null_ui == "tau", "trend", "level"), "\n", sep = "")
      cat(" • Ha: the series is non-stationary (has a unit root)\n")
      cat(" • α :", format(alpha, nsmall = 2), "\n")
      cat(" • Reject H0 if P-value < α.\n\n")
      
      cat(" STATISTICS:\n")
      cat("    - Test Statistic (η) : ", fmt_num(kpss_stat, 4), "\n", sep = "")
      cat("    - P-Value            : ", fmt_p(kpss_p), "\n\n", sep = "")
      
      cat(" DECISION:\n")
      if (isTRUE(kpss_reject)) {
        cat("  P-value < α ⇒ REJECT H0: Evidence of NON-STATIONARITY.\n")
      } else if (identical(kpss_reject, FALSE)) {
        cat("  P-value ≥ α ⇒ FAIL TO REJECT H0: Stationarity is plausible.\n")
      } else {
        cat("  KPSS p-value missing; decision is inconclusive.\n")
      }
    }
    cat("--------------------------------------------------------------------------\n")
    
    # ---- 4) Phillips–Perron ----
    cat(" 4. PHILLIPS–PERRON — tseries::pp.test\n")
    if (is.null(st$pp)) {
      cat("    Not available.\n")
    } else {
      cat(" DECISION RULE:\n")
      cat(" • Use case: unit-root test robust to serial correlation/heteroskedasticity\n")
      cat("   (nonparametric corrections); complements ADF.\n")
      cat(" • H0: the series has a unit root (non-stationary)\n")
      cat(" • Ha: the series is stationary\n")
      cat(" • α :", format(alpha, nsmall = 2), "\n")
      cat(" • Reject H0 if P-value < α.\n\n")
      
      cat(" STATISTICS:\n")
      cat("    - Test Statistic     : ", fmt_num(pp_stat, 4), "\n", sep = "")
      cat("    - P-Value            : ", fmt_p(pp_p), "\n\n", sep = "")
      
      cat(" DECISION:\n")
      if (isTRUE(pp_reject)) {
        cat("  P-value < α ⇒ REJECT H0: Stationarity suggested.\n")
      } else if (identical(pp_reject, FALSE)) {
        cat("  P-value ≥ α ⇒ FAIL TO REJECT H0: Possible unit root (non-stationary).\n")
      } else {
        cat("  PP p-value missing; decision is inconclusive.\n")
      }
    }
    cat("--------------------------------------------------------------------------\n")
    
    # ---- Final synthesis ----
    cat(" RECOMMENDATION:\n")
    cat(" • Combine test evidence with plots (Time, ACF/PACF) to choose differencing orders d and D.\n")
    if ( (isTRUE(adf_dec_ts) || isTRUE(adf_dec_ur) || isTRUE(pp_reject)) && !isTRUE(kpss_reject) ) {
      cat(" • Convergent evidence for STATIONARITY (ADF/PP agree; KPSS does not reject): consider d = 0 (check seasonality).\n")
    } else if ( (identical(adf_dec_ts, FALSE) || identical(adf_dec_ur, FALSE) || identical(pp_reject, FALSE)) && isTRUE(kpss_reject) ) {
      cat(" • Convergent evidence for NON-STATIONARITY: apply differencing (d and/or D), then re-test.\n")
    } else {
      cat(" • Mixed signals: be conservative — prefer differencing and validate via residual diagnostics/forecast performance.\n")
    }
    cat("\n==========================================================================\n\n")
  })
  
  
  
  # output$stationarity_interpretation <- renderPrint({
  #   req(stationarity())
  #   st <- stationarity()
  #   
  #   # ---- small helpers (local) ----
  #   `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
  #   to_num  <- function(x, d = NA_real_) { y <- suppressWarnings(as.numeric(x)); ifelse(is.finite(y), y, d) }
  #   fmt_p   <- get("fmt_p",   envir = parent.env(environment()))  # use your global helpers
  #   fmt_num <- get("fmt_num", envir = parent.env(environment()))
  #   
  #   # Inputs (with safe defaults)
  #   adf_type_ui  <- input$adf_type  %||% "trend"  # for ur.df
  #   kpss_null_ui <- input$kpss_type %||% "mu"
  #   k_lags       <- to_num(input$adf_lags, 10)
  #   alpha        <- to_num(input$alphaSt2, 0.05)  # will fall back to 0.05 if not present in this tab
  #   
  #   # Map for ur.df rows/cols (critical values)
  #   tau_row <- switch(adf_type_ui, "none" = "tau1", "drift" = "tau2", "trend" = "tau3", "tau3")
  #   alpha_col <- switch(as.character(alpha),
  #                       "0.01"="1pct","0.05"="5pct","0.1"="10pct","0.10"="10pct","5pct")
  #   
  #   # Try to extract ur.df tau stat and critical value
  #   tau_obs  <- NA_real_; tau_crit <- NA_real_
  #   if (!is.null(st$ur)) {
  #     # observed
  #     tau_obs <- suppressWarnings(as.numeric(st$ur@teststat[tau_row]))
  #     if (!is.finite(tau_obs)) {
  #       tau_obs <- suppressWarnings(as.numeric(st$ur@teststat[1]))
  #     }
  #     # critical
  #     cv <- tryCatch(st$ur@cval, error = function(e) NULL)
  #     if (!is.null(cv) && !is.null(dim(cv)) && tau_row %in% rownames(cv) && alpha_col %in% colnames(cv)) {
  #       tau_crit <- suppressWarnings(as.numeric(cv[tau_row, alpha_col]))
  #     } else if (!is.null(cv) && is.null(dim(cv)) && !is.null(names(cv)) && alpha_col %in% names(cv)) {
  #       tau_crit <- suppressWarnings(as.numeric(cv[[alpha_col]]))
  #     }
  #   }
  #   
  #   # Extract main test numbers
  #   adf_p   <- if (!is.null(st$adf)) to_num(st$adf$p.value) else NA_real_
  #   adf_stat<- if (!is.null(st$adf)) to_num(st$adf$statistic) else NA_real_
  #   
  #   kpss_p   <- if (!is.null(st$kpss)) to_num(st$kpss$p.value) else NA_real_
  #   kpss_stat<- if (!is.null(st$kpss)) to_num(st$kpss$statistic) else NA_real_
  #   
  #   pp_p     <- if (!is.null(st$pp)) to_num(st$pp$p.value) else NA_real_
  #   pp_stat  <- if (!is.null(st$pp)) to_num(st$pp$statistic) else NA_real_
  #   
  #   # Decisions
  #   # ADF (ur.df rule): reject H0 (unit root) if tau_obs < tau_crit
  #   adf_dec_ur <- if (is.finite(tau_obs) && is.finite(tau_crit)) (tau_obs < tau_crit) else NA
  #   adf_dec_ts <- if (is.finite(adf_p)) (adf_p < alpha) else NA  # tseries reference
  #   
  #   # KPSS: reject H0 (stationary) if p < alpha
  #   kpss_reject <- if (is.finite(kpss_p)) (kpss_p < alpha) else NA
  #   
  #   # PP: reject H0 (unit root) if p < alpha
  #   pp_reject <- if (is.finite(pp_p)) (pp_p < alpha) else NA
  #   
  #   # ---- HEADER ----
  #   cat("==========================================================================\n")
  #   cat("                 ACADEMIC REPORT: STATIONARITY ANALYSIS                   \n")
  #   cat("==========================================================================\n")
  #   cat(" DECISION RULES:\n")
  #   cat(" • ADF (tseries::adf.test): H0 = Unit Root (non-stationary).\n")
  #   cat("     Reject H0 if P-value < α.\n")
  #   cat(" • ADF (urca::ur.df): H0 = Unit Root; type =", adf_type_ui, ", lags k =", k_lags, "\n")
  #   cat("     Reject H0 if Tau-Observed < Tau-Critical (", alpha_col, ").\n", sep = "")
  #   cat(" • KPSS (tseries::kpss.test): H0 = Stationary (null =", kpss_null_ui, ").\n")
  #   cat("     Reject H0 if P-value < α.\n")
  #   cat(" • PP (tseries::pp.test): H0 = Unit Root (non-stationary).\n")
  #   cat("     Reject H0 if P-value < α.\n")
  #   cat(" • Significance Level (α): ", format(alpha, nsmall = 2), "\n", sep = "")
  #   cat("--------------------------------------------------------------------------\n")
  #   
  #   # ---- 1) ADF (tseries) ----
  #   cat(" 1. AUGMENTED DICKEY–FULLER (tseries)\n")
  #   if (is.null(st$adf)) {
  #     cat("    Not available.\n")
  #   } else {
  #     cat("    - Test Statistic (DF) :", fmt_num(adf_stat, 4), "\n")
  #     cat("    - P-Value             :", fmt_p(adf_p), "\n\n")
  #     cat(" RESULT:\n")
  #     cat("  The ADF test from tseries evaluates the unit-root hypothesis in a regression\n")
  #     cat("  with lagged differences (k =", k_lags, ").\n\n")
  #     cat(" DECISION:\n")
  #     if (isTRUE(adf_dec_ts)) {
  #       cat("  At α =", format(alpha, nsmall = 2), " the p-value indicates REJECTION of H0.\n")
  #       cat("  Evidence suggests the series is STATIONARY.\n")
  #     } else if (identical(adf_dec_ts, FALSE)) {
  #       cat("  At α =", format(alpha, nsmall = 2), " the p-value is not small enough.\n")
  #       cat("  FAIL TO REJECT H0: Evidence suggests NON-STATIONARITY (unit root).\n")
  #     } else {
  #       cat("  ADF p-value is NA/Inf; decision is inconclusive.\n")
  #     }
  #   }
  #   cat("--------------------------------------------------------------------------\n")
  #   
  #   # ---- 2) ADF (urca::ur.df) with critical value ----
  #   cat(" 2. ADF (urca::ur.df) — Tau vs Critical\n")
  #   if (is.null(st$ur)) {
  #     cat("    Not available.\n")
  #   } else {
  #     cat("    - Tau Observed        :", fmt_num(tau_obs, 4), "\n")
  #     cat("    - Tau Critical (", alpha_col, "): ", fmt_num(tau_crit, 4), "\n", sep = "")
  #     cat("    - Model type          : ", adf_type_ui, " | Lags k: ", k_lags, "\n", sep = "")
  #     cat("\n RESULT:\n")
  #     cat("  The ur.df framework compares the observed tau statistic to its critical value\n")
  #     cat("  under the specified deterministic component (none/drift/trend).\n\n")
  #     cat(" DECISION:\n")
  #     if (isTRUE(adf_dec_ur)) {
  #       cat("  Tau-Obs < Tau-Crit ⇒ REJECT H0: Stationarity supported.\n")
  #     } else if (identical(adf_dec_ur, FALSE)) {
  #       cat("  Tau-Obs ≥ Tau-Crit ⇒ FAIL TO REJECT H0: Possible unit root (non-stationary).\n")
  #     } else {
  #       cat("  Tau or critical value missing; cannot form a decision.\n")
  #     }
  #   }
  #   cat("--------------------------------------------------------------------------\n")
  #   
  #   # ---- 3) KPSS ----
  #   cat(" 3. KPSS (tseries::kpss.test)\n")
  #   if (is.null(st$kpss)) {
  #     cat("    Not available.\n")
  #   } else {
  #     cat("    - Test Statistic (η)  :", fmt_num(kpss_stat, 4), "\n")
  #     cat("    - P-Value             :", fmt_p(kpss_p), "\n")
  #     cat("    - Null hypothesis     : Stationary around ", ifelse(kpss_null_ui == "tau", "a trend", "a level"), "\n", sep = "")
  #     cat("\n RESULT:\n")
  #     cat("  Small p-values indicate deviations from stationarity with respect to the chosen null.\n\n")
  #     cat(" DECISION:\n")
  #     if (isTRUE(kpss_reject)) {
  #       cat("  P-value < α ⇒ REJECT H0: Evidence of NON-STATIONARITY.\n")
  #     } else if (identical(kpss_reject, FALSE)) {
  #       cat("  P-value ≥ α ⇒ FAIL TO REJECT H0: Stationarity is plausible.\n")
  #     } else {
  #       cat("  KPSS p-value missing; decision is inconclusive.\n")
  #     }
  #   }
  #   cat("--------------------------------------------------------------------------\n")
  #   
  #   # ---- 4) Phillips–Perron ----
  #   cat(" 4. PHILLIPS–PERRON (tseries::pp.test)\n")
  #   if (is.null(st$pp)) {
  #     cat("    Not available.\n")
  #   } else {
  #     cat("    - Test Statistic      :", fmt_num(pp_stat, 4), "\n")
  #     cat("    - P-Value             :", fmt_p(pp_p), "\n\n")
  #     cat(" RESULT:\n")
  #     cat("  The PP test addresses serial correlation/heteroskedasticity nonparametrically.\n\n")
  #     cat(" DECISION:\n")
  #     if (isTRUE(pp_reject)) {
  #       cat("  P-value < α ⇒ REJECT H0: Stationarity suggested.\n")
  #     } else if (identical(pp_reject, FALSE)) {
  #       cat("  P-value ≥ α ⇒ FAIL TO REJECT H0: Possible unit root (non-stationary).\n")
  #     } else {
  #       cat("  PP p-value missing; decision is inconclusive.\n")
  #     }
  #   }
  #   cat("--------------------------------------------------------------------------\n")
  #   
  #   # ---- Final Advice ----
  #   cat(" ADVICE:\n")
  #   if (isTRUE(adf_dec_ts) || isTRUE(adf_dec_ur)) {
  #     if (!isTRUE(kpss_reject) && isTRUE(pp_reject)) {
  #       cat(" • Convergent evidence (ADF/PP) with KPSS non-rejection ⇒ stationarity is well supported.\n")
  #     } else if (isTRUE(kpss_reject)) {
  #       cat(" • ADF/PP vs KPSS conflict: consider differencing (d or D), check model type (none/drift/trend),\n")
  #       cat("   and inspect for structural breaks/seasonality.\n")
  #     } else {
  #       cat(" • ADF indicates stationarity; verify visually and with residual diagnostics.\n")
  #     }
  #   } else if (identical(adf_dec_ts, FALSE) || identical(adf_dec_ur, FALSE) || identical(pp_reject, FALSE)) {
  #     cat(" • Evidence of a unit root: apply differencing (d and/or D), then re-test.\n")
  #   } else {
  #     cat(" • Mixed/inconclusive: rely on plots (time, ACF/PACF) and adopt a conservative differencing choice.\n")
  #   }
  #   cat("\n==========================================================================\n\n")
  # })
  
  
  
  

  # output$stationarity_interpretation <- renderPrint({
  #   req(stationarity())
  #   st <- stationarity()
  #   cat("Interpretation (heuristic):\n\n")
  #   cat("- ADF: H0 = unit root (non-stationary). Small p suggests stationarity.\n")
  #   cat("- KPSS: H0 = stationary. Small p suggests non-stationary.\n")
  #   cat("- PP: H0 = unit root (non-stationary). Small p suggests stationarity.\n\n")
  #   if (!is.null(st$adf)) cat("ADF:", fmt_p(st$adf$p.value), "\n")
  #   if (!is.null(st$kpss)) cat("KPSS:", fmt_p(st$kpss$p.value), "\n")
  #   if (!is.null(st$pp)) cat("PP:", fmt_p(st$pp$p.value), "\n\n")
  #   cat("Recommendation: use combined evidence (tests + plots) to decide d and D.\n")
  # })

  diff_preview <- eventReactive(input$preview_diff, {
    s <- ts_train_test()
    x <- ts(c(as.numeric(s$ts_train), as.numeric(s$ts_test)), start = 1, frequency = frequency(s$ts_train))
    d <- as.numeric(input$d_preview)
    D <- as.numeric(input$D_preview)
    x2 <- x
    if (d > 0) x2 <- diff(x2, differences = d)
    if (D > 0) x2 <- diff(x2, lag = frequency(x), differences = D)
    x2
  })

  output$diff_plot <- renderPlot({ req(diff_preview()); plot(diff_preview(), main = "Differenced series preview", ylab = "Value", xlab = "Time") })

  output$apa_stationarity_paragraph <- renderPrint({
    req(stationarity(), prepared())
    st <- stationarity()
    s_per <- prepared()$freq
    adf_txt <- if (!is.null(st$adf)) paste0("An Augmented Dickey–Fuller test indicated ", ifelse(st$adf$p.value < 0.05, "stationarity", "non-stationarity"), ", ", fmt_p(st$adf$p.value), ". ") else ""
    kpss_txt <- if (!is.null(st$kpss)) paste0("A KPSS test suggested ", ifelse(st$kpss$p.value < 0.05, "non-stationarity", "stationarity"), ", ", fmt_p(st$kpss$p.value), ". ") else ""
    pp_txt <- if (!is.null(st$pp)) paste0("A Phillips–Perron test suggested ", ifelse(st$pp$p.value < 0.05, "stationarity", "non-stationarity"), ", ", fmt_p(st$pp$p.value), ". ") else ""
    cat(
      "APA-ready paragraph (edit differencing decisions as needed):\n\n",
      adf_txt, kpss_txt, pp_txt,
      "Based on combined evidence and inspection of differenced series plots, appropriate non-seasonal (d) and seasonal (D) differencing ",
      "were selected prior to fitting SARIMA models (season length s = ", s_per, ").\n",
      sep = ""
    )
  })

  # ---- Step 5 Auto-ARIMA ----

  auto_fit <- eventReactive(input$fit_auto, {
    s <- ts_train_test()
    forecast::auto.arima(
      s$ts_train,
      seasonal = isTRUE(input$auto_seasonal),
      stepwise = isTRUE(input$auto_stepwise),
      approximation = isTRUE(input$auto_approx),
      allowmean = isTRUE(input$auto_allow_mean),
      allowdrift = isTRUE(input$auto_allow_mean),
      max.order = as.numeric(input$auto_max_order)
    )
  })

  auto_fc <- reactive({
    req(auto_fit(), ts_train_test(), prepared())
    s <- ts_train_test()
    p <- prepared()
    user_h <- if (!is.na(input$auto_h) && input$auto_h > 0) as.numeric(input$auto_h) else NA_real_
    if (s$test_n > 0) {
      h <- s$test_n
    } else {
      h <- if (!is.na(user_h)) user_h else max(1, frequency(s$ts_train))
    }
    fc <- forecast::forecast(auto_fit(), h = h)
    list(fc = fc, h = h, by = p$by)
  })

  output$auto_horizon_note <- renderPrint({
    req(ts_train_test(), auto_fc())
    s <- ts_train_test()
    if (s$test_n > 0) {
      cat("Validation mode: forecast horizon is forced to the test length (h =", s$test_n, ") to overlay the test period.\n")
    } else {
      cat("Future mode: no test set. Forecast horizon uses your 'Future horizon h' setting (or defaults to one seasonal cycle).\n")
    }
  })

  output$auto_model_spec <- renderPrint({
    req(auto_fit(), prepared())
    fit <- auto_fit()
    p <- prepared()
    cat("Auto-ARIMA selected model:\n")
    cat("- Specification:", as.character(fit), "\n")
    cat("- Seasonal period (s):", p$freq, "\n")
    cat("- Information criteria: AICc =", fmt_num(fit$aicc, 2), ", BIC =", fmt_num(fit$bic, 2), "\n")
  })

  output$auto_coef_table <- renderTable({ req(auto_fit()); coef_table(auto_fit()) }, rownames = FALSE)
  output$auto_resid_ts <- renderPlot({ req(auto_fit()); plot(residuals(auto_fit()), main = "Residuals (Auto-ARIMA)", ylab = "Residual", xlab = "Time") })
  output$auto_resid_acf <- renderPlot({ req(auto_fit()); plot(acf(residuals(auto_fit()), plot = FALSE), main = "Residual ACF (Auto-ARIMA)") })
  output$auto_resid_hist <- renderPlot({ req(auto_fit()); hist(residuals(auto_fit()), breaks = 30, main = "Residual histogram", xlab = "Residual") })
  output$auto_resid_qq <- renderPlot({ req(auto_fit()); qqnorm(residuals(auto_fit())); qqline(residuals(auto_fit())) })
  
  
  
  
  
  
  
  
  
  #===============================================================
  #===============================================================
  #===============================================================
  
  
  
  # ============================================================
  # --- FIXED: Auto-ARIMA equation renderer (MathJax-safe) ---
  # ============================================================
  
  auto_equations <- reactive({
    req(auto_fit())
    fit <- auto_fit()
    
    # --- Orders from forecast::Arima (arma = c(p, q, P, Q, m, d, D)) ---
    arma <- fit$arma
    p <- arma[1]; q <- arma[2]; P <- arma[3]; Q <- arma[4]
    s <- ifelse(isTRUE(arma[5] > 0), arma[5], 1L)
    d <- arma[6]; D <- arma[7]
    
    coefs <- stats::coef(fit)
    nm    <- names(coefs)
    
    ip <- if (p > 0) seq_len(p) else integer(0)
    iq <- if (q > 0) seq_len(q) else integer(0)
    iP <- if (P > 0) seq_len(P) else integer(0)
    iQ <- if (Q > 0) seq_len(Q) else integer(0)
    
    # Intercept/mean (show only when meaningful)
    intercept_name <- intersect(c("intercept", "mean"), nm)
    intercept_val  <- if (length(intercept_name) > 0) unname(coefs[intercept_name[1]]) else NA_real_
    show_intercept <- is.finite(intercept_val) && abs(intercept_val) > 1e-8 && d == 0 && D == 0
    intercept_num  <- if (show_intercept) sprintf("%.3f", intercept_val) else ""
    
    # Drift (present only when differencing is used and auto.arima included it)
    drift_val <- if ("drift" %in% nm) unname(coefs["drift"]) else NA_real_
    show_drift <- is.finite(drift_val) && abs(drift_val) > 1e-8 && (d > 0 || D > 0)
    
    # MathJax wrappers
    tex_display <- function(x) paste0("\\[", x, "\\]")
    
    # Cosmetic cleanup for numeric line
    simplify_tex <- function(x) {
      x <- gsub("\\(1\\)", "", x)
      x <- gsub("\\s+", " ", x)
      x <- gsub("\\+\\s*\\+", "+", x)
      x <- gsub("\\+\\s*-", "-", x)
      x <- gsub("-\\s*\\+", "-", x)
      x <- gsub("-\\s*-", "+", x)
      x <- gsub("\\s*\\+\\s*0\\.000\\b", "", x)
      trimws(x)
    }
    
    # ---- Parameter-polynomial strings (symbolic) ----
    poly_param_ar  <- function() if (p == 0) "1" else paste0("1", paste0(" - \\phi_{", ip, "}L^{", ip, "}", collapse = ""))
    poly_param_sar <- function() if (P == 0) "1" else paste0("1", paste0(" - \\Phi_{", iP, "}L^{", s * iP, "}", collapse = ""))
    poly_param_ma  <- function() if (q == 0) "1" else paste0("1", paste0(" + \\theta_{", iq, "}L^{", iq, "}", collapse = ""))
    poly_param_sma <- function() if (Q == 0) "1" else paste0("1", paste0(" + \\Theta_{", iQ, "}L^{", s * iQ, "}", collapse = ""))
    
    # ---- Numeric polynomials (use estimated coefficients) ----
    poly_num_ar <- function() {
      if (p == 0) return("1")
      v <- suppressWarnings(unname(coefs[paste0("ar", ip)]))
      keep <- is.finite(v)
      if (!any(keep)) return("1")
      paste0("1", paste(sprintf(" %+.3fL^{%d}", -v[keep], ip[keep]), collapse = ""))
    }
    poly_num_sar <- function() {
      if (P == 0) return("1")
      v <- suppressWarnings(unname(coefs[paste0("sar", iP)]))
      keep <- is.finite(v)
      if (!any(keep)) return("1")
      lags <- s * iP
      paste0("1", paste(sprintf(" %+.3fL^{%d}", -v[keep], lags[keep]), collapse = ""))
    }
    poly_num_ma <- function() {
      if (q == 0) return("1")
      v <- suppressWarnings(unname(coefs[paste0("ma", iq)]))
      keep <- is.finite(v)
      if (!any(keep)) return("1")
      paste0("1", paste(sprintf(" %+.3fL^{%d}", v[keep], iq[keep]), collapse = ""))
    }
    poly_num_sma <- function() {
      if (Q == 0) return("1")
      v <- suppressWarnings(unname(coefs[paste0("sma", iQ)]))
      keep <- is.finite(v)
      if (!any(keep)) return("1")
      lags <- s * iQ
      paste0("1", paste(sprintf(" %+.3fL^{%d}", v[keep], lags[keep]), collapse = ""))
    }
    
    # Differencing operators (omit if exponent is zero)
    diff_part  <- if (d > 0) paste0("(1-L)^{", d, "}") else ""
    sdiff_part <- if (D > 0) paste0("(1-L^{", s, "})^{", D, "}") else ""
    
    # ---------- Line 1: General operator form ----------
    # \phi_p(L)\Phi_P(L^S)(1-L)^d(1-L^S)^D Y_t = c + \theta_q(L)\Theta_Q(L^S)\varepsilon_t (+ drift when appropriate)
    line1 <- paste0(
      "\\phi_p(L)\\,\\Phi_P(L^{S})\\,(1-L)^{d}(1-L^{S})^{D}Y_t",
      " = c + \\theta_q(L)\\,\\Theta_Q(L^{S})\\varepsilon_t",
      if (show_drift) " + \\delta t" else ""
    )
    
    # ---------- Line 2: Expanded operator (summation) ----------
    line2 <- paste0(
      "\\left(1-\\sum_{i=1}^{p}\\phi_i L^{i}\\right)",
      "\\left(1-\\sum_{j=1}^{P}\\Phi_j L^{jS}\\right)",
      "(1-L)^{d}(1-L^{S})^{D}Y_t",
      " = c + ",
      "\\left(1+\\sum_{i=1}^{q}\\theta_i L^{i}\\right)",
      "\\left(1+\\sum_{j=1}^{Q}\\Theta_j L^{jS}\\right)",
      "\\varepsilon_t",
      if (show_drift) " + \\delta t" else ""
    )
    
    # ---------- Line 3: Parameter-expanded polynomials ----------
    line3 <- paste0(
      "(", poly_param_ar(), ")",
      "(", poly_param_sar(), ")",
      diff_part, sdiff_part,
      "Y_t = c + ",
      "(", poly_param_ma(), ")",
      "(", poly_param_sma(), ")\\varepsilon_t",
      if (show_drift) " + \\delta t" else ""
    )
    
    # ---------- Line 4: Numeric-expanded polynomials ----------
    rhs_intercept <- if (show_intercept) paste0(intercept_num, " + ") else ""
    line4 <- paste0(
      "(", poly_num_ar(), ")",
      "(", poly_num_sar(), ")",
      diff_part, sdiff_part,
      "Y_t = ",
      rhs_intercept,
      "(", poly_num_ma(), ")",
      "(", poly_num_sma(), ")\\varepsilon_t",
      if (show_drift) paste0(" + ", sprintf("%.3f", drift_val), "t") else ""
    )
    line4 <- simplify_tex(line4)
    
    # ---------- Coefficient legend ----------
    coef_lines <- c()
    if (show_intercept) coef_lines <- c(coef_lines, paste0("c (intercept/mean) = ", sprintf("%.4f", intercept_val)))
    if (show_drift)     coef_lines <- c(coef_lines, paste0("drift (\\(\\delta\\)) = ", sprintf("%.4f", drift_val)))
    
    if (p > 0) {
      for (i in ip) {
        nm_i <- paste0("ar", i)
        if (nm_i %in% nm && is.finite(coefs[nm_i])) {
          coef_lines <- c(coef_lines, paste0("\\(\\phi_", i, "\\) = ", sprintf("%.4f", coefs[nm_i])))
        }
      }
    }
    if (q > 0) {
      for (i in iq) {
        nm_i <- paste0("ma", i)
        if (nm_i %in% nm && is.finite(coefs[nm_i])) {
          coef_lines <- c(coef_lines, paste0("\\(\\theta_", i, "\\) = ", sprintf("%.4f", coefs[nm_i])))
        }
      }
    }
    if (P > 0) {
      for (i in iP) {
        nm_i <- paste0("sar", i)
        if (nm_i %in% nm && is.finite(coefs[nm_i])) {
          coef_lines <- c(coef_lines, paste0("\\(\\Phi_", i, "\\) = ", sprintf("%.4f", coefs[nm_i])))
        }
      }
    }
    if (Q > 0) {
      for (i in iQ) {
        nm_i <- paste0("sma", i)
        if (nm_i %in% nm && is.finite(coefs[nm_i])) {
          coef_lines <- c(coef_lines, paste0("\\(\\Theta_", i, "\\) = ", sprintf("%.4f", coefs[nm_i])))
        }
      }
    }
    if (length(coef_lines) == 0) coef_lines <- "No coefficients available."
    
    list(
      p = p, d = d, q = q, P = P, D = D, Q = Q, s = s,
      eq_general      = tex_display(line1),
      eq_expanded     = tex_display(line2),
      eq_line3        = tex_display(line3),
      eq_line4        = tex_display(line4),
      coef_lines      = coef_lines
    )
  })
  
  output$auto_model_equation <- renderUI({
    req(auto_equations())
    eq <- auto_equations()
    
    tagList(
      tags$div(
        style = "text-align:left;",
        tags$h4("Auto-ARIMA model"),
        tags$p(sprintf("ARIMA(%d,%d,%d)%s",
                       eq$p, eq$d, eq$q,
                       if (eq$s > 1) sprintf(" × (%d,%d,%d)[%d]", eq$P, eq$D, eq$Q, eq$s) else "")),
        
        tags$h5("Estimated coefficients"),
        tags$ul(lapply(eq$coef_lines, function(x) tags$li(HTML(x)))),
        
        tags$hr(),
        
        tags$h4("General SARIMA formulation"),
        HTML(eq$eq_general),
        
        tags$hr(),
        
        tags$h4("Expanded operator form"),
        HTML(eq$eq_expanded),
        
        tags$hr(),
        
        tags$h4("Numerical model"),
        HTML(eq$eq_line3),
        tags$hr(),
        
        # HTML("\\[\\text{------------}\\]"),
        HTML(eq$eq_line4),
        tags$hr(),
        
      ),
      
      # Force MathJax typesetting for dynamically injected content
      tags$script(HTML("
      if (window.MathJax) {
        if (window.MathJax.Hub) { MathJax.Hub.Queue(['Typeset', MathJax.Hub]); }
        else if (window.MathJax.typesetPromise) { MathJax.typesetPromise(); }
      }
    "))
    )
  })
  
  
  
  
  
  
  # # ---- Helpers for equation formatting (put once in server.R) ----
  # latex_poly <- function(k, greek, seasonal = FALSE, s = 1, type = c("AR","MA")) {
  #   # Returns "1 - phi_1 B - ... - phi_k B^k" (AR) or "1 + theta_1 B + ... + theta_k B^k" (MA)
  #   type <- match.arg(type)
  #   if (k <= 0) return("1")
  #   opsign <- if (type == "AR") "-" else "+"
  #   terms <- vapply(seq_len(k), function(i) {
  #     exp_txt <- if (seasonal) paste0(i, "\\,", s) else i
  #     paste0(" ", opsign, " \\\\", greek, "_", i, " B^{", exp_txt, "}")
  #   }, character(1))
  #   paste0("1", paste0(terms, collapse = ""))
  # }
  # 
  # latex_diff <- function(d, D, s) {
  #   # (1 - B)^d (1 - B^s)^D   with smart omission if exponent is 0
  #   parts <- c()
  #   if (d > 0) parts <- c(parts, paste0("(1 - B)^{", d, "}"))
  #   if (D > 0) parts <- c(parts, paste0("(1 - B^{", s, "})^{", D, "}"))
  #   if (length(parts) == 0) "1" else paste(parts, collapse = " ")
  # }
  # 
  # latex_coeff_list <- function(coefs, greek, label_prefix = "", wrap_dollars = TRUE) {
  #   # coefs is a numeric vector; greek like "phi"/"theta"/"Phi"/"Theta"
  #   if (length(coefs) == 0) return(NULL)
  #   items <- vapply(seq_along(coefs), function(i) {
  #     val <- sprintf("%.4f", coefs[i])
  #     sym <- paste0("\\\\", greek, "_", i)
  #     if (wrap_dollars) {
  #       paste0("<li>\\(", label_prefix, sym, " = ", val, "\\)</li>")
  #     } else {
  #       paste0("<li>", label_prefix, sym, " = ", val, "</li>")
  #     }
  #   }, character(1))
  #   paste0("<ul style='margin:4px 0 10px 20px;'>", paste(items, collapse = ""), "</ul>")
  # }
  # 
  # # ---- NEW: Auto-ARIMA model equation output ----
  # output$auto_model_equation <- renderUI({
  #   req(auto_fit())
  #   fit <- auto_fit()
  #   
  #   # Extract ARIMA structure from forecast::Arima
  #   # arma = c(p, q, P, Q, m, d, D)
  #   arma <- fit$arma
  #   p <- arma[1]; q <- arma[2]; P <- arma[3]; Q <- arma[4]
  #   s <- arma[5]; d <- arma[6]; D <- arma[7]
  #   
  #   # Coefficient lookup by names (forecast uses "ar", "ma", "sar", "sma", "intercept", "drift")
  #   cf   <- stats::coef(fit)
  #   nm   <- names(cf)
  #   phi  <- unname(cf[grep("^ar\\d+$", nm)])          # non-seasonal AR
  #   theta<- unname(cf[grep("^ma\\d+$", nm)])          # non-seasonal MA
  #   PhiS <- unname(cf[grep("^sar\\d+$", nm)])         # seasonal AR
  #   TheS <- unname(cf[grep("^sma\\d+$", nm)])         # seasonal MA
  #   mu   <- if (any(nm == "intercept")) unname(cf["intercept"]) else NA_real_
  #   drift<- if (any(nm == "drift"))     unname(cf["drift"])     else NA_real_
  #   sig2 <- if (!is.null(fit$sigma2)) fit$sigma2 else NA_real_
  #   
  #   # Build symbolic polynomials
  #   AR_ns   <- latex_poly(p,  "phi", seasonal = FALSE, s = s, type = "AR")
  #   MA_ns   <- latex_poly(q,  "theta", seasonal = FALSE, s = s, type = "MA")
  #   AR_seas <- latex_poly(P,  "Phi", seasonal = TRUE,  s = s, type = "AR")
  #   MA_seas <- latex_poly(Q,  "Theta", seasonal = TRUE, s = s, type = "MA")
  #   DIFF    <- latex_diff(d, D, s)
  #   
  #   # Constant / drift term string
  #   rhs_const <- NULL
  #   if (is.finite(mu) && d == 0 && D == 0) {
  #     rhs_const <- "\\; + \\; \\mu"
  #   } else if (is.finite(drift) && (d > 0 || D > 0)) {
  #     rhs_const <- "\\; + \\; \\delta \\, t"
  #   } else {
  #     rhs_const <- ""
  #   }
  #   
  #   # Main symbolic equation
  #   eq <- paste0(
  #     "\\[",
  #     "\\underbrace{(", AR_ns, ")}_{\\text{AR(", p, ")}}\\,",
  #     "\\underbrace{(", AR_seas, ")}_{\\text{Seasonal AR(", P, ")}}\\,",
  #     "\\underbrace{", DIFF, "}_{\\text{Differencing}}\\, y_t",
  #     " \\;=\\; ",
  #     "\\underbrace{(", MA_ns, ")}_{\\text{MA(", q, ")}}\\,",
  #     "\\underbrace{(", MA_seas, ")}_{\\text{Seasonal MA(", Q, ")}}\\, \\varepsilon_t",
  #     rhs_const,
  #     "\\]"
  #   )
  #   
  #   # Coefficient legend (only those that exist)
  #   legend_html <- c(
  #     if (length(phi))  latex_coeff_list(phi,  "phi"),
  #     if (length(theta))latex_coeff_list(theta,"theta"),
  #     if (length(PhiS)) latex_coeff_list(PhiS,"Phi"),
  #     if (length(TheS)) latex_coeff_list(TheS,"Theta"),
  #     if (is.finite(mu) && d == 0 && D == 0)
  #       sprintf("<div>\\(\\mu\\) (intercept) = <b>%.4f</b></div>", mu),
  #     if (is.finite(drift) && (d > 0 || D > 0))
  #       sprintf("<div>\\(\\delta\\) (drift) = <b>%.4f</b></div>", drift),
  #     if (is.finite(sig2))
  #       sprintf("<div>\\(\\sigma^2_{\\varepsilon}\\) (innovation variance) = <b>%.4f</b></div>", sig2)
  #   )
  #   legend_html <- paste(Filter(Negate(is.null), legend_html), collapse = "\n")
  #   
  #   # Header with compact summary of orders
  #   head_html <- sprintf(
  #     "<div style='margin-bottom:8px;'>
  #      <b>Selected model:</b> ARIMA(%d,%d,%d)%s
  #      <br/><b>Seasonal period (s)</b> = %d
  #    </div>",
  #     p, d, q,
  #     if (s > 1) sprintf(" × (%d,%d,%d)[%d]", P, D, Q, s) else "",
  #     s
  #   )
  #   
  #   HTML(paste0(head_html, "<div>", eq, "</div>",
  #               "<div style='margin-top:6px;'><b>Estimated coefficients</b>:</div>",
  #               legend_html))
  # })
  # 
  # 
  
  #===============================================================
  #===============================================================
  #===============================================================
  
  
  # output$auto_diag_tests <- renderText({ req(auto_fit()); diag_tests_text(residuals(auto_fit()), lag = as.numeric(input$diag_lag), fitdf = length(coef(auto_fit()))) })

  
  output$auto_diag_tests <- renderText({
    req(auto_fit())
    res   <- residuals(auto_fit())
    L     <- as.integer(input$diag_lag %||% 12L)
    fitdf <- length(coef(auto_fit()))
    alpha <- to_num_safe(input$alphaSt2 %||% 0.05, 0.05)
    residual_diagnostics_report(res, L = L, fitdf = fitdf, alpha = alpha)
  })
  
  
  # output$auto_diag_tests <- renderText({
  #   req(auto_fit())
  #   residual_diagnostics_report(
  #     res   = residuals(auto_fit()),
  #     L     = as.numeric(input$diag_lag),
  #     fitdf = length(coef(auto_fit())),
  #     alpha = to_num_safe(input$alphaSt2 %||% 0.05, 0.05)
  #   )
  # })
  
  
  
  residual_diagnostics_report <- function(res, L = 12L, fitdf = 0L, alpha = 0.05) {
    # ---- Safe helpers (use local fallbacks if not already defined globally) ----
    if (!exists("fmt_p", inherits = TRUE)) {
      fmt_p <- function(p) {
        if (!is.finite(p)) return("NA")
        if (p < .001) "<0.001" else sprintf("%.6f", p)
      }
    }
    if (!exists("fmt_num", inherits = TRUE)) {
      fmt_num <- function(z, d = 6) ifelse(is.finite(z), sprintf(paste0("%.", d, "f"), z), "NA")
    }
    
    # ---- Sanitize inputs ----
    res   <- as.numeric(res)
    res   <- res[is.finite(res)]
    N     <- length(res)
    L     <- max(1L, as.integer(L))
    fitdf <- max(0L, as.integer(fitdf))
    alpha <- suppressWarnings(as.numeric(alpha)); if (!is.finite(alpha)) alpha <- 0.05
    if (N < 8) {
      return(
        "==========================================================================\n" %+%
          "                     RESIDUAL DIAGNOSTIC BATTERY                          \n" %+%
          "==========================================================================\n" %+%
          "Too few residuals (N < 8) to run the requested tests."
      )
    }
    
    # Convenience
    `%+%` <- function(a, b) paste0(a, b)
    df_lb <- max(L - fitdf, 1L)             # Ljung–Box / Box–Pierce df after parameter adjustment
    cv_lb <- stats::qchisq(1 - alpha, df = df_lb)
    
    # ---- 1) Ljung–Box (portmanteau for autocorrelation) ----
    lb <- tryCatch(stats::Box.test(res, lag = L, type = "Ljung-Box", fitdf = fitdf), error = function(e) NULL)
    lb_stat <- if (!is.null(lb)) as.numeric(lb$statistic) else NA_real_
    lb_p    <- if (!is.null(lb)) as.numeric(lb$p.value)   else NA_real_
    
    # ---- 2) Box–Pierce (classic portmanteau) ----
    bp <- tryCatch(stats::Box.test(res, lag = L, type = "Box-Pierce", fitdf = fitdf), error = function(e) NULL)
    bp_stat <- if (!is.null(bp)) as.numeric(bp$statistic) else NA_real_
    bp_p    <- if (!is.null(bp)) as.numeric(bp$p.value)   else NA_real_
    
    # ---- 3) Jarque–Bera (normality of residuals; large-sample χ^2_2) ----
    jb <- if (requireNamespace("tseries", quietly = TRUE))
      tryCatch(tseries::jarque.bera.test(res), error = function(e) NULL) else NULL
    jb_stat <- if (!is.null(jb)) as.numeric(jb$statistic) else NA_real_
    jb_p    <- if (!is.null(jb)) as.numeric(jb$p.value)   else NA_real_
    cv_jb   <- stats::qchisq(1 - alpha, df = 2)  # asymptotic
    
    # ---- 4) Shapiro–Wilk (normality; exact test for N in [3, 5000]) ----
    sw <- if (N >= 3 && N <= 5000)
      tryCatch(stats::shapiro.test(res), error = function(e) NULL) else NULL
    sw_W <- if (!is.null(sw)) as.numeric(sw$statistic) else NA_real_
    sw_p <- if (!is.null(sw)) as.numeric(sw$p.value)   else NA_real_
    
    # ---- 5) Engle’s ARCH LM (conditional heteroskedasticity) ----
    arch_lags <- min(L, max(1L, floor(N/10)))
    arch <- if (requireNamespace("FinTS", quietly = TRUE))
      tryCatch(FinTS::ArchTest(res, lags = arch_lags), error = function(e) NULL) else NULL
    arch_stat <- if (!is.null(arch)) as.numeric(arch$statistic) else NA_real_
    arch_p    <- if (!is.null(arch)) as.numeric(arch$p.value)   else NA_real_
    cv_arch   <- stats::qchisq(1 - alpha, df = arch_lags)
    
    # ---- 6) Runs test (randomness / independence of signs) ----
    run <- if (requireNamespace("tseries", quietly = TRUE))
      tryCatch(tseries::runs.test(res), error = function(e) NULL) else NULL
    run_Z <- if (!is.null(run)) as.numeric(run$statistic) else NA_real_
    run_p <- if (!is.null(run)) as.numeric(run$p.value)   else NA_real_
    zcrit <- stats::qnorm(1 - alpha/2)  # two-sided
    
    # ---- 7) Anderson–Darling for normality (tail-sensitive; optional) ----
    ad <- if (requireNamespace("nortest", quietly = TRUE))
      tryCatch(nortest::ad.test(res), error = function(e) NULL) else NULL
    ad_A2 <- if (!is.null(ad)) as.numeric(ad$statistic) else NA_real_
    ad_p  <- if (!is.null(ad)) as.numeric(ad$p.value)   else NA_real_
    
    # ---- 8) KPSS on residuals (should be stationary if model is adequate; optional) ----
    kpss_res <- if (requireNamespace("tseries", quietly = TRUE))
      tryCatch(tseries::kpss.test(res, null = "Level"), error = function(e) NULL) else NULL
    kpss_stat <- if (!is.null(kpss_res)) as.numeric(kpss_res$statistic) else NA_real_
    kpss_p    <- if (!is.null(kpss_res)) as.numeric(kpss_res$p.value)   else NA_real_
    
    # ---- Compose the academic-friendly text ----
    out <- c(
      "==========================================================================",
      "                     RESIDUAL DIAGNOSTIC BATTERY                          ",
      "==========================================================================",
      sprintf(" SAMPLE SIZE (residuals used): %d   |   α: %s   |   Lag (L): %d   |   fitdf: %d",
              N, fmt_num(alpha, 4), L, fitdf),
      "--------------------------------------------------------------------------",
      
      # Ljung–Box
      "TEST 1: Ljung–Box Portmanteau (autocorrelation)",
      " Purpose     : To determine whether any linear autocorrelation remains in the residuals up to lag L after fitting the SARIMA model.",
      "               Passing this test supports the idea that the model has successfully captured the serial dependence in the training data.",
      " Description : The Ljung–Box statistic aggregates squared sample autocorrelations of the residuals across lags 1..L,",
      "               with a small-sample correction. Under the null hypothesis that residuals are white noise, Q follows",
      "               approximately a chi-square distribution with degrees of freedom equal to max(L − fitdf, 1).",
      sprintf(" Statistic   : Q(LB) = %s  |  df = %d  |  p-value = %s", fmt_num(lb_stat, 4), df_lb, fmt_p(lb_p)),
      sprintf(" Critical    : χ^2_(%d, 1-α) = %s", df_lb, fmt_num(cv_lb, 4)),
      sprintf(" Decision    : Reject H0 (white noise) if Q(LB) > χ^2_(%d, 1-α) (equivalently, p < α).", df_lb),
      sprintf(" Result      : %s",
              if (!is.finite(lb_p)) {
                "Inconclusive (statistic/p-value unavailable)."
              } else if (lb_p < alpha) {
                "Reject H0 → residual autocorrelation detected."
              } else {
                "Fail to reject H0 → residuals consistent with white noise."
              }),
      sprintf(" Conclusion  : %s",
              if (!is.finite(lb_p)) {
                "Because the test could not be evaluated, do not draw conclusions about left-over autocorrelation. Recheck model estimation and sample size."
              } else if (lb_p < alpha) {
                "There is statistical evidence of remaining serial correlation up to the chosen lag. This suggests the SARIMA orders may be underspecified (e.g., too few AR/MA or seasonal terms) or that differencing/seasonality was not fully addressed. Consider revisiting (p,d,q)(P,D,Q)m, increasing L modestly, or examining ACF/PACF of residuals to guide refinement."
              } else {
                "No material linear autocorrelation is detected at the examined lags. This supports the adequacy of the fitted orders for the training data and increases confidence that forecast errors are not biased by serial dependence left in the residuals."
              }),
      "--------------------------------------------------------------------------",
      
      # Box–Pierce
      "TEST 2: Box–Pierce Portmanteau (autocorrelation, classic)",
      " Purpose     : Same goal as Ljung–Box—screen for any left-over autocorrelation up to lag L—using the original Box–Pierce statistic.",
      " Description : The Box–Pierce statistic is the uncorrected sum of squared residual autocorrelations over lags 1..L.",
      "               Under the white-noise null, it is approximately χ^2 with df = max(L − fitdf, 1).",
      sprintf(" Statistic   : Q(BP) = %s  |  df = %d  |  p-value = %s", fmt_num(bp_stat, 4), df_lb, fmt_p(bp_p)),
      sprintf(" Critical    : χ^2_(%d, 1-α) = %s", df_lb, fmt_num(cv_lb, 4)),
      sprintf(" Decision    : Reject H0 if Q(BP) > χ^2_(%d, 1-α) (equivalently, p < α).", df_lb),
      sprintf(" Result      : %s",
              if (!is.finite(bp_p)) {
                "Inconclusive (statistic/p-value unavailable)."
              } else if (bp_p < alpha) {
                "Reject H0 → residual autocorrelation detected (classic statistic)."
              } else {
                "Fail to reject H0 → no strong evidence of residual autocorrelation."
              }),
      sprintf(" Conclusion  : %s",
              if (!is.finite(bp_p)) {
                "Since the test could not be computed, rely on Ljung–Box and graphical diagnostics. If both are inconclusive, try a different L or expand the sample."
              } else if (bp_p < alpha) {
                "Autocorrelation persists by the Box–Pierce criterion. In practice, prioritize Ljung–Box; if both agree, refine the model orders or differencing. If they disagree, prefer Ljung–Box due to its small-sample correction."
              } else {
                "Results align with white-noise residuals by the Box–Pierce criterion. Together with Ljung–Box, this strengthens the case that the fitted SARIMA captured the main serial structure."
              }),
      "--------------------------------------------------------------------------",
      
      # Jarque–Bera
      "TEST 3: Jarque–Bera (normality, large-sample χ^2)",
      " Purpose     : Assess whether residuals are approximately normal—important when using normal-based prediction intervals and likelihood-based inference.",
      " Description : The statistic combines squared residual skewness and excess kurtosis. Under normality (large N), JB ~ χ^2 with 2 degrees of freedom.",
      sprintf(" Statistic   : JB = %s  |  df = 2  |  p-value = %s", fmt_num(jb_stat, 4), fmt_p(jb_p)),
      sprintf(" Critical    : χ^2_(2, 1-α) = %s", fmt_num(cv_jb, 4)),
      " Decision    : Reject H0 (normality) if JB > χ^2_(2, 1-α) (equivalently, p < α).",
      sprintf(" Result      : %s",
              if (!is.finite(jb_p)) {
                "Inconclusive (statistic/p-value unavailable)."
              } else if (jb_p < alpha) {
                "Reject H0 → residuals deviate from normality."
              } else {
                "Fail to reject H0 → residual normality is plausible."
              }),
      sprintf(" Conclusion  : %s",
              if (!is.finite(jb_p)) {
                "Because the test could not be evaluated, rely on the QQ-plot and Shapiro–Wilk (if available) to judge normality. Normality mainly affects interval calibration rather than point forecasts."
              } else if (jb_p < alpha) {
                "Non-normal residuals indicate skewness and/or heavy tails. Point forecasts remain unbiased if the mean is correctly specified, but nominal prediction intervals may undercover or overcover. Consider variance-stabilizing transforms (e.g., log/Box–Cox), robust inference, or heavier-tailed error models if tails matter."
              } else {
                "Residuals are consistent with normality at the chosen α. This supports the use of standard Gaussian prediction intervals and likelihood-based comparisons for competing SARIMA specifications."
              }),
      "--------------------------------------------------------------------------",
      
      # Shapiro–Wilk
      "TEST 4: Shapiro–Wilk (normality, small/medium samples)",
      " Purpose     : Provide a powerful small-sample test of normality when N ≤ 5000.",
      " Description : The W statistic compares ordered residuals to expected normal order statistics. Lower W indicates departure from normality.",
      sprintf(" Statistic   : W = %s  |  p-value = %s  |  Range: N ∈ [3, 5000]",
              fmt_num(sw_W, 4), if (N > 5000) "n/a (N>5000)" else fmt_p(sw_p)),
      " Critical    : R relies on p-value; explicit critical values are not printed.",
      " Decision    : Reject H0 (normality) if p < α.",
      sprintf(" Result      : %s",
              if (N > 5000) {
                "Omitted (Shapiro–Wilk is defined for N ≤ 5000)."
              } else if (!is.finite(sw_p)) {
                "Inconclusive (statistic/p-value unavailable)."
              } else if (sw_p < alpha) {
                "Reject H0 → residuals deviate from normality."
              } else {
                "Fail to reject H0 → residual normality is plausible."
              }),
      sprintf(" Conclusion  : %s",
              if (N > 5000) {
                "For large samples, rely on Jarque–Bera and visual diagnostics (histogram, QQ-plot)."
              } else if (!is.finite(sw_p)) {
                "Unable to conclude about normality from Shapiro–Wilk. Cross-check with QQ-plot and Jarque–Bera before changing the model."
              } else if (sw_p < alpha) {
                "Evidence against normality in smaller samples warrants caution with Gaussian intervals. Inspect QQ-plots for systematic S-shapes (tails) or asymmetry (skew). Consider transformations or alternative error distributions if interval accuracy is a goal."
              } else {
                "Normality appears adequate by Shapiro–Wilk. Combined with QQ-plot and JB, this supports the standard Gaussian assumptions used in SARIMA teaching examples."
              }),
      "--------------------------------------------------------------------------",
      
      # ARCH LM
      "TEST 5: Engle’s ARCH LM (conditional heteroskedasticity)",
      sprintf(" Purpose     : Check whether residual variance changes over time (ARCH effects) up to %d lags—important when volatility clustering is present.", arch_lags),
      " Description : Regress squared residuals on their own lags and test whether lag coefficients are jointly zero. Under H0 (no ARCH), LM ~ χ^2 with df equal to the lag count.",
      sprintf(" Statistic   : LM = %s  |  df = %d  |  p-value = %s", fmt_num(arch_stat, 4), arch_lags, fmt_p(arch_p)),
      sprintf(" Critical    : χ^2_(%d, 1-α) = %s", arch_lags, fmt_num(cv_arch, 4)),
      " Decision    : Reject H0 (no ARCH) if LM > χ^2_(lags, 1-α) (equivalently, p < α).",
      sprintf(" Result      : %s",
              if (is.null(arch)) {
                "Skipped (package 'FinTS' not installed)."
              } else if (!is.finite(arch_p)) {
                "Inconclusive (statistic/p-value unavailable)."
              } else if (arch_p < alpha) {
                "Reject H0 → ARCH present."
              } else {
                "Fail to reject H0 → no strong evidence of ARCH."
              }),
      sprintf(" Conclusion  : %s",
              if (is.null(arch)) {
                "Install the 'FinTS' package to assess time-varying volatility rigorously. In the meantime, inspect a rolling variance plot; volatility clustering can bias interval calibration."
              } else if (!is.finite(arch_p)) {
                "ARCH could not be assessed; consider re-running with a different lag count or larger sample. Visual checks of squared residuals may still reveal volatility clusters."
              } else if (arch_p < alpha) {
                "Time-varying variance is indicated. If volatility matters (finance, energy), consider augmenting the mean model with a GARCH-type variance model or using heteroskedasticity-robust intervals. If mean forecasts are the sole target, point forecasts can still be useful but intervals should be treated with caution."
              } else {
                "No strong evidence of ARCH effects. Constant-variance assumptions used in basic SARIMA are reasonable for these residuals, improving confidence in prediction interval calibration."
              }),
      "--------------------------------------------------------------------------",
      
      # Runs test
      "TEST 6: Runs Test (randomness of signs)",
      " Purpose     : Evaluate whether positive and negative residuals occur in a random order.",
      "               Systematic alternation or clustering of signs may indicate remaining structure (e.g., nonlinearity or omitted seasonal effects).",
      " Description : Counts runs of consecutive positive/negative residuals. Under independence, the standardized count Z is approximately N(0,1).",
      sprintf(" Statistic   : Z = %s  |  p-value = %s  |  Two-sided z-crit = ±%s",
              fmt_num(run_Z, 4), fmt_p(run_p), fmt_num(zcrit, 3)),
      " Critical    : Reject H0 if |Z| > z_crit (equivalently, p < α).",
      sprintf(" Result      : %s",
              if (is.null(run)) {
                "Skipped (package 'tseries' not installed)."
              } else if (!is.finite(run_p)) {
                "Inconclusive (statistic/p-value unavailable)."
              } else if (run_p < alpha) {
                "Reject H0 → residual signs are not random."
              } else {
                "Fail to reject H0 → residual signs appear random."
              }),
      sprintf(" Conclusion  : %s",
              if (is.null(run)) {
                "Install 'tseries' to include the runs test in your diagnostic battery. Without it, rely on the residual time plot to spot sign clustering."
              } else if (!is.finite(run_p)) {
                "Because the statistic could not be computed, defer to visual checks. Persistent runs may point to nonlinearity or missing seasonal terms."
              } else if (run_p < alpha) {
                "Non-random sign patterns suggest unresolved structure (e.g., threshold dynamics or seasonal mismatches). Inspect seasonal plots and consider adding nonlinearity or revisiting seasonal differencing."
              } else {
                "Signs occur in a pattern consistent with randomness, which complements portmanteau tests by checking a simple, intuitive notion of independence."
              }),
      "--------------------------------------------------------------------------",
      
      # Anderson–Darling (optional)
      "TEST 7: Anderson–Darling (normality, tail-sensitive) [optional]",
      " Purpose     : Detect departures from normality with particular emphasis on the tails—useful when extremes matter for risk or service-level planning.",
      " Description : The A² statistic compares the empirical CDF of residuals to the normal CDF, weighting discrepancies more heavily in the tails.",
      sprintf(" Statistic   : A² = %s  |  p-value = %s", fmt_num(ad_A2, 4), if (!is.null(ad)) fmt_p(ad_p) else "n/a (package 'nortest' not installed)"),
      " Critical    : Decision by p-value (package provides the calibration).",
      sprintf(" Result      : %s",
              if (is.null(ad)) {
                "Skipped (package 'nortest' not installed)."
              } else if (!is.finite(ad_p)) {
                "Inconclusive (statistic/p-value unavailable)."
              } else if (ad_p < alpha) {
                "Reject H0 → tail behavior not normal."
              } else {
                "Fail to reject H0 → tails are consistent with normality."
              }),
      sprintf(" Conclusion  : %s",
              if (is.null(ad)) {
                "If tail accuracy matters, consider installing 'nortest' or inspecting extreme quantiles in the QQ-plot. Heavy tails call for robust or heavy-tailed error models."
              } else if (!is.finite(ad_p)) {
                "Unable to assess tail behavior quantitatively. Rely on QQ-plots focusing on the extremes before modifying the model."
              } else if (ad_p < alpha) {
                "Evidence of tail non-normality suggests Gaussian prediction intervals may under- or over-cover in the extremes. Consider transformations, t-errors, or bootstrapped intervals if extreme events are decision-critical."
              } else {
                "Tail behavior looks compatible with normality, supporting standard Gaussian intervals for planning purposes."
              }),
      "--------------------------------------------------------------------------",
      
      # KPSS on residuals (optional)
      "TEST 8: KPSS on Residuals (level-stationarity) [optional]",
      " Purpose     : Verify that residuals are stationary around a fixed level (no unit root), as expected from a well-specified SARIMA model.",
      " Description : The KPSS statistic accumulates partial sums of residuals; large values indicate non-stationarity. The test is run with the 'Level' null.",
      sprintf(" Statistic   : KPSS = %s  |  p-value = %s", fmt_num(kpss_stat, 5),
              if (!is.null(kpss_res)) fmt_p(kpss_p) else "n/a (package 'tseries' not installed)"),
      " Critical    : Decision by p-value or package-reported thresholds.",
      sprintf(" Result      : %s",
              if (is.null(kpss_res)) {
                "Skipped (package 'tseries' not installed)."
              } else if (!is.finite(kpss_p)) {
                "Inconclusive (statistic/p-value unavailable)."
              } else if (kpss_p < alpha) {
                "Reject H0 → residuals may contain non-stationary structure."
              } else {
                "Fail to reject H0 → residuals behave as stationary noise."
              }),
      sprintf(" Conclusion  : %s",
              if (is.null(kpss_res)) {
                "Install 'tseries' to include KPSS on residuals. In its absence, rely on the residual ACF and time plot to confirm stationarity."
              } else if (!is.finite(kpss_p)) {
                "Stationarity of residuals could not be certified statistically. Check for slow drifts or level shifts; if present, revisit differencing or break handling."
              } else if (kpss_p < alpha) {
                "Some non-stationary behavior remains in the residuals, which can undermine the white-noise assumption and degrade forecast uncertainty estimates. Consider additional differencing, seasonal adjustments, or modeling identified breaks."
              } else {
                "Residuals appear stationary around a fixed mean, consistent with a correctly differenced and seasonally specified SARIMA model."
              }),
      
      "==========================================================================",
      "INTERPRETATION GUIDE",
      " • Good SARIMA residuals typically: (i) pass Ljung–Box/Box–Pierce (no linear autocorrelation),",
      "   (ii) show no strong ARCH unless volatility modeling is intended, and (iii) are roughly normal",
      "   if you rely on Gaussian prediction intervals. When diagnostics fail, prefer changing the model",
      "   (orders, differencing, seasonality) rather than over-fitting residuals with ad-hoc fixes."
    )
    
    paste(out, collapse = "\n")
  }
  
  
  # residual_diagnostics_report <- function(res, L = 12L, fitdf = 0L, alpha = 0.05) {
  #   # sanitize inputs
  #   res <- as.numeric(res)
  #   res <- res[is.finite(res)]
  #   N   <- length(res)
  #   L   <- max(1L, as.integer(L))
  #   fitdf <- max(0L, as.integer(fitdf))
  #   alpha <- to_num_safe(alpha, 0.05)
  #   if (N < 8) {
  #     return("Residual Diagnostic Battery\n----------------------------------------\nToo few residuals (N < 8) to run the requested tests.")
  #   }
  #   
  #   # convenience
  #   df_lb <- max(L - fitdf, 1L)  # Ljung–Box / Box–Pierce df after parameter adjustment
  #   cv_lb <- stats::qchisq(1 - alpha, df = df_lb)
  #   
  #   # 1) Ljung–Box (portmanteau for autocorrelation)
  #   lb <- tryCatch(stats::Box.test(res, lag = L, type = "Ljung-Box", fitdf = fitdf), error = function(e) NULL)
  #   lb_stat <- if (!is.null(lb)) as.numeric(lb$statistic) else NA_real_
  #   lb_p    <- if (!is.null(lb)) as.numeric(lb$p.value)   else NA_real_
  #   
  #   # 2) Box–Pierce (classic portmanteau)
  #   bp <- tryCatch(stats::Box.test(res, lag = L, type = "Box-Pierce", fitdf = fitdf), error = function(e) NULL)
  #   bp_stat <- if (!is.null(bp)) as.numeric(bp$statistic) else NA_real_
  #   bp_p    <- if (!is.null(bp)) as.numeric(bp$p.value)   else NA_real_
  #   
  #   # 3) Jarque–Bera (normality of residuals; large-sample χ^2_2)
  #   jb <- if (requireNamespace("tseries", quietly = TRUE))
  #     tryCatch(tseries::jarque.bera.test(res), error = function(e) NULL) else NULL
  #   jb_stat <- if (!is.null(jb)) as.numeric(jb$statistic) else NA_real_
  #   jb_p    <- if (!is.null(jb)) as.numeric(jb$p.value)   else NA_real_
  #   cv_jb   <- stats::qchisq(1 - alpha, df = 2)  # JB ~ χ^2(2) asymptotically
  #   
  #   # 4) Shapiro–Wilk (normality; exact W, but no simple CV to print)
  #   sw <- if (N >= 3 && N <= 5000)
  #     tryCatch(stats::shapiro.test(res), error = function(e) NULL) else NULL
  #   sw_W <- if (!is.null(sw)) as.numeric(sw$statistic) else NA_real_
  #   sw_p <- if (!is.null(sw)) as.numeric(sw$p.value)   else NA_real_
  #   
  #   # 5) Engle’s ARCH LM (conditional heteroskedasticity)
  #   # Uses FinTS::ArchTest if available; df = lags
  #   arch_lags <- min(L, max(1L, floor(N/10)))
  #   arch <- if (requireNamespace("FinTS", quietly = TRUE))
  #     tryCatch(FinTS::ArchTest(res, lags = arch_lags), error = function(e) NULL) else NULL
  #   arch_stat <- if (!is.null(arch)) as.numeric(arch$statistic) else NA_real_
  #   arch_p    <- if (!is.null(arch)) as.numeric(arch$p.value)   else NA_real_
  #   cv_arch   <- stats::qchisq(1 - alpha, df = arch_lags)
  #   
  #   # 6) Runs test (randomness / independence)
  #   run <- if (requireNamespace("tseries", quietly = TRUE))
  #     tryCatch(tseries::runs.test(res), error = function(e) NULL) else NULL
  #   run_Z <- if (!is.null(run)) as.numeric(run$statistic) else NA_real_
  #   run_p <- if (!is.null(run)) as.numeric(run$p.value)   else NA_real_
  #   zcrit <- stats::qnorm(1 - alpha/2)  # two-sided
  #   
  #   # 7) (Optional) Anderson–Darling for normality (if nortest installed)
  #   ad <- if (requireNamespace("nortest", quietly = TRUE))
  #     tryCatch(nortest::ad.test(res), error = function(e) NULL) else NULL
  #   ad_A2 <- if (!is.null(ad)) as.numeric(ad$statistic) else NA_real_
  #   ad_p  <- if (!is.null(ad)) as.numeric(ad$p.value)   else NA_real_
  #   
  #   # ---- Compose the academic-friendly text ----
  #   out <- c(
  #     "==========================================================================",
  #     "                     RESIDUAL DIAGNOSTIC BATTERY                          ",
  #     "==========================================================================",
  #     sprintf(" SAMPLE SIZE (residuals used): %d   |   α: %s   |   Lag (L): %d   |   fitdf: %d",
  #             N, fmt_num(alpha, 4), L, fitdf),
  #     "--------------------------------------------------------------------------",
  #     
  #     # Ljung–Box
  #     "TEST 1: Ljung–Box Portmanteau (autocorrelation)",
  #     " Purpose : Detect remaining serial correlation up to lag L in residuals.",
  #     sprintf(" Statistic : Q(LB) = %s  |  df = max(L - fitdf, 1) = %d  |  p-value = %s",
  #             fmt_num(lb_stat, 4), df_lb, fmt_p(lb_p)),
  #     sprintf(" Critical  : χ^2_(%d, 1-α) = %s", df_lb, fmt_num(cv_lb, 4)),
  #     sprintf(" Decision  : Reject H0 (white noise) if Q(LB) > χ^2_(%d,1-α) (equivalently p < α).", df_lb),
  #     sprintf(" Result    : %s",
  #             if (!is.finite(lb_p)) {
  #               "Inconclusive (statistic/p-value unavailable)."
  #             } else if (lb_p < alpha) {
  #               "Reject H0 → residuals show autocorrelation; model may be underspecified."
  #             } else {
  #               "Fail to reject H0 → residuals are consistent with white noise."
  #             }),
  #     "--------------------------------------------------------------------------",
  #     
  #     # Box–Pierce
  #     "TEST 2: Box–Pierce Portmanteau (autocorrelation, classic)",
  #     " Purpose : Historical portmanteau alternative to Ljung–Box.",
  #     sprintf(" Statistic : Q(BP) = %s  |  df = %d  |  p-value = %s",
  #             fmt_num(bp_stat, 4), df_lb, fmt_p(bp_p)),
  #     sprintf(" Critical  : χ^2_(%d, 1-α) = %s", df_lb, fmt_num(cv_lb, 4)),
  #     sprintf(" Decision  : Reject H0 if Q(BP) > χ^2_(%d,1-α) (equivalently p < α).", df_lb),
  #     sprintf(" Result    : %s",
  #             if (!is.finite(bp_p)) {
  #               "Inconclusive (statistic/p-value unavailable)."
  #             } else if (bp_p < alpha) {
  #               "Reject H0 → residuals show autocorrelation; consider increasing AR/MA orders or seasonal terms."
  #             } else {
  #               "Fail to reject H0 → no strong evidence of residual autocorrelation."
  #             }),
  #     "--------------------------------------------------------------------------",
  #     
  #     # Jarque–Bera
  #     "TEST 3: Jarque–Bera (normality, large-sample χ^2)",
  #     " Purpose : Assess normality via skewness and kurtosis of residuals.",
  #     sprintf(" Statistic : JB = %s  |  df = 2  |  p-value = %s", fmt_num(jb_stat, 4), fmt_p(jb_p)),
  #     sprintf(" Critical  : χ^2_(2, 1-α) = %s", fmt_num(cv_jb, 4)),
  #     " Decision  : Reject H0 (normality) if JB > χ^2_(2,1-α) (equivalently p < α).",
  #     sprintf(" Result    : %s",
  #             if (!is.finite(jb_p)) {
  #               "Inconclusive (statistic/p-value unavailable)."
  #             } else if (jb_p < alpha) {
  #               "Reject H0 → residuals deviate from normality; prediction intervals may be miscalibrated."
  #             } else {
  #               "Fail to reject H0 → residual normality is plausible."
  #             }),
  #     "--------------------------------------------------------------------------",
  #     
  #     # Shapiro–Wilk
  #     "TEST 4: Shapiro–Wilk (normality, small/medium samples)",
  #     " Purpose : Sensitive test for normality (recommended for N ≤ 5000).",
  #     sprintf(" Statistic : W = %s  |  p-value = %s  |  Note: critical values are not printed by R.",
  #             fmt_num(sw_W, 4), if (N > 5000) "n/a (N>5000)" else fmt_p(sw_p)),
  #     " Critical  : n/a (decision based on p-value).",
  #     " Decision  : Reject H0 (normality) if p < α.",
  #     sprintf(" Result    : %s",
  #             if (N > 5000) {
  #               "Omitted (Shapiro–Wilk only defined for N ≤ 5000)."
  #             } else if (!is.finite(sw_p)) {
  #               "Inconclusive (statistic/p-value unavailable)."
  #             } else if (sw_p < alpha) {
  #               "Reject H0 → residuals deviate from normality."
  #             } else {
  #               "Fail to reject H0 → residual normality is plausible."
  #             }),
  #     "--------------------------------------------------------------------------",
  #     
  #     # ARCH LM
  #     "TEST 5: Engle’s ARCH LM (conditional heteroskedasticity)",
  #     sprintf(" Purpose : Detect ARCH effects (time-varying variance) up to %d lags.", arch_lags),
  #     sprintf(" Statistic : LM = %s  |  df = %d  |  p-value = %s",
  #             fmt_num(arch_stat, 4), arch_lags, fmt_p(arch_p)),
  #     sprintf(" Critical  : χ^2_(%d, 1-α) = %s", arch_lags, fmt_num(cv_arch, 4)),
  #     " Decision  : Reject H0 (no ARCH) if LM > χ^2_(lags,1-α) (equivalently p < α).",
  #     sprintf(" Result    : %s",
  #             if (is.null(arch)) {
  #               "Skipped (package 'FinTS' not installed)."
  #             } else if (!is.finite(arch_p)) {
  #               "Inconclusive (statistic/p-value unavailable)."
  #             } else if (arch_p < alpha) {
  #               "Reject H0 → ARCH present; consider SARIMA + GARCH if volatility clustering matters."
  #             } else {
  #               "Fail to reject H0 → no strong evidence of ARCH."
  #             }),
  #     "--------------------------------------------------------------------------",
  #     
  #     # Runs test
  #     "TEST 6: Runs Test (randomness of signs)",
  #     " Purpose : Check independence/randomness of residual signs.",
  #     sprintf(" Statistic : Z = %s  |  p-value = %s  |  Two-sided z-crit = ±%s",
  #             fmt_num(run_Z, 4), fmt_p(run_p), fmt_num(zcrit, 3)),
  #     " Critical  : Reject H0 if |Z| > z_crit (equivalently p < α).",
  #     sprintf(" Result    : %s",
  #             if (is.null(run)) {
  #               "Skipped (package 'tseries' not installed)."
  #             } else if (!is.finite(run_p)) {
  #               "Inconclusive (statistic/p-value unavailable)."
  #             } else if (run_p < alpha) {
  #               "Reject H0 → residual signs are not random (possible structure left)."
  #             } else {
  #               "Fail to reject H0 → residual signs appear random."
  #             }),
  #     "==========================================================================",
  #     "INTERPRETATION GUIDE",
  #     " - Good SARIMA residuals typically: pass Ljung–Box (no autocorrelation),",
  #     "   show no strong ARCH (unless volatility modeling is intended), and are",
  #     "   roughly normal (for well-calibrated prediction intervals)."
  #   )
  #   
  #   paste(out, collapse = "\n")
  # }
  
  
  
  
  
  
  
  
  
  
  
  
  #===============================================================
  #===============================================================
  #===============================================================
  #===============================================================
  
  
  output$auto_forecast_plot <- renderPlot({
    req(auto_fc(), ts_train_test(), prepared())
    s <- ts_train_test()
    p <- prepared()
    fc <- auto_fc()$fc

    obs_df <- s$dfm[, c("x", "y_trans")]
    names(obs_df) <- c("x", "y")

    fc_df <- plot_forecast_df(obs_df, s$train_n, fc, by = p$by)
    gg_forecast_plot(obs_df, s$train_n, fc_df, title = "Auto-ARIMA forecast (train/test + intervals)")
  })

  output$auto_forecast_table <- renderTable({ req(auto_fc()); head(forecast_table(auto_fc()$fc), 25) }, rownames = FALSE)

  output$auto_accuracy_table <- renderTable({
    req(auto_fc(), ts_train_test())
    s <- ts_train_test()
    if (s$test_n == 0) return(data.frame(message = "No test set (training = 100%). Reduce training to compute accuracy."))
    accuracy_df(s$ts_test, auto_fc()$fc$mean)
  }, rownames = FALSE)

  output$apa_auto_paragraph <- renderPrint({
    req(auto_fit(), auto_fc(), ts_train_test())
    fit <- auto_fit()
    s <- ts_train_test()
    fc <- auto_fc()$fc
    lag <- as.numeric(input$diag_lag)
    lb <- tryCatch(Box.test(residuals(fit), lag = lag, type = "Ljung-Box", fitdf = length(coef(fit))), error = function(e) NULL)
    acc_line <- ""
    if (s$test_n > 0) {
      acc <- accuracy_df(s$ts_test, fc$mean)
      rmse <- acc$Value[acc$Metric == "RMSE"]
      mae <- acc$Value[acc$Metric == "MAE"]
      acc_line <- paste0("Forecast accuracy on the holdout set was RMSE = ", fmt_num(rmse, 2), " and MAE = ", fmt_num(mae, 2), ". ")
    }
    lb_line <- if (!is.null(lb)) paste0("The Ljung–Box test suggested ", ifelse(lb$p.value > 0.05, "no strong residual autocorrelation", "residual autocorrelation"), " (", fmt_p(lb$p.value), "). ") else ""
    cat(
      "APA-ready paragraph:\n\n",
      "An Auto-ARIMA procedure selected a seasonal ARIMA model (", as.character(fit), "). ",
      lb_line, acc_line,
      "Forecasts were generated with prediction intervals to quantify uncertainty.\n",
      sep = ""
    )
  })

  # ---- Step 6 Manual SARIMA ----

  output$manual_split_text <- renderPrint({
    req(ts_train_test(), prepared())
    s <- ts_train_test()
    df <- s$dfm
    train_n <- s$train_n
    test_n <- s$test_n
    x <- df$x

    cat("Train/Test split (Manual SARIMA)\n")
    cat("- Training proportion:", fmt_pct(as.numeric(input$train_prop)), "\n")
    cat("- Training size (n):", train_n, "\n")
    cat("- Test size (n):", test_n, "\n")
    if (test_n > 0) {
      cat("- Validation forecast horizon (h):", test_n, "(forced)\n\n")
    } else {
      cat("- Validation forecast horizon (h): none (future mode)\n\n")
    }

    cat("- Training range:", as.character(x[1]), "to", as.character(x[train_n]), "\n")
    if (test_n > 0) cat("- Test range:", as.character(x[train_n + 1]), "to", as.character(x[length(x)]), "\n")
  })

  output$manual_split_plot <- renderPlot({
    req(ts_train_test(), prepared())
    p <- prepared()
    s <- ts_train_test()
    df <- s$dfm
    df$set <- ifelse(seq_len(nrow(df)) <= s$train_n, "Train", "Test")
    split_df <- data.frame(xint = df$x[s$train_n])

    ggplot(df, aes(x = x, y = y_trans, color = set)) +
      geom_line(linewidth = 0.9) +
      geom_vline(data = split_df, aes(xintercept = xint), linetype = "dashed") +
      theme_minimal() +
      labs(title = "Train/Test split used for Manual SARIMA", x = p$x_label, y = "Value", color = NULL) +
      theme(legend.position = "bottom")
  })

  manual_fit <- eventReactive(input$fit_manual, {
    s <- ts_train_test()
    p <- prepared()
    period <- if (!is.na(input$s) && input$s >= 1) as.numeric(input$s) else p$freq
    forecast::Arima(
      s$ts_train,
      order = c(as.numeric(input$p), as.numeric(input$d), as.numeric(input$q)),
      seasonal = list(order = c(as.numeric(input$P), as.numeric(input$D), as.numeric(input$Q)), period = period),
      include.drift = isTRUE(input$manual_drift),
      include.mean = isTRUE(input$manual_drift)
    )
  })

  manual_fc <- reactive({
    req(manual_fit(), ts_train_test(), prepared())
    s <- ts_train_test()
    p <- prepared()

    # FIX v2 (core issue):
    # When a test set exists, we forecast exactly h = test_n so forecasts overlay the test period.
    user_h <- if (!is.na(input$manual_h) && input$manual_h > 0) as.numeric(input$manual_h) else NA_real_
    if (s$test_n > 0) {
      h <- s$test_n
    } else {
      h <- if (!is.na(user_h)) user_h else max(1, frequency(s$ts_train))
    }

    fc <- forecast::forecast(manual_fit(), h = h)
    list(fc = fc, h = h, by = p$by)
  })

  output$manual_horizon_note <- renderPrint({
    req(ts_train_test(), manual_fc())
    s <- ts_train_test()
    if (s$test_n > 0) {
      cat("Validation mode: forecast horizon is forced to the test length (h =", s$test_n, ") so the forecast overlays the test period.\n")
    } else {
      cat("Future mode: no test set. Horizon uses your 'Future horizon h' input (or defaults to one seasonal cycle).\n")
    }
  })

  output$manual_model_spec <- renderPrint({
    req(manual_fit(), prepared())
    p <- prepared()
    cat("Manual SARIMA model:\n")
    cat("- Non-seasonal (p,d,q): (", input$p, ",", input$d, ",", input$q, ")\n", sep = "")
    cat("- Seasonal (P,D,Q)[s]: (", input$P, ",", input$D, ",", input$Q, ")[", ifelse(is.na(input$s), p$freq, input$s), "]\n", sep = "")
    cat("- Drift/mean included:", isTRUE(input$manual_drift), "\n")
    cat("- AICc =", fmt_num(manual_fit()$aicc, 2), ", BIC =", fmt_num(manual_fit()$bic, 2), "\n")
  })

  output$manual_coef_table <- renderTable({ req(manual_fit()); coef_table(manual_fit()) }, rownames = FALSE)
  output$manual_resid_ts <- renderPlot({ req(manual_fit()); plot(residuals(manual_fit()), main = "Residuals (Manual SARIMA)", ylab = "Residual", xlab = "Time") })
  output$manual_resid_acf <- renderPlot({ req(manual_fit()); plot(acf(residuals(manual_fit()), plot = FALSE), main = "Residual ACF (Manual SARIMA)") })
  output$manual_resid_hist <- renderPlot({ req(manual_fit()); hist(residuals(manual_fit()), breaks = 30, main = "Residual histogram", xlab = "Residual") })
  output$manual_resid_qq <- renderPlot({ req(manual_fit()); qqnorm(residuals(manual_fit())); qqline(residuals(manual_fit())) })
  
  
  
  # --------------------------------------------- 
  # --------------------------------------------- 
  
  # output$manual_diag_tests <- renderText({ req(manual_fit()); diag_tests_text(residuals(manual_fit()), lag = as.numeric(input$diag_lag), fitdf = length(coef(manual_fit()))) })

  # Manual SARIMA residual tests (formatted)
  output$manual_diag_tests <- renderPrint({
    req(manual_fit())
    cat(
      residual_diagnostics_report_full(
        res   = residuals(manual_fit()),
        L     = suppressWarnings(as.integer(input$diag_lag %||% 12)),
        fitdf = length(coef(manual_fit())),
        alpha = to_num_safe(input$alphaSt2 %||% 0.05, 0.05)
      )
    )
  })
  
  
  residual_diagnostics_report_full <- function(res, L = 12L, fitdf = 0L, alpha = 0.05) {
    # sanitize
    res   <- as.numeric(res); res <- res[is.finite(res)]
    N     <- length(res)
    L     <- max(1L, as.integer(L))
    fitdf <- max(0L, as.integer(fitdf))
    alpha <- to_num_safe(alpha, 0.05)
    
    if (N < 8) {
      return(paste(
        "==========================================================================",
        "                    RESIDUAL DIAGNOSTIC BATTERY                           ",
        "==========================================================================",
        " Too few residuals (N < 8) to run the requested tests.",
        sep = "\n"
      ))
    }
    
    # convenience numbers
    df_lb  <- max(L - fitdf, 1L)
    cv_lb  <- stats::qchisq(1 - alpha, df = df_lb)
    cv_jb  <- stats::qchisq(1 - alpha, df = 2)
    arch_m <- min(L, max(1L, floor(N / 10)))
    zcrit  <- stats::qnorm(1 - alpha/2)
    
    # tests
    lb   <- tryCatch(stats::Box.test(res, lag = L, type = "Ljung-Box", fitdf = fitdf), error = function(e) NULL)
    bp   <- tryCatch(stats::Box.test(res, lag = L, type = "Box-Pierce", fitdf = fitdf), error = function(e) NULL)
    jb   <- if (requireNamespace("tseries", quietly = TRUE)) tryCatch(tseries::jarque.bera.test(res), error = function(e) NULL) else NULL
    sw   <- if (N >= 3 && N <= 5000) tryCatch(stats::shapiro.test(res), error = function(e) NULL) else NULL
    arch <- if (requireNamespace("FinTS", quietly = TRUE)) tryCatch(FinTS::ArchTest(res, lags = arch_m), error = function(e) NULL) else NULL
    run  <- if (requireNamespace("tseries", quietly = TRUE)) tryCatch(tseries::runs.test(res), error = function(e) NULL) else NULL
    ad   <- if (requireNamespace("nortest", quietly = TRUE)) tryCatch(nortest::ad.test(res), error = function(e) NULL) else NULL
    
    # pull numbers safely
    lb_stat  <- if (!is.null(lb))  as.numeric(lb$statistic) else NA_real_
    lb_p     <- if (!is.null(lb))  as.numeric(lb$p.value)   else NA_real_
    bp_stat  <- if (!is.null(bp))  as.numeric(bp$statistic) else NA_real_
    bp_p     <- if (!is.null(bp))  as.numeric(bp$p.value)   else NA_real_
    jb_stat  <- if (!is.null(jb))  as.numeric(jb$statistic) else NA_real_
    jb_p     <- if (!is.null(jb))  as.numeric(jb$p.value)   else NA_real_
    sw_W     <- if (!is.null(sw))  as.numeric(sw$statistic) else NA_real_
    sw_p     <- if (!is.null(sw))  as.numeric(sw$p.value)   else NA_real_
    arch_stat<- if (!is.null(arch))as.numeric(arch$statistic) else NA_real_
    arch_p   <- if (!is.null(arch))as.numeric(arch$p.value)    else NA_real_
    run_Z    <- if (!is.null(run)) as.numeric(run$statistic)  else NA_real_
    run_p    <- if (!is.null(run)) as.numeric(run$p.value)    else NA_real_
    ad_A2    <- if (!is.null(ad))  as.numeric(ad$statistic)   else NA_real_
    ad_p     <- if (!is.null(ad))  as.numeric(ad$p.value)     else NA_real_
    
    # flags for the overall conclusion
    lb_ok   <- is.finite(lb_p)   && lb_p   >= alpha
    bp_ok   <- is.finite(bp_p)   && bp_p   >= alpha
    norm_ok <- (is.finite(jb_p) && jb_p >= alpha) || (is.finite(sw_p) && sw_p >= alpha) || (is.finite(ad_p) && ad_p >= alpha)
    arch_ok <- is.null(arch) || (is.finite(arch_p) && arch_p >= alpha)
    runs_ok <- is.null(run)  || (is.finite(run_p)  && run_p  >= alpha)
    
    out <- c(
      "==========================================================================",
      "                     RESIDUAL DIAGNOSTIC BATTERY                          ",
      "==========================================================================",
      sprintf(" SAMPLE SIZE (residuals): %d   |   α: %s   |   Lag (L): %d   |   fitdf: %d",
              N, fmt_num(alpha, 4), L, fitdf),
      "--------------------------------------------------------------------------",
      
      # TEST 1: Ljung–Box
      "TEST 1 — LJUNG–BOX PORTMANTEAU (AUTOCORRELATION)",
      " Purpose: Detect remaining serial correlation up to lag L in the residuals of the fitted model.",
      " Description: The Ljung–Box statistic sums squared sample autocorrelations with a small-sample",
      "  correction. After estimating ARMA/SARIMA parameters, the degrees of freedom are reduced by the",
      "  number of fitted coefficients (fitdf). Well-specified residuals should resemble white noise.",
      " • H0: Residuals are white noise (no serial correlation up to lag L).",
      " • Ha: Residuals are autocorrelated (model may be underspecified).",
      sprintf(" → CRITERIA: Reject H0 if Q(LB) > χ^2_(%d,1-α)  (equivalently p-value < α).", df_lb),
      " RESULT:",
      sprintf("  - Q(LB)          : %s", fmt_num(lb_stat, 4)),
      sprintf("  - df (L - fitdf) : %d", df_lb),
      sprintf("  - χ^2 crit       : %s", fmt_num(cv_lb, 4)),
      sprintf("  - p-value        : %s", fmt_p(lb_p)),
      " DECISION & INTERPRETATION:",
      if (!is.finite(lb_p)) {
        "  Inconclusive: statistic or p-value unavailable."
      } else if (lb_p < alpha) {
        "  Reject H0. The residuals still contain autocorrelation up to lag L. This suggests the model may be\n  missing AR/MA or seasonal terms, or that differencing is insufficient. Review ACF/PACF of residuals and\n  consider adjusting orders, seasonal components, or transformations."
      } else {
        "  Fail to reject H0. Residuals behave like white noise up to lag L. This supports the adequacy of the model’s\n  dynamic specification (orders), which is desirable before interpreting parameters or forecasting."
      },
      "--------------------------------------------------------------------------",
      
      # TEST 2: Box–Pierce
      "TEST 2 — BOX–PIERCE PORTMANTEAU (CLASSIC AUTOCORRELATION)",
      " Purpose: Older portmanteau test for residual autocorrelation; conceptually similar to Ljung–Box.",
      " Description: Uses a simpler large-sample approximation without the Ljung–Box small-sample correction.",
      "  It is less accurate in small samples but should broadly agree with Ljung–Box when N is moderate/large.",
      " • H0: Residuals are white noise.",
      " • Ha: Residuals are autocorrelated.",
      sprintf(" → CRITERIA: Reject H0 if Q(BP) > χ^2_(%d,1-α)  (equivalently p-value < α).", df_lb),
      " RESULT:",
      sprintf("  - Q(BP)    : %s", fmt_num(bp_stat, 4)),
      sprintf("  - df       : %d", df_lb),
      sprintf("  - χ^2 crit : %s", fmt_num(cv_lb, 4)),
      sprintf("  - p-value  : %s", fmt_p(bp_p)),
      " DECISION & INTERPRETATION:",
      if (!is.finite(bp_p)) {
        "  Inconclusive: statistic or p-value unavailable."
      } else if (bp_p < alpha) {
        "  Reject H0. Autocorrelation remains. Combined with Ljung–Box, this strengthens the case for revising the model."
      } else {
        "  Fail to reject H0. No strong evidence of residual autocorrelation by Box–Pierce."
      },
      "--------------------------------------------------------------------------",
      
      # TEST 3: Jarque–Bera
      "TEST 3 — JARQUE–BERA (NORMALITY: SKEWNESS & KURTOSIS)",
      " Purpose: Evaluate whether residuals are approximately normally distributed by combining deviations in skewness",
      "  and kurtosis. Normal residuals help ensure well-calibrated prediction intervals and valid t-statistics in",
      "  regression-type outputs.",
      " Description: Asymptotically follows χ^2 with 2 df. Sensitive to heavy tails and skew.",
      " • H0: Residuals are normally distributed.",
      " • Ha: Residuals are not normal (skewed and/or heavy/light tails).",
      " → CRITERIA: Reject H0 if JB > χ^2_(2,1-α) (equivalently p-value < α).",
      " RESULT:",
      sprintf("  - JB statistic : %s", fmt_num(jb_stat, 4)),
      sprintf("  - χ^2 crit    : %s", fmt_num(cv_jb, 4)),
      sprintf("  - p-value     : %s", fmt_p(jb_p)),
      " DECISION & INTERPRETATION:",
      if (is.null(jb)) {
        "  Skipped: package 'tseries' not installed."
      } else if (!is.finite(jb_p)) {
        "  Inconclusive: statistic or p-value unavailable."
      } else if (jb_p < alpha) {
        "  Reject H0. Residuals deviate from normality. Forecast means are still unbiased if the model is correct, but\n  interval forecasts may be miscalibrated. Consider transformations (e.g., log/Box–Cox), robust modeling, or\n  heavy-tailed error models (e.g., t innovations) if this materially affects your goals."
      } else {
        "  Fail to reject H0. Normality is plausible by JB, supporting standard interval calibration."
      },
      "--------------------------------------------------------------------------",
      
      # TEST 4: Shapiro–Wilk
      "TEST 4 — SHAPIRO–WILK (NORMALITY: SMALL/MEDIUM N)",
      " Purpose: Sensitive test for normality, recommended when sample size is small to moderate.",
      " Description: Based on correlation between ordered sample values and corresponding normal scores. Does not print a simple",
      "  critical value in base R; decisions are p-value based. Defined for 3 ≤ N ≤ 5000.",
      " • H0: Residuals come from a normal distribution.",
      " • Ha: Residuals are non-normal.",
      " → CRITERIA: Reject H0 if p-value < α.",
      " RESULT:",
      sprintf("  - W statistic : %s", fmt_num(sw_W, 4)),
      sprintf("  - p-value     : %s", if (N > 5000) "n/a (N > 5000)" else fmt_p(sw_p)),
      " DECISION & INTERPRETATION:",
      if (N > 5000) {
        "  Omitted: Shapiro–Wilk is defined only up to N = 5000. For large N, prefer visual tools (QQ plot) and JB."
      } else if (is.null(sw)) {
        "  Inconclusive: Shapiro–Wilk did not run."
      } else if (!is.finite(sw_p)) {
        "  Inconclusive: statistic or p-value unavailable."
      } else if (sw_p < alpha) {
        "  Reject H0. Residuals depart from normality. Inspect QQ plot to see whether tails or skew drive the result; the remedy\n  depends on whether the distribution is heavy-tailed, skewed, or affected by outliers/level shifts."
      } else {
        "  Fail to reject H0. Shapiro–Wilk supports approximate normality."
      },
      "--------------------------------------------------------------------------",
      
      # TEST 5: Engle ARCH LM
      "TEST 5 — ENGLE'S ARCH LM (TIME-VARYING VARIANCE)",
      sprintf(" Purpose: Detect ARCH effects (conditional heteroskedasticity) up to %d lags. If present, variance clusters over time,", arch_m),
      "  which violates the constant-variance assumption and can distort interval forecasts.",
      " Description: Regress squared residuals on their own lags; under H0, the LM statistic ~ χ^2 with degrees of freedom equal",
      "  to the number of lags. Often used to motivate GARCH-type extensions.",
      " • H0: No ARCH effects (variance is constant over time).",
      " • Ha: ARCH effects present (variance changes with time).",
      sprintf(" → CRITERIA: Reject H0 if LM > χ^2_(%d,1-α) (equivalently p-value < α).", arch_m),
      " RESULT:",
      sprintf("  - LM statistic : %s", fmt_num(arch_stat, 4)),
      sprintf("  - df (lags)    : %d", arch_m),
      sprintf("  - χ^2 crit     : %s", fmt_num(stats::qchisq(1 - alpha, df = arch_m), 4)),
      sprintf("  - p-value      : %s", fmt_p(arch_p)),
      " DECISION & INTERPRETATION:",
      if (is.null(arch)) {
        "  Skipped: package 'FinTS' not installed. If volatility clustering is suspected (e.g., financial data), consider installing it."
      } else if (!is.finite(arch_p)) {
        "  Inconclusive: statistic or p-value unavailable."
      } else if (arch_p < alpha) {
        "  Reject H0. Evidence of time-varying variance. Consider SARIMA + GARCH (or other conditional variance models) if volatility matters\n  for your application; otherwise, robust intervals or transformations can help."
      } else {
        "  Fail to reject H0. No strong evidence of ARCH effects; the constant variance assumption is reasonable."
      },
      "--------------------------------------------------------------------------",
      
      # TEST 6: Runs test
      "TEST 6 — RUNS TEST (RANDOMNESS OF SIGNS)",
      " Purpose: Check whether residual signs (+/−) appear in random order. Non-random runs can indicate leftover structure",
      "  (e.g., bias, level shifts) even when autocorrelations are small.",
      " Description: Based on the number of sign runs compared with its expectation under randomness; large |Z| rejects randomness.",
      " • H0: Residual signs occur in random order (independent signs).",
      " • Ha: Residual signs are not random (patterns/clustering of signs).",
      sprintf(" → CRITERIA: Reject H0 if |Z| > %s  (equivalently p-value < α).", fmt_num(zcrit, 3)),
      " RESULT:",
      sprintf("  - Z statistic : %s", fmt_num(run_Z, 4)),
      sprintf("  - p-value     : %s", fmt_p(run_p)),
      " DECISION & INTERPRETATION:",
      if (is.null(run)) {
        "  Skipped: package 'tseries' not installed."
      } else if (!is.finite(run_p)) {
        "  Inconclusive: statistic or p-value unavailable."
      } else if (run_p < alpha) {
        "  Reject H0. Residual signs are not random, hinting at systematic under- or over-prediction or structural change.\n  Inspect fitted vs. actual plots and consider adding trend/seasonal terms or handling breaks/outliers."
      } else {
        "  Fail to reject H0. Residual signs appear random, which is consistent with an unbiased model."
      },
      "--------------------------------------------------------------------------",
      
      # TEST 7: Anderson–Darling
      "TEST 7 — ANDERSON–DARLING (NORMALITY, TAIL-SENSITIVE)",
      " Purpose: Additional normality check that gives more weight to the tails than Shapiro–Wilk/Jarque–Bera.",
      " Description: Often more sensitive to deviations in the extremes; useful when tail behavior matters for prediction intervals.",
      " • H0: Residuals are normally distributed.",
      " • Ha: Residuals are not normal.",
      " → CRITERIA: Reject H0 if p-value < α.",
      " RESULT:",
      sprintf("  - A^2 statistic : %s", fmt_num(ad_A2, 4)),
      sprintf("  - p-value       : %s", fmt_p(ad_p)),
      " DECISION & INTERPRETATION:",
      if (is.null(ad)) {
        "  Skipped: package 'nortest' not installed."
      } else if (!is.finite(ad_p)) {
        "  Inconclusive: statistic or p-value unavailable."
      } else if (ad_p < alpha) {
        "  Reject H0. Tail behavior deviates from normality; extreme forecast errors may be more common than normal theory predicts."
      } else {
        "  Fail to reject H0. Normal tails are plausible by Anderson–Darling."
      },
      
      "==========================================================================",
      "OVERALL CONCLUSION (HOW TO READ THESE TESTS TOGETHER)",
      if (!lb_ok || !bp_ok) {
        " Serial correlation is still present (one or both portmanteau tests rejected). Before relying on forecasts or coefficients,\n refine the SARIMA specification: check ACF/PACF of residuals, reconsider differencing (d/D), and try alternative AR/MA and\n seasonal orders until residual autocorrelation is removed."
      } else if (!arch_ok) {
        " Autocorrelation is under control, but variance likely varies over time (ARCH). For applications where interval accuracy matters,\n consider augmenting SARIMA with a volatility model (e.g., GARCH)."
      } else if (!norm_ok) {
        " Dynamics look adequate, but residuals deviate from normality. Forecast means remain useful; however, prediction intervals may\n need robust/bootstrapped methods, transformations, or heavy-tailed innovations."
      } else if (!runs_ok) {
        " No strong autocorrelation and variance looks stable, yet residual signs are not random. This can indicate bias or a structural\n feature not captured by the model. Inspect time plots/level shifts and consider adding trend/seasonal regressors or break handling."
      } else {
        " All core checks pass (no autocorrelation, no clear ARCH, normality plausible, signs random). Residual behavior is consistent with\n a well-specified SARIMA. Proceed to forecasting and keep monitoring residuals after re-estimation on new data."
      }
    )
    
    paste(out, collapse = "\n")
  }
  
  
  # --------------------------------------------- 
  # --------------------------------------------- 
  
  
  
  # ============================================================
  # --- MOD: Manual SARIMA equation renderer (FULL code) ---
  # ============================================================
  
  
  
  manual_equations <- reactive({
    req(manual_fit(), ts_train_test())
    
    fit   <- manual_fit()
    coefs <- coef(fit)
    
    # Orders
    p <- as.integer(input$p); d <- as.integer(input$d); q <- as.integer(input$q)
    P <- as.integer(input$P); D <- as.integer(input$D); Q <- as.integer(input$Q)
    
    # Seasonal period (never NA)
    s <- if (!is.na(input$s) && input$s >= 1) {
      as.integer(input$s)
    } else {
      f <- frequency(ts_train_test()$ts_train)
      if (is.na(f) || f < 1) 1L else as.integer(f)
    }
    
    ip <- if (p > 0) seq_len(p) else integer(0)
    iq <- if (q > 0) seq_len(q) else integer(0)
    iP <- if (P > 0) seq_len(P) else integer(0)
    iQ <- if (Q > 0) seq_len(Q) else integer(0)
    
    # Intercept/mean (show only if not ~0)
    intercept_name <- intersect(c("intercept", "mean"), names(coefs))
    intercept_val  <- if (length(intercept_name) > 0) unname(coefs[intercept_name[1]]) else 0
    show_intercept <- is.finite(intercept_val) && abs(intercept_val) > 1e-8
    intercept_num  <- if (show_intercept) sprintf("%.3f", intercept_val) else ""
    
    # Drift
    drift_sym <- if (isTRUE(input$manual_drift)) " + \\delta t" else ""
    drift_val <- if (isTRUE(input$manual_drift) && "drift" %in% names(coefs)) unname(coefs["drift"]) else NA_real_
    drift_num <- if (isTRUE(input$manual_drift) && is.finite(drift_val) && abs(drift_val) > 1e-8) {
      paste0(" + ", sprintf("%.3f", drift_val), "t")
    } else if (isTRUE(input$manual_drift)) {
      " + \\delta t"
    } else ""
    
    # MathJax display wrapper
    tex_display <- function(x) paste0("\\[", x, "\\]")
    
    # Cleanup for numeric line
    simplify_tex <- function(x) {
      x <- gsub("\\(1\\)", "", x)
      x <- gsub("\\s+", " ", x)
      x <- gsub("\\+\\s*\\+", "+", x)
      x <- gsub("\\+\\s*-", "-", x)
      x <- gsub("-\\s*\\+", "-", x)
      x <- gsub("-\\s*-", "+", x)
      x <- gsub("\\s*\\+\\s*0\\.000\\b", "", x)
      trimws(x)
    }
    
    # Parameter-polynomial strings (symbolic)
    poly_param_ar  <- function() if (p == 0) "1" else paste0("1", paste0(" - \\phi_{", ip, "}L^{", ip, "}", collapse = ""))
    poly_param_sar <- function() if (P == 0) "1" else paste0("1", paste0(" - \\Phi_{", iP, "}L^{", s * iP, "}", collapse = ""))
    poly_param_ma  <- function() if (q == 0) "1" else paste0("1", paste0(" + \\theta_{", iq, "}L^{", iq, "}", collapse = ""))
    poly_param_sma <- function() if (Q == 0) "1" else paste0("1", paste0(" + \\Theta_{", iQ, "}L^{", s * iQ, "}", collapse = ""))
    
    # Numeric-polynomial strings (from estimates)
    poly_num_ar <- function() {
      if (p == 0) return("1")
      nms <- paste0("ar", ip)
      v <- suppressWarnings(unname(coefs[nms]))
      keep <- is.finite(v)
      if (!any(keep)) return("1")
      paste0("1", paste(sprintf(" %+.3fL^{%d}", -v[keep], ip[keep]), collapse = ""))
    }
    poly_num_sar <- function() {
      if (P == 0) return("1")
      nms <- paste0("sar", iP)
      v <- suppressWarnings(unname(coefs[nms]))
      keep <- is.finite(v)
      if (!any(keep)) return("1")
      lags <- s * iP
      paste0("1", paste(sprintf(" %+.3fL^{%d}", -v[keep], lags[keep]), collapse = ""))
    }
    poly_num_ma <- function() {
      if (q == 0) return("1")
      nms <- paste0("ma", iq)
      v <- suppressWarnings(unname(coefs[nms]))
      keep <- is.finite(v)
      if (!any(keep)) return("1")
      paste0("1", paste(sprintf(" %+.3fL^{%d}", v[keep], iq[keep]), collapse = ""))
    }
    poly_num_sma <- function() {
      if (Q == 0) return("1")
      nms <- paste0("sma", iQ)
      v <- suppressWarnings(unname(coefs[nms]))
      keep <- is.finite(v)
      if (!any(keep)) return("1")
      lags <- s * iQ
      paste0("1", paste(sprintf(" %+.3fL^{%d}", v[keep], lags[keep]), collapse = ""))
    }
    
    # Differencing operators
    diff_part  <- if (d > 0) paste0("(1-L)^{", d, "}") else ""
    sdiff_part <- if (D > 0) paste0("(1-L^{", s, "})^{", D, "}") else ""
    
    # Line 1: General operator form
    line1 <- paste0(
      "\\phi_p(L)\\,\\Phi_P(L^{S})\\,(1-L)^{d}(1-L^{S})^{D}Y_t",
      " = c + \\theta_q(L)\\,\\Theta_Q(L^{S})\\varepsilon_t", drift_sym
    )
    
    # Line 2: Expanded operator (summation)
    line2 <- paste0(
      "\\left(1-\\sum_{i=1}^{p}\\phi_i L^{i}\\right)",
      "\\left(1-\\sum_{j=1}^{P}\\Phi_j L^{jS}\\right)",
      "(1-L)^{d}(1-L^{S})^{D}Y_t",
      " = c + ",
      "\\left(1+\\sum_{i=1}^{q}\\theta_i L^{i}\\right)",
      "\\left(1+\\sum_{j=1}^{Q}\\Theta_j L^{jS}\\right)",
      "\\varepsilon_t", drift_sym
    )
    
    # Line 3: parameter-expanded polynomials
    line3 <- paste0(
      "(", poly_param_ar(), ")",
      "(", poly_param_sar(), ")",
      diff_part, sdiff_part,
      "Y_t = c + ",
      "(", poly_param_ma(), ")",
      "(", poly_param_sma(), ")\\varepsilon_t",
      drift_sym
    )
    
    # Line 4: numeric-expanded polynomials
    rhs_intercept <- if (show_intercept) paste0(intercept_num, " + ") else ""
    line4 <- paste0(
      "(", poly_num_ar(), ")",
      "(", poly_num_sar(), ")",
      diff_part, sdiff_part,
      "Y_t = ",
      rhs_intercept,
      "(", poly_num_ma(), ")",
      "(", poly_num_sma(), ")\\varepsilon_t",
      drift_num
    )
    line4 <- simplify_tex(line4)
    
    # Time-domain (teaching form)
    ar_vals  <- if (p > 0) suppressWarnings(unname(coefs[paste0("ar", ip)])) else numeric(0)
    sar_vals <- if (P > 0) suppressWarnings(unname(coefs[paste0("sar", iP)])) else numeric(0)
    ma_vals  <- if (q > 0) suppressWarnings(unname(coefs[paste0("ma", iq)])) else numeric(0)
    sma_vals <- if (Q > 0) suppressWarnings(unname(coefs[paste0("sma", iQ)])) else numeric(0)
    
    keep_ar  <- is.finite(ar_vals)
    keep_sar <- is.finite(sar_vals)
    keep_ma  <- is.finite(ma_vals)
    keep_sma <- is.finite(sma_vals)
    
    td <- paste0(
      "Y_t = ",
      if (show_intercept) intercept_num else "0",
      if (p > 0 && any(keep_ar))  paste0(paste(sprintf(" %+.3fY_{t-%d}", ar_vals[keep_ar], ip[keep_ar]), collapse = "")) else "",
      if (P > 0 && any(keep_sar)) paste0(paste(sprintf(" %+.3fY_{t-%d}", sar_vals[keep_sar], (s * iP)[keep_sar]), collapse = "")) else "",
      if (q > 0 && any(keep_ma))  paste0(paste(sprintf(" %+.3f\\varepsilon_{t-%d}", ma_vals[keep_ma], iq[keep_ma]), collapse = "")) else "",
      if (Q > 0 && any(keep_sma)) paste0(paste(sprintf(" %+.3f\\varepsilon_{t-%d}", sma_vals[keep_sma], (s * iQ)[keep_sma]), collapse = "")) else "",
      " + \\varepsilon_t",
      if (d > 0 || D > 0) paste0(" \\\\ \\text{(with differencing: }(1-L)^{", d, "}(1-L^{", s, "})^{", D, "}\\text{)}") else ""
    )
    td <- simplify_tex(td)
    
    # ---------- Estimated coefficients (MathJax-friendly) ----------
    coef_lines <- c()
    
    if (show_intercept) {
      coef_lines <- c(coef_lines, paste0("\\(c\\) (intercept/mean) = ", sprintf("%.4f", intercept_val)))
    }
    if (isTRUE(input$manual_drift)) {
      if (is.finite(drift_val)) {
        coef_lines <- c(coef_lines, paste0("drift \\((\\delta)\\) = ", sprintf("%.4f", drift_val)))
      } else {
        coef_lines <- c(coef_lines, "drift \\((\\delta)\\) included (value not estimated explicitly)")
      }
    }
    
    if (p > 0) {
      for (i in ip) {
        nm <- paste0("ar", i)
        if (nm %in% names(coefs) && is.finite(coefs[nm])) {
          coef_lines <- c(coef_lines, paste0("ar", i, ": \\(\\phi_{", i, "}\\) = ", sprintf("%.4f", coefs[nm])))
        }
      }
    }
    if (q > 0) {
      for (i in iq) {
        nm <- paste0("ma", i)
        if (nm %in% names(coefs) && is.finite(coefs[nm])) {
          coef_lines <- c(coef_lines, paste0("ma", i, ": \\(\\theta_{", i, "}\\) = ", sprintf("%.4f", coefs[nm])))
        }
      }
    }
    if (P > 0) {
      for (i in iP) {
        nm <- paste0("sar", i)
        if (nm %in% names(coefs) && is.finite(coefs[nm])) {
          coef_lines <- c(coef_lines, paste0("sar", i, ": \\(\\Phi_{", i, "}\\) = ", sprintf("%.4f", coefs[nm])))
        }
      }
    }
    if (Q > 0) {
      for (i in iQ) {
        nm <- paste0("sma", i)
        if (nm %in% names(coefs) && is.finite(coefs[nm])) {
          coef_lines <- c(coef_lines, paste0("sma", i, ": \\(\\Theta_{", i, "}\\) = ", sprintf("%.4f", coefs[nm])))
        }
      }
    }
    
    if (length(coef_lines) == 0) coef_lines <- "No coefficients available."
    
    list(
      p = p, d = d, q = q, P = P, D = D, Q = Q, s = s,
      coef_lines   = coef_lines,
      eq_general   = tex_display(line1),
      eq_expanded  = tex_display(line2),
      eq_line3     = tex_display(line3),
      eq_line4     = tex_display(line4),
      eq_time_domain = tex_display(td)
    )
  })
  
  
  
  
  
  
  
  # --- MOD: Replace/Update your manual_equations() reactive with this FULL version ---
  # manual_equations <- reactive({
  #   req(manual_fit(), ts_train_test())
  #   
  #   fit <- manual_fit()
  #   coefs <- coef(fit)
  #   
  #   # Orders
  #   p <- as.integer(input$p); d <- as.integer(input$d); q <- as.integer(input$q)
  #   P <- as.integer(input$P); D <- as.integer(input$D); Q <- as.integer(input$Q)
  #   
  #   # Seasonal period (never NA)
  #   s <- if (!is.na(input$s) && input$s >= 1) {
  #     as.integer(input$s)
  #   } else {
  #     f <- frequency(ts_train_test()$ts_train)
  #     if (is.na(f) || f < 1) 1L else as.integer(f)
  #   }
  #   
  #   ip <- if (p > 0) seq_len(p) else integer(0)
  #   iq <- if (q > 0) seq_len(q) else integer(0)
  #   iP <- if (P > 0) seq_len(P) else integer(0)
  #   iQ <- if (Q > 0) seq_len(Q) else integer(0)
  #   
  #   # Intercept/mean (show only if not ~0)
  #   intercept_name <- intersect(c("intercept", "mean"), names(coefs))
  #   intercept_val <- if (length(intercept_name) > 0) unname(coefs[intercept_name[1]]) else 0
  #   show_intercept <- is.finite(intercept_val) && abs(intercept_val) > 1e-8
  #   intercept_num <- if (show_intercept) sprintf("%.3f", intercept_val) else ""
  #   
  #   # Drift
  #   drift_sym <- if (isTRUE(input$manual_drift)) " + \\delta t" else ""
  #   drift_val <- if (isTRUE(input$manual_drift) && "drift" %in% names(coefs)) unname(coefs["drift"]) else NA_real_
  #   drift_num <- if (isTRUE(input$manual_drift) && is.finite(drift_val) && abs(drift_val) > 1e-8) {
  #     paste0(" + ", sprintf("%.3f", drift_val), "t")
  #   } else if (isTRUE(input$manual_drift)) {
  #     " + \\delta t"
  #   } else ""
  #   
  #   # --- MOD: MathJax-safe equation blocks ---
  #   # Use \\[ ... \\] and align with aligned/array (supported).
  #   tex_display <- function(x) paste0("\\[", x, "\\]")
  #   
  #   # --- MOD: clean cosmetics in numeric equation line ---
  #   simplify_tex <- function(x) {
  #     x <- gsub("\\(1\\)", "", x)                      # remove (1)
  #     x <- gsub("\\s+", " ", x)                        # normalize spaces
  #     x <- gsub("\\+\\s*\\+", "+", x)                  # ++ -> +
  #     x <- gsub("\\+\\s*-", "-", x)                    # +- -> -
  #     x <- gsub("-\\s*\\+", "-", x)                    # -+ -> -
  #     x <- gsub("-\\s*-", "+", x)                      # -- -> +
  #     x <- gsub("\\s*\\+\\s*0\\.000\\b", "", x)        # remove + 0.000
  #     # x <- gsub("\\b0\\.000\\s*\\+\\s*", "", x)        # remove leading 0.000 +
  #     # x <- gsub("\\s*\\+\\s*0\\b", "", x)              # remove + 0
  #     # x <- gsub("\\b0\\s*\\+\\s*", "", x)              # remove leading 0 +
  #     trimws(x)
  #   }
  #   
  #   # Parameter-polynomial strings (line 3)
  #   poly_param_ar  <- function() if (p == 0) "1" else paste0("1", paste0(" - \\phi_{", ip, "}L^{", ip, "}", collapse = ""))
  #   poly_param_sar <- function() if (P == 0) "1" else paste0("1", paste0(" - \\Phi_{", iP, "}L^{", s * iP, "}", collapse = ""))
  #   poly_param_ma  <- function() if (q == 0) "1" else paste0("1", paste0(" + \\theta_{", iq, "}L^{", iq, "}", collapse = ""))
  #   poly_param_sma <- function() if (Q == 0) "1" else paste0("1", paste0(" + \\Theta_{", iQ, "}L^{", s * iQ, "}", collapse = ""))
  #   
  #   # Numeric-polynomial strings (line 4) with NA-safe filtering and correct signs
  #   poly_num_ar <- function() {
  #     if (p == 0) return("1")
  #     nms <- paste0("ar", ip)
  #     v <- suppressWarnings(unname(coefs[nms]))
  #     keep <- is.finite(v)
  #     if (!any(keep)) return("1")
  #     paste0("1", paste(sprintf(" %+.3fL^{%d}", -v[keep], ip[keep]), collapse = ""))
  #   }
  #   poly_num_sar <- function() {
  #     if (P == 0) return("1")
  #     nms <- paste0("sar", iP)
  #     v <- suppressWarnings(unname(coefs[nms]))
  #     keep <- is.finite(v)
  #     if (!any(keep)) return("1")
  #     lags <- s * iP
  #     paste0("1", paste(sprintf(" %+.3fL^{%d}", -v[keep], lags[keep]), collapse = ""))
  #   }
  #   poly_num_ma <- function() {
  #     if (q == 0) return("1")
  #     nms <- paste0("ma", iq)
  #     v <- suppressWarnings(unname(coefs[nms]))
  #     keep <- is.finite(v)
  #     if (!any(keep)) return("1")
  #     paste0("1", paste(sprintf(" %+.3fL^{%d}", v[keep], iq[keep]), collapse = ""))
  #   }
  #   poly_num_sma <- function() {
  #     if (Q == 0) return("1")
  #     nms <- paste0("sma", iQ)
  #     v <- suppressWarnings(unname(coefs[nms]))
  #     keep <- is.finite(v)
  #     if (!any(keep)) return("1")
  #     lags <- s * iQ
  #     paste0("1", paste(sprintf(" %+.3fL^{%d}", v[keep], lags[keep]), collapse = ""))
  #   }
  #   
  #   # Differencing operators
  #   diff_part  <- if (d > 0) paste0("(1-L)^{", d, "}") else ""
  #   sdiff_part <- if (D > 0) paste0("(1-L^{", s, "})^{", D, "}") else ""
  #   
  #   # ---------- Line 1: General operator form ----------
  #   line1 <- paste0(
  #     "\\phi_p(L)\\,\\Phi_P(L^{S})\\,(1-L)^{d}(1-L^{S})^{D}Y_t",
  #     " = c + \\theta_q(L)\\,\\Theta_Q(L^{S})\\varepsilon_t", drift_sym
  #   )
  #   
  #   # ---------- Line 2: Expanded operator (summation) ----------
  #   line2 <- paste0(
  #     "\\left(1-\\sum_{i=1}^{p}\\phi_i L^{i}\\right)",
  #     "\\left(1-\\sum_{j=1}^{P}\\Phi_j L^{jS}\\right)",
  #     "(1-L)^{d}(1-L^{S})^{D}Y_t",
  #     " = c + ",
  #     "\\left(1+\\sum_{i=1}^{q}\\theta_i L^{i}\\right)",
  #     "\\left(1+\\sum_{j=1}^{Q}\\Theta_j L^{jS}\\right)",
  #     "\\varepsilon_t", drift_sym
  #   )
  #   
  #   # ---------- Line 3: parameter-expanded polynomials ----------
  #   line3 <- paste0(
  #     "(", poly_param_ar(), ")",
  #     "(", poly_param_sar(), ")",
  #     diff_part, sdiff_part,
  #     "Y_t = c + ",
  #     "(", poly_param_ma(), ")",
  #     "(", poly_param_sma(), ")\\varepsilon_t",
  #     drift_sym
  #   )
  #   
  #   # ---------- Line 4: numeric-expanded polynomials ----------
  #   rhs_intercept <- if (show_intercept) paste0(intercept_num, " + ") else ""
  #   line4 <- paste0(
  #     "(", poly_num_ar(), ")",
  #     "(", poly_num_sar(), ")",
  #     diff_part, sdiff_part,
  #     "Y_t = ",
  #     rhs_intercept,
  #     "(", poly_num_ma(), ")",
  #     "(", poly_num_sma(), ")\\varepsilon_t",
  #     drift_num
  #   )
  #   line4 <- simplify_tex(line4)
  #   
  #   # ---------- Equivalent time-domain (teaching form) ----------
  #   ar_vals  <- if (p > 0) suppressWarnings(unname(coefs[paste0("ar", ip)])) else numeric(0)
  #   sar_vals <- if (P > 0) suppressWarnings(unname(coefs[paste0("sar", iP)])) else numeric(0)
  #   ma_vals  <- if (q > 0) suppressWarnings(unname(coefs[paste0("ma", iq)])) else numeric(0)
  #   sma_vals <- if (Q > 0) suppressWarnings(unname(coefs[paste0("sma", iQ)])) else numeric(0)
  #   
  #   keep_ar  <- is.finite(ar_vals)
  #   keep_sar <- is.finite(sar_vals)
  #   keep_ma  <- is.finite(ma_vals)
  #   keep_sma <- is.finite(sma_vals)
  #   
  #   td <- paste0(
  #     "Y_t = ",
  #     if (show_intercept) intercept_num else "0",
  #     if (p > 0 && any(keep_ar))  paste0(paste(sprintf(" %+.3fY_{t-%d}", ar_vals[keep_ar], ip[keep_ar]), collapse = "")) else "",
  #     if (P > 0 && any(keep_sar)) paste0(paste(sprintf(" %+.3fY_{t-%d}", sar_vals[keep_sar], (s * iP)[keep_sar]), collapse = "")) else "",
  #     if (q > 0 && any(keep_ma))  paste0(paste(sprintf(" %+.3f\\varepsilon_{t-%d}", ma_vals[keep_ma], iq[keep_ma]), collapse = "")) else "",
  #     if (Q > 0 && any(keep_sma)) paste0(paste(sprintf(" %+.3f\\varepsilon_{t-%d}", sma_vals[keep_sma], (s * iQ)[keep_sma]), collapse = "")) else "",
  #     " + \\varepsilon_t",
  #     if (d > 0 || D > 0) paste0(" \\\\ \\text{(with differencing: }(1-L)^{", d, "}(1-L^{", s, "})^{", D, "}\\text{)}") else ""
  #   )
  #   td <- simplify_tex(td)
  #   
  #   # --- MOD: coefficient listing, including ar/ma/sar/sma ---
  #   coef_lines <- c()
  #   
  #   if (show_intercept) coef_lines <- c(coef_lines, paste0("c (intercept/mean) = ", sprintf("%.4f", intercept_val)))
  #   if (isTRUE(input$manual_drift)) {
  #     if (is.finite(drift_val)) coef_lines <- c(coef_lines, paste0("drift = ", sprintf("%.4f", drift_val)))
  #     else coef_lines <- c(coef_lines, "drift included (value not estimated explicitly)")
  #   }
  #   
  #   if (p > 0) {
  #     for (i in ip) {
  #       nm <- paste0("ar", i)
  #       if (nm %in% names(coefs) && is.finite(coefs[nm])) {
  #         coef_lines <- c(coef_lines, paste0("ar", i, " (\\phi_", i, ") = ", sprintf("%.4f", coefs[nm])))
  #       }
  #     }
  #   }
  #   if (q > 0) {
  #     for (i in iq) {
  #       nm <- paste0("ma", i)
  #       if (nm %in% names(coefs) && is.finite(coefs[nm])) {
  #         coef_lines <- c(coef_lines, paste0("ma", i, " (\\theta_", i, ") = ", sprintf("%.4f", coefs[nm])))
  #       }
  #     }
  #   }
  #   if (P > 0) {
  #     for (i in iP) {
  #       nm <- paste0("sar", i)
  #       if (nm %in% names(coefs) && is.finite(coefs[nm])) {
  #         coef_lines <- c(coef_lines, paste0("sar", i, " (\\Phi_", i, ") = ", sprintf("%.4f", coefs[nm])))
  #       }
  #     }
  #   }
  #   if (Q > 0) {
  #     for (i in iQ) {
  #       nm <- paste0("sma", i)
  #       if (nm %in% names(coefs) && is.finite(coefs[nm])) {
  #         coef_lines <- c(coef_lines, paste0("sma", i, " (\\Theta_", i, ") = ", sprintf("%.4f", coefs[nm])))
  #       }
  #     }
  #   }
  #   
  #   if (length(coef_lines) == 0) coef_lines <- c("No coefficients available.")
  #   
  #   # --- MOD: package sections as MathJax display blocks ---
  #   list(
  #     p = p, d = d, q = q, P = P, D = D, Q = Q, s = s,
  #     coef_lines = coef_lines,
  #     eq_general = tex_display(line1),
  #     eq_expanded = tex_display(line2),
  #     eq_line3 = tex_display(line3),
  #     eq_line4 = tex_display(line4),
  #     eq_time_domain = tex_display(td)
  #   )
  # })
  
  # ============================================================
  # --- MOD: Render the equation panel with LEFT alignment (CSS) and headings ---
  #   NOTE: Left alignment is done via an HTML wrapper div.
  # ============================================================
  output$manual_model_equation <- renderUI({
    req(manual_equations())
    eq <- manual_equations()
    
    tagList(
      # --- MOD: left alignment wrapper for all MathJax blocks ---
      tags$div(
        style = "text-align:left;",
        
        tags$h4("Manual SARIMA model"),
        tags$p(sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)[%d]", eq$p, eq$d, eq$q, eq$P, eq$D, eq$Q, eq$s)),
        
        tags$h5("Estimated coefficients"),
        tags$ul(lapply(eq$coef_lines, function(x) tags$li(HTML(x)))),
        
        tags$hr(),
        
        tags$h4("General SARIMA formulation"),
        HTML(eq$eq_general),
        
        tags$hr(),
        
        tags$h4("Expanded operator form"),
        HTML(eq$eq_expanded),
        
        tags$hr(),
        
        tags$h4("Numerical model"),
        # --- MOD: show line 3 then line 4 under the same heading ---
        HTML(eq$eq_line3),
        tags$hr(),
        
        # HTML(tex_display("\\text{------------}")),
        HTML(eq$eq_line4),
        
        tags$hr(),
        
        # tags$h4("Equivalent time-domain representation (teaching form)"),
        # HTML(eq$eq_time_domain)
      ),
      
      # --- MOD: Force MathJax typesetting for dynamically inserted content ---
      tags$script(HTML("if (window.MathJax && MathJax.Hub) MathJax.Hub.Queue(['Typeset', MathJax.Hub]);"))
    )
  })
  
  # --- MOD: helper used above inside renderUI (place near other helpers or above output block) ---
  tex_display <- function(x) paste0("\\[", x, "\\]")
  
  
  

  # --------------------------------------------- 
  # --------------------------------------------- 
  
  output$manual_forecast_plot <- renderPlot({
    req(manual_fc(), ts_train_test(), prepared())
    s <- ts_train_test()
    p <- prepared()
    fc <- manual_fc()$fc

    obs_df <- s$dfm[, c("x", "y_trans")]
    names(obs_df) <- c("x", "y")

    fc_df <- plot_forecast_df(obs_df, s$train_n, fc, by = p$by)
    gg_forecast_plot(obs_df, s$train_n, fc_df, title = "Manual SARIMA forecast (train/test + intervals)")
  })
  
  

  output$manual_forecast_table <- renderTable({ req(manual_fc()); head(forecast_table(manual_fc()$fc), 25) }, rownames = FALSE)

  
  
  output$manual_accuracy_table <- renderTable({
    req(manual_fc(), ts_train_test())
    s <- ts_train_test()
    if (s$test_n == 0) return(data.frame(message = "No test set (training = 100%). Reduce training to compute accuracy."))
    accuracy_df(s$ts_test, manual_fc()$fc$mean)
  }, rownames = FALSE)

  output$apa_manual_paragraph <- renderPrint({
    req(manual_fit(), manual_fc(), ts_train_test())
    fit <- manual_fit()
    s <- ts_train_test()
    fc <- manual_fc()$fc
    lag <- as.numeric(input$diag_lag)
    lb <- tryCatch(Box.test(residuals(fit), lag = lag, type = "Ljung-Box", fitdf = length(coef(fit))), error = function(e) NULL)
    acc_line <- ""
    if (s$test_n > 0) {
      acc <- accuracy_df(s$ts_test, fc$mean)
      rmse <- acc$Value[acc$Metric == "RMSE"]
      mae <- acc$Value[acc$Metric == "MAE"]
      acc_line <- paste0("Forecast accuracy on the holdout set was RMSE = ", fmt_num(rmse, 2), " and MAE = ", fmt_num(mae, 2), ". ")
    }
    lb_line <- if (!is.null(lb)) paste0("The Ljung–Box test suggested ", ifelse(lb$p.value > 0.05, "no strong residual autocorrelation", "residual autocorrelation"), " (", fmt_p(lb$p.value), "). ") else ""
    cat(
      "APA-ready paragraph:\n\n",
      "A manual seasonal ARIMA model was specified as (", input$p, ",", input$d, ",", input$q, ")(",
      input$P, ",", input$D, ",", input$Q, ")[", ifelse(is.na(input$s), frequency(s$ts_train), input$s), "]. ",
      lb_line, acc_line,
      "Forecasts were produced with prediction intervals to quantify uncertainty.\n",
      sep = ""
    )
  })

  # ---- Step 7 Comparison & Paper builder ----

  comparison <- reactive({
    s <- ts_train_test()
    out <- data.frame(Model = character(0), AICc = numeric(0), BIC = numeric(0), Test_RMSE = numeric(0), Test_MAE = numeric(0), stringsAsFactors = FALSE)

    if (!is.null(input$fit_auto)) {
      fit <- tryCatch(auto_fit(), error = function(e) NULL)
      if (!is.null(fit)) {
        rmse <- NA_real_; mae <- NA_real_
        if (s$test_n > 0) {
          acc <- accuracy_df(s$ts_test, auto_fc()$fc$mean)
          rmse <- acc$Value[acc$Metric == "RMSE"]
          mae <- acc$Value[acc$Metric == "MAE"]
        }
        out <- rbind(out, data.frame(Model = "Auto-ARIMA", AICc = fit$aicc, BIC = fit$bic, Test_RMSE = rmse, Test_MAE = mae))
      }
    }

    if (!is.null(input$fit_manual)) {
      fit <- tryCatch(manual_fit(), error = function(e) NULL)
      if (!is.null(fit)) {
        rmse <- NA_real_; mae <- NA_real_
        if (s$test_n > 0) {
          acc <- accuracy_df(s$ts_test, manual_fc()$fc$mean)
          rmse <- acc$Value[acc$Metric == "RMSE"]
          mae <- acc$Value[acc$Metric == "MAE"]
        }
        out <- rbind(out, data.frame(Model = "Manual SARIMA", AICc = fit$aicc, BIC = fit$bic, Test_RMSE = rmse, Test_MAE = mae))
      }
    }
    out
  })

  output$comparison_table <- renderTable({
    df <- comparison()
    if (nrow(df) == 0) return(data.frame(message = "Fit Auto-ARIMA and/or Manual SARIMA to compare models."))
    df
  }, rownames = FALSE)

  output$comparison_interpretation <- renderPrint({
    df <- comparison()
    if (nrow(df) == 0) { cat("Fit models first to generate comparison interpretation.\n"); return() }
    cat("Interpretation guide:\n")
    cat("- Lower AICc/BIC indicates better in-sample trade-off (fit vs complexity).\n")
    cat("- Lower RMSE/MAE indicates better holdout forecast accuracy.\n\n")
    best <- df[order(df$Test_RMSE, df$AICc, na.last = TRUE), , drop = FALSE]
    cat("Top-ranked (by RMSE then AICc):", best$Model[1], "\n")
  })

  output$apa_methods_draft <- renderPrint({
    req(prepared())
    p <- prepared()
    df <- p$df
    n <- nrow(df)
    cat(
      "Methods (APA draft — edit names/unit and add study context):\n\n",
      "Time series analyses were conducted using R. The series consisted of ", n, " observations spanning ",
      format(min(df$date)), " to ", format(max(df$date)), " with a seasonal period of s = ", p$freq, ". ",
      "Missing values were handled using the ", input$missing_policy, " method. ",
      "Stationarity was assessed using the Augmented Dickey–Fuller (ADF), KPSS, and Phillips–Perron tests, and differencing decisions were informed by these tests alongside visual inspection of differenced series plots. ",
      "Seasonal ARIMA models were estimated using maximum likelihood as implemented in the forecast package.\n",
      sep = ""
    )
  })

  output$apa_results_draft <- renderPrint({
    cat(
      "Results (APA draft — refine after you finalize the model):\n\n",
      "Exploratory analysis suggested trend and seasonal dynamics in the series, motivating the evaluation of seasonal ARIMA models. ",
      "Model adequacy was evaluated via residual diagnostics, including the Ljung–Box test for residual autocorrelation and additional normality and heteroskedasticity checks when available. ",
      "Forecast performance was evaluated using holdout accuracy metrics (e.g., RMSE and MAE) when a test set was available. ",
      "The final model specification, diagnostics, and forecast results are reported in the corresponding panels.\n",
      sep = ""
    )
  })

  observeEvent(input$refresh_all, { invisible(NULL) })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # ====================================================================
  # ====================================================================
  # ====================================================================
  #
  #                          - D?d?log? -
  #
  # ====================================================================
  # ====================================================================
  # ==================================================================== 
  
  
  
  # helper used below (safe defaulting)
  # `%||%` <- function(a, b) if (is.null(a) || is.na(a) || identical(a, "")) b else a
  
  
  # ============================================================================
  # 0) SMALL HELPERS (safe input fallback + safe numeric)
  # ============================================================================
  `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
  to_num_safe <- function(v, default = NA_real_) {
    out <- suppressWarnings(as.numeric(v))
    if (length(out) == 0 || all(is.na(out)) || !is.finite(out[1])) default else out[1]
  }
  
  getPlotDim <- function(x, default = "100%") {
    if (is.null(x) || is.na(x) || identical(x, "")) return(default)
    if (is.numeric(x)) return(paste0(x, "px"))
    # accept strings like "800" or "100%" or "800px"
    xs <- as.character(x)
    if (grepl("^[0-9]+$", xs)) return(paste0(xs, "px"))
    xs
  }
  
  
  # # ---- helpers.R (or at the top of server.R) ----
  # # Convert numeric -> "Npx"; pass through "100%" or "auto"; provide a default.
  # getPlotDim <- function(x, default = "100%") {
  #   if (is.null(x)) return(default)
  #   if (is.na(x))   return(default)
  #   if (is.numeric(x)) return(paste0(x, "px"))
  #   # if it's already a character like "100%" or "450px", just return it
  #   as.character(x)
  # }
  
  
  # helper: map input$plot_theme -> a ggplot2 theme object
  theme_picker <- function(key = "Minimal") {
    switch(key,
           "Minimal" = ggplot2::theme_minimal(),
           "Classic" = ggplot2::theme_classic(),
           "Light"   = ggplot2::theme_light(),
           "Dark"    = ggplot2::theme_dark(),
           "BW"      = ggplot2::theme_bw(),
           "Void"    = ggplot2::theme_void(),
           # default
           ggplot2::theme_gray()
    )
  }
  
  # common theming applied to ggplot objects
  add_theme <- function(g) {
    g +
      theme_picker(input$plot_theme) +
      ggplot2::theme(
        axis.text  = ggplot2::element_text(size = input$tickSize),
        axis.title = ggplot2::element_text(size = input$tickSize + 2),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )
  }
  
  
  
  
  # Map the log checkbox to a simple flag
  values <- reactiveValues(islog = "No")
  observe({
    values$islog <- if (isTRUE(input$check_box)) "Yes" else "No"
  })
  
  # Core transformation helper (log, then D, then d)
  getMyData <- function(tsData, frequency, islog = "No", d_n = 0, DS_n = 0) {
    shiny::req(tsData)
    freq <- as.numeric(frequency)
    d_n  <- as.numeric(d_n)
    DS_n <- as.numeric(DS_n)
    
    working_ts <- if (identical(islog, "Yes")) log(tsData) else tsData
    
    # Seasonal differencing (D)
    if (!is.na(DS_n) && DS_n > 0 && !is.na(freq) && freq > 1) {
      working_ts <- diff(working_ts, lag = freq, differences = DS_n)
    }
    
    # Non-seasonal differencing (d)
    if (!is.na(d_n) && d_n > 0) {
      working_ts <- diff(working_ts, lag = 1, differences = d_n)
    }
    
    working_ts
  }
  
  # Base ts comes from your prepared() reactive
  # prepared() returns list(df=..., freq=..., by=..., x_label=...)
  # We'll use the filled values (before transform) to let the playground own the transforms
  ts_base <- reactive({
    p <- prepared()                                   # df/freq are defined here
    ts(p$df$y_filled, start = 1, frequency = p$freq)  # build a ts from the cleaned series
  })
  
  # Central reactive used by all (*) panels in this tab
  # ========= replace the whole myData_Choice() reactive with this =========
  myData_Choice <- reactive({
    req(prepared())
    base_ts <- ts_base()                 # full, cleaned & regularly spaced series (y_filled)
    
    # Figure out how many observations belong to the training window
    s <- ts_train_test()                 # has $train_n computed from input$train_prop
    train_n <- s$train_n
    
    # Choose full vs. training window based on the new checkbox
    if (isTRUE(input$use_train_explore) && is.finite(train_n) && train_n >= 2) {
      # keep start and frequency; truncate by index
      base_vec <- as.numeric(base_ts)
      base_vec <- base_vec[seq_len(min(train_n, length(base_vec)))]
      base_ts  <- ts(base_vec, start = start(base_ts), frequency = frequency(base_ts))
    }
    
    # Apply the Exploration controls (log / d / D) consistently to whatever window we chose
    getMyData(
      tsData    = base_ts,
      frequency = frequency(base_ts),
      islog     = values$islog,          # this already mirrors input$check_box via your observeEvent
      d_n       = input$d_n,
      DS_n      = input$DS_n
    )
  })
  
  
  
  
  
  
  # myData_Choice <- reactive({
  #   getMyData(
  #     tsData    = ts_base(),
  #     frequency = frequency(ts_base()),
  #     islog     = values$islog,
  #     d_n       = input$d_n,
  #     DS_n      = input$DS_n
  #   )
  # })
  
  # # --- UI wrapper for the combined diagnostics plot ---
  # output$d_D_Log_ts_Choice_UI <- renderUI({
  #   plotOutput("d_D_Log_ts_Choice", width = "100%", height = 420)
  # })
  
  # UI wrapper: uses textInput("plot_width"), textInput("plot_height")
  output$d_D_Log_ts_Choice_UI <- renderUI({
    plotOutput(
      "d_D_Log_ts_Choice",
      width  = getPlotDim(input$plot_width  %||% "800"),
      height = getPlotDim(input$plot_height %||% "500")
    )
  })
  
  
  # --- ggtsdisplay of the transformed series (time plot + ACF + PACF) ---
  output$d_D_Log_ts_Choice <- renderPlot({
    req(myData_Choice())
    p <- prepared()
    forecast::ggtsdisplay(
      myData_Choice(),
      main = "Diagnostics of transformed series (d/D/log)",
      xlab = p$x_label,
      ylab = "Transformed value"
    )
  }, res = 96)



  
  
  
  # ====================================================================
  # ====================================================================
  # ====================================================================
  
  
  
  
  
  # map input$plot_theme to a ggplot theme
  theme_picker <- function(key = "Minimal") {
    switch(key,
           "Minimal" = ggplot2::theme_minimal(),
           "Classic" = ggplot2::theme_classic(),
           "Light"   = ggplot2::theme_light(),
           "Dark"    = ggplot2::theme_dark(),
           "BW"      = ggplot2::theme_bw(),
           "Void"    = ggplot2::theme_void(),
           ggplot2::theme_gray()  # default "Gray"
    )
  }
  
  # apply selected theme + common text sizes
  add_theme <- function(g) {
    g +
      theme_picker(input$plot_theme) +
      ggplot2::theme(
        axis.text  = ggplot2::element_text(size = input$tickSize),
        axis.title = ggplot2::element_text(size = input$tickSize + 2),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )
  }
  
  # ---- UI wrapper so the plot can size dynamically ----
  output$tsPlot_Choice_UI <- renderUI({
    plotOutput(
      "tsPlot_Choice",
      width  = getPlotDim(input$plot_width  %||% "800"),
      height = getPlotDim(input$plot_height %||% "500")
    )
  })
  
  # ---- Color picker UI (already in your sidebar) ----
  output$ts_color_ui <- renderUI({
    tagList(
      tags$input(type = "color", id = "ts_line_color", value = "#2C7FB8"),
      tags$label("Series color", style = "color: #2C7FB8;"),
      br(), br(),
      tags$script(HTML("
      $(document).ready(function() {
        var el = document.getElementById('ts_line_color');
        if (el) Shiny.setInputValue('ts_line_color', el.value);
        $(document).on('input', '#ts_line_color', function() {
          Shiny.setInputValue(this.id, this.value);
        });
      });
    "))
    )
  })
  
  # ---- Plot router for the Plot (*) tab (uses theme) ----
  output$tsPlot_Choice <- renderPlot({
    req(myData_Choice())
    req(input$plot_type_choice)
    
    ts_obj <- myData_Choice()
    p      <- prepared()  # for x-axis label
    freq   <- tryCatch(stats::frequency(ts_obj), error = function(e) NA_real_)
    
    df_ts <- function(z) {
      if (inherits(z, "ts")) {
        data.frame(t = as.numeric(stats::time(z)), y = as.numeric(z))
      } else {
        data.frame(t = seq_along(z), y = as.numeric(z))
      }
    }
    
    k_ma  <- max(2L, as.integer(input$ma_k %||% 5))
    lag_m <- max(1L, as.integer(input$lag_m %||% 12))
    
    plt <- switch(
      input$plot_type_choice,
      
      "Line" = {
        forecast::autoplot(
          ts_obj,
          size   = 1,
          colour = input$ts_line_color %||% "#2C7FB8"
        ) +
          ggplot2::labs(title = "Transformed series", x = p$x_label, y = "Transformed value")
      },
      
      "Points" = {
        d <- df_ts(ts_obj)
        ggplot2::ggplot(d, ggplot2::aes(t, y)) +
          ggplot2::geom_point(size = 1) +
          ggplot2::labs(title = "Points", x = p$x_label, y = "Transformed value")
      },
      
      "Line + Points" = {
        d <- df_ts(ts_obj)
        ggplot2::ggplot(d, ggplot2::aes(t, y)) +
          ggplot2::geom_line() +
          ggplot2::geom_point(size = 0.9, alpha = 0.8) +
          ggplot2::labs(title = "Line + Points", x = p$x_label, y = "Transformed value")
      },
      
      "Smoothed (LOESS)" = {
        d <- df_ts(ts_obj)
        ggplot2::ggplot(d, ggplot2::aes(t, y)) +
          ggplot2::geom_line(alpha = 0.4) +
          ggplot2::geom_smooth(method = "loess", se = FALSE, span = 0.2) +
          ggplot2::labs(title = "LOESS smooth", x = p$x_label, y = "Transformed value")
      },
      
      "Moving average" = {
        d <- df_ts(ts_obj)
        ma <- stats::filter(d$y, rep(1/k_ma, k_ma), sides = 2)
        d$ma <- as.numeric(ma)
        ggplot2::ggplot(d, ggplot2::aes(t, y)) +
          ggplot2::geom_line(alpha = 0.4) +
          ggplot2::geom_line(ggplot2::aes(y = ma), size = 1) +
          ggplot2::labs(title = sprintf("Moving average (k = %d)", k_ma), x = p$x_label, y = "Transformed value")
      },
      
      "Cumulative sum" = {
        d <- df_ts(ts_obj); d$cum <- cumsum(d$y)
        ggplot2::ggplot(d, ggplot2::aes(t, cum)) +
          ggplot2::geom_line() +
          ggplot2::labs(title = "Cumulative sum", x = p$x_label, y = "Cumulative value")
      },
      
      "Seasonal plot" = {
        validate(need(is.finite(freq) && freq > 1, "Seasonal plot requires frequency > 1."))
        forecast::ggseasonplot(ts_obj, year.labels = TRUE) +
          ggplot2::labs(title = "Seasonal plot", x = p$x_label, y = "Value")
      },
      
      "Seasonal subseries" = {
        validate(need(is.finite(freq) && freq > 1, "Seasonal subseries requires frequency > 1."))
        forecast::ggsubseriesplot(ts_obj) +
          ggplot2::labs(title = "Seasonal subseries", x = p$x_label, y = "Value")
      },
      
      "Polar seasonal" = {
        validate(need(is.finite(freq) && freq > 1, "Polar seasonal requires frequency > 1."))
        forecast::ggseasonplot(ts_obj, polar = TRUE) +
          ggplot2::labs(title = "Polar seasonal plot", x = p$x_label, y = "Value")
      },
      
      "Seasonal boxplot" = {
        validate(need(is.finite(freq) && freq > 1, "Seasonal boxplot requires frequency > 1."))
        d <- if (inherits(ts_obj, "ts")) data.frame(season = stats::cycle(ts_obj), y = as.numeric(ts_obj))
        else data.frame(season = factor(1), y = as.numeric(ts_obj))
        ggplot2::ggplot(d, ggplot2::aes(x = factor(season), y = y)) +
          ggplot2::geom_boxplot() +
          ggplot2::labs(title = "Seasonal boxplot", x = "Season", y = "Value")
      },
      
      "Classical decomposition (additive)" = {
        validate(need(is.finite(freq) && freq > 1, "Classical decomposition requires frequency > 1."))
        dc <- stats::decompose(ts_obj, type = "additive")
        forecast::autoplot(dc) + ggplot2::labs(title = "Classical decomposition (additive)")
      },
      
      "Classical decomposition (multiplicative)" = {
        validate(need(is.finite(freq) && freq > 1, "Classical decomposition requires frequency > 1."))
        dc <- stats::decompose(ts_obj, type = "multiplicative")
        forecast::autoplot(dc) + ggplot2::labs(title = "Classical decomposition (multiplicative)")
      },
      
      "STL decomposition" = {
        validate(need(is.finite(freq) && freq > 1, "STL decomposition requires frequency > 1."))
        decomp <- stats::stl(ts_obj, s.window = "periodic")
        forecast::autoplot(decomp) + ggplot2::labs(title = "STL decomposition")
      },
      
      "Histogram" = {
        xx <- as.numeric(stats::na.omit(ts_obj))
        ggplot2::ggplot(data.frame(x = xx), ggplot2::aes(x)) +
          ggplot2::geom_histogram(bins = 30) +
          ggplot2::labs(title = "Histogram", x = "Value", y = "Count")
      },
      
      "Density" = {
        xx <- as.numeric(stats::na.omit(ts_obj))
        ggplot2::ggplot(data.frame(x = xx), ggplot2::aes(x)) +
          ggplot2::geom_density() +
          ggplot2::labs(title = "Density", x = "Value", y = "Density")
      },
      
      "QQ plot" = {
        xx <- as.numeric(stats::na.omit(ts_obj))
        ggplot2::ggplot(data.frame(x = xx), ggplot2::aes(sample = x)) +
          ggplot2::stat_qq() +
          ggplot2::stat_qq_line() +
          ggplot2::labs(title = "Normal Q-Q plot", x = "Theoretical quantiles", y = "Sample quantiles")
      },
      
      "Lag-1 scatter" = {
        xx <- as.numeric(stats::na.omit(ts_obj))
        validate(need(length(xx) >= 2, "Not enough data for lag-1 scatter."))
        d <- data.frame(x = xx[-length(xx)], y = xx[-1])
        ggplot2::ggplot(d, ggplot2::aes(x, y)) +
          ggplot2::geom_point() +
          ggplot2::geom_smooth(method = "lm", se = FALSE) +
          ggplot2::labs(title = "Lag-1 scatter", x = "y(t-1)", y = "y(t)")
      },
      
      "Lag plot (1..m)" = {
        xx <- as.numeric(stats::na.omit(ts_obj))
        validate(need(length(xx) > (lag_m + 1), "Increase data or reduce m for lag plots."))
        forecast::gglagplot(ts_obj, lags = lag_m) +
          ggplot2::labs(title = sprintf("Lag plot (1..%d)", lag_m))
      },
      
      "ACF" = {
        forecast::ggAcf(ts_obj) + ggplot2::labs(title = "ACF")
      },
      
      "PACF" = {
        forecast::ggPacf(ts_obj) + ggplot2::labs(title = "PACF")
      },
      

      
      
      # inside your switch(input$plot_type_choice, ...)
      "ACF+PACF" = {
        # ACF (top) & PACF (bottom)
        p_acf  <- forecast::ggAcf(ts_obj)  + ggplot2::labs(title = "ACF")
        p_pacf <- forecast::ggPacf(ts_obj) + ggplot2::labs(title = "PACF")
        
        # Apply your theme helper
        p_acf  <- add_theme(p_acf)
        p_pacf <- add_theme(p_pacf)
        
        # Vertical layout: ACF on top, PACF below
        gridExtra::grid.arrange(p_acf, p_pacf, ncol = 1, heights = c(1, 1))
      },
      
      
      
      # # inside your switch(input$plot_type_choice, ...)
      # "ACF+PACF" = {
      #   # ACF / PACF only
      #   p_acf  <- forecast::ggAcf(ts_obj)  + ggplot2::labs(title = "ACF")
      #   p_pacf <- forecast::ggPacf(ts_obj) + ggplot2::labs(title = "PACF")
      #   
      #   # Apply theme
      #   p_acf  <- add_theme(p_acf)
      #   p_pacf <- add_theme(p_pacf)
      #   
      #   # Side-by-side layout (1 row, 2 columns)
      #   gridExtra::grid.arrange(p_acf, p_pacf, ncol = 2)
      # },
      
      
       
      # "ACF+PACF" = {
      #   # Time plot
      #   p_time <- forecast::autoplot(
      #     ts_obj,
      #     size   = 1,
      #     colour = input$ts_line_color %||% "#2C7FB8"
      #   ) +
      #     ggplot2::labs(
      #       title = "Time plot",
      #       x = p$x_label,
      #       y = "Transformed value"
      #     )
      #   
      #   # ACF / PACF
      #   p_acf  <- forecast::ggAcf(ts_obj)  + ggplot2::labs(title = "ACF")
      #   p_pacf <- forecast::ggPacf(ts_obj) + ggplot2::labs(title = "PACF")
      #   
      #   # Apply your theme to each ggplot
      #   p_time <- add_theme(p_time)
      #   p_acf  <- add_theme(p_acf)
      #   p_pacf <- add_theme(p_pacf)
      #   
      #   # Stack them (1 column)
      #   gridExtra::grid.arrange(p_time, p_acf, p_pacf, ncol = 1,
      #                           heights = c(2, 1, 1))
      # },
      
      
      
      # at top of server.R (if not already)
      # library(gridExtra)
      
      "Time + ACF+PACF" = {
        # Time plot (top)
        p_time <- forecast::autoplot(
          ts_obj,
          size   = 1,
          colour = input$ts_line_color %||% "#2C7FB8"
        ) +
          ggplot2::labs(
            title = "Time plot",
            x = p$x_label,
            y = "Transformed value"
          )
        
        # ACF (bottom-left) & PACF (bottom-right)
        p_acf  <- forecast::ggAcf(ts_obj)  + ggplot2::labs(title = "ACF")
        p_pacf <- forecast::ggPacf(ts_obj) + ggplot2::labs(title = "PACF")
        
        # Apply your selected theme to each plot
        p_time <- add_theme(p_time)
        p_acf  <- add_theme(p_acf)
        p_pacf <- add_theme(p_pacf)
        
        # Bottom row: ACF (left) | PACF (right)
        bottom_row <- gridExtra::arrangeGrob(p_acf, p_pacf, ncol = 2)
        
        # Final layout: Time on top; ACF+PACF in one row at bottom
        gridExtra::grid.arrange(p_time, bottom_row, ncol = 1, heights = c(1.3, 1))
      },
      
      
      
      # "Time + ACF+PACF" = {
      #   forecast::ggtsdisplay(
      #     ts_obj,
      #     main = "Diagnostics (Time, ACF, PACF)",
      #     xlab = p$x_label, ylab = "Transformed value"
      #   )
      # },
      
      "Periodogram" = {
        xx <- as.numeric(stats::na.omit(ts_obj))
        validate(need(length(xx) > 5, "More data needed for periodogram."))
        sp <- stats::spec.pgram(xx, detrend = TRUE, taper = 0.1, plot = FALSE)
        d  <- data.frame(freq = sp$freq, spec = sp$spec)
        ggplot2::ggplot(d, ggplot2::aes(freq, spec)) +
          ggplot2::geom_line() +
          ggplot2::labs(title = "Periodogram", x = "Frequency", y = "Spectral density")
      }
    )
    
    # apply theme (for ggplot outputs)
    if (inherits(plt, "ggplot")) {
      plt <- add_theme(plt)
    }
    plt
  }, res = 96)
  
  
  
  
  
  # # --- UI wrapper so the plot can size dynamically (same style as your other UI wrappers)
  # output$tsPlot_Choice_UI <- renderUI({
  #   plotOutput(
  #     "tsPlot_Choice",
  #     width  = getPlotDim(input$plot_width  %||% 800),
  #     height = getPlotDim(input$plot_height %||% 500)
  #   )
  # })
  # 
  # # --- Color picker UI (used in this tab’s sidebar)
  # output$ts_color_ui <- renderUI({
  #   tagList(
  #     tags$input(type = "color", id = "ts_line_color", value = "#2C7FB8"),
  #     tags$label("Series color", style = "color: #2C7FB8;"),
  #     br(), br(),
  #     tags$script(HTML("
  #     $(document).ready(function() {
  #       var el = document.getElementById('ts_line_color');
  #       if (el) Shiny.setInputValue('ts_line_color', el.value);
  #       $(document).on('input', '#ts_line_color', function() {
  #         Shiny.setInputValue(this.id, this.value);
  #       });
  #     });
  #   "))
  #   )
  # })
  # 
  # # --- Main time plot of the transformed series (uses myData_Choice())
  # output$tsPlot_Choice <- renderPlot({
  #   req(myData_Choice())
  #   p <- prepared()  # for axis label (date vs index)
  #   
  #   forecast::autoplot(
  #     myData_Choice(),
  #     size   = 1,
  #     colour = input$ts_line_color %||% "#2C7FB8"
  #   ) +
  #     ggplot2::labs(
  #       title = "Transformed series",
  #       x = p$x_label,
  #       y = "Transformed value"
  #     ) +
  #     ggplot2::theme_minimal() +
  #     ggplot2::theme(
  #       axis.text  = ggplot2::element_text(size = input$tickSize),
  #       axis.title = ggplot2::element_text(size = input$tickSize + 2)
  #     )
  # }, res = 96)
  
  
  
  
  
  
  
  
  # ====================================================================
  # ====================================================================
  # ====================================================================
  
  
  
  
  # ---- UI wrapper (dynamic size) ----
  output$difference2ACFPACF_UI <- renderUI({
    plotOutput(
      "difference2ACFPACF",
      width  = getPlotDim(input$plot_width  %||% 800),
      height = getPlotDim(input$plot_height %||% 600)
    )
  })
  
  # ---- Combined ACF + PACF (stacked) ----
  output$difference2ACFPACF <- renderPlot({
    req(myData_Choice())
    
    p1 <- forecast::ggAcf(myData_Choice()) +
      ggplot2::labs(title = "ACF") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text  = ggplot2::element_text(size = input$tickSize),
        axis.title = ggplot2::element_text(size = input$tickSize + 2),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )
    
    p2 <- forecast::ggPacf(myData_Choice()) +
      ggplot2::labs(title = "PACF") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text  = ggplot2::element_text(size = input$tickSize),
        axis.title = ggplot2::element_text(size = input$tickSize + 2),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )
    
    gridExtra::grid.arrange(p1, p2, ncol = 1, top = "ACF & PACF of transformed series")
  }, res = 96)
  
  
  
  
  
  # ====================================================================
  # ====================================================================
  # ====================================================================
  
  
  output$teststationarited3St <- renderPrint({
    # ---------- Helpers ----------
    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
    to_num_safe <- function(v, default = NA_real_) {
      out <- suppressWarnings(as.numeric(v))
      if (length(out) == 0 || all(is.na(out)) || !is.finite(out[1])) default else out[1]
    }
    to_int_safe <- function(v, default = 0L) {
      out <- suppressWarnings(as.integer(v))
      if (length(out) == 0 || is.na(out[1]) || !is.finite(out[1])) default else out[1]
    }
    safe_head_tail <- function(x, n = 5) {
      x <- as.numeric(x); x <- x[is.finite(x)]
      if (length(x) == 0) return(list(head = numeric(0), tail = numeric(0)))
      list(head = head(x, n), tail = tail(x, n))
    }
    tau_row_for <- function(type_in) switch(type_in, "none"="tau1", "drift"="tau2", "trend"="tau3", "tau3")
    cval_pick_safe <- function(cval_obj, row, col) {
      if (is.null(cval_obj)) return(NA_real_)
      ans <- NA_real_
      if (is.matrix(cval_obj)) {
        if (!missing(row) && !missing(col) && row %in% rownames(cval_obj) && col %in% colnames(cval_obj)) {
          ans <- suppressWarnings(as.numeric(cval_obj[row, col]))
        } else if (!missing(col) && col %in% colnames(cval_obj)) {
          ans <- suppressWarnings(as.numeric(cval_obj[1, col]))
        } else {
          ans <- suppressWarnings(as.numeric(cval_obj[1, 1]))
        }
      } else {
        nm <- names(cval_obj)
        if (!is.null(nm) && col %in% nm) ans <- suppressWarnings(as.numeric(cval_obj[[col]]))
        if (!is.finite(ans) && length(cval_obj) >= 3) {
          if (identical(col, "10pct")) ans <- suppressWarnings(as.numeric(cval_obj[1]))
          if (identical(col, "5pct"))  ans <- suppressWarnings(as.numeric(cval_obj[2]))
          if (identical(col, "1pct"))  ans <- suppressWarnings(as.numeric(cval_obj[3]))
        }
        if (!is.finite(ans)) ans <- suppressWarnings(as.numeric(cval_obj[1]))
      }
      ans
    }
    fmt_p <- function(p) {
      if (!is.finite(p)) return("NA")
      if (p < .001) "<0.001" else sprintf("%.6f", p)
    }
    fmt_num <- function(z, d=6) ifelse(is.finite(z), sprintf(paste0("%.", d ,"f"), z), "NA")
    tick <- function(ok) if (isTRUE(ok)) "[✓]" else "[X]"
    warn <- function(ok) if (isTRUE(ok)) "[✓]" else "[!]"
    qmark <- function(ok) if (isTRUE(ok)) "[✓]" else "[?]"
    
    # ---------- Inputs (from your UI) ----------
    alt_in  <- input$alternd2St %||% input$alternSt       # "stationary", "explosive", "regression"
    lag_in  <- input$LagOrderADFd2St %||% input$LagOrderADFSt
    a_in    <- input$alphaSt2
    type_in <- input$adfTypeSt2                            # "none", "drift", "trend"
    
    # Transformation flags (provenance)
    d_in   <- input$d_n  %||% NA
    D_in   <- input$DS_n %||% NA
    log_in <- input$check_box %||% FALSE
    
    # ---------- Data ----------
    req(myData_Choice())
    x_raw <- myData_Choice()
    na_before <- sum(is.na(x_raw))
    x_class <- paste(class(x_raw), collapse = ", ")
    x_freq  <- if (inherits(x_raw, "ts")) tryCatch(stats::frequency(x_raw), error = function(e) NA_integer_) else NA_integer_
    x <- as.numeric(stats::na.omit(x_raw))
    N <- length(x)
    
    # ---------- Hyper-parameters ----------
    k <- to_int_safe(lag_in, default = 0L); if (!is.finite(k) || k < 0) k <- 0L
    alpha_raw <- as.character(a_in)
    alpha_val <- to_num_safe(alpha_raw, default = 0.05)
    alpha_col <- switch(alpha_raw, "0.01"="1pct", "0.05"="5pct", "0.1"="10pct", "0.10"="10pct",
                        "1pct"="1pct", "5pct"="5pct", "10pct"="10pct", "5pct")
    tau_row <- tau_row_for(type_in)
    
    # ---------- Sanity checks ----------
    if (N < 5) {
      cat("==========================================================================\n")
      cat("                         ADF UNIT ROOT DIAGNOSTIC                         \n")
      cat("==========================================================================\n")
      cat("CRITICAL ERROR: Too few observations (N < 5). Provide more data.\n")
      return(invisible(NULL))
    }
    if (!is.finite(stats::sd(x)) || stats::sd(x) == 0) {
      cat("==========================================================================\n")
      cat("                         ADF UNIT ROOT DIAGNOSTIC                         \n")
      cat("==========================================================================\n")
      cat("CRITICAL ERROR: Series is constant/invalid (sd = 0 or NA). Check transforms.\n")
      return(invisible(NULL))
    }
    
    # ---------- Package check ----------
    if (!requireNamespace("urca", quietly = TRUE) || !requireNamespace("tseries", quietly = TRUE)) {
      cat("==========================================================================\n")
      cat("                         ADF UNIT ROOT DIAGNOSTIC                         \n")
      cat("==========================================================================\n")
      cat("ERROR: Please install.packages(c('urca','tseries')) to run ADF/KPSS.\n")
      return(invisible(NULL))
    }
    
    # ---------- Common quantities ----------
    m_mean <- mean(x, na.rm = TRUE)
    m_sd   <- stats::sd(x, na.rm = TRUE)
    lb_lag <- max(1L, min(10L, floor(N / 5)))
    
    # ---------- Core ADF (urca) for chosen type & k ----------
    ur  <- tryCatch(urca::ur.df(x, type = type_in, lags = k), error = function(e) NULL)
    tau_obs  <- if (!is.null(ur)) suppressWarnings(as.numeric(ur@teststat[tau_row])) else NA_real_
    if (!is.finite(tau_obs) && !is.null(ur)) tau_obs <- suppressWarnings(as.numeric(ur@teststat[1]))
    tau_crit <- if (!is.null(ur)) cval_pick_safe(ur@cval, tau_row, alpha_col) else NA_real_
    adf_ts <- if (N > (k + 10)) tryCatch(tseries::adf.test(x, alternative = alt_in, k = k),
                                         error = function(e) NULL) else NULL
    adf_p  <- if (!is.null(adf_ts)) to_num_safe(adf_ts$p.value) else NA_real_
    lb_main <- if (!is.null(ur)) stats::Box.test(ur@res, lag = lb_lag, type = "Ljung-Box") else NULL
    lb_stat<- if (!is.null(lb_main)) to_num_safe(lb_main$statistic) else NA_real_
    lb_p   <- if (!is.null(lb_main)) to_num_safe(lb_main$p.value) else NA_real_
    
    adf_ok_vals <- is.finite(tau_obs) && is.finite(tau_crit)
    adf_stationary <- adf_ok_vals && (tau_obs < tau_crit)
    
    # ---------- KPSS (full sample) ----------
    kpss_type <- if (type_in == "trend") "Trend" else "Level"
    kpss_ts   <- tryCatch(tseries::kpss.test(x, null = kpss_type), error = function(e) NULL)
    kpss_p    <- if (!is.null(kpss_ts)) to_num_safe(kpss_ts$p.value) else NA_real_
    kpss_uc   <- tryCatch(urca::ur.kpss(x, type = if (type_in == "trend") "tau" else "mu"),
                          error = function(e) NULL)
    eta_obs_uc  <- if (!is.null(kpss_uc)) to_num_safe(kpss_uc@teststat) else NA_real_
    eta_col     <- if (alpha_val <= 0.01) "1pct" else if (alpha_val <= 0.05) "5pct" else "10pct"
    eta_crit_uc <- if (!is.null(kpss_uc)) cval_pick_safe(kpss_uc@cval, 1, eta_col) else NA_real_
    kpss_by_p   <- is.finite(kpss_p)     && (kpss_p < alpha_val)
    kpss_by_eta <- is.finite(eta_obs_uc) && is.finite(eta_crit_uc) && (eta_obs_uc > eta_crit_uc)
    kpss_stationary <- !(kpss_by_p || kpss_by_eta)
    
    # ---------- Structural break (Pettitt) ----------
    pett_U <- pett_p <- NA_real_; pett_cp <- NA_integer_
    if (requireNamespace("trend", quietly = TRUE)) {
      pet <- tryCatch(trend::pettitt.test(x), error = function(e) NULL)
      if (!is.null(pet)) {
        pett_U <- to_num_safe(pet$statistic)
        pett_p <- to_num_safe(pet$p.value)
        # trend::pettitt.test typically returns estimate (break index); guard for NA
        pett_cp <- to_int_safe(pet$estimate[1], default = NA_integer_)
        if (!is.finite(pett_cp)) {
          # conservative fallback: midpoint (only if needed)
          pett_cp <- as.integer(floor(N/2))
        }
        if (pett_cp < 2L || pett_cp > (N - 2L)) pett_cp <- NA_integer_  # unusable splits
      }
    }
    
    # ---------- KPSS by segments (if break usable) ----------
    seg1_p <- seg2_p <- NA_real_
    if (is.finite(pett_cp)) {
      if (pett_cp >= 8 && (N - pett_cp) >= 8) {
        seg1 <- tryCatch(tseries::kpss.test(x[seq_len(pett_cp)], null = kpss_type), error = function(e) NULL)
        seg2 <- tryCatch(tseries::kpss.test(x[(pett_cp + 1L):N], null = kpss_type), error = function(e) NULL)
        if (!is.null(seg1)) seg1_p <- to_num_safe(seg1$p.value)
        if (!is.null(seg2)) seg2_p <- to_num_safe(seg2$p.value)
      }
    }
    
    # ---------- Phase 0: k* scan & spec sensitivity ----------
    # k* = smallest k such that Ljung–Box p-value on ADF residuals > alpha, for the user-chosen type_in
    k_grid <- 0L:min(12L, max(1L, floor(N/10)))
    lbp_by_k <- rep(NA_real_, length(k_grid))
    for (i in seq_along(k_grid)) {
      kk <- k_grid[i]
      ur_i <- tryCatch(urca::ur.df(x, type = type_in, lags = kk), error = function(e) NULL)
      lb_i <- if (!is.null(ur_i)) tryCatch(stats::Box.test(ur_i@res, lag = lb_lag, type = "Ljung-Box"), error = function(e) NULL) else NULL
      lbp_by_k[i] <- if (!is.null(lb_i)) to_num_safe(lb_i$p.value) else NA_real_
    }
    k_suggest <- NA_integer_
    idx_ok <- which(is.finite(lbp_by_k) & lbp_by_k > alpha_val)
    if (length(idx_ok) > 0) k_suggest <- k_grid[min(idx_ok)] else k_suggest <- k
    
    # Spec sensitivity across types for k = selected and k = k_suggest
    types <- c("none","drift","trend")
    spec_line <- function(kk) {
      for (tp in types) {
        tau_row_tp <- tau_row_for(tp)
        ur_tp <- tryCatch(urca::ur.df(x, type = tp, lags = kk), error = function(e) NULL)
        tau_tp  <- if (!is.null(ur_tp)) suppressWarnings(as.numeric(ur_tp@teststat[tau_row_tp])) else NA_real_
        if (!is.finite(tau_tp) && !is.null(ur_tp)) tau_tp <- suppressWarnings(as.numeric(ur_tp@teststat[1]))
        crit_tp <- if (!is.null(ur_tp)) cval_pick_safe(ur_tp@cval, tau_row_tp, alpha_col) else NA_real_
        lb_tp   <- if (!is.null(ur_tp)) tryCatch(stats::Box.test(ur_tp@res, lag = lb_lag, type = "Ljung-Box"), error = function(e) NULL) else NULL
        lbp_tp  <- if (!is.null(lb_tp)) to_num_safe(lb_tp$p.value) else NA_real_
        adf_dec <- if (is.finite(tau_tp) && is.finite(crit_tp) && (tau_tp < crit_tp)) "ADF=STATIONARY" else "ADF=NON-STATIONARY"
        cat(sprintf("     %s type=%-5s | tau=%s | crit=%s | %s | LB p=%s\n",
                    tick(is.finite(lbp_tp) && lbp_tp > alpha_val),
                    tp,
                    fmt_num(tau_tp, 4),
                    fmt_num(crit_tp, 4),
                    adf_dec,
                    fmt_num(lbp_tp, 4)))
      }
    }
    
    # ---------- Agreement & flags ----------
    agreement <- (adf_stationary && kpss_stationary) || (!adf_stationary && !kpss_stationary)
    
    # ---------- HEADER ----------
    cat("==========================================================================\n")
    cat("                        ADF UNIT ROOT DIAGNOSTIC                          \n")
    cat("==========================================================================\n")
    cat(sprintf(" MODEL TYPE : %-10s | SAMPLE SIZE (N) : %s\n", toupper(type_in), N))
    cat(sprintf(" LAG ORDER  : %-10s | SIGNIFICANCE (α) : %s\n", k, fmt_num(alpha_val, 4)))
    cat(sprintf(" MOMENTS    : Mean: %s | Std.Dev: %s\n",
                fmt_num(m_mean, 4), fmt_num(m_sd, 4)))
    cat("--------------------------------------------------------------------------\n")
    cat(" TRANSFORMATION PROVENANCE (what ALL tests use):\n")
    cat(sprintf(" [ ] Source object class           : %s\n", x_class))
    cat(sprintf(" [ ] Frequency (if ts)             : %s\n", ifelse(is.finite(x_freq), x_freq, "NA")))
    cat(sprintf(" [ ] NA count before na.omit       : %s\n", na_before))
    cat(sprintf(" [ ] Transformation inputs (UI)    : log=%s | d=%s | D=%s\n",
                ifelse(isTRUE(log_in),"ON","OFF"), as.character(d_in), as.character(D_in)))
    ht <- safe_head_tail(x, 5)
    cat(sprintf(" [ ] First 5 values used in tests  : %s\n", paste(round(ht$head, 4), collapse = ", ")))
    cat(sprintf(" [ ] Last  5 values used in tests  : %s\n\n", paste(round(ht$tail, 4), collapse = ", ")))
    
    # ---------- PHASE 0 ----------
    cat("==========================================================================\n")
    cat("PHASE 0: DECISION AID (Lag + Spec Sensitivity)\n")
    cat("==========================================================================\n")
    cat(sprintf(" • Ljung-Box reference lag used in scan : %d  ; (LB lag = min(10, floor(N/5))\n", lb_lag))
    cat(sprintf(" [!] Suggested k* (scan)                : %s  (smallest k with Ljung-Box p-value > alpha (whiter residuals))\n",
                ifelse(is.finite(k_suggest), k_suggest, "NA")))
    if (is.finite(k_suggest) && k_suggest != k) {
      cat(sprintf(" [!] You selected k=%d. Consider trying k=%d and re-running.\n", k, k_suggest))
    } else {
      cat(sprintf(" [✓] You selected k=%d. This already meets the LB>α rule-of-thumb.\n", k))
    }
    cat("--------------------------------------------------------------------------\n")
    cat(" ADF SPEC SENSITIVITY (ur.df):\n")
    cat("   (Prefer: LB ok + stable decision across types)\n")
    cat(sprintf("  • k=%d\n", k)); spec_line(k)
    if (is.finite(k_suggest) && k_suggest != k) {
      cat(sprintf("  • k=%d\n", k_suggest)); spec_line(k_suggest)
    }
    cat("--------------------------------------------------------------------------\n")
    
    # ---------- PHASE 1 ----------
    cat("==========================================================================\n")
    cat("PHASE 1: ADF UNIT ROOT TEST\n")
    cat("==========================================================================\n")
    cat(" • H0: The series has a Unit Root (Non-Stationary).\n")
    cat(" • Ha: The series is Stationary (Mean Reverting).\n")
    cat(sprintf(" -> CRITERIA: Reject H0 if Tau-Obs (%s) < Tau-Crit (%s)\n",
                fmt_num(tau_obs, 4), fmt_num(tau_crit, 4)))
    cat("\n RESULT:\n")
    cat(sprintf("  - Tau Observed : %s \n", fmt_num(tau_obs, 6)))
    cat(sprintf("  - Tau Critical : %s \n", fmt_num(tau_crit, 2)))
    cat(sprintf("  - P-Value (Ref): %s \n", ifelse(is.finite(adf_p), fmt_p(adf_p), "NA")))
    cat("\n DECISION:\n")
    if (adf_ok_vals && adf_stationary) {
      cat("  -> REJECT H0: Evidence suggests the series is STATIONARY.\n")
    } else if (adf_ok_vals && !adf_stationary) {
      cat("  -> FAIL TO REJECT H0: Evidence suggests the series is NON-STATIONARY.\n")
    } else {
      cat("  -> INCONCLUSIVE: Missing statistic/critical value; regression may be ill-conditioned.\n")
    }
    
    # ---------- PHASE 2 ----------
    cat("\n==========================================================================\n")
    cat("PHASE 2: RESIDUAL DIAGNOSTICS (LJUNG-BOX)\n")
    cat("==========================================================================\n")
    cat(" • H0: Residuals are White Noise (No Autocorrelation).\n")
    cat(" • Ha: Residuals are Correlated (Lags are insufficient).\n")
    cat(sprintf(" -> CRITERIA: Reject H0 if P-Value (%s) < α (%s)\n", fmt_num(lb_p, 6), fmt_num(alpha_val, 4)))
    cat("\n RESULT:\n")
    cat(sprintf("  - LB Statistic : %s \n", fmt_num(lb_stat, 6)))
    cat(sprintf("  - LB P-Value   : %s \n", fmt_num(lb_p, 6)))
    cat(sprintf("  - LB Lag used  : %d \n", lb_lag))
    cat("    it’s normal to choose LB lag by a rule-of-thumb (like min(10, floor(N/5)) \n")
    cat("\n DECISION:\n")
    if (is.finite(lb_p)) {
      if (lb_p > alpha_val) cat("  -> FAIL TO REJECT H0: Residuals are White Noise. [ADF more reliable]\n")
      else                  cat("  -> REJECT H0: Residual autocorrelation remains; increase k or difference.\n")
    } else {
      cat("  -> INCONCLUSIVE: Ljung–Box p-value is NA.\n")
    }
    
    # ---------- PHASE 3A ----------
    cat("\n==========================================================================\n")
    cat("PHASE 3: KPSS + STRUCTURAL BREAK (Pettitt) + KPSS SEGMENTS\n")
    cat("==========================================================================\n\n")
    cat("PHASE 3A: KPSS (Stationarity Confirmation)\n")
    cat(sprintf(" • H0: The series is Stationary around a %s.\n", if (kpss_type=="Trend") "Trend" else "Level"))
    cat(" • Ha: The series is Non-Stationary.\n")
    cat(sprintf(" • CRITERIA (p-value) : Reject H0 if p-value (%s) < α (%s)\n",
                ifelse(is.finite(kpss_p), fmt_num(kpss_p, 6), "NA"), fmt_num(alpha_val,4)))
    cat(sprintf(" • CRITERIA (eta)     : Reject H0 if Eta-Obs (%s) > Eta-Crit (%s)  [urca]\n",
                fmt_num(eta_obs_uc, 6), fmt_num(eta_crit_uc, 6)))
    cat("\n RESULT:\n")
    cat(sprintf("  - Eta (Observed value) [urca]   : %s \n", fmt_num(eta_obs_uc, 5)))
    cat(sprintf("  - Eta (Critical value) [urca]   : %s \n", fmt_num(eta_crit_uc, 3)))
    cat(sprintf("  - p-value (one-tailed) [tseries]: %s \n", ifelse(is.finite(kpss_p), fmt_num(kpss_p, 6), "NA")))
    cat(sprintf("  - Eta (Observed) [tseries, FYI] : %s \n", fmt_num(eta_obs_uc, 5)))
    cat("\n DECISION:\n")
    if (kpss_stationary) cat("  -> FAIL TO REJECT H0: Stationarity supported by KPSS.\n")
    else                 cat("  -> REJECT H0: Non-stationarity indicated by KPSS.\n")
    
    # ---------- PHASE 3B ----------
    cat("--------------------------------------------------------------------------\n")
    cat("PHASE 3B: STRUCTURAL BREAK CHECK (Pettitt) + KPSS SEGMENT CHECK\n")
    cat(" • Goal: Detect a single change-point (median shift) that can distort KPSS.\n")
    cat(" • If a break exists, we re-run KPSS before/after the break.\n")
    cat(sprintf(" • Pettitt p-value: %s | α: %s\n\n", ifelse(is.finite(pett_p), fmt_num(pett_p, 6), "NA"), fmt_num(alpha_val, 4)))
    cat(" RESULT:\n")
    cat(sprintf("  - Break index (estimate): %s \n", ifelse(is.finite(pett_cp), pett_cp, "NA")))
    cat(sprintf("  - Pettitt statistic     : %s \n", fmt_num(pett_U, 0)))
    cat(sprintf("  - Pettitt p-value       : %s \n", ifelse(is.finite(pett_p), fmt_num(pett_p, 6), "NA")))
    cat("\n DECISION:\n")
    if (is.finite(pett_p)) {
      if (pett_p < alpha_val) cat("  -> REJECT H0: Break detected. Full-sample KPSS/ADF may be contaminated.\n")
      else                    cat("  -> FAIL TO REJECT H0: No strong single break detected.\n")
    } else {
      cat("  -> INCONCLUSIVE: Pettitt p-value NA.\n")
    }
    
    # ---------- PHASE 3C ----------
    cat("--------------------------------------------------------------------------\n")
    cat("PHASE 3C: KPSS RE-CHECK BY SEGMENTS\n")
    cat("  • Segment 1 = [1 .. break]\n")
    cat("  • Segment 2 = [break+1 .. N]\n")
    cat(sprintf("  - KPSS p-value (Segment 1): %s\n", ifelse(is.finite(seg1_p), fmt_num(seg1_p, 6), "NA")))
    cat(sprintf("  - KPSS p-value (Segment 2): %s\n", ifelse(is.finite(seg2_p), fmt_num(seg2_p, 6), "NA")))
    cat("\n INTERPRETATION:\n")
    if (is.finite(seg1_p) && is.finite(seg2_p)) {
      if (seg1_p >= alpha_val && seg2_p >= alpha_val)
        cat("  [✓] Both segments look stationary by KPSS.\n      -> Full-sample non-stationarity may be break-driven.\n")
      else if (seg1_p < alpha_val && seg2_p < alpha_val)
        cat("  [X] Both segments look non-stationary by KPSS.\n")
      else
        cat("  [?] Mixed evidence: one segment stationary, the other not.\n")
    } else {
      cat("  [?] Segment KPSS not available (short segments or missing package).\n")
    }
    
    # ---------- PHASE 4: Final verdict & advice ----------
    cat("\n==========================================================================\n")
    cat("PHASE 4: FINAL ACADEMIC VERDICT & ADVICE\n")
    cat("==========================================================================\n")
    if (adf_stationary && kpss_stationary) {
      cat(" [✓] VERDICT: CONVERGENT STATIONARITY (ADF & KPSS agree; residuals likely white).\n")
      cat("     ADVICE: Proceed with SARIMA identification with d=0 (choose D via seasonality), then residual checks.\n")
    } else if (!adf_stationary && !kpss_stationary) {
      cat(" [X] VERDICT: CONVERGENT NON-STATIONARITY (ADF & KPSS agree).\n")
      cat("     ADVICE: Difference the series (d=1). If seasonal (m≥2), consider D=1, then re-run tests.\n")
    } else {
      cat(" [?] VERDICT: CONFLICTING RESULTS (ADF vs KPSS).\n")
      cat("     ADVICE: Near-unit-root, trend-vs-drift mismatch, or break contamination are common causes.\n")
      cat("             Use PHASE 0 to pick k/type with LB ok and stable decisions; consider seasonal differencing and break handling.\n")
    }
    
    # ---------- Technical appendix: ADF regression coefficients ----------
    cat("\n TECHNICAL APPENDIX (ADF Regression Coefficients):\n")
    if (!is.null(ur) && !is.null(ur@testreg)) {
      cf <- tryCatch(coef(summary(ur@testreg)), error = function(e) NULL)
      if (!is.null(cf)) {
        printCoefmat(cf, digits = 7, signif.stars = FALSE, P.values = TRUE, has.Pvalue = TRUE)
      } else {
        cat("  (coefficients unavailable)\n")
      }
    } else {
      cat("  (ADF regression unavailable)\n")
    }
    
    # ---------- PHASE 5: Evidence snapshot & checklist ----------
    cat("\n==========================================================================\n")
    cat("PHASE 5: POST-SUMMARY (Academic-quality snapshot)\n")
    cat("==========================================================================\n")
    cat(" EVIDENCE SNAPSHOT (All key outcomes in one place):\n")
    cat(sprintf(" [ ] N (effective sample size)              : %s\n", N))
    cat(sprintf(" [ ] Model type (ADF)                       : %s  (tau row: %s)\n", tolower(type_in), tau_row))
    cat(sprintf(" [ ] Lag order (k)                          : %s\n", k))
    cat(sprintf(" [ ] Alpha (α)                              : %s\n", fmt_num(alpha_val, 4)))
    cat(sprintf(" [ ] Tau-Observed (urca)                    : %s\n", fmt_num(tau_obs, 6)))
    cat(sprintf(" [ ] Tau-Critical (urca, %s)                  : %s\n", alpha_col, fmt_num(tau_crit, 6)))
    cat(sprintf(" [ ] ADF p-value (tseries reference)        : %s\n", ifelse(is.finite(adf_p), fmt_num(adf_p, 6), "NA")))
    cat(sprintf(" [ ] Ljung-Box p-value (residuals)          : %s\n", fmt_num(lb_p, 6)))
    cat(sprintf(" [ ] KPSS Eta observed (urca)               : %s\n", fmt_num(eta_obs_uc, 6)))
    cat(sprintf(" [ ] KPSS Eta critical (urca, %s)             : %s\n", eta_col, fmt_num(eta_crit_uc, 6)))
    cat(sprintf(" [ ] KPSS p-value (one-tailed, tseries)     : %s\n", ifelse(is.finite(kpss_p), fmt_num(kpss_p, 6), "NA")))
    cat(sprintf(" [!] Suggested k* (scan)                     : %s\n", ifelse(is.finite(k_suggest), k_suggest, "NA")))
    cat("--------------------------------------------------------------------------\n")
    cat(" CHECKLIST (Academic-quality acceptance criteria):\n")
    cat(sprintf(" [ ] Tests run on transformed series (x=myData_Choice) : %s YES (same x used everywhere)\n", tick(TRUE)))
    cat(sprintf(" [ ] Variance usable (sd>0)                          : %s SATISFIED\n", tick(TRUE)))
    cat(sprintf(" [ ] Tau is finite and usable                        : %s %s\n", tick(is.finite(tau_obs)), ifelse(is.finite(tau_obs),"SATISFIED","CHECK")))
    cat(sprintf(" [ ] Tau critical value extracted                    : %s %s\n", tick(is.finite(tau_crit)), ifelse(is.finite(tau_crit),"SATISFIED","CHECK")))
    cat(sprintf(" [ ] Residuals pass Ljung-Box (white noise)          : %s %s\n", tick(is.finite(lb_p) && lb_p > alpha_val), ifelse(is.finite(lb_p) && lb_p > alpha_val,"SATISFIED","CHECK")))
    cat(sprintf(" [ ] KPSS Eta observed & critical are usable         : %s %s\n", tick(is.finite(eta_obs_uc) && is.finite(eta_crit_uc)), ifelse(is.finite(eta_obs_uc) && is.finite(eta_crit_uc),"SATISFIED","CHECK")))
    cat(sprintf(" [ ] KPSS p-value is usable                          : %s %s\n", tick(is.finite(kpss_p)), ifelse(is.finite(kpss_p),"SATISFIED","CHECK")))
    cat(sprintf(" [ ] Sample size adequacy                            : %s STRONG\n", tick(N >= 50)))
    cat(sprintf(" [ ] Lag order reasonable relative to N              : %s %s\n", tick(k <= max(12L, floor(N/10))), ifelse(k <= max(12L, floor(N/10)),"OK","LARGE")))
    cat(sprintf(" [ ] Seasonal differencing indicated (UI D>0)         : %s %s\n", warn(isTRUE(to_int_safe(D_in,0L) > 0)), ifelse(isTRUE(to_int_safe(D_in,0L) > 0),"YES","NO / UNKNOWN")))
    cat(sprintf(" [ ] ADF alternative mode (Stationary/Explosive)      : %s\n", paste0(ifelse(is.null(alt_in),"NA",alt_in))))
    cat(sprintf(" [ ] ADF & KPSS agreement                             : %s %s\n",
                qmark(agreement), ifelse(agreement,"AGREEMENT","CONFLICT")))
    cat("     [?] NOTE: conflicts are common with near-unit-root series, trend vs drift mismatch,\n")
    cat("         structural breaks (Pettitt), or missing seasonal differencing.\n")
    cat("--------------------------------------------------------------------------\n")
    cat(" STRUCTURAL BREAK SUMMARY (Pettitt):\n")
    cat(sprintf(" [ ] Pettitt p-value : %s  (Reject H0 if < α)\n", ifelse(is.finite(pett_p), fmt_num(pett_p, 6), "NA")))
    cat("--------------------------------------------------------------------------\n")
    cat("==========================================================================\n")
    cat(" ACTIONABLE NEXT STEPS (What to do now):\n")
    cat("==========================================================================\n")
    cat(sprintf(" %s [1] RESIDUAL AUTOCORRELATION CHECK\n", tick(is.finite(lb_p) && lb_p > alpha_val)))
    cat("     • Ljung-Box is acceptable → ADF regression is less likely biased.\n")
    cat(" [?] [2] RESOLVE ADF vs KPSS CONFLICT\n")
    cat("     • Use PHASE 0 'ADF SPEC SENSITIVITY' to pick the type with LB ok and stable decision.\n")
    cat("     • Try ADF model type variants: none / drift / trend (match KPSS Level vs Trend).\n")
    cat("     • If Pettitt indicates a break: split sample and re-test.\n")
    cat("     • If series is seasonal: test after seasonal differencing (D=1) + maybe log.\n")
    cat("     • Consider variance stabilization: log or Box-Cox (if positive data).\n")
    cat(" [!] [3] SEASONALITY SANITY (especially for AirPassengers-like series)\n")
    if (is.finite(x_freq) && x_freq >= 2) {
      cat(sprintf("     • Detected frequency=%d → seasonality is plausible.\n", x_freq))
    } else {
      cat("     • Frequency unknown → inspect ACF/PACF for seasonal spikes.\n")
    }
    cat(" [✓] [4] EXPLOSIVE MODE NOTE\n")
    cat(sprintf("     • %s\n", ifelse(identical(alt_in,"explosive"),
                                      "Explosive alternative was requested; interpret ADF accordingly.",
                                      "Not in explosive mode → standard stationarity workflow applies.")))
    cat("\n PRACTICAL MODELING PATH (for your Shiny workflow):\n")
    if (!adf_stationary || !kpss_stationary) {
      cat(" [X] Apply differencing (d and/or D) → re-run ADF/KPSS → then identify ARMA.\n")
    } else {
      cat(" [✓] Keep d=0 (and decide D via seasonality) → SARIMA identification and residual diagnostics.\n")
    }
    cat("--------------------------------------------------------------------------\n")
  })
  
  
  
  # output$teststationarited3St <- renderPrint({
  #   # ---------- Helpers ----------
  #   `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
  #   to_num_safe <- function(v, default = NA_real_) {
  #     out <- suppressWarnings(as.numeric(v))
  #     if (length(out) == 0 || all(is.na(out)) || !is.finite(out[1])) default else out[1]
  #   }
  #   to_int_safe <- function(v, default = 0L) {
  #     out <- suppressWarnings(as.integer(v))
  #     if (length(out) == 0 || is.na(out[1]) || !is.finite(out[1])) default else out[1]
  #   }
  #   safe_head_tail <- function(x, n = 5) {
  #     x <- as.numeric(x); x <- x[is.finite(x)]
  #     if (length(x) == 0) return(list(head = numeric(0), tail = numeric(0)))
  #     list(head = head(x, n), tail = tail(x, n))
  #   }
  #   tau_row_for <- function(type_in) switch(type_in, "none"="tau1", "drift"="tau2", "trend"="tau3", "tau3")
  #   cval_pick_safe <- function(cval_obj, row, col) {
  #     if (is.null(cval_obj)) return(NA_real_)
  #     ans <- NA_real_
  #     if (is.matrix(cval_obj)) {
  #       if (!missing(row) && !missing(col) && row %in% rownames(cval_obj) && col %in% colnames(cval_obj)) {
  #         ans <- suppressWarnings(as.numeric(cval_obj[row, col]))
  #       } else if (!missing(col) && col %in% colnames(cval_obj)) {
  #         ans <- suppressWarnings(as.numeric(cval_obj[1, col]))
  #       } else {
  #         ans <- suppressWarnings(as.numeric(cval_obj[1, 1]))
  #       }
  #     } else {
  #       nm <- names(cval_obj)
  #       if (!is.null(nm) && col %in% nm) ans <- suppressWarnings(as.numeric(cval_obj[[col]]))
  #       if (!is.finite(ans) && length(cval_obj) >= 3) {
  #         if (identical(col, "10pct")) ans <- suppressWarnings(as.numeric(cval_obj[1]))
  #         if (identical(col, "5pct"))  ans <- suppressWarnings(as.numeric(cval_obj[2]))
  #         if (identical(col, "1pct"))  ans <- suppressWarnings(as.numeric(cval_obj[3]))
  #       }
  #       if (!is.finite(ans)) ans <- suppressWarnings(as.numeric(cval_obj[1]))
  #     }
  #     ans
  #   }
  #   fmt_p   <- function(p) { if (!is.finite(p)) "NA" else if (p < .001) "<0.001" else sprintf("%.6f", p) }
  #   fmt_num <- function(z, d=6) ifelse(is.finite(z), sprintf(paste0("%.", d ,"f"), z), "NA")
  #   tick    <- function(ok) if (isTRUE(ok)) "[✓]" else "[X]"
  #   warn    <- function(ok) if (isTRUE(ok)) "[✓]" else "[!]"
  #   qmark   <- function(ok) if (isTRUE(ok)) "[✓]" else "[?]"
  #   
  #   # ---------- Inputs (from your UI) ----------
  #   alt_in  <- input$alternd2St %||% input$alternSt   # "stationary", "explosive", "regression"
  #   lag_in  <- input$LagOrderADFd2St %||% input$LagOrderADFSt
  #   a_in    <- input$alphaSt2
  #   type_in <- input$adfTypeSt2                        # "none", "drift", "trend"
  #   
  #   # Transformation flags (provenance only; we assume your upstream pipeline already applies them)
  #   d_in   <- input$d_n  %||% NA
  #   D_in   <- input$DS_n %||% NA
  #   log_in <- input$check_box %||% FALSE
  #   
  #   # ---------- DATA = TRAINING-ONLY if available ----------
  #   # Try to pull the training split; fallback to full series if helper not present
  #   x_source <- NULL
  #   if (exists("ts_train_test", mode = "function")) {
  #     s <- tryCatch(ts_train_test(), error = function(e) NULL)
  #     if (!is.null(s) && !is.null(s$ts_train)) x_source <- s$ts_train
  #   }
  #   if (is.null(x_source)) {
  #     req(myData_Choice())
  #     x_source <- myData_Choice()
  #   }
  #   
  #   x_class <- paste(class(x_source), collapse = ", ")
  #   x_freq  <- if (inherits(x_source, "ts")) tryCatch(stats::frequency(x_source), error = function(e) NA_integer_) else NA_integer_
  #   na_before <- sum(is.na(x_source))
  #   
  #   # The actual sample that ALL tests use (training window, post any upstream transforms)
  #   x_vec <- as.numeric(stats::na.omit(x_source))
  #   N     <- length(x_vec)
  #   
  #   # ---------- Hyper-parameters ----------
  #   k <- to_int_safe(lag_in, default = 0L); if (!is.finite(k) || k < 0) k <- 0L
  #   alpha_raw <- as.character(a_in)
  #   alpha_val <- to_num_safe(alpha_raw, default = 0.05)
  #   alpha_col <- switch(alpha_raw, "0.01"="1pct", "0.05"="5pct", "0.1"="10pct", "0.10"="10pct",
  #                       "1pct"="1pct", "5pct"="5pct", "10pct"="10pct", "5pct")
  #   tau_row <- tau_row_for(type_in)
  #   
  #   # ---------- Sanity checks ----------
  #   if (N < 5) {
  #     cat("==========================================================================\n")
  #     cat("                         ADF UNIT ROOT DIAGNOSTIC                         \n")
  #     cat("==========================================================================\n")
  #     cat("CRITICAL ERROR: Too few observations in TRAINING window (N < 5). Provide more data or reduce differencing.\n")
  #     return(invisible(NULL))
  #   }
  #   if (!is.finite(stats::sd(x_vec)) || stats::sd(x_vec) == 0) {
  #     cat("==========================================================================\n")
  #     cat("                         ADF UNIT ROOT DIAGNOSTIC                         \n")
  #     cat("==========================================================================\n")
  #     cat("CRITICAL ERROR: Training series is constant/invalid (sd = 0 or NA). Check transforms.\n")
  #     return(invisible(NULL))
  #   }
  #   
  #   # ---------- Package check ----------
  #   if (!requireNamespace("urca", quietly = TRUE) || !requireNamespace("tseries", quietly = TRUE)) {
  #     cat("==========================================================================\n")
  #     cat("                         ADF UNIT ROOT DIAGNOSTIC                         \n")
  #     cat("==========================================================================\n")
  #     cat("ERROR: Please install.packages(c('urca','tseries')) to run ADF/KPSS.\n")
  #     return(invisible(NULL))
  #   }
  #   
  #   # ---------- Common quantities ----------
  #   m_mean <- mean(x_vec, na.rm = TRUE)
  #   m_sd   <- stats::sd(x_vec, na.rm = TRUE)
  #   lb_lag <- max(1L, min(10L, floor(N / 5)))
  #   
  #   # ---------- Core ADF (urca) for chosen type & k ----------
  #   ur  <- tryCatch(urca::ur.df(x_vec, type = type_in, lags = k), error = function(e) NULL)
  #   tau_obs  <- if (!is.null(ur)) suppressWarnings(as.numeric(ur@teststat[tau_row])) else NA_real_
  #   if (!is.finite(tau_obs) && !is.null(ur)) tau_obs <- suppressWarnings(as.numeric(ur@teststat[1]))
  #   tau_crit <- if (!is.null(ur)) cval_pick_safe(ur@cval, tau_row, alpha_col) else NA_real_
  #   adf_ts   <- if (N > (k + 10)) tryCatch(tseries::adf.test(x_vec, alternative = alt_in, k = k),
  #                                          error = function(e) NULL) else NULL
  #   adf_p    <- if (!is.null(adf_ts)) to_num_safe(adf_ts$p.value) else NA_real_
  #   lb_main  <- if (!is.null(ur)) stats::Box.test(ur@res, lag = lb_lag, type = "Ljung-Box") else NULL
  #   lb_stat  <- if (!is.null(lb_main)) to_num_safe(lb_main$statistic) else NA_real_
  #   lb_p     <- if (!is.null(lb_main)) to_num_safe(lb_main$p.value) else NA_real_
  #   
  #   # ---------- KPSS ----------
  #   kpss_type <- if (type_in == "trend") "Trend" else "Level"
  #   kpss_ts   <- tryCatch(tseries::kpss.test(x_vec, null = kpss_type), error = function(e) NULL)
  #   kpss_p    <- if (!is.null(kpss_ts)) to_num_safe(kpss_ts$p.value) else NA_real_
  #   kpss_uc   <- tryCatch(urca::ur.kpss(x_vec, type = if (type_in == "trend") "tau" else "mu"),
  #                         error = function(e) NULL)
  #   eta_obs_uc  <- if (!is.null(kpss_uc)) to_num_safe(kpss_uc@teststat) else NA_real_
  #   eta_col     <- if (alpha_val <= 0.01) "1pct" else if (alpha_val <= 0.05) "5pct" else "10pct"
  #   eta_crit_uc <- if (!is.null(kpss_uc)) cval_pick_safe(kpss_uc@cval, 1, eta_col) else NA_real_
  #   kpss_by_p   <- is.finite(kpss_p)     && (kpss_p < alpha_val)
  #   kpss_by_eta <- is.finite(eta_obs_uc) && is.finite(eta_crit_uc) && (eta_obs_uc > eta_crit_uc)
  #   kpss_stationary <- !(kpss_by_p || kpss_by_eta)
  #   
  #   # ---------- Structural break (Pettitt) ----------
  #   pett_U <- pett_p <- NA_real_; pett_cp <- NA_integer_
  #   if (requireNamespace("trend", quietly = TRUE)) {
  #     pet <- tryCatch(trend::pettitt.test(x_vec), error = function(e) NULL)
  #     if (!is.null(pet)) {
  #       pett_U <- to_num_safe(pet$statistic)
  #       pett_p <- to_num_safe(pet$p.value)
  #       pett_cp <- to_int_safe(pet$estimate[1], default = NA_integer_)
  #       if (!is.finite(pett_cp)) pett_cp <- as.integer(floor(N/2))
  #       if (pett_cp < 2L || pett_cp > (N - 2L)) pett_cp <- NA_integer_
  #     }
  #   }
  #   
  #   # ---------- KPSS by segments (if break usable) ----------
  #   seg1_p <- seg2_p <- NA_real_
  #   if (is.finite(pett_cp)) {
  #     if (pett_cp >= 8 && (N - pett_cp) >= 8) {
  #       seg1 <- tryCatch(tseries::kpss.test(x_vec[seq_len(pett_cp)], null = kpss_type), error = function(e) NULL)
  #       seg2 <- tryCatch(tseries::kpss.test(x_vec[(pett_cp + 1L):N], null = kpss_type), error = function(e) NULL)
  #       if (!is.null(seg1)) seg1_p <- to_num_safe(seg1$p.value)
  #       if (!is.null(seg2)) seg2_p <- to_num_safe(seg2$p.value)
  #     }
  #   }
  #   
  #   # ---------- Phase 0: suggested k (LB scan on TRAINING window) ----------
  #   k_grid <- 0L:min(12L, max(1L, floor(N/10)))
  #   lbp_by_k <- rep(NA_real_, length(k_grid))
  #   for (i in seq_along(k_grid)) {
  #     kk <- k_grid[i]
  #     ur_i <- tryCatch(urca::ur.df(x_vec, type = type_in, lags = kk), error = function(e) NULL)
  #     lb_i <- if (!is.null(ur_i)) tryCatch(stats::Box.test(ur_i@res, lag = lb_lag, type = "Ljung-Box"), error = function(e) NULL) else NULL
  #     lbp_by_k[i] <- if (!is.null(lb_i)) to_num_safe(lb_i$p.value) else NA_real_
  #   }
  #   k_suggest <- NA_integer_
  #   idx_ok <- which(is.finite(lbp_by_k) & lbp_by_k > alpha_val)
  #   if (length(idx_ok) > 0) k_suggest <- k_grid[min(idx_ok)] else k_suggest <- k
  #   
  #   # Spec sensitivity across types for k = selected and k = k_suggest
  #   types <- c("none","drift","trend")
  #   spec_line <- function(kk) {
  #     for (tp in types) {
  #       tau_row_tp <- tau_row_for(tp)
  #       ur_tp <- tryCatch(urca::ur.df(x_vec, type = tp, lags = kk), error = function(e) NULL)
  #       tau_tp  <- if (!is.null(ur_tp)) suppressWarnings(as.numeric(ur_tp@teststat[tau_row_tp])) else NA_real_
  #       if (!is.finite(tau_tp) && !is.null(ur_tp)) tau_tp <- suppressWarnings(as.numeric(ur_tp@teststat[1]))
  #       crit_tp <- if (!is.null(ur_tp)) cval_pick_safe(ur_tp@cval, tau_row_tp, alpha_col) else NA_real_
  #       lb_tp   <- if (!is.null(ur_tp)) tryCatch(stats::Box.test(ur_tp@res, lag = lb_lag, type = "Ljung-Box"), error = function(e) NULL) else NULL
  #       lbp_tp  <- if (!is.null(lb_tp)) to_num_safe(lb_tp$p.value) else NA_real_
  #       adf_dec <- if (is.finite(tau_tp) && is.finite(crit_tp) && (tau_tp < crit_tp)) "ADF=STATIONARY" else "ADF=NON-STATIONARY"
  #       cat(sprintf("     %s type=%-5s | tau=%s | crit=%s | %s | LB p=%s\n",
  #                   tick(is.finite(lbp_tp) && lbp_tp > alpha_val),
  #                   tp,
  #                   fmt_num(tau_tp, 4),
  #                   fmt_num(crit_tp, 4),
  #                   adf_dec,
  #                   fmt_num(lbp_tp, 4)))
  #     }
  #   }
  #   
  #   # ---------- Agreement ----------
  #   adf_ok_vals <- is.finite(tau_obs) && is.finite(tau_crit)
  #   adf_stationary <- adf_ok_vals && (tau_obs < tau_crit)
  #   agreement <- (adf_stationary && kpss_stationary) || (!adf_stationary && !kpss_stationary)
  #   
  #   # ---------- HEADER ----------
  #   cat("==========================================================================\n")
  #   cat("                        ADF UNIT ROOT DIAGNOSTIC                          \n")
  #   cat("==========================================================================\n")
  #   cat(sprintf(" MODEL TYPE : %-10s | SAMPLE SIZE (N) : %s (TRAINING)\n", toupper(type_in), N))
  #   cat(sprintf(" LAG ORDER  : %-10s | SIGNIFICANCE (α) : %s\n", k, fmt_num(alpha_val, 4)))
  #   cat(sprintf(" MOMENTS    : Mean: %s | Std.Dev: %s\n",
  #               fmt_num(m_mean, 4), fmt_num(m_sd, 4)))
  #   cat("--------------------------------------------------------------------------\n")
  #   cat(" TRANSFORMATION PROVENANCE (what ALL tests use):\n")
  #   cat(sprintf(" [ ] Source object class           : %s\n", x_class))
  #   cat(sprintf(" [ ] Frequency (if ts)             : %s\n", ifelse(is.finite(x_freq), x_freq, "NA")))
  #   cat(sprintf(" [ ] NA count before na.omit       : %s\n", na_before))
  #   cat(sprintf(" [ ] Transformation inputs (UI)    : log=%s | d=%s | D=%s\n",
  #               ifelse(isTRUE(log_in),"ON","OFF"), as.character(d_in), as.character(D_in)))
  #   ht <- safe_head_tail(x_vec, 5)
  #   cat(sprintf(" [ ] First 5 values used in tests  : %s\n", paste(round(ht$head, 4), collapse = ", ")))
  #   cat(sprintf(" [ ] Last  5 values used in tests  : %s\n\n", paste(round(ht$tail, 4), collapse = ", ")))
  #   
  #   # ---------- PHASE 0 ----------
  #   cat("==========================================================================\n")
  #   cat("PHASE 0: DECISION AID (Lag + Spec Sensitivity)\n")
  #   cat("==========================================================================\n")
  #   cat(sprintf(" • Ljung-Box reference lag used in scan : %d  ; (LB lag = min(10, floor(N/5))\n", lb_lag))
  #   cat(sprintf(" [!] Suggested k* (scan)                : %s  (smallest k with Ljung-Box p-value > alpha (whiter residuals))\n",
  #               ifelse(is.finite(k_suggest), k_suggest, "NA")))
  #   if (is.finite(k_suggest) && k_suggest != k) {
  #     cat(sprintf(" [!] You selected k=%d. Consider trying k=%d and re-running.\n", k, k_suggest))
  #   } else {
  #     cat(sprintf(" [✓] You selected k=%d. This already meets the LB>α rule-of-thumb.\n", k))
  #   }
  #   cat("--------------------------------------------------------------------------\n")
  #   cat(" ADF SPEC SENSITIVITY (ur.df):\n")
  #   cat("   (Prefer: LB ok + stable decision across types)\n")
  #   cat(sprintf("  • k=%d\n", k)); spec_line(k)
  #   if (is.finite(k_suggest) && k_suggest != k) {
  #     cat(sprintf("  • k=%d\n", k_suggest)); spec_line(k_suggest)
  #   }
  #   cat("--------------------------------------------------------------------------\n")
  #   
  #   # ---------- PHASE 1 ----------
  #   cat("==========================================================================\n")
  #   cat("PHASE 1: ADF UNIT ROOT TEST\n")
  #   cat("==========================================================================\n")
  #   cat(" • H0: The series has a Unit Root (Non-Stationary).\n")
  #   cat(" • Ha: The series is Stationary (Mean Reverting).\n")
  #   cat(sprintf(" -> CRITERIA: Reject H0 if Tau-Obs (%s) < Tau-Crit (%s)\n",
  #               fmt_num(tau_obs, 4), fmt_num(tau_crit, 4)))
  #   cat("\n RESULT:\n")
  #   cat(sprintf("  - Tau Observed : %s \n", fmt_num(tau_obs, 6)))
  #   cat(sprintf("  - Tau Critical : %s \n", fmt_num(tau_crit, 2)))
  #   cat(sprintf("  - P-Value (Ref): %s \n", ifelse(is.finite(adf_p), fmt_p(adf_p), "NA")))
  #   cat("\n DECISION:\n")
  #   if (adf_ok_vals && adf_stationary) {
  #     cat("  -> REJECT H0: Evidence suggests the series is STATIONARY (on training window).\n")
  #   } else if (adf_ok_vals && !adf_stationary) {
  #     cat("  -> FAIL TO REJECT H0: Evidence suggests the series is NON-STATIONARY (on training window).\n")
  #   } else {
  #     cat("  -> INCONCLUSIVE: Missing statistic/critical value; regression may be ill-conditioned.\n")
  #   }
  #   
  #   # ---------- PHASE 2 ----------
  #   cat("\n==========================================================================\n")
  #   cat("PHASE 2: RESIDUAL DIAGNOSTICS (LJUNG-BOX)\n")
  #   cat("==========================================================================\n")
  #   cat(" • H0: Residuals are White Noise (No Autocorrelation).\n")
  #   cat(" • Ha: Residuals are Correlated (Lags are insufficient).\n")
  #   cat(sprintf(" -> CRITERIA: Reject H0 if P-Value (%s) < α (%s)\n", fmt_num(lb_p, 6), fmt_num(alpha_val, 4)))
  #   cat("\n RESULT:\n")
  #   cat(sprintf("  - LB Statistic : %s \n", fmt_num(lb_stat, 6)))
  #   cat(sprintf("  - LB P-Value   : %s \n", fmt_num(lb_p, 6)))
  #   cat(sprintf("  - LB Lag used  : %d \n", lb_lag))
  #   cat("    it’s normal to choose LB lag by a rule-of-thumb (like min(10, floor(N/5)) \n")
  #   cat("\n DECISION:\n")
  #   if (is.finite(lb_p)) {
  #     if (lb_p > alpha_val) cat("  -> FAIL TO REJECT H0: Residuals are White Noise. [ADF more reliable]\n")
  #     else                  cat("  -> REJECT H0: Residual autocorrelation remains; increase k or difference.\n")
  #   } else {
  #     cat("  -> INCONCLUSIVE: Ljung–Box p-value is NA.\n")
  #   }
  #   
  #   # ---------- PHASE 3A ----------
  #   cat("\n==========================================================================\n")
  #   cat("PHASE 3: KPSS + STRUCTURAL BREAK (Pettitt) + KPSS SEGMENTS\n")
  #   cat("==========================================================================\n\n")
  #   cat("PHASE 3A: KPSS (Stationarity Confirmation)\n")
  #   cat(sprintf(" • H0: The series is Stationary around a %s.\n", if (kpss_type=="Trend") "Trend" else "Level"))
  #   cat(" • Ha: The series is Non-Stationary.\n")
  #   cat(sprintf(" • CRITERIA (p-value) : Reject H0 if p-value (%s) < α (%s)\n", 
  #               ifelse(is.finite(kpss_p), fmt_num(kpss_p, 6), "NA"), fmt_num(alpha_val,4)))
  #   cat(sprintf(" • CRITERIA (eta)     : Reject H0 if Eta-Obs (%s) > Eta-Crit (%s)  [urca]\n", 
  #               fmt_num(eta_obs_uc, 6), fmt_num(eta_crit_uc, 6)))
  #   cat("\n RESULT:\n")
  #   cat(sprintf("  - Eta (Observed value) [urca]   : %s \n", fmt_num(eta_obs_uc, 5)))
  #   cat(sprintf("  - Eta (Critical value) [urca]   : %s \n", fmt_num(eta_crit_uc, 3)))
  #   cat(sprintf("  - p-value (one-tailed) [tseries]: %s \n", ifelse(is.finite(kpss_p), fmt_num(kpss_p, 6), "NA")))
  #   cat(sprintf("  - Eta (Observed) [tseries, FYI] : %s \n", fmt_num(eta_obs_uc, 5)))
  #   cat("\n DECISION:\n")
  #   if (kpss_stationary) cat("  -> FAIL TO REJECT H0: Stationarity supported by KPSS.\n")
  #   else                 cat("  -> REJECT H0: Non-stationarity indicated by KPSS.\n")
  #   
  #   # ---------- PHASE 3B ----------
  #   cat("--------------------------------------------------------------------------\n")
  #   cat("PHASE 3B: STRUCTURAL BREAK CHECK (Pettitt) + KPSS SEGMENT CHECK\n")
  #   cat(" • Goal: Detect a single change-point (median shift) that can distort KPSS.\n")
  #   cat(" • If a break exists, we re-run KPSS before/after the break.\n")
  #   cat(sprintf(" • Pettitt p-value: %s | α: %s\n\n", ifelse(is.finite(pett_p), fmt_num(pett_p, 6), "NA"), fmt_num(alpha_val, 4)))
  #   cat(" RESULT:\n")
  #   cat(sprintf("  - Break index (estimate): %s \n", ifelse(is.finite(pett_cp), pett_cp, "NA")))
  #   cat(sprintf("  - Pettitt statistic     : %s \n", fmt_num(pett_U, 0)))
  #   cat(sprintf("  - Pettitt p-value       : %s \n", ifelse(is.finite(pett_p), fmt_num(pett_p, 6), "NA")))
  #   cat("\n DECISION:\n")
  #   if (is.finite(pett_p)) {
  #     if (pett_p < alpha_val) cat("  -> REJECT H0: Break detected. Full-sample KPSS/ADF may be contaminated.\n")
  #     else                    cat("  -> FAIL TO REJECT H0: No strong single break detected.\n")
  #   } else {
  #     cat("  -> INCONCLUSIVE: Pettitt p-value NA.\n")
  #   }
  #   
  #   # ---------- PHASE 3C ----------
  #   cat("--------------------------------------------------------------------------\n")
  #   cat("PHASE 3C: KPSS RE-CHECK BY SEGMENTS\n")
  #   cat("  • Segment 1 = [1 .. break]\n")
  #   cat("  • Segment 2 = [break+1 .. N]\n")
  #   cat(sprintf("  - KPSS p-value (Segment 1): %s\n", ifelse(is.finite(seg1_p), fmt_num(seg1_p, 6), "NA")))
  #   cat(sprintf("  - KPSS p-value (Segment 2): %s\n", ifelse(is.finite(seg2_p), fmt_num(seg2_p, 6), "NA")))
  #   cat("\n INTERPRETATION:\n")
  #   if (is.finite(seg1_p) && is.finite(seg2_p)) {
  #     if (seg1_p >= alpha_val && seg2_p >= alpha_val)
  #       cat("  [✓] Both segments look stationary by KPSS.\n      -> Full-sample non-stationarity may be break-driven.\n")
  #     else if (seg1_p < alpha_val && seg2_p < alpha_val)
  #       cat("  [X] Both segments look non-stationary by KPSS.\n")
  #     else
  #       cat("  [?] Mixed evidence: one segment stationary, the other not.\n")
  #   } else {
  #     cat("  [?] Segment KPSS not available (short segments or missing package).\n")
  #   }
  #   
  #   # ---------- PHASE 4: Final verdict & advice ----------
  #   cat("\n==========================================================================\n")
  #   cat("PHASE 4: FINAL ACADEMIC VERDICT & ADVICE\n")
  #   cat("==========================================================================\n")
  #   if (adf_stationary && kpss_stationary) {
  #     cat(" [✓] VERDICT: CONVERGENT STATIONARITY (ADF & KPSS agree; residuals likely white).\n")
  #     cat("     ADVICE: Proceed with SARIMA identification with d=0 (choose D via seasonality), then residual checks.\n")
  #   } else if (!adf_stationary && !kpss_stationary) {
  #     cat(" [X] VERDICT: CONVERGENT NON-STATIONARITY (ADF & KPSS agree).\n")
  #     cat("     ADVICE: Difference the series (d=1). If seasonal (m≥2), consider D=1, then re-run tests.\n")
  #   } else {
  #     cat(" [?] VERDICT: CONFLICTING RESULTS (ADF vs KPSS).\n")
  #     cat("     ADVICE: Near-unit-root, trend-vs-drift mismatch, or break contamination are common causes.\n")
  #     cat("             Use PHASE 0 to pick k/type with LB ok and stable decisions; consider seasonal differencing and break handling.\n")
  #   }
  #   
  #   # ---------- Technical appendix: ADF regression coefficients ----------
  #   cat("\n TECHNICAL APPENDIX (ADF Regression Coefficients):\n")
  #   if (!is.null(ur) && !is.null(ur@testreg)) {
  #     cf <- tryCatch(coef(summary(ur@testreg)), error = function(e) NULL)
  #     if (!is.null(cf)) {
  #       printCoefmat(cf, digits = 7, signif.stars = FALSE, P.values = TRUE, has.Pvalue = TRUE)
  #     } else {
  #       cat("  (coefficients unavailable)\n")
  #     }
  #   } else {
  #     cat("  (ADF regression unavailable)\n")
  #   }
  #   
  #   # ---------- PHASE 5: Evidence snapshot & checklist ----------
  #   cat("\n==========================================================================\n")
  #   cat("PHASE 5: POST-SUMMARY (Academic-quality snapshot)\n")
  #   cat("==========================================================================\n")
  #   cat(" EVIDENCE SNAPSHOT (All key outcomes in one place):\n")
  #   cat(sprintf(" [ ] N (effective sample size)              : %s\n", N))
  #   cat(sprintf(" [ ] Model type (ADF)                       : %s  (tau row: %s)\n", tolower(type_in), tau_row))
  #   cat(sprintf(" [ ] Lag order (k)                          : %s\n", k))
  #   cat(sprintf(" [ ] Alpha (α)                              : %s\n", fmt_num(alpha_val, 4)))
  #   cat(sprintf(" [ ] Tau-Observed (urca)                    : %s\n", fmt_num(tau_obs, 6)))
  #   cat(sprintf(" [ ] Tau-Critical (urca, %s)                  : %s\n", alpha_col, fmt_num(tau_crit, 6)))
  #   cat(sprintf(" [ ] ADF p-value (tseries reference)        : %s\n", ifelse(is.finite(adf_p), fmt_num(adf_p, 6), "NA")))
  #   cat(sprintf(" [ ] Ljung-Box p-value (residuals)          : %s\n", fmt_num(lb_p, 6)))
  #   cat(sprintf(" [ ] KPSS Eta observed (urca)               : %s\n", fmt_num(eta_obs_uc, 6)))
  #   cat(sprintf(" [ ] KPSS Eta critical (urca, %s)             : %s\n", eta_col, fmt_num(eta_crit_uc, 6)))
  #   cat(sprintf(" [ ] KPSS p-value (one-tailed, tseries)     : %s\n", ifelse(is.finite(kpss_p), fmt_num(kpss_p, 6), "NA")))
  #   cat(sprintf(" [!] Suggested k* (scan)                     : %s\n", ifelse(is.finite(k_suggest), k_suggest, "NA")))
  #   cat("--------------------------------------------------------------------------\n")
  #   cat(" CHECKLIST (Academic-quality acceptance criteria):\n")
  #   cat(sprintf(" [ ] Tests run on TRAINING window only            : %s YES\n", tick(TRUE)))
  #   cat(sprintf(" [ ] Variance usable (sd>0)                       : %s SATISFIED\n", tick(TRUE)))
  #   cat(sprintf(" [ ] Tau is finite and usable                     : %s %s\n", tick(is.finite(tau_obs)), ifelse(is.finite(tau_obs),"SATISFIED","CHECK")))
  #   cat(sprintf(" [ ] Tau critical value extracted                 : %s %s\n", tick(is.finite(tau_crit)), ifelse(is.finite(tau_crit),"SATISFIED","CHECK")))
  #   cat(sprintf(" [ ] Residuals pass Ljung-Box (white noise)       : %s %s\n", tick(is.finite(lb_p) && lb_p > alpha_val), ifelse(is.finite(lb_p) && lb_p > alpha_val,"SATISFIED","CHECK")))
  #   cat(sprintf(" [ ] KPSS eta & critical are usable               : %s %s\n", tick(is.finite(eta_obs_uc) && is.finite(eta_crit_uc)), ifelse(is.finite(eta_obs_uc) && is.finite(eta_crit_uc),"SATISFIED","CHECK")))
  #   cat(sprintf(" [ ] KPSS p-value is usable                       : %s %s\n", tick(is.finite(kpss_p)), ifelse(is.finite(kpss_p),"SATISFIED","CHECK")))
  #   cat(sprintf(" [ ] Sample size adequacy                         : %s STRONG\n", tick(N >= 50)))
  #   cat(sprintf(" [ ] Lag order reasonable relative to N           : %s %s\n", tick(k <= max(12L, floor(N/10))), ifelse(k <= max(12L, floor(N/10)),"OK","LARGE")))
  #   cat(sprintf(" [ ] Seasonal differencing indicated (UI D>0)      : %s %s\n", warn(isTRUE(to_int_safe(D_in,0L) > 0)), ifelse(isTRUE(to_int_safe(D_in,0L) > 0),"YES","NO / UNKNOWN")))
  #   cat(sprintf(" [ ] ADF alternative mode (Stationary/Explosive)   : %s\n", paste0(ifelse(is.null(alt_in),"NA",alt_in))))
  #   cat(sprintf(" [ ] ADF & KPSS agreement                          : %s %s\n",
  #               qmark(agreement), ifelse(agreement,"AGREEMENT","CONFLICT")))
  #   cat("     [?] NOTE: conflicts are common with near-unit-root series, trend vs drift mismatch,\n")
  #   cat("         structural breaks (Pettitt), or missing seasonal differencing.\n")
  #   cat("--------------------------------------------------------------------------\n")
  #   cat(" STRUCTURAL BREAK SUMMARY (Pettitt):\n")
  #   cat(sprintf(" [ ] Pettitt p-value : %s  (Reject H0 if < α)\n", ifelse(is.finite(pett_p), fmt_num(pett_p, 6), "NA")))
  #   cat("--------------------------------------------------------------------------\n")
  #   cat("==========================================================================\n")
  #   cat(" ACTIONABLE NEXT STEPS (What to do now):\n")
  #   cat("==========================================================================\n")
  #   cat(sprintf(" %s [1] RESIDUAL AUTOCORRELATION CHECK\n", tick(is.finite(lb_p) && lb_p > alpha_val)))
  #   cat("     • Ljung-Box is acceptable → ADF regression is less likely biased.\n")
  #   cat(" [?] [2] RESOLVE ADF vs KPSS CONFLICT\n")
  #   cat("     • Use PHASE 0 'ADF SPEC SENSITIVITY' to pick the type with LB ok and stable decision.\n")
  #   cat("     • Try ADF model type variants: none / drift / trend (match KPSS Level vs Trend).\n")
  #   cat("     • If Pettitt indicates a break: split sample and re-test.\n")
  #   cat("     • If series is seasonal: test after seasonal differencing (D=1) + maybe log.\n")
  #   cat("     • Consider variance stabilization: log or Box-Cox (if positive data).\n")
  #   cat(" [!] [3] SEASONALITY SANITY (especially for AirPassengers-like series)\n")
  #   if (is.finite(x_freq) && x_freq >= 2) {
  #     cat(sprintf("     • Detected frequency=%d → seasonality is plausible.\n", x_freq))
  #   } else {
  #     cat("     • Frequency unknown → inspect ACF/PACF for seasonal spikes.\n")
  #   }
  #   cat(" [✓] [4] EXPLOSIVE MODE NOTE\n")
  #   cat(sprintf("     • %s\n", ifelse(identical(alt_in,"explosive"),
  #                                     "Explosive alternative was requested; interpret ADF accordingly.",
  #                                     "Not in explosive mode → standard stationarity workflow applies.")))
  #   cat("\n PRACTICAL MODELING PATH (for your Shiny workflow):\n")
  #   if (!adf_stationary || !kpss_stationary) {
  #     cat(" [X] Apply differencing (d and/or D) → re-run ADF/KPSS → then identify ARMA.\n")
  #   } else {
  #     cat(" [✓] Keep d=0 (and decide D via seasonality) → SARIMA identification and residual diagnostics.\n")
  #   }
  #   cat("--------------------------------------------------------------------------\n")
  # })
  
  
  
  

  
  # output$teststationarited3St <- renderPrint({
  #   # ---------- Helpers ----------
  #   `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
  #   to_num_safe <- function(v, default = NA_real_) {
  #     out <- suppressWarnings(as.numeric(v))
  #     if (length(out) == 0 || all(is.na(out)) || !is.finite(out[1])) default else out[1]
  #   }
  #   to_int_safe <- function(v, default = 0L) {
  #     out <- suppressWarnings(as.integer(v))
  #     if (length(out) == 0 || is.na(out[1]) || !is.finite(out[1])) default else out[1]
  #   }
  #   safe_head_tail <- function(x, n = 5) {
  #     x <- as.numeric(x); x <- x[is.finite(x)]
  #     if (length(x) == 0) return(list(head = numeric(0), tail = numeric(0)))
  #     list(head = head(x, n), tail = tail(x, n))
  #   }
  #   tau_row_for <- function(type_in) switch(type_in, "none"="tau1", "drift"="tau2", "trend"="tau3", "tau3")
  #   title <- function(txt) cat(sprintf("\n%s\n", txt))
  #   bullet <- function(txt) cat(sprintf(" • %s\n", txt))
  #   line <- function(k, v) cat(sprintf("   - %-22s %s\n", paste0(k, ":"), v))
  #   
  #   # Critical value extractor (robust to naming)
  #   cval_pick_safe <- function(cval_obj, row, col) {
  #     if (is.null(cval_obj)) return(NA_real_)
  #     ans <- NA_real_
  #     if (is.matrix(cval_obj)) {
  #       if (!missing(row) && !missing(col) && row %in% rownames(cval_obj) && col %in% colnames(cval_obj)) {
  #         ans <- suppressWarnings(as.numeric(cval_obj[row, col]))
  #       } else if (!missing(col) && col %in% colnames(cval_obj)) {
  #         ans <- suppressWarnings(as.numeric(cval_obj[1, col]))
  #       } else {
  #         ans <- suppressWarnings(as.numeric(cval_obj[1, 1]))
  #       }
  #     } else {
  #       nm <- names(cval_obj)
  #       if (!is.null(nm) && col %in% nm) ans <- suppressWarnings(as.numeric(cval_obj[[col]]))
  #       if (!is.finite(ans) && length(cval_obj) >= 3) {
  #         if (identical(col, "10pct")) ans <- suppressWarnings(as.numeric(cval_obj[1]))
  #         if (identical(col, "5pct"))  ans <- suppressWarnings(as.numeric(cval_obj[2]))
  #         if (identical(col, "1pct"))  ans <- suppressWarnings(as.numeric(cval_obj[3]))
  #       }
  #       if (!is.finite(ans)) ans <- suppressWarnings(as.numeric(cval_obj[1]))
  #     }
  #     ans
  #   }
  #   
  #   # APA helpers
  #   fmt_p <- function(p) {
  #     if (!is.finite(p)) return("p = NA")
  #     if (p < .001) "p < .001" else paste0("p = ", sub("^0\\.", ".", sprintf("%.3f", p)))
  #   }
  #   apa_line_adf  <- function(tau, p, k, N, type)
  #     paste0("According to the ADF test (", type, ", k = ", k, ", N = ", N, "), the observed statistic was τ = ",
  #            sprintf("%.3f", tau), " with ", fmt_p(p), ".")
  #   apa_line_pp   <- function(stat, p)
  #     paste0("According to the Phillips–Perron test, the observed statistic was Z_α = ",
  #            ifelse(is.finite(stat), sprintf("%.3f", stat), "NA"), " with ", fmt_p(p), ".")
  #   apa_line_dfgls<- function(stat, crit)
  #     paste0("The DF-GLS statistic was ", ifelse(is.finite(stat), sprintf("%.3f", stat), "NA"),
  #            if (is.finite(crit)) paste0(", compared to a 5% critical value of ", sprintf("%.3f", crit), ".") else ".")
  #   apa_line_kpss <- function(eta, p, mode)
  #     paste0("According to the KPSS test (", mode, "), the observed statistic was η = ",
  #            sprintf("%.3f", eta), " with ", fmt_p(p), ".")
  #   apa_line_lb   <- function(stat, p, L)
  #     paste0("Using a Ljung–Box test at lag ", L, ", the statistic was Q = ",
  #            ifelse(is.finite(stat), sprintf("%.3f", stat), "NA"), " with ", fmt_p(p), ".")
  #   apa_line_pett <- function(U, p)
  #     paste0("The Pettitt test yielded U = ", ifelse(is.finite(U), sprintf("%.3f", U), "NA"),
  #            " with ", fmt_p(p), ".")
  #   apa_line_jb   <- function(stat, p)
  #     paste0("The Jarque–Bera normality test returned JB = ",
  #            ifelse(is.finite(stat), sprintf("%.3f", stat), "NA"), " with ", fmt_p(p), ".")
  #   apa_line_arch <- function(stat, p, L)
  #     paste0("The ARCH LM test (lags = ", L, ") returned LM = ",
  #            ifelse(is.finite(stat), sprintf("%.3f", stat), "NA"), " with ", fmt_p(p), ".")
  #   
  #   # ---------- Inputs (from your UI) ----------
  #   alt_in  <- input$alternd2St %||% input$alternSt       # "stationary", "explosive", "regression" (ADF/PP)
  #   lag_in  <- input$LagOrderADFd2St %||% input$LagOrderADFSt
  #   a_in    <- input$alphaSt2
  #   type_in <- input$adfTypeSt2                            # "none", "drift", "trend"
  #   
  #   # Transformation flags (provenance)
  #   d_in   <- input$d_n  %||% NA
  #   D_in   <- input$DS_n %||% NA
  #   log_in <- input$check_box %||% FALSE
  #   
  #   # ---------- Data ----------
  #   req(myData_Choice())
  #   x_raw <- myData_Choice()
  #   na_before <- sum(is.na(x_raw))
  #   x_class <- paste(class(x_raw), collapse = ", ")
  #   x_freq  <- if (inherits(x_raw, "ts")) tryCatch(stats::frequency(x_raw), error = function(e) NA_integer_) else NA_integer_
  #   x <- as.numeric(stats::na.omit(x_raw))
  #   N <- length(x)
  #   
  #   # ---------- Hyper-parameters ----------
  #   k <- to_int_safe(lag_in, default = 0L); if (!is.finite(k) || k < 0) k <- 0L
  #   alpha_raw <- as.character(a_in)
  #   alpha_val <- to_num_safe(alpha_raw, default = 0.05)
  #   alpha_col <- switch(alpha_raw, "0.01"="1pct", "0.05"="5pct", "0.1"="10pct", "0.10"="10pct",
  #                       "1pct"="1pct", "5pct"="5pct", "10pct"="10pct", "5pct")
  #   tau_row <- tau_row_for(type_in)
  #   seasonality_resolved <- isTRUE(to_int_safe(D_in, 0L) > 0L)
  #   
  #   # ---------- Sanity checks ----------
  #   if (N < 5) { title("Critical error"); bullet("Too few observations (N < 5). Provide more data."); return(invisible(NULL)) }
  #   if (!is.finite(stats::sd(x)) || stats::sd(x) == 0) {
  #     title("Critical error"); bullet("Series is constant/invalid (sd = 0 or NA). Check transforms."); return(invisible(NULL))
  #   }
  #   
  #   # ---------- Package check ----------
  #   if (!requireNamespace("urca", quietly = TRUE) || !requireNamespace("tseries", quietly = TRUE)) {
  #     title("Missing packages")
  #     bullet("Please install.packages(c('urca','tseries')) to run ADF/KPSS/PP.")
  #     return(invisible(NULL))
  #   }
  #   
  #   # ---------- ADF (urca + tseries p-value) ----------
  #   ur  <- tryCatch(urca::ur.df(x, type = type_in, lags = k), error = function(e) NULL)
  #   tau_obs  <- if (!is.null(ur)) suppressWarnings(as.numeric(ur@teststat[tau_row])) else NA_real_
  #   if (!is.finite(tau_obs) && !is.null(ur)) tau_obs <- suppressWarnings(as.numeric(ur@teststat[1]))
  #   tau_crit <- if (!is.null(ur)) cval_pick_safe(ur@cval, tau_row, alpha_col) else NA_real_
  #   
  #   adf_ts <- if (N > (k + 10)) tryCatch(tseries::adf.test(x, alternative = alt_in, k = k), error = function(e) NULL) else NULL
  #   adf_p  <- if (!is.null(adf_ts)) to_num_safe(adf_ts$p.value) else NA_real_
  #   
  #   # ---------- Phillips–Perron (robust to serial corr.) ----------
  #   alt_pp <- if (alt_in %in% c("stationary","explosive")) alt_in else "stationary"
  #   pp_res <- tryCatch(tseries::pp.test(x, alternative = alt_pp), error = function(e) NULL)
  #   pp_stat<- if (!is.null(pp_res)) to_num_safe(pp_res$statistic[1]) else NA_real_
  #   pp_p   <- if (!is.null(pp_res)) to_num_safe(pp_res$p.value) else NA_real_
  #   
  #   # ---------- DF-GLS / ERS (more power) ----------
  #   ers_res <- tryCatch(urca::ur.ers(x, type = "DF-GLS", model = if (type_in == "trend") "trend" else "constant",
  #                                    lag.max = max(k, 1L)), error = function(e) NULL)
  #   ers_stat <- if (!is.null(ers_res)) to_num_safe(ers_res@teststat) else NA_real_
  #   ers_crit5<- if (!is.null(ers_res)) cval_pick_safe(ers_res@cval, 1, "5pct") else NA_real_
  #   ers_stationary <- is.finite(ers_stat) && is.finite(ers_crit5) && (ers_stat < ers_crit5)
  #   
  #   # ---------- Ljung–Box on ADF residuals ----------
  #   lb_lag <- max(1L, min(10L, floor(N / 5)))
  #   lb     <- if (!is.null(ur)) stats::Box.test(ur@res, lag = lb_lag, type = "Ljung-Box") else NULL
  #   lb_stat<- if (!is.null(lb)) to_num_safe(lb$statistic) else NA_real_
  #   lb_p   <- if (!is.null(lb)) to_num_safe(lb$p.value) else NA_real_
  #   
  #   # ---------- KPSS ----------
  #   kpss_type <- if (type_in == "trend") "Trend" else "Level"
  #   kpss_ts   <- tryCatch(tseries::kpss.test(x, null = kpss_type), error = function(e) NULL)
  #   kpss_p    <- if (!is.null(kpss_ts)) to_num_safe(kpss_ts$p.value) else NA_real_
  #   
  #   kpss_uc   <- tryCatch(urca::ur.kpss(x, type = if (type_in == "trend") "tau" else "mu"), error = function(e) NULL)
  #   eta_obs_uc  <- if (!is.null(kpss_uc)) to_num_safe(kpss_uc@teststat) else NA_real_
  #   eta_col     <- if (alpha_val <= 0.01) "1pct" else if (alpha_val <= 0.05) "5pct" else "10pct"
  #   eta_crit_uc <- if (!is.null(kpss_uc)) cval_pick_safe(kpss_uc@cval, 1, eta_col) else NA_real_
  #   
  #   # ---------- Pettitt (optional structural break) ----------
  #   pett_U <- pett_p <- NA_real_
  #   if (requireNamespace("trend", quietly = TRUE)) {
  #     pet <- tryCatch(trend::pettitt.test(x), error = function(e) NULL)
  #     if (!is.null(pet)) { pett_U <- to_num_safe(pet$statistic); pett_p <- to_num_safe(pet$p.value) }
  #   }
  #   
  #   # ---------- Seasonal ACF summary (for SARIMA students) ----------
  #   seasonal_summary <- NULL
  #   if (is.finite(x_freq) && x_freq >= 2) {
  #     ac <- tryCatch(stats::acf(x, lag.max = min(2L * x_freq, max(10L, x_freq)), plot = FALSE), error = function(e) NULL)
  #     if (!is.null(ac) && length(ac$acf) > 1) {
  #       # acf includes lag 0; build a small summary
  #       lags <- as.integer(round(ac$lag * x_freq))  # convert fractional lags if ts
  #       vals <- as.numeric(ac$acf)
  #       idx  <- lags >= 1
  #       lags <- lags[idx]; vals <- vals[idx]
  #       seas_lags <- lags[lags %in% c(x_freq, 2L * x_freq)]
  #       seas_vals <- vals[lags %in% c(x_freq, 2L * x_freq)]
  #       thr <- 1.96 / sqrt(N)
  #       seasonal_summary <- list(threshold = thr,
  #                                s = x_freq,
  #                                pairs = if (length(seas_lags) > 0) data.frame(lag = seas_lags, acf = round(seas_vals, 3)) else NULL,
  #                                top = { o <- order(-abs(vals)); head(data.frame(lag = lags[o], acf = round(vals[o], 3)), 5) })
  #     }
  #   }
  #   
  #   # ---------- Optional: recommended differencing (forecast) ----------
  #   nd <- nsd <- NA_integer_
  #   if (requireNamespace("forecast", quietly = TRUE)) {
  #     nd  <- tryCatch(forecast::ndiffs(x, alpha = alpha_val, test = "kpss"), error = function(e) NA_integer_)
  #     nsd <- if (is.finite(x_freq) && x_freq >= 2)
  #       tryCatch(forecast::nsdiffs(x, m = x_freq, test = "ch"), error = function(e) NA_integer_) else NA_integer_
  #   }
  #   
  #   # ---------- Optional: residual distribution/variance checks ----------
  #   jb_stat <- jb_p <- NA_real_
  #   arch_stat <- arch_p <- NA_real_
  #   if (is.finite(lb_p)) { # only if we ran an ADF regression and got residuals whiteness test
  #     # Use the same residuals (ur@res) for these checks if available
  #     if (!is.null(ur)) {
  #       # Normality of the series (simple) — many syllabi still show JB on raw data
  #       jbr <- tryCatch(tseries::jarque.bera.test(x), error = function(e) NULL)
  #       if (!is.null(jbr)) { jb_stat <- to_num_safe(jbr$statistic[1]); jb_p <- to_num_safe(jbr$p.value) }
  #       # ARCH LM on series (pedagogical; strictly, residuals from a fitted model are better)
  #       if (requireNamespace("FinTS", quietly = TRUE)) {
  #         ar <- tryCatch(FinTS::ArchTest(x, lags = lb_lag), error = function(e) NULL)
  #         if (!is.null(ar)) { arch_stat <- to_num_safe(ar$statistic[1]); arch_p <- to_num_safe(ar$p.value) }
  #       }
  #     }
  #   }
  #   
  #   # ---------- Decisions ----------
  #   adf_ok_vals <- is.finite(tau_obs) && is.finite(tau_crit)
  #   adf_stationary <- adf_ok_vals && (tau_obs < tau_crit)
  #   
  #   kpss_by_p   <- is.finite(kpss_p)     && (kpss_p < alpha_val)
  #   kpss_by_eta <- is.finite(eta_obs_uc) && is.finite(eta_crit_uc) && (eta_obs_uc > eta_crit_uc)
  #   kpss_stationary <- !(kpss_by_p || kpss_by_eta)
  #   
  #   agreement <- (adf_stationary && kpss_stationary) || (!adf_stationary && !kpss_stationary)
  #   lb_pass   <- is.finite(lb_p) && (lb_p > alpha_val)
  #   
  #   # ---------- CHECKLIST ----------
  #   title("Checklist")
  #   line("ADF decision", if (adf_ok_vals && adf_stationary) "STATIONARY" else "NON-STATIONARY")
  #   line("KPSS decision", if (kpss_stationary) "STATIONARY" else "NON-STATIONARY")
  #   line("ADF vs KPSS", if (agreement) "AGREEMENT" else "CONFLICT")
  #   line("Residual whiteness (LB)", if (lb_pass) "PASS" else if (is.finite(lb_p)) "FAIL" else "NA")
  #   line("Seasonal differencing in UI (D>0)", if (seasonality_resolved) "YES" else "NO/UNKNOWN")
  #   
  #   # ---------- SPECIFICATION & PROVENANCE ----------
  #   title("Specification & provenance")
  #   line("Model type (ADF)", toupper(type_in))
  #   line("Lag k", k)
  #   line("Alpha (α)", alpha_raw)
  #   line("Series class", x_class)
  #   line("Frequency (m)", ifelse(is.finite(x_freq), x_freq, "NA"))
  #   line("NA before omit", na_before)
  #   ht <- safe_head_tail(x, 5)
  #   line("Transforms", paste0("log=", ifelse(isTRUE(log_in), "ON", "OFF"),
  #                             ", d=", as.character(d_in), ", D=", as.character(D_in)))
  #   line("First 5 used", paste(round(ht$head, 4), collapse = ", "))
  #   line("Last 5 used",  paste(round(ht$tail, 4), collapse = ", "))
  #   
  #   # ---------- PHASE 1: ADF ----------
  #   title("ADF: Augmented Dickey–Fuller")
  #   bullet("Purpose: Detect a nonseasonal unit root (random-walk behavior). This helps decide if nonseasonal differencing (d) is required.")
  #   bullet("Hypotheses:")
  #   bullet("  H0 — The series has a unit root (non-stationary).")
  #   bullet("  Ha — The series is stationary under the chosen deterministic terms.")
  #   bullet(sprintf("Decision rule: at α = %s, reject H0 if τ_obs < τ_crit.", alpha_raw))
  #   bullet(sprintf("Statistics: τ_obs = %.6f; τ_crit(%s) = %.6f; p(ref, tseries) = %s",
  #                  tau_obs, alpha_col, tau_crit,
  #                  if (is.finite(adf_p)) format.pval(adf_p, digits = 4) else "NA"))
  #   cat(" ", apa_line_adf(tau_obs, adf_p, k, N, type_in), "\n", sep = "")
  #   if (adf_ok_vals && adf_stationary) {
  #     bullet(sprintf("Decision: REJECT H0 (stationary) because τ_obs (%.6f) < τ_crit(%s) (%.6f).", tau_obs, alpha_col, tau_crit))
  #   } else if (adf_ok_vals && !adf_stationary) {
  #     bullet(sprintf("Decision: FAIL TO REJECT H0 (non-stationary) because τ_obs (%.6f) ≥ τ_crit(%s) (%.6f).", tau_obs, alpha_col, tau_crit))
  #   } else {
  #     bullet("Decision: INCONCLUSIVE — missing statistic/critical value (often due to large k vs N or constant series).")
  #   }
  #   
  #   # ---------- PHASE 1b: Phillips–Perron ----------
  #   title("Phillips–Perron (complement to ADF)")
  #   bullet("Purpose: Same null as ADF but corrects for serial correlation and heteroskedasticity in the errors via nonparametric adjustments.")
  #   bullet("Hypotheses: H0 — unit root; Ha — stationary (same as ADF).")
  #   bullet(sprintf("Statistics: Z_α = %s; p = %s",
  #                  ifelse(is.finite(pp_stat), sprintf("%.6f", pp_stat), "NA"),
  #                  if (is.finite(pp_p)) format.pval(pp_p, digits = 4) else "NA"))
  #   cat(" ", apa_line_pp(pp_stat, pp_p), "\n", sep = "")
  #   if (is.finite(pp_p)) {
  #     bullet(if (pp_p < alpha_val) "Decision: REJECT H0 — PP supports stationarity."
  #            else "Decision: FAIL TO REJECT H0 — PP supports a unit root.")
  #   } else bullet("Decision: INCONCLUSIVE — PP p-value not available.")
  #   
  #   # ---------- PHASE 1c: DF-GLS / ERS ----------
  #   title("DF-GLS (ERS) — higher power near unit roots")
  #   bullet("Purpose: A more powerful unit-root test obtained by GLS detrending before the ADF regression.")
  #   bullet(sprintf("Statistic: DF-GLS = %s; 5%% critical ≈ %s",
  #                  ifelse(is.finite(ers_stat), sprintf("%.6f", ers_stat), "NA"),
  #                  ifelse(is.finite(ers_crit5), sprintf("%.6f", ers_crit5), "NA")))
  #   cat(" ", apa_line_dfgls(ers_stat, ers_crit5), "\n", sep = "")
  #   if (is.finite(ers_stat) && is.finite(ers_crit5)) {
  #     bullet(if (ers_stationary) "Decision: REJECT H0 — DF-GLS indicates stationarity."
  #            else "Decision: FAIL TO REJECT H0 — DF-GLS indicates a unit root.")
  #   }
  #   
  #   # ---------- PHASE 2: Residual whiteness (LB) ----------
  #   title("Residual whiteness (Ljung–Box)")
  #   bullet("Purpose: Check that ADF regression residuals behave like white noise (helps validate the regression).")
  #   bullet("Hypotheses: H0 — residuals are white noise; Ha — residuals are autocorrelated.")
  #   bullet(sprintf("Statistics: lag = %d; Q = %s; p = %s",
  #                  lb_lag,
  #                  ifelse(is.finite(lb_stat), sprintf("%.6f", lb_stat), "NA"),
  #                  if (is.finite(lb_p)) format.pval(lb_p, digits = 4) else "NA"))
  #   cat(" ", apa_line_lb(lb_stat, lb_p, lb_lag), "\n", sep = "")
  #   if (is.finite(lb_p)) {
  #     bullet(if (lb_p > alpha_val) "Decision: FAIL TO REJECT H0 — residuals look white (favorable)."
  #            else "Decision: REJECT H0 — residual autocorrelation remains; increase k or difference.")
  #   } else bullet("Decision: INCONCLUSIVE — Ljung–Box p-value is NA.")
  #   
  #   # ---------- PHASE 3: KPSS ----------
  #   title("KPSS — confirm stationarity (opposite null)")
  #   bullet(sprintf("Purpose: Cross-check stationarity; H0 — stationary around a %s; Ha — non-stationary.",
  #                  if (kpss_type == "Trend") "trend" else "level"))
  #   bullet(sprintf("Statistics: η_obs = %.6f; η_crit(%s) = %.6f; p = %s",
  #                  eta_obs_uc, eta_col, eta_crit_uc,
  #                  if (is.finite(kpss_p)) format.pval(kpss_p, digits = 4) else "NA"))
  #   cat(" ", apa_line_kpss(eta_obs_uc, kpss_p, if (kpss_type == "Trend") "Trend" else "Level"), "\n", sep = "")
  #   bullet(if ( (is.finite(kpss_p) && kpss_p >= alpha_val) ||
  #               (is.finite(eta_obs_uc) && is.finite(eta_crit_uc) && eta_obs_uc <= eta_crit_uc)) {
  #     "Decision: FAIL TO REJECT H0 — stationarity supported."
  #   } else {
  #     "Decision: REJECT H0 — non-stationarity indicated."
  #   })
  #   
  #   # ---------- PHASE 4: Structural break (Pettitt) ----------
  #   title("Structural break (Pettitt)")
  #   bullet("Purpose: Single change-point can distort unit-root tests and suggest higher differencing than needed.")
  #   if (requireNamespace("trend", quietly = TRUE)) {
  #     bullet(sprintf("Statistics: U = %s; p = %s",
  #                    ifelse(is.finite(pett_U), sprintf("%.6f", pett_U), "NA"),
  #                    if (is.finite(pett_p)) format.pval(pett_p, digits = 4) else "NA"))
  #     cat(" ", apa_line_pett(pett_U, pett_p), "\n", sep = "")
  #     if (is.finite(pett_p)) {
  #       bullet(if (pett_p < alpha_val) "Decision: REJECT H0 — evidence of a structural break; consider segment-wise testing."
  #              else "Decision: FAIL TO REJECT H0 — no strong single break detected.")
  #     } else bullet("Decision: INCONCLUSIVE — Pettitt p-value not available.")
  #   } else bullet("Package 'trend' not installed — Pettitt skipped. (install.packages('trend'))")
  #   
  #   # ---------- SARIMA-oriented helpers ----------
  #   title("SARIMA helpers for students")
  #   if (is.finite(x_freq) && x_freq >= 2 && !is.null(seasonal_summary)) {
  #     bullet(sprintf("Seasonal ACF (m = %d). 95%% no-corr threshold ≈ ±%.3f.", seasonal_summary$s, seasonal_summary$threshold))
  #     if (!is.null(seasonal_summary$pairs) && nrow(seasonal_summary$pairs) > 0) {
  #       for (i in seq_len(nrow(seasonal_summary$pairs))) {
  #         bullet(sprintf("ACF at lag %d: %.3f", seasonal_summary$pairs$lag[i], seasonal_summary$pairs$acf[i]))
  #       }
  #     }
  #     bullet("Largest absolute autocorrelations (first five):")
  #     if (!is.null(seasonal_summary$top)) {
  #       apply(seasonal_summary$top, 1, function(r) bullet(sprintf("lag %s: %s", r[1], r[2])))
  #     }
  #   } else {
  #     bullet("Seasonal ACF: frequency not available or < 2.")
  #   }
  #   
  #   if (requireNamespace("forecast", quietly = TRUE)) {
  #     bullet(sprintf("Recommended nonseasonal differencing d (forecast::ndiffs) = %s", ifelse(is.finite(nd), nd, "NA")))
  #     if (is.finite(x_freq) && x_freq >= 2)
  #       bullet(sprintf("Recommended seasonal differencing D (forecast::nsdiffs) = %s", ifelse(is.finite(nsd), nsd, "NA")))
  #   } else {
  #     bullet("Install 'forecast' to see recommended d and D via (n)sdiffs.")
  #   }
  #   
  #   # ---------- Distribution/variance checks (pedagogical) ----------
  #   title("Distribution & variance checks (pedagogical)")
  #   bullet("Purpose: These do not decide differencing, but help explain residual diagnostics students will use after fitting SARIMA.")
  #   bullet(sprintf("Jarque–Bera normality: %s", if (is.finite(jb_stat)) apa_line_jb(jb_stat, jb_p) else "not available"))
  #   if (requireNamespace("FinTS", quietly = TRUE))
  #     bullet(sprintf("ARCH LM (volatility clustering): %s",
  #                    if (is.finite(arch_stat)) apa_line_arch(arch_stat, arch_p, lb_lag) else "not available"))
  #   else bullet("Install 'FinTS' to run an ARCH LM test (volatility).")
  #   
  #   # ---------- Final academic advice ----------
  #   title("Final academic advice")
  #   if (adf_stationary && kpss_stationary && lb_pass) {
  #     bullet("Convergent evidence of stationarity with white residuals. For SARIMA, set d = 0 (and consider D based on seasonality) and proceed to ACF/PACF model identification.")
  #   } else if (!adf_stationary && !kpss_stationary) {
  #     bullet("Both ADF and KPSS indicate non-stationarity. Difference the series (try d = 1); if m ≥ 2 and seasonal ACF is strong at lag m, consider D = 1. Re-check tests before fitting SARIMA.")
  #   } else {
  #     bullet("ADF and KPSS conflict. Try: (i) a different ADF deterministic term (none/drift/trend), (ii) addressing seasonality (D = 1), (iii) checking for breaks (Pettitt) and re-testing, (iv) DF-GLS/PP for robustness.")
  #   }
  #   
  #   # ---------- Snapshot ----------
  #   title("Snapshot")
  #   line("N", N)
  #   line("type", type_in)
  #   line("k", k)
  #   line("α", alpha_raw)
  #   line("freq (m)", ifelse(is.finite(x_freq), x_freq, "NA"))
  #   line("log", ifelse(isTRUE(log_in), "ON", "OFF"))
  #   line("d (UI)", as.character(d_in))
  #   line("D (UI)", as.character(D_in))
  #   line("ADF", sprintf("τ=%.4f, τ_crit(%s)=%.4f, %s", tau_obs, alpha_col, tau_crit, fmt_p(adf_p)))
  #   line("PP",  sprintf("%s, %s",
  #                       ifelse(is.finite(pp_stat), paste0("Z_α=", sprintf('%.4f', pp_stat)), "Z_α=NA"),
  #                       fmt_p(pp_p)))
  #   line("DF-GLS", if (is.finite(ers_stat)) sprintf("stat=%.4f; 5%% crit=%s", ers_stat,
  #                                                   ifelse(is.finite(ers_crit5), sprintf('%.4f', ers_crit5), "NA")) else "NA")
  #   line("KPSS", sprintf("η=%.4f, η_crit(%s)=%.4f, %s", eta_obs_uc, eta_col, eta_crit_uc, fmt_p(kpss_p)))
  #   line("LB", sprintf("%s (lag %d)", fmt_p(lb_p), lb_lag))
  #   if (requireNamespace("trend", quietly = TRUE)) line("Pettitt", fmt_p(pett_p))
  #   if (requireNamespace("forecast", quietly = TRUE)) {
  #     line("ndiffs (d)", ifelse(is.finite(nd), nd, "NA"))
  #     if (is.finite(x_freq) && x_freq >= 2) line("nsdiffs (D)", ifelse(is.finite(nsd), nsd, "NA"))
  #   }
  #   
  #   # ---------- Learning notes (concise) ----------
  #   title("Learning notes")
  #   bullet("Why multiple tests? ADF can lose power with autocorrelated errors; PP adjusts for this; DF-GLS gains power near a unit root; KPSS flips the null to guard against Type II errors.")
  #   bullet("Seasonal structure: large ACF at lag m (and possibly 2m) suggests seasonal AR or seasonal differencing; compare ACF/PACF patterns after differencing.")
  #   bullet("Breaks: a single level/trend shift can mimic a unit root. Always inspect plots and consider break tests before differencing too aggressively.")
  #   bullet("After differencing: fit candidate SARIMA(p,d,q)(P,D,Q)[m], then check residuals with Ljung–Box, normality, and (optionally) ARCH LM; compare models by AIC/AICc/BIC.")
  #   
  # })
  
  
  
  
  # output$teststationarited3St <- renderPrint({
  #   # ---------- Helpers ----------
  #   `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
  #   to_num_safe <- function(v, default = NA_real_) {
  #     out <- suppressWarnings(as.numeric(v))
  #     if (length(out) == 0 || all(is.na(out)) || !is.finite(out[1])) default else out[1]
  #   }
  #   to_int_safe <- function(v, default = 0L) {
  #     out <- suppressWarnings(as.integer(v))
  #     if (length(out) == 0 || is.na(out[1]) || !is.finite(out[1])) default else out[1]
  #   }
  #   safe_head_tail <- function(x, n = 5) {
  #     x <- as.numeric(x); x <- x[is.finite(x)]
  #     if (length(x) == 0) return(list(head = numeric(0), tail = numeric(0)))
  #     list(head = head(x, n), tail = tail(x, n))
  #   }
  #   tau_row_for <- function(type_in) switch(type_in, "none"="tau1", "drift"="tau2", "trend"="tau3", "tau3")
  #   
  #   # --- Pretty printers ---
  #   hr <- function() cat(paste0(strrep("=", 74), "\n"))
  #   subhr <- function() cat(paste0(strrep("-", 74), "\n"))
  #   h <- function(title) { hr(); cat(" ", title, "\n"); hr() }
  #   sec <- function(title) { subhr(); cat(" ", title, "\n"); subhr() }
  #   blt <- function(text) cat(sprintf("   • %s\n", text))
  #   kv  <- function(key, val) cat(sprintf("   - %-20s: %s\n", key, val))
  #   
  #   # Critical value extractor for urca objects (robust to naming)
  #   cval_pick_safe <- function(cval_obj, row, col) {
  #     if (is.null(cval_obj)) return(NA_real_)
  #     ans <- NA_real_
  #     if (is.matrix(cval_obj)) {
  #       if (!missing(row) && !missing(col) && row %in% rownames(cval_obj) && col %in% colnames(cval_obj)) {
  #         ans <- suppressWarnings(as.numeric(cval_obj[row, col]))
  #       } else if (!missing(col) && col %in% colnames(cval_obj)) {
  #         ans <- suppressWarnings(as.numeric(cval_obj[1, col]))
  #       } else {
  #         ans <- suppressWarnings(as.numeric(cval_obj[1, 1]))
  #       }
  #     } else {
  #       nm <- names(cval_obj)
  #       if (!is.null(nm) && col %in% nm) ans <- suppressWarnings(as.numeric(cval_obj[[col]]))
  #       if (!is.finite(ans) && length(cval_obj) >= 3) {
  #         if (identical(col, "10pct")) ans <- suppressWarnings(as.numeric(cval_obj[1]))
  #         if (identical(col, "5pct"))  ans <- suppressWarnings(as.numeric(cval_obj[2]))
  #         if (identical(col, "1pct"))  ans <- suppressWarnings(as.numeric(cval_obj[3]))
  #       }
  #       if (!is.finite(ans)) ans <- suppressWarnings(as.numeric(cval_obj[1]))
  #     }
  #     ans
  #   }
  #   
  #   # APA helpers
  #   fmt_p <- function(p) {
  #     if (!is.finite(p)) return("p = NA")
  #     if (p < .001) "p < .001" else paste0("p = ", sub("^0\\.", ".", sprintf("%.3f", p)))
  #   }
  #   apa_line_adf <- function(tau, p, k, N, type) {
  #     paste0("According to the ADF test (", type, ", k = ", k, ", N = ", N, "), the observed statistic was τ = ",
  #            sprintf("%.3f", tau), " with ", fmt_p(p), ".")
  #   }
  #   apa_line_kpss <- function(eta, p, mode) {
  #     paste0("According to the KPSS test (", mode, "), the observed statistic was η = ",
  #            sprintf("%.3f", eta), " with ", fmt_p(p), ".")
  #   }
  #   apa_line_lb <- function(stat, p, L) {
  #     paste0("Using a Ljung–Box test at lag ", L, ", the statistic was Q = ",
  #            ifelse(is.finite(stat), sprintf("%.3f", stat), "NA"), " with ", fmt_p(p), ".")
  #   }
  #   apa_line_pettitt <- function(U, p) {
  #     paste0("The Pettitt test yielded U = ", ifelse(is.finite(U), sprintf("%.3f", U), "NA"),
  #            " with ", fmt_p(p), ".")
  #   }
  #   
  #   # ---------- Inputs (from your UI) ----------
  #   alt_in  <- input$alternd2St %||% input$alternSt       # "stationary", "explosive", "regression"
  #   lag_in  <- input$LagOrderADFd2St %||% input$LagOrderADFSt
  #   a_in    <- input$alphaSt2
  #   type_in <- input$adfTypeSt2                            # "none", "drift", "trend"
  #   
  #   # Transformation flags (for provenance only; actual transform is in myData_Choice())
  #   d_in   <- input$d_n  %||% NA
  #   D_in   <- input$DS_n %||% NA
  #   log_in <- input$check_box %||% FALSE
  #   
  #   # ---------- Data ----------
  #   req(myData_Choice())
  #   x_raw <- myData_Choice()            # already transformed upstream in your app
  #   na_before <- sum(is.na(x_raw))
  #   x_class <- paste(class(x_raw), collapse = ", ")
  #   x_freq  <- if (inherits(x_raw, "ts")) tryCatch(stats::frequency(x_raw), error = function(e) NA_integer_) else NA_integer_
  #   x <- as.numeric(stats::na.omit(x_raw))
  #   N <- length(x)
  #   
  #   # ---------- Hyper-parameters ----------
  #   k <- to_int_safe(lag_in, default = 0L)
  #   if (!is.finite(k) || k < 0) k <- 0L
  #   alpha_raw <- as.character(a_in)
  #   alpha_val <- to_num_safe(alpha_raw, default = 0.05)
  #   alpha_col <- switch(alpha_raw, "0.01"="1pct", "0.05"="5pct", "0.1"="10pct", "0.10"="10pct",
  #                       "1pct"="1pct", "5pct"="5pct", "10pct"="10pct", "5pct")
  #   tau_row <- tau_row_for(type_in)
  #   seasonality_resolved <- isTRUE(to_int_safe(D_in, 0L) > 0L)
  #   
  #   # ---------- Sanity checks ----------
  #   if (N < 5) {
  #     h("CRITICAL ERROR")
  #     blt("Too few observations (N < 5). Provide more data.")
  #     return(invisible(NULL))
  #   }
  #   if (!is.finite(stats::sd(x)) || stats::sd(x) == 0) {
  #     h("CRITICAL ERROR")
  #     blt("Series is constant/invalid (sd = 0 or NA). Check transforms.")
  #     return(invisible(NULL))
  #   }
  #   
  #   # ---------- Packages ----------
  #   if (!requireNamespace("urca", quietly = TRUE) || !requireNamespace("tseries", quietly = TRUE)) {
  #     h("MISSING PACKAGES")
  #     blt("Please install.packages(c('urca','tseries')) to run ADF/KPSS.")
  #     return(invisible(NULL))
  #   }
  #   
  #   # ---------- ADF (urca + tseries p-value) ----------
  #   ur  <- tryCatch(urca::ur.df(x, type = type_in, lags = k), error = function(e) NULL)
  #   tau_obs  <- if (!is.null(ur)) suppressWarnings(as.numeric(ur@teststat[tau_row])) else NA_real_
  #   if (!is.finite(tau_obs) && !is.null(ur)) tau_obs <- suppressWarnings(as.numeric(ur@teststat[1]))
  #   tau_crit <- if (!is.null(ur)) cval_pick_safe(ur@cval, tau_row, alpha_col) else NA_real_
  #   
  #   adf_ts <- if (N > (k + 10)) tryCatch(tseries::adf.test(x, alternative = alt_in, k = k),
  #                                        error = function(e) NULL) else NULL
  #   adf_p  <- if (!is.null(adf_ts)) to_num_safe(adf_ts$p.value) else NA_real_
  #   
  #   # ---------- Ljung–Box on ADF residuals ----------
  #   lb_lag <- max(1L, min(10L, floor(N / 5)))
  #   lb     <- if (!is.null(ur)) stats::Box.test(ur@res, lag = lb_lag, type = "Ljung-Box") else NULL
  #   lb_stat<- if (!is.null(lb)) to_num_safe(lb$statistic) else NA_real_
  #   lb_p   <- if (!is.null(lb)) to_num_safe(lb$p.value) else NA_real_
  #   
  #   # ---------- KPSS (tseries + urca) ----------
  #   kpss_type <- if (type_in == "trend") "Trend" else "Level"
  #   kpss_ts   <- tryCatch(tseries::kpss.test(x, null = kpss_type), error = function(e) NULL)
  #   kpss_p    <- if (!is.null(kpss_ts)) to_num_safe(kpss_ts$p.value) else NA_real_
  #   
  #   kpss_uc   <- tryCatch(urca::ur.kpss(x, type = if (type_in == "trend") "tau" else "mu"),
  #                         error = function(e) NULL)
  #   eta_obs_uc  <- if (!is.null(kpss_uc)) to_num_safe(kpss_uc@teststat) else NA_real_
  #   eta_col     <- if (alpha_val <= 0.01) "1pct" else if (alpha_val <= 0.05) "5pct" else "10pct"
  #   eta_crit_uc <- if (!is.null(kpss_uc)) cval_pick_safe(kpss_uc@cval, 1, eta_col) else NA_real_
  #   
  #   # ---------- Pettitt (optional) ----------
  #   pett_U <- pett_p <- NA_real_
  #   if (requireNamespace("trend", quietly = TRUE)) {
  #     pet <- tryCatch(trend::pettitt.test(x), error = function(e) NULL)
  #     if (!is.null(pet)) {
  #       pett_U <- to_num_safe(pet$statistic); pett_p <- to_num_safe(pet$p.value)
  #     }
  #   }
  #   
  #   # ---------- Decisions ----------
  #   adf_ok_vals <- is.finite(tau_obs) && is.finite(tau_crit)
  #   adf_stationary <- adf_ok_vals && (tau_obs < tau_crit)
  #   
  #   kpss_by_p   <- is.finite(kpss_p)     && (kpss_p < alpha_val)
  #   kpss_by_eta <- is.finite(eta_obs_uc) && is.finite(eta_crit_uc) && (eta_obs_uc > eta_crit_uc)
  #   kpss_stationary <- !(kpss_by_p || kpss_by_eta)
  #   
  #   agreement <- (adf_stationary && kpss_stationary) || (!adf_stationary && !kpss_stationary)
  #   lb_pass   <- is.finite(lb_p) && (lb_p > alpha_val)
  #   
  #   # ---------- CHECKLIST ----------
  #   h("CHECKLIST")
  #   kv("ADF decision",        if (adf_ok_vals && adf_stationary) "[✓] STATIONARY" else "[X] NON-STATIONARY")
  #   kv("KPSS decision",       if (kpss_stationary) "[✓] STATIONARY" else "[X] NON-STATIONARY")
  #   kv("ADF vs KPSS",         if (agreement) "[✓] AGREEMENT" else "[?] CONFLICT")
  #   kv("Residual whiteness",  if (lb_pass) "[✓] PASS" else if (is.finite(lb_p)) "[X] FAIL" else "[?] NA")
  #   kv("Seasonal D>0 (UI)",   if (seasonality_resolved) "[✓] YES" else "[!] NO / UNKNOWN")
  #   kv("Pettitt available",   if (requireNamespace("trend", quietly = TRUE)) "[✓] YES" else "[!] NO")
  #   
  #   # ---------- SPECIFICATION & PROVENANCE ----------
  #   h("SPECIFICATION & PROVENANCE")
  #   kv("Model type (ADF)", toupper(type_in))
  #   kv("Lag k", k)
  #   kv("Alpha (α)", alpha_raw)
  #   kv("Series class", x_class)
  #   kv("Frequency", ifelse(is.finite(x_freq), x_freq, "NA"))
  #   kv("NA before omit", na_before)
  #   kv("Transforms (UI)", paste0("log=", ifelse(isTRUE(log_in), "ON", "OFF"),
  #                                ", d=", as.character(d_in),
  #                                ", D=", as.character(D_in)))
  #   ht <- safe_head_tail(x, 5)
  #   kv("First 5 used", paste(round(ht$head, 4), collapse = ", "))
  #   kv("Last 5 used",  paste(round(ht$tail, 4), collapse = ", "))
  #   
  #   # ---------- PHASE 1: ADF ----------
  #   h("PHASE 1 — AUGMENTED DICKEY–FULLER (ADF)")
  #   sec("Purpose")
  #   blt("Detect a nonseasonal unit root (random-walk behavior).")
  #   
  #   sec("Hypotheses")
  #   blt("H0: The series has a unit root (non-stationary).")
  #   blt("Ha: The series is stationary under the chosen deterministic terms.")
  #   
  #   sec("Decision rule")
  #   blt(sprintf("At α = %s, reject H0 if τ_obs < τ_crit.", alpha_raw))
  #   
  #   sec("Statistics")
  #   blt(sprintf("τ_obs = %.6f", tau_obs))
  #   blt(sprintf("τ_crit(%s) = %.6f", alpha_col, tau_crit))
  #   blt(sprintf("Reference p-value (tseries) = %s", if (is.finite(adf_p)) format.pval(adf_p, digits = 4) else "NA"))
  #   
  #   sec("Decision (APA)")
  #   cat(" ", apa_line_adf(tau_obs, adf_p, k, N, type_in), "\n", sep = "")
  #   if (adf_ok_vals && adf_stationary) {
  #     blt(sprintf("Decision: REJECT H0 (stationary) because τ_obs (%.6f) < τ_crit(%s) (%.6f).",
  #                 tau_obs, alpha_col, tau_crit))
  #   } else if (adf_ok_vals && !adf_stationary) {
  #     blt(sprintf("Decision: FAIL TO REJECT H0 (non-stationary) because τ_obs (%.6f) ≥ τ_crit(%s) (%.6f).",
  #                 tau_obs, alpha_col, tau_crit))
  #   } else {
  #     blt("Decision: INCONCLUSIVE — missing statistic/critical value.")
  #   }
  #   
  #   # ---------- PHASE 2: Residual whiteness (LB) ----------
  #   h("PHASE 2 — RESIDUAL WHITENESS (LJUNG–BOX)")
  #   sec("Hypotheses")
  #   blt("H0: ADF residuals are white noise (no autocorrelation).")
  #   blt("Ha: ADF residuals are autocorrelated.")
  #   
  #   sec("Statistics")
  #   blt(sprintf("Lag (L) = %d", lb_lag))
  #   blt(sprintf("Q = %s", ifelse(is.finite(lb_stat), sprintf("%.6f", lb_stat), "NA")))
  #   blt(sprintf("p-value = %s", if (is.finite(lb_p)) format.pval(lb_p, digits = 4) else "NA"))
  #   
  #   sec("Decision (APA)")
  #   cat(" ", apa_line_lb(lb_stat, lb_p, lb_lag), "\n", sep = "")
  #   if (is.finite(lb_p)) {
  #     if (lb_p > alpha_val) {
  #       blt("Decision: FAIL TO REJECT H0 → residuals are consistent with white noise (favorable).")
  #     } else {
  #       blt("Decision: REJECT H0 → residual autocorrelation remains (unfavorable).")
  #     }
  #   } else {
  #     blt("Decision: INCONCLUSIVE — Ljung–Box p-value is NA.")
  #   }
  #   
  #   # ---------- PHASE 3: KPSS ----------
  #   h("PHASE 3 — KPSS (STATIONARITY CONFIRMATION)")
  #   sec("Hypotheses")
  #   blt(sprintf("H0: Series is stationary around a %s.", if (kpss_type == "Trend") "trend" else "level"))
  #   blt("Ha: Series is non-stationary.")
  #   
  #   sec("Decision rule")
  #   blt(sprintf("At α = %s, reject H0 if p < α (tseries) OR if η_obs > η_crit (urca).", alpha_raw))
  #   
  #   sec("Statistics")
  #   blt(sprintf("η_obs (urca) = %.6f", eta_obs_uc))
  #   blt(sprintf("η_crit(%s) = %.6f", eta_col, eta_crit_uc))
  #   blt(sprintf("p-value (tseries) = %s", if (is.finite(kpss_p)) format.pval(kpss_p, digits = 4) else "NA"))
  #   
  #   sec("Decision (APA)")
  #   cat(" ", apa_line_kpss(eta_obs_uc, kpss_p, if (kpss_type == "Trend") "Trend" else "Level"), "\n", sep = "")
  #   if (kpss_stationary) {
  #     blt("Decision: FAIL TO REJECT H0 — stationarity supported.")
  #   } else {
  #     blt("Decision: REJECT H0 — non-stationarity indicated.")
  #   }
  #   
  #   # ---------- PHASE 4: Pettitt (optional) ----------
  #   h("PHASE 4 — STRUCTURAL BREAK (Pettitt)")
  #   if (requireNamespace("trend", quietly = TRUE)) {
  #     sec("Hypotheses")
  #     blt("H0: No change-point in location.")
  #     blt("Ha: A single change-point is present.")
  #     
  #     sec("Statistics")
  #     blt(sprintf("U = %s", ifelse(is.finite(pett_U), sprintf("%.6f", pett_U), "NA")))
  #     blt(sprintf("p-value = %s", if (is.finite(pett_p)) format.pval(pett_p, digits = 4) else "NA"))
  #     
  #     sec("Decision (APA)")
  #     cat(" ", apa_line_pettitt(pett_U, pett_p), "\n", sep = "")
  #     if (is.finite(pett_p)) {
  #       if (pett_p < alpha_val) {
  #         blt("Decision: REJECT H0 — evidence of a structural break; consider segment-wise testing.")
  #       } else {
  #         blt("Decision: FAIL TO REJECT H0 — no strong single break detected.")
  #       }
  #     } else {
  #       blt("Decision: INCONCLUSIVE — Pettitt p-value is NA.")
  #     }
  #   } else {
  #     blt("Package 'trend' not installed — Pettitt test skipped. (install.packages('trend'))")
  #   }
  #   
  #   # ---------- FINAL ACADEMIC ADVICE ----------
  #   h("FINAL ACADEMIC ADVICE")
  #   if (adf_stationary && kpss_stationary && lb_pass) {
  #     blt("Convergent evidence of stationarity (ADF & KPSS) with white residuals — proceed with ARMA/ARIMA on current series (d = 0).")
  #   } else if (!adf_stationary && !kpss_stationary) {
  #     blt("Both ADF and KPSS indicate non-stationarity — difference (d = 1); if seasonal, consider D = 1; then re-test.")
  #   } else {
  #     blt("ADF and KPSS conflict — try a different ADF type (none/drift/trend), resolve seasonality (D = 1), check for breaks, or stabilize variance (log/Box–Cox).")
  #   }
  #   
  #   # ---------- ACTIONABLE NEXT STEPS ----------
  #   h("ACTIONABLE NEXT STEPS")
  #   if (is.finite(lb_p) && lb_p <= alpha_val) {
  #     blt("Increase ADF lag k gradually until Ljung–Box p > α, or add differencing to reduce autocorrelation.")
  #   } else if (!is.finite(lb_p)) {
  #     blt("Ljung–Box p is NA — try a smaller k and ensure enough observations after transformations.")
  #   } else {
  #     blt("Residuals look white at current k — no immediate action needed.")
  #   }
  #   if (!adf_stationary || !kpss_stationary) {
  #     blt("For monthly multiplicative seasonality: try log + diff(1) + seasonal diff(12), then re-test ADF/KPSS.")
  #   } else {
  #     blt("Proceed to ARMA/ARIMA identification and compare by information criteria.")
  #   }
  #   if (requireNamespace("trend", quietly = TRUE) && is.finite(pett_p) && pett_p < alpha_val) {
  #     blt("Break detected by Pettitt — split at the change-point and re-run KPSS/ADF on segments.")
  #   } else {
  #     blt("If plots suggest shifts, consider a structural-break check (Perron, Bai–Perron) in addition to Pettitt.")
  #   }
  #   
  #   # ---------- SNAPSHOT ----------
  #   h("SNAPSHOT")
  #   kv("N", N)
  #   kv("type", type_in)
  #   kv("k", k)
  #   kv("α", alpha_raw)
  #   kv("freq", ifelse(is.finite(x_freq), x_freq, "NA"))
  #   kv("log", ifelse(isTRUE(log_in), "ON", "OFF"))
  #   kv("d", as.character(d_in))
  #   kv("D", as.character(D_in))
  #   kv("ADF", sprintf("τ=%.4f, τ_crit(%s)=%.4f, %s", tau_obs, alpha_col, tau_crit, fmt_p(adf_p)))
  #   kv("KPSS", sprintf("η=%.4f, η_crit(%s)=%.4f, %s", eta_obs_uc, eta_col, eta_crit_uc, fmt_p(kpss_p)))
  #   kv("LB", sprintf("%s (lag %d)", fmt_p(lb_p), lb_lag))
  #   if (requireNamespace("trend", quietly = TRUE)) kv("Pettitt", fmt_p(pett_p))
  #   hr()
  # })
  
  
  
  
  # output$teststationarited3St <- renderPrint({
  #   # ---------- Helpers ----------
  #   `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
  #   to_num_safe <- function(v, default = NA_real_) {
  #     out <- suppressWarnings(as.numeric(v))
  #     if (length(out) == 0 || all(is.na(out)) || !is.finite(out[1])) default else out[1]
  #   }
  #   to_int_safe <- function(v, default = 0L) {
  #     out <- suppressWarnings(as.integer(v))
  #     if (length(out) == 0 || is.na(out[1]) || !is.finite(out[1])) default else out[1]
  #   }
  #   safe_head_tail <- function(x, n = 5) {
  #     x <- as.numeric(x); x <- x[is.finite(x)]
  #     if (length(x) == 0) return(list(head = numeric(0), tail = numeric(0)))
  #     list(head = head(x, n), tail = tail(x, n))
  #   }
  #   tau_row_for <- function(type_in) switch(type_in, "none" = "tau1", "drift" = "tau2", "trend" = "tau3", "tau3")
  #   # Critical value extractor for urca objects (robust to naming)
  #   cval_pick_safe <- function(cval_obj, row, col) {
  #     if (is.null(cval_obj)) return(NA_real_)
  #     ans <- NA_real_
  #     if (is.matrix(cval_obj)) {
  #       if (!missing(row) && !missing(col) && row %in% rownames(cval_obj) && col %in% colnames(cval_obj)) {
  #         ans <- suppressWarnings(as.numeric(cval_obj[row, col]))
  #       } else if (!missing(col) && col %in% colnames(cval_obj)) {
  #         ans <- suppressWarnings(as.numeric(cval_obj[1, col]))
  #       } else {
  #         ans <- suppressWarnings(as.numeric(cval_obj[1, 1]))
  #       }
  #     } else {
  #       nm <- names(cval_obj)
  #       if (!is.null(nm) && col %in% nm) ans <- suppressWarnings(as.numeric(cval_obj[[col]]))
  #       if (!is.finite(ans) && length(cval_obj) >= 3) {
  #         if (identical(col, "10pct")) ans <- suppressWarnings(as.numeric(cval_obj[1]))
  #         if (identical(col, "5pct"))  ans <- suppressWarnings(as.numeric(cval_obj[2]))
  #         if (identical(col, "1pct"))  ans <- suppressWarnings(as.numeric(cval_obj[3]))
  #       }
  #       if (!is.finite(ans)) ans <- suppressWarnings(as.numeric(cval_obj[1]))
  #     }
  #     ans
  #   }
  #   # APA helpers
  #   fmt_p <- function(p) {
  #     if (!is.finite(p)) return("p = NA")
  #     if (p < .001) "p < .001" else paste0("p = ", sub("^0\\.", ".", sprintf("%.3f", p)))
  #   }
  #   apa_line_adf <- function(tau, p, k, N, type) {
  #     paste0("ADF (", type, ", k = ", k, ", N = ", N, "): τ = ", sprintf("%.3f", tau), ", ", fmt_p(p), ".")
  #   }
  #   apa_line_kpss <- function(eta, p, mode) {
  #     paste0("KPSS (", mode, "): η = ", sprintf("%.3f", eta), ", ", fmt_p(p), ".")
  #   }
  #   apa_line_lb <- function(stat, p, L) {
  #     paste0("Ljung–Box (lag = ", L, "): Q = ", ifelse(is.finite(stat), sprintf("%.3f", stat), "NA"), ", ", fmt_p(p), ".")
  #   }
  #   apa_line_pettitt <- function(U, p) {
  #     paste0("Pettitt: U = ", ifelse(is.finite(U), sprintf("%.3f", U), "NA"), ", ", fmt_p(p), ".")
  #   }
  #   
  #   # ---------- Inputs (from your UI) ----------
  #   alt_in  <- input$alternd2St %||% input$alternSt       # "stationary", "explosive", "regression"
  #   lag_in  <- input$LagOrderADFd2St %||% input$LagOrderADFSt
  #   a_in    <- input$alphaSt2
  #   type_in <- input$adfTypeSt2                            # "none", "drift", "trend"
  #   
  #   # Transformation flags (for provenance only; actual transform is in myData_Choice())
  #   d_in   <- input$d_n  %||% NA
  #   D_in   <- input$DS_n %||% NA
  #   log_in <- input$check_box %||% FALSE
  #   
  #   # ---------- Data ----------
  #   req(myData_Choice())
  #   x_raw <- myData_Choice()            # already transformed upstream in your app
  #   na_before <- sum(is.na(x_raw))
  #   x_class <- paste(class(x_raw), collapse = ", ")
  #   x_freq  <- if (inherits(x_raw, "ts")) tryCatch(stats::frequency(x_raw), error = function(e) NA_integer_) else NA_integer_
  #   x <- as.numeric(stats::na.omit(x_raw))
  #   N <- length(x)
  #   
  #   # ---------- Hyper-parameters ----------
  #   k <- to_int_safe(lag_in, default = 0L)
  #   if (!is.finite(k) || k < 0) k <- 0L
  #   alpha_raw <- as.character(a_in)
  #   alpha_val <- to_num_safe(alpha_raw, default = 0.05)
  #   alpha_col <- switch(alpha_raw, "0.01" = "1pct", "0.05" = "5pct", "0.1" = "10pct", "0.10" = "10pct",
  #                       "1pct" = "1pct", "5pct" = "5pct", "10pct" = "10pct", "5pct")
  #   tau_row <- tau_row_for(type_in)
  #   seasonality_resolved <- isTRUE(to_int_safe(D_in, 0L) > 0L)
  #   
  #   # ---------- Sanity checks ----------
  #   if (N < 5) {
  #     cat("==========================================================================\n")
  #     cat("CRITICAL ERROR: Too few observations (N < 5). Provide more data.\n")
  #     cat("==========================================================================\n"); return(invisible(NULL))
  #   }
  #   if (!is.finite(stats::sd(x)) || stats::sd(x) == 0) {
  #     cat("==========================================================================\n")
  #     cat("CRITICAL ERROR: Series is constant/invalid (sd = 0 or NA). Check transforms.\n")
  #     cat("==========================================================================\n"); return(invisible(NULL))
  #   }
  #   
  #   # ---------- Packages ----------
  #   if (!requireNamespace("urca", quietly = TRUE) || !requireNamespace("tseries", quietly = TRUE)) {
  #     cat("==========================================================================\n")
  #     cat("ERROR: Please install.packages(c('urca','tseries')) to run ADF/KPSS.\n")
  #     cat("==========================================================================\n"); return(invisible(NULL))
  #   }
  #   
  #   # ---------- ADF (urca + tseries p-value) ----------
  #   ur  <- tryCatch(urca::ur.df(x, type = type_in, lags = k), error = function(e) NULL)
  #   tau_obs  <- if (!is.null(ur)) suppressWarnings(as.numeric(ur@teststat[tau_row])) else NA_real_
  #   if (!is.finite(tau_obs) && !is.null(ur)) tau_obs <- suppressWarnings(as.numeric(ur@teststat[1]))
  #   tau_crit <- if (!is.null(ur)) cval_pick_safe(ur@cval, tau_row, alpha_col) else NA_real_
  #   
  #   adf_ts <- if (N > (k + 10)) tryCatch(tseries::adf.test(x, alternative = alt_in, k = k),
  #                                        error = function(e) NULL) else NULL
  #   adf_p  <- if (!is.null(adf_ts)) to_num_safe(adf_ts$p.value) else NA_real_
  #   
  #   # ---------- Ljung–Box on ADF residuals ----------
  #   lb_lag <- max(1L, min(10L, floor(N / 5)))
  #   lb     <- if (!is.null(ur)) stats::Box.test(ur@res, lag = lb_lag, type = "Ljung-Box") else NULL
  #   lb_stat<- if (!is.null(lb)) to_num_safe(lb$statistic) else NA_real_
  #   lb_p   <- if (!is.null(lb)) to_num_safe(lb$p.value) else NA_real_
  #   
  #   # ---------- KPSS (tseries + urca) ----------
  #   kpss_type <- if (type_in == "trend") "Trend" else "Level"
  #   kpss_ts   <- tryCatch(tseries::kpss.test(x, null = kpss_type), error = function(e) NULL)
  #   kpss_p    <- if (!is.null(kpss_ts)) to_num_safe(kpss_ts$p.value) else NA_real_
  #   
  #   kpss_uc   <- tryCatch(urca::ur.kpss(x, type = if (type_in == "trend") "tau" else "mu"),
  #                         error = function(e) NULL)
  #   eta_obs_uc  <- if (!is.null(kpss_uc)) to_num_safe(kpss_uc@teststat) else NA_real_
  #   eta_col     <- if (alpha_val <= 0.01) "1pct" else if (alpha_val <= 0.05) "5pct" else "10pct"
  #   eta_crit_uc <- if (!is.null(kpss_uc)) cval_pick_safe(kpss_uc@cval, 1, eta_col) else NA_real_
  #   
  #   # ---------- Pettitt (optional) ----------
  #   pett_U <- pett_p <- NA_real_
  #   if (requireNamespace("trend", quietly = TRUE)) {
  #     pet <- tryCatch(trend::pettitt.test(x), error = function(e) NULL)
  #     if (!is.null(pet)) {
  #       pett_U <- to_num_safe(pet$statistic); pett_p <- to_num_safe(pet$p.value)
  #     }
  #   }
  #   
  #   # ---------- Decisions ----------
  #   # ADF decision (left-tail): reject H0 if tau_obs < tau_crit
  #   adf_ok_vals <- is.finite(tau_obs) && is.finite(tau_crit)
  #   adf_stationary <- adf_ok_vals && (tau_obs < tau_crit)
  #   
  #   # KPSS decision (H0 = stationary): reject H0 if p < alpha OR eta_obs > eta_crit
  #   kpss_by_p   <- is.finite(kpss_p)    && (kpss_p < alpha_val)
  #   kpss_by_eta <- is.finite(eta_obs_uc) && is.finite(eta_crit_uc) && (eta_obs_uc > eta_crit_uc)
  #   kpss_stationary <- !(kpss_by_p || kpss_by_eta)
  #   
  #   agreement <- (adf_stationary && kpss_stationary) || (!adf_stationary && !kpss_stationary)
  #   lb_pass   <- is.finite(lb_p) && (lb_p > alpha_val)
  #   
  #   # ---------- PRINT: CHECKLIST first ----------
  #   cat("==========================================================================\n")
  #   cat(" CHECKLIST\n")
  #   cat("==========================================================================\n")
  #   cat(sprintf(" [ ] ADF decision (reject unit root => stationary) : %s\n",
  #               if (adf_ok_vals && adf_stationary) "[✓] STATIONARY" else "[X] NON-STATIONARY"))
  #   cat(sprintf(" [ ] KPSS decision (fail reject => stationary)     : %s\n",
  #               if (kpss_stationary) "[✓] STATIONARY" else "[X] NON-STATIONARY"))
  #   cat(sprintf(" [ ] ADF vs KPSS agreement                         : %s\n",
  #               if (agreement) "[✓] AGREEMENT" else "[?] CONFLICT"))
  #   cat(sprintf(" [ ] Residual whiteness (LB (Ljung–Box) p>α)       : %s\n",
  #               if (lb_pass) "[✓] PASS" else if (is.finite(lb_p)) "[X] FAIL" else "[?] NA"))
  #   cat(sprintf(" [ ] Seasonal differencing indicated (UI D>0)      : %s\n",
  #               if (seasonality_resolved) "[✓] YES" else "[!] NO / UNKNOWN"))
  #   cat(sprintf(" [ ] Pettitt break check available                 : %s\n\n",
  #               if (requireNamespace("trend", quietly = TRUE)) "[✓] YES" else "[!] NO"))
  #   
  #   # ---------- SPECIFICATION & PROVENANCE ----------
  #   cat("SPECIFICATION & PROVENANCE\n")
  #   cat("--------------------------------------------------------------------------\n")
  #   cat(sprintf(" Model type (ADF): %s | Lag k: %d | α: %s\n", toupper(type_in), k, alpha_raw))
  #   cat(sprintf(" Series class: %s | frequency: %s | NA before omit: %d\n",
  #               x_class, ifelse(is.finite(x_freq), x_freq, "NA"), na_before))
  #   cat(sprintf(" Transform flags (UI): log=%s, d=%s, D=%s\n",
  #               ifelse(isTRUE(log_in), "ON", "OFF"),
  #               as.character(d_in), as.character(D_in)))
  #   ht <- safe_head_tail(x, 5)
  #   cat(sprintf(" First 5 values used: %s\n", paste(round(ht$head, 4), collapse = ", ")))
  #   cat(sprintf(" Last  5 values used: %s\n\n", paste(round(ht$tail, 4), collapse = ", ")))
  #   
  #   # ---------- PHASE 1: ADF ----------
  #   cat("==========================================================================\n")
  #   cat(" PHASE 1 — AUGMENTED DICKEY–FULLER (ADF) TEST\n")
  #   cat("==========================================================================\n")
  #   cat(" Use case: detect a nonseasonal unit root (random walk behavior).\n")
  #   cat(" H0 (unit root): the series is non-stationary.\n")
  #   cat(" Ha (stationary): the series is stationary (mean-reverting) under the chosen deterministic terms.\n")
  #   cat(sprintf(" Statistic: τ (tau). Decision rule at α=%s: Reject H0 if τ_obs < τ_crit.\n", alpha_raw))
  #   cat(sprintf(" Observed: τ_obs = %.6f | Critical: τ_crit(%s) = %.6f | p(ref, tseries) = %s\n",
  #               tau_obs, alpha_col, tau_crit, if (is.finite(adf_p)) format.pval(adf_p, digits = 4) else "NA"))
  #   
  #   if (adf_ok_vals && adf_stationary) {
  #     cat(" Decision: REJECT H0 (stationary).\n")
  #     cat(sprintf(" Justification: τ_obs (%.6f) is more negative than τ_crit(%s) (%.6f), providing evidence that the series is stationary under the %s specification.\n",
  #                 tau_obs, alpha_col, tau_crit, type_in))
  #   } else if (adf_ok_vals && !adf_stationary) {
  #     cat(" Decision: FAIL TO REJECT H0 (non-stationary).\n")
  #     cat(sprintf(" Justification: τ_obs (%.6f) is not more negative than τ_crit(%s) (%.6f), so we lack evidence to reject the unit-root null under the %s specification.\n",
  #                 tau_obs, alpha_col, tau_crit, type_in))
  #   } else {
  #     cat(" Decision: INCONCLUSIVE (missing statistic/critical value).\n")
  #     cat(" Justification: τ or τ_crit is NA/Inf; the regression could not produce valid inference (often due to k too large vs N, or constant series).\n")
  #   }
  #   cat(" APA: ", apa_line_adf(tau_obs, adf_p, k, N, type_in), "\n\n", sep = "")
  #   
  #   # ---------- PHASE 2: Residual diagnostics ----------
  #   cat("==========================================================================\n")
  #   cat(" PHASE 2 — RESIDUAL WHITENESS (LJUNG–BOX)\n")
  #   cat("==========================================================================\n")
  #   cat(" Use case: verify that the ADF regression residuals are approximately white noise.\n")
  #   cat(" H0: residuals are white noise (no autocorrelation).  Ha: residuals are autocorrelated.\n")
  #   cat(sprintf(" Observed: Q = %s | p-value = %s | lag = %d\n",
  #               ifelse(is.finite(lb_stat), sprintf("%.6f", lb_stat), "NA"),
  #               ifelse(is.finite(lb_p), format.pval(lb_p, digits = 4), "NA"),
  #               lb_lag))
  #   if (is.finite(lb_p)) {
  #     if (lb_p > alpha_val) {
  #       cat(" Decision: FAIL TO REJECT H0 → residuals are consistent with white noise (favorable).\n")
  #       cat(" Justification: p-value exceeds α, indicating insufficient evidence of residual autocorrelation.\n")
  #     } else {
  #       cat(" Decision: REJECT H0 → residual autocorrelation remains (unfavorable).\n")
  #       cat(" Justification: p-value below α suggests remaining autocorrelation; increase k or difference appropriately.\n")
  #     }
  #   } else {
  #     cat(" Decision: INCONCLUSIVE (LB p-value NA). Consider smaller k or more data.\n")
  #   }
  #   cat(" APA: ", apa_line_lb(lb_stat, lb_p, lb_lag), "\n\n", sep = "")
  #   
  #   # ---------- PHASE 3: KPSS ----------
  #   cat("==========================================================================\n")
  #   cat(" PHASE 3 — KPSS TEST (STATIONARITY CONFIRMATION)\n")
  #   cat("==========================================================================\n")
  #   cat(" Use case: confirm stationarity (complements ADF with opposite null).\n")
  #   cat(sprintf(" H0: series is stationary around a %s.  Ha: series is non-stationary.\n", if (kpss_type == "Trend") "trend" else "level"))
  #   cat(sprintf(" Statistics: η (eta). Decision at α=%s: reject H0 if p < α (tseries) OR η_obs > η_crit (urca).\n", alpha_raw))
  #   cat(sprintf(" Observed: η_obs (urca) = %.6f | η_crit(%s) = %.6f | p(tseries) = %s\n",
  #               eta_obs_uc, eta_col, eta_crit_uc,
  #               if (is.finite(kpss_p)) format.pval(kpss_p, digits = 4) else "NA"))
  #   
  #   if (kpss_stationary) {
  #     cat(" Decision: FAIL TO REJECT H0 (stationarity supported).\n")
  #     cat(" Justification: KPSS does not provide sufficient evidence against stationarity under the tested specification.\n")
  #   } else {
  #     cat(" Decision: REJECT H0 (non-stationary).\n")
  #     if (is.finite(kpss_p) && kpss_p < alpha_val) {
  #       cat(" Justification: The one-sided p-value is below α, indicating stationarity is unlikely.\n")
  #     } else if (is.finite(eta_obs_uc) && is.finite(eta_crit_uc) && eta_obs_uc > eta_crit_uc) {
  #       cat(" Justification: The observed η exceeds its critical value, indicating non-stationarity.\n")
  #     } else {
  #       cat(" Justification: KPSS indicates departure from stationarity via at least one decision criterion.\n")
  #     }
  #   }
  #   cat(" APA: ", apa_line_kpss(eta_obs_uc, kpss_p, if (kpss_type == "Trend") "Trend" else "Level"), "\n\n", sep = "")
  #   
  #   # ---------- PHASE 4: Pettitt (optional) ----------
  #   cat("==========================================================================\n")
  #   cat(" PHASE 4 — STRUCTURAL BREAK (Pettitt)\n")
  #   cat("==========================================================================\n")
  #   if (requireNamespace("trend", quietly = TRUE)) {
  #     cat(" Use case: detect a single change-point that can contaminate ADF/KPSS decisions.\n")
  #     cat(sprintf(" Observed: U = %s | p-value = %s\n",
  #                 ifelse(is.finite(pett_U), sprintf("%.6f", pett_U), "NA"),
  #                 ifelse(is.finite(pett_p), format.pval(pett_p, digits = 4), "NA")))
  #     if (is.finite(pett_p)) {
  #       if (pett_p < alpha_val) {
  #         cat(" Decision: REJECT H0 → evidence of a structural break; consider segment-wise testing.\n")
  #       } else {
  #         cat(" Decision: FAIL TO REJECT H0 → no strong single break detected.\n")
  #       }
  #     } else {
  #       cat(" Decision: INCONCLUSIVE (Pettitt p-value NA).\n")
  #     }
  #     cat(" APA: ", apa_line_pettitt(pett_U, pett_p), "\n\n", sep = "")
  #   } else {
  #     cat(" Package 'trend' not installed — Pettitt test skipped. (install.packages('trend'))\n\n")
  #   }
  #   
  #   # ---------- FINAL ACADEMIC ADVICE ----------
  #   cat("==========================================================================\n")
  #   cat(" FINAL ACADEMIC ADVICE\n")
  #   cat("==========================================================================\n")
  #   if (adf_stationary && kpss_stationary && lb_pass) {
  #     cat(" The series exhibits convergent evidence of stationarity (ADF, KPSS) with white residuals.\n")
  #     cat(" Proceed with ARMA/ARIMA identification on the current (transformed) series (d=0), and validate with residual diagnostics.\n\n")
  #   } else if (!adf_stationary && !kpss_stationary) {
  #     cat(" Both ADF and KPSS point to non-stationarity. Apply differencing (d=1) and, if seasonal, D=1.\n")
  #     cat(" After transforming (e.g., log + diff(1) + seasonal diff), re-run diagnostics before modeling.\n\n")
  #   } else {
  #     cat(" ADF and KPSS are conflicting. This is common near a unit root or under mis-specified deterministic terms.\n")
  #     cat(" Consider: (i) trying a different ADF type (none/drift/trend), (ii) resolving seasonality (D=1),\n")
  #     cat(" (iii) checking for structural breaks (Pettitt) and testing segments, and (iv) variance stabilization (log/Box–Cox).\n\n")
  #   }
  #   
  #   # ---------- ACTIONABLE NEXT STEPS ----------
  #   cat("==========================================================================\n")
  #   cat(" ACTIONABLE NEXT STEPS\n")
  #   cat("==========================================================================\n")
  #   if (is.finite(lb_p) && lb_p <= alpha_val) {
  #     cat(" [1] Increase ADF lag k gradually until Ljung–Box p > α, or apply additional differencing to reduce autocorrelation.\n")
  #   } else if (!is.finite(lb_p)) {
  #     cat(" [1] Ljung–Box p is NA: try a smaller k and ensure enough observations after transformations.\n")
  #   } else {
  #     cat(" [1] Residuals are acceptably white at current k; no immediate action needed.\n")
  #   }
  #   if (!adf_stationary || !kpss_stationary) {
  #     cat(" [2] For monthly data with multiplicative seasonality: try log + diff(1) + seasonal diff(12); then re-test ADF/KPSS.\n")
  #   } else {
  #     cat(" [2] Proceed to model identification (ARMA/ARIMA) and information-criteria comparison.\n")
  #   }
  #   if (requireNamespace("trend", quietly = TRUE) && is.finite(pett_p) && pett_p < alpha_val) {
  #     cat(" [3] Break detected by Pettitt: split at the estimated change-point and re-run KPSS/ADF on each segment.\n")
  #   } else {
  #     cat(" [3] Consider a structural-break check (Pettitt) if plots suggest level shifts.\n")
  #   }
  #   cat("\n--------------------------------------------------------------------------\n")
  #   cat(sprintf(" SNAPSHOT — N=%d | type=%s | k=%d | α=%s | freq=%s | log=%s | d=%s | D=%s\n",
  #               N, type_in, k, alpha_raw, ifelse(is.finite(x_freq), x_freq, "NA"),
  #               ifelse(isTRUE(log_in), "ON", "OFF"), as.character(d_in), as.character(D_in)))
  #   cat(sprintf(" ADF: τ=%.4f, τ_crit(%s)=%.4f, %s\n", tau_obs, alpha_col, tau_crit, fmt_p(adf_p)))
  #   cat(sprintf(" KPSS: η=%.4f, η_crit(%s)=%.4f, %s\n", eta_obs_uc, eta_col, eta_crit_uc, fmt_p(kpss_p)))
  #   cat(sprintf(" LB: %s (lag %d)\n", fmt_p(lb_p), lb_lag))
  #   if (requireNamespace("trend", quietly = TRUE)) {
  #     cat(sprintf(" Pettitt: %s\n", fmt_p(pett_p)))
  #   }
  #   cat("==========================================================================\n")
  # })
  
  
  
  
  
  # --- Replace existing output$teststationarited3St with this ---
  # output$teststationarited3St <- renderPrint({
  #   `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
  #   to_num_safe <- function(v, default = NA_real_) { out <- suppressWarnings(as.numeric(v)); if (length(out)==0 || all(is.na(out)) || !is.finite(out[1])) default else out[1] }
  #   to_int_safe <- function(v, default = 0L) { out <- suppressWarnings(as.integer(v)); if (length(out)==0 || is.na(out[1]) || !is.finite(out[1])) default else out[1] }
  #   safe_head_tail <- function(x, n = 5) { x <- as.numeric(x); x <- x[is.finite(x)]; if (!length(x)) return(list(head=numeric(0),tail=numeric(0))); list(head=head(x,n), tail=tail(x,n)) }
  #   tau_row_for <- function(type_in) switch(type_in, "none"="tau1","drift"="tau2","trend"="tau3","tau3")
  #   
  #   # ---- Inputs (using your ids) ----
  #   alt_in  <- input$alternd2St %||% input$alternSt
  #   lag_in  <- input$LagOrderADFd2St %||% input$LagOrderADFSt
  #   a_in    <- input$alphaSt2
  #   type_in <- input$adfTypeSt2
  #   
  #   # (Only for reporting provenance)
  #   d_in   <- input$d_n  %||% NA
  #   D_in   <- input$DS_n %||% NA
  #   log_in <- input$check_box %||% FALSE
  #   
  #   # ---- Data ----
  #   req(myData_Choice())
  #   x_raw <- myData_Choice()
  #   na_before <- sum(is.na(x_raw))
  #   x_class <- paste(class(x_raw), collapse = ", ")
  #   x_freq  <- if (inherits(x_raw, "ts")) tryCatch(stats::frequency(x_raw), error = function(e) NA_integer_) else NA_integer_
  #   x <- as.numeric(stats::na.omit(x_raw))
  #   N <- length(x)
  #   
  #   # ---- Hyper-params ----
  #   k <- to_int_safe(lag_in, default = 0L)
  #   alpha_raw <- as.character(a_in)
  #   alpha_val <- to_num_safe(alpha_raw, default = 0.05)
  #   alpha_col <- switch(alpha_raw,
  #                       "0.01"="1pct","0.05"="5pct","0.1"="10pct","0.10"="10pct",
  #                       "1pct"="1pct","5pct"="5pct","10pct"="10pct","5pct")
  #   tau_row <- tau_row_for(type_in)
  #   seasonality_resolved <- to_int_safe(D_in,0L) > 0L
  #   
  #   # ---- Quick sanity ----
  #   cat("==========================================================================\n")
  #   cat("                        ADF UNIT ROOT DIAGNOSTIC                          \n")
  #   cat("==========================================================================\n")
  #   cat(sprintf(" MODEL TYPE : %-10s | SAMPLE SIZE (N) : %d\n", toupper(type_in), N))
  #   cat(sprintf(" LAG ORDER  : %-10d | SIGNIFICANCE (α) : %s\n", k, alpha_raw))
  #   cat(sprintf(" MOMENTS    : Mean: %.4f | Std.Dev: %.4f\n", mean(x), stats::sd(x)))
  #   cat("--------------------------------------------------------------------------\n")
  #   cat(" TRANSFORMATION PROVENANCE (common to ALL tests):\n")
  #   cat(sprintf(" [ ] Source class: %s | Frequency: %s | NA before omit: %d\n",
  #               x_class, ifelse(is.finite(x_freq), x_freq, NA), na_before))
  #   cat(sprintf(" [ ] UI transforms: log=%s | d=%s | D=%s\n",
  #               ifelse(isTRUE(log_in),"ON","OFF"), as.character(d_in), as.character(D_in)))
  #   ht <- safe_head_tail(x,5)
  #   cat(sprintf(" [ ] First 5 used: %s\n", paste(round(ht$head,4), collapse=", ")))
  #   cat(sprintf(" [ ] Last  5 used: %s\n\n", paste(round(ht$tail,4), collapse=", ")))
  #   
  #   if (N < 5) { cat(" [!] Too few observations (N < 5). Abort.\n==========================================================================\n"); return(invisible(NULL)) }
  #   if (!is.finite(stats::sd(x)) || stats::sd(x) == 0) { cat(" [!] Series constant/invalid (sd=0 or NA). Abort.\n==========================================================================\n"); return(invisible(NULL)) }
  #   if (N <= (k + 10)) cat(" [!] WARNING: high lag vs sample; ADF may be weak. Consider smaller k.\n\n")
  #   
  #   # ---- Compute tests ----
  #   suppressMessages({
  #     need_urca   <- requireNamespace("urca", quietly = TRUE)
  #     need_tseries<- requireNamespace("tseries", quietly = TRUE)
  #   })
  #   if (!need_urca || !need_tseries) {
  #     cat(" [!] Please install.packages(c('urca','tseries')) to run ADF/KPSS.\n==========================================================================\n")
  #     return(invisible(NULL))
  #   }
  #   
  #   # Main ADF via urca (for tau + criticals)
  #   ur  <- tryCatch(urca::ur.df(x, type = type_in, lags = k), error = function(e) NULL)
  #   tau_obs  <- if (!is.null(ur)) suppressWarnings(as.numeric(ur@teststat[tau_row])) else NA_real_
  #   if (!is.finite(tau_obs) && !is.null(ur)) tau_obs <- suppressWarnings(as.numeric(ur@teststat[1]))
  #   tau_crit <- if (!is.null(ur)) suppressWarnings(as.numeric(ur@cval[tau_row, alpha_col])) else NA_real_
  #   
  #   # p-value ref from tseries (only if df ok)
  #   adf_ts <- if (N > (k + 10)) tryCatch(tseries::adf.test(x, alternative = alt_in, k = k),
  #                                        error = function(e) NULL) else NULL
  #   
  #   # Residual whiteness from ur.df regression
  #   lb_lag <- max(1L, min(10L, floor(N/5)))
  #   lb     <- if (!is.null(ur)) Box.test(ur@res, lag = lb_lag, type = "Ljung-Box") else NULL
  #   
  #   # KPSS (both flavors: tseries p and urca eta)
  #   kpss_type <- if (type_in == "trend") "Trend" else "Level"
  #   kpss_ts   <- tryCatch(tseries::kpss.test(x, null = kpss_type), error = function(e) NULL)
  #   kpss_uc   <- tryCatch(urca::ur.kpss(x, type = if (type_in=="trend") "tau" else "mu"),
  #                         error = function(e) NULL)
  #   eta_obs_uc  <- if (!is.null(kpss_uc)) suppressWarnings(as.numeric(kpss_uc@teststat)) else NA_real_
  #   eta_col     <- if (alpha_val <= .01) "1pct" else if (alpha_val <= .05) "5pct" else "10pct"
  #   eta_crit_uc <- if (!is.null(kpss_uc)) {
  #     cv <- kpss_uc@cval
  #     if (is.matrix(cv) && eta_col %in% colnames(cv)) as.numeric(cv[1, eta_col]) else
  #       if (!is.matrix(cv) && !is.null(names(cv)) && eta_col %in% names(cv)) as.numeric(cv[[eta_col]]) else NA_real_
  #   } else NA_real_
  #   
  #   # ---- PHASE 1: ADF ----
  #   cat("==========================================================================\n")
  #   cat("PHASE 1: ADF (urca::ur.df + tseries::adf.test)\n")
  #   cat("==========================================================================\n")
  #   cat(sprintf(" Tau-Obs = %.4f | Tau-Crit(%s) = %.4f | p(ref) = %s\n",
  #               tau_obs, alpha_col, tau_crit,
  #               if (!is.null(adf_ts)) format.pval(adf_ts$p.value, digits = 4) else "NA"))
  #   adf_stationary <- is.finite(tau_obs) && is.finite(tau_crit) && (tau_obs < tau_crit)
  #   cat(" Decision: ", if (isTRUE(adf_stationary)) "REJECT H0 (stationary)\n" else "Fail to reject H0 (non-stationary)\n")
  #   
  #   # ---- PHASE 2: Residual whiteness ----
  #   cat("\n==========================================================================\n")
  #   cat("PHASE 2: Ljung–Box on ADF regression residuals\n")
  #   cat("==========================================================================\n")
  #   lb_p <- if (!is.null(lb)) suppressWarnings(as.numeric(lb$p.value)) else NA_real_
  #   lb_st<- if (!is.null(lb)) suppressWarnings(as.numeric(lb$statistic)) else NA_real_
  #   cat(sprintf(" LB lag=%d | LB stat=%.4f | LB p=%.4f  → %s\n",
  #               lb_lag, ifelse(is.finite(lb_st), lb_st, NA), ifelse(is.finite(lb_p), lb_p, NA),
  #               if (is.finite(lb_p) && lb_p > alpha_val) "white noise OK" else "autocorr remains"))
  #   
  #   # ---- PHASE 3: KPSS ----
  #   cat("\n==========================================================================\n")
  #   cat("PHASE 3: KPSS (tseries p + urca eta)\n")
  #   cat("==========================================================================\n")
  #   kpss_p <- if (!is.null(kpss_ts)) suppressWarnings(as.numeric(kpss_ts$p.value)) else NA_real_
  #   cat(sprintf(" KPSS p (tseries) = %s | Eta(obs) (urca) = %.4f | Eta(crit %s) = %.4f\n",
  #               if (is.finite(kpss_p)) format.pval(kpss_p, digits = 4) else "NA",
  #               eta_obs_uc, eta_col, eta_crit_uc))
  #   kpss_stationary <- !( (is.finite(kpss_p) && kpss_p < alpha_val) ||
  #                           (is.finite(eta_obs_uc) && is.finite(eta_crit_uc) && eta_obs_uc > eta_crit_uc) )
  #   cat(" Decision: ", if (kpss_stationary) "Fail to reject H0 (stationary)\n" else "REJECT H0 (non-stationary)\n")
  #   
  #   # ---- FINAL VERDICT ----
  #   cat("\n==========================================================================\n")
  #   cat("PHASE 4: FINAL VERDICT & ACTIONS\n")
  #   cat("==========================================================================\n")
  #   agree <- (adf_stationary && kpss_stationary) || (!adf_stationary && !kpss_stationary)
  #   if (isTRUE(lb_p <= alpha_val)) cat(" [!] Residual autocorrelation → consider increasing k or differencing.\n")
  #   if (adf_stationary && kpss_stationary) {
  #     cat(" [✓] Strong evidence of STATIONARITY. Proceed with ARMA/ARIMA at d=0.\n")
  #   } else if (!adf_stationary && !kpss_stationary) {
  #     cat(" [X] Clear UNIT ROOT. Apply differencing (d=1) and/or seasonal differencing (D=1 if seasonal).\n")
  #   } else {
  #     cat(" [?] ADF vs KPSS conflict. Check trend spec (none/drift/trend), seasonality (D), and breaks.\n")
  #   }
  #   
  #   # ---- Snapshot ----
  #   cat("\n--------------------------------------------------------------------------\n")
  #   cat("SNAPSHOT\n")
  #   cat("--------------------------------------------------------------------------\n")
  #   cat(sprintf(" N=%d | type=%s | k=%d | α=%s | freq=%s | log=%s | d=%s | D=%s\n",
  #               N, type_in, k, alpha_raw, ifelse(is.finite(x_freq), x_freq, "NA"),
  #               ifelse(isTRUE(log_in),"ON","OFF"), as.character(d_in), as.character(D_in)))
  #   cat(sprintf(" ADF: tau=%.4f, crit(%s)=%.4f, p(ref)=%s\n",
  #               tau_obs, alpha_col, tau_crit,
  #               if (!is.null(adf_ts)) format.pval(adf_ts$p.value, digits = 4) else "NA"))
  #   cat(sprintf(" LB:  p=%.4f (lag %d)\n", ifelse(is.finite(lb_p), lb_p, NA), lb_lag))
  #   cat(sprintf(" KPSS: p=%.4f, eta(obs)=%.4f vs eta(crit %s)=%.4f\n",
  #               ifelse(is.finite(kpss_p), kpss_p, NA), eta_obs_uc, eta_col, eta_crit_uc))
  #   cat("==========================================================================\n")
  # })
  
  
  
  # output$teststationarited3St <- renderPrint({
  #   
  #   # ============================================================================
  #   # 0) SMALL HELPERS (safe input fallback + safe numeric)
  #   # ============================================================================
  #   `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
  #   to_num_safe <- function(v, default = NA_real_) {
  #     out <- suppressWarnings(as.numeric(v))
  #     if (length(out) == 0 || all(is.na(out)) || !is.finite(out[1])) default else out[1]
  #   }
  #   
  #   # ============================================================================
  #   # 1) INPUT COLLECTION (supports either naming convention in your UI)
  #   #    - This prevents "silent req stops" if your UI ids differ.
  #   # ============================================================================
  #   alt_in  <- input$alternd2St %||% input$alternSt
  #   lag_in  <- input$LagOrderADFd2St %||% input$LagOrderADFSt
  #   a_in    <- input$alphaSt2
  #   type_in <- input$adfTypeSt2
  #   
  #   # Required data
  #   req(myData_Choice())
  #   
  #   # Required inputs (validated safely)
  #   if (is.null(alt_in) || is.null(lag_in) || is.null(a_in) || is.null(type_in)) {
  #     cat("==========================================================================\n")
  #     cat("                STATE-OF-THE-ART ADF UNIT ROOT DIAGNOSTIC                 \n")
  #     cat("==========================================================================\n")
  #     cat(" [!] INPUT ERROR: One or more inputs are NULL.\n")
  #     cat("     This usually means a UI/server ID mismatch.\n")
  #     cat("     Needed inputs (either naming is OK):\n")
  #     cat("       - alternative : input$alternd2St OR input$alternSt\n")
  #     cat("       - lag         : input$LagOrderADFd2St OR input$LagOrderADFSt\n")
  #     cat("       - alpha       : input$alphaSt2\n")
  #     cat("       - adf type    : input$adfTypeSt2\n")
  #     cat("==========================================================================\n")
  #     return(invisible(NULL))
  #   }
  #   
  #   # ============================================================================
  #   # 2) DATA PREP
  #   # ============================================================================
  #   x <- as.numeric(stats::na.omit(myData_Choice()))
  #   valid_N <- length(x)
  #   
  #   # lag (k) must be integer >= 0
  #   k <- suppressWarnings(as.integer(lag_in))
  #   if (length(k) == 0 || is.na(k) || k < 0) k <- 0L
  #   
  #   # alpha mapping (robust)
  #   alpha_raw <- as.character(a_in)
  #   alpha_val <- to_num_safe(alpha_raw, default = 0.05)
  #   
  #   alpha_col <- switch(alpha_raw,
  #                       "0.01" = "1pct",
  #                       "0.05" = "5pct",
  #                       "0.1"  = "10pct"
  #   )
  #   
  #   # tau row mapping
  #   tau_row <- switch(type_in,
  #                     "none"  = "tau1",
  #                     "drift" = "tau2",
  #                     "trend" = "tau3",
  #                     "tau3"
  #   )
  #   
  #   # ============================================================================
  #   # SECTION 1: DATA CONTEXT & SPECIFICATION
  #   # ============================================================================
  #   cat("==========================================================================\n")
  #   cat("                        ADF UNIT ROOT DIAGNOSTIC                          \n")
  #   cat("==========================================================================\n")
  #   cat(sprintf(" MODEL TYPE : %-10s | SAMPLE SIZE (N) : %d\n", toupper(type_in), valid_N))
  #   cat(sprintf(" LAG ORDER  : %-10d | SIGNIFICANCE (α) : %s\n", k, alpha_raw))
  #   cat(sprintf(" MOMENTS    : Mean: %.4f | Std.Dev: %.4f\n", mean(x), stats::sd(x)))
  #   # cat("--------------------------------------------------------------------------\n")
  #   cat("\n")
  #   # Basic sanity
  #   if (valid_N < 5) {
  #     cat(" [!] CRITICAL ERROR: Too few observations (N < 5).\n")
  #     cat("     DIRECTIVE: Provide more data points.\n")
  #     cat("==========================================================================\n")
  #     return(invisible(NULL))
  #   }
  #   
  #   if (!is.finite(stats::sd(x)) || stats::sd(x) == 0) {
  #     cat(" [!] CRITICAL ERROR: Series is constant or invalid (sd = 0 / NA).\n")
  #     cat("     DIRECTIVE: Check transformation (log/diff) or data quality.\n")
  #     cat("==========================================================================\n")
  #     return(invisible(NULL))
  #   }
  #   
  #   # Degrees of Freedom Safety Check (conservative)
  #   if (valid_N <= (k + 10)) {
  #     cat(" [!] WARNING: Sample size is small relative to selected lag.\n")
  #     cat("     This can break tseries::adf.test and weaken inference.\n")
  #     cat("     DIRECTIVE: Decrease Lag (k) or increase data points.\n")
  #     # cat("--------------------------------------------------------------------------\n")
  #     cat("\n")
  #   }
  #   
  #   # ============================================================================
  #   # 3) PREDECLARE OBJECTS (prevents crashes in post-summary)
  #   # ============================================================================
  #   res_urca <- NULL
  #   
  #   tau_obs  <- NA_real_
  #   tau_crit <- NA_real_
  #   
  #   res_tseries <- list(p.value = NA_real_)
  #   lb_test     <- list(statistic = NA_real_, p.value = NA_real_)
  #   kpss_type   <- if (type_in == "trend") "Trend" else "Level"
  #   res_kpss    <- list(statistic = NA_real_, p.value = NA_real_)
  #   
  #   is_stationary   <- FALSE
  #   kpss_stationary <- FALSE
  #   
  #   # Pettitt placeholders (ONLY — Buishand removed)
  #   pettitt_res <- list(statistic = NA_real_, p.value = NA_real_, estimate = NULL)
  #   
  #   # ============================================================================
  #   # 4) COMPUTATIONS (tryCatch keeps UI alive, cats remain)
  #   # ============================================================================
  #   tryCatch({
  #     
  #     # --- Package checks (explicit + friendly) ---
  #     if (!requireNamespace("urca", quietly = TRUE)) {
  #       stop("Package 'urca' is not installed. ACTION: install.packages('urca')")
  #     }
  #     if (!requireNamespace("tseries", quietly = TRUE)) {
  #       stop("Package 'tseries' is not installed. ACTION: install.packages('tseries')")
  #     }
  #     
  #     # ==========================================================================
  #     # PRIMARY TEST: ADF via urca::ur.df
  #     # ==========================================================================
  #     res_urca <- urca::ur.df(x, type = type_in, lags = k)
  #     
  #     # Tau statistic: extract by mapped row safely
  #     tau_obs <- suppressWarnings(as.numeric(res_urca@teststat[tau_row]))
  #     if (!is.finite(tau_obs)) {
  #       # fallback attempt (some objects behave differently in edge cases)
  #       tau_obs_fallback <- suppressWarnings(as.numeric(res_urca@teststat[1]))
  #       if (is.finite(tau_obs_fallback)) tau_obs <- tau_obs_fallback
  #     }
  #     
  #     # Critical value
  #     tau_crit <- suppressWarnings(as.numeric(res_urca@cval[tau_row, alpha_col]))
  #     
  #     # Secondary ADF (p-value reference) — only if feasible
  #     if (valid_N > (k + 10)) {
  #       res_tseries <- tseries::adf.test(x, alternative = alt_in, k = k)
  #     } else {
  #       res_tseries <- list(p.value = NA_real_)
  #     }
  #     
  #     # ==========================================================================
  #     # QUALITY CHECK: Ljung-Box on ur.df residuals
  #     # ==========================================================================
  #     lb_lag <- max(1L, min(10L, floor(valid_N / 5)))
  #     lb_test <- Box.test(res_urca@res, lag = lb_lag, type = "Ljung-Box")
  #     
  #     # ==========================================================================
  #     # CONFIRMATORY: KPSS (H0 = stationary)
  #     # ==========================================================================
  #     kpss_type <- if (type_in == "trend") "Trend" else "Level"
  #     res_kpss  <- tseries::kpss.test(x, null = kpss_type)
  #     
  #     # ==========================================================================
  #     # SAFE FLAGS (avoid NA inside if-conditions)
  #     # ==========================================================================
  #     tau_ok  <- is.finite(tau_obs) && is.finite(tau_crit)
  #     lb_ok   <- is.finite(to_num_safe(lb_test$p.value))
  #     kpss_ok <- is.finite(to_num_safe(res_kpss$p.value))
  #     
  #     is_stationary <- isTRUE(tau_ok) && (tau_obs < tau_crit)
  #     
  #     if (kpss_ok && (to_num_safe(res_kpss$p.value) > alpha_val)) {
  #       kpss_stationary <- TRUE
  #     } else {
  #       kpss_stationary <- FALSE
  #     }
  #     
  #     lb_white_safe <- lb_ok && (to_num_safe(lb_test$p.value) > alpha_val)
  #     
  #     # ==========================================================================
  #     # SECTION 2: PRIMARY TEST - AUGMENTED DICKEY-FULLER (ADF)
  #     # ==========================================================================
  #     cat("==========================================================================\n")
  #     cat("PHASE 1: ADF UNIT ROOT TEST\n")
  #     cat("==========================================================================\n")
  #     
  #     cat(" • H0: The series has a Unit Root (Non-Stationary).\n")
  #     cat(" • Ha: The series is Stationary (Mean Reverting).\n")
  #     cat(sprintf(" -> CRITERIA: Reject H0 if Tau-Obs (%.4f) < Tau-Crit (%.4f)\n",
  #                 tau_obs, tau_crit))
  #     
  #     cat("\n RESULT:\n")
  #     cat(paste("  - Tau Observed :", round(tau_obs, 4), "\n"))
  #     cat(paste("  - Tau Critical :", round(tau_crit, 4), "\n"))
  #     cat(paste("  - P-Value (Ref):", round(to_num_safe(res_tseries$p.value), 4), "\n"))
  #     
  #     cat("\n DECISION:\n")
  #     if (isTRUE(tau_ok) && isTRUE(is_stationary)) {
  #       cat("  -> REJECT H0: Evidence suggests the series is STATIONARY.\n")
  #     } else if (isTRUE(tau_ok) && !isTRUE(is_stationary)) {
  #       cat("  -> FAIL TO REJECT H0: Evidence suggests the series is NON-STATIONARY.\n")
  #     } else {
  #       cat("  -> [!] WARNING: Tau/Critical value is NA/Inf. ADF inference is not valid.\n")
  #     }
  #     # cat("--------------------------------------------------------------------------\n")
  #     cat("\n")
  #     
  #     # ==========================================================================
  #     # SECTION 3: QUALITY CHECK - LJUNG-BOX TEST
  #     # ==========================================================================
  #     cat("==========================================================================\n")
  #     cat("PHASE 2: RESIDUAL DIAGNOSTICS (LJUNG-BOX)\n")
  #     cat("==========================================================================\n")
  #     
  #     cat(" • H0: Residuals are White Noise (No Autocorrelation).\n")
  #     cat(" • Ha: Residuals are Correlated (Lags are insufficient).\n")
  #     cat(sprintf(" -> CRITERIA: Reject H0 if P-Value (%.4f) < α (%.2f)\n",
  #                 to_num_safe(lb_test$p.value), alpha_val))
  #     
  #     cat("\n RESULT:\n")
  #     cat(paste("  - LB Statistic :", round(to_num_safe(lb_test$statistic), 4), "\n"))
  #     cat(paste("  - LB P-Value   :", round(to_num_safe(lb_test$p.value), 4), "\n"))
  #     
  #     cat("\n DECISION:\n")
  #     if (isTRUE(lb_white_safe)) {
  #       cat("  -> FAIL TO REJECT H0: Residuals are White Noise. [ADF VALID]\n")
  #     } else if (isTRUE(lb_ok)) {
  #       cat("  -> REJECT H0: Residuals are Correlated. [ADF BIASED]\n")
  #     } else {
  #       cat("  -> [!] WARNING: Ljung-Box P-Value is NA/Inf. Residual diagnosis failed.\n")
  #     }
  #     # cat("--------------------------------------------------------------------------\n")
  #     cat("\n")
  #     
  #     
  #     
  #     # ==========================================================================
  #     # PHASE 3: KPSS + STRUCTURAL BREAK CHECK (Pettitt) + Segment KPSS
  #     # ==========================================================================
  #     
  #     cat("==========================================================================\n")
  #     cat("PHASE 3: CONFIRMATORY ANALYSIS (KPSS)+ STRUCTURAL BREAK CHECK (Pettitt) + Segment KPSS\n")
  #     cat("==========================================================================\n\n")
  #     cat("PHASE 3A: CONFIRMATORY ANALYSIS (KPSS)\n")
  #     cat(paste(" • H0: The series is Stationary around a", kpss_type, ".\n"))
  #     cat(" • Ha: The series is Non-Stationary (Unit Root).\n")
  #     cat(sprintf(" • CRITERIA: Reject H0 if P-Value (%.4f) < α (%.2f)\n",
  #                 to_num_safe(res_kpss$p.value), alpha_val))
  #     
  #     cat("\n RESULT:\n")
  #     cat(paste("  - KPSS P-Value :", round(to_num_safe(res_kpss$p.value), 4), "\n"))
  #     
  #     cat("\n DECISION:\n")
  #     if (isTRUE(kpss_ok) && isTRUE(kpss_stationary)) {
  #       cat("  -> FAIL TO REJECT H0: The series is STATIONARY. [KPSS]\n")
  #     } else if (isTRUE(kpss_ok) && !isTRUE(kpss_stationary)) {
  #       cat("  -> REJECT H0: The series is NON-STATIONARY. [KPSS]\n")
  #     } else {
  #       cat("  -> [!] WARNING: KPSS P-Value is NA/Inf. KPSS inference is not valid.\n")
  #     }
  #     cat("--------------------------------------------------------------------------\n")
  #     
  #     # ==========================================================================
  #     # PHASE 3B: STRUCTURAL BREAK CHECK (Pettitt) + KPSS re-check by segments
  #     # ==========================================================================
  #     cat("PHASE 3B: STRUCTURAL BREAK CHECK (Pettitt) + KPSS SEGMENT CHECK\n")
  #     cat(" • Goal: Detect a single change-point (shift) that can contaminate KPSS.\n")
  #     cat(" • If a break exists, we re-run KPSS before/after the break.\n")
  #     
  #     pettitt_break <- FALSE
  #     break_idx <- NA_integer_
  #     
  #     if (requireNamespace("trend", quietly = TRUE)) {
  #       
  #       # Pettitt safely (never kill the rest of the report)
  #       tryCatch({
  #         pettitt_res <- trend::pettitt.test(x)
  #       }, error = function(e) {
  #         pettitt_res <- list(statistic = NA_real_, p.value = NA_real_, estimate = NULL)
  #         cat(" [!] WARNING: Pettitt test failed.\n")
  #         cat("     ERROR: ", e$message, "\n", sep = "")
  #       })
  #       
  #       pett_p <- to_num_safe(pettitt_res$p.value)
  #       cat(sprintf(" • Pettitt P-Value: %.6f | α: %.4f\n", pett_p, alpha_val))
  #       
  #       # estimate is often the location (index)
  #       if (!is.null(pettitt_res$estimate) && is.finite(as.numeric(pettitt_res$estimate))) {
  #         break_idx <- as.integer(as.numeric(pettitt_res$estimate))
  #       }
  #       
  #       pettitt_break <- is.finite(pett_p) && (pett_p < alpha_val) && is.finite(break_idx)
  #       
  #       cat("\n RESULT:\n")
  #       cat(paste("  - Break index (estimate):", ifelse(is.finite(break_idx), break_idx, "NA"), "\n"))
  #       cat(paste("  - Pettitt p-value       :", round(pett_p, 6), "\n"))
  #       
  #       cat("\n DECISION:\n")
  #       if (pettitt_break) {
  #         cat("  -> REJECT H0: Break detected. KPSS may be distorted on full sample.\n")
  #       } else {
  #         cat("  -> FAIL TO REJECT H0: No strong evidence of a single break.\n")
  #       }
  #       
  #       # ------------------------------------------------------------
  #       # Segment KPSS (only if break is detected and segments are big enough)
  #       # ------------------------------------------------------------
  #       if (pettitt_break) {
  #         
  #         # Split around break_idx
  #         x1 <- x[1:break_idx]
  #         x2 <- x[(break_idx + 1):valid_N]
  #         
  #         cat("--------------------------------------------------------------------------\n")
  #         
  #         cat("PHASE 3C: KPSS RE-CHECK BY SEGMENTS:\n")
  #         
  #         cat("  • Segment 1 = [1 .. break]\n")
  #         cat("  • Segment 2 = [break+1 .. N]\n")
  #         
  #         # Minimum length guard (KPSS needs enough points)
  #         min_seg <- 12L
  #         if (length(x1) < min_seg || length(x2) < min_seg) {
  #           cat("  [!] WARNING: One segment is too short for reliable KPSS.\n")
  #           cat(sprintf("     Segment lengths: n1=%d, n2=%d (min recommended=%d)\n",
  #                       length(x1), length(x2), min_seg))
  #         } else {
  #           
  #           # Run KPSS on each segment safely
  #           kpss1 <- tryCatch(tseries::kpss.test(x1, null = kpss_type),
  #                             error = function(e) list(p.value = NA_real_))
  #           kpss2 <- tryCatch(tseries::kpss.test(x2, null = kpss_type),
  #                             error = function(e) list(p.value = NA_real_))
  #           
  #           p1 <- to_num_safe(kpss1$p.value)
  #           p2 <- to_num_safe(kpss2$p.value)
  #           
  #           cat(sprintf("  - KPSS p-value (Segment 1): %.6f\n", p1))
  #           cat(sprintf("  - KPSS p-value (Segment 2): %.6f\n", p2))
  #           
  #           # Interpret segment KPSS
  #           s1_stat <- is.finite(p1) && (p1 > alpha_val)
  #           s2_stat <- is.finite(p2) && (p2 > alpha_val)
  #           
  #           cat("\n INTERPRETATION:\n")
  #           if (s1_stat && s2_stat) {
  #             cat("  [✓] Both segments look stationary by KPSS.\n")
  #             cat("      -> Full-sample KPSS non-stationarity may be break-driven.\n")
  #           } else if (!s1_stat && !s2_stat) {
  #             cat("  [X] Both segments look non-stationary by KPSS.\n")
  #             cat("      -> Non-stationarity is not only due to a break.\n")
  #           } else {
  #             cat("  [?] Mixed: one segment stationary, the other not.\n")
  #             cat("      -> Consider regime modeling or re-check transformation.\n")
  #           }
  #         }
  #       }
  #       
  #     } else {
  #       cat(" [!] WARNING: Package 'trend' not installed → Pettitt break check skipped.\n")
  #       cat("     ACTION: install.packages('trend')\n")
  #     }
  #     
  #     
  #     
  #     # cat("--------------------------------------------------------------------------\n")
  #     
  #     
  #     
  #     
  #     # ==========================================================================
  #     # SECTION 5: FINAL VERDICT & STUDENT DIRECTIVES
  #     # ==========================================================================
  #     cat("==========================================================================\n")
  #     cat("PHASE 4: FINAL ACADEMIC VERDICT & ADVICE\n")
  #     cat("==========================================================================\n")
  #     
  #     
  #     if (isTRUE(lb_ok) && !isTRUE(lb_white_safe)) {
  #       cat(" [!] WARNING: Your ADF 'Tau' statistic is technically biased.\n")
  #       cat("     ADVICE: Increase 'Lag Order' until Ljung-Box P-Value > α.\n\n")
  #     }
  #     
  #     if (isTRUE(is_stationary) && isTRUE(kpss_stationary)) {
  #       cat(" [✓] VERDICT: STRONG STATIONARITY confirmed by both ADF and KPSS.\n")
  #       cat("     ACTION: Proceed to ARMA/ARIMA modeling with d=0.\n")
  #     } else if (isTRUE(!is_stationary) && isTRUE(!kpss_stationary)) {
  #       cat(" [X] VERDICT: CLEAR UNIT ROOT. Both tests confirm Non-Stationarity.\n")
  #       cat("     ACTION: Apply differencing (d=1) to stabilize the mean.\n")
  #     } else {
  #       cat(" [?] VERDICT: CONFLICTING RESULTS (ADF vs KPSS).\n")
  #       cat("     ADVICE: The series may be highly persistent or trend-stationary.\n")
  #       cat("     Double-check plots for breaks or changing variance.\n")
  #     }
  #     
  #     cat("\n TECHNICAL APPENDIX (ADF Regression Coefficients):\n")
  #     if (!is.null(res_urca)) {
  #       print(stats::coef(res_urca@testreg))
  #     } else {
  #       cat("  [!] Not available (ur.df did not run).\n")
  #     }
  #     
  #   }, error = function(e) {
  #     cat(" EXECUTION ERROR: ", e$message, "\n")
  #   })
  #   
  #   # cat("--------------------------------------------------------------------------\n")
  #   cat("\n")
  #   
  #   # ============================================================================
  #   # 5) POST-SUMMARY (always runs; safe against NA objects)
  #   # ============================================================================
  #   cat("==========================================================================\n")
  #   cat("PHASE 5: POST-SUMMARY (Academic-quality snapshot)\n")
  #   cat("==========================================================================\n")
  #   
  #   cat(" EVIDENCE SNAPSHOT (All key outcomes in one place):\n")
  #   cat(sprintf(" [ ] N (effective sample size)        : %d\n", valid_N))
  #   cat(sprintf(" [ ] Model type (ADF)                 : %s  (tau row: %s)\n", type_in, tau_row))
  #   cat(sprintf(" [ ] Lag order (k)                    : %d\n", k))
  #   cat(sprintf(" [ ] Alpha (α)                        : %.4f\n", alpha_val))
  #   cat(sprintf(" [ ] Tau-Observed (urca)              : %.6f\n", tau_obs))
  #   cat(sprintf(" [ ] Tau-Critical (urca, %s)          : %.6f\n", alpha_col, tau_crit))
  #   cat(sprintf(" [ ] ADF p-value (tseries reference)  : %.6f\n", to_num_safe(res_tseries$p.value)))
  #   cat(sprintf(" [ ] Ljung-Box p-value (residuals)    : %.6f\n", to_num_safe(lb_test$p.value)))
  #   cat(sprintf(" [ ] KPSS p-value (%s null)           : %.6f\n", kpss_type, to_num_safe(res_kpss$p.value)))
  #   cat("--------------------------------------------------------------------------\n")
  #   
  #   cat(" CHECKLIST (Academic-quality acceptance criteria):\n")
  #   cat(sprintf(" [ ] Tau is finite and usable                 : %s\n",
  #               ifelse(is.finite(tau_obs), "[✓] SATISFIED", "[!] NOT SATISFIED")))
  #   if (!is.finite(tau_obs)) {
  #     cat("     [!] DANGER: Tau is NA/Inf → ADF inference cannot be trusted.\n")
  #     cat("     [!] Likely causes: constant/near-constant series, k too high,\n")
  #     cat("         wrong extraction index, too small effective N.\n")
  #   }
  #   
  #   cat(sprintf(" [ ] Tau critical value extracted             : %s\n",
  #               ifelse(is.finite(tau_crit), "[✓] SATISFIED", "[!] NOT SATISFIED")))
  #   if (!is.finite(tau_crit)) {
  #     cat("     [!] DANGER: Critical value missing → alpha/type mapping issue.\n")
  #     cat("     [!] ACTION: verify alpha_col and tau_row exist in urca cval table.\n")
  #   }
  #   
  #   lb_p <- to_num_safe(lb_test$p.value)
  #   cat(sprintf(" [ ] Residuals pass Ljung-Box (white noise)   : %s\n",
  #               ifelse(is.finite(lb_p) && lb_p > alpha_val, "[✓] SATISFIED", "[X] NOT SATISFIED")))
  #   if (is.finite(lb_p) && lb_p <= alpha_val) {
  #     cat("     [!] WARNING: residual autocorrelation detected.\n")
  #     cat("     [!] Interpretation: ADF may be biased (insufficient k).\n")
  #     cat("     [!] Suggested fix: increase k cautiously OR re-check transformation.\n")
  #   }
  #   
  #   cat(sprintf(" [ ] Sample size adequacy                     : %s\n",
  #               ifelse(valid_N >= 30, "[✓] STRONG", ifelse(valid_N >= 15, "[?] BORDERLINE", "[X] WEAK"))))
  #   if (valid_N < 15) {
  #     cat("     [!] WARNING: very low power. ADF/KPSS decisions may be unstable.\n")
  #     cat("     [!] Advice: combine with plots (time plot + ACF/PACF) and be conservative.\n")
  #   }
  #   
  #   
  #   
  #   #  lag sanity
  #   cat(sprintf(" [ ] Lag order reasonable relative to N       : %s\n",
  #               ifelse(valid_N > (k + 5), "[✓] OK", "[!] TOO HIGH / RISKY")))
  #   if (valid_N <= (k + 5)) {
  #     cat("     [!] DANGER: k too high for N → regression can break or give NA stats.\n")
  #     cat("     [!] ACTION: reduce k immediately.\n")
  #   }
  #   
  #   
  #   #  cross-test agreement (ADF vs KPSS)
  #   agreement_safe <- (isTRUE(is_stationary) && isTRUE(kpss_stationary)) ||
  #     (isTRUE(!is_stationary) && isTRUE(!kpss_stationary))
  #   
  #   cat(sprintf(" [ ] ADF & KPSS agreement                      : %s\n",
  #               ifelse(agreement_safe, "[✓] AGREEMENT", "[?] CONFLICT")))
  #   
  #   if (!agreement_safe) {
  #     cat("     [?] NOTE: conflict can happen for near-unit-root series, structural breaks,\n")
  #     cat("         or when model type (none/drift/trend) is mis-specified.\n")
  #   }
  #   
  #   
  #   # Structural break summary (Pettitt only)
  #   if (requireNamespace("trend", quietly = TRUE)) {
  #     cat("--------------------------------------------------------------------------\n")
  #     cat(" STRUCTURAL BREAK SUMMARY (Pettitt):\n")
  #     cat(sprintf(" [ ] Pettitt p-value   : %.6f  (Reject H0 if < α)\n", to_num_safe(pettitt_res$p.value)))
  #     
  #     pett_break2 <- is.finite(to_num_safe(pettitt_res$p.value)) &&
  #       (to_num_safe(pettitt_res$p.value) < alpha_val)
  #     
  #     cat("\n INTERPRETATION:\n")
  #     if (pett_break2) {
  #       cat(" [!] WARNING: Evidence of a break detected (Pettitt).\n")
  #       cat("     ADVICE: Consider split-sample testing or break-aware modeling.\n")
  #       cat("     NOTE: Breaks can cause ADF/KPSS instability or conflicts.\n")
  #     } else {
  #       cat(" [✓] No strong evidence of a major break at α (Pettitt).\n")
  #       cat("     ADVICE: Stationarity conclusions are less likely to be break-contaminated.\n")
  #     }
  #   } else {
  #     cat("--------------------------------------------------------------------------\n")
  #     cat(" STRUCTURAL BREAK SUMMARY (Pettitt):\n")
  #     cat(" [!] WARNING: Package 'trend' is not installed. Pettitt summary unavailable.\n")
  #     cat("     ACTION: install.packages('trend')\n")
  #   }
  #   
  #   cat("--------------------------------------------------------------------------\n")
  #   
  #   if (isTRUE(is_stationary) && isTRUE(kpss_stationary)) {
  #     cat(" [✓] VERDICT: CONSISTENT STATIONARITY. Both tests agree the series is I(0).\n")
  #     cat(" ADVICE: You may proceed to model identification using the original series.\n")
  #   } else if (isTRUE(!is_stationary) && isTRUE(!kpss_stationary)) {
  #     cat(" [X] VERDICT: CONSISTENT UNIT ROOT. Both tests confirm Non-Stationarity.\n")
  #     cat(" ADVICE: Apply First Differencing (d=1) to achieve stationarity.\n")
  #   } else {
  #     cat(" [?] VERDICT: CONFLICTING RESULTS. One test suggests a Unit Root while the other doesn't.\n")
  #     cat(" ADVICE: This often happens with near-integrated or break-contaminated series. \n")
  #     cat("         Use differencing for safety.\n")
  #   }
  #   
  #   cat("\n ACTIONABLE NEXT STEPS (What to do now):\n")
  #   
  #   if (!is.finite(tau_obs) || !is.finite(tau_crit)) {
  #     cat(" [!] PRIORITY 1: Fix the ADF computation first.\n")
  #     cat("     • Reduce k, verify variance is not zero, and confirm urca outputs.\n")
  #     cat("     • If Tau remains NA, your series may be too short after transformation.\n")
  #   } else {
  #     cat(" [✓] PRIORITY 1: Tau and critical values are usable.\n")
  #   }
  #   
  #   if (is.finite(lb_p) && lb_p <= alpha_val) {
  #     cat(" [X] PRIORITY 2: Residual autocorrelation detected.\n")
  #     cat("     • Increase the Lag k (carefully) OR re-check transformation (d/D/log).\n")
  #   } else {
  #     cat(" [✓] PRIORITY 2: Residual whiteness acceptable → ADF inference more reliable.\n")
  #   }
  #   
  #   if (!agreement_safe) {
  #     cat(" [?] PRIORITY 3: Resolve ADF vs KPSS conflict.\n")
  #     cat("     • Try alternative model types (none/drift/trend).\n")
  #     cat("     • Consider break-aware modeling if Pettitt flags a break.\n")
  #     cat("     • Use visual diagnostics + conservative differencing.\n")
  #   } else {
  #     cat(" [✓] PRIORITY 3: Tests agree → proceed with confidence.\n")
  #   }
  #   
  #   cat("\n PRACTICAL MODELING PATH (for your Shiny workflow):\n")
  #   if (isTRUE(is_stationary) && is.finite(lb_p) && (lb_p > alpha_val)) {
  #     cat(" [✓] Use ARMA identification on current series → fit → residual analysis.\n")
  #   } else if (!isTRUE(is_stationary)) {
  #     cat(" [X] Apply differencing (d and/or D) → re-run ADF/KPSS → then identify ARMA.\n")
  #   } else {
  #     cat(" [?] Mixed evidence: consider differencing for safety and validate visually.\n")
  #   }
  #   
  #   cat("--------------------------------------------------------------------------\n\n")
  # })
  
  
  
  
  
  
  
  
  # ====================================================================
  # ====================================================================
  # ====================================================================
  
  
  
  
  output$CHECKLIST <- renderPrint({
    
    # ============================================================================
    # 0) SMALL HELPERS
    # ============================================================================
    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
    
    to_num_safe <- function(v, default = NA_real_) {
      out <- suppressWarnings(as.numeric(v))
      if (length(out) == 0 || all(is.na(out)) || !is.finite(out[1])) default else out[1]
    }
    
    to_int_safe <- function(v, default = 0L) {
      out <- suppressWarnings(as.integer(v))
      if (length(out) == 0 || is.na(out[1]) || !is.finite(out[1])) default else out[1]
    }
    
    safe_head_tail <- function(x, n = 5) {
      x <- as.numeric(x)
      x <- x[is.finite(x)]
      if (length(x) == 0) return(list(head = numeric(0), tail = numeric(0)))
      list(head = head(x, n), tail = tail(x, n))
    }
    
    # Robust extractor for urca @cval (vector OR matrix; row/col name variations)
    get_uc_cval <- function(cv, key) {
      if (is.null(cv) || is.null(key) || !nzchar(key)) return(NA_real_)
      
      if (is.matrix(cv)) {
        if (!is.null(colnames(cv)) && key %in% colnames(cv)) return(as.numeric(cv[1, key]))
        if (!is.null(rownames(cv)) && key %in% rownames(cv)) return(as.numeric(cv[key, 1]))
        
        if (!is.null(colnames(cv))) {
          m <- which(grepl(key, colnames(cv), fixed = TRUE))
          if (length(m) > 0) return(as.numeric(cv[1, m[1]]))
        }
        if (!is.null(rownames(cv))) {
          m <- which(grepl(key, rownames(cv), fixed = TRUE))
          if (length(m) > 0) return(as.numeric(cv[m[1], 1]))
        }
        return(NA_real_)
      }
      
      nm <- names(cv)
      if (!is.null(nm) && key %in% nm) return(as.numeric(cv[[key]]))
      if (!is.null(nm)) {
        m <- which(grepl(key, nm, fixed = TRUE))
        if (length(m) > 0) return(as.numeric(cv[[m[1]]]))
      }
      NA_real_
    }
    
    # ============================================================================
    # 1) INPUT COLLECTION (supports either naming convention in your UI)
    # ============================================================================
    alt_in  <- input$alternd2St %||% input$alternSt
    lag_in  <- input$LagOrderADFd2St %||% input$LagOrderADFSt
    a_in    <- input$alphaSt2
    type_in <- input$adfTypeSt2
    
    # Optional transformation UI inputs (informational only)
    d_in   <- input$d_n  %||% input$d  %||% NA
    D_in   <- input$DS_n %||% input$D  %||% NA
    log_in <- input$check_box %||% input$islog %||% FALSE
    
    req(myData_Choice())
    
    if (is.null(alt_in) || is.null(lag_in) || is.null(a_in) || is.null(type_in)) {
      cat("==========================================================================\n")
      cat(" [!] INPUT ERROR: One or more inputs are NULL (likely UI/server ID mismatch)\n")
      cat("     Needed inputs:\n")
      cat("       - alternative : input$alternd2St OR input$alternSt\n")
      cat("       - lag         : input$LagOrderADFd2St OR input$LagOrderADFSt\n")
      cat("       - alpha       : input$alphaSt2\n")
      cat("       - adf type    : input$adfTypeSt2\n")
      cat("==========================================================================\n")
      return(invisible(NULL))
    }
    
    # ============================================================================
    # 2) DATA PREP
    # ============================================================================
    x_raw <- myData_Choice()
    na_before <- sum(is.na(x_raw))
    
    x_class <- paste(class(x_raw), collapse = ", ")
    x_freq  <- NA_integer_
    if (inherits(x_raw, "ts")) x_freq <- tryCatch(stats::frequency(x_raw), error = function(e) NA_integer_)
    
    x <- as.numeric(stats::na.omit(x_raw))
    valid_N <- length(x)
    
    k <- to_int_safe(lag_in, default = 0L)
    if (!is.finite(k) || k < 0) k <- 0L
    
    alpha_raw <- as.character(a_in)
    alpha_val <- to_num_safe(alpha_raw, default = 0.05)
    
    alpha_col <- switch(alpha_raw,
                        "0.01" = "1pct",
                        "0.05" = "5pct",
                        "0.1"  = "10pct",
                        "0.10" = "10pct",
                        "1pct" = "1pct",
                        "5pct" = "5pct",
                        "10pct"= "10pct",
                        "5pct")
    
    tau_row <- switch(type_in,
                      "none"  = "tau1",
                      "drift" = "tau2",
                      "trend" = "tau3",
                      "tau3")
    
    D_ui <- to_int_safe(D_in, default = 0L)
    seasonality_resolved <- isTRUE(D_ui > 0L)
    
    # Quick sanity exits (still only printing the 3 sections)
    if (valid_N < 5 || !is.finite(stats::sd(x)) || stats::sd(x) == 0) {
      cat("==========================================================================\n")
      cat(" CHECKLIST\n")
      cat("==========================================================================\n\n")
      cat(sprintf(" [ ] N (effective)                    : %d\n", valid_N))
      cat(sprintf(" [ ] Variance usable (sd>0)           : %s\n",
                  ifelse(is.finite(stats::sd(x)) && stats::sd(x) > 0, "[✓]", "[!]")))
      cat("--------------------------------------------------------------------------\n")
      cat(" FINAL ACADEMIC ADVICE\n")
      cat("--------------------------------------------------------------------------\n\n")
      if (valid_N < 5) {
        cat(" [!] Too few observations after transformation. Provide more data points.\n")
      } else {
        cat(" [!] Series is constant/invalid after transformation (sd=0/NA).\n")
        cat("     Fix data quality or transformation pipeline in myData_Choice().\n")
      }
      cat("--------------------------------------------------------------------------\n")
      cat(" ACTIONABLE NEXT STEPS\n")
      cat("--------------------------------------------------------------------------\n")
      cat(" [!] [1] Verify myData_Choice() returns the transformed series (log/d/D).\n")
      cat(" [!] [2] Print/plot the transformed series to confirm it is not constant.\n")
      cat(" [!] [3] Re-run tests after fixing the above.\n")
      cat("==========================================================================\n")
      return(invisible(NULL))
    }
    
    # ============================================================================
    # 3) PREDECLARE OUTPUTS
    # ============================================================================
    res_urca <- NULL
    res_tseries <- list(p.value = NA_real_)
    lb_test <- list(statistic = NA_real_, p.value = NA_real_)
    lb_lag <- NA_integer_
    
    res_kpss_ts <- list(statistic = NA_real_, p.value = NA_real_)
    res_kpss_uc <- NULL
    
    tau_obs <- NA_real_
    tau_crit <- NA_real_
    
    eta_obs_uc <- NA_real_
    eta_crit_uc <- NA_real_
    eta_obs_ts <- NA_real_
    eta_p_one <- NA_real_
    eta_col <- NA_character_
    
    pettitt_res <- list(statistic = NA_real_, p.value = NA_real_, estimate = NULL)
    
    is_stationary <- FALSE
    kpss_stationary <- FALSE
    agreement_safe <- FALSE
    lb_white_safe <- FALSE
    
    # ============================================================================
    # 4) COMPUTE TESTS (no verbose phases; we only use results for the 3 outputs)
    # ============================================================================
    tryCatch({
      
      if (!requireNamespace("urca", quietly = TRUE)) stop("Package 'urca' missing. install.packages('urca')")
      if (!requireNamespace("tseries", quietly = TRUE)) stop("Package 'tseries' missing. install.packages('tseries')")
      
      # ADF (urca)
      res_urca <- urca::ur.df(x, type = type_in, lags = k)
      
      tau_obs <- suppressWarnings(as.numeric(res_urca@teststat[tau_row]))
      if (!is.finite(tau_obs)) {
        tau_obs_fallback <- suppressWarnings(as.numeric(res_urca@teststat[1]))
        if (is.finite(tau_obs_fallback)) tau_obs <- tau_obs_fallback
      }
      tau_crit <- suppressWarnings(as.numeric(res_urca@cval[tau_row, alpha_col]))
      
      # ADF p-value reference (tseries)
      if (valid_N > (k + 10)) {
        res_tseries <- tseries::adf.test(x, alternative = alt_in, k = k)
      }
      
      # Ljung-Box on ADF residuals
      lb_lag <- max(1L, min(10L, floor(valid_N / 5)))
      lb_test <- Box.test(res_urca@res, lag = lb_lag, type = "Ljung-Box")
      
      # KPSS (tseries p-value, urca cvals)
      kpss_type <- if (type_in == "trend") "Trend" else "Level"
      res_kpss_ts <- tseries::kpss.test(x, null = kpss_type)
      eta_obs_ts  <- to_num_safe(res_kpss_ts$statistic)
      eta_p_one   <- to_num_safe(res_kpss_ts$p.value)
      
      kpss_type_urca <- if (type_in == "trend") "tau" else "mu"
      res_kpss_uc    <- urca::ur.kpss(x, type = kpss_type_urca)
      
      eta_obs_uc <- suppressWarnings(as.numeric(res_kpss_uc@teststat))
      
      eta_col <- switch(alpha_raw,
                        "0.01" = "1pct",
                        "0.05" = "5pct",
                        "0.1"  = "10pct",
                        "0.10" = "10pct",
                        "1pct" = "1pct",
                        "5pct" = "5pct",
                        "10pct"= "10pct",
                        "5pct")
      
      eta_crit_uc <- suppressWarnings(get_uc_cval(res_kpss_uc@cval, eta_col))
      
      # Optional Pettitt
      if (requireNamespace("trend", quietly = TRUE)) {
        pettitt_res <- tryCatch(trend::pettitt.test(x),
                                error = function(e) list(statistic = NA_real_, p.value = NA_real_, estimate = NULL))
      }
      
      # Decisions
      tau_ok <- is.finite(tau_obs) && is.finite(tau_crit)
      lb_p   <- to_num_safe(lb_test$p.value)
      lb_ok  <- is.finite(lb_p)
      
      is_stationary <- isTRUE(tau_ok) && (tau_obs < tau_crit)
      
      # KPSS rejects stationarity if p<alpha OR eta_obs>eta_crit (when available)
      kpss_reject_by_p   <- is.finite(eta_p_one) && (eta_p_one < alpha_val)
      kpss_reject_by_eta <- is.finite(eta_obs_uc) && is.finite(eta_crit_uc) && (eta_obs_uc > eta_crit_uc)
      kpss_stationary <- !(kpss_reject_by_p || kpss_reject_by_eta)
      
      lb_white_safe <- lb_ok && (lb_p > alpha_val)
      
      agreement_safe <- (isTRUE(is_stationary) && isTRUE(kpss_stationary)) ||
        (isTRUE(!is_stationary) && isTRUE(!kpss_stationary))
      
    }, error = function(e) {
      cat("==========================================================================\n")
      cat(" CHECKLIST\n")
      cat("==========================================================================\n")
      cat(" [!] Execution error while computing tests:\n")
      cat("     ", e$message, "\n", sep = "")
      cat("--------------------------------------------------------------------------\n")
      cat(" FINAL ACADEMIC ADVICE\n")
      cat("--------------------------------------------------------------------------\n")
      cat(" [!] Install/enable required packages and re-run.\n")
      cat("--------------------------------------------------------------------------\n")
      cat(" ACTIONABLE NEXT STEPS\n")
      cat("--------------------------------------------------------------------------\n")
      cat(" [!] [1] install.packages('urca') and install.packages('tseries')\n")
      cat(" [!] [2] Re-run after verifying myData_Choice() returns a numeric/ts vector.\n")
      cat("==========================================================================\n")
      return(invisible(NULL))
    })
    
    # ============================================================================
    # 5) OUTPUT ONLY: CHECKLIST, FINAL ACADEMIC ADVICE, ACTIONABLE NEXT STEPS
    # ============================================================================
    lb_p <- to_num_safe(lb_test$p.value)
    adf_p <- to_num_safe(res_tseries$p.value)
    pett_p <- to_num_safe(pettitt_res$p.value)
    
    # --------------------
    # CHECKLIST
    # --------------------
    cat("==========================================================================\n")
    cat(" CHECKLIST\n")
    cat("==========================================================================\n")
    
    # ht <- safe_head_tail(x, 5)
    
    # cat(sprintf(" [ ] N (effective)                              : %d\n", valid_N))
    # cat(sprintf(" [ ] ADF model type (none/drift/trend)          : %s\n", as.character(type_in)))
    # cat(sprintf(" [ ] Lag k                                      : %d\n", k))
    # cat(sprintf(" [ ] Alpha (α)                                  : %.4f\n", alpha_val))
    # cat(sprintf(" [ ] Data class                                 : %s\n", x_class))
    # cat(sprintf(" [ ] Frequency (if ts)                          : %s\n", ifelse(is.finite(x_freq), as.character(x_freq), "NA / not ts")))
    # cat(sprintf(" [ ] NA count before na.omit                     : %d\n", na_before))
    # cat(sprintf(" [ ] UI transform flags                          : log=%s | d=%s | D=%s\n",
    #             ifelse(isTRUE(log_in), "ON", "OFF"),
    #             ifelse(is.na(d_in), "NA", as.character(d_in)),
    #             ifelse(is.na(D_in), "NA", as.character(D_in))))
    # cat(sprintf(" [ ] First 5 values tested                       : %s\n", paste(round(ht$head, 4), collapse = ", ")))
    # cat(sprintf(" [ ] Last  5 values tested                       : %s\n", paste(round(ht$tail, 4), collapse = ", ")))
    # cat("--------------------------------------------------------------------------\n")
    
    # cat(sprintf(" [ ] ADF Tau observed is finite                  : %s\n", ifelse(is.finite(tau_obs), "[✓]", "[!]")))
    # cat(sprintf(" [ ] ADF Tau critical is finite                  : %s\n", ifelse(is.finite(tau_crit), "[✓]", "[!]")))
    # cat(sprintf(" [ ] Ljung-Box p-value is finite                 : %s\n", ifelse(is.finite(lb_p), "[✓]", "[!]")))
    # cat(sprintf(" [ ] KPSS Eta observed (urca) is finite          : %s\n", ifelse(is.finite(eta_obs_uc), "[✓]", "[!]")))
    # cat(sprintf(" [ ] KPSS Eta critical (urca) is finite          : %s\n", ifelse(is.finite(eta_crit_uc), "[✓]", "[!]" )))
    # cat(sprintf(" [ ] KPSS p-value (tseries) is finite            : %s\n", ifelse(is.finite(eta_p_one), "[✓]", "[!]" )))
    # cat("--------------------------------------------------------------------------\n")
    
    # Key decisions summary (compact)
    cat(sprintf(" [ ] ADF decision (reject unit root => stationary) : %s\n",
                ifelse(isTRUE(is_stationary), "[✓] STATIONARY", "[X] NON-STATIONARY")))
    cat(sprintf(" [ ] KPSS decision (fail reject => stationary)     : %s\n",
                ifelse(isTRUE(kpss_stationary), "[✓] STATIONARY", "[X] NON-STATIONARY")))
    cat(sprintf(" [ ] ADF vs KPSS agreement                         : %s\n",
                ifelse(isTRUE(agreement_safe), "[✓] AGREEMENT", "[?] CONFLICT")))
    cat(sprintf(" [ ] Residual whiteness (LB (Ljung–Box) p>α)       : %s\n",
                ifelse(is.finite(lb_p) && lb_p > alpha_val, "[✓] OK", ifelse(is.finite(lb_p), "[X] FAIL", "[?] UNKNOWN"))))
    cat(sprintf(" [ ] Seasonal differencing indicated (UI D>0)      : %s\n",
                ifelse(isTRUE(seasonality_resolved), "[✓] YES", "[!] NO / UNKNOWN")))
    cat(sprintf(" [ ] Pettitt break check available                 : %s\n",
                ifelse(requireNamespace("trend", quietly = TRUE), "[✓] YES", "[!] NO (trend pkg missing)")))
    # if (requireNamespace("trend", quietly = TRUE)) {
    #   cat(sprintf(" [ ] Pettitt p-value                               : %.6f\n", pett_p))
    # }
    
    cat("\n")
    
    # --------------------
    # FINAL ACADEMIC ADVICE
    # --------------------
    cat("==========================================================================\n")
    cat(" FINAL ACADEMIC ADVICE\n")
    cat("==========================================================================\n")
    cat("\n")
    # Explain conflict in a compact, decision-useful way
    if (isTRUE(agreement_safe) && isTRUE(is_stationary) && isTRUE(kpss_stationary)) {
      cat(" [✓] Strong evidence of stationarity (I(0)) from BOTH ADF and KPSS.\n")
      if (!(is.finite(lb_p) && lb_p > alpha_val)) {
        cat(" [!] But residual autocorrelation suggests your ADF lag may be too small.\n")
        cat("     Treat the ADF conclusion as less reliable until LB passes.\n")
      } else {
        cat(" [✓] Residuals are consistent with a well-specified ADF regression.\n")
      }
      cat("     Academic implication: proceed with ARMA/ARIMA using d=0 (and D as already applied).\n")
      
    } else if (isTRUE(agreement_safe) && !isTRUE(is_stationary) && !isTRUE(kpss_stationary)) {
      cat(" [X] Strong evidence of a unit root (non-stationarity) from BOTH tests.\n")
      cat("     Academic implication: apply differencing (d=1) and re-test.\n")
      if (is.finite(x_freq) && x_freq > 1 && !isTRUE(seasonality_resolved)) {
        cat(" [!] Seasonal frequency detected; consider seasonal differencing (D=1) before increasing k.\n")
      }
      
    } else {
      cat(" [?] ADF and KPSS are in conflict.\n")
      cat("     This is common when the series is near-unit-root, trend-stationary vs difference-stationary,\n")
      cat("     has structural breaks, seasonality not properly removed, or lag k is mis-specified.\n")
      
      if (requireNamespace("trend", quietly = TRUE) && is.finite(pett_p) && pett_p < alpha_val) {
        cat(" [!] Pettitt suggests a structural break (p < α). Breaks often cause ADF/KPSS disagreement.\n")
      }
      
      if (is.finite(x_freq) && x_freq > 1 && !isTRUE(seasonality_resolved)) {
        cat(" [!] Seasonal frequency detected but seasonal differencing is NOT indicated (D=0).\n")
        cat("     Unremoved seasonality often triggers KPSS non-stationarity while ADF looks borderline.\n")
      }
      
      if (is.finite(lb_p) && lb_p <= alpha_val) {
        cat(" [!] Ljung-Box fails => ADF regression residuals are autocorrelated; your ADF decision is less trustworthy.\n")
      }
      
      cat("     Academic implication: be conservative—prefer differencing (and/or D=1 if seasonal), then re-test.\n")
    }
    
    cat(sprintf("\n (Numbers) ADF tau=%.6f | tau_crit=%.6f | ADF p(ref)=%.6f\n", tau_obs, tau_crit, adf_p))
    cat(sprintf("           KPSS eta_obs=%.6f | eta_crit=%.6f | KPSS p=%.6f\n", eta_obs_uc, eta_crit_uc, eta_p_one))
    cat(sprintf("           LB (the Ljung–Box Q test) p=%.6f (lag=%s)\n", lb_p, ifelse(is.finite(lb_lag), as.character(lb_lag), "NA")))
    cat("            |__ it’s normal to choose LB lag by a rule-of-thumb (like min(10, floor(N/5)) \n")
    cat("\n")
    
    # --------------------
    # ACTIONABLE NEXT STEPS
    # --------------------
    cat("==========================================================================\n")
    cat(" ACTIONABLE NEXT STEPS\n")
    cat("==========================================================================\n")
    
    
    # [1] Residual autocorrelation
    cat("\n")
    if (is.finite(lb_p) && lb_p <= alpha_val) {
      cat(" [X] [1] FIX RESIDUAL AUTOCORRELATION (LB failed)\n")
      cat("     • Increase k gradually (k+1, k+2) and re-run.\n")
      cat("     • If k becomes large vs N, prefer (d=1) and/or (D=1 if seasonal) instead of pushing k.\n")
    } else if (is.finite(lb_p) && lb_p > alpha_val) {
      cat(" [✓] [1] RESIDUALS LOOK OK (LB passed)\n")
      cat("     • Your ADF regression is less likely biased by autocorrelation.\n")
    } else {
      cat(" [?] [1] RESIDUAL DIAGNOSTICS INCONCLUSIVE\n")
      cat("     • LB p-value is NA/Inf. Reduce k and re-run; verify x is not pathological.\n")
    }
    
    # [2] Agreement vs conflict
    cat("\n")
    if (!isTRUE(agreement_safe)) {
      cat(" [?] [2] RESOLVE THE ADF vs KPSS CONFLICT\n")
      cat("     • Re-test ADF with model types: none / drift / trend (choose what your plot suggests).\n")
      cat("     • If seasonal frequency exists, try D=1 (seasonal differencing) before debating k.\n")
      cat("     • If series is strictly positive, try log or Box-Cox for variance stabilization.\n")
      if (requireNamespace("trend", quietly = TRUE)) {
        if (is.finite(pett_p) && pett_p < alpha_val) {
          cat("     • Break detected: split sample around the break and re-run KPSS/ADF on each segment.\n")
        } else {
          cat("     • No strong break evidence at α: focus on model type (trend/drift) + seasonality + k.\n")
        }
      } else {
        cat("     • Install 'trend' to check breaks: install.packages('trend')\n")
      }
    } else {
      cat(" [✓] [2] TESTS AGREE — MOVE FORWARD\n")
      cat("     • Use this decision to pick d (and D if seasonal) for ARIMA/SARIMA modeling.\n")
    }
    
    # [3] Seasonality sanity with required rule: [✓] if resolved, else [!]
    cat("\n")
    season_tag <- if (isTRUE(seasonality_resolved)) "[✓]" else "[!]"
    cat(sprintf(" %s [3] SEASONALITY SANITY\n", season_tag))
    if (is.finite(x_freq) && x_freq > 1) {
      cat(sprintf("     • Frequency=%d detected.\n", x_freq))
      if (isTRUE(seasonality_resolved)) {
        cat("     • D>0 in UI indicates seasonality treatment is ON (as long as myData_Choice() applies it).\n")
      } else {
        cat("     • Consider D=1 if ACF shows seasonal spikes at lag = frequency, 2*frequency, ...\n")
      }
    } else {
      cat("     • Frequency not available: rely on plot/ACF to decide if seasonal differencing is needed.\n")
    }
    
    # [4] Explosive note
    cat("\n")
    if (identical(as.character(alt_in), "explosive")) {
      cat(" [!] [4] EXPLOSIVE MODE NOTE\n")
      cat("     • KPSS is not an explosive test. Use tseries ADF p-value + plots.\n")
    } else {
      cat(" [✓] [4] EXPLOSIVE MODE NOTE\n")
      cat("     • Standard stationarity workflow applies.\n")
    }
    
    
    
    # --------------------------------------------------------------------------
    # EXTRA ACTIONABLE NEXT STEPS (more complete workflow)
    # --------------------------------------------------------------------------
    
    # [5] Choose ADF model type systematically (avoid random picking)
    cat("\n")
    cat(" [!] [5] CHOOSE ADF MODEL TYPE SYSTEMATICALLY (avoid mis-specification)\n")
    cat("     • Use your time plot:\n")
    cat("       - Clear deterministic trend  → set type='trend'\n")
    cat("       - No clear trend, non-zero mean → set type='drift'\n")
    cat("       - Mean around ~0 (rare)      → set type='none'\n")
    cat("     • Wrong type_in is a top cause of ADF vs KPSS conflict.\n")
    
    # [6] Lag strategy: use Ljung-Box as your guardrail
    cat("\n")
    cat(" [!] [6] USE LB AS A LAG-SELECTION GUARDRAIL\n")
    cat("     • Increase k until LB p-value > α (residuals approx. white).\n")
    cat("     • Stop increasing k if N becomes too small relative to k (risk: N <= k+10).\n")
    cat("     • If you hit that risk, prefer differencing (d or D) instead of more k.\n")
    
    # [7] Seasonality protocol (if frequency known)
    cat("\n")
    if (is.finite(x_freq) && x_freq > 1) {
      cat(" [!] [7] SEASONALITY PROTOCOL (freq detected)\n")
      cat("     • Check ACF for spikes at seasonal lags: freq, 2*freq, ...\n")
      if (isTRUE(seasonality_resolved)) {
        cat("     [✓] D>0 indicated → ensure myData_Choice() truly applied seasonal differencing.\n")
      } else {
        cat("     [X] D=0 indicated → if seasonal spikes exist, set D=1 and re-test.\n")
      }
    }
    
    # [8] Break protocol (if Pettitt available)
    cat("\n")
    if (requireNamespace("trend", quietly = TRUE)) {
      cat(" [!] [8] BREAK PROTOCOL (Pettitt)\n")
      if (is.finite(pett_p) && (pett_p < alpha_val)) {
        cat("     [!] Break detected → do this:\n")
        cat("       1) Split the series around the estimated break index.\n")
        cat("       2) Re-run KPSS/ADF on each segment.\n")
        cat("       3) If segments are stationary but full sample is not → break-driven non-stationarity.\n")
      } else {
        cat("     [✓] No strong single-break evidence at α.\n")
        cat("     • If conflict persists, consider multiple breaks or gradual regime changes.\n")
      }
    } else {
      cat(" [!] [8] BREAK PROTOCOL (Pettitt)\n")
      cat("     [!] trend package not installed → install.packages('trend') to enable break diagnostics.\n")
    }
    
    # [9] Conservative “default safe choice” rule (useful in teaching apps)
    cat("\n")
    cat(" [!] [9] DEFAULT SAFE CHOICE (when unsure)\n")
    cat("     • If ADF/KPSS conflict persists after fixing seasonality + lag:\n")
    cat("       - Prefer differencing (d=1) (and D=1 if seasonal) then re-test.\n")
    cat("     • This reduces the chance of building ARMA on a near-integrated series.\n")
    
    # [10] Sanity check transformations inside myData_Choice()
    cat("\n")
    cat(" [!] [10] TRANSFORMATION SANITY INSIDE myData_Choice()\n")
    cat("     • Ensure the same transformed object is returned for ALL downstream tests.\n")
    cat("     • Avoid mixing ts and numeric conversions before applying frequency-based operations.\n")
    cat("     • After log, verify positivity and handle zeros (e.g., log1p) if needed.\n")
    
    
    cat("==========================================================================\n\n")
    
    # --------------------------------------------------------------------------
    # EXTRA FINAL ACADEMIC ADVICE (more complete decision logic)
    # --------------------------------------------------------------------------
    
    # Power / sample size warnings (high impact)
    if (valid_N < 30) {
      cat(" [!] Power warning: with N < 30, ADF/KPSS can be unstable (low power / size distortions).\n")
      if (valid_N < 15) {
        cat("     [X] N < 15 is very weak for reliable inference; treat conclusions as provisional.\n")
      } else {
        cat("     [?] N is borderline; combine with plots + conservative transformations.\n")
      }
    } else {
      cat(" [✓] Sample size is generally adequate for classical stationarity tests.\n")
    }
    
    # Lag-vs-N safety
    if (valid_N <= (k + 10)) {
      cat(" [!] Lag risk: k is large relative to N (N <= k+10). Regression-based tests may misbehave.\n")
      cat("     Action: reduce k OR increase N, and prefer differencing/seasonal differencing over huge k.\n")
    } else {
      cat(" [✓] Lag order looks safe relative to N.\n")
    }
    
    # Model-type specification advice (none/drift/trend)
    cat(" [!] Model-type (none/drift/trend) matters:\n")
    cat("     • If the series has a visible non-zero mean but no deterministic trend → prefer 'drift'.\n")
    cat("     • If the series has a clear deterministic trend → prefer 'trend'.\n")
    cat("     • If the series oscillates around zero (rare in real data) → 'none'.\n")
    
    # Seasonality: warn if frequency exists but D not indicated
    if (is.finite(x_freq) && x_freq > 1 && !isTRUE(seasonality_resolved)) {
      cat(" [!] Seasonality risk: frequency suggests seasonality, but D=0 in UI.\n")
      cat("     Missing seasonal differencing can cause KPSS to reject stationarity and/or ADF to look borderline.\n")
    } else if (is.finite(x_freq) && x_freq > 1 && isTRUE(seasonality_resolved)) {
      cat(" [✓] Seasonality flag: D>0 in UI indicates seasonal treatment is intended.\n")
    }
    
    # Structural break contamination (Pettitt)
    if (requireNamespace("trend", quietly = TRUE) && is.finite(pett_p) && (pett_p < alpha_val)) {
      cat(" [!] Structural break contamination: Pettitt indicates a change-point (p < α).\n")
      cat("     ADF/KPSS disagreements are common under breaks; consider segment tests or break-aware models.\n")
    }
    
    # Near-unit-root / borderline zone heuristic
    # (We don't have direct 'borderline' threshold universally, so we use agreement + p-values)
    if (!isTRUE(agreement_safe) && is.finite(adf_p) && is.finite(eta_p_one)) {
      if (adf_p > 0.01 && adf_p < 0.10 && eta_p_one > 0.01 && eta_p_one < 0.10) {
        cat(" [?] Near-unit-root zone: both tests are near typical cutoffs.\n")
        cat("     Treat the series as highly persistent; prefer conservative differencing and validate by forecasting performance.\n")
      }
    }
    
    # Explosive alternative caveat (important correctness note)
    if (identical(as.character(alt_in), "explosive")) {
      cat(" [!] Explosive caveat: urca tau critical values are for the usual left-tail unit-root framework.\n")
      cat("     For explosive detection, rely primarily on tseries::adf.test p-value (right-tail) + plots.\n")
    }
    
    
    cat("==========================================================================\n\n")
    
    # Practical path
    cat(" PRACTICAL MODELING PATH:\n")
    if (isTRUE(is_stationary) && is.finite(lb_p) && (lb_p > alpha_val)) {
      cat(" [✓] Treat as I(0): identify ARMA on current series → fit → residual analysis.\n")
    } else if (!isTRUE(is_stationary)) {
      cat(" [X] Treat as I(1): apply differencing (d and/or D) → re-test → then identify ARMA.\n")
    } else {
      cat(" [?] Borderline/mixed: prefer conservative differencing and validate with residual diagnostics.\n")
    }
    
    cat("==========================================================================\n\n")
    
    
    
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # ====================================================================
  # ====================================================================
  # ====================================================================
  #
  #                          - HELP -
  #
  # ====================================================================
  # ====================================================================
  # ====================================================================
  
  
  
  
  # ---- Roadmap Detailed & teaching notes ----
  
  
  output$roadmap_Detailed_ui <- renderUI({
    tags$div(
      style = "background:#f7f7f7;padding:14px;border-radius:8px;",
      
      tags$hr(),
      
      tags$p(
        "Below is a ",
        tags$b("practical, SARIMA Modeling roadmap"),
      ),
      
      tags$hr(),
      
      tags$h4("[0] - Set the stage: define the modeling problem"),
      
      tags$h5("What students do"),
      tags$ul(
        tags$li("Define the ", tags$b("response series"), " (y_t) (what you forecast)."),
        tags$li("Define the ", tags$b("time index"), " (daily/weekly/monthly), and confirm it’s consistent."),
        tags$li(
          "Define the ", tags$b("forecast task"), ":",
          tags$ul(
            tags$li("horizon (e.g., 12 months ahead),"),
            tags$li("evaluation scheme (rolling-origin or simple train/test split),"),
            tags$li("loss metric (MAE/RMSE/MAPE/sMAPE).")
          )
        ),
        tags$li(
          "Decide whether you’ll model in:",
          tags$ul(
            tags$li(tags$b("levels"), " (raw data),"),
            tags$li(tags$b("log-levels"), " (common if variance grows with level),"),
            tags$li(tags$b("Box–Cox"), " transformed space (more general).")
          )
        )
      ),
      
      tags$h5("What they write (paper)"),
      tags$p(
        tags$b("Methods (Data & Objective). "),
        "“We modeled the univariate time series (y_t) observed at a [monthly] frequency from [start] to [end] (n=...). ",
        "The objective was to forecast (h=...) steps ahead. Model performance was evaluated using [metric(s)] under a ",
        "[train/test or rolling-origin] evaluation design.”"
      ),
      
      tags$h5("Pitfalls"),
      tags$ul(
        tags$li("SARIMA assumes ", tags$b("regular spacing"), "; irregular timestamps need fixing before anything else."),
        tags$li("SARIMA models ", tags$b("one series"), " (no predictors). If you have external regressors, that’s SARIMAX.")
      ),
      
      tags$hr(),
      
      tags$h4("[1] - Describe the data: sample size, missing values, descriptive statistics"),
      
      tags$h5("What students do"),
      tags$ol(
        tags$li(
          "Report:",
          tags$ul(
            tags$li("sample size (n),"),
            tags$li("start/end dates,"),
            tags$li("frequency,"),
            tags$li("number/percent missing values.")
          )
        ),
        tags$li(
          "Handle missingness:",
          tags$ul(
            tags$li("If rare and random: impute (linear interpolation, seasonal interpolation)."),
            tags$li("If many: reconsider the series, frequency, or data source.")
          )
        ),
        tags$li(
          "Descriptive stats:",
          tags$ul(
            tags$li("mean, median, sd, min/max,"),
            tags$li("maybe skewness/kurtosis,"),
            tags$li("and seasonal summaries (e.g., average by month).")
          )
        )
      ),
      
      tags$h5("What they write (APA-style)"),
      tags$p(
        tags$b("Results (Data description). "),
        "“The series contained (n=...) observations spanning [dates] at a [frequency] frequency. ",
        "Missing values accounted for (...%) of observations (k=... points). Missing observations were handled using ",
        "[method], selected because [reason]. The distribution of (y_t) showed a mean of (...) (SD=(...)), median (...), ",
        "and range ([...,...]).”"
      ),
      
      tags$h5("Pitfalls"),
      tags$ul(
        tags$li("Don’t “silently” impute—", tags$b("always justify"), " it."),
        tags$li("If you log-transform, describe the transformed series stats too.")
      ),
      
      tags$hr(),
      
      tags$h4("[2] - Explore visually: trend/seasonality/outliers; report observations"),
      
      tags$h5("What students do"),
      tags$p("Make plots and annotate:"),
      tags$ul(
        tags$li(tags$b("Line plot"), " of (y_t)."),
        tags$li(tags$b("Seasonal plot"), " (e.g., month-of-year lines)."),
        tags$li(tags$b("Boxplot by season"), " (month/quarter/week)."),
        tags$li(
          tags$b("Outlier check"), ":",
          tags$ul(
            tags$li("z-scores, IQR rule, or robust methods,"),
            tags$li("but also context (holidays, policy changes, measurement glitches).")
          )
        )
      ),
      
      tags$h5("What they write"),
      tags$p(
        tags$b("Results (Exploratory analysis). "),
        "“Visual inspection indicated [an upward/downward] trend and recurring seasonal fluctuations with period (s=...). ",
        "Variability appeared [constant/increasing with level], suggesting [no transformation / log transformation]. ",
        "Several potential outliers were observed around [dates], likely associated with [context], and were ",
        "[retained/adjusted] because [reason].”"
      ),
      
      tags$h5("Pitfalls"),
      tags$ul(
        tags$li("Outliers aren’t automatically “bad”—they might be real events that your forecast must respect."),
        tags$li("If variance grows with level, SARIMA often behaves better after a ", tags$b("log"), " or ", tags$b("Box–Cox"), " transform.")
      ),
      
      tags$hr(),
      
      tags$h4("[3] - Decompose: justify additive vs multiplicative; use STL when robust needed"),
      
      tags$h5("What students do"),
      tags$p("Perform decomposition to separate:"),
      tags$ul(
        tags$li("trend,"),
        tags$li("seasonality,"),
        tags$li("remainder.")
      ),
      
      tags$p(tags$b("Choose model form:")),
      tags$ul(
        tags$li(tags$b("Additive:"), " (y_t = T_t + S_t + e_t). Use when seasonal amplitude is roughly constant."),
        tags$li(tags$b("Multiplicative:"), " (y_t = T_t × S_t × e_t). Use when seasonal amplitude grows with the level (often solved by log transform → additive in log space).")
      ),
      
      tags$p(tags$b("Use STL decomposition when:")),
      tags$ul(
        tags$li("seasonality changes slowly over time,"),
        tags$li("you want robustness to outliers.")
      ),
      
      tags$h5("What they write"),
      tags$p(
        tags$b("Methods (Decomposition). "),
        "“We assessed additive versus multiplicative structure by examining whether seasonal amplitude scaled with the series level. ",
        "Because [seasonal variation was approximately constant / increased with level], we used an [additive model / log transformation] ",
        "and decomposed the series using [classical decomposition / STL]. STL was selected due to its robustness to outliers and its flexibility ",
        "in modeling evolving seasonality.”"
      ),
      
      tags$h5("Pitfalls"),
      tags$ul(
        tags$li("“Multiplicative seasonality” and “log transform” are basically best friends."),
        tags$li("STL decomposition is descriptive; SARIMA fitting still needs stationarity checks.")
      ),
      
      tags$hr(),
      
      tags$h4("[4] - Check stationarity: ADF/KPSS/PP; justify differencing (d and D)"),
      tags$p("SARIMA needs stationarity ", tags$b("after differencing"), "."),
      
      tags$h5("What students do"),
      tags$ol(
        tags$li("Define seasonal period (s) (e.g., 12 for monthly, 7 for daily-with-weekly seasonality)."),
        tags$li(
          "Test stationarity on:",
          tags$ul(
            tags$li("original series,"),
            tags$li("after ", tags$b("regular differencing"), " ((1-B)^d),"),
            tags$li("after ", tags$b("seasonal differencing"), " ((1-B^s)^D),"),
            tags$li("and sometimes after both.")
          )
        )
      ),
      
      tags$h5("Stationarity tests: what they do, H0/Ha, and how to conclude"),
      tags$ul(
        
        tags$li(
          tags$b("ADF test (Augmented Dickey–Fuller)"),
          tags$ul(
            tags$li(tags$b("What it does:"), " tests whether the series behaves like it has a ", tags$b("unit root"),
                    " (a stochastic trend), which implies non-stationarity; the test estimates a regression where lagged differences are added to handle autocorrelation."),
            tags$li(tags$b("H0:"), " the series has a unit root (non-stationary; shocks have permanent effects)."),
            tags$li(tags$b("Ha:"), " the series does not have a unit root (stationary around a mean or around a deterministic trend, depending on the ADF specification)."),
            tags$li(tags$b("Conclusion sentence template:"), " “The ADF test yielded p = [p-value]; therefore, at α = [alpha] we ",
                    tags$b("[reject / fail to reject]"), " H0 that the series contains a unit root. This implies the series is ",
                    tags$b("[stationary / non-stationary]"), " under the ADF framework, so we ",
                    tags$b("[did not apply additional differencing / applied]"), " [d=…] regular and/or [D=…] seasonal differencing to obtain an approximately stationary series suitable for SARIMA estimation.”")
          )
        ),
        
        tags$li(
          tags$b("KPSS test (Kwiatkowski–Phillips–Schmidt–Shin)"),
          tags$ul(
            tags$li(tags$b("What it does:"), " tests stationarity by examining whether the cumulative sum of residuals (from a level or trend regression) is too large; it is designed as a complement to ADF by flipping the null hypothesis."),
            tags$li(tags$b("H0:"), " the series is stationary (level-stationary, or trend-stationary if a trend is included)."),
            tags$li(tags$b("Ha:"), " the series is non-stationary (contains a unit root or otherwise violates stationarity)."),
            tags$li(tags$b("Conclusion sentence template:"), " “The KPSS test produced p = [p-value]; thus, at α = [alpha] we ",
                    tags$b("[reject / fail to reject]"), " H0 of stationarity. This indicates the series is ",
                    tags$b("[not stationary / consistent with stationarity]"), " in the KPSS sense, which ",
                    tags$b("[supports applying / does not require]"), " additional differencing; we therefore selected differencing orders [d=…] and [D=…] and rechecked stationarity on the transformed series.”")
          )
        ),
        
        tags$li(
          tags$b("PP test (Phillips–Perron)"),
          tags$ul(
            tags$li(tags$b("What it does:"), " tests for a unit root like ADF, but uses a nonparametric correction for autocorrelation and heteroskedasticity in the errors (instead of adding many lagged difference terms)."),
            tags$li(tags$b("H0:"), " the series has a unit root (non-stationary)."),
            tags$li(tags$b("Ha:"), " the series does not have a unit root (stationary)."),
            tags$li(tags$b("Conclusion sentence template:"), " “The PP test returned p = [p-value]; accordingly, at α = [alpha] we ",
                    tags$b("[reject / fail to reject]"), " H0 of a unit root. Interpreted alongside ADF and KPSS results, this suggests the series is ",
                    tags$b("[stationary / non-stationary]"), " after applying [d=…] regular and [D=…] seasonal differences, supporting the use of SARIMA on the differenced series.”")
          )
        )
      ),
      
      tags$p(
        tags$b("How to interpret ADF/KPSS/PP together (the logic students should write).")
      ),
      tags$ul(
        tags$li(tags$b("Best-case agreement:"), " ADF/PP reject unit root (small p) and KPSS fails to reject stationarity (large p) → strong evidence of stationarity."),
        tags$li(tags$b("Clear non-stationarity:"), " ADF/PP fail to reject unit root (large p) and KPSS rejects stationarity (small p) → strong evidence you need differencing."),
        tags$li(tags$b("Conflicts happen:"), " when tests disagree, prioritize the combination of evidence: plots + ACF behavior + results after differencing, and report that conclusions were based on the converging pattern rather than a single p-value.")
      ),
      
      tags$p(tags$b("Differencing logic:")),
      tags$ul(
        tags$li("Choose ", tags$b("d"), " to remove trend / unit root."),
        tags$li("Choose ", tags$b("D"), " to remove seasonal unit root."),
        tags$li("Stop as soon as stationarity is reasonable; avoid over-differencing.")
      ),
      
      tags$h5("What they write"),
      tags$p(
        tags$b("Methods (Stationarity and differencing). "),
        "“Stationarity was assessed using ADF, KPSS, and PP tests to triangulate evidence because the tests use different null hypotheses. ",
        "Based on the combined results and visual diagnostics, we selected [d=...] regular differences and [D=...] seasonal differences with seasonal period (s=...). ",
        "This differencing order was chosen to remove [trend/seasonal unit root] while avoiding over-differencing, and stationarity was re-evaluated on the transformed series before fitting SARIMA models.”"
      ),
      
      tags$h5("Pitfalls (classic)"),
      tags$ul(
        tags$li(
          tags$b("Over-differencing"),
          " causes:",
          tags$ul(
            tags$li("strong negative lag-1 autocorrelation,"),
            tags$li("inflated variance,"),
            tags$li("messier forecasts.")
          )
        ),
        tags$li("D is usually ", tags$b("0 or 1"), " in real life. If you need (D=2), the series may be weird or the seasonal period is wrong.")
      ),
      
      tags$hr(),
      
      tags$h4("[5] - Fit a baseline model: Auto-ARIMA to obtain a strong starting SARIMA"),
      
      tags$h5("What students do"),
      tags$ul(
        tags$li("Use an ", tags$b("auto-ARIMA"), " procedure (AICc-based or similar) to propose: ((p,d,q)(P,D,Q)_s)."),
        tags$li(
          "Keep track of:",
          tags$ul(
            tags$li("transformations used,"),
            tags$li("constraints (max p/q etc.),"),
            tags$li("whether stepwise search was used.")
          )
        ),
        tags$li(tags$b("Important:"), " Auto-ARIMA gives a baseline, not truth carved into granite.")
      ),
      
      tags$h5("What they write"),
      tags$p(
        tags$b("Methods (Baseline model). "),
        "“A baseline SARIMA model was selected using an automated information-criterion approach (minimizing AICc) over candidate orders ((p,q,P,Q)) ",
        "subject to [bounds]. The chosen baseline specification was SARIMA((p,d,q)(P,D,Q)_s), which served as the reference model for subsequent theory-driven refinement.”"
      ),
      
      tags$h5("Pitfalls"),
      tags$ul(
        tags$li("Auto-ARIMA can pick models that are statistically fine but ", tags$b("hard to interpret"), " or slightly unstable."),
        tags$li("If your evaluation is forecast-focused, it’s okay to prefer ", tags$b("simpler"), " models with similar accuracy.")
      ),
      
      tags$hr(),
      
      tags$h4("[6] - Fit a theory-driven model: Manual SARIMA using ACF/PACF + tests"),
      
      tags$h5("What students do"),
      tags$p("Using the differenced series (after chosen (d, D)):"),
      tags$ol(
        tags$li("Plot ", tags$b("ACF/PACF"), "."),
        tags$li(
          "Propose candidate structures:",
          tags$ul(
            tags$li(
              "Nonseasonal:",
              tags$ul(
                tags$li("AR(p): PACF cuts off around p; ACF tails."),
                tags$li("MA(q): ACF cuts off around q; PACF tails.")
              )
            ),
            tags$li(
              "Seasonal:",
              tags$ul(
                tags$li("seasonal AR(P): PACF spikes at lags (s, 2s, ...)."),
                tags$li("seasonal MA(Q): ACF spikes at lags (s, 2s, ...).")
              )
            )
          )
        ),
        tags$li("Fit a ", tags$b("small set"), " of plausible candidates (e.g., 3–8 models)."),
        tags$li(
          "Use tests/criteria:",
          tags$ul(
            tags$li("AICc/BIC,"),
            tags$li("parameter significance (with caution),"),
            tags$li("stability/invertibility checks.")
          )
        )
      ),
      
      tags$h5("What they write"),
      tags$p(
        tags$b("Methods (Theory-driven model building). "),
        "“Candidate SARIMA structures were proposed based on ACF/PACF behavior of the differenced series. Spikes at [lags] suggested nonseasonal [AR/MA] components, ",
        "while prominent autocorrelation at multiples of (s) indicated seasonal [AR/MA] terms. Several candidate models were fitted and compared using [AICc/BIC], with final selection also considering parsimony and diagnostic adequacy.”"
      ),
      
      tags$h5("Pitfalls"),
      tags$ul(
        tags$li("ACF/PACF heuristics are ", tags$b("guides"), ", not commandments."),
        tags$li("Don’t brute-force 200 models and pretend it’s “theory-driven.” Pick a ", tags$b("small, reasoned set"), ".")
      ),
      
      tags$hr(),
      
      tags$h4("[7] - Diagnose & compare: residual tests + forecast accuracy; choose final model"),
      
      tags$h5("What students do"),
      tags$p(tags$b("Residual diagnostics (must-do):")),
      tags$ul(
        tags$li("Residual time plot (should look like noise)."),
        tags$li("Residual ACF (no big spikes)."),
        tags$li(tags$b("Ljung–Box test"), " for residual autocorrelation."),
        tags$li("Normality checks (QQ plot; Shapiro-Wilk is too sensitive for big n)."),
        tags$li("Check heteroskedasticity (variance changing over time).")
      ),
      
      tags$p(tags$b("Forecast evaluation (must-do):")),
      tags$ul(
        tags$li("Holdout or rolling cross-validation."),
        tags$li("Metrics: MAE/RMSE; MAPE only if data never near zero."),
        tags$li(
          "Compare:",
          tags$ul(
            tags$li("baseline auto-ARIMA,"),
            tags$li("manual SARIMA candidates,"),
            tags$li("maybe a simple benchmark (seasonal naive).")
          )
        )
      ),
      
      tags$p(tags$b("Model choice rule (healthy):")),
      tags$ul(
        tags$li("Must pass diagnostics reasonably well."),
        tags$li("Must beat naive benchmark."),
        tags$li("Prefer simpler model if accuracy is essentially tied.")
      ),
      
      tags$h5("What they write"),
      tags$p(
        tags$b("Results (Model diagnostics and performance). "),
        "“Residual diagnostics indicated approximate white-noise behavior: residual autocorrelations were small and the Ljung–Box test was [non-significant/significant] at (α=...). ",
        "Forecast performance over the evaluation window showed MAE=(...) and RMSE=(...), outperforming the baseline and benchmark models. Based on diagnostic adequacy and predictive performance, the final selected model was SARIMA((p,d,q)(P,D,Q)_s).”"
      ),
      
      tags$h5("Pitfalls"),
      tags$ul(
        tags$li("A model with great AIC but autocorrelated residuals is basically ", tags$i("lying to you politely"), "."),
        tags$li("If residuals are non-normal, forecasts can still be good; the bigger issue is ", tags$b("autocorrelation"), " left in residuals.")
      ),
      
      tags$hr(),
      
      tags$h4("[8] - Write your paper: use APA paragraphs in each step; assemble Methods/Results"),
      
      tags$h5("What students do (assembly checklist)"),
      tags$p(tags$b("Methods section")),
      tags$ul(
        tags$li("Data (source, frequency, missing handling, transformation)."),
        tags$li("Exploratory approach (plots used, decomposition method)."),
        tags$li("Stationarity tests and differencing decisions."),
        tags$li("Baseline (auto-ARIMA settings)."),
        tags$li("Manual model selection rationale (ACF/PACF + candidates)."),
        tags$li("Diagnostics and evaluation scheme.")
      ),
      
      tags$p(tags$b("Results section")),
      tags$ul(
        tags$li("Data summary + key visual observations."),
        tags$li("Decomposition findings (trend/seasonality statements)."),
        tags$li("Stationarity test outcomes and chosen (d, D)."),
        tags$li("Final model parameters."),
        tags$li("Diagnostics results and accuracy results."),
        tags$li("Forecast plot + table of errors.")
      ),
      
      tags$h5("What they write (APA structure guidance)"),
      tags$ul(
        tags$li("Keep each subsection as: ", tags$b("What we did → Why → What we found → What we concluded"), "."),
        tags$li("Use past tense for Methods, results-oriented past tense for Results."),
        tags$li("Put the math in-line sparingly; put full model spec once, clearly.")
      ),
      
      tags$hr(),
      
      tags$h4("A clean “deliverable package” students should submit"),
      tags$ul(
        tags$li(
          "A notebook/script that:",
          tags$ul(
            tags$li("loads data,"),
            tags$li("handles missingness,"),
            tags$li("performs EDA plots,"),
            tags$li("decomposition,"),
            tags$li("stationarity tests,"),
            tags$li("baseline auto-ARIMA,"),
            tags$li("manual candidates,"),
            tags$li("diagnostics,"),
            tags$li("evaluation,"),
            tags$li("and final forecast.")
          )
        ),
        tags$li(
          "A short paper with:",
          tags$ul(
            tags$li("Methods + Results sections aligned to steps 1–7,"),
            tags$li("figures: time plot, decomposition, ACF/PACF, residual ACF, forecast plot,"),
            tags$li("a table comparing candidate models (AICc + metrics).")
          )
        )
      ),
      
      tags$hr(),
      
    )
  })
  
  
  
  
  
  
  
}
