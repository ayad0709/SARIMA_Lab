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
    s <- ts_train_test()
    x <- s$ts_train
    k <- as.numeric(input$adf_lags)
    list(
      adf = tryCatch(tseries::adf.test(x, k = k), error = function(e) NULL),
      kpss = tryCatch(tseries::kpss.test(x, null = input$kpss_type), error = function(e) NULL),
      pp = tryCatch(tseries::pp.test(x), error = function(e) NULL),
      ur = tryCatch(urca::ur.df(x, type = input$adf_type, lags = k), error = function(e) NULL)
    )
  })

  output$stationarity_results <- renderPrint({
    req(stationarity())
    st <- stationarity()
    cat("Stationarity tests (training series)\n\n")
    cat("ADF (tseries::adf.test)\n"); if (is.null(st$adf)) cat("  Not available.\n") else print(st$adf)
    cat("\nKPSS (tseries::kpss.test)\n"); if (is.null(st$kpss)) cat("  Not available.\n") else print(st$kpss)
    cat("\nPhillips–Perron (tseries::pp.test)\n"); if (is.null(st$pp)) cat("  Not available.\n") else print(st$pp)
    cat("\nADF (urca::ur.df) summary\n"); if (is.null(st$ur)) cat("  Not available.\n") else print(summary(st$ur))
  })

  output$stationarity_interpretation <- renderPrint({
    req(stationarity())
    st <- stationarity()
    cat("Interpretation (heuristic):\n\n")
    cat("- ADF: H0 = unit root (non-stationary). Small p suggests stationarity.\n")
    cat("- KPSS: H0 = stationary. Small p suggests non-stationary.\n")
    cat("- PP: H0 = unit root (non-stationary). Small p suggests stationarity.\n\n")
    if (!is.null(st$adf)) cat("ADF:", fmt_p(st$adf$p.value), "\n")
    if (!is.null(st$kpss)) cat("KPSS:", fmt_p(st$kpss$p.value), "\n")
    if (!is.null(st$pp)) cat("PP:", fmt_p(st$pp$p.value), "\n\n")
    cat("Recommendation: use combined evidence (tests + plots) to decide d and D.\n")
  })

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
  output$auto_diag_tests <- renderText({ req(auto_fit()); diag_tests_text(residuals(auto_fit()), lag = as.numeric(input$diag_lag), fitdf = length(coef(auto_fit()))) })

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
  output$manual_diag_tests <- renderText({ req(manual_fit()); diag_tests_text(residuals(manual_fit()), lag = as.numeric(input$diag_lag), fitdf = length(coef(manual_fit()))) })

  # --------------------------------------------- 
  # --------------------------------------------- 
  
  
  
  # ============================================================
  # --- MOD: Manual SARIMA equation renderer (FULL code) ---
  # ============================================================
  
  # --- MOD: Replace/Update your manual_equations() reactive with this FULL version ---
  manual_equations <- reactive({
    req(manual_fit(), ts_train_test())
    
    fit <- manual_fit()
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
    intercept_val <- if (length(intercept_name) > 0) unname(coefs[intercept_name[1]]) else 0
    show_intercept <- is.finite(intercept_val) && abs(intercept_val) > 1e-8
    intercept_num <- if (show_intercept) sprintf("%.3f", intercept_val) else ""
    
    # Drift
    drift_sym <- if (isTRUE(input$manual_drift)) " + \\delta t" else ""
    drift_val <- if (isTRUE(input$manual_drift) && "drift" %in% names(coefs)) unname(coefs["drift"]) else NA_real_
    drift_num <- if (isTRUE(input$manual_drift) && is.finite(drift_val) && abs(drift_val) > 1e-8) {
      paste0(" + ", sprintf("%.3f", drift_val), "t")
    } else if (isTRUE(input$manual_drift)) {
      " + \\delta t"
    } else ""
    
    # --- MOD: MathJax-safe equation blocks ---
    # Use \\[ ... \\] and align with aligned/array (supported).
    tex_display <- function(x) paste0("\\[", x, "\\]")
    
    # --- MOD: clean cosmetics in numeric equation line ---
    simplify_tex <- function(x) {
      x <- gsub("\\(1\\)", "", x)                      # remove (1)
      x <- gsub("\\s+", " ", x)                        # normalize spaces
      x <- gsub("\\+\\s*\\+", "+", x)                  # ++ -> +
      x <- gsub("\\+\\s*-", "-", x)                    # +- -> -
      x <- gsub("-\\s*\\+", "-", x)                    # -+ -> -
      x <- gsub("-\\s*-", "+", x)                      # -- -> +
      x <- gsub("\\s*\\+\\s*0\\.000\\b", "", x)        # remove + 0.000
      # x <- gsub("\\b0\\.000\\s*\\+\\s*", "", x)        # remove leading 0.000 +
      # x <- gsub("\\s*\\+\\s*0\\b", "", x)              # remove + 0
      # x <- gsub("\\b0\\s*\\+\\s*", "", x)              # remove leading 0 +
      trimws(x)
    }
    
    # Parameter-polynomial strings (line 3)
    poly_param_ar  <- function() if (p == 0) "1" else paste0("1", paste0(" - \\phi_{", ip, "}L^{", ip, "}", collapse = ""))
    poly_param_sar <- function() if (P == 0) "1" else paste0("1", paste0(" - \\Phi_{", iP, "}L^{", s * iP, "}", collapse = ""))
    poly_param_ma  <- function() if (q == 0) "1" else paste0("1", paste0(" + \\theta_{", iq, "}L^{", iq, "}", collapse = ""))
    poly_param_sma <- function() if (Q == 0) "1" else paste0("1", paste0(" + \\Theta_{", iQ, "}L^{", s * iQ, "}", collapse = ""))
    
    # Numeric-polynomial strings (line 4) with NA-safe filtering and correct signs
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
    
    # ---------- Line 1: General operator form ----------
    line1 <- paste0(
      "\\phi_p(L)\\,\\Phi_P(L^{S})\\,(1-L)^{d}(1-L^{S})^{D}Y_t",
      " = c + \\theta_q(L)\\,\\Theta_Q(L^{S})\\varepsilon_t", drift_sym
    )
    
    # ---------- Line 2: Expanded operator (summation) ----------
    line2 <- paste0(
      "\\left(1-\\sum_{i=1}^{p}\\phi_i L^{i}\\right)",
      "\\left(1-\\sum_{j=1}^{P}\\Phi_j L^{jS}\\right)",
      "(1-L)^{d}(1-L^{S})^{D}Y_t",
      " = c + ",
      "\\left(1+\\sum_{i=1}^{q}\\theta_i L^{i}\\right)",
      "\\left(1+\\sum_{j=1}^{Q}\\Theta_j L^{jS}\\right)",
      "\\varepsilon_t", drift_sym
    )
    
    # ---------- Line 3: parameter-expanded polynomials ----------
    line3 <- paste0(
      "(", poly_param_ar(), ")",
      "(", poly_param_sar(), ")",
      diff_part, sdiff_part,
      "Y_t = c + ",
      "(", poly_param_ma(), ")",
      "(", poly_param_sma(), ")\\varepsilon_t",
      drift_sym
    )
    
    # ---------- Line 4: numeric-expanded polynomials ----------
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
    
    # ---------- Equivalent time-domain (teaching form) ----------
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
    
    # --- MOD: coefficient listing, including ar/ma/sar/sma ---
    coef_lines <- c()
    
    if (show_intercept) coef_lines <- c(coef_lines, paste0("c (intercept/mean) = ", sprintf("%.4f", intercept_val)))
    if (isTRUE(input$manual_drift)) {
      if (is.finite(drift_val)) coef_lines <- c(coef_lines, paste0("drift = ", sprintf("%.4f", drift_val)))
      else coef_lines <- c(coef_lines, "drift included (value not estimated explicitly)")
    }
    
    if (p > 0) {
      for (i in ip) {
        nm <- paste0("ar", i)
        if (nm %in% names(coefs) && is.finite(coefs[nm])) {
          coef_lines <- c(coef_lines, paste0("ar", i, " (\\phi_", i, ") = ", sprintf("%.4f", coefs[nm])))
        }
      }
    }
    if (q > 0) {
      for (i in iq) {
        nm <- paste0("ma", i)
        if (nm %in% names(coefs) && is.finite(coefs[nm])) {
          coef_lines <- c(coef_lines, paste0("ma", i, " (\\theta_", i, ") = ", sprintf("%.4f", coefs[nm])))
        }
      }
    }
    if (P > 0) {
      for (i in iP) {
        nm <- paste0("sar", i)
        if (nm %in% names(coefs) && is.finite(coefs[nm])) {
          coef_lines <- c(coef_lines, paste0("sar", i, " (\\Phi_", i, ") = ", sprintf("%.4f", coefs[nm])))
        }
      }
    }
    if (Q > 0) {
      for (i in iQ) {
        nm <- paste0("sma", i)
        if (nm %in% names(coefs) && is.finite(coefs[nm])) {
          coef_lines <- c(coef_lines, paste0("sma", i, " (\\Theta_", i, ") = ", sprintf("%.4f", coefs[nm])))
        }
      }
    }
    
    if (length(coef_lines) == 0) coef_lines <- c("No coefficients available.")
    
    # --- MOD: package sections as MathJax display blocks ---
    list(
      p = p, d = d, q = q, P = P, D = D, Q = Q, s = s,
      coef_lines = coef_lines,
      eq_general = tex_display(line1),
      eq_expanded = tex_display(line2),
      eq_line3 = tex_display(line3),
      eq_line4 = tex_display(line4),
      eq_time_domain = tex_display(td)
    )
  })
  
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
        HTML(tex_display("\\text{------------}")),
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
  myData_Choice <- reactive({
    getMyData(
      tsData    = ts_base(),
      frequency = frequency(ts_base()),
      islog     = values$islog,
      d_n       = input$d_n,
      DS_n      = input$DS_n
    )
  })
  
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
    
    # ============================================================================
    # 0) SMALL HELPERS (safe input fallback + safe numeric)
    # ============================================================================
    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b
    to_num_safe <- function(v, default = NA_real_) {
      out <- suppressWarnings(as.numeric(v))
      if (length(out) == 0 || all(is.na(out)) || !is.finite(out[1])) default else out[1]
    }
    
    # ============================================================================
    # 1) INPUT COLLECTION (supports either naming convention in your UI)
    #    - This prevents "silent req stops" if your UI ids differ.
    # ============================================================================
    alt_in  <- input$alternd2St %||% input$alternSt
    lag_in  <- input$LagOrderADFd2St %||% input$LagOrderADFSt
    a_in    <- input$alphaSt2
    type_in <- input$adfTypeSt2
    
    # Required data
    req(myData_Choice())
    
    # Required inputs (validated safely)
    if (is.null(alt_in) || is.null(lag_in) || is.null(a_in) || is.null(type_in)) {
      cat("==========================================================================\n")
      cat("                STATE-OF-THE-ART ADF UNIT ROOT DIAGNOSTIC                 \n")
      cat("==========================================================================\n")
      cat(" [!] INPUT ERROR: One or more inputs are NULL.\n")
      cat("     This usually means a UI/server ID mismatch.\n")
      cat("     Needed inputs (either naming is OK):\n")
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
    x <- as.numeric(stats::na.omit(myData_Choice()))
    valid_N <- length(x)
    
    # lag (k) must be integer >= 0
    k <- suppressWarnings(as.integer(lag_in))
    if (length(k) == 0 || is.na(k) || k < 0) k <- 0L
    
    # alpha mapping (robust)
    alpha_raw <- as.character(a_in)
    alpha_val <- to_num_safe(alpha_raw, default = 0.05)
    
    alpha_col <- switch(alpha_raw,
                        "0.01" = "1pct",
                        "0.05" = "5pct",
                        "0.1"  = "10pct"
    )
    
    # tau row mapping
    tau_row <- switch(type_in,
                      "none"  = "tau1",
                      "drift" = "tau2",
                      "trend" = "tau3",
                      "tau3"
    )
    
    # ============================================================================
    # SECTION 1: DATA CONTEXT & SPECIFICATION
    # ============================================================================
    cat("==========================================================================\n")
    cat("                        ADF UNIT ROOT DIAGNOSTIC                          \n")
    cat("==========================================================================\n")
    cat(sprintf(" MODEL TYPE : %-10s | SAMPLE SIZE (N) : %d\n", toupper(type_in), valid_N))
    cat(sprintf(" LAG ORDER  : %-10d | SIGNIFICANCE (α) : %s\n", k, alpha_raw))
    cat(sprintf(" MOMENTS    : Mean: %.4f | Std.Dev: %.4f\n", mean(x), stats::sd(x)))
    # cat("--------------------------------------------------------------------------\n")
    cat("\n")
    # Basic sanity
    if (valid_N < 5) {
      cat(" [!] CRITICAL ERROR: Too few observations (N < 5).\n")
      cat("     DIRECTIVE: Provide more data points.\n")
      cat("==========================================================================\n")
      return(invisible(NULL))
    }
    
    if (!is.finite(stats::sd(x)) || stats::sd(x) == 0) {
      cat(" [!] CRITICAL ERROR: Series is constant or invalid (sd = 0 / NA).\n")
      cat("     DIRECTIVE: Check transformation (log/diff) or data quality.\n")
      cat("==========================================================================\n")
      return(invisible(NULL))
    }
    
    # Degrees of Freedom Safety Check (conservative)
    if (valid_N <= (k + 10)) {
      cat(" [!] WARNING: Sample size is small relative to selected lag.\n")
      cat("     This can break tseries::adf.test and weaken inference.\n")
      cat("     DIRECTIVE: Decrease Lag (k) or increase data points.\n")
      # cat("--------------------------------------------------------------------------\n")
      cat("\n")
    }
    
    # ============================================================================
    # 3) PREDECLARE OBJECTS (prevents crashes in post-summary)
    # ============================================================================
    res_urca <- NULL
    
    tau_obs  <- NA_real_
    tau_crit <- NA_real_
    
    res_tseries <- list(p.value = NA_real_)
    lb_test     <- list(statistic = NA_real_, p.value = NA_real_)
    kpss_type   <- if (type_in == "trend") "Trend" else "Level"
    res_kpss    <- list(statistic = NA_real_, p.value = NA_real_)
    
    is_stationary   <- FALSE
    kpss_stationary <- FALSE
    
    # Pettitt placeholders (ONLY — Buishand removed)
    pettitt_res <- list(statistic = NA_real_, p.value = NA_real_, estimate = NULL)
    
    # ============================================================================
    # 4) COMPUTATIONS (tryCatch keeps UI alive, cats remain)
    # ============================================================================
    tryCatch({
      
      # --- Package checks (explicit + friendly) ---
      if (!requireNamespace("urca", quietly = TRUE)) {
        stop("Package 'urca' is not installed. ACTION: install.packages('urca')")
      }
      if (!requireNamespace("tseries", quietly = TRUE)) {
        stop("Package 'tseries' is not installed. ACTION: install.packages('tseries')")
      }
      
      # ==========================================================================
      # PRIMARY TEST: ADF via urca::ur.df
      # ==========================================================================
      res_urca <- urca::ur.df(x, type = type_in, lags = k)
      
      # Tau statistic: extract by mapped row safely
      tau_obs <- suppressWarnings(as.numeric(res_urca@teststat[tau_row]))
      if (!is.finite(tau_obs)) {
        # fallback attempt (some objects behave differently in edge cases)
        tau_obs_fallback <- suppressWarnings(as.numeric(res_urca@teststat[1]))
        if (is.finite(tau_obs_fallback)) tau_obs <- tau_obs_fallback
      }
      
      # Critical value
      tau_crit <- suppressWarnings(as.numeric(res_urca@cval[tau_row, alpha_col]))
      
      # Secondary ADF (p-value reference) — only if feasible
      if (valid_N > (k + 10)) {
        res_tseries <- tseries::adf.test(x, alternative = alt_in, k = k)
      } else {
        res_tseries <- list(p.value = NA_real_)
      }
      
      # ==========================================================================
      # QUALITY CHECK: Ljung-Box on ur.df residuals
      # ==========================================================================
      lb_lag <- max(1L, min(10L, floor(valid_N / 5)))
      lb_test <- Box.test(res_urca@res, lag = lb_lag, type = "Ljung-Box")
      
      # ==========================================================================
      # CONFIRMATORY: KPSS (H0 = stationary)
      # ==========================================================================
      kpss_type <- if (type_in == "trend") "Trend" else "Level"
      res_kpss  <- tseries::kpss.test(x, null = kpss_type)
      
      # ==========================================================================
      # SAFE FLAGS (avoid NA inside if-conditions)
      # ==========================================================================
      tau_ok  <- is.finite(tau_obs) && is.finite(tau_crit)
      lb_ok   <- is.finite(to_num_safe(lb_test$p.value))
      kpss_ok <- is.finite(to_num_safe(res_kpss$p.value))
      
      is_stationary <- isTRUE(tau_ok) && (tau_obs < tau_crit)
      
      if (kpss_ok && (to_num_safe(res_kpss$p.value) > alpha_val)) {
        kpss_stationary <- TRUE
      } else {
        kpss_stationary <- FALSE
      }
      
      lb_white_safe <- lb_ok && (to_num_safe(lb_test$p.value) > alpha_val)
      
      # ==========================================================================
      # SECTION 2: PRIMARY TEST - AUGMENTED DICKEY-FULLER (ADF)
      # ==========================================================================
      cat("==========================================================================\n")
      cat("PHASE 1: ADF UNIT ROOT TEST\n")
      cat("==========================================================================\n")
      
      cat(" • H0: The series has a Unit Root (Non-Stationary).\n")
      cat(" • Ha: The series is Stationary (Mean Reverting).\n")
      cat(sprintf(" -> CRITERIA: Reject H0 if Tau-Obs (%.4f) < Tau-Crit (%.4f)\n",
                  tau_obs, tau_crit))
      
      cat("\n RESULT:\n")
      cat(paste("  - Tau Observed :", round(tau_obs, 4), "\n"))
      cat(paste("  - Tau Critical :", round(tau_crit, 4), "\n"))
      cat(paste("  - P-Value (Ref):", round(to_num_safe(res_tseries$p.value), 4), "\n"))
      
      cat("\n DECISION:\n")
      if (isTRUE(tau_ok) && isTRUE(is_stationary)) {
        cat("  -> REJECT H0: Evidence suggests the series is STATIONARY.\n")
      } else if (isTRUE(tau_ok) && !isTRUE(is_stationary)) {
        cat("  -> FAIL TO REJECT H0: Evidence suggests the series is NON-STATIONARY.\n")
      } else {
        cat("  -> [!] WARNING: Tau/Critical value is NA/Inf. ADF inference is not valid.\n")
      }
      # cat("--------------------------------------------------------------------------\n")
      cat("\n")
      
      # ==========================================================================
      # SECTION 3: QUALITY CHECK - LJUNG-BOX TEST
      # ==========================================================================
      cat("==========================================================================\n")
      cat("PHASE 2: RESIDUAL DIAGNOSTICS (LJUNG-BOX)\n")
      cat("==========================================================================\n")
      
      cat(" • H0: Residuals are White Noise (No Autocorrelation).\n")
      cat(" • Ha: Residuals are Correlated (Lags are insufficient).\n")
      cat(sprintf(" -> CRITERIA: Reject H0 if P-Value (%.4f) < α (%.2f)\n",
                  to_num_safe(lb_test$p.value), alpha_val))
      
      cat("\n RESULT:\n")
      cat(paste("  - LB Statistic :", round(to_num_safe(lb_test$statistic), 4), "\n"))
      cat(paste("  - LB P-Value   :", round(to_num_safe(lb_test$p.value), 4), "\n"))
      
      cat("\n DECISION:\n")
      if (isTRUE(lb_white_safe)) {
        cat("  -> FAIL TO REJECT H0: Residuals are White Noise. [ADF VALID]\n")
      } else if (isTRUE(lb_ok)) {
        cat("  -> REJECT H0: Residuals are Correlated. [ADF BIASED]\n")
      } else {
        cat("  -> [!] WARNING: Ljung-Box P-Value is NA/Inf. Residual diagnosis failed.\n")
      }
      # cat("--------------------------------------------------------------------------\n")
      cat("\n")
      
      
      
      # ==========================================================================
      # PHASE 3: KPSS + STRUCTURAL BREAK CHECK (Pettitt) + Segment KPSS
      # ==========================================================================
      
      cat("==========================================================================\n")
      cat("PHASE 3: CONFIRMATORY ANALYSIS (KPSS)+ STRUCTURAL BREAK CHECK (Pettitt) + Segment KPSS\n")
      cat("==========================================================================\n\n")
      cat("PHASE 3A: CONFIRMATORY ANALYSIS (KPSS)\n")
      cat(paste(" • H0: The series is Stationary around a", kpss_type, ".\n"))
      cat(" • Ha: The series is Non-Stationary (Unit Root).\n")
      cat(sprintf(" • CRITERIA: Reject H0 if P-Value (%.4f) < α (%.2f)\n",
                  to_num_safe(res_kpss$p.value), alpha_val))
      
      cat("\n RESULT:\n")
      cat(paste("  - KPSS P-Value :", round(to_num_safe(res_kpss$p.value), 4), "\n"))
      
      cat("\n DECISION:\n")
      if (isTRUE(kpss_ok) && isTRUE(kpss_stationary)) {
        cat("  -> FAIL TO REJECT H0: The series is STATIONARY. [KPSS]\n")
      } else if (isTRUE(kpss_ok) && !isTRUE(kpss_stationary)) {
        cat("  -> REJECT H0: The series is NON-STATIONARY. [KPSS]\n")
      } else {
        cat("  -> [!] WARNING: KPSS P-Value is NA/Inf. KPSS inference is not valid.\n")
      }
      cat("--------------------------------------------------------------------------\n")
      
      # ==========================================================================
      # PHASE 3B: STRUCTURAL BREAK CHECK (Pettitt) + KPSS re-check by segments
      # ==========================================================================
      cat("PHASE 3B: STRUCTURAL BREAK CHECK (Pettitt) + KPSS SEGMENT CHECK\n")
      cat(" • Goal: Detect a single change-point (shift) that can contaminate KPSS.\n")
      cat(" • If a break exists, we re-run KPSS before/after the break.\n")
      
      pettitt_break <- FALSE
      break_idx <- NA_integer_
      
      if (requireNamespace("trend", quietly = TRUE)) {
        
        # Pettitt safely (never kill the rest of the report)
        tryCatch({
          pettitt_res <- trend::pettitt.test(x)
        }, error = function(e) {
          pettitt_res <- list(statistic = NA_real_, p.value = NA_real_, estimate = NULL)
          cat(" [!] WARNING: Pettitt test failed.\n")
          cat("     ERROR: ", e$message, "\n", sep = "")
        })
        
        pett_p <- to_num_safe(pettitt_res$p.value)
        cat(sprintf(" • Pettitt P-Value: %.6f | α: %.4f\n", pett_p, alpha_val))
        
        # estimate is often the location (index)
        if (!is.null(pettitt_res$estimate) && is.finite(as.numeric(pettitt_res$estimate))) {
          break_idx <- as.integer(as.numeric(pettitt_res$estimate))
        }
        
        pettitt_break <- is.finite(pett_p) && (pett_p < alpha_val) && is.finite(break_idx)
        
        cat("\n RESULT:\n")
        cat(paste("  - Break index (estimate):", ifelse(is.finite(break_idx), break_idx, "NA"), "\n"))
        cat(paste("  - Pettitt p-value       :", round(pett_p, 6), "\n"))
        
        cat("\n DECISION:\n")
        if (pettitt_break) {
          cat("  -> REJECT H0: Break detected. KPSS may be distorted on full sample.\n")
        } else {
          cat("  -> FAIL TO REJECT H0: No strong evidence of a single break.\n")
        }
        
        # ------------------------------------------------------------
        # Segment KPSS (only if break is detected and segments are big enough)
        # ------------------------------------------------------------
        if (pettitt_break) {
          
          # Split around break_idx
          x1 <- x[1:break_idx]
          x2 <- x[(break_idx + 1):valid_N]
          
          cat("--------------------------------------------------------------------------\n")
          
          cat("PHASE 3C: KPSS RE-CHECK BY SEGMENTS:\n")
          
          cat("  • Segment 1 = [1 .. break]\n")
          cat("  • Segment 2 = [break+1 .. N]\n")
          
          # Minimum length guard (KPSS needs enough points)
          min_seg <- 12L
          if (length(x1) < min_seg || length(x2) < min_seg) {
            cat("  [!] WARNING: One segment is too short for reliable KPSS.\n")
            cat(sprintf("     Segment lengths: n1=%d, n2=%d (min recommended=%d)\n",
                        length(x1), length(x2), min_seg))
          } else {
            
            # Run KPSS on each segment safely
            kpss1 <- tryCatch(tseries::kpss.test(x1, null = kpss_type),
                              error = function(e) list(p.value = NA_real_))
            kpss2 <- tryCatch(tseries::kpss.test(x2, null = kpss_type),
                              error = function(e) list(p.value = NA_real_))
            
            p1 <- to_num_safe(kpss1$p.value)
            p2 <- to_num_safe(kpss2$p.value)
            
            cat(sprintf("  - KPSS p-value (Segment 1): %.6f\n", p1))
            cat(sprintf("  - KPSS p-value (Segment 2): %.6f\n", p2))
            
            # Interpret segment KPSS
            s1_stat <- is.finite(p1) && (p1 > alpha_val)
            s2_stat <- is.finite(p2) && (p2 > alpha_val)
            
            cat("\n INTERPRETATION:\n")
            if (s1_stat && s2_stat) {
              cat("  [✓] Both segments look stationary by KPSS.\n")
              cat("      -> Full-sample KPSS non-stationarity may be break-driven.\n")
            } else if (!s1_stat && !s2_stat) {
              cat("  [X] Both segments look non-stationary by KPSS.\n")
              cat("      -> Non-stationarity is not only due to a break.\n")
            } else {
              cat("  [?] Mixed: one segment stationary, the other not.\n")
              cat("      -> Consider regime modeling or re-check transformation.\n")
            }
          }
        }
        
      } else {
        cat(" [!] WARNING: Package 'trend' not installed → Pettitt break check skipped.\n")
        cat("     ACTION: install.packages('trend')\n")
      }
      
      
      
      # cat("--------------------------------------------------------------------------\n")
      
      
      
      
      # ==========================================================================
      # SECTION 5: FINAL VERDICT & STUDENT DIRECTIVES
      # ==========================================================================
      cat("==========================================================================\n")
      cat("PHASE 4: FINAL ACADEMIC VERDICT & ADVICE\n")
      cat("==========================================================================\n")
      
      
      if (isTRUE(lb_ok) && !isTRUE(lb_white_safe)) {
        cat(" [!] WARNING: Your ADF 'Tau' statistic is technically biased.\n")
        cat("     ADVICE: Increase 'Lag Order' until Ljung-Box P-Value > α.\n\n")
      }
      
      if (isTRUE(is_stationary) && isTRUE(kpss_stationary)) {
        cat(" [✓] VERDICT: STRONG STATIONARITY confirmed by both ADF and KPSS.\n")
        cat("     ACTION: Proceed to ARMA/ARIMA modeling with d=0.\n")
      } else if (isTRUE(!is_stationary) && isTRUE(!kpss_stationary)) {
        cat(" [X] VERDICT: CLEAR UNIT ROOT. Both tests confirm Non-Stationarity.\n")
        cat("     ACTION: Apply differencing (d=1) to stabilize the mean.\n")
      } else {
        cat(" [?] VERDICT: CONFLICTING RESULTS (ADF vs KPSS).\n")
        cat("     ADVICE: The series may be highly persistent or trend-stationary.\n")
        cat("     Double-check plots for breaks or changing variance.\n")
      }
      
      cat("\n TECHNICAL APPENDIX (ADF Regression Coefficients):\n")
      if (!is.null(res_urca)) {
        print(stats::coef(res_urca@testreg))
      } else {
        cat("  [!] Not available (ur.df did not run).\n")
      }
      
    }, error = function(e) {
      cat(" EXECUTION ERROR: ", e$message, "\n")
    })
    
    # cat("--------------------------------------------------------------------------\n")
    cat("\n")
    
    # ============================================================================
    # 5) POST-SUMMARY (always runs; safe against NA objects)
    # ============================================================================
    cat("==========================================================================\n")
    cat("PHASE 5: POST-SUMMARY (Academic-quality snapshot)\n")
    cat("==========================================================================\n")
    
    cat(" EVIDENCE SNAPSHOT (All key outcomes in one place):\n")
    cat(sprintf(" [ ] N (effective sample size)        : %d\n", valid_N))
    cat(sprintf(" [ ] Model type (ADF)                 : %s  (tau row: %s)\n", type_in, tau_row))
    cat(sprintf(" [ ] Lag order (k)                    : %d\n", k))
    cat(sprintf(" [ ] Alpha (α)                        : %.4f\n", alpha_val))
    cat(sprintf(" [ ] Tau-Observed (urca)              : %.6f\n", tau_obs))
    cat(sprintf(" [ ] Tau-Critical (urca, %s)          : %.6f\n", alpha_col, tau_crit))
    cat(sprintf(" [ ] ADF p-value (tseries reference)  : %.6f\n", to_num_safe(res_tseries$p.value)))
    cat(sprintf(" [ ] Ljung-Box p-value (residuals)    : %.6f\n", to_num_safe(lb_test$p.value)))
    cat(sprintf(" [ ] KPSS p-value (%s null)           : %.6f\n", kpss_type, to_num_safe(res_kpss$p.value)))
    cat("--------------------------------------------------------------------------\n")
    
    cat(" CHECKLIST (Academic-quality acceptance criteria):\n")
    cat(sprintf(" [ ] Tau is finite and usable                 : %s\n",
                ifelse(is.finite(tau_obs), "[✓] SATISFIED", "[!] NOT SATISFIED")))
    if (!is.finite(tau_obs)) {
      cat("     [!] DANGER: Tau is NA/Inf → ADF inference cannot be trusted.\n")
      cat("     [!] Likely causes: constant/near-constant series, k too high,\n")
      cat("         wrong extraction index, too small effective N.\n")
    }
    
    cat(sprintf(" [ ] Tau critical value extracted             : %s\n",
                ifelse(is.finite(tau_crit), "[✓] SATISFIED", "[!] NOT SATISFIED")))
    if (!is.finite(tau_crit)) {
      cat("     [!] DANGER: Critical value missing → alpha/type mapping issue.\n")
      cat("     [!] ACTION: verify alpha_col and tau_row exist in urca cval table.\n")
    }
    
    lb_p <- to_num_safe(lb_test$p.value)
    cat(sprintf(" [ ] Residuals pass Ljung-Box (white noise)   : %s\n",
                ifelse(is.finite(lb_p) && lb_p > alpha_val, "[✓] SATISFIED", "[X] NOT SATISFIED")))
    if (is.finite(lb_p) && lb_p <= alpha_val) {
      cat("     [!] WARNING: residual autocorrelation detected.\n")
      cat("     [!] Interpretation: ADF may be biased (insufficient k).\n")
      cat("     [!] Suggested fix: increase k cautiously OR re-check transformation.\n")
    }
    
    cat(sprintf(" [ ] Sample size adequacy                     : %s\n",
                ifelse(valid_N >= 30, "[✓] STRONG", ifelse(valid_N >= 15, "[?] BORDERLINE", "[X] WEAK"))))
    if (valid_N < 15) {
      cat("     [!] WARNING: very low power. ADF/KPSS decisions may be unstable.\n")
      cat("     [!] Advice: combine with plots (time plot + ACF/PACF) and be conservative.\n")
    }
    
    
    
    #  lag sanity
    cat(sprintf(" [ ] Lag order reasonable relative to N       : %s\n",
                ifelse(valid_N > (k + 5), "[✓] OK", "[!] TOO HIGH / RISKY")))
    if (valid_N <= (k + 5)) {
      cat("     [!] DANGER: k too high for N → regression can break or give NA stats.\n")
      cat("     [!] ACTION: reduce k immediately.\n")
    }
    
    
    #  cross-test agreement (ADF vs KPSS)
    agreement_safe <- (isTRUE(is_stationary) && isTRUE(kpss_stationary)) ||
      (isTRUE(!is_stationary) && isTRUE(!kpss_stationary))
    
    cat(sprintf(" [ ] ADF & KPSS agreement                      : %s\n",
                ifelse(agreement_safe, "[✓] AGREEMENT", "[?] CONFLICT")))
    
    if (!agreement_safe) {
      cat("     [?] NOTE: conflict can happen for near-unit-root series, structural breaks,\n")
      cat("         or when model type (none/drift/trend) is mis-specified.\n")
    }
    
    
    # Structural break summary (Pettitt only)
    if (requireNamespace("trend", quietly = TRUE)) {
      cat("--------------------------------------------------------------------------\n")
      cat(" STRUCTURAL BREAK SUMMARY (Pettitt):\n")
      cat(sprintf(" [ ] Pettitt p-value   : %.6f  (Reject H0 if < α)\n", to_num_safe(pettitt_res$p.value)))
      
      pett_break2 <- is.finite(to_num_safe(pettitt_res$p.value)) &&
        (to_num_safe(pettitt_res$p.value) < alpha_val)
      
      cat("\n INTERPRETATION:\n")
      if (pett_break2) {
        cat(" [!] WARNING: Evidence of a break detected (Pettitt).\n")
        cat("     ADVICE: Consider split-sample testing or break-aware modeling.\n")
        cat("     NOTE: Breaks can cause ADF/KPSS instability or conflicts.\n")
      } else {
        cat(" [✓] No strong evidence of a major break at α (Pettitt).\n")
        cat("     ADVICE: Stationarity conclusions are less likely to be break-contaminated.\n")
      }
    } else {
      cat("--------------------------------------------------------------------------\n")
      cat(" STRUCTURAL BREAK SUMMARY (Pettitt):\n")
      cat(" [!] WARNING: Package 'trend' is not installed. Pettitt summary unavailable.\n")
      cat("     ACTION: install.packages('trend')\n")
    }
    
    cat("--------------------------------------------------------------------------\n")
    
    if (isTRUE(is_stationary) && isTRUE(kpss_stationary)) {
      cat(" [✓] VERDICT: CONSISTENT STATIONARITY. Both tests agree the series is I(0).\n")
      cat(" ADVICE: You may proceed to model identification using the original series.\n")
    } else if (isTRUE(!is_stationary) && isTRUE(!kpss_stationary)) {
      cat(" [X] VERDICT: CONSISTENT UNIT ROOT. Both tests confirm Non-Stationarity.\n")
      cat(" ADVICE: Apply First Differencing (d=1) to achieve stationarity.\n")
    } else {
      cat(" [?] VERDICT: CONFLICTING RESULTS. One test suggests a Unit Root while the other doesn't.\n")
      cat(" ADVICE: This often happens with near-integrated or break-contaminated series. \n")
      cat("         Use differencing for safety.\n")
    }
    
    cat("\n ACTIONABLE NEXT STEPS (What to do now):\n")
    
    if (!is.finite(tau_obs) || !is.finite(tau_crit)) {
      cat(" [!] PRIORITY 1: Fix the ADF computation first.\n")
      cat("     • Reduce k, verify variance is not zero, and confirm urca outputs.\n")
      cat("     • If Tau remains NA, your series may be too short after transformation.\n")
    } else {
      cat(" [✓] PRIORITY 1: Tau and critical values are usable.\n")
    }
    
    if (is.finite(lb_p) && lb_p <= alpha_val) {
      cat(" [X] PRIORITY 2: Residual autocorrelation detected.\n")
      cat("     • Increase the Lag k (carefully) OR re-check transformation (d/D/log).\n")
    } else {
      cat(" [✓] PRIORITY 2: Residual whiteness acceptable → ADF inference more reliable.\n")
    }
    
    if (!agreement_safe) {
      cat(" [?] PRIORITY 3: Resolve ADF vs KPSS conflict.\n")
      cat("     • Try alternative model types (none/drift/trend).\n")
      cat("     • Consider break-aware modeling if Pettitt flags a break.\n")
      cat("     • Use visual diagnostics + conservative differencing.\n")
    } else {
      cat(" [✓] PRIORITY 3: Tests agree → proceed with confidence.\n")
    }
    
    cat("\n PRACTICAL MODELING PATH (for your Shiny workflow):\n")
    if (isTRUE(is_stationary) && is.finite(lb_p) && (lb_p > alpha_val)) {
      cat(" [✓] Use ARMA identification on current series → fit → residual analysis.\n")
    } else if (!isTRUE(is_stationary)) {
      cat(" [X] Apply differencing (d and/or D) → re-run ADF/KPSS → then identify ARMA.\n")
    } else {
      cat(" [?] Mixed evidence: consider differencing for safety and validate visually.\n")
    }
    
    cat("--------------------------------------------------------------------------\n\n")
  })
  
  
  
  
  
  
  
  
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
