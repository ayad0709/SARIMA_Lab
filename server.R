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
  
  
  # prepared <- reactive({
  #   req(raw_data(), input$dateCol, input$valueCol)
  #   f <- freq_value(input)
  #   by <- freq_to_by(f)
  # 
  #   df <- raw_data()
  #   d <- parse_dates(df[[input$dateCol]])
  #   y <- suppressWarnings(as.numeric(df[[input$valueCol]]))
  # 
  #   keep <- !is.na(d)
  #   df2 <- data.frame(date = as.Date(d[keep]), y_raw = y[keep])
  #   df2 <- df2[order(df2$date), , drop = FALSE]
  # 
  #   if (isTRUE(input$align_regular) && !is.null(by)) {
  #     grid <- make_regular_grid(df2$date, by = by)
  #     df2 <- merge(data.frame(date = grid), df2, by = "date", all.x = TRUE, sort = TRUE)
  #   }
  # 
  #   df2$y_filled <- fill_missing(df2$y_raw, input$missing_policy, f)
  # 
  #   df2$y_trans <- tryCatch(
  #     apply_transform(df2$y_filled, input$transform, input$lambda),
  #     error = function(e) { validate(e$message); df2$y_filled }
  #   )
  # 
  #   if (!is.null(by)) {
  #     df2$x <- df2$date
  #     x_label <- "Date"
  #   } else {
  #     df2$x <- seq_len(nrow(df2))
  #     x_label <- "Index"
  #   }
  # 
  #   list(df = df2, freq = f, by = by, x_label = x_label)
  # })

  ts_train_test <- reactive({
    p <- prepared()
    df <- p$df

    ok <- is.finite(df$y_trans)
    dfm <- df[ok, , drop = FALSE]
    validate(need(nrow(dfm) >= 10, "Not enough valid observations after cleaning."))

    train_n <- max(2, floor(nrow(dfm) * as.numeric(input$train_prop)))
    y_tr <- dfm$y_trans[seq_len(train_n)]
    y_te <- if (train_n < nrow(dfm)) dfm$y_trans[(train_n + 1):nrow(dfm)] else numeric(0)

    list(
      dfm = dfm,
      train_n = train_n,
      test_n = length(y_te),
      ts_train = ts(y_tr, start = 1, frequency = p$freq),
      ts_test = if (length(y_te) > 0) ts(y_te, start = train_n + 1, frequency = p$freq) else ts(numeric(0), frequency = p$freq)
    )
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
  #   Fixes:
  #   1) Add AR/MA/SAR/SMA coefficient listings (as requested)
  #   2) Fix MathJax rendering: remove unsupported LaTeX env "flushleft"
  #      (MathJax often prints it as raw text). Use a MathJax-safe block:
  #      \\[ ... \\] with aligned/array, and LEFT-ALIGN via HTML/CSS wrapper.
  #   3) Print sections in this order:
  #      - Manual SARIMA model + coefficients (including drift if any)
  #      - "General SARIMA formulation" -> equation line 1
  #      - "Expanded operator form"     -> equation line 2
  #      - "Numerical model (one-line)" -> equation lines 3 and 4
  #      - "Equivalent time-domain representation (teaching form)" -> equation
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
  
  
  
  
  
  
  
  
  # # ============================================================
  # # --- MOD: Fix 4th-line cleanup + remove NA terms in numeric polynomials ---
  # #   Fixes requested:
  # #   1) Remove useless "(1)" factors (multiplying by 1)
  # #   2) Remove "0.000" intercept when ~0 (adding 0)
  # #   3) Prevent NA coefficients / NA lag powers from appearing (e.g., "1NAL12")
  # # ============================================================
  # 
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
  #   # --- MOD: safer seasonal period (never allow NA) ---
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
  #   # --- MOD: intercept shown only if meaningfully non-zero ---
  #   intercept_name <- intersect(c("intercept", "mean"), names(coefs))
  #   intercept_val <- if (length(intercept_name) > 0) unname(coefs[intercept_name[1]]) else 0
  #   show_intercept <- is.finite(intercept_val) && abs(intercept_val) > 1e-8
  #   intercept_num <- if (show_intercept) sprintf("%.3f", intercept_val) else ""
  # 
  #   # Drift flag (symbolic vs numeric)
  #   drift_sym <- if (isTRUE(input$manual_drift)) " + \\delta t" else ""
  #   drift_num <- if (isTRUE(input$manual_drift) && "drift" %in% names(coefs) && is.finite(coefs["drift"])) {
  #     paste0(" + ", sprintf("%.3f", unname(coefs["drift"])), "t")
  #   } else if (isTRUE(input$manual_drift)) {
  #     " + \\delta t"
  #   } else ""
  # 
  #   # Helpers
  #   tex_block <- function(x) paste0("\\[", x, "\\]")
  # 
  #   # --- MOD: helper to safely fetch coefficients and drop NA values ---
  #   get_coefs_safe <- function(nms) {
  #     if (length(nms) == 0) return(numeric(0))
  #     v <- suppressWarnings(unname(coefs[nms]))
  #     v <- v[is.finite(v)]  # drop NA/Inf
  #     v
  #   }
  # 
  #   # --- MOD: Build numeric polynomials without NA and with correct sign logic ---
  #   # AR/SAR polynomials: (1 - ar_i L^i) => we add term with coefficient = -ar_i
  #   # MA/SMA polynomials: (1 + ma_i L^i) => we add term with coefficient = +ma_i
  #   poly_num_ar <- function() {
  #     if (p == 0) return("1")
  #     nms <- paste0("ar", ip)
  #     v <- suppressWarnings(unname(coefs[nms]))
  #     keep <- is.finite(v)
  #     if (!any(keep)) return("1")  # MOD: if missing, treat as 1
  #     paste0("1", paste(sprintf(" %+.3fL^{%d}", -v[keep], ip[keep]), collapse = ""))
  #   }
  #   poly_num_sar <- function() {
  #     if (P == 0) return("1")
  #     nms <- paste0("sar", iP)
  #     v <- suppressWarnings(unname(coefs[nms]))
  #     keep <- is.finite(v)
  #     if (!any(keep)) return("1")  # MOD: if missing, treat as 1
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
  #   # Parameter polynomials (for line 3)
  #   poly_param_ar <- function() {
  #     if (p == 0) return("1")
  #     paste0("1", paste0(" - \\phi_{", ip, "}L^{", ip, "}", collapse = ""))
  #   }
  #   poly_param_sar <- function() {
  #     if (P == 0) return("1")
  #     paste0("1", paste0(" - \\Phi_{", iP, "}L^{", s * iP, "}", collapse = ""))
  #   }
  #   poly_param_ma <- function() {
  #     if (q == 0) return("1")
  #     paste0("1", paste0(" + \\theta_{", iq, "}L^{", iq, "}", collapse = ""))
  #   }
  #   poly_param_sma <- function() {
  #     if (Q == 0) return("1")
  #     paste0("1", paste0(" + \\Theta_{", iQ, "}L^{", s * iQ, "}", collapse = ""))
  #   }
  # 
  #   # Differencing operators
  #   diff_part  <- if (d > 0) paste0("(1-L)^{", d, "}") else ""
  #   sdiff_part <- if (D > 0) paste0("(1-L^{", s, "})^{", D, "}") else ""
  # 
  #   # --- MOD: cleanup function for line 4 (remove (1), remove +0, handle spacing/signs) ---
  #   simplify_tex <- function(x) {
  #     x <- gsub("\\(1\\)", "", x)                      # remove (1)
  #     x <- gsub("\\s+", " ", x)                        # normalize spaces
  #     x <- gsub("\\+\\s*\\+", "+", x)                  # ++ -> +
  #     x <- gsub("\\+\\s*-", "-", x)                    # +- -> -
  #     x <- gsub("-\\s*\\+", "-", x)                    # -+ -> -
  #     x <- gsub("-\\s*-", "+", x)                      # -- -> +
  #     x <- gsub("\\s*\\+\\s*0\\.000\\b", "", x)        # remove + 0.000
  #     x <- gsub("\\b0\\.000\\s*\\+\\s*", "", x)        # remove leading 0.000 +
  #     # x <- gsub("\\s*\\+\\s*0\\b", "", x)              # remove + 0
  #     # x <- gsub("\\b0\\s*\\+\\s*", "", x)              # remove leading 0 +
  #     trimws(x)
  #   }
  # 
  #   # ---------- Line 1: General operator form ----------
  #   line1 <- paste0(
  #     "\\phi_p(L)\\,\\Phi_P(L^{S})\\,(1-L)^{d}(1-L^{S})^{D}Y_t",
  #     " = c + \\theta_q(L)\\,\\Theta_Q(L^{S})\\,\\varepsilon_t", drift_sym
  #   )
  # 
  #   # ---------- Line 2: Expanded summation form ----------
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
  #   # ---------- Line 3: Parameter-expanded polynomial form ----------
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
  #   # ---------- Line 4: Numeric polynomial form ----------
  #   # --- MOD: build RHS without injecting "0.000 +" when intercept is ~0 ---
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
  #   # --- MOD: remove (1), remove 0.000, fix signs/spaces ---
  #   line4 <- simplify_tex(line4)
  # 
  #   # Build aligned 4-line display
  #   aligned <- paste0(
  #     "\\begin{aligned}",
  #     line1, " \\\\[-2pt]",
  #     "\\text{------------} \\\\[-2pt]",
  #     line2, " \\\\[-2pt]",
  #     "\\text{------------} \\\\[-2pt]",
  #     line3, " \\\\[-2pt]",
  #     "\\text{------------} \\\\[-2pt]",
  #     line4,
  #     "\\end{aligned}"
  #   )
  # 
  #   list(
  #     four_lines = tex_block(aligned),
  #     line1 = tex_block(line1),
  #     line2 = tex_block(line2),
  #     line3 = tex_block(line3),
  #     line4 = tex_block(line4)
  #   )
  # })
  # 
  # # ============================================================
  # # --- MOD: Render the 4-line equation (unchanged UI output structure) ---
  # # ============================================================
  # output$manual_model_equation <- renderUI({
  #   req(manual_equations())
  #   eq <- manual_equations()
  # 
  #   tagList(
  #     HTML(eq$four_lines),
  #     tags$script(HTML("if (window.MathJax && MathJax.Hub) MathJax.Hub.Queue(['Typeset', MathJax.Hub]);"))
  #   )
  # })

  
  
  
  
  # # ============================================================
  # # --- MOD: Rebuild SARIMA equation output to match the 4-line layout in the example image ---
  # #   Line 1: General operator form (with φp, ΦP, θq, ΘQ + differencing + drift)
  # #   Line 2: Expanded summation form (Σ notation)
  # #   Line 3: Same model with parameter subscripts expanded (φ1, φ2, … ; Θ1, …) using p,q,P,Q,s
  # #   Line 4: Same as line 3 but coefficients replaced by estimated numeric values
  # #   Notes:
  # #   - Uses \\[ ... \\] + aligned block for reliable MathJax rendering in Shiny.
  # #   - Avoids 1:0 bug by using seq_len() only when order > 0.
  # # ============================================================
  # manual_equations <- reactive({
  #   req(manual_fit(), ts_train_test())
  # 
  #   fit <- manual_fit()
  #   coefs <- coef(fit)
  # 
  #   # Orders
  #   p <- as.integer(input$p); d <- as.integer(input$d); q <- as.integer(input$q)
  #   P <- as.integer(input$P); D <- as.integer(input$D); Q <- as.integer(input$Q)
  #   s <- if (!is.na(input$s) && input$s >= 1) as.integer(input$s) else frequency(ts_train_test()$ts_train)
  # 
  #   ip <- if (p > 0) seq_len(p) else integer(0)
  #   iq <- if (q > 0) seq_len(q) else integer(0)
  #   iP <- if (P > 0) seq_len(P) else integer(0)
  #   iQ <- if (Q > 0) seq_len(Q) else integer(0)
  # 
  #   # Intercept/mean
  #   intercept_name <- intersect(c("intercept", "mean"), names(coefs))
  #   intercept_val <- if (length(intercept_name) > 0) unname(coefs[intercept_name[1]]) else 0
  #   intercept_num <- sprintf("%.3f", intercept_val)
  # 
  #   # Drift flag
  #   drift_sym <- if (isTRUE(input$manual_drift)) " + \\delta t" else ""
  #   drift_num <- if (isTRUE(input$manual_drift) && "drift" %in% names(coefs)) {
  #     paste0(" + ", sprintf("%.3f", unname(coefs["drift"])), "t")
  #   } else if (isTRUE(input$manual_drift)) {
  #     " + \\delta t"
  #   } else ""
  # 
  #   # --- helpers ---
  #   tex_block <- function(x) paste0("\\[", x, "\\]")
  # 
  #   # Build param-polynomial strings like: (1 - φ1 L^1 - φ2 L^2)
  #   poly_param_ar <- function() {
  #     if (p == 0) return("1")
  #     paste0("1", paste0(" - \\phi_{", ip, "}L^{", ip, "}", collapse = ""))
  #   }
  #   poly_param_sar <- function() {
  #     if (P == 0) return("1")
  #     paste0("1", paste0(" - \\Phi_{", iP, "}L^{", s * iP, "}", collapse = ""))
  #   }
  #   poly_param_ma <- function() {
  #     if (q == 0) return("1")
  #     paste0("1", paste0(" + \\theta_{", iq, "}L^{", iq, "}", collapse = ""))
  #   }
  #   poly_param_sma <- function() {
  #     if (Q == 0) return("1")
  #     paste0("1", paste0(" + \\Theta_{", iQ, "}L^{", s * iQ, "}", collapse = ""))
  #   }
  # 
  #   # Build numeric-polynomial strings with correct signs:
  #   # AR/SAR appear as "- ar_i L^i" inside (1 - ar_i L^i) — if ar_i is negative, becomes +.
  #   # MA/SMA appear as "+ ma_i L^i" inside (1 + ma_i L^i) — if ma_i is negative, becomes -.
  #   poly_num_ar <- function() {
  #     if (p == 0) return("1")
  #     v <- if (length(ip) > 0) unname(coefs[paste0("ar", ip)]) else numeric(0)
  #     paste0("1", paste(sprintf(" %+.3fL^{%d}", -v, ip), collapse = ""))
  #   }
  #   poly_num_sar <- function() {
  #     if (P == 0) return("1")
  #     v <- if (length(iP) > 0) unname(coefs[paste0("sar", iP)]) else numeric(0)
  #     paste0("1", paste(sprintf(" %+.3fL^{%d}", -v, s * iP), collapse = ""))
  #   }
  #   poly_num_ma <- function() {
  #     if (q == 0) return("1")
  #     v <- if (length(iq) > 0) unname(coefs[paste0("ma", iq)]) else numeric(0)
  #     paste0("1", paste(sprintf(" %+.3fL^{%d}", v, iq), collapse = ""))
  #   }
  #   poly_num_sma <- function() {
  #     if (Q == 0) return("1")
  #     v <- if (length(iQ) > 0) unname(coefs[paste0("sma", iQ)]) else numeric(0)
  #     paste0("1", paste(sprintf(" %+.3fL^{%d}", v, s * iQ), collapse = ""))
  #   }
  # 
  #   # Differencing operators
  #   diff_part  <- if (d > 0) paste0("(1-L)^{", d, "}") else ""
  #   sdiff_part <- if (D > 0) paste0("(1-L^{", s, "})^{", D, "}") else ""
  # 
  #   # ---------- Line 1: General operator form ----------
  #   line1 <- paste0(
  #     "\\phi_p(L)\\,\\Phi_P(L^{S})\\,(1-L)^{d}(1-L^{S})^{D}Y_t",
  #     " = c + \\theta_q(L)\\,\\Theta_Q(L^{S})\\,\\varepsilon_t", drift_sym
  #   )
  # 
  #   # ---------- Line 2: Expanded summation form ----------
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
  #   # ---------- Line 3: Parameter-expanded polynomial form ----------
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
  #   # ---------- Line 4: Numeric polynomial form ----------
  #   line4 <- paste0(
  #     "(", poly_num_ar(), ")",
  #     "(", poly_num_sar(), ")",
  #     diff_part, sdiff_part,
  #     "Y_t = ", intercept_num, " + ",
  #     "(", poly_num_ma(), ")",
  #     "(", poly_num_sma(), ")\\varepsilon_t",
  #     drift_num
  #   )
  # 
  #   # Build aligned 4-line display with separators like your image
  #   aligned <- paste0(
  #     "\\begin{aligned}",
  #     line1, " \\\\[-2pt]",
  #     "\\text{------------} \\\\[-2pt]",
  #     line2, " \\\\[-2pt]",
  #     "\\text{------------} \\\\[-2pt]",
  #     line3, " \\\\[-2pt]",
  #     "\\text{------------} \\\\[-2pt]",
  #     line4,
  #     "\\end{aligned}"
  #   )
  # 
  #   list(
  #     # This is the one you want to show in the tab
  #     four_lines = tex_block(aligned),
  # 
  #     # Keep these available if you want separate sub-panels later
  #     line1 = tex_block(line1),
  #     line2 = tex_block(line2),
  #     line3 = tex_block(line3),
  #     line4 = tex_block(line4)
  #   )
  # })
  # 
  # # ============================================================
  # # --- MOD: Render the 4-line equation exactly like the reference image ---
  # # ============================================================
  # output$manual_model_equation <- renderUI({
  #   req(manual_equations())
  #   eq <- manual_equations()
  # 
  #   tagList(
  #     HTML(eq$four_lines),
  #     # MOD: force MathJax to typeset dynamically inserted content
  #     tags$script(HTML("if (window.MathJax && MathJax.Hub) MathJax.Hub.Queue(['Typeset', MathJax.Hub]);"))
  #   )
  # })
  
  
  
  
  
  
  # # ============================================================
  # # --- MOD: Build ALL SARIMA equations (symbolic + numerical) safely ---
  # #      Fixes:
  # #      1) Avoid 1:0 -> c(1,0) bug (creates L^0 terms) by using seq_len() only if >0
  # #      2) Handle coefficient signs correctly (no "1 - -0.57L")
  # #      3) Use \\[ ... \\] instead of $$ ... $$ (more reliable with MathJax in Shiny)
  # #      4) Force MathJax re-typeset for dynamically-rendered UI
  # # ============================================================
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
  #   # Seasonal period
  #   s <- if (!is.na(input$s) && input$s >= 1) as.integer(input$s) else frequency(ts_train_test()$ts_train)
  #   
  #   # Safe index vectors (prevents 1:0 = c(1,0))
  #   ip <- if (p > 0) seq_len(p) else integer(0)
  #   iq <- if (q > 0) seq_len(q) else integer(0)
  #   iP <- if (P > 0) seq_len(P) else integer(0)
  #   iQ <- if (Q > 0) seq_len(Q) else integer(0)
  #   
  #   # Intercept/mean handling
  #   intercept_name <- intersect(c("intercept", "mean"), names(coefs))
  #   intercept_val <- if (length(intercept_name) > 0) unname(coefs[intercept_name[1]]) else 0
  #   intercept <- sprintf("%.3f", intercept_val)
  #   
  #   # Drift (kept symbolic; only shown if requested)
  #   drift <- if (isTRUE(input$manual_drift)) " + \\delta t" else ""
  #   
  #   # Helpers
  #   tex_block <- function(x) paste0("\\[", x, "\\]")
  #   tex_poly_ar <- function(phi, lags) {
  #     # φ(L) = 1 - φ1 L - ... ; terms are " %+ .3f L^{k}" with coefficient = -phi_k
  #     if (length(lags) == 0) return("1")
  #     paste0("1", paste(sprintf(" %+.3fL^{%d}", -phi, lags), collapse = ""))
  #   }
  #   tex_poly_ma <- function(theta, lags) {
  #     # θ(L) = 1 + θ1 L + ... ; terms are " %+ .3f L^{k}" with coefficient = +theta_k
  #     if (length(lags) == 0) return("1")
  #     paste0("1", paste(sprintf(" %+.3fL^{%d}", theta, lags), collapse = ""))
  #   }
  #   tex_poly_sar <- function(Phi, lagsS) {
  #     if (length(lagsS) == 0) return("1")
  #     paste0("1", paste(sprintf(" %+.3fL^{%d}", -Phi, lagsS), collapse = ""))
  #   }
  #   tex_poly_sma <- function(Theta, lagsS) {
  #     if (length(lagsS) == 0) return("1")
  #     paste0("1", paste(sprintf(" %+.3fL^{%d}", Theta, lagsS), collapse = ""))
  #   }
  #   
  #   # ---------------- SYMBOLIC (operator form + definitions) ----------------
  #   sym_ar  <- if (p > 0) paste0(" - \\phi_{", ip, "}L^{", ip, "}") else ""
  #   sym_ma  <- if (q > 0) paste0(" + \\theta_{", iq, "}L^{", iq, "}") else ""
  #   sym_sar <- if (P > 0) paste0(" - \\Phi_{", iP, "}L^{", s * iP, "}") else ""
  #   sym_sma <- if (Q > 0) paste0(" + \\Theta_{", iQ, "}L^{", s * iQ, "}") else ""
  #   
  #   symbolic <- paste0(
  #     "\\phi_p(L)\\,\\Phi_P(L^{", s, "})(1-L)^{", d, "}(1-L^{", s, "})^{", D, "}Y_t",
  #     " = c + \\theta_q(L)\\,\\Theta_Q(L^{", s, "})\\varepsilon_t", drift,
  #     " \\\\ \\text{where} \\\\ ",
  #     "\\phi_p(L)=1", sym_ar, " \\\\ ",
  #     "\\Phi_P(L^{", s, "})=1", sym_sar, " \\\\ ",
  #     "\\theta_q(L)=1", sym_ma, " \\\\ ",
  #     "\\Theta_Q(L^{", s, "})=1", sym_sma
  #   )
  #   symbolic <- tex_block(symbolic)
  #   
  #   # ---------------- NUMERICAL (operator polynomial form) ----------------
  #   ar_names  <- if (p > 0) paste0("ar", ip) else character(0)
  #   ma_names  <- if (q > 0) paste0("ma", iq) else character(0)
  #   sar_names <- if (P > 0) paste0("sar", iP) else character(0)
  #   sma_names <- if (Q > 0) paste0("sma", iQ) else character(0)
  #   
  #   ar_vals  <- if (length(ar_names) > 0) unname(coefs[ar_names]) else numeric(0)
  #   ma_vals  <- if (length(ma_names) > 0) unname(coefs[ma_names]) else numeric(0)
  #   sar_vals <- if (length(sar_names) > 0) unname(coefs[sar_names]) else numeric(0)
  #   sma_vals <- if (length(sma_names) > 0) unname(coefs[sma_names]) else numeric(0)
  #   
  #   poly_ar  <- tex_poly_ar(ar_vals, ip)
  #   poly_sar <- tex_poly_sar(sar_vals, s * iP)
  #   poly_ma  <- tex_poly_ma(ma_vals, iq)
  #   poly_sma <- tex_poly_sma(sma_vals, s * iQ)
  #   
  #   diff_part  <- if (d > 0) paste0("(1-L)^{", d, "}") else ""
  #   sdiff_part <- if (D > 0) paste0("(1-L^{", s, "})^{", D, "}") else ""
  #   
  #   numerical <- paste0(
  #     "(", poly_ar, ")(", poly_sar, ")", diff_part, sdiff_part, "Y_t",
  #     " = ", intercept,
  #     "(", poly_ma, ")(", poly_sma, ")\\varepsilon_t", drift
  #   )
  #   numerical <- tex_block(numerical)
  #   
  #   # ---------------- NUMERICAL (one-line compact) ----------------
  #   numerical_one_line <- tex_block(
  #     paste0(
  #       "(", poly_ar, ")(", poly_sar, ")", diff_part, sdiff_part, "Y_t",
  #       " = ", intercept,
  #       "(", poly_ma, ")(", poly_sma, ")\\varepsilon_t", drift
  #     )
  #   )
  #   
  #   # ---------------- TIME-DOMAIN (direct form; partial, teaching-friendly) ----------------
  #   # Note: A full time-domain expansion with d/D generally introduces many lagged terms.
  #   # Here we provide the "AR + MA + intercept" representation for intuition.
  #   # (We label it explicitly to avoid implying it's the fully expanded differenced form.)
  #   td <- paste0(
  #     "Y_t = ", intercept,
  #     if (p > 0) paste0(paste(sprintf(" %+.3fY_{t-%d}", ar_vals, ip), collapse = "")) else "",
  #     if (P > 0) paste0(paste(sprintf(" %+.3fY_{t-%d}", sar_vals, s * iP), collapse = "")) else "",
  #     if (q > 0) paste0(paste(sprintf(" %+.3f\\varepsilon_{t-%d}", ma_vals, iq), collapse = "")) else "",
  #     if (Q > 0) paste0(paste(sprintf(" %+.3f\\varepsilon_{t-%d}", sma_vals, s * iQ), collapse = "")) else "",
  #     " + \\varepsilon_t",
  #     if (d > 0 || D > 0) paste0(" \\\\ \\text{(with differencing: }(1-L)^{", d, "}(1-L^{", s, "})^{", D, "}\\text{)}") else ""
  #   )
  #   numerical_one_line_Y_t <- tex_block(td)
  #   
  #   list(
  #     symbolic = symbolic,
  #     numerical = numerical,
  #     numerical_one_line = numerical_one_line,
  #     numerical_one_line_Y_t = numerical_one_line_Y_t
  #   )
  # })
  # 
  # # ============================================================
  # # --- MOD: Render all SARIMA equations in the UI + retrigger MathJax typeset ---
  # # ============================================================
  # output$manual_model_equation <- renderUI({
  #   req(manual_equations())
  #   eq <- manual_equations()
  #   
  #   tagList(
  #     tags$h4("General SARIMA formulation"),
  #     HTML(eq$symbolic),
  #     tags$hr(),
  #     tags$h4("Expanded operator form"),
  #     HTML(eq$numerical),
  #     tags$hr(),
  #     tags$h4("Numerical model (one-line)"),
  #     HTML(eq$numerical_one_line),
  #     tags$hr(),
  #     tags$h4("Equivalent time-domain representation (teaching form)"),
  #     HTML(eq$numerical_one_line_Y_t),
  #     
  #     # MOD: force MathJax to typeset dynamically inserted content
  #     tags$script(HTML("if (window.MathJax && MathJax.Hub) MathJax.Hub.Queue(['Typeset', MathJax.Hub]);"))
  #   )
  # })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # # ============================================================
  # # --- MOD: Build ALL SARIMA equations (symbolic + numerical) ---
  # # ============================================================
  # manual_equations <- reactive({
  #   req(manual_fit())
  #   
  #   fit <- manual_fit()
  #   coefs <- coef(fit)
  #   
  #   p <- input$p; d <- input$d; q <- input$q
  #   P <- input$P; D <- input$D; Q <- input$Q
  #   s <- ifelse(is.na(input$s), frequency(ts_train_test()$ts_train), input$s)
  #   
  #   intercept <- if ("intercept" %in% names(coefs)) {
  #     sprintf("%.3f", coefs["intercept"])
  #   } else {
  #     "0"
  #   }
  #   
  #   drift <- if (isTRUE(input$manual_drift)) " + \\delta t" else ""
  #   
  #   # ---------- SYMBOLIC COMPONENTS ----------
  #   symbolic_ar  <- paste0(" - \\phi_", 1:p, "L^{", 1:p, "}")
  #   symbolic_ma  <- paste0(" + \\theta_", 1:q, "L^{", 1:q, "}")
  #   symbolic_sar <- paste0(" - \\Phi_", 1:P, "L^{", s * (1:P), "}")
  #   symbolic_sma <- paste0(" + \\Theta_", 1:Q, "L^{", s * (1:Q), "}")
  #   
  #   symbolic_eq <- paste0(
  #     "\\phi_p(L)\\Phi_P(L^{", s, "})(1-L)^{", d, "}(1-L^{", s, "})^{", D, "}Y_t = ",
  #     "c + \\theta_q(L)\\Theta_Q(L^{", s, "})\\varepsilon_t", drift,
  #     " \\\\ \\text{where} \\\\ ",
  #     "\\phi_p(L)=1", if (p > 0) paste0(symbolic_ar, collapse = ""), " \\\\ ",
  #     "\\Phi_P(L^{", s, "})=1", if (P > 0) paste0(symbolic_sar, collapse = ""), " \\\\ ",
  #     "\\theta_q(L)=1", if (q > 0) paste0(symbolic_ma, collapse = ""), " \\\\ ",
  #     "\\Theta_Q(L^{", s, "})=1", if (Q > 0) paste0(symbolic_sma, collapse = "")
  #   )
  #   symbolic_eq <- paste0("$$", symbolic_eq, "$$")
  #   
  #   # ---------- NUMERICAL POLYNOMIAL FORM ----------
  #   numerical_ar  <- paste0(" - ", sprintf("%.3f", coefs[names(coefs) %in% paste0("ar", 1:p)]), "L^{", 1:p, "}")
  #   numerical_ma  <- paste0(" + ", sprintf("%.3f", coefs[names(coefs) %in% paste0("ma", 1:q)]), "L^{", 1:q, "}")
  #   numerical_sar <- paste0(" - ", sprintf("%.3f", coefs[names(coefs) %in% paste0("sar", 1:P)]), "L^{", s * (1:P), "}")
  #   numerical_sma <- paste0(" + ", sprintf("%.3f", coefs[names(coefs) %in% paste0("sma", 1:Q)]), "L^{", s * (1:Q), "}")
  #   
  #   numerical_eq <- paste0(
  #     "(1", paste0(numerical_ar, collapse = ""), ")",
  #     "(1", paste0(numerical_sar, collapse = ""), ")",
  #     "(1-L)^{", d, "}(1-L^{", s, "})^{", D, "}Y_t = ",
  #     intercept,
  #     "(1", paste0(numerical_ma, collapse = ""), ")",
  #     "(1", paste0(numerical_sma, collapse = ""), ")\\varepsilon_t",
  #     drift
  #   )
  #   numerical_eq <- gsub("\\(1\\)", "", numerical_eq)
  #   numerical_eq <- paste0("$$", numerical_eq, "$$")
  #   
  #   # ---------- ONE-LINE OPERATOR FORM ----------
  #   numerical_one_line <- paste0(
  #     "(1", paste0(numerical_ar, collapse = ""), ")",
  #     "(1", paste0(numerical_sar, collapse = ""), ")",
  #     "(1-L)^{", d, "}(1-L^{", s, "})^{", D, "}Y_t = ",
  #     intercept,
  #     "(1", paste0(numerical_ma, collapse = ""), ")",
  #     "(1", paste0(numerical_sma, collapse = ""), ")\\varepsilon_t",
  #     drift
  #   )
  #   numerical_one_line <- gsub("\\(1\\)", "", numerical_one_line)
  #   numerical_one_line <- paste0("$$", numerical_one_line, "$$")
  #   
  #   # ---------- TIME-DOMAIN REPRESENTATION ----------
  #   numerical_one_line_Y_t <- paste0(
  #     "Y_t = ", intercept,
  #     if (p > 0) paste0(" + ", paste0(sprintf("%.3f", coefs[names(coefs) %in% paste0("ar", 1:p)]), "Y_{t-", 1:p, "}"), collapse = ""),
  #     if (P > 0) paste0(" + ", paste0(sprintf("%.3f", coefs[names(coefs) %in% paste0("sar", 1:P)]), "Y_{t-", s * (1:P), "}"), collapse = ""),
  #     if (q > 0) paste0(" + ", paste0(sprintf("%.3f", coefs[names(coefs) %in% paste0("ma", 1:q)]), "\\varepsilon_{t-", 1:q, "}"), collapse = ""),
  #     if (Q > 0) paste0(" + ", paste0(sprintf("%.3f", coefs[names(coefs) %in% paste0("sma", 1:Q)]), "\\varepsilon_{t-", s * (1:Q), "}"), collapse = ""),
  #     " + \\varepsilon_t"
  #   )
  #   numerical_one_line_Y_t <- paste0("$$", numerical_one_line_Y_t, "$$")
  #   
  #   list(
  #     symbolic = symbolic_eq,
  #     numerical = numerical_eq,
  #     numerical_one_line = numerical_one_line,
  #     numerical_one_line_Y_t = numerical_one_line_Y_t
  #   )
  # })
  # 
  # # ============================================================
  # # --- MOD: Render all SARIMA equations in the UI (MathJax) ---
  # # ============================================================
  # output$manual_model_equation <- renderUI({
  #   req(manual_equations())
  #   eq <- manual_equations()
  #   
  #   tagList(
  #     tags$h4("General SARIMA formulation"),
  #     HTML(eq$symbolic),
  #     tags$hr(),
  #     tags$h4("Expanded numerical operator form"),
  #     HTML(eq$numerical),
  #     tags$hr(),
  #     tags$h4("Compact numerical representation"),
  #     HTML(eq$numerical_one_line),
  #     tags$hr(),
  #     tags$h4("Equivalent time-domain equation"),
  #     HTML(eq$numerical_one_line_Y_t)
  #   )
  # })
  # 
  # 
  # 
  # 
  # # ============================================================
  # # --- MOD: Render full SARIMA model equation (MathJax) ---
  # # ============================================================
  # output$manual_model_equation <- renderUI({
  #   req(manual_fit(), manual_equations())
  #   
  #   eq <- manual_equations()
  #   
  #   tagList(
  #     tags$h4("General SARIMA formulation"),
  #     HTML(eq$symbolic),
  #     
  #     tags$hr(),
  #     
  #     tags$h4("Expanded operator form"),
  #     HTML(eq$numerical),
  #     
  #     tags$hr(),
  #     
  #     tags$h4("Numerical model (one-line)"),
  #     HTML(eq$numerical_one_line),
  #     
  #     tags$hr(),
  #     
  #     tags$h4("Equivalent time-domain representation"),
  #     HTML(eq$numerical_one_line_Y_t)
  #   )
  # })
  
  
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
}
