# server.R — SARIMA Scientific Writing Lab (v2)
# FIX v2: Manual SARIMA validation forecast aligns with test horizon (h = test length).

library(shiny)
library(ggplot2)
library(forecast)
library(lubridate)
library(zoo)
library(tseries)
library(urca)
library(gridExtra)
library(colourpicker)
library(patchwork)
library(scales)
library(DT)


has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)





# ---------- Safety helpers (add once) ----------

is_finite_vec <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  length(x) > 0 && any(is.finite(x))
}

safe_range <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0) return(NULL)
  r <- suppressWarnings(range(x, na.rm = TRUE, finite = TRUE))
  if (length(r) != 2 || any(!is.finite(r))) return(NULL)
  if (diff(r) == 0) return(NULL)
  r
}

# Safely extract and standardize coefficient table from MANY model types
coef_table_safe <- function(fit) {
  sm <- tryCatch(summary(fit), error = function(e) NULL)
  if (is.null(sm)) return(data.frame())
  
  # try common places:
  co <- NULL
  if (!is.null(sm$coef)) co <- sm$coef
  if (is.null(co) && !is.null(sm$coefficients)) co <- sm$coefficients
  if (is.null(co)) return(data.frame())
  
  co <- as.data.frame(co)
  co$term <- rownames(co)
  rownames(co) <- NULL
  
  # normalize column names
  nm <- names(co)
  nm0 <- tolower(gsub("[^a-z]+", "", nm))
  
  # Attempt to find estimate/se/t/p columns flexibly
  pick <- function(keys) {
    idx <- which(nm0 %in% keys)
    if (length(idx) == 0) NA_integer_ else idx[1]
  }
  
  i_est <- pick(c("estimate", "est", "coef", "value"))
  i_se  <- pick(c("se", "stderror", "stderr", "sestd"))
  i_t   <- pick(c("tvalue", "tstat", "zvalue", "zstat", "statistic"))
  i_p   <- pick(c("prtz", "prgtz", "prgt", "prtt", "prgtz", "pvalue", "pr", "prt"))
  
  out <- data.frame(
    term = co$term,
    est  = if (!is.na(i_est)) suppressWarnings(as.numeric(co[[i_est]])) else NA_real_,
    se   = if (!is.na(i_se))  suppressWarnings(as.numeric(co[[i_se]]))  else NA_real_,
    stat = if (!is.na(i_t))   suppressWarnings(as.numeric(co[[i_t]]))   else NA_real_,
    p    = if (!is.na(i_p))   suppressWarnings(as.numeric(co[[i_p]]))   else NA_real_,
    stringsAsFactors = FALSE
  )
  
  out
}

sig_stars <- function(p) {
  if (!is.finite(p)) return("")
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else if (p < 0.1) "." else ""
}






build_xreg_split <- function(df, cols, train_n, test_n, scale_x = FALSE) {
  if (length(cols) == 0) return(list(x_train = NULL, x_test = NULL, x_all = NULL))
  
  X <- df[, cols, drop = FALSE]
  
  # keep numeric columns only (convert logical to numeric)
  for (nm in names(X)) {
    if (is.logical(X[[nm]])) X[[nm]] <- as.numeric(X[[nm]])
  }
  # force numeric
  X <- data.frame(lapply(X, function(z) suppressWarnings(as.numeric(z))), check.names = FALSE)
  
  # handle missing values in xreg (simple: linear interpolation + LOCF fallback)
  for (j in seq_along(X)) {
    v <- X[[j]]
    idx <- which(is.finite(v))
    if (length(idx) >= 2) {
      v[!is.finite(v)] <- approx(x = idx, y = v[idx], xout = which(!is.finite(v)), method = "linear", rule = 2)$y
    }
    # if still NA (all missing), set 0
    v[!is.finite(v)] <- 0
    X[[j]] <- v
  }
  
  if (isTRUE(scale_x)) {
    X <- as.data.frame(scale(X), check.names = FALSE)
  }
  
  X_mat <- as.matrix(X)
  x_train <- X_mat[seq_len(train_n), , drop = FALSE]
  x_test  <- if (test_n > 0) X_mat[(train_n + 1):(train_n + test_n), , drop = FALSE] else NULL
  
  list(x_train = x_train, x_test = x_test, x_all = X_mat)
}



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



fmt_num <- function(x, digits = 4, trim = TRUE) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || !is.finite(x)) return("NA")
  format(round(x, digits), nsmall = digits, trim = trim)
}

# fmt_num <- function(x, digits = 2) {
#   if (is.na(x)) return("NA")
#   format(round(x, digits), nsmall = digits, trim = TRUE)
# }




fmt_pct <- function(x, digits = 1) {
  if (is.na(x)) return("NA")
  paste0(fmt_num(100 * x, digits), "%")
}

# ---------------- Dates & time grid ----------------

# parse_dates <- function(x) {
#   if (inherits(x, "Date")) return(x)
#   if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
#   x_chr <- as.character(x)
# 
#   d <- suppressWarnings(as.Date(zoo::as.yearmon(x_chr)))
#   if (all(is.na(d))) {
#     d <- suppressWarnings(lubridate::parse_date_time(
#       x_chr,
#       orders = c("ymd", "dmy", "mdy", "Ymd", "Y-m-d", "d-m-Y", "m/d/Y", "Y", "ym", "my", "bY", "Y-b", "any")
#     ))
#     d <- as.Date(d)
#   }
#   d
# }

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




















#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================
#=========================================================================================================


pkg_status_li <- function(pkg) {
  ok <- has_pkg(pkg)
  tags$li(
    if (ok) "✅" else "❌",
    tags$span(
      paste0(" ", pkg),
      style = if (ok) "color: #2e7d32;" else "color: #b71c1c;"
    )
  )
}


# pkg_status_li <- function(pkg) {
#   ok <- has_pkg(pkg)
#   tags$li(
#     tags$span(if (ok) "✅" else "❌", style = "display:inline-block; width: 1.4em;"),
#     tags$span(
#       pkg,
#       style = if (ok) "color:#2e7d32;" else "color:#b71c1c;"
#     )
#   )
# }


# ---------------- Shiny server ----------------

server <- function(input, output, session) {

  # ---- Roadmap & teaching notes ----

  # output$roadmap_ui <- renderUI({
  #   tags$div(
  #     style = "background:#f7f7f7;padding:12px;border-radius:8px;",
  #     tags$h4("Roadmap (what students do, what they write)"),
  #     tags$ol(
  #       tags$li(tags$b("Describe the data"), ": sample size, missing values, descriptive statistics."),
  #       tags$li(tags$b("Explore visually"), ": trend/seasonality/outliers; report observations."),
  #       tags$li(tags$b("Decompose"), ": justify additive vs multiplicative; use STL when robust needed."),
  #       tags$li(tags$b("Check stationarity"), ": ADF/KPSS/PP; justify differencing (d and D)."),
  #       tags$li(tags$b("Fit a baseline model"), ": Auto-ARIMA to obtain a strong starting SARIMA."),
  #       tags$li(tags$b("Fit a theory-driven model"), ": Manual SARIMA using ACF/PACF + tests."),
  #       tags$li(tags$b("Diagnose & compare"), ": residual tests + forecast accuracy; choose final model."),
  #       tags$li(tags$b("Write your paper"), ": use APA paragraphs in each step; assemble Methods/Results.")
  #     )
  #   )
  # })
  

  
  
  
  # Required packages (app will fail without these)
  output$package_status <- renderPrint({
    
    required_pkgs <- c(
      "shiny",        # app framework
      "ggplot2",      # all plotting
      "forecast",     # SARIMA, ACF/PACF, decomposition, BoxCox
      "lubridate",    # date parsing
      "zoo",          # na.locf, na.approx, as.yearmon
      "tseries",      # ADF, Jarque-Bera
      "urca",         # unit root tests
      "DT",           # data tables
      "scales"        # plot scales
    )
    
    optional_pkgs <- c(
      "readxl",       # only needed for XLS/XLSX upload
      "colourpicker", # S(t) plot color pickers
      "patchwork",    # combined plots
      "gridExtra",    # plot layouts
      "shinythemes",  # UI theme
      "shinyjs",      # JS helpers (disable/enable UI)
      "nortest",      # Anderson–Darling normality test
      "FinTS"         # ARCH LM test
    )
    
    
    cat("Package status\n")
    cat("==============\n\n")
    
    cat("Required packages:\n")
    for (p in required_pkgs) {
      cat(" -", p, ":", has_pkg(p), "\n")
    }
    
    cat("\nOptional packages:\n")
    for (p in optional_pkgs) {
      cat(" -", p, ":", has_pkg(p), "\n")
    }
    
    if (any(!sapply(required_pkgs, has_pkg))) {
      cat("\n⚠️  Some required packages are missing.\n")
      cat("Install them with:\n")
      cat("install.packages(c(",
          paste0('"', required_pkgs, '"', collapse = ", "),
          "))\n")
    } else {
      cat("\n✅ All required packages are installed.\n")
    }
  })
  
  
  

  
  
  # output$package_status <- renderPrint({
  #   
  #   required_pkgs <- c(
  #     "shiny", "ggplot2", "forecast", "lubridate", "zoo",
  #     "tseries", "urca", "DT", "scales"
  #   )
  #   
  #   optional_pkgs <- c(
  #     "gridExtra", "patchwork", "colourpicker",
  #     "nortest", "FinTS", "readxl",
  #     "shinythemes", "shinyjs"
  #   )
  #   
  #   cat("Package status\n")
  #   cat("==============\n\n")
  #   
  #   cat("Required packages:\n")
  #   for (p in required_pkgs) {
  #     cat(" -", p, ":", has_pkg(p), "\n")
  #   }
  #   
  #   cat("\nOptional packages:\n")
  #   for (p in optional_pkgs) {
  #     cat(" -", p, ":", has_pkg(p), "\n")
  #   }
  #   
  #   if (any(!sapply(required_pkgs, has_pkg))) {
  #     cat("\n⚠️  Some required packages are missing.\n")
  #     cat("Install them with:\n")
  #     cat("install.packages(c(",
  #         paste0('"', required_pkgs, '"', collapse = ", "),
  #         "))\n")
  #   } else {
  #     cat("\n✅ All required packages are installed.\n")
  #   }
  # })
  
  
  
  
  
  # output$roadmap_ui <- renderUI({
  #   
  #   required_pkgs <- c(
  #     "shiny", "ggplot2", "forecast", "lubridate", "zoo",
  #     "tseries", "urca", "DT", "scales"
  #   )
  #   
  #   optional_pkgs <- c(
  #     "readxl", "colourpicker", "patchwork", "gridExtra",
  #     "shinythemes", "shinyjs", "nortest", "FinTS"
  #   )
  #   
  #   tags$div(
  #     style = "background:#f7f7f7;padding:12px;border-radius:8px;",
  #     
  #     tags$h4("\n💡 Roadmap (what students do, what they write)"),
  #     
  #     # 🔹 NEW: Package status block
  #     tags$div(
  #       style = "margin-bottom:12px;",
  #       
  #       tags$b("R environment check"),
  #       
  #       tags$ul(
  #         style = "margin-top:6px;",
  #         tags$li(tags$b("Required packages 📦")),
  #         tags$ul(lapply(required_pkgs, pkg_status_li)),
  #         # tags$li(tags$b("Optional packages ")),
  #         tags$ul(lapply(optional_pkgs, pkg_status_li))
  #       )
  #     ),
  #     
  #     tags$hr(),
  #     
  #     # 🔹 Existing roadmap content
  #     tags$ol(
  #       tags$li(tags$b("Describe the data"), ": sample size, missing values, descriptive statistics."),
  #       tags$li(tags$b("Explore visually"), ": trend/seasonality/outliers; report observations."),
  #       tags$li(tags$b("Decompose"), ": justify additive vs multiplicative; use STL when robust needed."),
  #       tags$li(tags$b("Check stationarity"), ": ADF/KPSS/PP; justify differencing (d and D)."),
  #       tags$li(tags$b("Fit a baseline model"), ": Auto-ARIMA to obtain a strong starting SARIMA."),
  #       tags$li(tags$b("Fit a theory-driven model"), ": Manual SARIMA using ACF/PACF + tests."),
  #       tags$li(tags$b("Diagnose & compare"), ": residual tests + forecast accuracy; choose final model."),
  #       tags$li(tags$b("Write your paper"), ": use APA paragraphs in each step; assemble Methods/Results.")
  #     )
  #   )
  # })
  
  
  output$roadmap_ui <- renderUI({
    
    required_pkgs <- c(
      "shiny", "ggplot2", "forecast", "lubridate", "zoo",
      "tseries", "urca", "DT", "scales"
    )
    
    optional_pkgs <- c(
      "readxl", "colourpicker", "patchwork", "gridExtra",
      "shinythemes", "shinyjs", "nortest", "FinTS"
    )
    
    # Helper for roadmap items with icons
    step_li <- function(ic, title_bold, rest_text) {
      tags$li(
        tags$span(icon(ic), style = "margin-right:8px; color:#2c3e50;"),
        tags$b(title_bold),
        HTML(paste0(": ", rest_text))
      )
    }
    
    tags$div(
      
      # =========================
      # Box 1: Roadmap steps
      # =========================
      tags$div(
        style = "background:#f7f7f7;padding:12px;border-radius:8px;margin-bottom:10px;",
        tags$h4("Roadmap (what students do, what they write)"),
        tags$ol(
          step_li("database",      "Describe the data",   "sample size, missing values, descriptive statistics."),
          step_li("chart-line",    "Explore visually",    "trend/seasonality/outliers; report observations."),
          step_li("layer-group",   "Decompose",           "justify additive vs multiplicative; use STL when robust needed."),
          step_li("check-circle",  "Check stationarity",  "ADF/KPSS/PP; justify differencing (d and D)."),
          step_li("robot",         "Fit a baseline model","Auto-ARIMA to obtain a strong starting SARIMA."),
          step_li("sliders-h",     "Fit a theory-driven model","Manual SARIMA using ACF/PACF + tests."),
          step_li("stethoscope",   "Diagnose & compare",  "residual tests + forecast accuracy; choose final model."),
          step_li("file-alt",      "Write your paper",    "use APA paragraphs in each step; assemble Methods/Results.")
        ),
        tags$br(),
      ),
      
      # =========================
      # Box 2: Package status (below roadmap)
      # =========================
      tags$div(
        style = "background:#eef5ff;padding:12px;border-radius:8px;",
        tags$h4("R environment check"),
        tags$p(
          style = "margin-top:-6px; font-size: 13px; color:#34495e;",
          "✅ installed  •  ❌ missing"
        ),
        
        tags$br(),
        
        fluidRow(
          
          # =========================
          # LEFT: Required packages
          # =========================
          column(
            width = 2,
            tags$b("Required packages"),
            tags$ul(lapply(required_pkgs, pkg_status_li))
          ),
          
          # =========================
          # RIGHT: Optional packages
          # =========================
          column(
            width = 2,
            tags$b("Optional packages"),
            tags$ul(lapply(optional_pkgs, pkg_status_li))
          )
        ),
        
        # Optional warning if required packages missing
        if (any(!vapply(required_pkgs, has_pkg, logical(1)))) {
          tags$div(
            style = "margin-top:10px; color:#b71c1c; font-size:13px;",
            tags$b("Some required packages are missing."),
            tags$div("Install with:"),
            tags$pre(
              style = "background:white; padding:8px; border-radius:6px;",
              paste0(
                "install.packages(c(",
                paste0('"', required_pkgs, '"', collapse = ", "),
                "))"
              )
            )
          )
        } else NULL
      )
      
      
      
      # tags$div(
      #   style = "background:#eef5ff;padding:12px;border-radius:8px;",
      #   tags$h4("R environment check"),
      #   tags$p(
      #     style = "margin-top:-6px; font-size: 13px; color:#34495e;",
      #     "✅ installed  •  ❌ missing"
      #   ),
      #   
      #   tags$b("Required packages"),
      #   tags$ul(lapply(required_pkgs, pkg_status_li)),
      #   
      #   tags$b("Optional packages"),
      #   tags$ul(lapply(optional_pkgs, pkg_status_li)),
      #   
      #   # Optional install hint (only if required missing)
      #   if (any(!vapply(required_pkgs, has_pkg, logical(1)))) {
      #     tags$div(
      #       style = "margin-top:10px; color:#b71c1c; font-size: 13px;",
      #       tags$b("Some required packages are missing."),
      #       tags$div("Install with:"),
      #       tags$pre(
      #         style = "background:white; padding:8px; border-radius:6px;",
      #         paste0(
      #           "install.packages(c(",
      #           paste0('"', required_pkgs, '"', collapse = ", "),
      #           "))"
      #         )
      #       )
      #     )
      #   } else NULL
      # )
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
  
  # output$full_data_table <- DT::renderDataTable({
  #   req(raw_data())
  #   
  #   df <- raw_data()
  #   
  #   # optional: make Date/POSIX columns print nicely
  #   for (nm in names(df)) {
  #     if (inherits(df[[nm]], "Date")) df[[nm]] <- format(df[[nm]], "%Y-%m-%d")
  #     if (inherits(df[[nm]], "POSIXt")) df[[nm]] <- format(df[[nm]], "%Y-%m-%d %H:%M:%S")
  #   }
  #   
  #   DT::datatable(
  #     df,
  #     rownames = FALSE,
  #     filter = "top",
  #     extensions = c("Scroller"),
  #     options = list(
  #       deferRender = TRUE,
  #       scrollX = TRUE,
  #       scrollY = 520,
  #       scroller = TRUE,
  #       pageLength = 25,
  #       lengthMenu = list(c(10, 25, 50, 100, -1), c("10", "25", "50", "100", "All"))
  #     )
  #   )
  # })
  
  output$full_data_table <- DT::renderDataTable({
    req(raw_data())
    
    df <- raw_data()
    
    # ✅ Force Date & POSIX columns to Day-Month-Year
    for (nm in names(df)) {
      if (inherits(df[[nm]], "Date")) {
        df[[nm]] <- format(df[[nm]], "%d-%m-%Y")
      }
      if (inherits(df[[nm]], "POSIXt")) {
        df[[nm]] <- format(df[[nm]], "%d-%m-%Y %H:%M:%S")
      }
    }
    
    DT::datatable(
      df,
      rownames = FALSE,
      filter = "top",
      extensions = "Scroller",
      options = list(
        deferRender = TRUE,
        scrollX = TRUE,
        scrollY = 520,
        scroller = TRUE,
        pageLength = 25,
        lengthMenu = list(
          c(10, 25, 50, 100, -1),
          c("10", "25", "50", "100", "All")
        )
      )
    )
  })
  
  
  
  output$data_preview <- renderTable({
    req(prepared())
    # df <- head(prepared()$df, 12)
    df <- prepared()$df
    
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
  
  
  # output$basic_stats <- renderTable({ req(prepared()); basic_stats_df(prepared()$df$y_filled) }, rownames = FALSE)

  # output$hist_plot <- renderPlot({
  #   req(prepared())
  #   df <- prepared()$df
  #   ggplot(df, aes(x = y_filled)) +
  #     geom_histogram(bins = 30) +
  #     theme_minimal() +
  #     labs(title = "Distribution (filled values)", x = "Value", y = "Count")
  # })

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
  
  output$plot_series <- renderPlot(
    {
      req(prepared(), ts_train_test())
      p <- prepared()
      s <- ts_train_test()
      df <- s$dfm
      df$set <- ifelse(seq_len(nrow(df)) <= s$train_n, "Train", "Test/Future")
      
      ggplot(df, aes(x = x, y = y_trans, color = set)) +
        geom_line(linewidth = 0.9) +
        theme_minimal() +
        labs(
          title = "Time series (transformed)",
          x = p$x_label,
          y = "Value",
          color = NULL
        ) +
        theme(legend.position = "bottom")
    },
    width = 1000,
    height = 650
  )

  # output$plot_series <- renderPlot({
  #   req(prepared(), ts_train_test())
  #   p <- prepared()
  #   s <- ts_train_test()
  #   df <- s$dfm
  #   df$set <- ifelse(seq_len(nrow(df)) <= s$train_n, "Train", "Test/Future")
  #   ggplot(df, aes(x = x, y = y_trans, color = set)) +
  #     geom_line(linewidth = 0.9) +
  #     theme_minimal() +
  #     labs(title = "Time series (transformed)", x = p$x_label, y = "Value", color = NULL) +
  #     theme(legend.position = "bottom")
  # })

  # output$season_plot <- renderPlot({
  #   req(ts_train_test())
  #   s <- ts_train_test()
  #   x <- ts(c(as.numeric(s$ts_train), as.numeric(s$ts_test)), start = 1, frequency = frequency(s$ts_train))
  #   validate(need(frequency(x) >= 2, "Seasonal plots need frequency >= 2."))
  #   forecast::seasonplot(x, s = frequency(x))
  # })

  # output$subseries_plot <- renderPlot({
  #   req(ts_train_test())
  #   s <- ts_train_test()
  #   x <- ts(c(as.numeric(s$ts_train), as.numeric(s$ts_test)), start = 1, frequency = frequency(s$ts_train))
  #   validate(need(frequency(x) >= 2, "Subseries plot needs frequency >= 2."))
  #   forecast::ggsubseriesplot(x) + theme_minimal()
  # })

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

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  
  
  
  
  
  # ============================================================
  # S(t) plots tab (SERVER) — paste INSIDE server()
  # ============================================================
  
  # ---- helpers ----
  to_num <- function(x, d = NA_real_) {
    y <- suppressWarnings(as.numeric(x))
    ifelse(is.finite(y), y, d)
  }
  
  # fmt_num <- function(x, d = 4) {
  #   if (!is.finite(x)) return("NA")
  #   format(round(x, d), nsmall = d)
  # }
  
  # ============================================================
  # 0) Scope control: prevent "Training only" when train_prop == 1.00
  # ============================================================
  output$stp_scope_ui <- renderUI({
    tp <- to_num(input$train_prop, 1)
    has_test <- isTRUE(tp < 1)
    
    if (!has_test) {
      radioButtons(
        "stp_scope",
        label = NULL,
        choices = c("Full series" = "full"),
        selected = "full"
      )
    } else {
      radioButtons(
        "stp_scope",
        label = NULL,
        choices = c("Full series" = "full", "Training only" = "train"),
        selected = input$stp_scope %||% "full"
      )
    }
  })
  
  output$stp_scope_warning <- renderUI({
    tp <- to_num(input$train_prop, 1)
    if (isTRUE(tp >= 1)) {
      tags$div(
        style = "color:#b22222; font-size:12px; margin-top:-6px;",
        "No test set → training = full series"
      )
    } else NULL
  })
  
  observeEvent(input$train_prop, {
    tp <- to_num(input$train_prop, 1)
    if (isTRUE(tp >= 1) && identical(input$stp_scope, "train")) {
      updateRadioButtons(session, "stp_scope", selected = "full")
    }
  }, ignoreInit = TRUE)
  
  # ============================================================
  # 1) Global transform info (read-only)
  # ============================================================
  output$stp_transform_info <- renderPrint({
    tr <- input$transform %||% "none"
    if (tr == "boxcox") {
      lam <- input$lambda
      cat("Transformation:", "Box-Cox", "\n")
      cat("λ:", if (is.null(lam) || (length(lam) == 1 && is.na(lam))) "auto (estimated)" else as.character(lam), "\n")
    } else if (tr == "log") {
      cat("Transformation:", "Log (ln y)", "\n")
    } else {
      cat("Transformation:", "None", "\n")
    }
  })
  
  output$stp_transform_note <- renderUI({
    tr <- input$transform %||% "none"
    lam <- input$lambda
    
    tr_lbl <- switch(tr,
                     "none" = "None",
                     "log" = "Log (ln y)",
                     "boxcox" = "Box-Cox",
                     tr)
    
    if (identical(tr, "boxcox")) {
      tags$div(
        style = "font-size:12px;",
        tags$b("Transformation (global): "), tr_lbl, " — ",
        tags$b("λ: "), if (is.null(lam) || is.na(lam)) "auto" else as.character(lam)
      )
    } else {
      tags$div(
        style = "font-size:12px;",
        tags$b("Transformation (global): "), tr_lbl
      )
    }
  })
  
  # ============================================================
  # 2) Theme/palette helpers  (IMPORTANT: rotation is controlled here)
  # ============================================================
  stp_theme_picker <- function(key) {
    switch(
      key,
      "Minimal" = ggplot2::theme_minimal(),
      "Classic" = ggplot2::theme_classic(),
      "Light"   = ggplot2::theme_light(),
      "Dark"    = ggplot2::theme_dark(),
      "BW"      = ggplot2::theme_bw(),
      "Void"    = ggplot2::theme_void(),
      ggplot2::theme_gray()
    )
  }
  
  stp_apply_theme <- function(g) {
    ang <- to_num(input$stp_x_angle, 30)
    rot <- isTRUE(input$stp_x_rotate %||% TRUE)
    
    g +
      stp_theme_picker(input$stp_theme %||% "Minimal") +
      ggplot2::theme(
        text = ggplot2::element_text(size = to_num(input$stp_base_size, 12)),
        plot.title = ggplot2::element_text(hjust = 0.5),
        
        # ✅ rotation is applied HERE so it always wins over the theme preset
        axis.text.x = ggplot2::element_text(
          angle = if (rot) ang else 0,
          hjust = if (!rot || ang == 0) 0.5 else 1,
          vjust = if (!rot || ang == 0) 0.5 else 1
        )
      )
  }
  
  stp_apply_palette <- function(g) {
    pal <- input$stp_palette %||% "Paired (brewer)"
    rev <- isTRUE(input$stp_palette_rev)
    
    if (pal == "Viridis") {
      g +
        ggplot2::scale_color_viridis_d(direction = ifelse(rev, -1, 1)) +
        ggplot2::scale_fill_viridis_d(direction = ifelse(rev, -1, 1))
    } else {
      brewer_name <- sub(" \\(brewer\\)$", "", pal)
      g +
        ggplot2::scale_color_brewer(palette = brewer_name, direction = ifelse(rev, -1, 1)) +
        ggplot2::scale_fill_brewer(palette = brewer_name, direction = ifelse(rev, -1, 1))
    }
  }
  
  # ============================================================
  # 3) Conditional style UI (COLOUR PICKERS)
  # ============================================================
  output$stp_style_ui <- renderUI({
    req(input$stp_plot_type)
    pt <- input$stp_plot_type
    
    if (!requireNamespace("colourpicker", quietly = TRUE)) {
      return(tags$div(
        style = "color:#b22222; font-size:12px;",
        "Package 'colourpicker' is required for color pickers. Install it: install.packages('colourpicker')"
      ))
    }
    
    is_liney <- pt %in% c("Line", "Line + Points", "Smoothed (LOESS)", "Moving average",
                          "Cumulative sum", "Seasonal plot", "Seasonal subseries",
                          "Polar seasonal", "Periodogram",
                          "Classical decomposition (additive)", "Classical decomposition (multiplicative)",
                          "STL decomposition")
    
    is_pointy <- pt %in% c("Points", "Line + Points", "Smoothed (LOESS)",
                           "Lag-1 scatter", "Lag plot (1..m)", "QQ plot")
    
    is_filly  <- pt %in% c("Histogram", "Density", "Seasonal boxplot")
    
    tagList(
      if (is_liney) tagList(
        colourpicker::colourInput("stp_line_color", "Line color", value = "#2C7FB8", allowTransparent = FALSE),
        sliderInput("stp_line_width", "Line width", min = 0.1, max = 5, value = 0.9, step = 0.1)
      ),
      if (is_pointy) tagList(
        colourpicker::colourInput("stp_point_color", "Point color", value = "#2C7FB8", allowTransparent = FALSE),
        sliderInput("stp_point_size", "Point size", min = 0.5, max = 8, value = 2, step = 0.5)
      ),
      if (is_filly) tagList(
        colourpicker::colourInput("stp_fill_color", "Fill color", value = "#2C7FB8", allowTransparent = FALSE)
      )
    )
  })
  
  # ============================================================
  # 4) Data reactive: choose scope + apply GLOBAL transform
  # ============================================================
  stp_data <- reactive({
    req(prepared(), ts_train_test())
    p <- prepared()
    s <- ts_train_test()
    
    df_all <- p$df
    req(df_all)
    
    validate(need("y_filled" %in% names(df_all), "prepared()$df must contain column y_filled."))
    
    if (!("x" %in% names(df_all))) df_all$x <- seq_len(nrow(df_all))
    
    df_all <- df_all[is.finite(df_all$y_filled), , drop = FALSE]
    validate(need(nrow(df_all) >= 5, "Not enough data to plot."))
    
    train_n <- s$train_n %||% floor(nrow(df_all) * to_num(input$train_prop, 1))
    train_n <- max(2, min(as.integer(train_n), nrow(df_all)))
    
    has_test <- isTRUE(to_num(input$train_prop, 1) < 1)
    scope_in <- input$stp_scope %||% "full"
    scope <- if (!has_test) "full" else scope_in
    
    df_use <- if (identical(scope, "train")) df_all[seq_len(train_n), , drop = FALSE] else df_all
    
    tr <- input$transform %||% "none"
    y <- df_use$y_filled
    y_plot <- y
    lambda_used <- NA_real_
    
    if (tr == "log") {
      validate(need(all(y > 0, na.rm = TRUE), "Log transform requires strictly positive values."))
      y_plot <- log(y)
    } else if (tr == "boxcox") {
      validate(need(all(y > 0, na.rm = TRUE), "Box-Cox transform requires strictly positive values."))
      lam <- input$lambda
      if (is.null(lam) || (length(lam) == 1 && is.na(lam))) {
        lam <- forecast::BoxCox.lambda(y, method = "guerrero")
      } else {
        lam <- as.numeric(lam)
      }
      lambda_used <- lam
      y_plot <- forecast::BoxCox(y, lam)
    }
    
    df_use$y_plot <- as.numeric(y_plot)
    
    freq_use <- p$freq %||% 1
    ts_use <- stats::ts(df_use$y_plot, start = 1, frequency = freq_use)
    
    list(
      df = df_use,
      ts = ts_use,
      freq = freq_use,
      scope = scope,
      train_n = train_n,
      n_total = nrow(df_all),
      transform = tr,
      lambda_used = lambda_used,
      x_is_date = inherits(df_use$x, "Date"),
      x_is_dt   = inherits(df_use$x, "POSIXct") || inherits(df_use$x, "POSIXt")
    )
  })
  
  # ============================================================
  # 5) Labels helper (teaching subtitle)
  # ============================================================
  stp_labels <- function(d, default_title = "S(t)", default_y = "Value") {
    title_in <- input$stp_title %||% ""
    sub_in   <- input$stp_subtitle %||% ""
    x_in     <- input$stp_xlab %||% ""
    y_in     <- input$stp_ylab %||% ""
    
    scope_txt <- if (identical(d$scope, "train")) "training" else "full"
    n_txt <- paste0("n=", nrow(d$df))
    s_txt <- paste0("s=", d$freq)
    
    tr_txt <- switch(
      d$transform,
      "log" = "transform=log",
      "boxcox" = paste0("transform=BoxCox(λ=", ifelse(is.finite(d$lambda_used), fmt_num(d$lambda_used, 3), "auto"), ")"),
      "none" = "transform=none",
      "transform=none"
    )
    
    teaching_sub <- paste(scope_txt, n_txt, s_txt, tr_txt, sep = " • ")
    
    list(
      title = if (nzchar(title_in)) title_in else default_title,
      subtitle = if (nzchar(sub_in)) sub_in else teaching_sub,
      x = if (nzchar(x_in)) x_in else "Time",
      y = if (nzchar(y_in)) y_in else default_y
    )
  }
  
  # ============================================================
  # 6) Smart date axis (NO ROTATION here)
  # ============================================================
  # stp_apply_x_scale <- function(g, d) {
  #   if (!isTRUE(d$x_is_date) && !isTRUE(d$x_is_dt)) return(g)
  #   
  #   n_ticks <- as.integer(to_num(input$stp_x_ticks, 8))
  #   
  #   fmt_choice <- input$stp_date_format %||% "auto"
  #   fmt_custom <- input$stp_date_format_custom %||% "%Y-%m"
  #   fmt <- if (identical(fmt_choice, "custom")) fmt_custom else fmt_choice
  #   
  #   if (identical(fmt, "auto")) {
  #     n <- nrow(d$df)
  #     if (n <= 24) fmt <- "%b %Y"
  #     else if (n <= 120) fmt <- "%Y-%m"
  #     else fmt <- "%Y"
  #   }
  #   
  #   if (isTRUE(d$x_is_date)) {
  #     g + ggplot2::scale_x_date(
  #       labels = scales::date_format(fmt),
  #       breaks = scales::pretty_breaks(n = n_ticks)
  #     )
  #   } else {
  #     g + ggplot2::scale_x_datetime(
  #       labels = scales::date_format(fmt),
  #       breaks = scales::pretty_breaks(n = n_ticks)
  #     )
  #   }
  # }
  
  # stp_apply_x_scale <- function(g, d) {
  #   if (!isTRUE(d$x_is_date) && !isTRUE(d$x_is_dt)) return(g)
  #   
  #   n_ticks <- as.integer(to_num(input$stp_x_ticks, 8))
  #   
  #   fmt_choice <- input$stp_date_format %||% "auto"
  #   fmt_custom <- input$stp_date_format_custom %||% "%Y-%m"
  #   lang <- input$stp_date_lang %||% "en"
  #   
  #   # ---- format selection ----
  #   fmt <- if (identical(fmt_choice, "custom")) fmt_custom else fmt_choice
  #   
  #   if (identical(fmt, "auto")) {
  #     n <- nrow(d$df)
  #     if (n <= 24) fmt <- "%b %Y"
  #     else if (n <= 120) fmt <- "%Y-%m"
  #     else fmt <- "%Y"
  #   }
  #   
  #   # ---- locale mapping ----
  #   locale_map <- c(
  #     "en" = "en",
  #     "fr" = "fr",
  #     "ar" = "ar"
  #   )
  #   locale_use <- locale_map[[lang]] %||% "en"
  #   
  #   label_fun <- scales::label_date(
  #     format = fmt,
  #     locale = locale_use
  #   )
  #   
  #   if (isTRUE(d$x_is_date)) {
  #     g + ggplot2::scale_x_date(
  #       labels = label_fun,
  #       breaks = scales::pretty_breaks(n = n_ticks)
  #     )
  #   } else {
  #     g + ggplot2::scale_x_datetime(
  #       labels = label_fun,
  #       breaks = scales::pretty_breaks(n = n_ticks)
  #     )
  #   }
  # }
  
  
  stp_apply_x_scale <- function(g, d) {
    if (!isTRUE(d$x_is_date) && !isTRUE(d$x_is_dt)) return(g)
    
    n_ticks <- as.integer(to_num(input$stp_x_ticks, 8))
    
    fmt_choice <- input$stp_date_format %||% "auto"
    fmt_custom <- input$stp_date_format_custom %||% "%Y-%m"
    lang <- input$stp_date_lang %||% "en"
    
    # ---- format selection ----
    fmt <- if (identical(fmt_choice, "custom")) fmt_custom else fmt_choice
    
    if (identical(fmt, "auto")) {
      n <- nrow(d$df)
      if (n <= 24) fmt <- "%b %Y"
      else if (n <= 120) fmt <- "%Y-%m"
      else fmt <- "%Y"
    }
    
    # ---- locale mapping ----
    locale_map <- c(
      "en" = "en",
      "fr" = "fr",
      "ar" = "ar"
    )
    locale_use <- locale_map[[lang]] %||% "en"
    
    # ---- label function (locale-aware) ----
    base_label_fun <- scales::label_date(
      format = fmt,
      locale = locale_use
    )
    
    # ---- force Western digits for Arabic labels only ----
    label_fun <- function(x) {
      lab <- base_label_fun(x)
      
      if (identical(lang, "ar")) {
        # Arabic-Indic: ٠١٢٣٤٥٦٧٨٩
        arabic_indic <- c("٠","١","٢","٣","٤","٥","٦","٧","٨","٩")
        # Eastern Arabic-Indic (Persian): ۰۱۲۳۴۵۶۷۸۹
        eastern_arabic_indic <- c("۰","۱","۲","۳","۴","۵","۶","۷","۸","۹")
        western <- as.character(0:9)
        
        lab <- as.character(lab)
        for (i in 0:9) {
          lab <- gsub(arabic_indic[i + 1], western[i + 1], lab, fixed = TRUE)
          lab <- gsub(eastern_arabic_indic[i + 1], western[i + 1], lab, fixed = TRUE)
        }
      }
      
      lab
    }
    
    if (isTRUE(d$x_is_date)) {
      g + ggplot2::scale_x_date(
        labels = label_fun,
        breaks = scales::pretty_breaks(n = n_ticks)
      )
    } else {
      g + ggplot2::scale_x_datetime(
        labels = label_fun,
        breaks = scales::pretty_breaks(n = n_ticks)
      )
    }
  }
  
  
  
  # ============================================================
  # 7) UI dispatcher for multi-panel plot types
  # ============================================================
  # output$stp_plot_ui <- renderUI({
  #   req(input$stp_plot_type)
  #   pt <- input$stp_plot_type
  #   h <- to_num(input$stp_plot_height_px, 520)
  #   
  #   if (pt == "ACF+PACF") {
  #     fluidRow(
  #       column(6, plotOutput("stp_acf",  width = "100%", height = h)),
  #       column(6, plotOutput("stp_pacf", width = "100%", height = h))
  #     )
  #   } else if (pt == "Time + ACF+PACF") {
  #     tagList(
  #       plotOutput("stp_main", width = "100%", height = round(h * 0.9)),
  #       fluidRow(
  #         column(6, plotOutput("stp_acf",  width = "100%", height = round(h * 0.8))),
  #         column(6, plotOutput("stp_pacf", width = "100%", height = round(h * 0.8)))
  #       )
  #     )
  #   } else if (pt == "Lag plot (1..m)") {
  #     plotOutput("stp_lag_grid", width = "100%", height = round(h * 1.2))
  #   } else if (pt == "ACF") {
  #     plotOutput("stp_acf", width = "100%", height = h)
  #   } else if (pt == "PACF") {
  #     plotOutput("stp_pacf", width = "100%", height = h)
  #   } else {
  #     plotOutput("stp_main", width = "100%", height = h)
  #   }
  # })
  
  # output$stp_plot_ui <- renderUI({
  #   req(input$stp_plot_type)
  #   pt <- input$stp_plot_type
  #   h  <- to_num(input$stp_plot_height_px, 520)
  #   
  #   if (pt == "ACF+PACF") {
  #     
  #     fluidRow(
  #       column(6, plotOutput("stp_acf",  width = "100%", height = h)),
  #       column(6, plotOutput("stp_pacf", width = "100%", height = h))
  #     )
  #     
  #   } else if (pt == "Time + ACF+PACF") {
  #     
  #     fluidRow(
  #       column(
  #         12,
  #         
  #         # TOP: same container width as bottom
  #         plotOutput("stp_main", width = "100%", height = round(h * 0.9)),
  #         
  #         # BOTTOM: ACF + PACF inside SAME column(12)
  #         fluidRow(
  #           column(6, plotOutput("stp_acf",  width = "100%", height = round(h * 0.8))),
  #           column(6, plotOutput("stp_pacf", width = "100%", height = round(h * 0.8)))
  #         )
  #       )
  #     )
  #     
  #   } else if (pt == "Lag plot (1..m)") {
  #     
  #     plotOutput("stp_lag_grid", width = "100%", height = round(h * 1.2))
  #     
  #   } else if (pt == "ACF") {
  #     
  #     plotOutput("stp_acf", width = "100%", height = h)
  #     
  #   } else if (pt == "PACF") {
  #     
  #     plotOutput("stp_pacf", width = "100%", height = h)
  #     
  #   } else {
  #     
  #     plotOutput("stp_main", width = "100%", height = h)
  #     
  #   }
  # })
  
  
  
  output$stp_plot_ui <- renderUI({
    req(input$stp_plot_type)
    pt <- input$stp_plot_type
    h  <- to_num(input$stp_plot_height_px, 520)
    
    if (pt == "ACF+PACF") {
      
      fluidRow(
        column(6, plotOutput("stp_acf",  width = "100%", height = h)),
        column(6, plotOutput("stp_pacf", width = "100%", height = h))
      )
      
    } else if (pt == "Time + ACF+PACF") {
      
      gap_px <- round(h * 0.12)  # height of empty middle row
      
      fluidRow(
        column(
          12,
          
          # ── Row 1: Time plot ──
          plotOutput("stp_main", width = "100%", height = round(h * 0.9)),
          
          # ── Row 2: empty spacer row ──
          fluidRow(
            column(
              12,
              div(style = paste0("height:", gap_px, "px;"))
            )
          ),
          
          # ── Row 3: ACF (left-indented) + PACF (right-indented) ──
          fluidRow(
            column(
              6,
              div(
                style = "padding-left: 24px;",
                plotOutput("stp_acf", width = "100%", height = round(h * 0.6))
              )
            ),
            column(
              6,
              div(
                style = "padding-right: 24px;",
                plotOutput("stp_pacf", width = "100%", height = round(h * 0.6))
              )
            )
          )
        )
      )
      
    } else if (pt == "Lag plot (1..m)") {
      
      plotOutput("stp_lag_grid", width = "100%", height = round(h * 1.2))
      
    } else if (pt == "ACF") {
      
      plotOutput("stp_acf", width = "100%", height = h)
      
    } else if (pt == "PACF") {
      
      plotOutput("stp_pacf", width = "100%", height = h)
      
    } else {
      
      plotOutput("stp_main", width = "100%", height = h)
      
    }
  })
  
  
  
  
  # ============================================================
  # 8) Main plot  (FIXED ORDER: theme first, x-scale last)
  # ============================================================
  output$stp_main <- renderPlot({
    d <- stp_data()
    df <- d$df
    pt <- input$stp_plot_type %||% "Line"
    
    line_col <- input$stp_line_color %||% "#2C7FB8"
    lw       <- to_num(input$stp_line_width, 1)
    pt_col   <- input$stp_point_color %||% line_col
    ps       <- to_num(input$stp_point_size, 2)
    fill_col <- input$stp_fill_color %||% line_col
    a        <- to_num(input$stp_alpha, 1)
    
    labs0 <- stp_labels(d, default_title = "S(t)", default_y = "Value")
    
    base <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y_plot)) +
      ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = labs0$x, y = labs0$y)
    
    # ---- plot switch ----
    if (pt == "Line") {
      g <- base + ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a)
      g <- stp_apply_theme(g)
      g <- stp_apply_x_scale(g, d)
      return(g)
    }
    
    if (pt == "Points") {
      g <- base + ggplot2::geom_point(color = pt_col, size = ps, alpha = a)
      g <- stp_apply_theme(g)
      g <- stp_apply_x_scale(g, d)
      return(g)
    }
    
    if (pt == "Line + Points") {
      g <- base +
        ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a) +
        ggplot2::geom_point(color = pt_col, size = ps, alpha = a)
      g <- stp_apply_theme(g)
      g <- stp_apply_x_scale(g, d)
      return(g)
    }
    
    if (pt == "Smoothed (LOESS)") {
      span <- to_num(input$stp_loess_span, 0.4)
      g <- base +
        ggplot2::geom_point(color = pt_col, size = ps, alpha = a) +
        ggplot2::geom_smooth(method = "loess", span = span, se = TRUE,
                             color = line_col, linewidth = lw, alpha = 0.2)
      g <- stp_apply_theme(g)
      g <- stp_apply_x_scale(g, d)
      return(g)
    }
    
    if (pt == "Moving average") {
      k <- as.integer(to_num(input$stp_ma_k, 5))
      show_raw <- isTRUE(input$stp_ma_show_raw)
      
      df2 <- df
      df2$ma <- zoo::rollmean(df2$y_plot, k = k, fill = NA, align = "center")
      
      g <- ggplot2::ggplot(df2, ggplot2::aes(x = x)) +
        ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = labs0$x, y = labs0$y)
      
      if (show_raw) {
        g <- g + ggplot2::geom_line(ggplot2::aes(y = y_plot), color = "gray60", linewidth = 0.6, alpha = 0.7)
      }
      g <- g + ggplot2::geom_line(ggplot2::aes(y = ma), color = line_col, linewidth = lw, alpha = a)
      
      g <- stp_apply_theme(g)
      g <- stp_apply_x_scale(g, d)
      return(g)
    }
    
    if (pt == "Cumulative sum") {
      df2 <- df
      df2$cs <- cumsum(df2$y_plot)
      
      g <- ggplot2::ggplot(df2, ggplot2::aes(x = x, y = cs)) +
        ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a) +
        ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = labs0$x, y = "Cumsum")
      
      g <- stp_apply_theme(g)
      g <- stp_apply_x_scale(g, d)
      return(g)
    }
    
    
    # ---- Seasonal plots ----
    if (pt %in% c("Seasonal plot", "Seasonal subseries", "Polar seasonal", "Seasonal boxplot")) {
      validate(need(is.finite(d$freq) && d$freq >= 2, "Seasonal plots require frequency >= 2."))
      
      x_ts <- d$ts  # ts object built from y_plot + frequency
      
      if (pt == "Seasonal plot") {
        g <- forecast::ggseasonplot(x_ts, polar = FALSE) +
          ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle)
        g <- stp_apply_palette(g)
        g <- stp_apply_theme(g)
        return(g)
      }
      
      if (pt == "Polar seasonal") {
        g <- forecast::ggseasonplot(x_ts, polar = TRUE) +
          ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle)
        g <- stp_apply_palette(g)
        g <- stp_apply_theme(g)
        return(g)
      }
      
      if (pt == "Seasonal subseries") {
        g <- forecast::ggsubseriesplot(x_ts) +
          ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle)
        g <- stp_apply_theme(g)
        return(g)
      }
      
      if (pt == "Seasonal boxplot") {
        df2 <- data.frame(season = factor(stats::cycle(x_ts)), y = as.numeric(x_ts))
        g <- ggplot2::ggplot(df2, ggplot2::aes(x = season, y = y, fill = season)) +
          ggplot2::geom_boxplot(alpha = 0.6) +
          ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = "Season", y = labs0$y)
        g <- stp_apply_palette(g)
        g <- stp_apply_theme(g)
        return(g)
      }
    }
    
    # ---- Decomposition plots ----
    if (pt %in% c("Classical decomposition (additive)", "Classical decomposition (multiplicative)")) {
      validate(need(is.finite(d$freq) && d$freq >= 2, "Decomposition requires frequency >= 2."))
      
      x_ts <- d$ts
      type <- if (pt == "Classical decomposition (multiplicative)") "multiplicative" else "additive"
      
      # multiplicative requires strictly positive values
      if (type == "multiplicative") {
        validate(need(all(as.numeric(x_ts) > 0), "Multiplicative decomposition requires strictly positive values."))
      }
      
      dc <- stats::decompose(x_ts, type = type)
      g <- forecast::autoplot(dc) +
        ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle)
      g <- stp_apply_theme(g)
      return(g)
    }
    
    if (pt == "STL decomposition") {
      validate(need(is.finite(d$freq) && d$freq >= 2, "STL requires frequency >= 2."))
      
      x_ts <- d$ts
      fit <- stats::stl(x_ts, s.window = "periodic", robust = TRUE)
      
      g <- forecast::autoplot(fit) +
        ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle)
      g <- stp_apply_theme(g)
      return(g)
    }
    
    
    
    
    if (pt == "Histogram") {
      bins <- as.integer(to_num(input$stp_hist_bins, 30))
      g <- ggplot2::ggplot(df, ggplot2::aes(x = y_plot)) +
        ggplot2::geom_histogram(bins = bins, fill = fill_col, alpha = a) +
        ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = labs0$y, y = "Count")
      return(stp_apply_theme(g))
    }
    
    if (pt == "Density") {
      bw_adj <- to_num(input$stp_bw_adj, 1)
      g <- ggplot2::ggplot(df, ggplot2::aes(x = y_plot)) +
        ggplot2::geom_density(adjust = bw_adj, fill = fill_col, alpha = 0.25) +
        ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a) +
        ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = labs0$y, y = "Density")
      return(stp_apply_theme(g))
    }
    
    if (pt == "QQ plot") {
      g <- ggplot2::ggplot(df, ggplot2::aes(sample = y_plot)) +
        ggplot2::stat_qq(color = pt_col, alpha = a) +
        ggplot2::stat_qq_line(color = line_col, linewidth = lw) +
        ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = "Theoretical", y = "Sample")
      return(stp_apply_theme(g))
    }
    
    if (pt == "Lag-1 scatter") {
      y <- df$y_plot
      validate(need(length(y) >= 3, "Not enough observations for lag scatter."))
      dfl <- data.frame(x = y[-length(y)], y = y[-1])
      g <- ggplot2::ggplot(dfl, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_point(color = pt_col, size = ps, alpha = a) +
        ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = "S(t-1)", y = "S(t)")
      return(stp_apply_theme(g))
    }
    
    # (seasonal/decomposition/periodogram parts unchanged...)
    stp_apply_theme(base + ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a))
    
  }, width  = function() to_num(input$stp_plot_width_px, 980),
  height = function() to_num(input$stp_plot_height_px, 520))
  
  # ============================================================
  # 9) ACF/PACF panels
  # ============================================================
  output$stp_acf <- renderPlot({
    d <- stp_data()
    x_ts <- d$ts
    L <- min(60, length(x_ts) - 1)
    g <- forecast::ggAcf(x_ts, lag.max = L) +
      ggplot2::labs(title = "ACF", subtitle = stp_labels(d)$subtitle)
    stp_apply_theme(g)
  })
  
  output$stp_pacf <- renderPlot({
    d <- stp_data()
    x_ts <- d$ts
    L <- min(60, length(x_ts) - 1)
    g <- forecast::ggPacf(x_ts, lag.max = L) +
      ggplot2::labs(title = "PACF", subtitle = stp_labels(d)$subtitle)
    stp_apply_theme(g)
  })
  
  # ============================================================
  # 10) Lag grid (1..m)
  # ============================================================
  output$stp_lag_grid <- renderPlot({
    d <- stp_data()
    df <- d$df
    
    m <- as.integer(to_num(input$stp_lag_m, 12))
    validate(need(m >= 1, "m must be >= 1."))
    
    y <- df$y_plot
    validate(need(length(y) >= (m + 2), "Not enough observations for requested lag grid."))
    
    out <- lapply(1:m, function(k) {
      data.frame(
        lag = paste0("Lag ", k),
        x = y[seq_len(length(y) - k)],
        y = y[(k + 1):length(y)]
      )
    })
    dfl <- do.call(rbind, out)
    
    pt_col <- input$stp_point_color %||% "#2C7FB8"
    ps <- to_num(input$stp_point_size, 2)
    a  <- to_num(input$stp_alpha, 1)
    
    g <- ggplot2::ggplot(dfl, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(color = pt_col, size = ps, alpha = a) +
      ggplot2::facet_wrap(~ lag, scales = "free", ncol = 3) +
      ggplot2::labs(title = "Lag plots (1..m)", subtitle = stp_labels(d)$subtitle, x = "S(t-k)", y = "S(t)")
    
    stp_apply_theme(g)
    
  }, width  = function() to_num(input$stp_plot_width_px, 980),
  height = function() round(to_num(input$stp_plot_height_px, 520) * 1.2))
  
  
  
  
  # 
  # 
  # # ============================================================
  # # S(t) plots tab (SERVER) — paste INSIDE server()
  # # ============================================================
  # 
  # # ---- helpers ----
  # to_num <- function(x, d = NA_real_) {
  #   y <- suppressWarnings(as.numeric(x))
  #   ifelse(is.finite(y), y, d)
  # }
  # 
  # fmt_num <- function(x, d = 4) {
  #   if (!is.finite(x)) return("NA")
  #   format(round(x, d), nsmall = d)
  # }
  # 
  # # ============================================================
  # # 0) Scope control: prevent "Training only" when train_prop == 1.00
  # # ============================================================
  # output$stp_scope_ui <- renderUI({
  #   tp <- to_num(input$train_prop, 1)
  #   has_test <- isTRUE(tp < 1)
  #   
  #   if (!has_test) {
  #     radioButtons(
  #       "stp_scope",
  #       label = NULL,
  #       choices = c("Full series" = "full"),
  #       selected = "full"
  #     )
  #   } else {
  #     radioButtons(
  #       "stp_scope",
  #       label = NULL,
  #       choices = c("Full series" = "full", "Training only" = "train"),
  #       selected = input$stp_scope %||% "full"
  #     )
  #   }
  # })
  # 
  # output$stp_scope_warning <- renderUI({
  #   tp <- to_num(input$train_prop, 1)
  #   if (isTRUE(tp >= 1)) {
  #     tags$div(
  #       style = "color:#b22222; font-size:12px; margin-top:-6px;",
  #       "No test set → training = full series (Training-only is disabled)."
  #     )
  #   } else NULL
  # })
  # 
  # observeEvent(input$train_prop, {
  #   tp <- to_num(input$train_prop, 1)
  #   if (isTRUE(tp >= 1) && identical(input$stp_scope, "train")) {
  #     updateRadioButtons(session, "stp_scope", selected = "full")
  #   }
  # }, ignoreInit = TRUE)
  # 
  # # ============================================================
  # # 1) Global transform info (read-only)
  # # ============================================================
  # output$stp_transform_info <- renderPrint({
  #   tr <- input$transform %||% "none"
  #   if (tr == "boxcox") {
  #     lam <- input$lambda
  #     cat("Transformation:", "Box-Cox", "\n")
  #     cat("λ:", if (is.null(lam) || (length(lam) == 1 && is.na(lam))) "auto (estimated)" else as.character(lam), "\n")
  #   } else if (tr == "log") {
  #     cat("Transformation:", "Log (ln y)", "\n")
  #   } else {
  #     cat("Transformation:", "None", "\n")
  #   }
  # })
  # 
  # output$stp_transform_note <- renderUI({
  #   tr <- input$transform %||% "none"
  #   lam <- input$lambda
  #   
  #   tr_lbl <- switch(tr,
  #                    "none" = "None",
  #                    "log" = "Log (ln y)",
  #                    "boxcox" = "Box-Cox",
  #                    tr)
  #   
  #   if (identical(tr, "boxcox")) {
  #     tags$div(
  #       style = "font-size:12px;",
  #       tags$b("Transformation (global): "), tr_lbl, " — ",
  #       tags$b("λ: "), if (is.null(lam) || is.na(lam)) "auto" else as.character(lam)
  #     )
  #   } else {
  #     tags$div(
  #       style = "font-size:12px;",
  #       tags$b("Transformation (global): "), tr_lbl
  #     )
  #   }
  # })
  # 
  # # ============================================================
  # # 2) Theme/palette helpers  (IMPORTANT: no forced axis rotation here)
  # # ============================================================
  # stp_theme_picker <- function(key) {
  #   switch(
  #     key,
  #     "Minimal" = ggplot2::theme_minimal(),
  #     "Classic" = ggplot2::theme_classic(),
  #     "Light"   = ggplot2::theme_light(),
  #     "Dark"    = ggplot2::theme_dark(),
  #     "BW"      = ggplot2::theme_bw(),
  #     "Void"    = ggplot2::theme_void(),
  #     ggplot2::theme_gray()
  #   )
  # }
  # 
  # 
  # stp_apply_theme <- function(g) {
  #   ang <- to_num(input$stp_x_angle, 30)
  #   rot <- isTRUE(input$stp_x_rotate %||% TRUE)
  #   
  #   g +
  #     stp_theme_picker(input$stp_theme %||% "Minimal") +
  #     ggplot2::theme(
  #       text = ggplot2::element_text(size = to_num(input$stp_base_size, 12)),
  #       plot.title = ggplot2::element_text(hjust = 0.5),
  #       axis.text.x = ggplot2::element_text(
  #         angle = if (rot) ang else 0,
  #         hjust = if (!rot || ang == 0) 0.5 else 1,
  #         vjust = if (!rot || ang == 0) 0.5 else 1
  #       )
  #     )
  # }
  # 
  # 
  # 
  # # stp_apply_theme <- function(g) {
  # #   g +
  # #     stp_theme_picker(input$stp_theme %||% "Minimal") +
  # #     ggplot2::theme(
  # #       text = ggplot2::element_text(size = to_num(input$stp_base_size, 12)),
  # #       plot.title = ggplot2::element_text(hjust = 0.5)
  # #       # NOTE: do NOT set axis.text.x angle here
  # #     )
  # # }
  # 
  # stp_apply_palette <- function(g) {
  #   pal <- input$stp_palette %||% "Paired (brewer)"
  #   rev <- isTRUE(input$stp_palette_rev)
  #   
  #   if (pal == "Viridis") {
  #     g +
  #       ggplot2::scale_color_viridis_d(direction = ifelse(rev, -1, 1)) +
  #       ggplot2::scale_fill_viridis_d(direction = ifelse(rev, -1, 1))
  #   } else {
  #     brewer_name <- sub(" \\(brewer\\)$", "", pal)
  #     g +
  #       ggplot2::scale_color_brewer(palette = brewer_name, direction = ifelse(rev, -1, 1)) +
  #       ggplot2::scale_fill_brewer(palette = brewer_name, direction = ifelse(rev, -1, 1))
  #   }
  # }
  # 
  # # ============================================================
  # # 3) Conditional style UI (COLOUR PICKERS)
  # # ============================================================
  # output$stp_style_ui <- renderUI({
  #   req(input$stp_plot_type)
  #   pt <- input$stp_plot_type
  #   
  #   if (!requireNamespace("colourpicker", quietly = TRUE)) {
  #     return(tags$div(
  #       style = "color:#b22222; font-size:12px;",
  #       "Package 'colourpicker' is required for color pickers. Install it: install.packages('colourpicker')"
  #     ))
  #   }
  #   
  #   is_liney <- pt %in% c("Line", "Line + Points", "Smoothed (LOESS)", "Moving average",
  #                         "Cumulative sum", "Seasonal plot", "Seasonal subseries",
  #                         "Polar seasonal", "Periodogram",
  #                         "Classical decomposition (additive)", "Classical decomposition (multiplicative)",
  #                         "STL decomposition")
  #   
  #   is_pointy <- pt %in% c("Points", "Line + Points", "Smoothed (LOESS)",
  #                          "Lag-1 scatter", "Lag plot (1..m)", "QQ plot")
  #   
  #   is_filly  <- pt %in% c("Histogram", "Density", "Seasonal boxplot")
  #   
  #   tagList(
  #     if (is_liney) tagList(
  #       colourpicker::colourInput("stp_line_color", "Line color", value = "#2C7FB8", allowTransparent = FALSE),
  #       sliderInput("stp_line_width", "Line width", min = 0.2, max = 4, value = 1, step = 0.1)
  #     ),
  #     if (is_pointy) tagList(
  #       colourpicker::colourInput("stp_point_color", "Point color", value = "#2C7FB8", allowTransparent = FALSE),
  #       sliderInput("stp_point_size", "Point size", min = 0.5, max = 8, value = 2, step = 0.5)
  #     ),
  #     if (is_filly) tagList(
  #       colourpicker::colourInput("stp_fill_color", "Fill color", value = "#2C7FB8", allowTransparent = FALSE)
  #     )
  #   )
  # })
  # 
  # # ============================================================
  # # 4) Data reactive: choose scope + apply GLOBAL transform
  # # ============================================================
  # stp_data <- reactive({
  #   req(prepared(), ts_train_test())
  #   p <- prepared()
  #   s <- ts_train_test()
  #   
  #   df_all <- p$df
  #   req(df_all)
  #   
  #   validate(need("y_filled" %in% names(df_all), "prepared()$df must contain column y_filled."))
  #   
  #   if (!("x" %in% names(df_all))) df_all$x <- seq_len(nrow(df_all))
  #   
  #   df_all <- df_all[is.finite(df_all$y_filled), , drop = FALSE]
  #   validate(need(nrow(df_all) >= 5, "Not enough data to plot."))
  #   
  #   train_n <- s$train_n %||% floor(nrow(df_all) * to_num(input$train_prop, 1))
  #   train_n <- max(2, min(as.integer(train_n), nrow(df_all)))
  #   
  #   has_test <- isTRUE(to_num(input$train_prop, 1) < 1)
  #   scope_in <- input$stp_scope %||% "full"
  #   scope <- if (!has_test) "full" else scope_in
  #   
  #   df_use <- if (identical(scope, "train")) df_all[seq_len(train_n), , drop = FALSE] else df_all
  #   
  #   tr <- input$transform %||% "none"
  #   y <- df_use$y_filled
  #   y_plot <- y
  #   lambda_used <- NA_real_
  #   
  #   if (tr == "log") {
  #     validate(need(all(y > 0, na.rm = TRUE), "Log transform requires strictly positive values."))
  #     y_plot <- log(y)
  #   } else if (tr == "boxcox") {
  #     validate(need(all(y > 0, na.rm = TRUE), "Box-Cox transform requires strictly positive values."))
  #     lam <- input$lambda
  #     if (is.null(lam) || (length(lam) == 1 && is.na(lam))) {
  #       lam <- forecast::BoxCox.lambda(y, method = "guerrero")
  #     } else {
  #       lam <- as.numeric(lam)
  #     }
  #     lambda_used <- lam
  #     y_plot <- forecast::BoxCox(y, lam)
  #   }
  #   
  #   df_use$y_plot <- as.numeric(y_plot)
  #   
  #   freq_use <- p$freq %||% 1
  #   ts_use <- stats::ts(df_use$y_plot, start = 1, frequency = freq_use)
  #   
  #   list(
  #     df = df_use,
  #     ts = ts_use,
  #     freq = freq_use,
  #     scope = scope,
  #     train_n = train_n,
  #     n_total = nrow(df_all),
  #     transform = tr,
  #     lambda_used = lambda_used,
  #     x_is_date = inherits(df_use$x, "Date"),
  #     x_is_dt   = inherits(df_use$x, "POSIXct") || inherits(df_use$x, "POSIXt")
  #   )
  # })
  # 
  # # ============================================================
  # # 5) Labels helper (teaching subtitle)
  # # ============================================================
  # stp_labels <- function(d, default_title = "S(t)", default_y = "Value") {
  #   title_in <- input$stp_title %||% ""
  #   sub_in   <- input$stp_subtitle %||% ""
  #   x_in     <- input$stp_xlab %||% ""
  #   y_in     <- input$stp_ylab %||% ""
  #   
  #   scope_txt <- if (identical(d$scope, "train")) "training" else "full"
  #   n_txt <- paste0("n=", nrow(d$df))
  #   s_txt <- paste0("s=", d$freq)
  #   
  #   tr_txt <- switch(
  #     d$transform,
  #     "log" = "transform=log",
  #     "boxcox" = paste0("transform=BoxCox(λ=", ifelse(is.finite(d$lambda_used), fmt_num(d$lambda_used, 3), "auto"), ")"),
  #     "none" = "transform=none",
  #     "transform=none"
  #   )
  #   
  #   teaching_sub <- paste(scope_txt, n_txt, s_txt, tr_txt, sep = " • ")
  #   
  #   list(
  #     title = if (nzchar(title_in)) title_in else default_title,
  #     subtitle = if (nzchar(sub_in)) sub_in else teaching_sub,
  #     x = if (nzchar(x_in)) x_in else "Time",
  #     y = if (nzchar(y_in)) y_in else default_y
  #   )
  # }
  # 
  # # ============================================================
  # # 6) Smart date axis + ROTATION (THIS is where rotation belongs)
  # # ============================================================
  # 
  # stp_apply_x_scale <- function(g, d) {
  #   if (!isTRUE(d$x_is_date) && !isTRUE(d$x_is_dt)) return(g)
  #   
  #   n_ticks <- as.integer(to_num(input$stp_x_ticks, 8))
  #   
  #   fmt_choice <- input$stp_date_format %||% "auto"
  #   fmt_custom <- input$stp_date_format_custom %||% "%Y-%m"
  #   fmt <- if (identical(fmt_choice, "custom")) fmt_custom else fmt_choice
  #   
  #   if (identical(fmt, "auto")) {
  #     n <- nrow(d$df)
  #     if (n <= 24) fmt <- "%b %Y"
  #     else if (n <= 120) fmt <- "%Y-%m"
  #     else fmt <- "%Y"
  #   }
  #   
  #   if (isTRUE(d$x_is_date)) {
  #     g + ggplot2::scale_x_date(
  #       labels = scales::date_format(fmt),
  #       breaks = scales::pretty_breaks(n = n_ticks)
  #     )
  #   } else {
  #     g + ggplot2::scale_x_datetime(
  #       labels = scales::date_format(fmt),
  #       breaks = scales::pretty_breaks(n = n_ticks)
  #     )
  #   }
  # }
  # 
  # 
  # 
  # # stp_apply_x_scale <- function(g, d) {
  # #   if (!isTRUE(d$x_is_date) && !isTRUE(d$x_is_dt)) return(g)
  # #   
  # #   smart   <- isTRUE(input$stp_date_smart %||% TRUE)
  # #   n_ticks <- as.integer(to_num(input$stp_x_ticks, 8))
  # #   
  # #   fmt_choice <- input$stp_date_format %||% "auto"
  # #   fmt_custom <- input$stp_date_format_custom %||% "%Y-%m"
  # #   fmt <- if (identical(fmt_choice, "custom")) fmt_custom else fmt_choice
  # #   
  # #   if (identical(fmt, "auto")) {
  # #     n <- nrow(d$df)
  # #     if (n <= 24) fmt <- "%b %Y"
  # #     else if (n <= 120) fmt <- "%Y-%m"
  # #     else fmt <- "%Y"
  # #   }
  # #   
  # #   rot <- isTRUE(input$stp_x_rotate %||% TRUE)
  # #   ang <- to_num(input$stp_x_angle, 30)
  # #   
  # #   g2 <- g
  # #   
  # #   # (smart vs non-smart kept for future extension; both limit ticks now)
  # #   if (isTRUE(d$x_is_date)) {
  # #     g2 <- g2 + ggplot2::scale_x_date(
  # #       labels = scales::date_format(fmt),
  # #       breaks = scales::pretty_breaks(n = n_ticks)
  # #     )
  # #   } else {
  # #     g2 <- g2 + ggplot2::scale_x_datetime(
  # #       labels = scales::date_format(fmt),
  # #       breaks = scales::pretty_breaks(n = n_ticks)
  # #     )
  # #   }
  # #   
  # #   # ✅ rotation with correct justification
  # #   if (rot) {
  # #     g2 <- g2 + ggplot2::theme(
  # #       axis.text.x = ggplot2::element_text(
  # #         angle = ang,
  # #         hjust = ifelse(ang == 0, 0.5, 1),
  # #         vjust = ifelse(ang == 0, 0.5, 1)
  # #       )
  # #     )
  # #   } else {
  # #     g2 <- g2 + ggplot2::theme(
  # #       axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5)
  # #     )
  # #   }
  # #   
  # #   g2
  # # }
  # 
  # # ============================================================
  # # 7) UI dispatcher for multi-panel plot types
  # # ============================================================
  # output$stp_plot_ui <- renderUI({
  #   req(input$stp_plot_type)
  #   pt <- input$stp_plot_type
  #   h <- to_num(input$stp_plot_height_px, 520)
  #   
  #   if (pt == "ACF+PACF") {
  #     fluidRow(
  #       column(6, plotOutput("stp_acf",  width = "100%", height = h)),
  #       column(6, plotOutput("stp_pacf", width = "100%", height = h))
  #     )
  #   } else if (pt == "Time + ACF+PACF") {
  #     tagList(
  #       plotOutput("stp_main", width = "100%", height = round(h * 0.9)),
  #       fluidRow(
  #         column(6, plotOutput("stp_acf",  width = "100%", height = round(h * 0.8))),
  #         column(6, plotOutput("stp_pacf", width = "100%", height = round(h * 0.8)))
  #       )
  #     )
  #   } else if (pt == "Lag plot (1..m)") {
  #     plotOutput("stp_lag_grid", width = "100%", height = round(h * 1.2))
  #   } else if (pt == "ACF") {
  #     plotOutput("stp_acf", width = "100%", height = h)
  #   } else if (pt == "PACF") {
  #     plotOutput("stp_pacf", width = "100%", height = h)
  #   } else {
  #     plotOutput("stp_main", width = "100%", height = h)
  #   }
  # })
  # 
  # # ============================================================
  # # 8) Main plot
  # # ============================================================
  # output$stp_main <- renderPlot({
  #   d <- stp_data()
  #   df <- d$df
  #   pt <- input$stp_plot_type %||% "Line"
  #   
  #   line_col <- input$stp_line_color %||% "#2C7FB8"
  #   lw       <- to_num(input$stp_line_width, 1)
  #   pt_col   <- input$stp_point_color %||% line_col
  #   ps       <- to_num(input$stp_point_size, 2)
  #   fill_col <- input$stp_fill_color %||% line_col
  #   a        <- to_num(input$stp_alpha, 1)
  #   
  #   labs0 <- stp_labels(d, default_title = "S(t)", default_y = "Value")
  #   
  #   base <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y_plot)) +
  #     ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = labs0$x, y = labs0$y)
  #   
  #   # ✅ apply date scale + rotation BEFORE theme
  #   base <- stp_apply_x_scale(base, d)
  #   
  #   if (pt == "Line") {
  #     return(stp_apply_theme(base + ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a)))
  #   }
  #   
  #   if (pt == "Points") {
  #     return(stp_apply_theme(base + ggplot2::geom_point(color = pt_col, size = ps, alpha = a)))
  #   }
  #   
  #   if (pt == "Line + Points") {
  #     g <- base +
  #       ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a) +
  #       ggplot2::geom_point(color = pt_col, size = ps, alpha = a)
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Smoothed (LOESS)") {
  #     span <- to_num(input$stp_loess_span, 0.4)
  #     g <- base +
  #       ggplot2::geom_point(color = pt_col, size = ps, alpha = a) +
  #       ggplot2::geom_smooth(method = "loess", span = span, se = TRUE,
  #                            color = line_col, linewidth = lw, alpha = 0.2)
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Moving average") {
  #     k <- as.integer(to_num(input$stp_ma_k, 5))
  #     show_raw <- isTRUE(input$stp_ma_show_raw)
  #     
  #     df2 <- df
  #     df2$ma <- zoo::rollmean(df2$y_plot, k = k, fill = NA, align = "center")
  #     
  #     g <- ggplot2::ggplot(df2, ggplot2::aes(x = x)) +
  #       ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = labs0$x, y = labs0$y)
  #     
  #     g <- stp_apply_x_scale(g, d)
  #     
  #     if (show_raw) g <- g + ggplot2::geom_line(ggplot2::aes(y = y_plot), color = "gray60", linewidth = 0.6, alpha = 0.7)
  #     g <- g + ggplot2::geom_line(ggplot2::aes(y = ma), color = line_col, linewidth = lw, alpha = a)
  #     
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Cumulative sum") {
  #     df2 <- df
  #     df2$cs <- cumsum(df2$y_plot)
  #     g <- ggplot2::ggplot(df2, ggplot2::aes(x = x, y = cs)) +
  #       ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a) +
  #       ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = labs0$x, y = "Cumsum")
  #     g <- stp_apply_x_scale(g, d)
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Histogram") {
  #     bins <- as.integer(to_num(input$stp_hist_bins, 30))
  #     g <- ggplot2::ggplot(df, ggplot2::aes(x = y_plot)) +
  #       ggplot2::geom_histogram(bins = bins, fill = fill_col, alpha = a) +
  #       ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = labs0$y, y = "Count")
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Density") {
  #     bw_adj <- to_num(input$stp_bw_adj, 1)
  #     g <- ggplot2::ggplot(df, ggplot2::aes(x = y_plot)) +
  #       ggplot2::geom_density(adjust = bw_adj, fill = fill_col, alpha = 0.25) +
  #       ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a) +
  #       ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = labs0$y, y = "Density")
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "QQ plot") {
  #     g <- ggplot2::ggplot(df, ggplot2::aes(sample = y_plot)) +
  #       ggplot2::stat_qq(color = pt_col, alpha = a) +
  #       ggplot2::stat_qq_line(color = line_col, linewidth = lw) +
  #       ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = "Theoretical", y = "Sample")
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Lag-1 scatter") {
  #     y <- df$y_plot
  #     validate(need(length(y) >= 3, "Not enough observations for lag scatter."))
  #     dfl <- data.frame(x = y[-length(y)], y = y[-1])
  #     g <- ggplot2::ggplot(dfl, ggplot2::aes(x = x, y = y)) +
  #       ggplot2::geom_point(color = pt_col, size = ps, alpha = a) +
  #       ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = "S(t-1)", y = "S(t)")
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   # (seasonal/decomposition/periodogram parts unchanged...)
  #   stp_apply_theme(base + ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a))
  #   
  # }, width  = function() to_num(input$stp_plot_width_px, 980),
  # height = function() to_num(input$stp_plot_height_px, 520))
  # 
  # # ============================================================
  # # 9) ACF/PACF panels
  # # ============================================================
  # output$stp_acf <- renderPlot({
  #   d <- stp_data()
  #   x_ts <- d$ts
  #   L <- min(60, length(x_ts) - 1)
  #   g <- forecast::ggAcf(x_ts, lag.max = L) + ggplot2::labs(title = "ACF", subtitle = stp_labels(d)$subtitle)
  #   stp_apply_theme(g)
  # })
  # 
  # output$stp_pacf <- renderPlot({
  #   d <- stp_data()
  #   x_ts <- d$ts
  #   L <- min(60, length(x_ts) - 1)
  #   g <- forecast::ggPacf(x_ts, lag.max = L) + ggplot2::labs(title = "PACF", subtitle = stp_labels(d)$subtitle)
  #   stp_apply_theme(g)
  # })
  # 
  # # ============================================================
  # # 10) Lag grid (1..m)
  # # ============================================================
  # output$stp_lag_grid <- renderPlot({
  #   d <- stp_data()
  #   df <- d$df
  #   
  #   m <- as.integer(to_num(input$stp_lag_m, 12))
  #   validate(need(m >= 1, "m must be >= 1."))
  #   
  #   y <- df$y_plot
  #   validate(need(length(y) >= (m + 2), "Not enough observations for requested lag grid."))
  #   
  #   out <- lapply(1:m, function(k) {
  #     data.frame(
  #       lag = paste0("Lag ", k),
  #       x = y[seq_len(length(y) - k)],
  #       y = y[(k + 1):length(y)]
  #     )
  #   })
  #   dfl <- do.call(rbind, out)
  #   
  #   pt_col <- input$stp_point_color %||% "#2C7FB8"
  #   ps <- to_num(input$stp_point_size, 2)
  #   a  <- to_num(input$stp_alpha, 1)
  #   
  #   g <- ggplot2::ggplot(dfl, ggplot2::aes(x = x, y = y)) +
  #     ggplot2::geom_point(color = pt_col, size = ps, alpha = a) +
  #     ggplot2::facet_wrap(~ lag, scales = "free", ncol = 3) +
  #     ggplot2::labs(title = "Lag plots (1..m)", subtitle = stp_labels(d)$subtitle, x = "S(t-k)", y = "S(t)")
  #   
  #   stp_apply_theme(g)
  #   
  # }, width  = function() to_num(input$stp_plot_width_px, 980),
  # height = function() round(to_num(input$stp_plot_height_px, 520) * 1.2))
  # 
  # 
  # 
  
  
  
  

  # # ============================================================
  # # S(t) plots tab (SERVER) — paste INSIDE server()
  # # ============================================================
  # 
  # # ---- helpers ----
  # # `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
  # 
  # to_num <- function(x, d = NA_real_) {
  #   y <- suppressWarnings(as.numeric(x))
  #   ifelse(is.finite(y), y, d)
  # }
  # 
  # fmt_num <- function(x, d = 4) {
  #   if (!is.finite(x)) return("NA")
  #   format(round(x, d), nsmall = d)
  # }
  # 
  # # ============================================================
  # # 0) Scope control: prevent "Training only" when train_prop == 1.00
  # #    (Use this if your UI uses uiOutput("stp_scope_ui") + uiOutput("stp_scope_warning"))
  # # ============================================================
  # output$stp_scope_ui <- renderUI({
  #   tp <- to_num(input$train_prop, 1)
  #   has_test <- isTRUE(tp < 1)
  #   
  #   if (!has_test) {
  #     radioButtons(
  #       "stp_scope",
  #       label = NULL,
  #       choices = c("Full series" = "full"),
  #       selected = "full"
  #     )
  #   } else {
  #     radioButtons(
  #       "stp_scope",
  #       label = NULL,
  #       choices = c("Full series" = "full", "Training only" = "train"),
  #       selected = input$stp_scope %||% "full"
  #     )
  #   }
  # })
  # 
  # output$stp_scope_warning <- renderUI({
  #   tp <- to_num(input$train_prop, 1)
  #   if (isTRUE(tp >= 1)) {
  #     tags$div(
  #       style = "color:#b22222; font-size:12px; margin-top:-6px;",
  #       "No test set → training = full series (Training-only is disabled)."
  #     )
  #   } else NULL
  # })
  # 
  # observeEvent(input$train_prop, {
  #   tp <- to_num(input$train_prop, 1)
  #   if (isTRUE(tp >= 1) && identical(input$stp_scope, "train")) {
  #     updateRadioButtons(session, "stp_scope", selected = "full")
  #   }
  # }, ignoreInit = TRUE)
  # 
  # # ============================================================
  # # 1) Global transform info (read-only)
  # #    - Keep renderPrint for compatibility if you still use verbatimTextOutput("stp_transform_info")
  # #    - Add renderUI note for the newer UI (uiOutput("stp_transform_note"))
  # # ============================================================
  # output$stp_transform_info <- renderPrint({
  #   tr <- input$transform %||% "none"
  #   if (tr == "boxcox") {
  #     lam <- input$lambda
  #     cat("Transformation:", "Box-Cox", "\n")
  #     cat("λ:", if (is.null(lam) || (length(lam) == 1 && is.na(lam))) "auto (estimated)" else as.character(lam), "\n")
  #   } else if (tr == "log") {
  #     cat("Transformation:", "Log (ln y)", "\n")
  #   } else {
  #     cat("Transformation:", "None", "\n")
  #   }
  # })
  # 
  # output$stp_transform_note <- renderUI({
  #   tr <- input$transform %||% "none"
  #   lam <- input$lambda
  #   
  #   tr_lbl <- switch(tr,
  #                    "none" = "None",
  #                    "log" = "Log (ln y)",
  #                    "boxcox" = "Box-Cox",
  #                    tr)
  #   
  #   if (identical(tr, "boxcox")) {
  #     tags$div(
  #       style = "font-size:12px;",
  #       tags$b("Transformation (global): "), tr_lbl, " — ",
  #       tags$b("λ: "), if (is.null(lam) || is.na(lam)) "auto" else as.character(lam)
  #     )
  #   } else {
  #     tags$div(
  #       style = "font-size:12px;",
  #       tags$b("Transformation (global): "), tr_lbl
  #     )
  #   }
  # })
  # 
  # # ============================================================
  # # 2) Theme/palette helpers
  # # ============================================================
  # stp_theme_picker <- function(key) {
  #   switch(
  #     key,
  #     "Minimal" = ggplot2::theme_minimal(),
  #     "Classic" = ggplot2::theme_classic(),
  #     "Light"   = ggplot2::theme_light(),
  #     "Dark"    = ggplot2::theme_dark(),
  #     "BW"      = ggplot2::theme_bw(),
  #     "Void"    = ggplot2::theme_void(),
  #     ggplot2::theme_gray()
  #   )
  # }
  # 
  # stp_apply_theme <- function(g) {
  #   g +
  #     stp_theme_picker(input$stp_theme %||% "Minimal") +
  #     ggplot2::theme(
  #       text = ggplot2::element_text(size = to_num(input$stp_base_size, 12)),
  #       plot.title = ggplot2::element_text(hjust = 0.5),
  #       axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1)
  #     )
  # }
  # 
  # stp_apply_palette <- function(g) {
  #   pal <- input$stp_palette %||% "Paired (brewer)"
  #   rev <- isTRUE(input$stp_palette_rev)
  #   
  #   if (pal == "Viridis") {
  #     g +
  #       ggplot2::scale_color_viridis_d(direction = ifelse(rev, -1, 1)) +
  #       ggplot2::scale_fill_viridis_d(direction = ifelse(rev, -1, 1))
  #   } else {
  #     brewer_name <- sub(" \\(brewer\\)$", "", pal)
  #     g +
  #       ggplot2::scale_color_brewer(palette = brewer_name, direction = ifelse(rev, -1, 1)) +
  #       ggplot2::scale_fill_brewer(palette = brewer_name, direction = ifelse(rev, -1, 1))
  #   }
  # }
  # 
  # # ============================================================
  # # 3) Conditional style UI (COLOUR PICKERS)
  # # ============================================================
  # output$stp_style_ui <- renderUI({
  #   req(input$stp_plot_type)
  #   pt <- input$stp_plot_type
  #   
  #   if (!requireNamespace("colourpicker", quietly = TRUE)) {
  #     return(tags$div(
  #       style = "color:#b22222; font-size:12px;",
  #       "Package 'colourpicker' is required for color pickers. Install it: install.packages('colourpicker')"
  #     ))
  #   }
  #   
  #   is_liney <- pt %in% c("Line", "Line + Points", "Smoothed (LOESS)", "Moving average",
  #                         "Cumulative sum", "Seasonal plot", "Seasonal subseries",
  #                         "Polar seasonal", "Periodogram",
  #                         "Classical decomposition (additive)", "Classical decomposition (multiplicative)",
  #                         "STL decomposition")
  #   
  #   is_pointy <- pt %in% c("Points", "Line + Points", "Smoothed (LOESS)",
  #                          "Lag-1 scatter", "Lag plot (1..m)", "QQ plot")
  #   
  #   is_filly  <- pt %in% c("Histogram", "Density", "Seasonal boxplot")
  #   
  #   tagList(
  #     if (is_liney) tagList(
  #       colourpicker::colourInput("stp_line_color", "Line color", value = "#2C7FB8", allowTransparent = FALSE),
  #       sliderInput("stp_line_width", "Line width", min = 0.2, max = 4, value = 1, step = 0.1)
  #     ),
  #     if (is_pointy) tagList(
  #       colourpicker::colourInput("stp_point_color", "Point color", value = "#2C7FB8", allowTransparent = FALSE),
  #       sliderInput("stp_point_size", "Point size", min = 0.5, max = 8, value = 2, step = 0.5)
  #     ),
  #     if (is_filly) tagList(
  #       colourpicker::colourInput("stp_fill_color", "Fill color", value = "#2C7FB8", allowTransparent = FALSE)
  #     )
  #   )
  # })
  # 
  # # ============================================================
  # # 4) Data reactive: choose scope + apply GLOBAL transform (none/log/boxcox)
  # #    Also keeps Date axis correctly if df$x is Date/POSIXct.
  # # ============================================================
  # stp_data <- reactive({
  #   req(prepared(), ts_train_test())
  #   p <- prepared()
  #   s <- ts_train_test()
  #   
  #   df_all <- p$df
  #   req(df_all)
  #   
  #   validate(need("y_filled" %in% names(df_all), "prepared()$df must contain column y_filled."))
  #   
  #   # Ensure x exists
  #   if (!("x" %in% names(df_all))) df_all$x <- seq_len(nrow(df_all))
  #   
  #   df_all <- df_all[is.finite(df_all$y_filled), , drop = FALSE]
  #   validate(need(nrow(df_all) >= 5, "Not enough data to plot."))
  #   
  #   # train_n
  #   train_n <- s$train_n %||% floor(nrow(df_all) * to_num(input$train_prop, 1))
  #   train_n <- max(2, min(as.integer(train_n), nrow(df_all)))
  #   
  #   # scope (force full if no split)
  #   has_test <- isTRUE(to_num(input$train_prop, 1) < 1)
  #   scope_in <- input$stp_scope %||% "full"
  #   scope <- if (!has_test) "full" else scope_in
  #   
  #   df_use <- if (identical(scope, "train")) df_all[seq_len(train_n), , drop = FALSE] else df_all
  #   
  #   # Apply GLOBAL transformation to plotted y
  #   tr <- input$transform %||% "none"
  #   y <- df_use$y_filled
  #   y_plot <- y
  #   lambda_used <- NA_real_
  #   
  #   if (tr == "log") {
  #     validate(need(all(y > 0, na.rm = TRUE), "Log transform requires strictly positive values."))
  #     y_plot <- log(y)
  #   } else if (tr == "boxcox") {
  #     validate(need(all(y > 0, na.rm = TRUE), "Box-Cox transform requires strictly positive values."))
  #     lam <- input$lambda
  #     if (is.null(lam) || (length(lam) == 1 && is.na(lam))) {
  #       lam <- forecast::BoxCox.lambda(y, method = "guerrero")
  #     } else {
  #       lam <- as.numeric(lam)
  #     }
  #     lambda_used <- lam
  #     y_plot <- forecast::BoxCox(y, lam)
  #   }
  #   
  #   df_use$y_plot <- as.numeric(y_plot)
  #   
  #   freq_use <- p$freq %||% 1
  #   ts_use <- stats::ts(df_use$y_plot, start = 1, frequency = freq_use)
  #   
  #   list(
  #     df = df_use,
  #     ts = ts_use,
  #     freq = freq_use,
  #     scope = scope,
  #     train_n = train_n,
  #     n_total = nrow(df_all),
  #     transform = tr,
  #     lambda_used = lambda_used,
  #     x_is_date = inherits(df_use$x, "Date"),
  #     x_is_dt   = inherits(df_use$x, "POSIXct") || inherits(df_use$x, "POSIXt")
  #   )
  # })
  # 
  # # ============================================================
  # # 5) Labels helper + teaching subtitle (n, s, scope, transform)
  # # ============================================================
  # stp_labels <- function(d, default_title = "S(t)", default_y = "Value") {
  #   title_in <- input$stp_title %||% ""
  #   sub_in   <- input$stp_subtitle %||% ""
  #   x_in     <- input$stp_xlab %||% ""
  #   y_in     <- input$stp_ylab %||% ""
  #   
  #   scope_txt <- if (identical(d$scope, "train")) "training" else "full"
  #   n_txt <- paste0("n=", nrow(d$df))
  #   s_txt <- paste0("s=", d$freq)
  #   
  #   tr_txt <- switch(
  #     d$transform,
  #     "log" = "transform=log",
  #     "boxcox" = paste0("transform=BoxCox(λ=", ifelse(is.finite(d$lambda_used), fmt_num(d$lambda_used, 3), "auto"), ")"),
  #     "none" = "transform=none",
  #     "transform=none"
  #   )
  #   
  #   teaching_sub <- paste(scope_txt, n_txt, s_txt, tr_txt, sep = " • ")
  #   
  #   list(
  #     title = if (nzchar(title_in)) title_in else default_title,
  #     subtitle = if (nzchar(sub_in)) sub_in else teaching_sub,
  #     x = if (nzchar(x_in)) x_in else "Time",
  #     y = if (nzchar(y_in)) y_in else default_y
  #   )
  # }
  # 
  # # ============================================================
  # # 6) Smart date axis formatting
  # # ============================================================
  # 
  # stp_apply_x_scale <- function(g, d) {
  #   # If not date/datetime, do nothing
  #   if (!isTRUE(d$x_is_date) && !isTRUE(d$x_is_dt)) return(g)
  #   
  #   # user controls
  #   smart <- isTRUE(input$stp_date_smart %||% TRUE)
  #   n_ticks <- as.integer(to_num(input$stp_x_ticks, 8))
  #   
  #   # label format
  #   fmt_choice <- input$stp_date_format %||% "auto"
  #   fmt_custom <- input$stp_date_format_custom %||% "%Y-%m"
  #   
  #   fmt <- if (identical(fmt_choice, "custom")) fmt_custom else fmt_choice
  #   
  #   # if Auto: pick something reasonable
  #   if (identical(fmt, "auto")) {
  #     # heuristic based on how many rows are shown
  #     n <- nrow(d$df)
  #     if (n <= 24) fmt <- "%b %Y"
  #     else if (n <= 120) fmt <- "%Y-%m"
  #     else fmt <- "%Y"
  #   }
  #   
  #   # rotation controls
  #   rot <- isTRUE(input$stp_x_rotate %||% TRUE)
  #   ang <- to_num(input$stp_x_angle, 30)
  #   
  #   g2 <- g
  #   
  #   # Use "smart" breaks = pretty breaks (gives ~n ticks)
  #   if (smart) {
  #     if (isTRUE(d$x_is_date)) {
  #       g2 <- g2 + ggplot2::scale_x_date(
  #         labels = scales::date_format(fmt),
  #         breaks = scales::pretty_breaks(n = n_ticks)
  #       )
  #     } else {
  #       g2 <- g2 + ggplot2::scale_x_datetime(
  #         labels = scales::date_format(fmt),
  #         breaks = scales::pretty_breaks(n = n_ticks)
  #       )
  #     }
  #   } else {
  #     # Non-smart: still limit to pretty breaks, but fewer ticks
  #     if (isTRUE(d$x_is_date)) {
  #       g2 <- g2 + ggplot2::scale_x_date(
  #         labels = scales::date_format(fmt),
  #         breaks = scales::pretty_breaks(n = n_ticks)
  #       )
  #     } else {
  #       g2 <- g2 + ggplot2::scale_x_datetime(
  #         labels = scales::date_format(fmt),
  #         breaks = scales::pretty_breaks(n = n_ticks)
  #       )
  #     }
  #   }
  #   
  #   # apply rotation consistently (override theme angle if desired)
  #   if (rot) {
  #     g2 <- g2 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = ang, hjust = 1, vjust = 1))
  #   } else {
  #     g2 <- g2 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5))
  #   }
  #   
  #   g2
  # }
  # 
  # 
  # 
  # 
  # # stp_apply_x_scale <- function(g, d) {
  # #   if (isTRUE(d$x_is_date)) {
  # #     g + ggplot2::scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months")
  # #   } else if (isTRUE(d$x_is_dt)) {
  # #     g + ggplot2::scale_x_datetime(date_labels = "%Y-%m-%d", date_breaks = "1 month")
  # #   } else {
  # #     g
  # #   }
  # # }
  # 
  # # ============================================================
  # # 7) UI dispatcher for multi-panel plot types
  # # ============================================================
  # output$stp_plot_ui <- renderUI({
  #   req(input$stp_plot_type)
  #   pt <- input$stp_plot_type
  #   h <- to_num(input$stp_plot_height_px, 520)
  #   
  #   if (pt == "ACF+PACF") {
  #     fluidRow(
  #       column(6, plotOutput("stp_acf",  width = "100%", height = h)),
  #       column(6, plotOutput("stp_pacf", width = "100%", height = h))
  #     )
  #   } else if (pt == "Time + ACF+PACF") {
  #     tagList(
  #       plotOutput("stp_main", width = "100%", height = round(h * 0.9)),
  #       fluidRow(
  #         column(6, plotOutput("stp_acf",  width = "100%", height = round(h * 0.8))),
  #         column(6, plotOutput("stp_pacf", width = "100%", height = round(h * 0.8)))
  #       )
  #     )
  #   } else if (pt == "Lag plot (1..m)") {
  #     plotOutput("stp_lag_grid", width = "100%", height = round(h * 1.2))
  #   } else if (pt == "ACF") {
  #     plotOutput("stp_acf", width = "100%", height = h)
  #   } else if (pt == "PACF") {
  #     plotOutput("stp_pacf", width = "100%", height = h)
  #   } else {
  #     plotOutput("stp_main", width = "100%", height = h)
  #   }
  # })
  # 
  # # ============================================================
  # # 8) Main plot
  # # ============================================================
  # output$stp_main <- renderPlot({
  #   d <- stp_data()
  #   df <- d$df
  #   pt <- input$stp_plot_type %||% "Line"
  #   
  #   # Style (fallback-safe)
  #   line_col <- input$stp_line_color %||% "#2C7FB8"
  #   lw       <- to_num(input$stp_line_width, 1)
  #   pt_col   <- input$stp_point_color %||% line_col
  #   ps       <- to_num(input$stp_point_size, 2)
  #   fill_col <- input$stp_fill_color %||% line_col
  #   a        <- to_num(input$stp_alpha, 1)
  #   
  #   labs0 <- stp_labels(d,
  #                       default_title = "S(t)",
  #                       default_y = "Value")
  #   
  #   base <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y_plot)) +
  #     ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = labs0$x, y = labs0$y)
  #   
  #   base <- stp_apply_x_scale(base, d)
  #   
  #   # ---- plot switch ----
  #   if (pt == "Line") {
  #     g <- base + ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a)
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Points") {
  #     g <- base + ggplot2::geom_point(color = pt_col, size = ps, alpha = a)
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Line + Points") {
  #     g <- base +
  #       ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a) +
  #       ggplot2::geom_point(color = pt_col, size = ps, alpha = a)
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Smoothed (LOESS)") {
  #     span <- to_num(input$stp_loess_span, 0.4)
  #     g <- base +
  #       ggplot2::geom_point(color = pt_col, size = ps, alpha = a) +
  #       ggplot2::geom_smooth(method = "loess", span = span, se = TRUE,
  #                            color = line_col, linewidth = lw, alpha = 0.2)
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Moving average") {
  #     k <- as.integer(to_num(input$stp_ma_k, 5))
  #     show_raw <- isTRUE(input$stp_ma_show_raw)
  #     
  #     df2 <- df
  #     df2$ma <- zoo::rollmean(df2$y_plot, k = k, fill = NA, align = "center")
  #     
  #     g <- ggplot2::ggplot(df2, ggplot2::aes(x = x)) +
  #       ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = labs0$x, y = labs0$y)
  #     
  #     g <- stp_apply_x_scale(g, d)
  #     
  #     if (show_raw) g <- g + ggplot2::geom_line(ggplot2::aes(y = y_plot), color = "gray60", linewidth = 0.6, alpha = 0.7)
  #     g <- g + ggplot2::geom_line(ggplot2::aes(y = ma), color = line_col, linewidth = lw, alpha = a)
  #     
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Cumulative sum") {
  #     df2 <- df
  #     df2$cs <- cumsum(df2$y_plot)
  #     g <- ggplot2::ggplot(df2, ggplot2::aes(x = x, y = cs)) +
  #       ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a) +
  #       ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = labs0$x, y = "Cumsum")
  #     g <- stp_apply_x_scale(g, d)
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Histogram") {
  #     bins <- as.integer(to_num(input$stp_hist_bins, 30))
  #     g <- ggplot2::ggplot(df, ggplot2::aes(x = y_plot)) +
  #       ggplot2::geom_histogram(bins = bins, fill = fill_col, alpha = a) +
  #       ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = labs0$y, y = "Count")
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Density") {
  #     bw_adj <- to_num(input$stp_bw_adj, 1)
  #     g <- ggplot2::ggplot(df, ggplot2::aes(x = y_plot)) +
  #       ggplot2::geom_density(adjust = bw_adj, fill = fill_col, alpha = 0.25) +
  #       ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a) +
  #       ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = labs0$y, y = "Density")
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "QQ plot") {
  #     g <- ggplot2::ggplot(df, ggplot2::aes(sample = y_plot)) +
  #       ggplot2::stat_qq(color = pt_col, alpha = a) +
  #       ggplot2::stat_qq_line(color = line_col, linewidth = lw) +
  #       ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = "Theoretical", y = "Sample")
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Lag-1 scatter") {
  #     y <- df$y_plot
  #     validate(need(length(y) >= 3, "Not enough observations for lag scatter."))
  #     dfl <- data.frame(x = y[-length(y)], y = y[-1])
  #     g <- ggplot2::ggplot(dfl, ggplot2::aes(x = x, y = y)) +
  #       ggplot2::geom_point(color = pt_col, size = ps, alpha = a) +
  #       ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = "S(t-1)", y = "S(t)")
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt %in% c("Seasonal plot", "Seasonal subseries", "Polar seasonal", "Seasonal boxplot")) {
  #     validate(need(d$freq >= 2, "Seasonal plots require frequency >= 2."))
  #     
  #     x_ts <- d$ts
  #     
  #     if (pt == "Seasonal plot") {
  #       g <- forecast::ggseasonplot(x_ts, polar = FALSE) +
  #         ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle)
  #       return(stp_apply_theme(stp_apply_palette(g)))
  #     }
  #     
  #     if (pt == "Polar seasonal") {
  #       g <- forecast::ggseasonplot(x_ts, polar = TRUE) +
  #         ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle)
  #       return(stp_apply_theme(stp_apply_palette(g)))
  #     }
  #     
  #     if (pt == "Seasonal subseries") {
  #       g <- forecast::ggsubseriesplot(x_ts) +
  #         ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle)
  #       return(stp_apply_theme(g))
  #     }
  #     
  #     if (pt == "Seasonal boxplot") {
  #       df2 <- df
  #       df2$season <- as.factor(cycle(stats::ts(df2$y_plot, frequency = d$freq)))
  #       g <- ggplot2::ggplot(df2, ggplot2::aes(x = season, y = y_plot, fill = season)) +
  #         ggplot2::geom_boxplot(alpha = 0.6) +
  #         ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = "Season", y = labs0$y)
  #       return(stp_apply_theme(stp_apply_palette(g)))
  #     }
  #   }
  #   
  #   if (pt == "Classical decomposition (additive)" || pt == "Classical decomposition (multiplicative)") {
  #     validate(need(d$freq >= 2, "Decomposition requires frequency >= 2."))
  #     x_ts <- d$ts
  #     type <- if (pt == "Classical decomposition (multiplicative)") "multiplicative" else "additive"
  #     dc <- stats::decompose(x_ts, type = type)
  #     g <- forecast::autoplot(dc) + ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle)
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "STL decomposition") {
  #     validate(need(d$freq >= 2, "STL requires frequency >= 2."))
  #     x_ts <- d$ts
  #     fit <- stats::stl(x_ts, s.window = "periodic", robust = TRUE)
  #     g <- forecast::autoplot(fit) + ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle)
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Periodogram") {
  #     x_ts <- d$ts
  #     validate(need(length(x_ts) >= 8, "Not enough observations for periodogram."))
  #     taper <- to_num(input$stp_spec_taper, 0.1)
  #     sp <- stats::spec.pgram(x_ts, taper = taper, plot = FALSE)
  #     dfp <- data.frame(freq = sp$freq, spec = sp$spec)
  #     g <- ggplot2::ggplot(dfp, ggplot2::aes(x = freq, y = spec)) +
  #       ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a) +
  #       ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle, x = "Frequency", y = "Spectral density")
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   # fallback
  #   stp_apply_theme(base + ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a))
  #   
  # }, width  = function() to_num(input$stp_plot_width_px, 980),
  # height = function() to_num(input$stp_plot_height_px, 520))
  # 
  # # ============================================================
  # # 9) ACF/PACF panels
  # # ============================================================
  # output$stp_acf <- renderPlot({
  #   d <- stp_data()
  #   x_ts <- d$ts
  #   L <- min(60, length(x_ts) - 1)
  #   g <- forecast::ggAcf(x_ts, lag.max = L) + ggplot2::labs(title = "ACF", subtitle = stp_labels(d)$subtitle)
  #   stp_apply_theme(g)
  # })
  # 
  # output$stp_pacf <- renderPlot({
  #   d <- stp_data()
  #   x_ts <- d$ts
  #   L <- min(60, length(x_ts) - 1)
  #   g <- forecast::ggPacf(x_ts, lag.max = L) + ggplot2::labs(title = "PACF", subtitle = stp_labels(d)$subtitle)
  #   stp_apply_theme(g)
  # })
  # 
  # # ============================================================
  # # 10) Lag grid (1..m)
  # # ============================================================
  # output$stp_lag_grid <- renderPlot({
  #   d <- stp_data()
  #   df <- d$df
  #   
  #   m <- as.integer(to_num(input$stp_lag_m, 12))
  #   validate(need(m >= 1, "m must be >= 1."))
  #   
  #   y <- df$y_plot
  #   validate(need(length(y) >= (m + 2), "Not enough observations for requested lag grid."))
  #   
  #   out <- lapply(1:m, function(k) {
  #     data.frame(
  #       lag = paste0("Lag ", k),
  #       x = y[seq_len(length(y) - k)],
  #       y = y[(k + 1):length(y)]
  #     )
  #   })
  #   dfl <- do.call(rbind, out)
  #   
  #   pt_col <- input$stp_point_color %||% "#2C7FB8"
  #   ps <- to_num(input$stp_point_size, 2)
  #   a  <- to_num(input$stp_alpha, 1)
  #   
  #   g <- ggplot2::ggplot(dfl, ggplot2::aes(x = x, y = y)) +
  #     ggplot2::geom_point(color = pt_col, size = ps, alpha = a) +
  #     ggplot2::facet_wrap(~ lag, scales = "free", ncol = 3) +
  #     ggplot2::labs(title = "Lag plots (1..m)", subtitle = stp_labels(d)$subtitle, x = "S(t-k)", y = "S(t)")
  #   
  #   stp_apply_theme(g)
  #   
  # }, width  = function() to_num(input$stp_plot_width_px, 980),
  # height = function() round(to_num(input$stp_plot_height_px, 520) * 1.2))
  # 
  
  
  
  
  
  
  
  
  
  
  
  # ============================================================
  # S(t) plots tab (SERVER) — MUST BE INSIDE server()
  # Uses GLOBAL input$transform and input$lambda
  # ============================================================
  
  # `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
  
  # to_num <- function(x, d = NA_real_) {
  #   y <- suppressWarnings(as.numeric(x))
  #   ifelse(is.finite(y), y, d)
  # }
  # 
  # fmt_num <- function(x, d = 4) {
  #   if (!is.finite(x)) return("NA")
  #   format(round(x, d), nsmall = d)
  # }
  # 
  # # ---- show global transform info (read-only) ----
  # output$stp_transform_info <- renderPrint({
  #   tr <- input$transform %||% "none"
  #   if (tr == "boxcox") {
  #     lam <- input$lambda
  #     cat("Transformation:", "Box-Cox", "\n")
  #     cat("λ:", if (is.null(lam) || (length(lam) == 1 && is.na(lam))) "auto (estimated)" else as.character(lam), "\n")
  #   } else if (tr == "log") {
  #     cat("Transformation:", "Log (ln y)", "\n")
  #   } else {
  #     cat("Transformation:", "None", "\n")
  #   }
  # })
  # 
  # # ---- Theme/palette helpers ----
  # stp_theme_picker <- function(key) {
  #   switch(
  #     key,
  #     "Minimal" = ggplot2::theme_minimal(),
  #     "Classic" = ggplot2::theme_classic(),
  #     "Light"   = ggplot2::theme_light(),
  #     "Dark"    = ggplot2::theme_dark(),
  #     "BW"      = ggplot2::theme_bw(),
  #     "Void"    = ggplot2::theme_void(),
  #     ggplot2::theme_gray()
  #   )
  # }
  # 
  # stp_apply_theme <- function(g) {
  #   g +
  #     stp_theme_picker(input$stp_theme %||% "Minimal") +
  #     ggplot2::theme(
  #       text = ggplot2::element_text(size = as.numeric(input$stp_base_size %||% 12)),
  #       plot.title = ggplot2::element_text(hjust = 0.5)
  #     )
  # }
  # 
  # stp_apply_palette <- function(g) {
  #   pal <- input$stp_palette %||% "Set1 (brewer)"
  #   rev <- isTRUE(input$stp_palette_rev)
  #   
  #   if (pal == "Viridis") {
  #     g +
  #       ggplot2::scale_color_viridis_d(direction = ifelse(rev, -1, 1)) +
  #       ggplot2::scale_fill_viridis_d(direction = ifelse(rev, -1, 1))
  #   } else {
  #     brewer_name <- sub(" \\(brewer\\)$", "", pal)
  #     g +
  #       ggplot2::scale_color_brewer(palette = brewer_name, direction = ifelse(rev, -1, 1)) +
  #       ggplot2::scale_fill_brewer(palette = brewer_name, direction = ifelse(rev, -1, 1))
  #   }
  # }
  # 
  # # ---- conditional style UI ----
  # output$stp_style_ui <- renderUI({
  #   req(input$stp_plot_type)
  #   
  #   pt <- input$stp_plot_type
  #   
  #   is_liney <- pt %in% c("Line", "Line + Points", "Smoothed (LOESS)", "Moving average",
  #                         "Cumulative sum", "Seasonal plot", "Seasonal subseries",
  #                         "Polar seasonal", "Periodogram",
  #                         "Classical decomposition (additive)", "Classical decomposition (multiplicative)",
  #                         "STL decomposition")
  #   is_pointy <- pt %in% c("Points", "Line + Points", "Smoothed (LOESS)",
  #                          "Lag-1 scatter", "Lag plot (1..m)", "QQ plot")
  #   is_filly  <- pt %in% c("Histogram", "Density", "Seasonal boxplot")
  #   
  #   tagList(
  #     # Line controls
  #     if (is_liney) tagList(
  #       colourpicker::colourInput("stp_line_color", "Line color", value = "#2C7FB8", allowTransparent = FALSE),
  #       sliderInput("stp_line_width", "Line width", min = 0.2, max = 4, value = 1, step = 0.1)
  #     ),
  #     
  #     # Point controls
  #     if (is_pointy) tagList(
  #       colourpicker::colourInput("stp_point_color", "Point color", value = "#2C7FB8", allowTransparent = FALSE),
  #       sliderInput("stp_point_size", "Point size", min = 0.5, max = 8, value = 2, step = 0.5)
  #     ),
  #     
  #     # Fill controls
  #     if (is_filly) tagList(
  #       colourpicker::colourInput("stp_fill_color", "Fill color", value = "#2C7FB8", allowTransparent = FALSE)
  #     )
  #   )
  # })
  # 
  # # ---- choose series scope + apply GLOBAL transform ----
  # stp_data <- reactive({
  #   req(prepared(), ts_train_test())
  #   p <- prepared()
  #   s <- ts_train_test()
  #   
  #   df_all <- p$df
  #   req(df_all)
  #   
  #   # Expect these columns: y_filled and x (if x missing, create index)
  #   validate(need("y_filled" %in% names(df_all), "prepared()$df must contain y_filled."))
  #   
  #   if (!("x" %in% names(df_all))) df_all$x <- seq_len(nrow(df_all))
  #   
  #   df_all <- df_all[is.finite(df_all$y_filled), , drop = FALSE]
  #   validate(need(nrow(df_all) >= 5, "Not enough data to plot."))
  #   
  #   train_n <- s$train_n %||% floor(nrow(df_all) * as.numeric(input$train_prop %||% 1))
  #   train_n <- max(2, min(as.integer(train_n), nrow(df_all)))
  #   
  #   scope <- input$stp_scope %||% "full"
  #   df_use <- if (identical(scope, "train")) df_all[seq_len(train_n), , drop = FALSE] else df_all
  #   
  #   # Apply GLOBAL transformation (None / Log / BoxCox) to the plotted y
  #   tr <- input$transform %||% "none"
  #   y <- df_use$y_filled
  #   
  #   y_plot <- y
  #   if (tr == "log") {
  #     # safeguard: log needs positive
  #     validate(need(all(y > 0, na.rm = TRUE), "Log transform requires strictly positive values."))
  #     y_plot <- log(y)
  #   } else if (tr == "boxcox") {
  #     # BoxCox requires positive
  #     validate(need(all(y > 0, na.rm = TRUE), "Box-Cox transform requires strictly positive values."))
  #     lam <- input$lambda
  #     if (is.null(lam) || (length(lam) == 1 && is.na(lam))) {
  #       lam <- forecast::BoxCox.lambda(y, method = "guerrero")
  #     } else {
  #       lam <- as.numeric(lam)
  #     }
  #     y_plot <- forecast::BoxCox(y, lam)
  #     attr(y_plot, "lambda_used") <- lam
  #   }
  #   
  #   df_use$y_plot <- as.numeric(y_plot)
  #   
  #   freq_use <- p$freq %||% 1
  #   ts_use <- stats::ts(df_use$y_plot, start = 1, frequency = freq_use)
  #   
  #   list(
  #     df = df_use,
  #     ts = ts_use,
  #     freq = freq_use,
  #     scope = scope,
  #     train_n = train_n,
  #     n_total = nrow(df_all),
  #     transform = tr,
  #     lambda_used = attr(y_plot, "lambda_used") %||% NA_real_
  #   )
  # })
  # 
  # # ---- label helper ----
  # stp_labels <- function(default_title = "S(t)", default_y = "Value") {
  #   title_in <- input$stp_title %||% ""
  #   sub_in   <- input$stp_subtitle %||% ""
  #   x_in     <- input$stp_xlab %||% ""
  #   y_in     <- input$stp_ylab %||% ""
  #   
  #   list(
  #     title = if (nzchar(title_in)) title_in else default_title,
  #     subtitle = if (nzchar(sub_in)) sub_in else NULL,
  #     x = if (nzchar(x_in)) x_in else "Time",
  #     y = if (nzchar(y_in)) y_in else default_y
  #   )
  # }
  # 
  # # ---- UI dispatcher for multi-panel plot types ----
  # output$stp_plot_ui <- renderUI({
  #   req(input$stp_plot_type)
  #   pt <- input$stp_plot_type
  #   h <- as.numeric(input$stp_plot_height_px %||% 520)
  #   
  #   if (pt == "ACF+PACF") {
  #     fluidRow(
  #       column(6, plotOutput("stp_acf",  width = "100%", height = h)),
  #       column(6, plotOutput("stp_pacf", width = "100%", height = h))
  #     )
  #   } else if (pt == "Time + ACF+PACF") {
  #     tagList(
  #       plotOutput("stp_main", width = "100%", height = round(h * 0.9)),
  #       fluidRow(
  #         column(6, plotOutput("stp_acf",  width = "100%", height = round(h * 0.8))),
  #         column(6, plotOutput("stp_pacf", width = "100%", height = round(h * 0.8)))
  #       )
  #     )
  #   } else if (pt == "Lag plot (1..m)") {
  #     plotOutput("stp_lag_grid", width = "100%", height = round(h * 1.2))
  #   } else if (pt == "ACF") {
  #     plotOutput("stp_acf", width = "100%", height = h)
  #   } else if (pt == "PACF") {
  #     plotOutput("stp_pacf", width = "100%", height = h)
  #   } else {
  #     plotOutput("stp_main", width = "100%", height = h)
  #   }
  # })
  # 
  # # ---- Main plot ----
  # output$stp_main <- renderPlot({
  #   d <- stp_data()
  #   df <- d$df
  #   scope_tag <- if (identical(d$scope, "train")) " (Training only)" else " (Full series)"
  #   
  #   pt <- input$stp_plot_type %||% "Line"
  #   
  #   # Style values (exist only if UI displayed them; use fallbacks)
  #   line_col <- input$stp_line_color %||% "#2C7FB8"
  #   lw <- as.numeric(input$stp_line_width %||% 1)
  #   pt_col <- input$stp_point_color %||% line_col
  #   ps <- as.numeric(input$stp_point_size %||% 2)
  #   fill_col <- input$stp_fill_color %||% line_col
  #   a  <- as.numeric(input$stp_alpha %||% 1)
  #   
  #   # teaching subtitle: transformation note (auto)
  #   tr_note <- switch(
  #     d$transform,
  #     "log" = "Log-transformed",
  #     "boxcox" = paste0("Box-Cox (λ=", ifelse(is.finite(d$lambda_used), fmt_num(d$lambda_used, 3), "auto"), ")"),
  #     "none" = "No transformation",
  #     "No transformation"
  #   )
  #   
  #   # labels
  #   labs0 <- stp_labels(
  #     default_title = paste0("S(t)", scope_tag),
  #     default_y = paste0("S(t) [", tr_note, "]")
  #   )
  #   
  #   base <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y_plot)) +
  #     ggplot2::labs(title = labs0$title, subtitle = labs0$subtitle %||% tr_note, x = labs0$x, y = labs0$y)
  #   
  #   # ---- plot switch ----
  #   if (pt == "Line") {
  #     return(stp_apply_theme(base + ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a)))
  #   }
  #   
  #   if (pt == "Points") {
  #     return(stp_apply_theme(base + ggplot2::geom_point(color = pt_col, size = ps, alpha = a)))
  #   }
  #   
  #   if (pt == "Line + Points") {
  #     return(stp_apply_theme(
  #       base +
  #         ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a) +
  #         ggplot2::geom_point(color = pt_col, size = ps, alpha = a)
  #     ))
  #   }
  #   
  #   if (pt == "Smoothed (LOESS)") {
  #     span <- as.numeric(input$stp_loess_span %||% 0.4)
  #     return(stp_apply_theme(
  #       base +
  #         ggplot2::geom_point(color = pt_col, size = ps, alpha = a) +
  #         ggplot2::geom_smooth(method = "loess", span = span, se = TRUE,
  #                              color = line_col, linewidth = lw, alpha = 0.2)
  #     ))
  #   }
  #   
  #   if (pt == "Moving average") {
  #     k <- as.integer(input$stp_ma_k %||% 5)
  #     show_raw <- isTRUE(input$stp_ma_show_raw)
  #     
  #     df2 <- df
  #     df2$ma <- zoo::rollmean(df2$y_plot, k = k, fill = NA, align = "center")
  #     
  #     g <- ggplot2::ggplot(df2, ggplot2::aes(x = x)) +
  #       ggplot2::labs(
  #         title = if (nzchar(input$stp_title %||% "")) labs0$title else paste0("Moving average (k=", k, ")", scope_tag),
  #         subtitle = labs0$subtitle %||% tr_note,
  #         x = labs0$x, y = labs0$y
  #       )
  #     
  #     if (show_raw) g <- g + ggplot2::geom_line(ggplot2::aes(y = y_plot), color = "gray60", linewidth = 0.6, alpha = 0.7)
  #     g <- g + ggplot2::geom_line(ggplot2::aes(y = ma), color = line_col, linewidth = lw, alpha = a)
  #     
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Cumulative sum") {
  #     df2 <- df
  #     df2$cs <- cumsum(df2$y_plot)
  #     g <- ggplot2::ggplot(df2, ggplot2::aes(x = x, y = cs)) +
  #       ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a) +
  #       ggplot2::labs(title = paste0("Cumulative sum", scope_tag), subtitle = tr_note, x = labs0$x, y = "Cumsum")
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Histogram") {
  #     bins <- as.integer(input$stp_hist_bins %||% 30)
  #     g <- ggplot2::ggplot(df, ggplot2::aes(x = y_plot)) +
  #       ggplot2::geom_histogram(bins = bins, fill = fill_col, alpha = a) +
  #       ggplot2::labs(title = paste0("Histogram", scope_tag), subtitle = tr_note, x = labs0$y, y = "Count")
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Density") {
  #     bw_adj <- as.numeric(input$stp_bw_adj %||% 1)
  #     g <- ggplot2::ggplot(df, ggplot2::aes(x = y_plot)) +
  #       ggplot2::geom_density(adjust = bw_adj, fill = fill_col, alpha = 0.25) +
  #       ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a) +
  #       ggplot2::labs(title = paste0("Density", scope_tag), subtitle = tr_note, x = labs0$y, y = "Density")
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "QQ plot") {
  #     g <- ggplot2::ggplot(df, ggplot2::aes(sample = y_plot)) +
  #       ggplot2::stat_qq(color = pt_col, alpha = a) +
  #       ggplot2::stat_qq_line(color = line_col, linewidth = lw) +
  #       ggplot2::labs(title = paste0("QQ plot", scope_tag), subtitle = tr_note, x = "Theoretical", y = "Sample")
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Lag-1 scatter") {
  #     y <- df$y_plot
  #     validate(need(length(y) >= 3, "Not enough observations for lag scatter."))
  #     dfl <- data.frame(x = y[-length(y)], y = y[-1])
  #     g <- ggplot2::ggplot(dfl, ggplot2::aes(x = x, y = y)) +
  #       ggplot2::geom_point(color = pt_col, size = ps, alpha = a) +
  #       ggplot2::labs(title = paste0("Lag-1 scatter", scope_tag), subtitle = tr_note, x = "S(t-1)", y = "S(t)")
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt %in% c("Seasonal plot", "Seasonal subseries", "Polar seasonal", "Seasonal boxplot")) {
  #     validate(need(d$freq >= 2, "Seasonal plots require frequency >= 2."))
  #     
  #     x_ts <- d$ts
  #     
  #     if (pt == "Seasonal plot") {
  #       g <- forecast::ggseasonplot(x_ts, polar = FALSE) +
  #         ggplot2::labs(title = paste0("Seasonal plot", scope_tag), subtitle = tr_note)
  #       return(stp_apply_theme(stp_apply_palette(g)))
  #     }
  #     
  #     if (pt == "Polar seasonal") {
  #       g <- forecast::ggseasonplot(x_ts, polar = TRUE) +
  #         ggplot2::labs(title = paste0("Polar seasonal plot", scope_tag), subtitle = tr_note)
  #       return(stp_apply_theme(stp_apply_palette(g)))
  #     }
  #     
  #     if (pt == "Seasonal subseries") {
  #       g <- forecast::ggsubseriesplot(x_ts) +
  #         ggplot2::labs(title = paste0("Seasonal subseries", scope_tag), subtitle = tr_note)
  #       return(stp_apply_theme(g))
  #     }
  #     
  #     if (pt == "Seasonal boxplot") {
  #       df2 <- df
  #       df2$season <- as.factor(cycle(stats::ts(df2$y_plot, frequency = d$freq)))
  #       g <- ggplot2::ggplot(df2, ggplot2::aes(x = season, y = y_plot, fill = season)) +
  #         ggplot2::geom_boxplot(alpha = 0.6) +
  #         ggplot2::labs(title = paste0("Seasonal boxplot", scope_tag), subtitle = tr_note, x = "Season", y = labs0$y)
  #       return(stp_apply_theme(stp_apply_palette(g)))
  #     }
  #   }
  #   
  #   if (pt == "Classical decomposition (additive)" || pt == "Classical decomposition (multiplicative)") {
  #     validate(need(d$freq >= 2, "Decomposition requires frequency >= 2."))
  #     
  #     x_ts <- d$ts
  #     type <- if (pt == "Classical decomposition (multiplicative)") "multiplicative" else "additive"
  #     dc <- stats::decompose(x_ts, type = type)
  #     g <- forecast::autoplot(dc) +
  #       ggplot2::labs(title = paste0("Classical decomposition (", type, ")", scope_tag), subtitle = tr_note)
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "STL decomposition") {
  #     validate(need(d$freq >= 2, "STL requires frequency >= 2."))
  #     x_ts <- d$ts
  #     fit <- stats::stl(x_ts, s.window = "periodic", robust = TRUE)
  #     g <- forecast::autoplot(fit) +
  #       ggplot2::labs(title = paste0("STL decomposition", scope_tag), subtitle = tr_note)
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   if (pt == "Periodogram") {
  #     x_ts <- d$ts
  #     validate(need(length(x_ts) >= 8, "Not enough observations for periodogram."))
  #     taper <- as.numeric(input$stp_spec_taper %||% 0.1)
  #     
  #     sp <- stats::spec.pgram(x_ts, taper = taper, plot = FALSE)
  #     dfp <- data.frame(freq = sp$freq, spec = sp$spec)
  #     
  #     g <- ggplot2::ggplot(dfp, ggplot2::aes(x = freq, y = spec)) +
  #       ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a) +
  #       ggplot2::labs(title = paste0("Periodogram", scope_tag), subtitle = tr_note, x = "Frequency", y = "Spectral density")
  #     return(stp_apply_theme(g))
  #   }
  #   
  #   # fallback
  #   stp_apply_theme(base + ggplot2::geom_line(color = line_col, linewidth = lw, alpha = a))
  #   
  # }, width = function() as.numeric(input$stp_plot_width_px %||% 980),
  # height = function() as.numeric(input$stp_plot_height_px %||% 520))
  # 
  # # ---- ACF/PACF panels ----
  # output$stp_acf <- renderPlot({
  #   d <- stp_data()
  #   x_ts <- d$ts
  #   scope_tag <- if (identical(d$scope, "train")) " (Training only)" else " (Full series)"
  #   L <- min(60, length(x_ts) - 1)
  #   stp_apply_theme(forecast::ggAcf(x_ts, lag.max = L) +
  #                     ggplot2::labs(title = paste0("ACF", scope_tag)))
  # })
  # 
  # output$stp_pacf <- renderPlot({
  #   d <- stp_data()
  #   x_ts <- d$ts
  #   scope_tag <- if (identical(d$scope, "train")) " (Training only)" else " (Full series)"
  #   L <- min(60, length(x_ts) - 1)
  #   stp_apply_theme(forecast::ggPacf(x_ts, lag.max = L) +
  #                     ggplot2::labs(title = paste0("PACF", scope_tag)))
  # })
  # 
  # # ---- Lag grid (1..m) ----
  # output$stp_lag_grid <- renderPlot({
  #   d <- stp_data()
  #   df <- d$df
  #   scope_tag <- if (identical(d$scope, "train")) " (Training only)" else " (Full series)"
  #   
  #   m <- as.integer(input$stp_lag_m %||% 12)
  #   validate(need(m >= 1, "m must be >= 1."))
  #   
  #   y <- df$y_plot
  #   validate(need(length(y) >= (m + 2), "Not enough observations for requested lag grid."))
  #   
  #   out <- lapply(1:m, function(k) {
  #     data.frame(
  #       lag = paste0("Lag ", k),
  #       x = y[seq_len(length(y) - k)],
  #       y = y[(k + 1):length(y)]
  #     )
  #   })
  #   dfl <- do.call(rbind, out)
  #   
  #   pt_col <- input$stp_point_color %||% "#2C7FB8"
  #   ps <- as.numeric(input$stp_point_size %||% 2)
  #   a  <- as.numeric(input$stp_alpha %||% 1)
  #   
  #   g <- ggplot2::ggplot(dfl, ggplot2::aes(x = x, y = y)) +
  #     ggplot2::geom_point(color = pt_col, size = ps, alpha = a) +
  #     ggplot2::facet_wrap(~ lag, scales = "free", ncol = 3) +
  #     ggplot2::labs(title = paste0("Lag plots (1..", m, ")", scope_tag), x = "S(t-k)", y = "S(t)")
  #   
  #   stp_apply_theme(g)
  #   
  # }, width = function() as.numeric(input$stp_plot_width_px %||% 980),
  # height = function() round(as.numeric(input$stp_plot_height_px %||% 520) * 1.2))
  # 
  # 
  # 
  # 
  # 
 

  
  
  
  
  
  
  
  
  
  
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  # ============================================================
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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
    x <- x[is.finite(x)]
    
    k <- max(0L, to_int(input$adf_lags, 10L))
    
    # single source of truth for the deterministic type
    type_used <- tolower(as.character(input$adf_type %||% input$adfTypeSt2 %||% "drift"))
    type_used <- if (type_used %in% c("none", "drift", "trend")) type_used else "drift"
    
    alpha_used <- suppressWarnings(as.numeric(input$alphaSt2 %||% 0.05))
    if (!is.finite(alpha_used)) alpha_used <- 0.05
    
    out <- list(
      type_used  = type_used,
      alpha_used = alpha_used,
      k_used     = k,
      series     = x
    )
    
    # run tests (keep errors if any)
    out$adf <- tryCatch(
      tseries::adf.test(x, k = k, alternative = input$alternd2St %||% "stationary"),
      error = function(e) { out$adf_error <<- e$message; NULL }
    )
    
    out$kpss <- tryCatch({
      null_kpss <- if (identical(type_used, "trend")) "Trend" else "Level"
      tseries::kpss.test(x, null = null_kpss)
    }, error = function(e) { out$kpss_error <<- e$message; NULL })
    
    out$pp <- tryCatch(
      tseries::pp.test(x),
      error = function(e) { out$pp_error <<- e$message; NULL }
    )
    
    out$ur <- tryCatch(
      urca::ur.df(x, type = type_used, lags = k),
      error = function(e) { out$ur_error <<- e$message; NULL }
    )
    
    # NEW: Auto-lag ADF (AIC)
    out$ur_auto <- tryCatch(
      urca::ur.df(x, type = type_used, selectlags = "AIC"),
      error = function(e) { out$ur_auto_error <<- e$message; NULL }
    )
    
    # Infer chosen k from regression terms: number of z.diff* terms ~= k
    out$k_auto <- tryCatch({
      if (is.null(out$ur_auto)) return(NA_integer_)
      sm <- summary(out$ur_auto)
      rn <- rownames(sm@testreg$coefficients)
      as.integer(sum(grepl("^z\\.diff", rn)))
    }, error = function(e) NA_integer_)
    
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
  
  
  
 
  
  
  
  
  
  
  
  output$stationarity_interpretation <- renderPrint({
    req(stationarity())
    st <- stationarity()
    
    # ---- helpers ----
    `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
    to_num <- function(x, d = NA_real_) {
      y <- suppressWarnings(as.numeric(x))
      ifelse(is.finite(y), y, d)
    }
    to_int <- function(x, d = NA_integer_) {
      y <- suppressWarnings(as.integer(x))
      ifelse(is.finite(y), y, d)
    }
    
    fmt_p <- if (exists("fmt_p", inherits = TRUE)) get("fmt_p") else function(p) {
      if (!is.finite(p)) return("NA")
      format.pval(p, digits = 4, eps = .Machine$double.eps)
    }
    fmt_num <- if (exists("fmt_num", inherits = TRUE)) get("fmt_num") else function(x, d = 4) {
      if (!is.finite(x)) return("NA")
      format(round(x, d), nsmall = d)
    }
    
    # ---- series used ----
    n_used <- if (!is.null(st$series)) length(st$series) else NA_integer_
    
    # Prefer stationarity() metadata when available
    adf_type_ui <- st$type_used %||% (input$adfTypeSt2 %||% "drift")
    alpha <- to_num(st$alpha_used, to_num(input$alphaSt2, 0.05))
    
    # fixed-lag (user-provided) metadata
    k_lags <- as.integer(to_num(st$k_used, to_num(input$LagOrderADFd2St, 10)))
    
    # auto-lag metadata (AIC among 0..Kmax)
    k_auto <- as.integer(to_num(st$k_auto, NA_real_))
    max_lag_auto <- as.integer(to_num(st$max_lag_auto, NA_real_))
    
    # KPSS null pairing (informational)
    kpss_null_ui <- if (adf_type_ui == "trend") "tau" else "mu"
    
    # Map for ur.df critical values
    tau_row <- switch(adf_type_ui, "none" = "tau1", "drift" = "tau2", "trend" = "tau3", "tau3")
    alpha_col <- switch(
      as.character(alpha),
      "0.01" = "1pct",
      "0.05" = "5pct",
      "0.1"  = "10pct",
      "0.10" = "10pct",
      "5pct"
    )
    
    # =========================
    # NEW: Python-like ADF (no k specified) + report chosen lag
    # tseries::adf.test default lag rule:
    #   k_default = trunc((n - 1)^(1/3))
    # =========================
    k_default_adf <- NA_integer_
    if (is.finite(n_used) && n_used >= 2) {
      k_default_adf <- as.integer(trunc((n_used - 1)^(1/3)))
      if (!is.finite(k_default_adf) || k_default_adf < 0L) k_default_adf <- NA_integer_
    }
    
    adf_auto_obj  <- NULL
    adf_auto_p    <- NA_real_
    adf_auto_stat <- NA_real_
    adf_auto_err  <- NULL
    
    # Run WITHOUT specifying k (Python-like)
    if (!is.null(st$series) && length(st$series) >= 10) {
      adf_auto_obj <- tryCatch(
        tseries::adf.test(st$series, alternative = input$alternd2St %||% "stationary"),
        error = function(e) { adf_auto_err <<- e$message; NULL }
      )
      if (!is.null(adf_auto_obj)) {
        adf_auto_p    <- to_num(adf_auto_obj$p.value)
        adf_auto_stat <- to_num(adf_auto_obj$statistic)
      }
    }
    
    # ---- fixed-lag ur.df tau & crit ----
    tau_obs <- NA_real_
    tau_crit <- NA_real_
    if (!is.null(st$ur)) {
      tau_obs <- suppressWarnings(as.numeric(st$ur@teststat[tau_row]))
      if (!is.finite(tau_obs)) tau_obs <- suppressWarnings(as.numeric(st$ur@teststat[1]))
      
      cv <- tryCatch(st$ur@cval, error = function(e) NULL)
      if (!is.null(cv) && is.matrix(cv) &&
          tau_row %in% rownames(cv) && alpha_col %in% colnames(cv)) {
        tau_crit <- suppressWarnings(as.numeric(cv[tau_row, alpha_col]))
      }
    }
    
    # ---- auto-lag ur.df tau & crit ----
    tau_obs_auto <- NA_real_
    tau_crit_auto <- NA_real_
    if (!is.null(st$ur_auto)) {
      tau_obs_auto <- suppressWarnings(as.numeric(st$ur_auto@teststat[tau_row]))
      if (!is.finite(tau_obs_auto)) tau_obs_auto <- suppressWarnings(as.numeric(st$ur_auto@teststat[1]))
      
      cv2 <- tryCatch(st$ur_auto@cval, error = function(e) NULL)
      if (!is.null(cv2) && is.matrix(cv2) &&
          tau_row %in% rownames(cv2) && alpha_col %in% colnames(cv2)) {
        tau_crit_auto <- suppressWarnings(as.numeric(cv2[tau_row, alpha_col]))
      }
    }
    
    # ---- tseries numbers (your existing fixed-k ADF) ----
    adf_p    <- if (!is.null(st$adf))  to_num(st$adf$p.value)    else NA_real_
    adf_stat <- if (!is.null(st$adf))  to_num(st$adf$statistic)  else NA_real_
    
    kpss_p    <- if (!is.null(st$kpss)) to_num(st$kpss$p.value)   else NA_real_
    kpss_stat <- if (!is.null(st$kpss)) to_num(st$kpss$statistic) else NA_real_
    
    pp_p    <- if (!is.null(st$pp)) to_num(st$pp$p.value) else NA_real_
    pp_stat <- if (!is.null(st$pp)) to_num(st$pp$statistic) else NA_real_
    
    # ---- lags used by tests (best effort) ----
    adf_lags_used <- k_lags
    kpss_lags_used <- if (!is.null(st$kpss) && !is.null(st$kpss$parameter)) as.integer(to_num(st$kpss$parameter)) else NA_integer_
    pp_lags_used   <- if (!is.null(st$pp)   && !is.null(st$pp$parameter))   as.integer(to_num(st$pp$parameter))   else NA_integer_
    
    # ---- decisions ----
    adf_auto_dec_ts <- if (is.finite(adf_auto_p)) (adf_auto_p < alpha) else NA
    
    adf_dec_ts      <- if (is.finite(adf_p)) (adf_p < alpha) else NA
    adf_dec_ur      <- if (is.finite(tau_obs) && is.finite(tau_crit)) (tau_obs < tau_crit) else NA
    adf_dec_ur_auto <- if (is.finite(tau_obs_auto) && is.finite(tau_crit_auto)) (tau_obs_auto < tau_crit_auto) else NA
    
    kpss_reject <- if (is.finite(kpss_p)) (kpss_p < alpha) else NA
    pp_reject   <- if (is.finite(pp_p)) (pp_p < alpha) else NA
    
    # ---- HEADER ----
    cat("==========================================================================\n")
    cat("                 ACADEMIC REPORT: STATIONARITY ANALYSIS                   \n")
    cat("==========================================================================\n")
    cat(" Number of Observations Used: ", ifelse(is.finite(n_used), n_used, "NA"), "\n", sep = "")
    cat(" Significance Level (α): ", format(alpha, nsmall = 2), "\n", sep = "")
    cat("--------------------------------------------------------------------------\n")
    
    # ---- TEACHING NOTES: overview ----
    cat(" STATIONARITY ANALYSIS (overview):\n")
    cat(" • Stationarity in SARIMA usually means: mean/variance roughly stable over time and autocorrelation decays.\n")
    cat(" • ADF / PP: H0 = unit root (non-stationary). KPSS: H0 = stationary. They complement each other.\n")
    cat(" • Differencing guidance:\n")
    cat("    - Use d for non-seasonal trend; use D for seasonal unit roots.\n")
    cat("    - Avoid over-differencing: it can induce negative autocorrelation and inflate forecast variance.\n")
    cat(" • Deterministic terms (ur.df type):\n")
    cat("    - none: no intercept/trend; drift: intercept; trend: intercept + deterministic trend.\n")
    cat(" • Caution: structural breaks/outliers can make unit-root tests misleading. Always check plots.\n")
    cat("--------------------------------------------------------------------------\n")
    
    # ---- 1A) Python-like ADF (tseries::adf.test with automatic k) ----
    cat(" 1A. ADF — tseries::adf.test (AUTO Lag : k)\n")
    cat("--------------------------------------------------------------------------\n")
    
    if (is.null(adf_auto_obj)) {
      cat("    Not available.\n")
      if (!is.null(adf_auto_err)) cat("    Reason: ", adf_auto_err, "\n", sep = "")
    } else {
      cat(" DECISION RULE:\n")
      cat(" • H0: the series has a unit root (non-stationary)\n")
      cat(" • Ha: the series is stationary\n")
      cat(" • Reject H0 if P-value < α.\n\n")
      
      cat(" SPECIFICATION:\n")
      cat("    - Lag selection       : automatic (no k provided)\n")
      cat("    - Default lag formula : k = trunc((n - 1)^(1/3))\n")
      cat("    - n used              : ", ifelse(is.finite(n_used), n_used, "NA"), "\n", sep = "")
      cat("    - #Lags chosen (k)    : ", ifelse(is.finite(k_default_adf), k_default_adf, "NA"), "\n", sep = "")
      cat("    - Critical Value      : not printed (p-value based)\n\n")
      
      cat(" STATISTICS:\n")
      cat("    - Test Statistic (DF) : ", fmt_num(adf_auto_stat, 4), "\n", sep = "")
      cat("    - P-Value             : ", fmt_p(adf_auto_p), "\n\n", sep = "")
      
      cat(" DECISION:\n")
      if (isTRUE(adf_auto_dec_ts)) {
        cat("  P-value < α ⇒ REJECT H0: Evidence suggests STATIONARITY.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  This result mirrors common Python workflows: the ADF lag length is chosen automatically based on sample size.\n")
        cat("  Rejecting H0 suggests the series does not behave like a unit-root process under this default lag choice.\n")
        cat("  In SARIMA practice, this reduces the need for non-seasonal differencing (d), but you should still check seasonal\n")
        cat("  patterns (possibly D) and confirm with KPSS/PP plus ACF/PACF.\n")
      } else if (identical(adf_auto_dec_ts, FALSE)) {
        cat("  P-value ≥ α ⇒ FAIL TO REJECT H0: Evidence suggests NON-STATIONARITY.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  With automatic lag choice, failing to reject H0 indicates a unit root remains plausible, so differencing is typically\n")
        cat("  justified (often d = 1). After differencing, re-run tests: if results flip to stationarity and residual diagnostics look\n")
        cat("  white-noise-like, you have a good basis for ARIMA modeling.\n")
      } else {
        cat("  ADF p-value missing; decision is inconclusive.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  When p-values are unavailable/non-finite, rely more on the ur.df tau/critical approach, KPSS/PP, and especially the time\n")
        cat("  plot + ACF/PACF. Also verify the series isn’t constant and has enough variation after cleaning.\n")
      }
    }
    cat("--------------------------------------------------------------------------\n")
    
    # ---- 1) ADF (tseries::adf.test) fixed k (your existing one) ----
    cat(" 1B. ADF — tseries::adf.test (FIXED k)\n")
    cat("--------------------------------------------------------------------------\n")
    
    if (is.null(st$adf)) {
      cat("    Not available.\n")
      if (!is.null(st$adf_error)) cat("    Reason: ", st$adf_error, "\n", sep = "")
    } else {
      cat(" DECISION RULE:\n")
      cat(" • H0: the series has a unit root (non-stationary)\n")
      cat(" • Ha: the series is stationary\n")
      cat(" • Reject H0 if P-value < α.\n\n")
      
      cat(" SPECIFICATION:\n")
      cat("    - #Lags Used (k)      : ", ifelse(is.finite(adf_lags_used), adf_lags_used, "NA"), "\n", sep = "")
      cat("    - Critical Value      : not printed (p-value based)\n\n")
      
      cat(" STATISTICS:\n")
      cat("    - Test Statistic (DF) : ", fmt_num(adf_stat, 4), "\n", sep = "")
      cat("    - P-Value             : ", fmt_p(adf_p), "\n\n", sep = "")
      
      cat(" DECISION:\n")
      if (isTRUE(adf_dec_ts)) {
        cat("  P-value < α ⇒ REJECT H0: Evidence suggests STATIONARITY.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  Under the chosen fixed lag k, the unit-root null is rejected. This supports treating the series as stationary after\n")
        cat("  accounting for autocorrelation up to k lags. If this disagrees with KPSS, treat evidence as mixed and validate with\n")
        cat("  residual diagnostics and forecast performance.\n")
      } else if (identical(adf_dec_ts, FALSE)) {
        cat("  P-value ≥ α ⇒ FAIL TO REJECT H0: Evidence suggests NON-STATIONARITY.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  With the selected fixed lag, a unit root cannot be ruled out. Consider differencing (d) and then re-check ACF/PACF.\n")
        cat("  If the ACF still shows seasonal spikes at multiples of s, consider seasonal differencing (D).\n")
      } else {
        cat("  ADF p-value missing; decision is inconclusive.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  When ADF is inconclusive, rely on ur.df tau/critical, PP, and KPSS along with plots. In teaching settings, it helps to\n")
        cat("  compare fixed-k vs auto-k to show sensitivity to lag choice.\n")
      }
    }
    cat("--------------------------------------------------------------------------\n")
    
    # ---- 2) ADF (urca::ur.df) fixed lag ----
    cat(" 2A. ADF — urca::ur.df (FIXED lags; tau vs. critical)\n")
    cat("--------------------------------------------------------------------------\n")
    
    if (is.null(st$ur)) {
      cat("    Not available.\n")
      if (!is.null(st$ur_error)) cat("    Reason: ", st$ur_error, "\n", sep = "")
    } else {
      cat(" DECISION RULE:\n")
      cat(" • H0: the series has a unit root (non-stationary)\n")
      cat(" • Ha: the series is stationary\n")
      cat(" • Reject H0 if Tau-Observed < Tau-Critical.\n\n")
      
      cat(" SPECIFICATION:\n")
      cat("    - Model type          : ", adf_type_ui, "  (none/drift/trend)\n", sep = "")
      cat("    - #Lags Used (k)      : ", ifelse(is.finite(k_lags), k_lags, "NA"), "\n", sep = "")
      
      cat(" CRITICAL VALUE:\n")
      cat("    - Tau Critical        : ", fmt_num(tau_crit, 4), " (", alpha_col, ")\n\n", sep = "")
      
      cat(" STATISTICS:\n")
      cat("    - Tau Observed        : ", fmt_num(tau_obs, 4), "\n\n", sep = "")
      
      cat(" DECISION:\n")
      if (isTRUE(adf_dec_ur)) {
        cat("  Tau-Obs < Tau-Crit ⇒ REJECT H0: Stationarity supported.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  The ur.df version is useful for teaching because it explicitly shows the deterministic component and uses a tau statistic\n")
        cat("  compared to critical values. Rejection strengthens the case for d = 0, but always check if seasonality implies D > 0.\n")
      } else if (identical(adf_dec_ur, FALSE)) {
        cat("  Tau-Obs ≥ Tau-Crit ⇒ FAIL TO REJECT H0: Possible unit root.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  Failing to reject here suggests non-stationarity under the chosen deterministic component and lag order. In practice, try\n")
        cat("  differencing and re-run tests; in teaching, compare type=drift vs type=trend to show how assumptions change conclusions.\n")
      } else {
        cat("  Tau or critical value missing; cannot form a decision.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  If critical values are missing (package/version issues), use PP + KPSS and emphasize plot-based diagnostics.\n")
      }
    }
    cat("--------------------------------------------------------------------------\n")
    
    # ---- 2A) ADF (urca::ur.df) automatic lag selection ----
    cat(" 2B. ADF — urca::ur.df (AUTOMATIC lags via AIC; tau vs. critical)\n")
    cat("--------------------------------------------------------------------------\n")
    
    if (is.null(st$ur_auto)) {
      cat("    Not available.\n")
      if (!is.null(st$ur_auto_error)) cat("    Reason: ", st$ur_auto_error, "\n", sep = "")
    } else {
      cat(" DECISION RULE:\n")
      cat(" • H0: the series has a unit root (non-stationary)\n")
      cat(" • Ha: the series is stationary\n")
      cat(" • Lags are selected automatically using AIC within 0..Kmax.\n")
      cat(" • Reject H0 if Tau-Observed < Tau-Critical.\n\n")
      
      cat(" SPECIFICATION:\n")
      cat("    - Model type          : ", adf_type_ui, "\n", sep = "")
      cat("    - Kmax searched       : ", ifelse(is.finite(max_lag_auto), max_lag_auto, "NA"), "\n", sep = "")
      cat("    - #Lags chosen (auto) : ", ifelse(is.finite(k_auto), k_auto, "NA"), "\n", sep = "")
      
      cat(" CRITICAL VALUE:\n")
      cat("    - Tau Critical        : ", fmt_num(tau_crit_auto, 4), " (", alpha_col, ")\n\n", sep = "")
      
      cat(" STATISTICS:\n")
      cat("    - Tau Observed        : ", fmt_num(tau_obs_auto, 4), "\n\n", sep = "")
      
      cat(" DECISION:\n")
      if (isTRUE(adf_dec_ur_auto)) {
        cat("  Tau-Obs < Tau-Crit ⇒ REJECT H0: Stationarity supported.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  This is the most teachable/defensible ADF workflow: you set a maximum lag (Kmax) and let AIC pick the lag that best\n")
        cat("  explains autocorrelation without overfitting. Agreement with PP (reject unit root) and KPSS (do not reject stationarity)\n")
        cat("  is strong evidence for d = 0, subject to seasonal structure.\n")
      } else if (identical(adf_dec_ur_auto, FALSE)) {
        cat("  Tau-Obs ≥ Tau-Crit ⇒ FAIL TO REJECT H0: Possible unit root.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  Even with AIC-selected lag, the unit-root null remains plausible. This supports differencing. In teaching, highlight that\n")
        cat("  lag choice matters because too few lags leave autocorrelation in errors; too many reduce power.\n")
      } else {
        cat("  Tau or critical value missing; cannot form a decision.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  If tau/critical values can’t be extracted, treat as supportive context and rely on the other tests + plots.\n")
      }
    }
    cat("--------------------------------------------------------------------------\n")
    
    # ---- 3) KPSS ----
    cat(" 3. KPSS — tseries::kpss.test\n")
    cat("--------------------------------------------------------------------------\n")
    
    if (is.null(st$kpss)) {
      cat("    Not available.\n")
      if (!is.null(st$kpss_error)) cat("    Reason: ", st$kpss_error, "\n", sep = "")
    } else {
      cat(" DECISION RULE:\n")
      cat(" • H0: the series is stationary around a ", ifelse(kpss_null_ui == "tau", "trend", "level"), "\n", sep = "")
      cat(" • Ha: the series is non-stationary\n")
      cat(" • Reject H0 if P-value < α.\n\n")
      
      cat(" SPECIFICATION:\n")
      cat("    - #Lags Used          : ", ifelse(is.finite(kpss_lags_used), kpss_lags_used, "NA"), "\n", sep = "")
      cat("    - Critical Value      : not printed (p-value based)\n\n")
      
      cat(" STATISTICS:\n")
      cat("    - Test Statistic (η)  : ", fmt_num(kpss_stat, 4), "\n", sep = "")
      cat("    - P-Value             : ", fmt_p(kpss_p), "\n\n", sep = "")
      
      cat(" DECISION:\n")
      if (isTRUE(kpss_reject)) {
        cat("  P-value < α ⇒ REJECT H0: Evidence of NON-STATIONARITY.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  KPSS is the mirror-image of ADF/PP: it rejects stationarity. If ADF/PP fail to reject unit root AND KPSS rejects\n")
        cat("  stationarity, that’s convergent evidence for differencing. If KPSS conflicts with ADF, emphasize mixed evidence and\n")
        cat("  validate with residual diagnostics and forecast accuracy.\n")
      } else if (identical(kpss_reject, FALSE)) {
        cat("  P-value ≥ α ⇒ FAIL TO REJECT H0: Stationarity is plausible.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  Non-rejection supports stationarity around level/trend. In teaching, KPSS is great to show why ‘fail to reject’ is not\n")
        cat("  ‘prove stationarity’, but combined with ADF/PP it builds a coherent story.\n")
      } else {
        cat("  KPSS p-value missing; decision is inconclusive.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  If KPSS is unavailable/inconclusive, rely more on ADF/PP + plots and consider using multiple transformations.\n")
      }
    }
    cat("--------------------------------------------------------------------------\n")
    
    # ---- 4) Phillips–Perron ----
    cat(" 4. PHILLIPS–PERRON — tseries::pp.test\n")
    cat("--------------------------------------------------------------------------\n")
    
    if (is.null(st$pp)) {
      cat("    Not available.\n")
      if (!is.null(st$pp_error)) cat("    Reason: ", st$pp_error, "\n", sep = "")
    } else {
      cat(" DECISION RULE:\n")
      cat(" • H0: the series has a unit root (non-stationary)\n")
      cat(" • Ha: the series is stationary\n")
      cat(" • Reject H0 if P-value < α.\n\n")
      
      cat(" SPECIFICATION:\n")
      cat("    - #Lags Used          : ", ifelse(is.finite(pp_lags_used), pp_lags_used, "NA"), "\n", sep = "")
      cat("    - Critical Value      : not printed (p-value based)\n\n")
      
      cat(" STATISTICS:\n")
      cat("    - Test Statistic      : ", fmt_num(pp_stat, 4), "\n", sep = "")
      cat("    - P-Value             : ", fmt_p(pp_p), "\n\n", sep = "")
      
      cat(" DECISION:\n")
      if (isTRUE(pp_reject)) {
        cat("  P-value < α ⇒ REJECT H0: Stationarity suggested.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  PP is a robust alternative to ADF that corrects for serial correlation/heteroskedasticity nonparametrically.\n")
        cat("  In teaching, it’s useful to show that agreement between PP and ADF strengthens inference, while disagreement suggests\n")
        cat("  sensitivity to assumptions and motivates a conservative workflow.\n")
      } else if (identical(pp_reject, FALSE)) {
        cat("  P-value ≥ α ⇒ FAIL TO REJECT H0: Possible unit root.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  PP does not reject a unit root, supporting differencing. Re-test after differencing and check residual whiteness.\n")
      } else {
        cat("  PP p-value missing; decision is inconclusive.\n\n")
        cat(" CONCLUSION (teaching paragraph):\n")
        cat("  If PP is inconclusive, rely on the other tests + plots, and consider transformations (log/Box-Cox) before differencing.\n")
      }
    }
    cat("--------------------------------------------------------------------------\n")
    
    # ---- Final synthesis ----
    cat(" RECOMMENDATION (practical SARIMA workflow):\n")
    cat("--------------------------------------------------------------------------\n")
    
    cat(" • Step 1: Inspect time plot + ACF/PACF (look for slow decay and seasonal spikes).\n")
    cat(" • Step 2: Use evidence from ADF/PP (unit root) + KPSS (stationarity null) to decide d and D.\n")
    cat(" • Step 3: Apply minimal differencing; re-check tests and plots.\n")
    cat(" • Step 4: Fit SARIMA; verify residuals ~ white noise (Ljung-Box, residual ACF).\n")
    cat(" • Step 5: Prefer the model that passes diagnostics and forecasts well (not just lowest AIC).\n\n")
    
    # teaching: show how to interpret convergence
    if ((isTRUE(adf_auto_dec_ts) || isTRUE(adf_dec_ts) || isTRUE(adf_dec_ur) || isTRUE(adf_dec_ur_auto) || isTRUE(pp_reject)) && !isTRUE(kpss_reject)) {
      cat(" • Convergent evidence: STATIONARITY likely (consider d = 0; then evaluate seasonal differencing D).\n")
    } else if ((identical(adf_auto_dec_ts, FALSE) || identical(adf_dec_ts, FALSE) || identical(adf_dec_ur, FALSE) || identical(adf_dec_ur_auto, FALSE) || identical(pp_reject, FALSE)) && isTRUE(kpss_reject)) {
      cat(" • Convergent evidence: NON-STATIONARITY likely (difference: d and/or D; then re-test).\n")
    } else {
      cat(" • Mixed evidence: be conservative — try minimal differencing and validate using residual diagnostics + forecast accuracy.\n")
    }
    
    cat("\n==========================================================================\n\n")
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

  
  
  
  output$diff_plot <- renderPlot(
    {
      # =======================
      # Layout knobs (EDIT ME)
      # =======================
      LAG_MAX <- 40
      
      # When SPLIT (top + middle + (acf|pacf))
      # ↓↓↓ Reduce TOP, increase BOTTOM so ACF/PACF are readable
      H_TOP <- 0.7   # top: original series with split line
      H_MID <- 1.0   # middle: train-differenced series
      H_BOT <- 1.6   # bottom: ACF/PACF row
      
      # Width split between ACF and PACF
      W_LEFT  <- 1.0
      W_RIGHT <- 1.0
      
      # Panel styling knobs (teaching-friendly)
      TOP_TITLE_SIZE <- 11
      MID_TITLE_SIZE <- 12
      BOT_TITLE_SIZE <- 12
      
      TOP_MARGIN <- ggplot2::margin(t = 2, r = 6, b = 2, l = 6)
      MID_MARGIN <- ggplot2::margin(t = 6, r = 6, b = 6, l = 6)
      BOT_MARGIN <- ggplot2::margin(t = 6, r = 6, b = 6, l = 6)
      
      s <- ts_train_test()
      req(s$ts_train)
      
      d <- as.numeric(input$d_preview)
      D <- as.numeric(input$D_preview)
      
      has_test <- isTRUE(input$train_prop < 1) && length(s$ts_test) > 0
      
      # Full ORIGINAL series (no differencing)
      x_full_orig <- ts(
        c(as.numeric(s$ts_train), as.numeric(s$ts_test)),
        start = 1,
        frequency = frequency(s$ts_train)
      )
      
      # Differencing helper
      do_diff <- function(x) {
        x2 <- x
        if (d > 0) x2 <- diff(x2, differences = d)
        if (D > 0) x2 <- diff(x2, lag = frequency(x), differences = D)
        x2
      }
      
      # =========================
      # NO SPLIT: top reflects d,D
      # =========================
      if (!has_test) {
        x_full_diff <- do_diff(x_full_orig)
        
        p_top <- autoplot(x_full_diff) +
          ggtitle(paste0("Series preview (full), differenced (d=", d, ", D=", D, ")")) +
          xlab("Time") + ylab("Value") +
          theme(
            plot.title = element_text(size = MID_TITLE_SIZE),
            plot.margin = MID_MARGIN
          )
        
        p_acf  <- ggAcf(x_full_diff, lag.max = LAG_MAX) +
          ggtitle("ACF (full series, after differencing)") +
          theme(
            plot.title = element_text(size = BOT_TITLE_SIZE),
            plot.margin = BOT_MARGIN
          )
        
        p_pacf <- ggPacf(x_full_diff, lag.max = LAG_MAX) +
          ggtitle("PACF (full series, after differencing)") +
          theme(
            plot.title = element_text(size = BOT_TITLE_SIZE),
            plot.margin = BOT_MARGIN
          )
        
        (p_top / (p_acf | p_pacf)) +
          patchwork::plot_layout(heights = c(1.1, 1.0), widths = c(W_LEFT, W_RIGHT))
      } else {
        
        # =========================
        # SPLIT: top original + split line
        # =========================
        p_top <- autoplot(x_full_orig) +
          ggtitle("Original time series (no transformation) with Train/Test split") +
          xlab("Time") + ylab("Value") +
          theme(
            plot.title = element_text(size = TOP_TITLE_SIZE),
            plot.margin = TOP_MARGIN
          )
        
        split_time <- time(x_full_orig)[length(s$ts_train)]
        p_top <- p_top +
          geom_vline(xintercept = split_time, linetype = "dashed", color = "red") +
          annotate(
            "text",
            x = split_time,
            y = max(as.numeric(x_full_orig), na.rm = TRUE),
            label = "Train/Test split",
            vjust = -0.5, hjust = 0,
            color = "red",
            size = 3
          )
        
        # Middle: TRAIN-only differenced
        x_train_diff <- do_diff(s$ts_train)
        
        p_mid <- autoplot(x_train_diff) +
          ggtitle(paste0("Training portion only, differenced (d=", d, ", D=", D, ")")) +
          xlab("Time") + ylab("Differenced value") +
          theme(
            plot.title = element_text(size = MID_TITLE_SIZE),
            plot.margin = MID_MARGIN
          )
        
        # Bottom: ACF/PACF on TRAIN-diff
        p_acf <- ggAcf(x_train_diff, lag.max = LAG_MAX) +
          ggtitle("ACF (train only, after differencing)") +
          theme(
            plot.title = element_text(size = BOT_TITLE_SIZE),
            plot.margin = BOT_MARGIN
          )
        
        p_pacf <- ggPacf(x_train_diff, lag.max = LAG_MAX) +
          ggtitle("PACF (train only, after differencing)") +
          theme(
            plot.title = element_text(size = BOT_TITLE_SIZE),
            plot.margin = BOT_MARGIN
          )
        
        bottom_row <- (p_acf | p_pacf) +
          patchwork::plot_layout(widths = c(W_LEFT, W_RIGHT))
        
        # IMPORTANT: Force row heights so top is NOT huge
        (p_top / p_mid / bottom_row) +
          patchwork::plot_layout(heights = c(H_TOP, H_MID, H_BOT))
      }
    },
    
    # =======================
    # Output device sizing
    # =======================
    width = local({
      W_PX <- 990   # overall width (EDIT ME)
      W_PX
    }),
    
    height = local({
      H_NO_SPLIT_PX <- 750   # overall height without split (EDIT ME)
      H_SPLIT_PX    <- 950   # overall height with split (EDIT ME)
      function() {
        s <- ts_train_test()
        has_test <- isTRUE(input$train_prop < 1) &&
          !is.null(s$ts_test) &&
          length(s$ts_test) > 0
        if (has_test) H_SPLIT_PX else H_NO_SPLIT_PX
      }
    })
  )
  
  
  
 
  
  output$diff_plot <- renderPlot(
    {
      s <- ts_train_test()
      req(s$ts_train)
      
      d <- as.numeric(input$d_preview)
      D <- as.numeric(input$D_preview)
      
      has_test <- isTRUE(input$train_prop < 1) && length(s$ts_test) > 0
      
      # Build full ORIGINAL series (no differencing)
      x_full_orig <- ts(
        c(as.numeric(s$ts_train), as.numeric(s$ts_test)),
        start = 1,
        frequency = frequency(s$ts_train)
      )
      
      # Differencing helper (applies to any ts)
      do_diff <- function(x) {
        x2 <- x
        if (d > 0) x2 <- diff(x2, differences = d)
        if (D > 0) x2 <- diff(x2, lag = frequency(x), differences = D)
        x2
      }
      
      # ===== If NO SPLIT: upper plot should reflect parameters (d, D) =====
      if (!has_test) {
        x_full_diff <- do_diff(x_full_orig)
        
        p_top <- autoplot(x_full_diff) +
          ggtitle(paste0("Series preview (full), differenced (d=", d, ", D=", D, ")")) +
          xlab("Time") + ylab("Value")
        
        p_acf  <- ggAcf(x_full_diff, lag.max = 40) +
          ggtitle("ACF (full series, after differencing)")
        
        p_pacf <- ggPacf(x_full_diff, lag.max = 40) +
          ggtitle("PACF (full series, after differencing)")
        
        return(p_top / (p_acf | p_pacf))
      }
      
      # ===== If SPLIT: keep top ORIGINAL + split, then show TRAIN-diff + ACF/PACF =====
      p_top <- autoplot(x_full_orig) +
        ggtitle("Original time series (no transformation) with Train/Test split") +
        xlab("Time") + ylab("Value")
      
      split_time <- time(x_full_orig)[length(s$ts_train)]
      p_top <- p_top +
        geom_vline(xintercept = split_time, linetype = "dashed", color = "red") +
        annotate(
          "text",
          x = split_time,
          y = max(as.numeric(x_full_orig), na.rm = TRUE),
          label = "Train/Test split",
          vjust = -0.5, hjust = 0,
          color = "red",
          size = 3
        )
      
      # TRAIN-only after differencing
      x_train_diff <- do_diff(s$ts_train)
      
      p_train <- autoplot(x_train_diff) +
        ggtitle(paste0("Training portion only, differenced (d=", d, ", D=", D, ")")) +
        xlab("Time") + ylab("Differenced value")
      
      p_acf <- ggAcf(x_train_diff, lag.max = 40) +
        ggtitle("ACF (train only, after differencing)")
      
      p_pacf <- ggPacf(x_train_diff, lag.max = 40) +
        ggtitle("PACF (train only, after differencing)")
      
      p_top / (p_train / (p_acf | p_pacf))
    },
    width  = 990,
    height = function() {
      s <- ts_train_test()
      has_test <- isTRUE(input$train_prop < 1) && !is.null(s$ts_test) && length(s$ts_test) > 0
      if (has_test) 1050 else 750
    }
  )
  
 
  
  
  
  
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
  
  # Ljung–Box p-values across lags (Auto-ARIMA)
  output$auto_resid_lb_pvals <- renderPlot({
    req(auto_fit())
    
    res <- residuals(auto_fit())
    res <- as.numeric(res)[is.finite(res)]
    
    # Use the same controls you already expose for diagnostics
    L_input <- suppressWarnings(as.integer(input$diag_lag))
    L       <- if (is.finite(L_input) && L_input > 0) L_input else 12L
    fitdf   <- length(coef(auto_fit()))
    alpha   <- suppressWarnings(as.numeric(input$alphaSt2)); if (!is.finite(alpha)) alpha <- 0.05
    
    N <- length(res)
    validate(need(N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
    L <- max(1L, min(L, N - 1L))
    
    # Compute p-values at cumulative lags 1..L (omit k <= fitdf where df<=0)
    pvals <- rep(NA_real_, L)
    for (k in seq_len(L)) {
      if (k > fitdf) {
        bt <- tryCatch(stats::Box.test(res, lag = k, type = "Ljung-Box", fitdf = fitdf),
                       error = function(e) NULL)
        pvals[k] <- if (!is.null(bt)) as.numeric(bt$p.value) else NA_real_
      }
    }
    
    # Base R plotting to match your existing diagnostics
    plot(seq_len(L), pvals,
         type = "h", lwd = 2,
         xlab = "Lag (k)",
         ylab = "p-value  (Ljung–Box up to lag k)",
         main = "Ljung–Box p-values by lag (Auto-ARIMA)",
         ylim = c(0, 1))
    points(seq_len(L), pvals, pch = 16)
    abline(h = alpha, lty = 2)
    mtext(sprintf("alpha = %.3f,   fitdf = %d", alpha, fitdf), side = 3, line = 0.2, cex = 0.8)
    
    if (fitdf >= 1 && fitdf < L) {
      rect(xleft = 0.5, ybottom = -0.02, xright = fitdf + 0.5, ytop = 1.02,
           border = NA, col = grDevices::adjustcolor("gray", alpha.f = 0.15))
      text(x = (fitdf + 1) / 2, y = 0.95, labels = "df ≤ 0 (omitted)", cex = 0.8)
    }
  })
  
  
  
  
  
  
  
  
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
        
        tags$br(),
        tags$hr(),
        tags$hr(),
        
        tags$h4("General SARIMA formulation"),
        HTML(eq$eq_general),
        
        tags$br(),
        tags$hr(),
        tags$hr(),
        
        tags$h4("Expanded operator form"),
        HTML(eq$eq_expanded),
        
        tags$br(),
        tags$hr(),
        tags$hr(),
        
        tags$h4("Numerical model"),
        HTML(eq$eq_line3),
        
        tags$br(),
        tags$hr(),
        tags$hr(),,
        
        # HTML("\\[\\text{------------}\\]"),
        HTML(eq$eq_line4),
        
        tags$br(),
        tags$hr(),
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
  
  
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  
  

  # ---- Auto-ARIMA: academic conclusion (cached on Fit click) ----
  auto_conclusion_full_obj <- eventReactive(input$fit_auto, {
    req(auto_fit(), auto_fc(), auto_equations(), ts_train_test(), prepared())
    
    fit <- auto_fit()
    fc0 <- auto_fc()
    fc  <- fc0$fc
    eq  <- auto_equations()
    s   <- ts_train_test()
    p0  <- prepared()
    
    # ---- safe values
    n_train <- suppressWarnings(as.integer(s$train_n)); if (!is.finite(n_train)) n_train <- length(residuals(fit))
    n_test  <- suppressWarnings(as.integer(s$test_n));  if (!is.finite(n_test))  n_test  <- 0L
    N <- n_train + n_test
    
    L_in <- suppressWarnings(as.integer(input$diag_lag))
    L <- if (is.finite(L_in) && L_in > 0) L_in else 12L
    fitdf <- length(coef(fit))
    
    # ---- key stats
    AICc_val <- suppressWarnings(as.numeric(fit$aicc))
    BIC_val  <- suppressWarnings(as.numeric(fit$bic))
    AIC_val  <- suppressWarnings(as.numeric(fit$aic))
    
    # ---- residual tests (minimal, fast)
    res <- as.numeric(residuals(fit))
    res <- res[is.finite(res)]
    lb <- tryCatch(Box.test(res, lag = min(L, max(1L, floor(length(res) / 3))), type = "Ljung-Box", fitdf = fitdf),
                   error = function(e) NULL)
    
    # ---- accuracy (only if test exists)
    acc_line <- NULL
    if (n_test > 0) {
      acc <- tryCatch(accuracy_df(s$ts_test, fc$mean), error = function(e) NULL)
      if (!is.null(acc) && all(c("Metric", "Value") %in% names(acc))) {
        rmse <- acc$Value[acc$Metric == "RMSE"][1]
        mae  <- acc$Value[acc$Metric == "MAE"][1]
        mape <- acc$Value[acc$Metric == "MAPE"][1]
        acc_line <- tags$p(
          tags$b("Forecast accuracy (test set). "),
          HTML(paste0(
            "Over the holdout period (n = ", n_test, "), performance was RMSE = ",
            fmt_num(rmse, 3), ", MAE = ", fmt_num(mae, 3),
            if (is.finite(mape)) paste0(", MAPE = ", fmt_num(100 * mape, 2), "%") else "",
            "."
          ))
        )
      }
    }
    if (is.null(acc_line)) {
      acc_line <- tags$p(tags$b("Forecast accuracy. "),
                         "No holdout test set was detected; therefore, out-of-sample accuracy was not computed.")
    }
    
    # ---- horizon narrative (align with your Step 5 logic)
    horizon_txt <- if (n_test > 0) {
      paste0("Validation mode was used: the forecast horizon was forced to match the test length (h = ", n_test, ").")
    } else {
      paste0("Future mode was used: forecasts were produced beyond the training sample (h = ", fc0$h, ").")
    }
    
    # ---- “model string” consistent with your equation panel
    s_term <- if (is.finite(eq$s) && eq$s > 1) sprintf(" × (%d,%d,%d)[%d]", eq$P, eq$D, eq$Q, eq$s) else ""
    model_str <- sprintf("ARIMA(%d,%d,%d)%s", eq$p, eq$d, eq$q, s_term)
    
    # ---- build tag UI (fast + correct MathJax)
    tagList(
      tags$h3("Auto-ARIMA: Full academic conclusion"),
      
      tags$h4("1. Objective and modelling context"),
      tags$p(
        "An automated ARIMA procedure was used to establish a reproducible baseline for linear time-series dynamics. ",
        "Auto-ARIMA searches across candidate ARIMA/SARIMA structures and selects a parsimonious specification based on information criteria, ",
        "subject to the user-defined search settings (e.g., stepwise/approximation and seasonal allowance)."
      ),
      
      tags$h4("2. Sample design"),
      tags$p(
        HTML(paste0(
          "The analysis used <b>N = ", N, "</b> observations (training <b>n = ", n_train,
          "</b>", if (n_test > 0) paste0(", test <b>n = ", n_test, "</b>") else "",
          "). The seasonal period used in the workflow was <b>s = ", p0$freq, "</b>."
        ))
      ),
      tags$p(tags$b("Forecast design. "), horizon_txt),
      
      tags$h4("3. Selected specification and fit"),
      tags$p(HTML(paste0("The selected model was <b>", as.character(fit), "</b> (reported as <b>", model_str, "</b> in order notation)."))),
      tags$ul(
        tags$li(HTML(paste0("AIC = <b>", fmt_num(AIC_val, 2), "</b>"))),
        tags$li(HTML(paste0("AICc = <b>", fmt_num(AICc_val, 2), "</b>"))),
        tags$li(HTML(paste0("BIC = <b>", fmt_num(BIC_val, 2), "</b>")))
      ),
      tags$p(
        "In academic reporting, these criteria support relative comparison among candidate models; lower values indicate improved parsimony-adjusted fit."
      ),
      
      tags$h4("4. Model equations (reporting-ready)"),
      tags$p(
        "For documentation and replication, the fitted model is expressed in standard SARIMA operator notation, followed by an expanded form ",
        "and a numerical representation using the estimated coefficients."
      ),
      tags$div(
        style = "padding:10px;border:1px solid #e5e5e5;border-radius:6px;background:#fcfcfc;text-align:left;",
        
        tags$br(),
        tags$hr(),
        tags$hr(),
        
        tags$h5("General SARIMA formulation:"),
        HTML(eq$eq_general),
        
        tags$br(),
        tags$hr(),
        tags$hr(),
        
        tags$h5("Expanded operator form:"),
        HTML(eq$eq_expanded),
        
        tags$br(),
        tags$hr(),
        tags$hr(),
        
        tags$h5("Numerical model:"),
        HTML(eq$eq_line3),
        
        tags$br(),
        tags$hr(),
        tags$hr(),
        
        HTML(eq$eq_line4),
        
        tags$br(),
        tags$hr(),
        tags$hr(),
      ),
      
      tags$h4("5. Residual diagnostics (adequacy of linear dynamics)"),
      tags$p(
        "Adequacy was evaluated using graphical diagnostics (residual series, ACF, histogram, Q–Q plot) and formal residual tests. ",
        "A key criterion for ARIMA adequacy is that residuals approximate white noise (no remaining autocorrelation)."
      ),
      if (!is.null(lb)) {
        tags$p(HTML(paste0(
          "<b>Ljung–Box test:</b> Q(", lb$parameter, ") = ", fmt_num(lb$statistic, 3),
          ", p ", fmt_p(lb$p.value), "."
        )))
      } else {
        tags$p(tags$b("Ljung–Box test:"), " unavailable (insufficient residuals or test error).")
      },
      
      tags$h4("6. Forecasting and predictive performance"),
      acc_line,
      
      tags$h4("7. Overall conclusion and recommended next steps"),
      tags$p(
        "Overall, the Auto-ARIMA baseline provides a defensible benchmark for forecasting and for comparison against theory-guided manual SARIMA candidates. ",
        "If diagnostics suggest residual autocorrelation, a refined manual specification (guided by ACF/PACF after differencing) is recommended. ",
        "If volatility clustering is present (e.g., significant ARCH effects), ARIMA may be complemented by conditional variance modelling (e.g., GARCH)."
      )
    )
  })
  
  output$auto_conclusion_full <- renderUI({
    # Show a helpful message if user hasn't clicked Fit yet
    validate(need(input$fit_auto > 0, "Click “Fit Auto-ARIMA” to generate the full academic conclusion."))
    req(auto_conclusion_full_obj())
    
    # Force MathJax typesetting in the conclusion container
    session$onFlushed(function() {
      session$sendCustomMessage("mathjax-typeset", "auto_conclusion_box")
    }, once = TRUE)
    
    auto_conclusion_full_obj()
  })
  
  
  
  
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  

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
  # ---------------------------------------------   # --------------------------------------------- 
  # ---------------------------------------------   # --------------------------------------------- 
  # ---------------------------------------------   # --------------------------------------------- 
  # ---------------------------------------------   # --------------------------------------------- 
  # ---------------------------------------------   # --------------------------------------------- 
  # ---------------------------------------------   # --------------------------------------------- 
  # ---------------------------------------------   # --------------------------------------------- 
  # ---------------------------------------------   # --------------------------------------------- 
  # --------------------------------------------- 
  
  
  
  # # =========================
  # # For Manual Diag
  # # =========================
  # 
  # nice_par <- function() {
  #   par(mar = c(4, 4, 2.2, 1), mgp = c(2.2, 0.7, 0), las = 1)
  # }
  # 
  # # =========================
  # # Manual SARIMA diagnostics tab (NEW independent outputs)
  # # =========================
  # 
  # output$manual_resid_ts_diag <- renderPlot({
  #   req(manual_fit())
  #   nice_par()
  #   res <- residuals(manual_fit())
  #   plot(res, type = "l",
  #        main = "Residuals over time (Diagnostics tab)",
  #        xlab = "Time", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  # })
  # 
  # output$manual_resid_acf_diag <- renderPlot({
  #   req(manual_fit())
  #   nice_par()
  #   res <- residuals(manual_fit())
  #   plot(acf(res, plot = FALSE),
  #        main = "Residual ACF (Diagnostics tab)")
  # })
  # 
  # output$manual_resid_hist_diag <- renderPlot({
  #   req(manual_fit())
  #   nice_par()
  #   res <- residuals(manual_fit())
  #   hist(res, breaks = 30, col = "gray85", border = "white",
  #        main = "Residual histogram (Diagnostics tab)",
  #        xlab = "Residual")
  # })
  # 
  # output$manual_resid_qq_diag <- renderPlot({
  #   req(manual_fit())
  #   nice_par()
  #   res <- residuals(manual_fit())
  #   qqnorm(res, main = "Normal Q–Q (Diagnostics tab)")
  #   qqline(res, col = "red", lwd = 2)
  # })
  # 
  # # NEW: Residuals vs fitted (mean equation adequacy / nonlinearity)
  # output$manual_resid_fitted_diag <- renderPlot({
  #   req(manual_fit())
  #   nice_par()
  #   fit <- manual_fit()
  #   res <- residuals(fit)
  #   fitted_vals <- fitted(fit)
  #   
  #   plot(fitted_vals, res,
  #        pch = 16, cex = 0.7, col = grDevices::adjustcolor("steelblue", alpha.f = 0.7),
  #        main = "Residuals vs fitted",
  #        xlab = "Fitted values", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  #   # light smooth trend for visual cue (safe)
  #   ok <- is.finite(fitted_vals) & is.finite(res)
  #   if (sum(ok) > 20) {
  #     lines(stats::lowess(fitted_vals[ok], res[ok], f = 2/3), col = "tomato", lwd = 2)
  #   }
  # })
  # 
  # # NEW: Scale-location proxy (variance stability)
  # output$manual_resid_scale_diag <- renderPlot({
  #   req(manual_fit())
  #   nice_par()
  #   fit <- manual_fit()
  #   res <- residuals(fit)
  #   fitted_vals <- fitted(fit)
  #   
  #   y <- sqrt(abs(res))
  #   plot(fitted_vals, y,
  #        pch = 16, cex = 0.7, col = grDevices::adjustcolor("darkgreen", alpha.f = 0.65),
  #        main = "Scale-location (sqrt(|res|) vs fitted)",
  #        xlab = "Fitted values", ylab = expression(sqrt("|Residual|")))
  #   ok <- is.finite(fitted_vals) & is.finite(y)
  #   if (sum(ok) > 20) {
  #     lines(stats::lowess(fitted_vals[ok], y[ok], f = 2/3), col = "tomato", lwd = 2)
  #   }
  # })
  # 
  # # NEW: Ljung–Box p-values by lag (Diagnostics tab version)
  # output$manual_resid_lb_pvals_diag <- renderPlot({
  #   req(manual_fit())
  #   nice_par()
  #   
  #   res <- residuals(manual_fit())
  #   res <- as.numeric(res)[is.finite(res)]
  #   
  #   L_input <- suppressWarnings(as.integer(input$diag_lag))
  #   L       <- if (is.finite(L_input) && L_input > 0) L_input else 12L
  #   fitdf   <- length(coef(manual_fit()))
  #   alpha   <- suppressWarnings(as.numeric(input$alphaSt2)); if (!is.finite(alpha)) alpha <- 0.05
  #   
  #   N <- length(res)
  #   validate(need(N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   L <- max(1L, min(L, N - 1L))
  #   
  #   pvals <- rep(NA_real_, L)
  #   for (k in seq_len(L)) {
  #     if (k > fitdf) {
  #       bt <- tryCatch(stats::Box.test(res, lag = k, type = "Ljung-Box", fitdf = fitdf),
  #                      error = function(e) NULL)
  #       pvals[k] <- if (!is.null(bt)) as.numeric(bt$p.value) else NA_real_
  #     }
  #   }
  #   
  #   plot(seq_len(L), pvals,
  #        type = "h", lwd = 2,
  #        xlab = "Lag (k)",
  #        ylab = "p-value  (Ljung–Box up to lag k)",
  #        main = "Ljung–Box p-values by lag (Diagnostics tab)",
  #        ylim = c(0, 1))
  #   points(seq_len(L), pvals, pch = 16)
  #   abline(h = alpha, lty = 2, col = "gray40")
  #   mtext(sprintf("alpha = %.3f, fitdf = %d", alpha, fitdf), side = 3, line = 0.2, cex = 0.85)
  #   
  #   if (fitdf >= 1 && fitdf < L) {
  #     rect(xleft = 0.5, ybottom = -0.02, xright = fitdf + 0.5, ytop = 1.02,
  #          border = NA, col = grDevices::adjustcolor("gray", alpha.f = 0.15))
  #     text(x = (fitdf + 1) / 2, y = 0.95, labels = "df ≤ 0 (omitted)", cex = 0.8)
  #   }
  # })
  # 
  # # NEW: Commentary + conclusions (Diagnostics tab)
  # output$manual_diag_commentary_diag <- renderPrint({
  #   req(manual_fit())
  #   fit   <- manual_fit()
  #   res   <- as.numeric(residuals(fit))
  #   res   <- res[is.finite(res)]
  #   N     <- length(res)
  #   
  #   L_input <- suppressWarnings(as.integer(input$diag_lag))
  #   L       <- if (is.finite(L_input) && L_input > 0) L_input else 12L
  #   fitdf   <- length(coef(fit))
  #   alpha   <- suppressWarnings(as.numeric(input$alphaSt2)); if (!is.finite(alpha)) alpha <- 0.05
  #   
  #   if (N < 8) {
  #     cat("Diagnostics summary\n",
  #         "- Too few residuals to run a stable diagnostic battery.\n",
  #         "- Fit the model on more observations or reduce differencing/orders.\n", sep = "")
  #     return(invisible())
  #   }
  #   
  #   # Ljung–Box at chosen lag
  #   L2 <- max(1L, min(L, N - 1L))
  #   lb <- tryCatch(stats::Box.test(res, lag = L2, type = "Ljung-Box", fitdf = fitdf),
  #                  error = function(e) NULL)
  #   
  #   # Simple normality flag (optional tests already exist elsewhere in your app)
  #   sw <- if (N >= 3 && N <= 5000) tryCatch(stats::shapiro.test(res), error = function(e) NULL) else NULL
  #   
  #   # Heuristic variance stability from scale-location trend:
  #   fitted_vals <- as.numeric(fitted(fit))
  #   ok <- is.finite(fitted_vals) & is.finite(res)
  #   slope_flag <- NA
  #   if (sum(ok) > 30) {
  #     # correlate fitted with abs(res) as a crude variance trend signal
  #     slope_flag <- suppressWarnings(stats::cor(fitted_vals[ok], abs(res[ok]), use = "complete.obs"))
  #   }
  #   
  #   cat("Diagnostics summary (Manual SARIMA)\n")
  #   cat("------------------------------------------------------------\n")
  #   cat("* White-noise / autocorrelation:\n")
  #   if (!is.null(lb) && is.finite(lb$p.value)) {
  #     cat("  - Ljung–Box Q(", lb$parameter, ") = ", round(as.numeric(lb$statistic), 3),
  #         ", p = ", signif(as.numeric(lb$p.value), 3),
  #         if (lb$p.value > alpha) "  -> No strong evidence of residual autocorrelation.\n"
  #         else "  -> Evidence of remaining autocorrelation (model may be under-specified).\n",
  #         sep = "")
  #   } else {
  #     cat("  - Ljung–Box unavailable (test error or insufficient df).\n")
  #   }
  #   
  #   cat("* Distribution / normality:\n")
  #   if (!is.null(sw) && is.finite(sw$p.value)) {
  #     cat("  - Shapiro–Wilk p = ", signif(sw$p.value, 3),
  #         if (sw$p.value > alpha) "  -> residuals not strongly inconsistent with normality.\n"
  #         else "  -> residuals deviate from normality (heavy tails/skew possible).\n",
  #         sep = "")
  #   } else {
  #     cat("  - Shapiro–Wilk not computed (N out of range or error).\n")
  #   }
  #   
  #   cat("* Variance stability:\n")
  #   if (is.finite(slope_flag)) {
  #     cat("  - Corr(fitted, |res|) = ", round(slope_flag, 3),
  #         if (abs(slope_flag) < 0.15) "  -> variance looks roughly stable.\n"
  #         else "  -> possible heteroskedasticity (consider variance modeling / transformation).\n",
  #         sep = "")
  #   } else {
  #     cat("  - Not enough usable points to assess variance trend.\n")
  #   }
  #   
  #   cat("\nConclusion / next steps:\n")
  #   cat("------------------------------------------------------------\n")
  #   cat("- If ACF shows significant spikes and Ljung–Box rejects: revisit (d, D) and AR/MA orders; use ACF/PACF guidance.\n")
  #   cat("- If residuals are heavy-tailed or variance changes with level: consider transformation and/or adding a volatility model (e.g., GARCH) for conditional variance.\n")
  #   cat("- If residual-vs-fitted shows curvature: the mean equation may need additional structure (regressors, intervention dummies, or nonlinear terms).\n")
  # })
  # 
  # 
  # 
  # 
  # 
  # # Ljung–Box p-values across lags (Manual SARIMA)
  # output$manual_resid_lb_pvals <- renderPlot({
  #   req(manual_fit())
  #   res <- residuals(manual_fit())
  #   res <- as.numeric(res)[is.finite(res)]
  #   
  #   # Use the same controls you already expose elsewhere
  #   L_input <- suppressWarnings(as.integer(input$diag_lag))
  #   L       <- if (is.finite(L_input) && L_input > 0) L_input else 12L
  #   fitdf   <- length(coef(manual_fit()))
  #   alpha   <- suppressWarnings(as.numeric(input$alphaSt2)); if (!is.finite(alpha)) alpha <- 0.05
  #   
  #   N <- length(res)
  #   validate(need(N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   L <- max(1L, min(L, N - 1L))  # sane cap
  #   
  #   # Compute p-values for lags 1..L. For k <= fitdf, df would be nonpositive;
  #   # we mark those as NA (not meaningful after parameter adjustment).
  #   pvals <- rep(NA_real_, L)
  #   for (k in seq_len(L)) {
  #     if (k > fitdf) {
  #       bt <- tryCatch(stats::Box.test(res, lag = k, type = "Ljung-Box", fitdf = fitdf),
  #                      error = function(e) NULL)
  #       pvals[k] <- if (!is.null(bt)) as.numeric(bt$p.value) else NA_real_
  #     }
  #   }
  #   
  #   # Base R plot (consistent with your other diagnostics)
  #   plot(seq_len(L), pvals,
  #        type = "h", lwd = 2,
  #        xlab = "Lag (k)",
  #        ylab = "p-value  (Ljung–Box up to lag k)",
  #        main = "Ljung–Box p-values by lag (Manual SARIMA)",
  #        ylim = c(0, 1))
  #   points(seq_len(L), pvals, pch = 16)
  #   abline(h = alpha, lty = 2)           # significance threshold
  #   mtext(sprintf("alpha = %.3f,   fitdf = %d", alpha, fitdf), side = 3, line = 0.2, cex = 0.8)
  #   # optional: annotate NA region when k <= fitdf
  #   if (fitdf >= 1 && fitdf < L) {
  #     rect(xleft = 0.5, ybottom = -0.02, xright = fitdf + 0.5, ytop = 1.02,
  #          border = NA, col = grDevices::adjustcolor("gray", alpha.f = 0.15))
  #     text(x = (fitdf + 1) / 2, y = 0.95, labels = "df ≤ 0 (omitted)", cex = 0.8)
  #   }
  # })
  # 
  
  
  
  
  # --------------------------------------------- 
  # --------------------------------------------- 
  
  # # =========================
  # # For Manual Diag (FULL REPLACEMENT)
  # # =========================
  # 
  # # `%||%` <- function(a, b) if (!is.null(a)) a else b
  # 
  # nice_par <- function() {
  #   par(mar = c(4, 4, 2.2, 1), mgp = c(2.2, 0.7, 0), las = 1)
  # }
  # 
  # # ---- Shared helpers ----
  # get_diag_controls <- function(input) {
  #   L_input <- suppressWarnings(as.integer(input$diag_lag))
  #   L <- if (is.finite(L_input) && L_input > 0) L_input else 12L
  #   
  #   alpha <- suppressWarnings(as.numeric(input$alphaSt2))
  #   if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05
  #   
  #   list(L = L, alpha = alpha)
  # }
  # 
  # compute_lb_pvals <- function(res, L, fitdf) {
  #   res <- as.numeric(res)
  #   res <- res[is.finite(res)]
  #   
  #   N <- length(res)
  #   if (N < 8) return(list(pvals = numeric(0), L = 0L, N = N))
  #   
  #   L <- as.integer(L)
  #   L <- max(1L, min(L, N - 1L))
  #   
  #   pvals <- rep(NA_real_, L)
  #   for (k in seq_len(L)) {
  #     if (k > fitdf) {
  #       bt <- tryCatch(
  #         stats::Box.test(res, lag = k, type = "Ljung-Box", fitdf = fitdf),
  #         error = function(e) NULL
  #       )
  #       pvals[k] <- if (!is.null(bt)) as.numeric(bt$p.value) else NA_real_
  #     }
  #   }
  #   
  #   list(pvals = pvals, L = L, N = N)
  # }
  # 
  # plot_lb_pvals <- function(pvals, L, alpha, fitdf, main_title) {
  #   plot(seq_len(L), pvals,
  #        type = "h", lwd = 2,
  #        xlab = "Lag (k)",
  #        ylab = "p-value  (Ljung–Box up to lag k)",
  #        main = main_title,
  #        ylim = c(0, 1))
  #   points(seq_len(L), pvals, pch = 16)
  #   abline(h = alpha, lty = 2, col = "gray40")
  #   mtext(sprintf("alpha = %.3f, fitdf = %d", alpha, fitdf), side = 3, line = 0.2, cex = 0.85)
  #   
  #   # Shade region where k <= fitdf (df <= 0 after adjustment -> omitted)
  #   if (is.finite(fitdf) && fitdf >= 1 && fitdf < L) {
  #     rect(xleft = 0.5, ybottom = -0.02, xright = fitdf + 0.5, ytop = 1.02,
  #          border = NA, col = grDevices::adjustcolor("gray", alpha.f = 0.15))
  #     text(x = (fitdf + 1) / 2, y = 0.95, labels = "df ≤ 0 (omitted)", cex = 0.8)
  #   }
  # }
  # 
  # # =========================
  # # Manual SARIMA diagnostics tab (NEW independent outputs)
  # # =========================
  # 
  # output$manual_resid_ts_diag <- renderPlot({
  #   req(manual_fit())
  #   nice_par()
  #   res <- as.numeric(residuals(manual_fit()))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 3, "Not enough residuals to plot."))
  #   
  #   plot(res, type = "l",
  #        main = "Residuals over time",
  #        xlab = "Time", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  # })
  # 
  # output$manual_resid_acf_diag <- renderPlot({
  #   req(manual_fit())
  #   nice_par()
  #   res <- as.numeric(residuals(manual_fit()))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for ACF."))
  #   
  #   plot(acf(res, plot = FALSE),
  #        main = "Residual ACF")
  # })
  # 
  # output$manual_resid_hist_diag <- renderPlot({
  #   req(manual_fit())
  #   nice_par()
  #   res <- as.numeric(residuals(manual_fit()))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for histogram."))
  #   
  #   hist(res, breaks = 30, col = "gray85", border = "white",
  #        main = "Residual histogram",
  #        xlab = "Residual")
  # })
  # 
  # output$manual_resid_qq_diag <- renderPlot({
  #   req(manual_fit())
  #   nice_par()
  #   res <- as.numeric(residuals(manual_fit()))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for Q-Q plot."))
  #   
  #   qqnorm(res, main = "Normal Q–Q")
  #   qqline(res, col = "red", lwd = 2)
  # })
  # 
  # # Residuals vs fitted (mean equation adequacy / nonlinearity)
  # output$manual_resid_fitted_diag <- renderPlot({
  #   req(manual_fit())
  #   nice_par()
  #   fit <- manual_fit()
  #   
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   validate(need(sum(ok) >= 10, "Not enough valid points for residuals vs fitted."))
  #   
  #   plot(fv[ok], res[ok],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("steelblue", alpha.f = 0.7),
  #        main = "Residuals vs fitted",
  #        xlab = "Fitted values", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  #   
  #   if (sum(ok) > 20) {
  #     lines(stats::lowess(fv[ok], res[ok], f = 2/3), col = "tomato", lwd = 2)
  #   }
  # })
  # 
  # # Scale-location proxy (variance stability)
  # output$manual_resid_scale_diag <- renderPlot({
  #   req(manual_fit())
  #   nice_par()
  #   fit <- manual_fit()
  #   
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   
  #   y <- sqrt(abs(res))
  #   ok2 <- ok & is.finite(y)
  #   validate(need(sum(ok2) >= 10, "Not enough valid points for scale-location plot."))
  #   
  #   plot(fv[ok2], y[ok2],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("darkgreen", alpha.f = 0.65),
  #        main = "Scale-location (sqrt(|res|) vs fitted)",
  #        xlab = "Fitted values", ylab = expression(sqrt("|Residual|")))
  #   if (sum(ok2) > 20) {
  #     lines(stats::lowess(fv[ok2], y[ok2], f = 2/3), col = "tomato", lwd = 2)
  #   }
  # })
  # 
  # # Ljung–Box p-values by lag (Diagnostics tab version)
  # output$manual_resid_lb_pvals_diag <- renderPlot({
  #   req(manual_fit())
  #   nice_par()
  #   
  #   fit <- manual_fit()
  #   res <- residuals(fit)
  #   fitdf <- length(coef(fit))
  #   
  #   ctrl <- get_diag_controls(input)
  #   out <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
  #   
  #   validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
  #   
  #   plot_lb_pvals(out$pvals, out$L, ctrl$alpha, fitdf,
  #                 main_title = "Ljung–Box p-values by lag (Diagnostics)")
  # })
  # 
  # # Commentary + conclusions (Diagnostics tab)
  # output$manual_diag_commentary_diag <- renderPrint({
  #   req(manual_fit())
  #   fit <- manual_fit()
  #   
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   N   <- length(res)
  #   
  #   ctrl  <- get_diag_controls(input)
  #   fitdf <- length(coef(fit))
  #   
  #   if (N < 8) {
  #     cat(
  #       "Diagnostics summary (Manual SARIMA)\n",
  #       "------------------------------------------------------------\n",
  #       "- Too few residuals to run a stable diagnostic battery.\n",
  #       "- Fit the model on more observations or reduce differencing/orders.\n",
  #       sep = ""
  #     )
  #     return(invisible())
  #   }
  #   
  #   # Ljung–Box at chosen lag
  #   L2 <- max(1L, min(as.integer(ctrl$L), N - 1L))
  #   lb <- tryCatch(
  #     stats::Box.test(res, lag = L2, type = "Ljung-Box", fitdf = fitdf),
  #     error = function(e) NULL
  #   )
  #   
  #   # Shapiro–Wilk (guarded)
  #   sw <- if (N >= 3 && N <= 5000) tryCatch(stats::shapiro.test(res), error = function(e) NULL) else NULL
  #   
  #   # Heteroskedasticity proxy: Corr(fitted, |res|)
  #   fv <- as.numeric(fitted(fit))
  #   ok <- is.finite(fv) & is.finite(res)
  #   rho <- if (sum(ok) > 30) suppressWarnings(stats::cor(fv[ok], abs(res[ok]), use = "complete.obs")) else NA_real_
  #   
  #   cat("Diagnostics summary (Manual SARIMA)\n")
  #   cat("------------------------------------------------------------\n")
  #   
  #   cat("* White-noise / autocorrelation:\n")
  #   if (!is.null(lb) && is.finite(lb$p.value)) {
  #     cat(
  #       "  - Ljung–Box Q(", lb$parameter, ") = ", round(as.numeric(lb$statistic), 3),
  #       ", p = ", signif(as.numeric(lb$p.value), 3),
  #       if (lb$p.value > ctrl$alpha) "  -> No strong evidence of residual autocorrelation.\n"
  #       else "  -> Evidence of remaining autocorrelation (consider revising AR/MA orders or differencing).\n",
  #       sep = ""
  #     )
  #   } else {
  #     cat("  - Ljung–Box unavailable (test error or insufficient df).\n")
  #   }
  #   
  #   cat("* Distribution / normality:\n")
  #   if (!is.null(sw) && is.finite(sw$p.value)) {
  #     cat(
  #       "  - Shapiro–Wilk p = ", signif(sw$p.value, 3),
  #       if (sw$p.value > ctrl$alpha) "  -> residuals not strongly inconsistent with normality.\n"
  #       else "  -> residuals deviate from normality (heavy tails/skew possible).\n",
  #       sep = ""
  #     )
  #   } else {
  #     cat("  - Shapiro–Wilk not computed (N out of range or error).\n")
  #   }
  #   
  #   cat("* Variance stability:\n")
  #   if (is.finite(rho)) {
  #     cat(
  #       "  - Corr(fitted, |res|) = ", round(rho, 3),
  #       if (abs(rho) < 0.15) "  -> variance looks roughly stable.\n"
  #       else "  -> possible heteroskedasticity (consider transformation or volatility modeling).\n",
  #       sep = ""
  #     )
  #   } else {
  #     cat("  - Not enough usable points to assess variance trend.\n")
  #   }
  #   
  #   cat("\nConclusion / next steps:\n")
  #   cat("------------------------------------------------------------\n")
  #   cat("- If ACF has spikes and Ljung–Box rejects: revisit (p,q)(P,Q), and ensure (d,D,s) are appropriate.\n")
  #   cat("- If Q–Q shows heavy tails: prediction intervals may be optimistic; consider robust alternatives or bootstrapping.\n")
  #   cat("- If residual-vs-fitted shows structure: consider exogenous regressors, interventions, or model refinement.\n")
  # })
  # 
  # # =========================
  # # Ljung–Box p-values across lags (Manual SARIMA)  [REPORT VERSION]
  # # Keep this output name because the Academic conclusion (full) may reference it.
  # # =========================
  # output$manual_resid_lb_pvals <- renderPlot({
  #   req(manual_fit())
  #   nice_par()
  #   
  #   fit <- manual_fit()
  #   res <- residuals(fit)
  #   fitdf <- length(coef(fit))
  #   
  #   ctrl <- get_diag_controls(input)
  #   out <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
  #   
  #   validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
  #   
  #   plot_lb_pvals(out$pvals, out$L, ctrl$alpha, fitdf,
  #                 main_title = "Ljung–Box p-values by lag (Manual SARIMA)")
  # })
  # 
  # # --------------------------------------------- 
  # # --------------------------------------------- 
  # 
  
  
  
  
  
  
  
  
  
  
  
  # # =========================
  # # For Manual Diag (FULL REPLACEMENT + UI BUILDER)
  # # =========================
  # 
  # # `%||%` <- function(a, b) if (!is.null(a)) a else b
  # 
  # nice_par <- function() {
  #   par(mar = c(4, 4, 2.2, 1), mgp = c(2.2, 0.7, 0), las = 1)
  # }
  # 
  # # ---- Shared helpers ----
  # get_diag_controls <- function(input) {
  #   L_input <- suppressWarnings(as.integer(input$diag_lag))
  #   L <- if (is.finite(L_input) && L_input > 0) L_input else 12L
  #   
  #   alpha <- suppressWarnings(as.numeric(input$alphaSt2))
  #   if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05
  #   
  #   list(L = L, alpha = alpha)
  # }
  # 
  # compute_lb_pvals <- function(res, L, fitdf) {
  #   res <- as.numeric(res)
  #   res <- res[is.finite(res)]
  #   
  #   N <- length(res)
  #   if (N < 8) return(list(pvals = numeric(0), L = 0L, N = N))
  #   
  #   L <- as.integer(L)
  #   L <- max(1L, min(L, N - 1L))
  #   
  #   pvals <- rep(NA_real_, L)
  #   for (k in seq_len(L)) {
  #     if (k > fitdf) {
  #       bt <- tryCatch(
  #         stats::Box.test(res, lag = k, type = "Ljung-Box", fitdf = fitdf),
  #         error = function(e) NULL
  #       )
  #       pvals[k] <- if (!is.null(bt)) as.numeric(bt$p.value) else NA_real_
  #     }
  #   }
  #   
  #   list(pvals = pvals, L = L, N = N)
  # }
  # 
  # plot_lb_pvals <- function(pvals, L, alpha, fitdf, main_title) {
  #   plot(seq_len(L), pvals,
  #        type = "h", lwd = 2,
  #        xlab = "Lag (k)",
  #        ylab = "p-value  (Ljung–Box up to lag k)",
  #        main = main_title,
  #        ylim = c(0, 1))
  #   points(seq_len(L), pvals, pch = 16)
  #   abline(h = alpha, lty = 2, col = "gray40")
  #   mtext(sprintf("alpha = %.3f, fitdf = %d", alpha, fitdf),
  #         side = 3, line = 0.2, cex = 0.85)
  #   
  #   if (is.finite(fitdf) && fitdf >= 1 && fitdf < L) {
  #     rect(xleft = 0.5, ybottom = -0.02, xright = fitdf + 0.5, ytop = 1.02,
  #          border = NA, col = grDevices::adjustcolor("gray", alpha.f = 0.15))
  #     text(x = (fitdf + 1) / 2, y = 0.95, labels = "df ≤ 0 (omitted)", cex = 0.8)
  #   }
  # }
  # 
  # # =========================
  # # Manual SARIMA diagnostics tab (NEW independent outputs)
  # # =========================
  # 
  # output$manual_resid_ts_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 3, "Not enough residuals to plot."))
  #   
  #   plot(res, type = "l",
  #        main = "Residuals over time",
  #        xlab = "Time", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  # })
  # 
  # output$manual_resid_acf_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for ACF."))
  #   
  #   plot(acf(res, plot = FALSE),
  #        main = "Residual ACF")
  # })
  # 
  # output$manual_resid_hist_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for histogram."))
  #   
  #   hist(res, breaks = 30, col = "gray85", border = "white",
  #        main = "Residual histogram",
  #        xlab = "Residual")
  # })
  # 
  # output$manual_resid_qq_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for Q-Q plot."))
  #   
  #   qqnorm(res, main = "Normal Q–Q")
  #   qqline(res, col = "red", lwd = 2)
  # })
  # 
  # output$manual_resid_fitted_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   validate(need(sum(ok) >= 10, "Not enough valid points for residuals vs fitted."))
  #   
  #   plot(fv[ok], res[ok],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("steelblue", alpha.f = 0.7),
  #        main = "Residuals vs fitted",
  #        xlab = "Fitted values", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  #   
  #   if (sum(ok) > 20) {
  #     lines(stats::lowess(fv[ok], res[ok], f = 2/3), col = "tomato", lwd = 2)
  #   }
  # })
  # 
  # output$manual_resid_scale_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   
  #   y <- sqrt(abs(res))
  #   ok2 <- ok & is.finite(y)
  #   validate(need(sum(ok2) >= 10, "Not enough valid points for scale-location plot."))
  #   
  #   plot(fv[ok2], y[ok2],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("darkgreen", alpha.f = 0.65),
  #        main = "Scale-location (sqrt(|res|) vs fitted)",
  #        xlab = "Fitted values", ylab = expression(sqrt("|Residual|")))
  #   if (sum(ok2) > 20) {
  #     lines(stats::lowess(fv[ok2], y[ok2], f = 2/3), col = "tomato", lwd = 2)
  #   }
  # })
  # 
  # output$manual_resid_lb_pvals_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res   <- residuals(fit)
  #   fitdf <- length(coef(fit))
  #   
  #   ctrl <- get_diag_controls(input)
  #   out  <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
  #   
  #   validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
  #   
  #   plot_lb_pvals(out$pvals, out$L, ctrl$alpha, fitdf,
  #                 main_title = "Ljung–Box p-values by lag (Diagnostics)")
  # })
  # 
  # output$manual_diag_commentary_diag <- renderPrint({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   N   <- length(res)
  #   
  #   ctrl  <- get_diag_controls(input)
  #   fitdf <- length(coef(fit))
  #   
  #   if (N < 8) {
  #     cat(
  #       "Diagnostics summary (Manual SARIMA)\n",
  #       "------------------------------------------------------------\n",
  #       "- Too few residuals to run a stable diagnostic battery.\n",
  #       "- Fit the model on more observations or reduce differencing/orders.\n",
  #       sep = ""
  #     )
  #     return(invisible())
  #   }
  #   
  #   L2 <- max(1L, min(as.integer(ctrl$L), N - 1L))
  #   lb <- tryCatch(
  #     stats::Box.test(res, lag = L2, type = "Ljung-Box", fitdf = fitdf),
  #     error = function(e) NULL
  #   )
  #   
  #   sw <- if (N >= 3 && N <= 5000) tryCatch(stats::shapiro.test(res), error = function(e) NULL) else NULL
  #   
  #   fv <- as.numeric(fitted(fit))
  #   ok <- is.finite(fv) & is.finite(res)
  #   rho <- if (sum(ok) > 30) suppressWarnings(stats::cor(fv[ok], abs(res[ok]), use = "complete.obs")) else NA_real_
  #   
  #   cat("Diagnostics summary (Manual SARIMA)\n")
  #   cat("------------------------------------------------------------\n")
  #   
  #   cat("* White-noise / autocorrelation:\n")
  #   if (!is.null(lb) && is.finite(lb$p.value)) {
  #     cat(
  #       "  - Ljung–Box Q(", lb$parameter, ") = ", round(as.numeric(lb$statistic), 3),
  #       ", p = ", signif(as.numeric(lb$p.value), 3),
  #       if (lb$p.value > ctrl$alpha) "  -> No strong evidence of residual autocorrelation.\n"
  #       else "  -> Evidence of remaining autocorrelation (revise p/q/P/Q or differencing).\n",
  #       sep = ""
  #     )
  #   } else {
  #     cat("  - Ljung–Box unavailable (test error or insufficient df).\n")
  #   }
  #   
  #   cat("* Distribution / normality:\n")
  #   if (!is.null(sw) && is.finite(sw$p.value)) {
  #     cat(
  #       "  - Shapiro–Wilk p = ", signif(sw$p.value, 3),
  #       if (sw$p.value > ctrl$alpha) "  -> residuals not strongly inconsistent with normality.\n"
  #       else "  -> residuals deviate from normality (heavy tails/skew possible).\n",
  #       sep = ""
  #     )
  #   } else {
  #     cat("  - Shapiro–Wilk not computed (N out of range or error).\n")
  #   }
  #   
  #   cat("* Variance stability:\n")
  #   if (is.finite(rho)) {
  #     cat(
  #       "  - Corr(fitted, |res|) = ", round(rho, 3),
  #       if (abs(rho) < 0.15) "  -> variance looks roughly stable.\n"
  #       else "  -> possible heteroskedasticity.\n",
  #       sep = ""
  #     )
  #   } else {
  #     cat("  - Not enough usable points to assess variance trend.\n")
  #   }
  #   
  #   cat("\nConclusion / next steps:\n")
  #   cat("------------------------------------------------------------\n")
  #   cat("- If ACF spikes and Ljung–Box rejects: revisit (p,q)(P,Q) and/or (d,D,s).\n")
  #   cat("- If Q–Q shows heavy tails: intervals may be optimistic; consider robust alternatives.\n")
  #   cat("- If residual-vs-fitted shows structure: consider regressors/interventions.\n")
  # })
  # 
  # # =========================
  # # Ljung–Box p-values across lags (Manual SARIMA)  [REPORT VERSION]
  # # (keep this output name because the Academic conclusion (full) may reference it)
  # # =========================
  # output$manual_resid_lb_pvals <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res   <- residuals(fit)
  #   fitdf <- length(coef(fit))
  #   
  #   ctrl <- get_diag_controls(input)
  #   out  <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
  #   
  #   validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
  #   
  #   plot_lb_pvals(out$pvals, out$L, ctrl$alpha, fitdf,
  #                 main_title = "Ljung–Box p-values by lag (Manual SARIMA)")
  # })
  # 
  # # =========================
  # # Diagnostics tab: dynamic plot layout builder (THIS WAS MISSING)
  # # =========================
  # 
  # manual_diag_plot_map <- list(
  #   ts     = list(id = "manual_resid_ts_diag",       title = "Residuals over time"),
  #   acf    = list(id = "manual_resid_acf_diag",      title = "Residual ACF"),
  #   hist   = list(id = "manual_resid_hist_diag",     title = "Residual histogram"),
  #   qq     = list(id = "manual_resid_qq_diag",       title = "Normal Q-Q"),
  #   fitted = list(id = "manual_resid_fitted_diag",   title = "Residuals vs fitted"),
  #   scale  = list(id = "manual_resid_scale_diag",    title = "Scale-location"),
  #   lb     = list(id = "manual_resid_lb_pvals_diag", title = "Ljung–Box p-values by lag")
  # )
  # 
  # make_diag_plot_output <- function(plot_key, height_px, width_px = 0) {
  #   info <- manual_diag_plot_map[[plot_key]]
  #   if (is.null(info)) return(NULL)
  #   
  #   w <- if (!is.finite(width_px) || width_px <= 0) "100%" else paste0(width_px, "px")
  #   
  #   tagList(
  #     tags$div(
  #       style = "margin-bottom:10px;",
  #       tags$h5(info$title),
  #       plotOutput(info$id, height = height_px, width = w)
  #     )
  #   )
  # }
  # 
  # get_diag_plot_order <- reactive({
  #   sel <- input$diag_plots_selected
  #   if (is.null(sel) || length(sel) == 0) return(character(0))
  #   
  #   preset <- input$diag_layout_preset %||% "preset_current"
  #   
  #   if (preset != "preset_custom") {
  #     default_order <- c("ts", "acf", "hist", "qq", "fitted", "scale", "lb")
  #     return(default_order[default_order %in% sel])
  #   }
  #   
  #   pos <- c(
  #     ts     = input$pos_ts,
  #     acf    = input$pos_acf,
  #     hist   = input$pos_hist,
  #     qq     = input$pos_qq,
  #     fitted = input$pos_fitted,
  #     scale  = input$pos_scale,
  #     lb     = input$pos_lb
  #   )
  #   
  #   pos <- pos[names(pos) %in% sel]
  #   pos_num <- suppressWarnings(as.numeric(pos))
  #   pos_num[!is.finite(pos_num)] <- 999
  #   
  #   ord <- order(pos_num, names(pos_num))
  #   names(pos_num)[ord]
  # })
  # 
  # output$manual_diag_plots_ui <- renderUI({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) first, then diagnostics will appear."))
  #   
  #   keys <- get_diag_plot_order()
  #   if (length(keys) == 0) {
  #     return(tags$div(
  #       style = "padding:10px;border:1px dashed #ccc;border-radius:10px;background:#fff;",
  #       tags$strong("No plots selected."),
  #       tags$p("Use the sidebar to select one or more plots.")
  #     ))
  #   }
  #   
  #   height_px <- input$diag_plot_height %||% 260
  #   width_px  <- input$diag_plot_width  %||% 0
  #   preset    <- input$diag_layout_preset %||% "preset_current"
  #   
  #   # Preset: current layout (2x2 + stacked extras)
  #   if (preset == "preset_current") {
  #     row1 <- intersect(c("ts", "acf"), keys)
  #     row2 <- intersect(c("hist", "qq"), keys)
  #     rest <- setdiff(keys, c("ts", "acf", "hist", "qq"))
  #     
  #     ui <- tagList()
  #     
  #     if (length(row1) > 0) {
  #       ui <- tagList(ui,
  #                     fluidRow(lapply(row1, function(k) column(6, make_diag_plot_output(k, height_px, width_px))))
  #       )
  #     }
  #     if (length(row2) > 0) {
  #       ui <- tagList(ui,
  #                     fluidRow(lapply(row2, function(k) column(6, make_diag_plot_output(k, height_px, width_px))))
  #       )
  #     }
  #     if (length(rest) > 0) {
  #       ui <- tagList(ui, lapply(rest, function(k) make_diag_plot_output(k, height_px, width_px)))
  #     }
  #     return(ui)
  #   }
  #   
  #   # Preset: two columns stacked
  #   if (preset == "preset_two_col") {
  #     left  <- keys[seq(1, length(keys), by = 2)]
  #     right <- keys[seq(2, length(keys), by = 2)]
  #     return(
  #       fluidRow(
  #         column(6, lapply(left,  function(k) make_diag_plot_output(k, height_px, width_px))),
  #         column(6, lapply(right, function(k) make_diag_plot_output(k, height_px, width_px)))
  #       )
  #     )
  #   }
  #   
  #   # Preset: one column stacked (and custom order)
  #   tagList(lapply(keys, function(k) make_diag_plot_output(k, height_px, width_px)))
  # })
  
  
  
  
  
  
  # # =========================
  # # For Manual Diag (FULL REPLACEMENT + UI BUILDER + WIDTH SLIDER)
  # # =========================
  # 
  # `%||%` <- function(a, b) if (!is.null(a)) a else b
  # 
  # nice_par <- function() {
  #   par(mar = c(4, 4, 2.2, 1), mgp = c(2.2, 0.7, 0), las = 1)
  # }
  # 
  # # ---- Shared helpers ----
  # get_diag_controls <- function(input) {
  #   L_input <- suppressWarnings(as.integer(input$diag_lag))
  #   L <- if (is.finite(L_input) && L_input > 0) L_input else 12L
  #   
  #   alpha <- suppressWarnings(as.numeric(input$alphaSt2))
  #   if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05
  #   
  #   list(L = L, alpha = alpha)
  # }
  # 
  # compute_lb_pvals <- function(res, L, fitdf) {
  #   res <- as.numeric(res)
  #   res <- res[is.finite(res)]
  #   
  #   N <- length(res)
  #   if (N < 8) return(list(pvals = numeric(0), L = 0L, N = N))
  #   
  #   L <- as.integer(L)
  #   L <- max(1L, min(L, N - 1L))
  #   
  #   pvals <- rep(NA_real_, L)
  #   for (k in seq_len(L)) {
  #     if (k > fitdf) {
  #       bt <- tryCatch(
  #         stats::Box.test(res, lag = k, type = "Ljung-Box", fitdf = fitdf),
  #         error = function(e) NULL
  #       )
  #       pvals[k] <- if (!is.null(bt)) as.numeric(bt$p.value) else NA_real_
  #     }
  #   }
  #   
  #   list(pvals = pvals, L = L, N = N)
  # }
  # 
  # plot_lb_pvals <- function(pvals, L, alpha, fitdf, main_title) {
  #   plot(seq_len(L), pvals,
  #        type = "h", lwd = 2,
  #        xlab = "Lag (k)",
  #        ylab = "p-value  (Ljung–Box up to lag k)",
  #        main = main_title,
  #        ylim = c(0, 1))
  #   points(seq_len(L), pvals, pch = 16)
  #   abline(h = alpha, lty = 2, col = "gray40")
  #   mtext(sprintf("alpha = %.3f, fitdf = %d", alpha, fitdf),
  #         side = 3, line = 0.2, cex = 0.85)
  #   
  #   if (is.finite(fitdf) && fitdf >= 1 && fitdf < L) {
  #     rect(xleft = 0.5, ybottom = -0.02, xright = fitdf + 0.5, ytop = 1.02,
  #          border = NA, col = grDevices::adjustcolor("gray", alpha.f = 0.15))
  #     text(x = (fitdf + 1) / 2, y = 0.95, labels = "df ≤ 0 (omitted)", cex = 0.8)
  #   }
  # }
  # 
  # # =========================
  # # Manual SARIMA diagnostics tab (NEW independent outputs)
  # # =========================
  # 
  # output$manual_resid_ts_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 3, "Not enough residuals to plot."))
  #   
  #   plot(res, type = "l",
  #        main = "Residuals over time",
  #        xlab = "Time", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  # })
  # 
  # output$manual_resid_acf_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for ACF."))
  #   
  #   plot(acf(res, plot = FALSE), main = "Residual ACF")
  # })
  # 
  # output$manual_resid_hist_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for histogram."))
  #   
  #   hist(res, breaks = 30, col = "gray85", border = "white",
  #        main = "Residual histogram", xlab = "Residual")
  # })
  # 
  # output$manual_resid_qq_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for Q-Q plot."))
  #   
  #   qqnorm(res, main = "Normal Q–Q")
  #   qqline(res, col = "red", lwd = 2)
  # })
  # 
  # output$manual_resid_fitted_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   validate(need(sum(ok) >= 10, "Not enough valid points for residuals vs fitted."))
  #   
  #   plot(fv[ok], res[ok],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("steelblue", alpha.f = 0.7),
  #        main = "Residuals vs fitted",
  #        xlab = "Fitted values", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  #   
  #   if (sum(ok) > 20) {
  #     lines(stats::lowess(fv[ok], res[ok], f = 2/3), col = "tomato", lwd = 2)
  #   }
  # })
  # 
  # output$manual_resid_scale_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   
  #   y <- sqrt(abs(res))
  #   ok2 <- ok & is.finite(y)
  #   validate(need(sum(ok2) >= 10, "Not enough valid points for scale-location plot."))
  #   
  #   plot(fv[ok2], y[ok2],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("darkgreen", alpha.f = 0.65),
  #        main = "Scale-location (sqrt(|res|) vs fitted)",
  #        xlab = "Fitted values", ylab = expression(sqrt("|Residual|")))
  #   if (sum(ok2) > 20) {
  #     lines(stats::lowess(fv[ok2], y[ok2], f = 2/3), col = "tomato", lwd = 2)
  #   }
  # })
  # 
  # output$manual_resid_lb_pvals_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res   <- residuals(fit)
  #   fitdf <- length(coef(fit))
  #   
  #   ctrl <- get_diag_controls(input)
  #   out  <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
  #   
  #   validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
  #   
  #   plot_lb_pvals(out$pvals, out$L, ctrl$alpha, fitdf,
  #                 main_title = "Ljung–Box p-values by lag (Diagnostics)")
  # })
  # 
  # output$manual_diag_commentary_diag <- renderPrint({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   N   <- length(res)
  #   
  #   ctrl  <- get_diag_controls(input)
  #   fitdf <- length(coef(fit))
  #   
  #   if (N < 8) {
  #     cat(
  #       "Diagnostics summary (Manual SARIMA)\n",
  #       "------------------------------------------------------------\n",
  #       "- Too few residuals to run a stable diagnostic battery.\n",
  #       "- Fit the model on more observations or reduce differencing/orders.\n",
  #       sep = ""
  #     )
  #     return(invisible())
  #   }
  #   
  #   L2 <- max(1L, min(as.integer(ctrl$L), N - 1L))
  #   lb <- tryCatch(
  #     stats::Box.test(res, lag = L2, type = "Ljung-Box", fitdf = fitdf),
  #     error = function(e) NULL
  #   )
  #   
  #   sw <- if (N >= 3 && N <= 5000) tryCatch(stats::shapiro.test(res), error = function(e) NULL) else NULL
  #   
  #   fv <- as.numeric(fitted(fit))
  #   ok <- is.finite(fv) & is.finite(res)
  #   rho <- if (sum(ok) > 30) suppressWarnings(stats::cor(fv[ok], abs(res[ok]), use = "complete.obs")) else NA_real_
  #   
  #   cat("Diagnostics summary (Manual SARIMA)\n")
  #   cat("------------------------------------------------------------\n")
  #   
  #   cat("* White-noise / autocorrelation:\n")
  #   if (!is.null(lb) && is.finite(lb$p.value)) {
  #     cat(
  #       "  - Ljung–Box Q(", lb$parameter, ") = ", round(as.numeric(lb$statistic), 3),
  #       ", p = ", signif(as.numeric(lb$p.value), 3),
  #       if (lb$p.value > ctrl$alpha) "  -> No strong evidence of residual autocorrelation.\n"
  #       else "  -> Evidence of remaining autocorrelation (revise p/q/P/Q or differencing).\n",
  #       sep = ""
  #     )
  #   } else {
  #     cat("  - Ljung–Box unavailable (test error or insufficient df).\n")
  #   }
  #   
  #   cat("* Distribution / normality:\n")
  #   if (!is.null(sw) && is.finite(sw$p.value)) {
  #     cat(
  #       "  - Shapiro–Wilk p = ", signif(sw$p.value, 3),
  #       if (sw$p.value > ctrl$alpha) "  -> residuals not strongly inconsistent with normality.\n"
  #       else "  -> residuals deviate from normality (heavy tails/skew possible).\n",
  #       sep = ""
  #     )
  #   } else {
  #     cat("  - Shapiro–Wilk not computed (N out of range or error).\n")
  #   }
  #   
  #   cat("* Variance stability:\n")
  #   if (is.finite(rho)) {
  #     cat(
  #       "  - Corr(fitted, |res|) = ", round(rho, 3),
  #       if (abs(rho) < 0.15) "  -> variance looks roughly stable.\n"
  #       else "  -> possible heteroskedasticity.\n",
  #       sep = ""
  #     )
  #   } else {
  #     cat("  - Not enough usable points to assess variance trend.\n")
  #   }
  #   
  #   cat("\nConclusion / next steps:\n")
  #   cat("------------------------------------------------------------\n")
  #   cat("- If ACF spikes and Ljung–Box rejects: revisit (p,q)(P,Q) and/or (d,D,s).\n")
  #   cat("- If Q–Q shows heavy tails: intervals may be optimistic; consider robust alternatives.\n")
  #   cat("- If residual-vs-fitted shows structure: consider regressors/interventions.\n")
  # })
  # 
  # # =========================
  # # Dynamic main panel (controls width of the plots panel)
  # # =========================
  # 
  # output$manual_diag_mainpanel_ui <- renderUI({
  #   w <- input$diag_main_width %||% 9
  #   w <- as.integer(w)
  #   if (!is.finite(w) || w < 6) w <- 6
  #   if (w > 11) w <- 11
  #   
  #   mainPanel(
  #     width = w,
  #     
  #     tags$div(
  #       style = "padding:10px;border:1px solid #e5e5e5;border-radius:10px;background:#fcfcfc;margin-bottom:12px;",
  #       tags$h4("Residual diagnostics (Manual SARIMA)"),
  #       tags$p("Use the controls on the left to choose plots, layout, sizing, and whether to show the academic conclusion.")
  #     ),
  #     
  #     uiOutput("manual_diag_plots_ui"),
  #     
  #     conditionalPanel(
  #       condition = "input.diag_show_conclusion == true",
  #       tags$h4("Academic conclusion (Diagnostics)"),
  #       tags$div(
  #         style = "padding:10px;border:1px solid #e5e5e5;border-radius:10px;background:#ffffff;",
  #         verbatimTextOutput("manual_diag_commentary_diag")
  #       )
  #     )
  #   )
  # })
  # 
  # # =========================
  # # Diagnostics tab: dynamic plot layout builder
  # # =========================
  # 
  # manual_diag_plot_map <- list(
  #   ts     = list(id = "manual_resid_ts_diag",       title = "Residuals over time"),
  #   acf    = list(id = "manual_resid_acf_diag",      title = "Residual ACF"),
  #   hist   = list(id = "manual_resid_hist_diag",     title = "Residual histogram"),
  #   qq     = list(id = "manual_resid_qq_diag",       title = "Normal Q-Q"),
  #   fitted = list(id = "manual_resid_fitted_diag",   title = "Residuals vs fitted"),
  #   scale  = list(id = "manual_resid_scale_diag",    title = "Scale-location"),
  #   lb     = list(id = "manual_resid_lb_pvals_diag", title = "Ljung–Box p-values by lag")
  # )
  # 
  # # Create a plot "card" that never overlaps other columns
  # make_diag_plot_output <- function(plot_key, height_px, width_px = 0) {
  #   info <- manual_diag_plot_map[[plot_key]]
  #   if (is.null(info)) return(NULL)
  #   
  #   # Prefer 100% width (responsive). If px requested, still constrain to container.
  #   w <- if (!is.finite(width_px) || width_px <= 0) "100%" else paste0(width_px, "px")
  #   
  #   tags$div(
  #     style = "margin-bottom:12px;",
  #     tags$h5(info$title),
  #     tags$div(
  #       style = "width:100%;max-width:100%;overflow:hidden;",  # prevents overlap
  #       plotOutput(info$id, height = height_px, width = w)
  #     )
  #   )
  # }
  # 
  # get_diag_plot_order <- reactive({
  #   sel <- input$diag_plots_selected
  #   if (is.null(sel) || length(sel) == 0) return(character(0))
  #   
  #   preset <- input$diag_layout_preset %||% "preset_current"
  #   
  #   if (preset != "preset_custom") {
  #     default_order <- c("ts", "acf", "hist", "qq", "fitted", "scale", "lb")
  #     return(default_order[default_order %in% sel])
  #   }
  #   
  #   pos <- c(
  #     ts     = input$pos_ts,
  #     acf    = input$pos_acf,
  #     hist   = input$pos_hist,
  #     qq     = input$pos_qq,
  #     fitted = input$pos_fitted,
  #     scale  = input$pos_scale,
  #     lb     = input$pos_lb
  #   )
  #   
  #   pos <- pos[names(pos) %in% sel]
  #   pos_num <- suppressWarnings(as.numeric(pos))
  #   pos_num[!is.finite(pos_num)] <- 999
  #   
  #   ord <- order(pos_num, names(pos_num))
  #   names(pos_num)[ord]
  # })
  # 
  # output$manual_diag_plots_ui <- renderUI({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) first, then diagnostics will appear."))
  #   
  #   keys <- get_diag_plot_order()
  #   if (length(keys) == 0) {
  #     return(tags$div(
  #       style = "padding:10px;border:1px dashed #ccc;border-radius:10px;background:#fff;",
  #       tags$strong("No plots selected."),
  #       tags$p("Use the sidebar to select one or more plots.")
  #     ))
  #   }
  #   
  #   height_px <- input$diag_plot_height %||% 260
  #   width_px  <- input$diag_plot_width  %||% 0
  #   preset    <- input$diag_layout_preset %||% "preset_current"
  #   
  #   # --- Preset: current layout (2x2 + stacked extras)
  #   if (preset == "preset_current") {
  #     row1 <- intersect(c("ts", "acf"), keys)
  #     row2 <- intersect(c("hist", "qq"), keys)
  #     rest <- setdiff(keys, c("ts", "acf", "hist", "qq"))
  #     
  #     ui <- tagList()
  #     
  #     if (length(row1) > 0) {
  #       ui <- tagList(ui,
  #                     fluidRow(lapply(row1, function(k) column(6, make_diag_plot_output(k, height_px, width_px))))
  #       )
  #     }
  #     if (length(row2) > 0) {
  #       ui <- tagList(ui,
  #                     fluidRow(lapply(row2, function(k) column(6, make_diag_plot_output(k, height_px, width_px))))
  #       )
  #     }
  #     if (length(rest) > 0) {
  #       ui <- tagList(ui, lapply(rest, function(k) make_diag_plot_output(k, height_px, width_px)))
  #     }
  #     return(ui)
  #   }
  #   
  #   # --- Preset: two columns stacked
  #   if (preset == "preset_two_col") {
  #     left  <- keys[seq(1, length(keys), by = 2)]
  #     right <- keys[seq(2, length(keys), by = 2)]
  #     return(
  #       fluidRow(
  #         column(6, lapply(left,  function(k) make_diag_plot_output(k, height_px, width_px))),
  #         column(6, lapply(right, function(k) make_diag_plot_output(k, height_px, width_px)))
  #       )
  #     )
  #   }
  #   
  #   # --- Preset: one column stacked (and custom order)
  #   tagList(lapply(keys, function(k) make_diag_plot_output(k, height_px, width_px)))
  # })
  # 
  
  
  
  # # =========================
  # # For Manual Diag (FULL REPLACEMENT + STAY-RIGHT WIDTHS)
  # # =========================
  # 
  # `%||%` <- function(a, b) if (!is.null(a)) a else b
  # 
  # nice_par <- function() {
  #   par(mar = c(4, 4, 2.2, 1), mgp = c(2.2, 0.7, 0), las = 1)
  # }
  # 
  # # ---- Shared helpers ----
  # get_diag_controls <- function(input) {
  #   L_input <- suppressWarnings(as.integer(input$diag_lag))
  #   L <- if (is.finite(L_input) && L_input > 0) L_input else 12L
  #   
  #   alpha <- suppressWarnings(as.numeric(input$alphaSt2))
  #   if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05
  #   
  #   list(L = L, alpha = alpha)
  # }
  # 
  # compute_lb_pvals <- function(res, L, fitdf) {
  #   res <- as.numeric(res)
  #   res <- res[is.finite(res)]
  #   
  #   N <- length(res)
  #   if (N < 8) return(list(pvals = numeric(0), L = 0L, N = N))
  #   
  #   L <- as.integer(L)
  #   L <- max(1L, min(L, N - 1L))
  #   
  #   pvals <- rep(NA_real_, L)
  #   for (k in seq_len(L)) {
  #     if (k > fitdf) {
  #       bt <- tryCatch(
  #         stats::Box.test(res, lag = k, type = "Ljung-Box", fitdf = fitdf),
  #         error = function(e) NULL
  #       )
  #       pvals[k] <- if (!is.null(bt)) as.numeric(bt$p.value) else NA_real_
  #     }
  #   }
  #   
  #   list(pvals = pvals, L = L, N = N)
  # }
  # 
  # plot_lb_pvals <- function(pvals, L, alpha, fitdf, main_title) {
  #   plot(seq_len(L), pvals,
  #        type = "h", lwd = 2,
  #        xlab = "Lag (k)",
  #        ylab = "p-value  (Ljung–Box up to lag k)",
  #        main = main_title,
  #        ylim = c(0, 1))
  #   points(seq_len(L), pvals, pch = 16)
  #   abline(h = alpha, lty = 2, col = "gray40")
  #   mtext(sprintf("alpha = %.3f, fitdf = %d", alpha, fitdf),
  #         side = 3, line = 0.2, cex = 0.85)
  #   
  #   if (is.finite(fitdf) && fitdf >= 1 && fitdf < L) {
  #     rect(xleft = 0.5, ybottom = -0.02, xright = fitdf + 0.5, ytop = 1.02,
  #          border = NA, col = grDevices::adjustcolor("gray", alpha.f = 0.15))
  #     text(x = (fitdf + 1) / 2, y = 0.95, labels = "df ≤ 0 (omitted)", cex = 0.8)
  #   }
  # }
  # 
  # # =========================
  # # Dynamic widths (KEEP ON SAME ROW)
  # # =========================
  # diag_main_w <- reactive({
  #   # main width slider in Bootstrap columns (1..11)
  #   w <- input$diag_main_width %||% 9
  #   w <- as.integer(w)
  #   if (!is.finite(w)) w <- 9
  #   w <- max(6L, min(11L, w))  # keep main panel usable
  #   w
  # })
  # 
  # diag_sidebar_w <- reactive({
  #   # Always keep sum = 12 so it never wraps
  #   12L - diag_main_w()
  # })
  # 
  # # =========================
  # # Sidebar UI (built in server so width can react)
  # # =========================
  # output$manual_diag_sidebar_ui <- renderUI({
  #   column(
  #     width = diag_sidebar_w(),
  #     
  #     tags$h4("Diagnostics layout builder"),
  #     
  #     selectInput(
  #       "diag_layout_preset",
  #       "Layout preset",
  #       choices = c(
  #         "Current layout (2x2 + Ljung-Box + Conclusion)" = "preset_current",
  #         "Two columns (all plots, stacked)"             = "preset_two_col",
  #         "Single column (all plots, stacked)"           = "preset_one_col",
  #         "Custom order (choose positions below)"        = "preset_custom"
  #       ),
  #       selected = "preset_current"
  #     ),
  #     
  #     tags$hr(),
  #     
  #     checkboxGroupInput(
  #       "diag_plots_selected",
  #       "Plots to display",
  #       choices = c(
  #         "Residuals over time"                          = "ts",
  #         "Residual ACF"                                 = "acf",
  #         "Residual histogram"                           = "hist",
  #         "Normal Q-Q"                                   = "qq",
  #         "Residuals vs fitted"                          = "fitted",
  #         "Scale-location (sqrt(|res|) vs fitted)"       = "scale",
  #         "Ljung–Box p-values by lag"                    = "lb"
  #       ),
  #       selected = c("ts", "acf", "hist", "qq", "lb", "fitted", "scale")
  #     ),
  #     
  #     tags$hr(),
  #     
  #     # Slider that controls MAIN PANEL width (sidebar adjusts automatically)
  #     sliderInput(
  #       "diag_main_width",
  #       "Plots panel width (Bootstrap columns)",
  #       min = 6, max = 11, value = 9, step = 1
  #     ),
  #     helpText("Main + sidebar always sum to 12, so the plots panel stays on the right."),
  #     
  #     tags$hr(),
  #     
  #     sliderInput(
  #       "diag_plot_height",
  #       "Plot height (px)",
  #       min = 180, max = 700, value = 260, step = 10
  #     ),
  #     
  #     sliderInput(
  #       "diag_plot_width",
  #       "Plot width (px, optional; constrained)",
  #       min = 0, max = 2000, value = 0, step = 25
  #     ),
  #     helpText("width = 0 uses full available width. Large widths are clipped to avoid overlap."),
  #     
  #     tags$hr(),
  #     
  #     checkboxInput("diag_show_conclusion", "Show academic conclusion panel", value = TRUE),
  #     
  #     conditionalPanel(
  #       condition = "input.diag_layout_preset == 'preset_custom'",
  #       tags$h5("Custom positions (1 = first)"),
  #       helpText("Give each enabled plot a position. Ties are broken by name."),
  #       
  #       numericInput("pos_ts",     "Residuals over time position", value = 1, min = 1, step = 1),
  #       numericInput("pos_acf",    "Residual ACF position",        value = 2, min = 1, step = 1),
  #       numericInput("pos_hist",   "Histogram position",           value = 3, min = 1, step = 1),
  #       numericInput("pos_qq",     "Q-Q position",                 value = 4, min = 1, step = 1),
  #       numericInput("pos_fitted", "Residuals vs fitted position", value = 5, min = 1, step = 1),
  #       numericInput("pos_scale",  "Scale-location position",      value = 6, min = 1, step = 1),
  #       numericInput("pos_lb",     "Ljung–Box p-values position",  value = 7, min = 1, step = 1)
  #     )
  #   )
  # })
  # 
  # # =========================
  # # Main panel UI (built in server so width can react)
  # # =========================
  # output$manual_diag_mainpanel_ui <- renderUI({
  #   column(
  #     width = diag_main_w(),
  #     
  #     tags$div(
  #       style = "padding:10px;border:1px solid #e5e5e5;border-radius:10px;background:#fcfcfc;margin-bottom:12px;",
  #       tags$h4("Residual diagnostics (Manual SARIMA)"),
  #       tags$p("Choose plots, layout, and sizing from the left panel.")
  #     ),
  #     
  #     uiOutput("manual_diag_plots_ui"),
  #     
  #     conditionalPanel(
  #       condition = "input.diag_show_conclusion == true",
  #       tags$h4("Academic conclusion (Diagnostics)"),
  #       tags$div(
  #         style = "padding:10px;border:1px solid #e5e5e5;border-radius:10px;background:#ffffff;",
  #         verbatimTextOutput("manual_diag_commentary_diag")
  #       )
  #     )
  #   )
  # })
  # 
  # # =========================
  # # Plot outputs (diag)
  # # =========================
  # 
  # output$manual_resid_ts_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 3, "Not enough residuals to plot."))
  #   
  #   plot(res, type = "l", main = "Residuals over time", xlab = "Time", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  # })
  # 
  # output$manual_resid_acf_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for ACF."))
  #   
  #   plot(acf(res, plot = FALSE), main = "Residual ACF")
  # })
  # 
  # output$manual_resid_hist_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for histogram."))
  #   
  #   hist(res, breaks = 30, col = "gray85", border = "white",
  #        main = "Residual histogram", xlab = "Residual")
  # })
  # 
  # output$manual_resid_qq_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for Q-Q plot."))
  #   
  #   qqnorm(res, main = "Normal Q–Q")
  #   qqline(res, col = "red", lwd = 2)
  # })
  # 
  # output$manual_resid_fitted_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   validate(need(sum(ok) >= 10, "Not enough valid points for residuals vs fitted."))
  #   
  #   plot(fv[ok], res[ok],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("steelblue", alpha.f = 0.7),
  #        main = "Residuals vs fitted",
  #        xlab = "Fitted values", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  #   if (sum(ok) > 20) lines(stats::lowess(fv[ok], res[ok], f = 2/3), col = "tomato", lwd = 2)
  # })
  # 
  # output$manual_resid_scale_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   
  #   y <- sqrt(abs(res))
  #   ok2 <- ok & is.finite(y)
  #   validate(need(sum(ok2) >= 10, "Not enough valid points for scale-location plot."))
  #   
  #   plot(fv[ok2], y[ok2],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("darkgreen", alpha.f = 0.65),
  #        main = "Scale-location (sqrt(|res|) vs fitted)",
  #        xlab = "Fitted values", ylab = expression(sqrt("|Residual|")))
  #   if (sum(ok2) > 20) lines(stats::lowess(fv[ok2], y[ok2], f = 2/3), col = "tomato", lwd = 2)
  # })
  # 
  # output$manual_resid_lb_pvals_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res   <- residuals(fit)
  #   fitdf <- length(coef(fit))
  #   
  #   ctrl <- get_diag_controls(input)
  #   out  <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
  #   
  #   validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
  #   
  #   plot_lb_pvals(out$pvals, out$L, ctrl$alpha, fitdf, "Ljung–Box p-values by lag (Diagnostics)")
  # })
  # 
  # output$manual_diag_commentary_diag <- renderPrint({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   N   <- length(res)
  #   
  #   ctrl  <- get_diag_controls(input)
  #   fitdf <- length(coef(fit))
  #   
  #   if (N < 8) {
  #     cat(
  #       "Diagnostics summary (Manual SARIMA)\n",
  #       "------------------------------------------------------------\n",
  #       "- Too few residuals to run a stable diagnostic battery.\n",
  #       "- Fit the model on more observations or reduce differencing/orders.\n",
  #       sep = ""
  #     )
  #     return(invisible())
  #   }
  #   
  #   L2 <- max(1L, min(as.integer(ctrl$L), N - 1L))
  #   lb <- tryCatch(stats::Box.test(res, lag = L2, type = "Ljung-Box", fitdf = fitdf),
  #                  error = function(e) NULL)
  #   
  #   sw <- if (N >= 3 && N <= 5000) tryCatch(stats::shapiro.test(res), error = function(e) NULL) else NULL
  #   
  #   fv <- as.numeric(fitted(fit))
  #   ok <- is.finite(fv) & is.finite(res)
  #   rho <- if (sum(ok) > 30) suppressWarnings(stats::cor(fv[ok], abs(res[ok]), use = "complete.obs")) else NA_real_
  #   
  #   cat("Diagnostics summary (Manual SARIMA)\n")
  #   cat("------------------------------------------------------------\n")
  #   
  #   cat("* White-noise / autocorrelation:\n")
  #   if (!is.null(lb) && is.finite(lb$p.value)) {
  #     cat("  - Ljung–Box p = ", signif(lb$p.value, 3),
  #         if (lb$p.value > ctrl$alpha) " -> no strong evidence of autocorrelation.\n"
  #         else " -> remaining autocorrelation; consider revising orders/differencing.\n", sep = "")
  #   } else {
  #     cat("  - Ljung–Box unavailable.\n")
  #   }
  #   
  #   cat("* Distribution / normality:\n")
  #   if (!is.null(sw) && is.finite(sw$p.value)) {
  #     cat("  - Shapiro–Wilk p = ", signif(sw$p.value, 3),
  #         if (sw$p.value > ctrl$alpha) " -> not strongly inconsistent with normality.\n"
  #         else " -> non-normal residuals likely.\n", sep = "")
  #   } else {
  #     cat("  - Shapiro–Wilk not computed.\n")
  #   }
  #   
  #   cat("* Variance stability:\n")
  #   if (is.finite(rho)) {
  #     cat("  - Corr(fitted, |res|) = ", round(rho, 3),
  #         if (abs(rho) < 0.15) " -> roughly stable.\n"
  #         else " -> possible heteroskedasticity.\n", sep = "")
  #   } else {
  #     cat("  - Not enough points to assess.\n")
  #   }
  # })
  # 
  # # =========================
  # # Diagnostics tab: dynamic plot layout builder
  # # =========================
  # 
  # manual_diag_plot_map <- list(
  #   ts     = list(id = "manual_resid_ts_diag",       title = "Residuals over time"),
  #   acf    = list(id = "manual_resid_acf_diag",      title = "Residual ACF"),
  #   hist   = list(id = "manual_resid_hist_diag",     title = "Residual histogram"),
  #   qq     = list(id = "manual_resid_qq_diag",       title = "Normal Q-Q"),
  #   fitted = list(id = "manual_resid_fitted_diag",   title = "Residuals vs fitted"),
  #   scale  = list(id = "manual_resid_scale_diag",    title = "Scale-location"),
  #   lb     = list(id = "manual_resid_lb_pvals_diag", title = "Ljung–Box p-values by lag")
  # )
  # 
  # make_diag_plot_output <- function(plot_key, height_px, width_px = 0) {
  #   info <- manual_diag_plot_map[[plot_key]]
  #   if (is.null(info)) return(NULL)
  #   
  #   # responsive by default
  #   w <- if (!is.finite(width_px) || width_px <= 0) "100%" else paste0(width_px, "px")
  #   
  #   tags$div(
  #     style = "margin-bottom:12px;",
  #     tags$h5(info$title),
  #     tags$div(
  #       style = "width:100%;max-width:100%;overflow:hidden;", # prevents overlap
  #       plotOutput(info$id, height = height_px, width = w)
  #     )
  #   )
  # }
  # 
  # get_diag_plot_order <- reactive({
  #   sel <- input$diag_plots_selected
  #   if (is.null(sel) || length(sel) == 0) return(character(0))
  #   
  #   preset <- input$diag_layout_preset %||% "preset_current"
  #   if (preset != "preset_custom") {
  #     default_order <- c("ts", "acf", "hist", "qq", "fitted", "scale", "lb")
  #     return(default_order[default_order %in% sel])
  #   }
  #   
  #   pos <- c(
  #     ts     = input$pos_ts,
  #     acf    = input$pos_acf,
  #     hist   = input$pos_hist,
  #     qq     = input$pos_qq,
  #     fitted = input$pos_fitted,
  #     scale  = input$pos_scale,
  #     lb     = input$pos_lb
  #   )
  #   pos <- pos[names(pos) %in% sel]
  #   pos_num <- suppressWarnings(as.numeric(pos))
  #   pos_num[!is.finite(pos_num)] <- 999
  #   ord <- order(pos_num, names(pos_num))
  #   names(pos_num)[ord]
  # })
  # 
  # output$manual_diag_plots_ui <- renderUI({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) first, then diagnostics will appear."))
  #   
  #   keys <- get_diag_plot_order()
  #   if (length(keys) == 0) {
  #     return(tags$div(
  #       style = "padding:10px;border:1px dashed #ccc;border-radius:10px;background:#fff;",
  #       tags$strong("No plots selected."),
  #       tags$p("Use the left panel to select plots.")
  #     ))
  #   }
  #   
  #   height_px <- input$diag_plot_height %||% 260
  #   width_px  <- input$diag_plot_width  %||% 0
  #   preset    <- input$diag_layout_preset %||% "preset_current"
  #   
  #   if (preset == "preset_current") {
  #     row1 <- intersect(c("ts", "acf"), keys)
  #     row2 <- intersect(c("hist", "qq"), keys)
  #     rest <- setdiff(keys, c("ts", "acf", "hist", "qq"))
  #     
  #     ui <- tagList()
  #     if (length(row1) > 0) ui <- tagList(ui, fluidRow(lapply(row1, function(k) column(6, make_diag_plot_output(k, height_px, width_px)))))
  #     if (length(row2) > 0) ui <- tagList(ui, fluidRow(lapply(row2, function(k) column(6, make_diag_plot_output(k, height_px, width_px)))))
  #     if (length(rest) > 0) ui <- tagList(ui, lapply(rest, function(k) make_diag_plot_output(k, height_px, width_px)))
  #     return(ui)
  #   }
  #   
  #   if (preset == "preset_two_col") {
  #     left  <- keys[seq(1, length(keys), by = 2)]
  #     right <- keys[seq(2, length(keys), by = 2)]
  #     return(
  #       fluidRow(
  #         column(6, lapply(left,  function(k) make_diag_plot_output(k, height_px, width_px))),
  #         column(6, lapply(right, function(k) make_diag_plot_output(k, height_px, width_px)))
  #       )
  #     )
  #   }
  #   
  #   tagList(lapply(keys, function(k) make_diag_plot_output(k, height_px, width_px)))
  # })
  
  
 
  
  
  
   
  # # =========================
  # # Manual SARIMA Diagnostics (FULL) - responsive panel + builder
  # # =========================
  # 
  # # `%||%` <- function(a, b) if (!is.null(a)) a else b
  # 
  # 
  # # main width from slider
  # diag_main_w <- reactive({
  #   w <- input$diag_main_width %||% 9
  #   w <- as.integer(w)
  #   if (!is.finite(w)) w <- 9
  #   w <- max(6L, min(11L, w))
  #   w
  # })
  # 
  # # sidebar width auto-adjust so main stays right (sum = 12)
  # diag_sidebar_w <- reactive({
  #   12L - diag_main_w()
  # })
  # 
  # output$manual_diag_sidebar_ui <- renderUI({
  #   column(
  #     width = diag_sidebar_w(),
  # 
  #     tags$h4("Diagnostics layout builder"),
  # 
  #     selectInput(
  #       "diag_layout_preset",
  #       "Layout preset",
  #       choices = c(
  #         "Current layout (2x2 + Ljung-Box + Conclusion)" = "preset_current",
  #         "Two columns (all plots, stacked)"             = "preset_two_col",
  #         "Single column (all plots, stacked)"           = "preset_one_col",
  #         "Custom order (choose positions below)"        = "preset_custom"
  #       ),
  #       selected = "preset_current"
  #     ),
  # 
  #     tags$hr(),
  # 
  #     checkboxGroupInput(
  #       "diag_plots_selected",
  #       "Plots to display",
  #       choices = c(
  #         "Residuals over time"                          = "ts",
  #         "Residual ACF"                                 = "acf",
  #         "Residual histogram"                           = "hist",
  #         "Normal Q-Q"                                   = "qq",
  #         "Residuals vs fitted"                          = "fitted",
  #         "Scale-location (sqrt(|res|) vs fitted)"       = "scale",
  #         "Ljung–Box p-values by lag"                    = "lb"
  #       ),
  #       selected = c("ts", "acf", "hist", "qq", "lb", "fitted", "scale")
  #     ),
  # 
  #     tags$hr(),
  # 
  #     sliderInput(
  #       "diag_main_width",
  #       "Plots panel width (Bootstrap columns)",
  #       min = 6, max = 11, value = 9, step = 1
  #     ),
  #     helpText("Main + sidebar always sum to 12, so the plots panel stays on the right."),
  # 
  #     tags$hr(),
  # 
  #     sliderInput("diag_plot_height", "Plot height (px)", min = 180, max = 700, value = 260, step = 10),
  #     sliderInput("diag_plot_width",  "Plot width (px, optional)", min = 0, max = 2000, value = 0, step = 25),
  #     helpText("width=0 uses full available width; large widths are constrained in the plot cards."),
  # 
  #     tags$hr(),
  # 
  #     checkboxInput("diag_show_conclusion", "Show academic conclusion panel", value = TRUE),
  # 
  #     conditionalPanel(
  #       condition = "input.diag_layout_preset == 'preset_custom'",
  #       tags$h5("Custom positions (1 = first)"),
  #       numericInput("pos_ts",     "Residuals over time position", value = 1, min = 1),
  #       numericInput("pos_acf",    "Residual ACF position",        value = 2, min = 1),
  #       numericInput("pos_hist",   "Histogram position",           value = 3, min = 1),
  #       numericInput("pos_qq",     "Q-Q position",                 value = 4, min = 1),
  #       numericInput("pos_fitted", "Residuals vs fitted position", value = 5, min = 1),
  #       numericInput("pos_scale",  "Scale-location position",      value = 6, min = 1),
  #       numericInput("pos_lb",     "Ljung–Box p-values position",  value = 7, min = 1)
  #     )
  #   )
  # })
  # 
  # output$manual_diag_mainpanel_ui <- renderUI({
  #   column(
  #     width = diag_main_w(),
  # 
  #     tags$div(
  #       style = "padding:10px;border:1px solid #e5e5e5;border-radius:10px;background:#fcfcfc;margin-bottom:12px;",
  #       tags$h4("Residual diagnostics (Manual SARIMA)"),
  #       tags$p("Use the controls on the left to choose plots, layout, size, and whether to show the academic conclusion.")
  #     ),
  # 
  #     uiOutput("manual_diag_plots_ui"),
  # 
  #     conditionalPanel(
  #       condition = "input.diag_show_conclusion == true",
  #       tags$h4("Academic conclusion (Diagnostics)"),
  #       tags$div(
  #         style = "padding:10px;border:1px solid #e5e5e5;border-radius:10px;background:#ffffff;",
  #         verbatimTextOutput("manual_diag_commentary_diag")
  #       )
  #     )
  #   )
  # })
  # 
  # 
  # 
  # 
  # 
  # 
  # nice_par <- function() {
  #   par(mar = c(4, 4, 2.2, 1), mgp = c(2.2, 0.7, 0), las = 1)
  # }
  # 
  # # ---- Shared helpers ----
  # get_diag_controls <- function(input) {
  #   L_input <- suppressWarnings(as.integer(input$diag_lag))
  #   L <- if (is.finite(L_input) && L_input > 0) L_input else 12L
  # 
  #   alpha <- suppressWarnings(as.numeric(input$alphaSt2))
  #   if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05
  # 
  #   list(L = L, alpha = alpha)
  # }
  # 
  # compute_lb_pvals <- function(res, L, fitdf) {
  #   res <- as.numeric(res)
  #   res <- res[is.finite(res)]
  # 
  #   N <- length(res)
  #   if (N < 8) return(list(pvals = numeric(0), L = 0L, N = N))
  # 
  #   L <- as.integer(L)
  #   L <- max(1L, min(L, N - 1L))
  # 
  #   pvals <- rep(NA_real_, L)
  #   for (k in seq_len(L)) {
  #     if (k > fitdf) {
  #       bt <- tryCatch(
  #         stats::Box.test(res, lag = k, type = "Ljung-Box", fitdf = fitdf),
  #         error = function(e) NULL
  #       )
  #       pvals[k] <- if (!is.null(bt)) as.numeric(bt$p.value) else NA_real_
  #     }
  #   }
  #   list(pvals = pvals, L = L, N = N)
  # }
  # 
  # plot_lb_pvals <- function(pvals, L, alpha, fitdf, main_title) {
  #   plot(seq_len(L), pvals,
  #        type = "h", lwd = 2,
  #        xlab = "Lag (k)",
  #        ylab = "p-value  (Ljung–Box up to lag k)",
  #        main = main_title,
  #        ylim = c(0, 1))
  #   points(seq_len(L), pvals, pch = 16)
  #   abline(h = alpha, lty = 2, col = "gray40")
  #   mtext(sprintf("alpha = %.3f, fitdf = %d", alpha, fitdf),
  #         side = 3, line = 0.2, cex = 0.85)
  # 
  #   if (is.finite(fitdf) && fitdf >= 1 && fitdf < L) {
  #     rect(xleft = 0.5, ybottom = -0.02, xright = fitdf + 0.5, ytop = 1.02,
  #          border = NA, col = grDevices::adjustcolor("gray", alpha.f = 0.15))
  #     text(x = (fitdf + 1) / 2, y = 0.95, labels = "df ≤ 0 (omitted)", cex = 0.8)
  #   }
  # }
  # 
  # # =========================
  # # Dynamic widths (KEEP MAIN ON RIGHT, NEVER WRAP)
  # # main + sidebar always sum to 12
  # # =========================
  # diag_main_w <- reactive({
  #   w <- input$diag_main_width %||% 9
  #   w <- as.integer(w)
  #   if (!is.finite(w)) w <- 9
  #   w <- max(6L, min(11L, w))
  #   w
  # })
  # 
  # diag_sidebar_w <- reactive({
  #   12L - diag_main_w()
  # })
  # 
  # # =========================
  # # Sidebar UI (dynamic width)
  # # =========================
  # output$manual_diag_sidebar_ui <- renderUI({
  #   column(
  #     width = diag_sidebar_w(),
  # 
  #     tags$h4("Diagnostics layout builder"),
  # 
  #     selectInput(
  #       "diag_layout_preset",
  #       "Layout preset",
  #       choices = c(
  #         "Current layout (2x2 + Ljung-Box + Conclusion)" = "preset_current",
  #         "Two columns (all plots, stacked)"             = "preset_two_col",
  #         "Single column (all plots, stacked)"           = "preset_one_col",
  #         "Custom order (choose positions below)"        = "preset_custom"
  #       ),
  #       selected = "preset_current"
  #     ),
  # 
  #     tags$hr(),
  # 
  #     checkboxGroupInput(
  #       "diag_plots_selected",
  #       "Plots to display",
  #       choices = c(
  #         "Residuals over time"                          = "ts",
  #         "Residual ACF"                                 = "acf",
  #         "Residual histogram"                           = "hist",
  #         "Normal Q-Q"                                   = "qq",
  #         "Residuals vs fitted"                          = "fitted",
  #         "Scale-location (sqrt(|res|) vs fitted)"       = "scale",
  #         "Ljung–Box p-values by lag"                    = "lb"
  #       ),
  #       selected = c("ts", "acf", "hist", "qq", "lb", "fitted", "scale")
  #     ),
  # 
  #     tags$hr(),
  # 
  #     # Controls MAIN panel width; sidebar auto-adjusts so it stays on the right
  #     sliderInput(
  #       "diag_main_width",
  #       "Plots panel width (Bootstrap columns)",
  #       min = 6, max = 11, value = 9, step = 1
  #     ),
  #     helpText("Main + sidebar always sum to 12, so plots panel stays on the right."),
  # 
  #     tags$hr(),
  # 
  #     sliderInput("diag_plot_height", "Plot height (px)", min = 180, max = 700, value = 260, step = 10),
  # 
  #     # Plot width request. Used to decide 1-col vs 2-col to avoid overlap.
  #     sliderInput(
  #       "diag_plot_width",
  #       "Preferred plot width (px, optional)",
  #       min = 0, max = 2000, value = 0, step = 25
  #     ),
  #     helpText("If a wide width is requested, 2-column layouts automatically switch to 1-column to prevent overlap."),
  # 
  #     tags$hr(),
  # 
  #     checkboxInput("diag_show_conclusion", "Show academic conclusion panel", value = TRUE),
  # 
  #     conditionalPanel(
  #       condition = "input.diag_layout_preset == 'preset_custom'",
  #       tags$h5("Custom positions (1 = first)"),
  #       helpText("Give each enabled plot a position. Ties are broken by name."),
  # 
  #       numericInput("pos_ts",     "Residuals over time position", value = 1, min = 1, step = 1),
  #       numericInput("pos_acf",    "Residual ACF position",        value = 2, min = 1, step = 1),
  #       numericInput("pos_hist",   "Histogram position",           value = 3, min = 1, step = 1),
  #       numericInput("pos_qq",     "Q-Q position",                 value = 4, min = 1, step = 1),
  #       numericInput("pos_fitted", "Residuals vs fitted position", value = 5, min = 1, step = 1),
  #       numericInput("pos_scale",  "Scale-location position",      value = 6, min = 1, step = 1),
  #       numericInput("pos_lb",     "Ljung–Box p-values position",  value = 7, min = 1, step = 1)
  #     )
  #   )
  # })
  # 
  # # =========================
  # # Main panel UI (dynamic width)
  # # =========================
  # output$manual_diag_mainpanel_ui <- renderUI({
  #   column(
  #     width = diag_main_w(),
  # 
  #     tags$div(
  #       style = "padding:10px;border:1px solid #e5e5e5;border-radius:10px;background:#fcfcfc;margin-bottom:12px;",
  #       tags$h4("Residual diagnostics (Manual SARIMA)"),
  #       tags$p("Choose plots, layout, and sizing from the left panel.")
  #     ),
  # 
  #     uiOutput("manual_diag_plots_ui"),
  # 
  #     conditionalPanel(
  #       condition = "input.diag_show_conclusion == true",
  #       tags$h4("Academic conclusion (Diagnostics)"),
  #       tags$div(
  #         style = "padding:10px;border:1px solid #e5e5e5;border-radius:10px;background:#ffffff;",
  #         verbatimTextOutput("manual_diag_commentary_diag")
  #       )
  #     )
  #   )
  # })
  # 
  # # =========================
  # # Plot outputs (diag)
  # # =========================
  # 
  # output$manual_resid_ts_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 3, "Not enough residuals to plot."))
  # 
  #   plot(res, type = "l", main = "Residuals over time", xlab = "Time", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  # })
  # 
  # output$manual_resid_acf_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for ACF."))
  # 
  #   plot(acf(res, plot = FALSE), main = "Residual ACF")
  # })
  # 
  # output$manual_resid_hist_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for histogram."))
  # 
  #   hist(res, breaks = 30, col = "gray85", border = "white",
  #        main = "Residual histogram", xlab = "Residual")
  # })
  # 
  # output$manual_resid_qq_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for Q-Q plot."))
  # 
  #   qqnorm(res, main = "Normal Q–Q")
  #   qqline(res, col = "red", lwd = 2)
  # })
  # 
  # output$manual_resid_fitted_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   validate(need(sum(ok) >= 10, "Not enough valid points for residuals vs fitted."))
  # 
  #   plot(fv[ok], res[ok],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("steelblue", alpha.f = 0.7),
  #        main = "Residuals vs fitted",
  #        xlab = "Fitted values", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  #   if (sum(ok) > 20) lines(stats::lowess(fv[ok], res[ok], f = 2/3), col = "tomato", lwd = 2)
  # })
  # 
  # output$manual_resid_scale_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  # 
  #   y <- sqrt(abs(res))
  #   ok2 <- ok & is.finite(y)
  #   validate(need(sum(ok2) >= 10, "Not enough valid points for scale-location plot."))
  # 
  #   plot(fv[ok2], y[ok2],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("darkgreen", alpha.f = 0.65),
  #        main = "Scale-location (sqrt(|res|) vs fitted)",
  #        xlab = "Fitted values", ylab = expression(sqrt("|Residual|")))
  #   if (sum(ok2) > 20) lines(stats::lowess(fv[ok2], y[ok2], f = 2/3), col = "tomato", lwd = 2)
  # })
  # 
  # output$manual_resid_lb_pvals_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   nice_par()
  #   res   <- residuals(fit)
  #   fitdf <- length(coef(fit))
  # 
  #   ctrl <- get_diag_controls(input)
  #   out  <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
  # 
  #   validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
  # 
  #   plot_lb_pvals(out$pvals, out$L, ctrl$alpha, fitdf, "Ljung–Box p-values by lag (Diagnostics)")
  # })
  # 
  # output$manual_diag_commentary_diag <- renderPrint({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   N   <- length(res)
  # 
  #   ctrl  <- get_diag_controls(input)
  #   fitdf <- length(coef(fit))
  # 
  #   if (N < 8) {
  #     cat(
  #       "Diagnostics summary (Manual SARIMA)\n",
  #       "------------------------------------------------------------\n",
  #       "- Too few residuals to run a stable diagnostic battery.\n",
  #       "- Fit the model on more observations or reduce differencing/orders.\n",
  #       sep = ""
  #     )
  #     return(invisible())
  #   }
  # 
  #   L2 <- max(1L, min(as.integer(ctrl$L), N - 1L))
  #   lb <- tryCatch(stats::Box.test(res, lag = L2, type = "Ljung-Box", fitdf = fitdf),
  #                  error = function(e) NULL)
  # 
  #   cat("Diagnostics summary (Manual SARIMA)\n")
  #   cat("------------------------------------------------------------\n")
  #   if (!is.null(lb) && is.finite(lb$p.value)) {
  #     cat("Ljung–Box p =", signif(lb$p.value, 3),
  #         if (lb$p.value > ctrl$alpha) " -> no strong evidence of autocorrelation.\n"
  #         else " -> remaining autocorrelation; revise orders/differencing.\n", sep = "")
  #   } else {
  #     cat("Ljung–Box unavailable.\n")
  #   }
  # })
  # 
  # # =========================
  # # Diagnostics plot layout builder (AUTO: 2-col -> 1-col if plots too wide)
  # # =========================
  # 
  # manual_diag_plot_map <- list(
  #   ts     = list(id = "manual_resid_ts_diag",       title = "Residuals over time"),
  #   acf    = list(id = "manual_resid_acf_diag",      title = "Residual ACF"),
  #   hist   = list(id = "manual_resid_hist_diag",     title = "Residual histogram"),
  #   qq     = list(id = "manual_resid_qq_diag",       title = "Normal Q-Q"),
  #   fitted = list(id = "manual_resid_fitted_diag",   title = "Residuals vs fitted"),
  #   scale  = list(id = "manual_resid_scale_diag",    title = "Scale-location"),
  #   lb     = list(id = "manual_resid_lb_pvals_diag", title = "Ljung–Box p-values by lag")
  # )
  # 
  # make_diag_plot_output <- function(plot_key, height_px, width_px = 0) {
  #   info <- manual_diag_plot_map[[plot_key]]
  #   if (is.null(info)) return(NULL)
  # 
  #   # Default responsive; if px requested, apply but constrain to container.
  #   w <- if (!is.finite(width_px) || width_px <= 0) "100%" else paste0(width_px, "px")
  # 
  #   tags$div(
  #     style = "margin-bottom:12px;",
  #     tags$h5(info$title),
  #     tags$div(
  #       style = "width:100%;max-width:100%;overflow:hidden;",
  #       plotOutput(info$id, height = height_px, width = w)
  #     )
  #   )
  # }
  # 
  # get_diag_plot_order <- reactive({
  #   sel <- input$diag_plots_selected
  #   if (is.null(sel) || length(sel) == 0) return(character(0))
  # 
  #   preset <- input$diag_layout_preset %||% "preset_current"
  #   if (preset != "preset_custom") {
  #     default_order <- c("ts", "acf", "hist", "qq", "fitted", "scale", "lb")
  #     return(default_order[default_order %in% sel])
  #   }
  # 
  #   pos <- c(
  #     ts     = input$pos_ts,
  #     acf    = input$pos_acf,
  #     hist   = input$pos_hist,
  #     qq     = input$pos_qq,
  #     fitted = input$pos_fitted,
  #     scale  = input$pos_scale,
  #     lb     = input$pos_lb
  #   )
  #   pos <- pos[names(pos) %in% sel]
  #   pos_num <- suppressWarnings(as.numeric(pos))
  #   pos_num[!is.finite(pos_num)] <- 999
  #   ord <- order(pos_num, names(pos_num))
  #   names(pos_num)[ord]
  # })
  # 
  # # Heuristic: when user requests a big plot width, switch to 1-column in 2-col presets.
  # # This avoids overlap and behaves like "column grows when plot grows".
  # use_one_col_due_to_width <- reactive({
  #   wpx <- input$diag_plot_width %||% 0
  #   wpx <- suppressWarnings(as.numeric(wpx))
  #   if (!is.finite(wpx) || wpx <= 0) return(FALSE)
  #   # threshold: if requested width is large, prefer 1 column
  #   wpx >= 700
  # })
  # 
  # output$manual_diag_plots_ui <- renderUI({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) first, then diagnostics will appear."))
  # 
  #   keys <- get_diag_plot_order()
  #   if (length(keys) == 0) {
  #     return(tags$div(
  #       style = "padding:10px;border:1px dashed #ccc;border-radius:10px;background:#fff;",
  #       tags$strong("No plots selected."),
  #       tags$p("Use the left panel to select plots.")
  #     ))
  #   }
  # 
  #   height_px <- input$diag_plot_height %||% 260
  #   width_px  <- input$diag_plot_width  %||% 0
  #   preset    <- input$diag_layout_preset %||% "preset_current"
  # 
  #   # Auto-switch: if plots are wide, force 1 column even in 2-col presets
  #   force_one_col <- isTRUE(use_one_col_due_to_width())
  # 
  #   if (preset == "preset_current" && !force_one_col) {
  #     row1 <- intersect(c("ts", "acf"), keys)
  #     row2 <- intersect(c("hist", "qq"), keys)
  #     rest <- setdiff(keys, c("ts", "acf", "hist", "qq"))
  # 
  #     ui <- tagList()
  #     if (length(row1) > 0) ui <- tagList(ui, fluidRow(lapply(row1, function(k) column(6, make_diag_plot_output(k, height_px, width_px)))))
  #     if (length(row2) > 0) ui <- tagList(ui, fluidRow(lapply(row2, function(k) column(6, make_diag_plot_output(k, height_px, width_px)))))
  #     if (length(rest) > 0) ui <- tagList(ui, lapply(rest, function(k) make_diag_plot_output(k, height_px, width_px)))
  #     return(ui)
  #   }
  # 
  #   if (preset == "preset_two_col" && !force_one_col) {
  #     left  <- keys[seq(1, length(keys), by = 2)]
  #     right <- keys[seq(2, length(keys), by = 2)]
  #     return(
  #       fluidRow(
  #         column(6, lapply(left,  function(k) make_diag_plot_output(k, height_px, width_px))),
  #         column(6, lapply(right, function(k) make_diag_plot_output(k, height_px, width_px)))
  #       )
  #     )
  #   }
  # 
  #   # preset_one_col, preset_custom, OR forced one-col due to wide plots:
  #   tagList(lapply(keys, function(k) make_diag_plot_output(k, height_px, width_px)))
  # })
  
  
  
  
  # # =========================
  # # Manual SARIMA Diagnostics (FULL) - responsive panel + builder
  # #   MOD: allow real enlargement via scrollable plots panel (px width)
  # # =========================
  # 
  # `%||%` <- function(a, b) if (!is.null(a)) a else b
  # 
  # nice_par <- function() {
  #   par(mar = c(4, 4, 2.2, 1), mgp = c(2.2, 0.7, 0), las = 1)
  # }
  # 
  # # ---- Shared helpers ----
  # get_diag_controls <- function(input) {
  #   L_input <- suppressWarnings(as.integer(input$diag_lag))
  #   L <- if (is.finite(L_input) && L_input > 0) L_input else 12L
  #   
  #   alpha <- suppressWarnings(as.numeric(input$alphaSt2))
  #   if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05
  #   
  #   list(L = L, alpha = alpha)
  # }
  # 
  # compute_lb_pvals <- function(res, L, fitdf) {
  #   res <- as.numeric(res)
  #   res <- res[is.finite(res)]
  #   
  #   N <- length(res)
  #   if (N < 8) return(list(pvals = numeric(0), L = 0L, N = N))
  #   
  #   L <- as.integer(L)
  #   L <- max(1L, min(L, N - 1L))
  #   
  #   pvals <- rep(NA_real_, L)
  #   for (k in seq_len(L)) {
  #     if (k > fitdf) {
  #       bt <- tryCatch(
  #         stats::Box.test(res, lag = k, type = "Ljung-Box", fitdf = fitdf),
  #         error = function(e) NULL
  #       )
  #       pvals[k] <- if (!is.null(bt)) as.numeric(bt$p.value) else NA_real_
  #     }
  #   }
  #   list(pvals = pvals, L = L, N = N)
  # }
  # 
  # plot_lb_pvals <- function(pvals, L, alpha, fitdf, main_title) {
  #   plot(seq_len(L), pvals,
  #        type = "h", lwd = 2,
  #        xlab = "Lag (k)",
  #        ylab = "p-value  (Ljung–Box up to lag k)",
  #        main = main_title,
  #        ylim = c(0, 1))
  #   points(seq_len(L), pvals, pch = 16)
  #   abline(h = alpha, lty = 2, col = "gray40")
  #   mtext(sprintf("alpha = %.3f, fitdf = %d", alpha, fitdf),
  #         side = 3, line = 0.2, cex = 0.85)
  #   
  #   if (is.finite(fitdf) && fitdf >= 1 && fitdf < L) {
  #     rect(xleft = 0.5, ybottom = -0.02, xright = fitdf + 0.5, ytop = 1.02,
  #          border = NA, col = grDevices::adjustcolor("gray", alpha.f = 0.15))
  #     text(x = (fitdf + 1) / 2, y = 0.95, labels = "df ≤ 0 (omitted)", cex = 0.8)
  #   }
  # }
  # 
  # # =========================
  # # IMPORTANT MOD:
  # # instead of "Bootstrap columns" resizing (which cannot exceed the screen),
  # # we allow a PX-width plots panel with horizontal scrolling.
  # # We'll still keep main+sidebar stable (sum=12) so it stays on the right.
  # # =========================
  # 
  # # main width from slider (keeps main on right; do NOT try to exceed screen with this)
  # diag_main_w <- reactive({
  #   w <- input$diag_main_width %||% 9
  #   w <- as.integer(w)
  #   if (!is.finite(w)) w <- 9
  #   w <- max(6L, min(11L, w))
  #   w
  # })
  # 
  # # sidebar width auto-adjust so main stays right (sum = 12)
  # diag_sidebar_w <- reactive({
  #   12L - diag_main_w()
  # })
  # 
  # # =========================
  # # Sidebar UI (dynamic width)
  # # NOTE: you must add diag_panel_px slider in your UI to control the plot panel width (px).
  # # =========================
  # output$manual_diag_sidebar_ui <- renderUI({
  #   column(
  #     width = diag_sidebar_w(),
  #     
  #     tags$h4("Diagnostics layout builder"),
  #     
  #     selectInput(
  #       "diag_layout_preset",
  #       "Layout preset",
  #       choices = c(
  #         "Current layout (2x2 + Ljung-Box + Conclusion)" = "preset_current",
  #         "Two columns (all plots, stacked)"             = "preset_two_col",
  #         "Single column (all plots, stacked)"           = "preset_one_col",
  #         "Custom order (choose positions below)"        = "preset_custom"
  #       ),
  #       selected = "preset_current"
  #     ),
  #     
  #     tags$hr(),
  #     
  #     checkboxGroupInput(
  #       "diag_plots_selected",
  #       "Plots to display",
  #       choices = c(
  #         "Residuals over time"                          = "ts",
  #         "Residual ACF"                                 = "acf",
  #         "Residual histogram"                           = "hist",
  #         "Normal Q-Q"                                   = "qq",
  #         "Residuals vs fitted"                          = "fitted",
  #         "Scale-location (sqrt(|res|) vs fitted)"       = "scale",
  #         "Ljung–Box p-values by lag"                    = "lb"
  #       ),
  #       selected = c("ts", "acf", "hist", "qq", "lb", "fitted", "scale")
  #     ),
  #     
  #     tags$hr(),
  #     
  #     sliderInput(
  #       "diag_main_width",
  #       "Main panel share (Bootstrap columns)",
  #       min = 6, max = 11, value = 9, step = 1
  #     ),
  #     helpText("This controls how much of the page the Diagnostics area uses. Real enlargement uses the px panel width below."),
  #     
  #     tags$hr(),
  #     
  #     # NEW (must exist in UI): plots panel width in pixels (can exceed viewport; scrolls)
  #     sliderInput(
  #       "diag_panel_px",
  #       "Plots panel width (px) — enables enlargement via scroll",
  #       min = 600, max = 4000, value = 1200, step = 50
  #     ),
  #     helpText("Increase this to make 2-column plots larger than your screen; a horizontal scrollbar will appear."),
  #     
  #     tags$hr(),
  #     
  #     sliderInput("diag_plot_height", "Plot height (px)", min = 180, max = 900, value = 260, step = 10),
  #     sliderInput("diag_plot_width",  "Plot width (px, optional)", min = 0, max = 3000, value = 0, step = 25),
  #     helpText("width=0 uses full column width inside the plots panel. Non-zero width uses a fixed px width."),
  #     
  #     tags$hr(),
  #     
  #     checkboxInput("diag_show_conclusion", "Show academic conclusion panel", value = TRUE),
  #     
  #     conditionalPanel(
  #       condition = "input.diag_layout_preset == 'preset_custom'",
  #       tags$h5("Custom positions (1 = first)"),
  #       numericInput("pos_ts",     "Residuals over time position", value = 1, min = 1),
  #       numericInput("pos_acf",    "Residual ACF position",        value = 2, min = 1),
  #       numericInput("pos_hist",   "Histogram position",           value = 3, min = 1),
  #       numericInput("pos_qq",     "Q-Q position",                 value = 4, min = 1),
  #       numericInput("pos_fitted", "Residuals vs fitted position", value = 5, min = 1),
  #       numericInput("pos_scale",  "Scale-location position",      value = 6, min = 1),
  #       numericInput("pos_lb",     "Ljung–Box p-values position",  value = 7, min = 1)
  #     )
  #   )
  # })
  # 
  # # =========================
  # # Main panel UI (dynamic width)
  # # MOD: wrap plots UI in a scroll container whose inner div has min-width = diag_panel_px
  # # =========================
  # output$manual_diag_mainpanel_ui <- renderUI({
  #   column(
  #     width = diag_main_w(),
  #     
  #     tags$div(
  #       style = "padding:10px;border:1px solid #e5e5e5;border-radius:10px;background:#fcfcfc;margin-bottom:12px;",
  #       tags$h4("Residual diagnostics (Manual SARIMA)"),
  #       tags$p("Increase 'Plots panel width (px)' to enlarge plots beyond the screen; scroll horizontally to view.")
  #     ),
  #     
  #     # Scroll wrapper + fixed/min px width panel
  #     tags$div(
  #       style = "width:100%; overflow-x:auto; padding-bottom:6px;",
  #       tags$div(
  #         style = paste0(
  #           "min-width:", (input$diag_panel_px %||% 1200), "px;",
  #           "max-width:none;",
  #           "background:#fff;",
  #           "border:1px solid #eee;",
  #           "border-radius:10px;",
  #           "padding:10px;"
  #         ),
  #         uiOutput("manual_diag_plots_ui")
  #       )
  #     ),
  #     
  #     conditionalPanel(
  #       condition = "input.diag_show_conclusion == true",
  #       tags$h4("Academic conclusion (Diagnostics)"),
  #       tags$div(
  #         style = "padding:10px;border:1px solid #e5e5e5;border-radius:10px;background:#ffffff;",
  #         verbatimTextOutput("manual_diag_commentary_diag")
  #       )
  #     )
  #   )
  # })
  # 
  # # =========================
  # # Plot outputs (diag)
  # # =========================
  # 
  # output$manual_resid_ts_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 3, "Not enough residuals to plot."))
  # 
  #   plot(res, type = "l", main = "Residuals over time", xlab = "Time", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  # })
  # 
  # output$manual_resid_acf_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for ACF."))
  # 
  #   plot(acf(res, plot = FALSE), main = "Residual ACF")
  # })
  # 
  # output$manual_resid_hist_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for histogram."))
  # 
  #   hist(res, breaks = 30, col = "gray85", border = "white",
  #        main = "Residual histogram", xlab = "Residual")
  # })
  # 
  # output$manual_resid_qq_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for Q-Q plot."))
  # 
  #   qqnorm(res, main = "Normal Q–Q")
  #   qqline(res, col = "red", lwd = 2)
  # })
  # 
  # output$manual_resid_fitted_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   validate(need(sum(ok) >= 10, "Not enough valid points for residuals vs fitted."))
  # 
  #   plot(fv[ok], res[ok],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("steelblue", alpha.f = 0.7),
  #        main = "Residuals vs fitted",
  #        xlab = "Fitted values", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  #   if (sum(ok) > 20) lines(stats::lowess(fv[ok], res[ok], f = 2/3), col = "tomato", lwd = 2)
  # })
  # 
  # output$manual_resid_scale_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  # 
  #   y <- sqrt(abs(res))
  #   ok2 <- ok & is.finite(y)
  #   validate(need(sum(ok2) >= 10, "Not enough valid points for scale-location plot."))
  # 
  #   plot(fv[ok2], y[ok2],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("darkgreen", alpha.f = 0.65),
  #        main = "Scale-location (sqrt(|res|) vs fitted)",
  #        xlab = "Fitted values", ylab = expression(sqrt("|Residual|")))
  #   if (sum(ok2) > 20) lines(stats::lowess(fv[ok2], y[ok2], f = 2/3), col = "tomato", lwd = 2)
  # })
  # 
  # output$manual_resid_lb_pvals_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   nice_par()
  #   res   <- residuals(fit)
  #   fitdf <- length(coef(fit))
  # 
  #   ctrl <- get_diag_controls(input)
  #   out  <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
  # 
  #   validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
  # 
  #   plot_lb_pvals(out$pvals, out$L, ctrl$alpha, fitdf, "Ljung–Box p-values by lag (Diagnostics)")
  # })
  # 
  # output$manual_diag_commentary_diag <- renderPrint({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  # 
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   N   <- length(res)
  # 
  #   ctrl  <- get_diag_controls(input)
  #   fitdf <- length(coef(fit))
  # 
  #   if (N < 8) {
  #     cat(
  #       "Diagnostics summary (Manual SARIMA)\n",
  #       "------------------------------------------------------------\n",
  #       "- Too few residuals to run a stable diagnostic battery.\n",
  #       "- Fit the model on more observations or reduce differencing/orders.\n",
  #       sep = ""
  #     )
  #     return(invisible())
  #   }
  # 
  #   L2 <- max(1L, min(as.integer(ctrl$L), N - 1L))
  #   lb <- tryCatch(stats::Box.test(res, lag = L2, type = "Ljung-Box", fitdf = fitdf),
  #                  error = function(e) NULL)
  # 
  #   cat("Diagnostics summary (Manual SARIMA)\n")
  #   cat("------------------------------------------------------------\n")
  #   if (!is.null(lb) && is.finite(lb$p.value)) {
  #     cat("Ljung–Box p =", signif(lb$p.value, 3),
  #         if (lb$p.value > ctrl$alpha) " -> no strong evidence of autocorrelation.\n"
  #         else " -> remaining autocorrelation; revise orders/differencing.\n", sep = "")
  #   } else {
  #     cat("Ljung–Box unavailable.\n")
  #   }
  # })
  # 
  # # =========================
  # # Diagnostics plot layout builder
  # # MOD: remove forced 1-col switch; 2-col can remain large because panel can expand (scroll)
  # # =========================
  # 
  # manual_diag_plot_map <- list(
  #   ts     = list(id = "manual_resid_ts_diag",       title = "Residuals over time"),
  #   acf    = list(id = "manual_resid_acf_diag",      title = "Residual ACF"),
  #   hist   = list(id = "manual_resid_hist_diag",     title = "Residual histogram"),
  #   qq     = list(id = "manual_resid_qq_diag",       title = "Normal Q-Q"),
  #   fitted = list(id = "manual_resid_fitted_diag",   title = "Residuals vs fitted"),
  #   scale  = list(id = "manual_resid_scale_diag",    title = "Scale-location"),
  #   lb     = list(id = "manual_resid_lb_pvals_diag", title = "Ljung–Box p-values by lag")
  # )
  # 
  # make_diag_plot_output <- function(plot_key, height_px, width_px = 0) {
  #   info <- manual_diag_plot_map[[plot_key]]
  #   if (is.null(info)) return(NULL)
  #   
  #   # NOTE: no max-width clipping; the panel scrolls if larger than viewport
  #   w <- if (!is.finite(width_px) || width_px <= 0) "100%" else paste0(width_px, "px")
  #   
  #   tags$div(
  #     style = "margin-bottom:12px;",
  #     tags$h5(info$title),
  #     plotOutput(info$id, height = height_px, width = w)
  #   )
  # }
  # 
  # get_diag_plot_order <- reactive({
  #   sel <- input$diag_plots_selected
  #   if (is.null(sel) || length(sel) == 0) return(character(0))
  #   
  #   preset <- input$diag_layout_preset %||% "preset_current"
  #   if (preset != "preset_custom") {
  #     default_order <- c("ts", "acf", "hist", "qq", "fitted", "scale", "lb")
  #     return(default_order[default_order %in% sel])
  #   }
  #   
  #   pos <- c(
  #     ts     = input$pos_ts,
  #     acf    = input$pos_acf,
  #     hist   = input$pos_hist,
  #     qq     = input$pos_qq,
  #     fitted = input$pos_fitted,
  #     scale  = input$pos_scale,
  #     lb     = input$pos_lb
  #   )
  #   pos <- pos[names(pos) %in% sel]
  #   pos_num <- suppressWarnings(as.numeric(pos))
  #   pos_num[!is.finite(pos_num)] <- 999
  #   ord <- order(pos_num, names(pos_num))
  #   names(pos_num)[ord]
  # })
  # 
  # output$manual_diag_plots_ui <- renderUI({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) first, then diagnostics will appear."))
  #   
  #   keys <- get_diag_plot_order()
  #   if (length(keys) == 0) {
  #     return(tags$div(
  #       style = "padding:10px;border:1px dashed #ccc;border-radius:10px;background:#fff;",
  #       tags$strong("No plots selected."),
  #       tags$p("Use the left panel to select plots.")
  #     ))
  #   }
  #   
  #   height_px <- input$diag_plot_height %||% 260
  #   width_px  <- input$diag_plot_width  %||% 0
  #   preset    <- input$diag_layout_preset %||% "preset_current"
  #   
  #   if (preset == "preset_current") {
  #     row1 <- intersect(c("ts", "acf"), keys)
  #     row2 <- intersect(c("hist", "qq"), keys)
  #     rest <- setdiff(keys, c("ts", "acf", "hist", "qq"))
  #     
  #     ui <- tagList()
  #     if (length(row1) > 0) ui <- tagList(ui, fluidRow(lapply(row1, function(k) column(6, make_diag_plot_output(k, height_px, width_px)))))
  #     if (length(row2) > 0) ui <- tagList(ui, fluidRow(lapply(row2, function(k) column(6, make_diag_plot_output(k, height_px, width_px)))))
  #     if (length(rest) > 0) ui <- tagList(ui, lapply(rest, function(k) make_diag_plot_output(k, height_px, width_px)))
  #     return(ui)
  #   }
  #   
  #   if (preset == "preset_two_col") {
  #     left  <- keys[seq(1, length(keys), by = 2)]
  #     right <- keys[seq(2, length(keys), by = 2)]
  #     return(
  #       fluidRow(
  #         column(6, lapply(left,  function(k) make_diag_plot_output(k, height_px, width_px))),
  #         column(6, lapply(right, function(k) make_diag_plot_output(k, height_px, width_px)))
  #       )
  #     )
  #   }
  #   
  #   # one column
  #   tagList(lapply(keys, function(k) make_diag_plot_output(k, height_px, width_px)))
  # })
  
  
  
 
  
  # # ============================================================
  # # Manual SARIMA Diagnostics (FULL) — scrollable canvas + builder
  # # + NEW Ljung–Box p-values by lag plot for Diagnostics tab ONLY
  # #
  # # IMPORTANT:
  # # - Keep your Academic conclusion plot output name (manual_resid_lb_pvals) unchanged.
  # # - Diagnostics tab uses a NEW output id: manual_resid_lb_pvals_diag2
  # # - UI must call uiOutput("manual_diag_canvas_ui") and uiOutput("manual_diag_plots_ui")
  # #   as in the flex/canvas layout we discussed.
  # # ============================================================
  # 
  # `%||%` <- function(a, b) if (!is.null(a)) a else b
  # 
  # nice_par <- function() {
  #   par(mar = c(4, 4, 2.2, 1), mgp = c(2.2, 0.7, 0), las = 1)
  # }
  # 
  # # ---- Shared helpers ----
  # get_diag_controls <- function(input) {
  #   L_input <- suppressWarnings(as.integer(input$diag_lag))
  #   L <- if (is.finite(L_input) && L_input > 0) L_input else 12L
  #   
  #   alpha <- suppressWarnings(as.numeric(input$alphaSt2))
  #   if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05
  #   
  #   list(L = L, alpha = alpha)
  # }
  # 
  # compute_lb_pvals <- function(res, L, fitdf) {
  #   res <- as.numeric(res)
  #   res <- res[is.finite(res)]
  #   
  #   N <- length(res)
  #   if (N < 8) return(list(pvals = numeric(0), L = 0L, N = N))
  #   
  #   L <- as.integer(L)
  #   L <- max(1L, min(L, N - 1L))
  #   
  #   pvals <- rep(NA_real_, L)
  #   for (k in seq_len(L)) {
  #     if (k > fitdf) {
  #       bt <- tryCatch(
  #         stats::Box.test(res, lag = k, type = "Ljung-Box", fitdf = fitdf),
  #         error = function(e) NULL
  #       )
  #       pvals[k] <- if (!is.null(bt)) as.numeric(bt$p.value) else NA_real_
  #     }
  #   }
  #   list(pvals = pvals, L = L, N = N)
  # }
  # 
  # plot_lb_pvals <- function(pvals, L, alpha, fitdf, main_title) {
  #   plot(seq_len(L), pvals,
  #        type = "h", lwd = 2,
  #        xlab = "Lag (k)",
  #        ylab = "p-value (Ljung–Box up to lag k)",
  #        main = main_title,
  #        ylim = c(0, 1))
  #   points(seq_len(L), pvals, pch = 16)
  #   abline(h = alpha, lty = 2, col = "gray40")
  #   mtext(sprintf("alpha = %.3f, fitdf = %d", alpha, fitdf),
  #         side = 3, line = 0.2, cex = 0.85)
  #   
  #   if (is.finite(fitdf) && fitdf >= 1 && fitdf < L) {
  #     rect(xleft = 0.5, ybottom = -0.02, xright = fitdf + 0.5, ytop = 1.02,
  #          border = NA, col = grDevices::adjustcolor("gray", alpha.f = 0.15))
  #     text(x = (fitdf + 1) / 2, y = 0.95, labels = "df ≤ 0 (omitted)", cex = 0.8)
  #   }
  # }
  # 
  # # ============================================================
  # # Plot outputs (Diagnostics tab versions)
  # # ============================================================
  # 
  # output$manual_resid_ts_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 3, "Not enough residuals to plot."))
  #   
  #   plot(res, type = "l", main = "Residuals over time", xlab = "Time", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  # })
  # 
  # output$manual_resid_acf_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for ACF."))
  #   
  #   plot(acf(res, plot = FALSE), main = "Residual ACF")
  # })
  # 
  # output$manual_resid_hist_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for histogram."))
  #   
  #   hist(res, breaks = 30, col = "gray85", border = "white",
  #        main = "Residual histogram", xlab = "Residual")
  # })
  # 
  # output$manual_resid_qq_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for Q–Q plot."))
  #   
  #   qqnorm(res, main = "Normal Q–Q")
  #   qqline(res, col = "red", lwd = 2)
  # })
  # 
  # output$manual_resid_fitted_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   validate(need(sum(ok) >= 10, "Not enough valid points for residuals vs fitted."))
  #   
  #   plot(fv[ok], res[ok],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("steelblue", alpha.f = 0.7),
  #        main = "Residuals vs fitted",
  #        xlab = "Fitted values", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  #   if (sum(ok) > 20) lines(stats::lowess(fv[ok], res[ok], f = 2/3), col = "tomato", lwd = 2)
  # })
  # 
  # output$manual_resid_scale_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   
  #   y <- sqrt(abs(res))
  #   ok2 <- ok & is.finite(y)
  #   validate(need(sum(ok2) >= 10, "Not enough valid points for scale-location plot."))
  #   
  #   plot(fv[ok2], y[ok2],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("darkgreen", alpha.f = 0.65),
  #        main = "Scale-location (sqrt(|res|) vs fitted)",
  #        xlab = "Fitted values", ylab = expression(sqrt("|Residual|")))
  #   if (sum(ok2) > 20) lines(stats::lowess(fv[ok2], y[ok2], f = 2/3), col = "tomato", lwd = 2)
  # })
  # 
  # # ------------------------------------------------------------
  # # NEW: Ljung–Box p-values by lag (Diagnostics tab ONLY)
  # # This is the one the Diagnostics UI should use.
  # # ------------------------------------------------------------
  # output$manual_resid_lb_pvals_diag2 <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res   <- residuals(fit)
  #   fitdf <- length(coef(fit))
  #   
  #   ctrl <- get_diag_controls(input)
  #   out  <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
  #   
  #   validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
  #   
  #   plot_lb_pvals(
  #     pvals = out$pvals,
  #     L = out$L,
  #     alpha = ctrl$alpha,
  #     fitdf = fitdf,
  #     main_title = "Ljung–Box p-values by lag (Diagnostics tab)"
  #   )
  # })
  # 
  # # ------------------------------------------------------------
  # # Commentary (Diagnostics tab)
  # # ------------------------------------------------------------
  # output$manual_diag_commentary_diag <- renderPrint({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   N   <- length(res)
  #   
  #   ctrl  <- get_diag_controls(input)
  #   fitdf <- length(coef(fit))
  #   
  #   if (N < 8) {
  #     cat(
  #       "Diagnostics summary (Manual SARIMA)\n",
  #       "------------------------------------------------------------\n",
  #       "- Too few residuals to run a stable diagnostic battery.\n",
  #       "- Fit the model on more observations or reduce differencing/orders.\n",
  #       sep = ""
  #     )
  #     return(invisible())
  #   }
  #   
  #   L2 <- max(1L, min(as.integer(ctrl$L), N - 1L))
  #   lb <- tryCatch(
  #     stats::Box.test(res, lag = L2, type = "Ljung-Box", fitdf = fitdf),
  #     error = function(e) NULL
  #   )
  #   
  #   cat("Diagnostics summary (Manual SARIMA)\n")
  #   cat("------------------------------------------------------------\n")
  #   
  #   if (!is.null(lb) && is.finite(lb$p.value)) {
  #     cat("Ljung–Box (lag ", L2, ") p = ", signif(lb$p.value, 3),
  #         if (lb$p.value > ctrl$alpha) " -> no strong evidence of residual autocorrelation.\n"
  #         else " -> residual autocorrelation remains; revise (p,q)(P,Q) and/or (d,D,s).\n",
  #         sep = "")
  #   } else {
  #     cat("Ljung–Box test unavailable.\n")
  #   }
  #   
  #   cat("\nNotes:\n")
  #   cat("- Use Residual ACF + Ljung–Box together: spikes + small p-values suggest underfitting.\n")
  #   cat("- If residual distribution is heavy-tailed (Q–Q), inference/intervals may be optimistic.\n")
  # })
  # 
  # # ============================================================
  # # KEEP: Academic conclusion plot output name unchanged
  # # (this is your REPORT / Academic conclusion version)
  # # ============================================================
  # output$manual_resid_lb_pvals <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   nice_par()
  #   res   <- residuals(fit)
  #   fitdf <- length(coef(fit))
  #   
  #   ctrl <- get_diag_controls(input)
  #   out  <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
  #   
  #   validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
  #   
  #   plot_lb_pvals(
  #     pvals = out$pvals,
  #     L = out$L,
  #     alpha = ctrl$alpha,
  #     fitdf = fitdf,
  #     main_title = "Ljung–Box p-values by lag (Manual SARIMA)"
  #   )
  # })
  # 
  # # ============================================================
  # # Diagnostics canvas + layout builder (NO size restrictions)
  # # - You can scroll right/down.
  # # - Canvas size controlled by diag_canvas_width / diag_canvas_height (UI sliders).
  # # ============================================================
  # 
  # manual_diag_plot_map <- list(
  #   ts     = list(id = "manual_resid_ts_diag",        title = "Residuals over time"),
  #   acf    = list(id = "manual_resid_acf_diag",       title = "Residual ACF"),
  #   hist   = list(id = "manual_resid_hist_diag",      title = "Residual histogram"),
  #   qq     = list(id = "manual_resid_qq_diag",        title = "Normal Q–Q"),
  #   fitted = list(id = "manual_resid_fitted_diag",    title = "Residuals vs fitted"),
  #   scale  = list(id = "manual_resid_scale_diag",     title = "Scale-location"),
  #   # NEW diagnostics-only LB plot:
  #   lb     = list(id = "manual_resid_lb_pvals_diag2", title = "Ljung–Box p-values by lag (Diagnostics)")
  # )
  # 
  # make_diag_plot_output <- function(plot_key, height_px, width_px = 0) {
  #   info <- manual_diag_plot_map[[plot_key]]
  #   if (is.null(info)) return(NULL)
  #   
  #   w <- if (!is.finite(width_px) || width_px <= 0) "100%" else paste0(width_px, "px")
  #   
  #   tags$div(
  #     style = "margin-bottom:14px; padding:10px; border:1px solid #eee; border-radius:10px; background:#fff;",
  #     tags$h5(info$title),
  #     plotOutput(info$id, height = height_px, width = w)
  #   )
  # }
  # 
  # get_diag_plot_order <- reactive({
  #   sel <- input$diag_plots_selected
  #   if (is.null(sel) || length(sel) == 0) return(character(0))
  #   
  #   preset <- input$diag_layout_preset %||% "preset_current"
  #   if (preset != "preset_custom") {
  #     default_order <- c("ts", "acf", "hist", "qq", "fitted", "scale", "lb")
  #     return(default_order[default_order %in% sel])
  #   }
  #   
  #   pos <- c(
  #     ts     = input$pos_ts,
  #     acf    = input$pos_acf,
  #     hist   = input$pos_hist,
  #     qq     = input$pos_qq,
  #     fitted = input$pos_fitted,
  #     scale  = input$pos_scale,
  #     lb     = input$pos_lb
  #   )
  #   pos <- pos[names(pos) %in% sel]
  #   pos_num <- suppressWarnings(as.numeric(pos))
  #   pos_num[!is.finite(pos_num)] <- 999
  #   ord <- order(pos_num, names(pos_num))
  #   names(pos_num)[ord]
  # })
  # 
  # # Outer scroll container should be in UI (flex main area with overflow:auto).
  # # This inner canvas sets a large area to place plots into.
  # output$manual_diag_canvas_ui <- renderUI({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) first, then diagnostics will appear."))
  #   
  #   canvas_w <- input$diag_canvas_width  %||% 1600
  #   canvas_h <- input$diag_canvas_height %||% 1800
  #   canvas_w <- suppressWarnings(as.integer(canvas_w)); if (!is.finite(canvas_w)) canvas_w <- 1600L
  #   canvas_h <- suppressWarnings(as.integer(canvas_h)); if (!is.finite(canvas_h)) canvas_h <- 1800L
  #   
  #   tags$div(
  #     style = paste0(
  #       "width:", canvas_w, "px;",
  #       "min-height:", canvas_h, "px;",
  #       "background:#fcfcfc;"
  #     ),
  #     uiOutput("manual_diag_plots_ui")
  #   )
  # })
  # 
  # output$manual_diag_plots_ui <- renderUI({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) first, then diagnostics will appear."))
  #   
  #   keys <- get_diag_plot_order()
  #   if (length(keys) == 0) {
  #     return(tags$div(
  #       style = "padding:10px;border:1px dashed #ccc;border-radius:10px;background:#fff;",
  #       tags$strong("No plots selected."),
  #       tags$p("Use the left panel to select plots.")
  #     ))
  #   }
  #   
  #   height_px <- input$diag_plot_height %||% 320
  #   width_px  <- input$diag_plot_width  %||% 0
  #   preset    <- input$diag_layout_preset %||% "preset_current"
  #   
  #   if (preset == "preset_current") {
  #     row1 <- intersect(c("ts", "acf"), keys)
  #     row2 <- intersect(c("hist", "qq"), keys)
  #     rest <- setdiff(keys, c("ts", "acf", "hist", "qq"))
  #     
  #     ui <- tagList()
  #     if (length(row1) > 0) {
  #       ui <- tagList(ui, fluidRow(lapply(row1, function(k) column(6, make_diag_plot_output(k, height_px, width_px)))))
  #     }
  #     if (length(row2) > 0) {
  #       ui <- tagList(ui, fluidRow(lapply(row2, function(k) column(6, make_diag_plot_output(k, height_px, width_px)))))
  #     }
  #     if (length(rest) > 0) {
  #       ui <- tagList(ui, lapply(rest, function(k) make_diag_plot_output(k, height_px, width_px)))
  #     }
  #     return(ui)
  #   }
  #   
  #   if (preset == "preset_two_col") {
  #     left  <- keys[seq(1, length(keys), by = 2)]
  #     right <- keys[seq(2, length(keys), by = 2)]
  #     return(
  #       fluidRow(
  #         column(6, lapply(left,  function(k) make_diag_plot_output(k, height_px, width_px))),
  #         column(6, lapply(right, function(k) make_diag_plot_output(k, height_px, width_px)))
  #       )
  #     )
  #   }
  #   
  #   # preset_one_col + preset_custom (one column)
  #   tagList(lapply(keys, function(k) make_diag_plot_output(k, height_px, width_px)))
  # })
  
  
  
  
  
  
  
  # # ============================================================
  # # Manual SARIMA Diagnostics — ONLY 2 sliders (plot width/height)
  # # Auto canvas sizing based on layout + number of plots selected
  # # Includes NEW diagnostics-only Ljung–Box p-values plot
  # # Keeps Academic plot output name unchanged
  # # ============================================================
  # 
  # `%||%` <- function(a, b) if (!is.null(a)) a else b
  # 
  # nice_par <- function() {
  #   par(mar = c(4, 4, 2.2, 1), mgp = c(2.2, 0.7, 0), las = 1)
  # }
  # 
  # # ---- Helpers ----
  # get_diag_controls <- function(input) {
  #   L_input <- suppressWarnings(as.integer(input$diag_lag))
  #   L <- if (is.finite(L_input) && L_input > 0) L_input else 12L
  #   
  #   alpha <- suppressWarnings(as.numeric(input$alphaSt2))
  #   if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05
  #   
  #   list(L = L, alpha = alpha)
  # }
  # 
  # compute_lb_pvals <- function(res, L, fitdf) {
  #   res <- as.numeric(res)
  #   res <- res[is.finite(res)]
  #   N <- length(res)
  #   if (N < 8) return(list(pvals = numeric(0), L = 0L, N = N))
  #   
  #   L <- as.integer(L)
  #   L <- max(1L, min(L, N - 1L))
  #   
  #   pvals <- rep(NA_real_, L)
  #   for (k in seq_len(L)) {
  #     if (k > fitdf) {
  #       bt <- tryCatch(stats::Box.test(res, lag = k, type = "Ljung-Box", fitdf = fitdf),
  #                      error = function(e) NULL)
  #       pvals[k] <- if (!is.null(bt)) as.numeric(bt$p.value) else NA_real_
  #     }
  #   }
  #   list(pvals = pvals, L = L, N = N)
  # }
  # 
  # plot_lb_pvals <- function(pvals, L, alpha, fitdf, main_title) {
  #   plot(seq_len(L), pvals,
  #        type = "h", lwd = 2,
  #        xlab = "Lag (k)",
  #        ylab = "p-value (Ljung–Box up to lag k)",
  #        main = main_title,
  #        ylim = c(0, 1))
  #   points(seq_len(L), pvals, pch = 16)
  #   abline(h = alpha, lty = 2, col = "gray40")
  #   mtext(sprintf("alpha = %.3f, fitdf = %d", alpha, fitdf),
  #         side = 3, line = 0.2, cex = 0.85)
  #   
  #   if (is.finite(fitdf) && fitdf >= 1 && fitdf < L) {
  #     rect(xleft = 0.5, ybottom = -0.02, xright = fitdf + 0.5, ytop = 1.02,
  #          border = NA, col = grDevices::adjustcolor("gray", alpha.f = 0.15))
  #     text(x = (fitdf + 1) / 2, y = 0.95, labels = "df ≤ 0 (omitted)", cex = 0.8)
  #   }
  # }
  # 
  # # ============================================================
  # # Diagnostics plot outputs (independent)
  # # ============================================================
  # 
  # output$manual_resid_ts_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 3, "Not enough residuals to plot."))
  #   
  #   plot(res, type = "l", main = "Residuals over time", xlab = "Time", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  # })
  # 
  # output$manual_resid_acf_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for ACF."))
  #   
  #   plot(acf(res, plot = FALSE), main = "Residual ACF")
  # })
  # 
  # output$manual_resid_hist_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for histogram."))
  #   
  #   hist(res, breaks = 30, col = "gray85", border = "white",
  #        main = "Residual histogram", xlab = "Residual")
  # })
  # 
  # output$manual_resid_qq_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for Q–Q plot."))
  #   
  #   qqnorm(res, main = "Normal Q–Q")
  #   qqline(res, col = "red", lwd = 2)
  # })
  # 
  # output$manual_resid_fitted_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   validate(need(sum(ok) >= 10, "Not enough valid points for residuals vs fitted."))
  #   
  #   plot(fv[ok], res[ok],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("steelblue", alpha.f = 0.7),
  #        main = "Residuals vs fitted",
  #        xlab = "Fitted values", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  #   if (sum(ok) > 20) lines(stats::lowess(fv[ok], res[ok], f = 2/3), col = "tomato", lwd = 2)
  # })
  # 
  # output$manual_resid_scale_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   y   <- sqrt(abs(res))
  #   ok2 <- ok & is.finite(y)
  #   validate(need(sum(ok2) >= 10, "Not enough valid points for scale-location plot."))
  #   
  #   plot(fv[ok2], y[ok2],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("darkgreen", alpha.f = 0.65),
  #        main = "Scale-location (sqrt(|res|) vs fitted)",
  #        xlab = "Fitted values", ylab = expression(sqrt("|Residual|")))
  #   if (sum(ok2) > 20) lines(stats::lowess(fv[ok2], y[ok2], f = 2/3), col = "tomato", lwd = 2)
  # })
  # 
  # # NEW diagnostics-only Ljung–Box plot (do NOT reuse academic one)
  # output$manual_resid_lb_pvals_diag2 <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res   <- residuals(fit)
  #   fitdf <- length(coef(fit))
  #   ctrl  <- get_diag_controls(input)
  #   out   <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
  #   
  #   validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
  #   
  #   plot_lb_pvals(out$pvals, out$L, ctrl$alpha, fitdf,
  #                 main_title = "Ljung–Box p-values by lag (Diagnostics)")
  # })
  # 
  # # Commentary (kept simple; you can expand later)
  # output$manual_diag_commentary_diag <- renderPrint({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   N   <- length(res)
  #   
  #   ctrl  <- get_diag_controls(input)
  #   fitdf <- length(coef(fit))
  #   
  #   if (N < 8) {
  #     cat("Diagnostics summary (Manual SARIMA)\n",
  #         "------------------------------------------------------------\n",
  #         "- Too few residuals to run diagnostics.\n", sep = "")
  #     return(invisible())
  #   }
  #   
  #   L2 <- max(1L, min(as.integer(ctrl$L), N - 1L))
  #   lb <- tryCatch(stats::Box.test(res, lag = L2, type = "Ljung-Box", fitdf = fitdf),
  #                  error = function(e) NULL)
  #   
  #   cat("Diagnostics summary (Manual SARIMA)\n")
  #   cat("------------------------------------------------------------\n")
  #   if (!is.null(lb) && is.finite(lb$p.value)) {
  #     cat("Ljung–Box (lag ", L2, ") p = ", signif(lb$p.value, 3),
  #         if (lb$p.value > ctrl$alpha) " -> no strong evidence of autocorrelation.\n"
  #         else " -> remaining autocorrelation; revise orders/differencing.\n",
  #         sep = "")
  #   } else {
  #     cat("Ljung–Box unavailable.\n")
  #   }
  # })
  # 
  # # ============================================================
  # # KEEP academic conclusion plot output unchanged
  # # ============================================================
  # output$manual_resid_lb_pvals <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res   <- residuals(fit)
  #   fitdf <- length(coef(fit))
  #   ctrl  <- get_diag_controls(input)
  #   out   <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
  #   
  #   validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
  #   
  #   plot_lb_pvals(out$pvals, out$L, ctrl$alpha, fitdf,
  #                 main_title = "Ljung–Box p-values by lag (Manual SARIMA)")
  # })
  # 
  # # ============================================================
  # # Diagnostics layout builder + AUTO canvas sizing
  # # ============================================================
  # 
  # manual_diag_plot_map <- list(
  #   ts     = list(id = "manual_resid_ts_diag",        title = "Residuals over time"),
  #   acf    = list(id = "manual_resid_acf_diag",       title = "Residual ACF"),
  #   hist   = list(id = "manual_resid_hist_diag",      title = "Residual histogram"),
  #   qq     = list(id = "manual_resid_qq_diag",        title = "Normal Q–Q"),
  #   fitted = list(id = "manual_resid_fitted_diag",    title = "Residuals vs fitted"),
  #   scale  = list(id = "manual_resid_scale_diag",     title = "Scale-location"),
  #   lb     = list(id = "manual_resid_lb_pvals_diag2", title = "Ljung–Box p-values by lag (Diagnostics)")
  # )
  # 
  # get_diag_plot_order <- reactive({
  #   sel <- input$diag_plots_selected
  #   if (is.null(sel) || length(sel) == 0) return(character(0))
  #   
  #   preset <- input$diag_layout_preset %||% "preset_current"
  #   if (preset != "preset_custom") {
  #     default_order <- c("ts", "acf", "hist", "qq", "fitted", "scale", "lb")
  #     return(default_order[default_order %in% sel])
  #   }
  #   
  #   pos <- c(
  #     ts     = input$pos_ts,
  #     acf    = input$pos_acf,
  #     hist   = input$pos_hist,
  #     qq     = input$pos_qq,
  #     fitted = input$pos_fitted,
  #     scale  = input$pos_scale,
  #     lb     = input$pos_lb
  #   )
  #   pos <- pos[names(pos) %in% sel]
  #   pos_num <- suppressWarnings(as.numeric(pos))
  #   pos_num[!is.finite(pos_num)] <- 999
  #   ord <- order(pos_num, names(pos_num))
  #   names(pos_num)[ord]
  # })
  # 
  # # Canvas sizing from plot size + layout + number of plots
  # diag_canvas_dims <- reactive({
  #   keys   <- get_diag_plot_order()
  #   preset <- input$diag_layout_preset %||% "preset_current"
  #   
  #   pw <- suppressWarnings(as.integer(input$diag_plot_width))
  #   ph <- suppressWarnings(as.integer(input$diag_plot_height))
  #   if (!is.finite(pw) || pw < 300) pw <- 650L
  #   if (!is.finite(ph) || ph < 200) ph <- 320L
  #   
  #   # layout columns
  #   ncol <- if (preset %in% c("preset_one_col", "preset_custom")) 1L else 2L
  #   if (preset == "preset_current") ncol <- 2L
  #   
  #   n <- length(keys)
  #   nrow <- if (n == 0) 1L else ceiling(n / ncol)
  #   
  #   # add padding/headers per plot "card"
  #   card_pad_w <- 60L
  #   card_pad_h <- 90L
  #   
  #   canvas_w <- ncol * (pw + card_pad_w) + 40L
  #   canvas_h <- nrow * (ph + card_pad_h) + 40L
  #   
  #   list(w = canvas_w, h = canvas_h, pw = pw, ph = ph, ncol = ncol)
  # })
  # 
  # make_diag_plot_output <- function(plot_key, pw, ph) {
  #   info <- manual_diag_plot_map[[plot_key]]
  #   if (is.null(info)) return(NULL)
  #   
  #   tags$div(
  #     style = sprintf(
  #       "margin-bottom:14px; padding:10px; border:1px solid #eee; border-radius:10px; background:#fff; width:%dpx;",
  #       pw + 20
  #     ),
  #     tags$h5(info$title),
  #     plotOutput(info$id, height = ph, width = pw)
  #   )
  # }
  # 
  # # Canvas UI: width/height computed automatically (no extra sliders)
  # output$manual_diag_canvas_ui <- renderUI({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) first, then diagnostics will appear."))
  #   
  #   dims <- diag_canvas_dims()
  #   keys <- get_diag_plot_order()
  #   if (length(keys) == 0) {
  #     return(tags$div(
  #       style = "padding:10px;border:1px dashed #ccc;border-radius:10px;background:#fff;",
  #       tags$strong("No plots selected."),
  #       tags$p("Use the left panel to select plots.")
  #     ))
  #   }
  #   
  #   tags$div(
  #     style = sprintf(
  #       "width:%dpx; min-height:%dpx; background:#fcfcfc;",
  #       dims$w, dims$h
  #     ),
  #     uiOutput("manual_diag_plots_ui")
  #   )
  # })
  # 
  # # Plot grid builder uses the same computed plot sizes
  # output$manual_diag_plots_ui <- renderUI({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) first, then diagnostics will appear."))
  #   
  #   keys <- get_diag_plot_order()
  #   dims <- diag_canvas_dims()
  #   pw <- dims$pw
  #   ph <- dims$ph
  #   preset <- input$diag_layout_preset %||% "preset_current"
  #   
  #   if (preset %in% c("preset_one_col", "preset_custom")) {
  #     return(tagList(lapply(keys, function(k) make_diag_plot_output(k, pw, ph))))
  #   }
  #   
  #   # Two-column layouts
  #   if (preset == "preset_current" || preset == "preset_two_col") {
  #     left  <- keys[seq(1, length(keys), by = 2)]
  #     right <- keys[seq(2, length(keys), by = 2)]
  #     return(
  #       fluidRow(
  #         column(6, lapply(left,  function(k) make_diag_plot_output(k, pw, ph))),
  #         column(6, lapply(right, function(k) make_diag_plot_output(k, pw, ph)))
  #       )
  #     )
  #   }
  #   
  #   # fallback
  #   tagList(lapply(keys, function(k) make_diag_plot_output(k, pw, ph)))
  # })
  
  
  # ============================================================
  # Diagnostics TAB (FULL) — container width slider + plot size sliders
  # - Slider controls the flex container width (replaces width:100%)
  # - Only TWO plot sliders: diag_plot_width + diag_plot_height
  # - NEW Ljung–Box plot for Diagnostics only: manual_resid_lb_pvals_diag2
  # - Keeps academic conclusion Ljung–Box output: manual_resid_lb_pvals
  # ============================================================
  
  # `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # nice_par <- function() {
  #   par(mar = c(4, 4, 2.2, 1), mgp = c(2.2, 0.7, 0), las = 1)
  # }
  # 
  # # ---- Shared helpers ----
  # get_diag_controls <- function(input) {
  #   L_input <- suppressWarnings(as.integer(input$diag_lag))
  #   L <- if (is.finite(L_input) && L_input > 0) L_input else 12L
  #   
  #   alpha <- suppressWarnings(as.numeric(input$alphaSt2))
  #   if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05
  #   
  #   list(L = L, alpha = alpha)
  # }
  # 
  # compute_lb_pvals <- function(res, L, fitdf) {
  #   res <- as.numeric(res)
  #   res <- res[is.finite(res)]
  #   
  #   N <- length(res)
  #   if (N < 8) return(list(pvals = numeric(0), L = 0L, N = N))
  #   
  #   L <- as.integer(L)
  #   L <- max(1L, min(L, N - 1L))
  #   
  #   pvals <- rep(NA_real_, L)
  #   for (k in seq_len(L)) {
  #     if (k > fitdf) {
  #       bt <- tryCatch(
  #         stats::Box.test(res, lag = k, type = "Ljung-Box", fitdf = fitdf),
  #         error = function(e) NULL
  #       )
  #       pvals[k] <- if (!is.null(bt)) as.numeric(bt$p.value) else NA_real_
  #     }
  #   }
  #   list(pvals = pvals, L = L, N = N)
  # }
  # 
  # plot_lb_pvals <- function(pvals, L, alpha, fitdf, main_title) {
  #   plot(seq_len(L), pvals,
  #        type = "h", lwd = 2,
  #        xlab = "Lag (k)",
  #        ylab = "p-value (Ljung–Box up to lag k)",
  #        main = main_title,
  #        ylim = c(0, 1))
  #   points(seq_len(L), pvals, pch = 16)
  #   abline(h = alpha, lty = 2, col = "gray40")
  #   mtext(sprintf("alpha = %.3f, fitdf = %d", alpha, fitdf),
  #         side = 3, line = 0.2, cex = 0.85)
  #   
  #   if (is.finite(fitdf) && fitdf >= 1 && fitdf < L) {
  #     rect(xleft = 0.5, ybottom = -0.02, xright = fitdf + 0.5, ytop = 1.02,
  #          border = NA, col = grDevices::adjustcolor("gray", alpha.f = 0.15))
  #     text(x = (fitdf + 1) / 2, y = 0.95, labels = "df ≤ 0 (omitted)", cex = 0.8)
  #   }
  # }
  # 
  # # ============================================================
  # # Diagnostics plots (independent outputs)
  # # ============================================================
  # 
  # output$manual_resid_ts_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 3, "Not enough residuals to plot."))
  #   
  #   plot(res, type = "l",
  #        main = "Residuals over time",
  #        xlab = "Time", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  # })
  # 
  # output$manual_resid_acf_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for ACF."))
  #   
  #   plot(acf(res, plot = FALSE), main = "Residual ACF")
  # })
  # 
  # output$manual_resid_hist_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for histogram."))
  #   
  #   hist(res, breaks = 30, col = "gray85", border = "white",
  #        main = "Residual histogram", xlab = "Residual")
  # })
  # 
  # output$manual_resid_qq_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   validate(need(length(res) >= 5, "Not enough residuals for Q–Q plot."))
  #   
  #   qqnorm(res, main = "Normal Q–Q")
  #   qqline(res, col = "red", lwd = 2)
  # })
  # 
  # output$manual_resid_fitted_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   validate(need(sum(ok) >= 10, "Not enough valid points for residuals vs fitted."))
  #   
  #   plot(fv[ok], res[ok],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("steelblue", alpha.f = 0.7),
  #        main = "Residuals vs fitted",
  #        xlab = "Fitted values", ylab = "Residual")
  #   abline(h = 0, col = "gray40", lty = 2)
  #   if (sum(ok) > 20) {
  #     lines(stats::lowess(fv[ok], res[ok], f = 2/3), col = "tomato", lwd = 2)
  #   }
  # })
  # 
  # output$manual_resid_scale_diag <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res <- as.numeric(residuals(fit))
  #   fv  <- as.numeric(fitted(fit))
  #   ok  <- is.finite(res) & is.finite(fv)
  #   y   <- sqrt(abs(res))
  #   ok2 <- ok & is.finite(y)
  #   validate(need(sum(ok2) >= 10, "Not enough valid points for scale-location plot."))
  #   
  #   plot(fv[ok2], y[ok2],
  #        pch = 16, cex = 0.7,
  #        col = grDevices::adjustcolor("darkgreen", alpha.f = 0.65),
  #        main = "Scale-location (sqrt(|res|) vs fitted)",
  #        xlab = "Fitted values", ylab = expression(sqrt("|Residual|")))
  #   if (sum(ok2) > 20) {
  #     lines(stats::lowess(fv[ok2], y[ok2], f = 2/3), col = "tomato", lwd = 2)
  #   }
  # })
  # 
  # # ------------------------------------------------------------
  # # NEW: Ljung–Box p-values by lag (Diagnostics ONLY)
  # # ------------------------------------------------------------
  # output$manual_resid_lb_pvals_diag2 <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res   <- residuals(fit)
  #   fitdf <- length(coef(fit))
  #   
  #   ctrl <- get_diag_controls(input)
  #   out  <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
  #   
  #   validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
  #   
  #   plot_lb_pvals(out$pvals, out$L, ctrl$alpha, fitdf,
  #                 main_title = "Ljung–Box p-values by lag (Diagnostics)")
  # })
  # 
  # # ------------------------------------------------------------
  # # Commentary (Diagnostics tab)
  # # ------------------------------------------------------------
  # output$manual_diag_commentary_diag <- renderPrint({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   N   <- length(res)
  #   
  #   ctrl  <- get_diag_controls(input)
  #   fitdf <- length(coef(fit))
  #   
  #   if (N < 8) {
  #     cat("Diagnostics summary (Manual SARIMA)\n",
  #         "------------------------------------------------------------\n",
  #         "- Too few residuals to run diagnostics.\n", sep = "")
  #     return(invisible())
  #   }
  #   
  #   L2 <- max(1L, min(as.integer(ctrl$L), N - 1L))
  #   lb <- tryCatch(
  #     stats::Box.test(res, lag = L2, type = "Ljung-Box", fitdf = fitdf),
  #     error = function(e) NULL
  #   )
  #   
  #   cat("Diagnostics summary (Manual SARIMA)\n")
  #   cat("------------------------------------------------------------\n")
  #   if (!is.null(lb) && is.finite(lb$p.value)) {
  #     cat("Ljung–Box (lag ", L2, ") p = ", signif(lb$p.value, 3),
  #         if (lb$p.value > ctrl$alpha) " -> no strong evidence of autocorrelation.\n"
  #         else " -> remaining autocorrelation; revise orders/differencing.\n",
  #         sep = "")
  #   } else {
  #     cat("Ljung–Box unavailable.\n")
  #   }
  #   
  #   cat("\nNotes:\n")
  #   cat("- Use Residual ACF + Ljung–Box together: spikes + small p-values suggest underfitting.\n")
  #   cat("- If Q–Q shows heavy tails: intervals/inference may be optimistic.\n")
  # })
  # 
  # # ============================================================
  # # KEEP: Academic conclusion plot output unchanged
  # # ============================================================
  # output$manual_resid_lb_pvals <- renderPlot({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
  #   fit <- manual_fit()
  #   nice_par()
  #   
  #   res   <- residuals(fit)
  #   fitdf <- length(coef(fit))
  #   
  #   ctrl <- get_diag_controls(input)
  #   out  <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
  #   
  #   validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
  #   validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
  #   
  #   plot_lb_pvals(out$pvals, out$L, ctrl$alpha, fitdf,
  #                 main_title = "Ljung–Box p-values by lag (Manual SARIMA)")
  # })
  # 
  # # ============================================================
  # # Layout builder + auto canvas sizing from ONLY plot width/height
  # # ============================================================
  # 
  # manual_diag_plot_map <- list(
  #   ts     = list(id = "manual_resid_ts_diag",        title = "Residuals over time"),
  #   acf    = list(id = "manual_resid_acf_diag",       title = "Residual ACF"),
  #   hist   = list(id = "manual_resid_hist_diag",      title = "Residual histogram"),
  #   qq     = list(id = "manual_resid_qq_diag",        title = "Normal Q–Q"),
  #   fitted = list(id = "manual_resid_fitted_diag",    title = "Residuals vs fitted"),
  #   scale  = list(id = "manual_resid_scale_diag",     title = "Scale-location"),
  #   lb     = list(id = "manual_resid_lb_pvals_diag2", title = "Ljung–Box p-values by lag (Diagnostics)")
  # )
  # 
  # get_diag_plot_order <- reactive({
  #   sel <- input$diag_plots_selected
  #   if (is.null(sel) || length(sel) == 0) return(character(0))
  #   
  #   preset <- input$diag_layout_preset %||% "preset_current"
  #   if (preset != "preset_custom") {
  #     default_order <- c("ts", "acf", "hist", "qq", "fitted", "scale", "lb")
  #     return(default_order[default_order %in% sel])
  #   }
  #   
  #   pos <- c(
  #     ts     = input$pos_ts,
  #     acf    = input$pos_acf,
  #     hist   = input$pos_hist,
  #     qq     = input$pos_qq,
  #     fitted = input$pos_fitted,
  #     scale  = input$pos_scale,
  #     lb     = input$pos_lb
  #   )
  #   pos <- pos[names(pos) %in% sel]
  #   pos_num <- suppressWarnings(as.numeric(pos))
  #   pos_num[!is.finite(pos_num)] <- 999
  #   ord <- order(pos_num, names(pos_num))
  #   names(pos_num)[ord]
  # })
  # 
  # # Canvas dimensions computed from plot size + how many plots + layout
  # diag_canvas_dims <- reactive({
  #   keys   <- get_diag_plot_order()
  #   preset <- input$diag_layout_preset %||% "preset_current"
  #   
  #   pw <- suppressWarnings(as.integer(input$diag_plot_width))
  #   ph <- suppressWarnings(as.integer(input$diag_plot_height))
  #   if (!is.finite(pw) || pw < 300) pw <- 650L
  #   if (!is.finite(ph) || ph < 200) ph <- 320L
  #   
  #   # Use 2 columns for current/two_col; 1 column for one_col/custom
  #   ncol <- if (preset %in% c("preset_one_col", "preset_custom")) 1L else 2L
  #   
  #   n <- length(keys)
  #   nrow <- if (n == 0) 1L else ceiling(n / ncol)
  #   
  #   # Card padding and heading space
  #   pad_w <- 60L
  #   pad_h <- 110L
  #   
  #   canvas_w <- ncol * (pw + pad_w) + 40L
  #   canvas_h <- nrow * (ph + pad_h) + 40L
  #   
  #   list(w = canvas_w, h = canvas_h, pw = pw, ph = ph)
  # })
  # 
  # make_diag_plot_output <- function(plot_key, pw, ph) {
  #   info <- manual_diag_plot_map[[plot_key]]
  #   if (is.null(info)) return(NULL)
  #   
  #   tags$div(
  #     style = sprintf(
  #       "margin-bottom:14px; padding:10px; border:1px solid #eee; border-radius:10px; background:#fff; width:%dpx;",
  #       pw + 30
  #     ),
  #     tags$h5(info$title),
  #     plotOutput(info$id, width = pw, height = ph)
  #   )
  # }
  # 
  # output$manual_diag_plots_ui <- renderUI({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) first, then diagnostics will appear."))
  #   
  #   keys <- get_diag_plot_order()
  #   if (length(keys) == 0) {
  #     return(tags$div(
  #       style = "padding:10px;border:1px dashed #ccc;border-radius:10px;background:#fff;",
  #       tags$strong("No plots selected."),
  #       tags$p("Use the left panel to select plots.")
  #     ))
  #   }
  #   
  #   dims <- diag_canvas_dims()
  #   pw <- dims$pw
  #   ph <- dims$ph
  #   preset <- input$diag_layout_preset %||% "preset_current"
  #   
  #   if (preset %in% c("preset_one_col", "preset_custom")) {
  #     return(tagList(lapply(keys, function(k) make_diag_plot_output(k, pw, ph))))
  #   }
  #   
  #   # Two-column
  #   left  <- keys[seq(1, length(keys), by = 2)]
  #   right <- keys[seq(2, length(keys), by = 2)]
  #   fluidRow(
  #     column(6, lapply(left,  function(k) make_diag_plot_output(k, pw, ph))),
  #     column(6, lapply(right, function(k) make_diag_plot_output(k, pw, ph)))
  #   )
  # })
  # 
  # output$manual_diag_canvas_ui <- renderUI({
  #   validate(need(input$fit_manual > 0, "Click 'Fit' (Manual SARIMA) first, then diagnostics will appear."))
  #   
  #   dims <- diag_canvas_dims()
  #   keys <- get_diag_plot_order()
  #   
  #   if (length(keys) == 0) return(NULL)
  #   
  #   tags$div(
  #     style = sprintf("width:%dpx; min-height:%dpx; background:#fcfcfc;", dims$w, dims$h),
  #     uiOutput("manual_diag_plots_ui")
  #   )
  # })
  # 
  # # ============================================================
  # # FULL Diagnostics container UI (THIS is where width:100% becomes slider-driven)
  # # ============================================================
  # 
  # output$diag_container_ui <- renderUI({
  #   # Container width slider controls what used to be width:100%
  #   container_w <- suppressWarnings(as.integer(input$diag_container_width_px))
  #   if (!is.finite(container_w) || container_w <= 0) container_w <- 1600L
  #   
  #   tags$div(
  #     style = sprintf("
  #     display:flex;
  #     gap:12px;
  #     align-items:flex-start;
  #     width:%dpx;          /* <-- THIS replaces width:100% */
  #     max-width:none;
  #   ", container_w),
  #     
  #     # Sidebar
  #     tags$div(
  #       style = "
  #       flex:0 0 320px;
  #       max-width:320px;
  #       border:1px solid #e5e5e5;
  #       border-radius:10px;
  #       padding:10px;
  #       background:#fff;
  #       max-height:85vh;
  #       overflow-y:auto;
  #     ",
  #       
  #       tags$h4("Diagnostics layout builder"),
  #       
  #       selectInput(
  #         "diag_layout_preset",
  #         "Layout preset",
  #         choices = c(
  #           "Current layout (2x2 + Ljung-Box + Conclusion)" = "preset_current",
  #           "Two columns (all plots, stacked)"             = "preset_two_col",
  #           "Single column (all plots, stacked)"           = "preset_one_col",
  #           "Custom order (choose positions below)"        = "preset_custom"
  #         ),
  #         selected = input$diag_layout_preset %||% "preset_current"
  #       ),
  #       
  #       tags$hr(),
  #       
  #       checkboxGroupInput(
  #         "diag_plots_selected",
  #         "Plots to display",
  #         choices = c(
  #           "Residuals over time"                          = "ts",
  #           "Residual ACF"                                 = "acf",
  #           "Residual histogram"                           = "hist",
  #           "Normal Q-Q"                                   = "qq",
  #           "Residuals vs fitted"                          = "fitted",
  #           "Scale-location (sqrt(|res|) vs fitted)"       = "scale",
  #           "Ljung–Box p-values by lag"                    = "lb"
  #         ),
  #         selected = input$diag_plots_selected %||% c("ts","acf","hist","qq","lb","fitted","scale")
  #       ),
  #       
  #       tags$hr(),
  #       
  #       # The slider you asked for: controls the flex container width
  #       sliderInput(
  #         "diag_container_width_px",
  #         "Diagnostics container width (px)",
  #         min = 900, max = 6000, value = container_w, step = 50
  #       ),
  #       helpText("This changes the flex container width that used to be width:100%. Scroll horizontally if it exceeds your screen."),
  #       
  #       tags$hr(),
  #       
  #       # ONLY TWO plot sliders
  #       sliderInput(
  #         "diag_plot_width",
  #         "Plot width (px)",
  #         min = 300, max = 2500, value = input$diag_plot_width %||% 650, step = 25
  #       ),
  #       sliderInput(
  #         "diag_plot_height",
  #         "Plot height (px)",
  #         min = 200, max = 1400, value = input$diag_plot_height %||% 320, step = 10
  #       ),
  #       
  #       tags$hr(),
  #       
  #       checkboxInput(
  #         "diag_show_conclusion",
  #         "Show academic conclusion panel",
  #         value = isTRUE(input$diag_show_conclusion %||% TRUE)
  #       ),
  #       
  #       conditionalPanel(
  #         condition = "input.diag_layout_preset == 'preset_custom'",
  #         tags$h5("Custom positions (1 = first)"),
  #         helpText("Give each enabled plot a position. Ties are broken by name."),
  #         numericInput("pos_ts",     "Residuals over time position", value = input$pos_ts %||% 1, min = 1, step = 1),
  #         numericInput("pos_acf",    "Residual ACF position",        value = input$pos_acf %||% 2, min = 1, step = 1),
  #         numericInput("pos_hist",   "Histogram position",           value = input$pos_hist %||% 3, min = 1, step = 1),
  #         numericInput("pos_qq",     "Q-Q position",                 value = input$pos_qq %||% 4, min = 1, step = 1),
  #         numericInput("pos_fitted", "Residuals vs fitted position", value = input$pos_fitted %||% 5, min = 1, step = 1),
  #         numericInput("pos_scale",  "Scale-location position",      value = input$pos_scale %||% 6, min = 1, step = 1),
  #         numericInput("pos_lb",     "Ljung–Box p-values position",  value = input$pos_lb %||% 7, min = 1, step = 1)
  #       )
  #     ),
  #     
  #     # Right panel (scrolls)
  #     tags$div(
  #       style = "
  #       flex:1 1 auto;
  #       border:1px solid #e5e5e5;
  #       border-radius:10px;
  #       background:#fcfcfc;
  #       padding:10px;
  #       overflow:auto;
  #       max-height:85vh;
  #     ",
  #       
  #       tags$div(
  #         style = "margin-bottom:12px;",
  #         tags$h4("Residual diagnostics (Manual SARIMA)"),
  #         tags$p("Use the container width + plot width/height sliders. Scroll right/down to see everything.")
  #       ),
  #       
  #       uiOutput("manual_diag_canvas_ui"),
  #       
  #       conditionalPanel(
  #         condition = "input.diag_show_conclusion == true",
  #         tags$h4("Academic conclusion (Diagnostics)"),
  #         tags$div(
  #           style = "padding:10px;border:1px solid #e5e5e5;border-radius:10px;background:#ffffff;",
  #           verbatimTextOutput("manual_diag_commentary_diag")
  #         )
  #       )
  #     )
  #   )
  # })
  
  
  
  
  # ============================================================
  # ✅ DIAGNOSTICS (COPY/PASTE COMPLETE SERVER CODE)
  # Fixes:
  # - NO more "valeur manquante là où TRUE/FALSE est requis"
  # - safe handling for NULL inputs (fit_manual, preset, sliders)
  # - Diagnostics-only Ljung–Box p-values plot: manual_resid_lb_pvals_diag2
  # - Keeps Academic conclusion plot name unchanged: manual_resid_lb_pvals
  # - Only TWO plot sliders used: diag_plot_width, diag_plot_height
  # - Container width slider: diag_container_width_px controls the flex width
  #
  # REQUIRED UI:
  # - tabPanel("Diagnostics", uiOutput("diag_container_ui"))
  # - plus: tags$head(tags$style(HTML("body { overflow-x: auto; }")))
  # ============================================================
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # ---------- SAFE helpers (prevent logical(0) / NA in if/need) ----------
  safe_chr1 <- function(x, default) {
    if (is.null(x) || length(x) == 0 || is.na(x[1])) return(default)
    as.character(x[1])
  }
  
  safe_int1 <- function(x, default) {
    y <- suppressWarnings(as.integer(x))
    if (is.null(y) || length(y) == 0 || !is.finite(y[1])) return(default)
    y[1]
  }
  
  fit_manual_clicked <- function(input) {
    fm <- input$fit_manual
    isTRUE(!is.null(fm) && length(fm) > 0 && is.finite(fm) && fm > 0)
  }
  
  nice_par <- function() {
    par(mar = c(4, 4, 2.2, 1), mgp = c(2.2, 0.7, 0), las = 1)
  }
  
  # ---- Shared helpers ----
  get_diag_controls <- function(input) {
    L_input <- suppressWarnings(as.integer(input$diag_lag))
    L <- if (is.finite(L_input) && L_input > 0) L_input else 12L
    
    alpha <- suppressWarnings(as.numeric(input$alphaSt2))
    if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05
    
    list(L = L, alpha = alpha)
  }
  
  compute_lb_pvals <- function(res, L, fitdf) {
    res <- as.numeric(res)
    res <- res[is.finite(res)]
    
    N <- length(res)
    if (N < 8) return(list(pvals = numeric(0), L = 0L, N = N))
    
    L <- as.integer(L)
    L <- max(1L, min(L, N - 1L))
    
    pvals <- rep(NA_real_, L)
    for (k in seq_len(L)) {
      if (k > fitdf) {
        bt <- tryCatch(
          stats::Box.test(res, lag = k, type = "Ljung-Box", fitdf = fitdf),
          error = function(e) NULL
        )
        pvals[k] <- if (!is.null(bt)) as.numeric(bt$p.value) else NA_real_
      }
    }
    list(pvals = pvals, L = L, N = N)
  }
  
  plot_lb_pvals <- function(pvals, L, alpha, fitdf, main_title) {
    plot(seq_len(L), pvals,
         type = "h", lwd = 2,
         xlab = "Lag (k)",
         ylab = "p-value (Ljung–Box up to lag k)",
         main = main_title,
         ylim = c(0, 1))
    points(seq_len(L), pvals, pch = 16)
    abline(h = alpha, lty = 2, col = "gray40")
    mtext(sprintf("alpha = %.3f, fitdf = %d", alpha, fitdf),
          side = 3, line = 0.2, cex = 0.85)
    
    if (is.finite(fitdf) && fitdf >= 1 && fitdf < L) {
      rect(xleft = 0.5, ybottom = -0.02, xright = fitdf + 0.5, ytop = 1.02,
           border = NA, col = grDevices::adjustcolor("gray", alpha.f = 0.15))
      text(x = (fitdf + 1) / 2, y = 0.95, labels = "df ≤ 0 (omitted)", cex = 0.8)
    }
  }
  
  # ============================================================
  # Diagnostics plots (independent outputs)
  # ============================================================
  
  output$manual_resid_ts_diag <- renderPlot({
    validate(need(fit_manual_clicked(input), "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
    req(manual_fit())
    fit <- manual_fit()
    nice_par()
    
    res <- as.numeric(residuals(fit))
    res <- res[is.finite(res)]
    validate(need(length(res) >= 3, "Not enough residuals to plot."))
    
    plot(res, type = "l",
         main = "Residuals over time",
         xlab = "Time", ylab = "Residual")
    abline(h = 0, col = "gray40", lty = 2)
  })
  
  output$manual_resid_acf_diag <- renderPlot({
    validate(need(fit_manual_clicked(input), "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
    req(manual_fit())
    fit <- manual_fit()
    nice_par()
    
    res <- as.numeric(residuals(fit))
    res <- res[is.finite(res)]
    validate(need(length(res) >= 5, "Not enough residuals for ACF."))
    
    plot(acf(res, plot = FALSE), main = "Residual ACF")
  })
  
  output$manual_resid_hist_diag <- renderPlot({
    validate(need(fit_manual_clicked(input), "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
    req(manual_fit())
    fit <- manual_fit()
    nice_par()
    
    res <- as.numeric(residuals(fit))
    res <- res[is.finite(res)]
    validate(need(length(res) >= 5, "Not enough residuals for histogram."))
    
    hist(res, breaks = 30, col = "gray85", border = "white",
         main = "Residual histogram", xlab = "Residual")
  })
  
  output$manual_resid_qq_diag <- renderPlot({
    validate(need(fit_manual_clicked(input), "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
    req(manual_fit())
    fit <- manual_fit()
    nice_par()
    
    res <- as.numeric(residuals(fit))
    res <- res[is.finite(res)]
    validate(need(length(res) >= 5, "Not enough residuals for Q–Q plot."))
    
    qqnorm(res, main = "Normal Q–Q")
    qqline(res, col = "red", lwd = 2)
  })
  
  output$manual_resid_fitted_diag <- renderPlot({
    validate(need(fit_manual_clicked(input), "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
    req(manual_fit())
    fit <- manual_fit()
    nice_par()
    
    res <- as.numeric(residuals(fit))
    fv  <- as.numeric(fitted(fit))
    ok  <- is.finite(res) & is.finite(fv)
    validate(need(sum(ok) >= 10, "Not enough valid points for residuals vs fitted."))
    
    plot(fv[ok], res[ok],
         pch = 16, cex = 0.7,
         col = grDevices::adjustcolor("steelblue", alpha.f = 0.7),
         main = "Residuals vs fitted",
         xlab = "Fitted values", ylab = "Residual")
    abline(h = 0, col = "gray40", lty = 2)
    if (sum(ok) > 20) {
      lines(stats::lowess(fv[ok], res[ok], f = 2/3), col = "tomato", lwd = 2)
    }
  })
  
  output$manual_resid_scale_diag <- renderPlot({
    validate(need(fit_manual_clicked(input), "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
    req(manual_fit())
    fit <- manual_fit()
    nice_par()
    
    res <- as.numeric(residuals(fit))
    fv  <- as.numeric(fitted(fit))
    ok  <- is.finite(res) & is.finite(fv)
    y   <- sqrt(abs(res))
    ok2 <- ok & is.finite(y)
    validate(need(sum(ok2) >= 10, "Not enough valid points for scale-location plot."))
    
    plot(fv[ok2], y[ok2],
         pch = 16, cex = 0.7,
         col = grDevices::adjustcolor("darkgreen", alpha.f = 0.65),
         main = "Scale-location (sqrt(|res|) vs fitted)",
         xlab = "Fitted values", ylab = expression(sqrt("|Residual|")))
    if (sum(ok2) > 20) {
      lines(stats::lowess(fv[ok2], y[ok2], f = 2/3), col = "tomato", lwd = 2)
    }
  })
  
  # ------------------------------------------------------------
  # NEW: Ljung–Box p-values by lag (Diagnostics ONLY)
  # ------------------------------------------------------------
  output$manual_resid_lb_pvals_diag2 <- renderPlot({
    validate(need(fit_manual_clicked(input), "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
    req(manual_fit())
    fit <- manual_fit()
    nice_par()
    
    res   <- residuals(fit)
    fitdf <- length(coef(fit))
    
    ctrl <- get_diag_controls(input)
    out  <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
    
    validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
    validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
    
    plot_lb_pvals(out$pvals, out$L, ctrl$alpha, fitdf,
                  main_title = "Ljung–Box p-values by lag (Diagnostics)")
  })
  
  # ------------------------------------------------------------
  # Commentary (Diagnostics tab)
  # ------------------------------------------------------------
  output$manual_diag_commentary_diag <- renderPrint({
    validate(need(fit_manual_clicked(input), "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
    req(manual_fit())
    fit <- manual_fit()
    
    res <- as.numeric(residuals(fit))
    res <- res[is.finite(res)]
    N   <- length(res)
    
    ctrl  <- get_diag_controls(input)
    fitdf <- length(coef(fit))
    
    if (N < 8) {
      cat("Diagnostics summary (Manual SARIMA)\n",
          "------------------------------------------------------------\n",
          "- Too few residuals to run diagnostics.\n", sep = "")
      return(invisible())
    }
    
    L2 <- max(1L, min(as.integer(ctrl$L), N - 1L))
    lb <- tryCatch(
      stats::Box.test(res, lag = L2, type = "Ljung-Box", fitdf = fitdf),
      error = function(e) NULL
    )
    
    cat("Diagnostics summary (Manual SARIMA)\n")
    cat("------------------------------------------------------------\n")
    if (!is.null(lb) && is.finite(lb$p.value)) {
      cat("Ljung–Box (lag ", L2, ") p = ", signif(lb$p.value, 3),
          if (lb$p.value > ctrl$alpha) " -> no strong evidence of autocorrelation.\n"
          else " -> remaining autocorrelation; revise orders/differencing.\n",
          sep = "")
    } else {
      cat("Ljung–Box unavailable.\n")
    }
    
    cat("\nNotes:\n")
    cat("- Use Residual ACF + Ljung–Box together: spikes + small p-values suggest underfitting.\n")
    cat("- If Q–Q shows heavy tails: intervals/inference may be optimistic.\n")
  })
  
  # ============================================================
  # KEEP: Academic conclusion plot output unchanged
  # ============================================================
  output$manual_resid_lb_pvals <- renderPlot({
    validate(need(fit_manual_clicked(input), "Click 'Fit' (Manual SARIMA) to generate diagnostics."))
    req(manual_fit())
    fit <- manual_fit()
    nice_par()
    
    res   <- residuals(fit)
    fitdf <- length(coef(fit))
    
    ctrl <- get_diag_controls(input)
    out  <- compute_lb_pvals(res, L = ctrl$L, fitdf = fitdf)
    
    validate(need(out$N >= 8, "Too few residuals (N < 8) to compute Ljung–Box p-values."))
    validate(need(out$L >= 1, "No valid lags available for Ljung–Box p-values."))
    
    plot_lb_pvals(out$pvals, out$L, ctrl$alpha, fitdf,
                  main_title = "Ljung–Box p-values by lag (Manual SARIMA)")
  })
  
  # ============================================================
  # Layout builder + auto canvas sizing from ONLY plot width/height
  # ============================================================
  
  manual_diag_plot_map <- list(
    ts     = list(id = "manual_resid_ts_diag",        title = "Residuals over time"),
    acf    = list(id = "manual_resid_acf_diag",       title = "Residual ACF"),
    hist   = list(id = "manual_resid_hist_diag",      title = "Residual histogram"),
    qq     = list(id = "manual_resid_qq_diag",        title = "Normal Q–Q"),
    fitted = list(id = "manual_resid_fitted_diag",    title = "Residuals vs fitted"),
    scale  = list(id = "manual_resid_scale_diag",     title = "Scale-location"),
    lb     = list(id = "manual_resid_lb_pvals_diag2", title = "Ljung–Box p-values by lag (Diagnostics)")
  )
  
  get_diag_plot_order <- reactive({
    sel <- input$diag_plots_selected
    if (is.null(sel) || length(sel) == 0) return(character(0))
    
    preset <- safe_chr1(input$diag_layout_preset, "preset_current")
    
    if (!identical(preset, "preset_custom")) {
      default_order <- c("ts", "acf", "hist", "qq", "fitted", "scale", "lb")
      return(default_order[default_order %in% sel])
    }
    
    pos <- c(
      ts     = input$pos_ts,
      acf    = input$pos_acf,
      hist   = input$pos_hist,
      qq     = input$pos_qq,
      fitted = input$pos_fitted,
      scale  = input$pos_scale,
      lb     = input$pos_lb
    )
    pos <- pos[names(pos) %in% sel]
    pos_num <- suppressWarnings(as.numeric(pos))
    pos_num[!is.finite(pos_num)] <- 999
    ord <- order(pos_num, names(pos_num))
    names(pos_num)[ord]
  })
  
  diag_canvas_dims <- reactive({
    keys   <- get_diag_plot_order()
    preset <- safe_chr1(input$diag_layout_preset, "preset_current")
    
    pw <- safe_int1(input$diag_plot_width, 650L)
    ph <- safe_int1(input$diag_plot_height, 320L)
    
    # If you want to allow "auto width", set slider min >=300; we keep px always.
    if (pw < 300) pw <- 300L
    if (ph < 200) ph <- 200L
    
    ncol <- if (identical(preset, "preset_one_col") || identical(preset, "preset_custom")) 1L else 2L
    n <- length(keys)
    nrow <- if (n == 0) 1L else ceiling(n / ncol)
    
    pad_w <- 60L
    pad_h <- 110L
    
    canvas_w <- ncol * (pw + pad_w) + 40L
    canvas_h <- nrow * (ph + pad_h) + 40L
    
    list(w = canvas_w, h = canvas_h, pw = pw, ph = ph)
  })
  
  make_diag_plot_output <- function(plot_key, pw, ph) {
    info <- manual_diag_plot_map[[plot_key]]
    if (is.null(info)) return(NULL)
    
    tags$div(
      style = sprintf(
        "margin-bottom:14px; padding:10px; border:1px solid #eee; border-radius:10px; background:#fff; width:%dpx;",
        pw + 30
      ),
      tags$h5(info$title),
      plotOutput(info$id, width = pw, height = ph)
    )
  }
  
  output$manual_diag_plots_ui <- renderUI({
    validate(need(fit_manual_clicked(input), "Click 'Fit' (Manual SARIMA) first, then diagnostics will appear."))
    
    keys <- get_diag_plot_order()
    if (length(keys) == 0) {
      return(tags$div(
        style = "padding:10px;border:1px dashed #ccc;border-radius:10px;background:#fff;",
        tags$strong("No plots selected."),
        tags$p("Use the left panel to select plots.")
      ))
    }
    
    dims <- diag_canvas_dims()
    pw <- dims$pw
    ph <- dims$ph
    preset <- safe_chr1(input$diag_layout_preset, "preset_current")
    
    if (identical(preset, "preset_one_col") || identical(preset, "preset_custom")) {
      return(tagList(lapply(keys, function(k) make_diag_plot_output(k, pw, ph))))
    }
    
    left  <- keys[seq(1, length(keys), by = 2)]
    right <- keys[seq(2, length(keys), by = 2)]
    fluidRow(
      column(6, lapply(left,  function(k) make_diag_plot_output(k, pw, ph))),
      column(6, lapply(right, function(k) make_diag_plot_output(k, pw, ph)))
    )
  })
  
  output$manual_diag_canvas_ui <- renderUI({
    validate(need(fit_manual_clicked(input), "Click 'Fit' (Manual SARIMA) first, then diagnostics will appear."))
    
    dims <- diag_canvas_dims()
    keys <- get_diag_plot_order()
    if (length(keys) == 0) return(NULL)
    
    tags$div(
      style = sprintf("width:%dpx; min-height:%dpx; background:#fcfcfc;", dims$w, dims$h),
      uiOutput("manual_diag_plots_ui")
    )
  })
  
  # ============================================================
  # FULL Diagnostics container UI (slider-driven width)
  # ============================================================
  
  output$diag_container_ui <- renderUI({
    container_w <- safe_int1(input$diag_container_width_px, 1600L)
    if (container_w < 900) container_w <- 900L
    
    # keep selected values stable even if UI rebuilds
    preset_now <- safe_chr1(input$diag_layout_preset, "preset_current")
    sel_now    <- input$diag_plots_selected %||% c("ts","acf","hist","qq","lb","fitted","scale")
    pw_now     <- safe_int1(input$diag_plot_width, 450L)
    ph_now     <- safe_int1(input$diag_plot_height, 320L)
    show_conc  <- isTRUE(input$diag_show_conclusion %||% TRUE)
    
    tags$div(
      style = sprintf("
      display:flex;
      gap:12px;
      align-items:flex-start;
      width:%dpx;
      max-width:none;
    ", container_w),
      
      # Sidebar
      tags$div(
        style = "
        flex:0 0 220px;
        max-width:320px;
        border:1px solid #e5e5e5;
        border-radius:10px;
        padding:10px;
        background:#fff;
        max-height:85vh;
        overflow-y:auto;
      ",
        
        tags$h4("Diagnostics layout builder"),
        
        selectInput(
          "diag_layout_preset",
          "Layout preset",
          choices = c(
            "Current layout" = "preset_current",
            "Two columns"    = "preset_two_col",
            "Single column"  = "preset_one_col",
            "Custom order"   = "preset_custom"
          ),
          selected = preset_now
        ),
        
        tags$hr(),
        
        checkboxGroupInput(
          "diag_plots_selected",
          "Plots to display",
          choices = c(
            "Residuals over time"                          = "ts",
            "Residual ACF"                                 = "acf",
            "Residual histogram"                           = "hist",
            "Normal Q-Q"                                   = "qq",
            "Residuals vs fitted"                          = "fitted",
            "Scale-location (sqrt(|res|) vs fitted)"       = "scale",
            "Ljung–Box p-values by lag"                    = "lb"
          ),
          selected = sel_now
        ),
        
        tags$hr(),
        
        sliderInput(
          "diag_container_width_px",
          "Container width (px)",
          min = 1000, max = 5000, value = container_w, step = 50
        ),
        # helpText("Controls the whole Diagnostics width (was width:100%). Scroll horizontally if needed."),
        
        tags$hr(),
        
        # ONLY TWO plot sliders
        sliderInput(
          "diag_plot_width",
          "Plot width (px)",
          min = 300, max = 2500, value = pw_now, step = 25
        ),
        sliderInput(
          "diag_plot_height",
          "Plot height (px)",
          min = 200, max = 1400, value = ph_now, step = 10
        ),
        
        tags$hr(),
        
        checkboxInput("diag_show_conclusion", "Show academic conclusion panel", value = show_conc),
        
        conditionalPanel(
          condition = "input.diag_layout_preset == 'preset_custom'",
          tags$h5("Custom positions (1 = first)"),
          helpText("Give each enabled plot a position. Ties are broken by name."),
          numericInput("pos_ts",     "Residuals over time position", value = safe_int1(input$pos_ts, 1L), min = 1, step = 1),
          numericInput("pos_acf",    "Residual ACF position",        value = safe_int1(input$pos_acf, 2L), min = 1, step = 1),
          numericInput("pos_hist",   "Histogram position",           value = safe_int1(input$pos_hist, 3L), min = 1, step = 1),
          numericInput("pos_qq",     "Q-Q position",                 value = safe_int1(input$pos_qq, 4L), min = 1, step = 1),
          numericInput("pos_fitted", "Residuals vs fitted position", value = safe_int1(input$pos_fitted, 5L), min = 1, step = 1),
          numericInput("pos_scale",  "Scale-location position",      value = safe_int1(input$pos_scale, 6L), min = 1, step = 1),
          numericInput("pos_lb",     "Ljung–Box p-values position",  value = safe_int1(input$pos_lb, 7L), min = 1, step = 1)
        )
      ),
      
      # Right panel
      tags$div(
        style = "
        flex:1 1 auto;
        border:1px solid #e5e5e5;
        border-radius:10px;
        background:#fcfcfc;
        padding:10px;
        overflow:auto;
        max-height:85vh;
      ",
        
        tags$div(
          style = "margin-bottom:12px;",
          tags$h4("Residual diagnostics (Manual SARIMA)"),
          tags$p("Use container width + plot width/height. Scroll right/down to see everything.")
        ),
        
        uiOutput("manual_diag_canvas_ui"),
        
        conditionalPanel(
          condition = "input.diag_show_conclusion == true",
          tags$h4("Academic conclusion (Diagnostics)"),
          tags$div(
            style = "padding:10px;border:1px solid #e5e5e5;border-radius:10px;background:#ffffff;",
            verbatimTextOutput("manual_diag_commentary_diag")
          )
        )
      )
    )
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------- 
  # ---------------------------------------------   # --------------------------------------------- 
  # ---------------------------------------------   # --------------------------------------------- 
  # ---------------------------------------------   # --------------------------------------------- 
  # ---------------------------------------------   # --------------------------------------------- 
  # ---------------------------------------------   # --------------------------------------------- 
  # ---------------------------------------------   # --------------------------------------------- 
  # ---------------------------------------------   # --------------------------------------------- 
  # ---------------------------------------------   # --------------------------------------------- 
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
  
  
  
  
  
 
  
  # ============================================================
  # --- MOD: Render the equation panel with LEFT alignment (CSS) and headings ---
  #   NOTE: Left alignment is done via an HTML wrapper div.
  # ============================================================
  
  # --- helper: scientific formatting (uppercase E), preserves NA ---
  fmt_sci <- function(x, digits = 4) {
    ifelse(is.na(x),
           NA_character_,
           toupper(formatC(x, format = "e", digits = digits)))
  }
  
  # --- helper: show scientific only when needed (tiny/huge), else fixed ---
  fmt_auto <- function(x, digits_fixed = 6, digits_sci = 4) {
    ifelse(
      is.na(x),
      NA_character_,
      ifelse(abs(x) > 0 & (abs(x) < 1e-4 | abs(x) >= 1e5),
             toupper(formatC(x, format = "e", digits = digits_sci)),
             formatC(x, format = "fg", digits = digits_fixed, flag = "#"))
    )
  }
  
  # Parameter significance table for Manual SARIMA (robust to tiny numbers)
  output$manual_param_table <- renderTable({
    req(manual_fit())
    fit <- manual_fit()
    
    # 1) Coefficients
    est <- tryCatch(stats::coef(fit), error = function(e) NULL)
    validate(need(!is.null(est) && length(est) > 0, "No estimated parameters available."))
    
    # 2) Covariance → SE (robust fallbacks)
    V <- tryCatch(stats::vcov(fit), error = function(e) NULL)
    if (is.null(V)) V <- tryCatch(fit$var.coef, error = function(e) NULL)  # forecast::Arima stores var.coef
    se <- if (!is.null(V)) sqrt(diag(V)) else rep(NA_real_, length(est))
    
    # 3) Test stats and p-values (normal/Z approx)
    tst <- est / se
    pvl <- 2 * stats::pnorm(abs(tst), lower.tail = FALSE)
    
    # 4) Nice names (optional)
    s <- suppressWarnings(tryCatch(manual_equations()$s, error = function(e) NA_integer_))
    map_name <- function(nm) {
      nm <- gsub("^ar(\\d+)$", "AR{\\1}", nm, ignore.case = TRUE)
      nm <- gsub("^ma(\\d+)$", "MA{\\1}", nm, ignore.case = TRUE)
      if (isTRUE(!is.na(s))) {
        nm <- gsub("^sar\\d+$", paste0("SAR{", s, "}"), nm, ignore.case = TRUE)
        nm <- gsub("^sma\\d+$", paste0("SMA{", s, "}"), nm, ignore.case = TRUE)
      } else {
        nm <- gsub("^sar\\d+$", "SAR", nm, ignore.case = TRUE)
        nm <- gsub("^sma\\d+$", "SMA", nm, ignore.case = TRUE)
      }
      nm <- gsub("^intercept$", "Constant", nm, ignore.case = TRUE)
      nm <- gsub("^mean$",      "Constant", nm, ignore.case = TRUE)
      nm <- gsub("^drift$",     "Drift",    nm, ignore.case = TRUE)
      nm
    }
    
    df_num <- data.frame(
      Parameter        = vapply(names(est), map_name, character(1)),
      Value            = as.numeric(est),
      `Standard Error` = as.numeric(se),
      `t Statistic`    = as.numeric(tst),
      `P-Value`        = as.numeric(pvl),
      check.names = FALSE
    )
    
    # 5) Append variance (no significance test)
    sigma2 <- suppressWarnings(as.numeric(fit$sigma2))
    if (is.finite(sigma2)) {
      df_num <- rbind(
        df_num,
        data.frame(Parameter = "Variance",
                   Value = sigma2,
                   `Standard Error` = NA_real_,
                   `t Statistic` = NA_real_,
                   `P-Value` = NA_real_,
                   check.names = FALSE)
      )
    }
    
    # 6) FORMAT: keep numbers that need scientific notation in E form (e.g., 5.2E-12)
    df_out <- transform(
      df_num,
      Value            = fmt_auto(Value),
      `Standard Error` = fmt_auto(`Standard Error`),
      `t Statistic`    = fmt_auto(`t Statistic`),
      `P-Value`        = fmt_sci(`P-Value`, digits = 3)  # always scientific for p-values
    )
    
    df_out
  }, rownames = FALSE)
  
  # Parameter significance table for Manual SARIMA
  # output$manual_param_table <- renderTable({
  #   req(manual_fit())
  #   
  #   fit <- manual_fit()
  #   
  #   # estimates and covariance (works for forecast::Arima or stats::arima)
  #   est <- tryCatch(coef(fit), error = function(e) NULL)
  #   V   <- tryCatch(vcov(fit), error = function(e) NULL)
  #   
  #   validate(need(!is.null(est) && length(est) > 0, "No estimated parameters available."))
  #   
  #   se <- if (!is.null(V)) sqrt(diag(V)) else rep(NA_real_, length(est))
  #   tst <- est / se
  #   pvl <- 2 * pnorm(abs(tst), lower.tail = FALSE)   # large-sample normal approx
  #   
  #   # friendly parameter names
  #   s <- tryCatch(manual_equations()$s, error = function(e) NA_integer_)
  #   map_name <- function(nm) {
  #     nm <- gsub("^ar(\\d+)$", "AR{\\1}", nm, ignore.case = TRUE)
  #     nm <- gsub("^ma(\\d+)$", "MA{\\1}", nm, ignore.case = TRUE)
  #     if (isTRUE(!is.na(s))) {
  #       nm <- gsub("^sar\\d+$", paste0("SAR{", s, "}"), nm, ignore.case = TRUE)
  #       nm <- gsub("^sma\\d+$", paste0("SMA{", s, "}"), nm, ignore.case = TRUE)
  #     } else {
  #       nm <- gsub("^sar\\d+$", "SAR", nm, ignore.case = TRUE)
  #       nm <- gsub("^sma\\d+$", "SMA", nm, ignore.case = TRUE)
  #     }
  #     nm <- gsub("^intercept$", "Constant", nm, ignore.case = TRUE)
  #     nm <- gsub("^mean$",      "Constant", nm, ignore.case = TRUE)
  #     nm <- gsub("^drift$",     "Drift",    nm, ignore.case = TRUE)
  #     nm
  #   }
  #   
  #   out <- data.frame(
  #     Parameter       = vapply(names(est), map_name, character(1)),
  #     Value           = as.numeric(est),
  #     `Standard Error`= as.numeric(se),
  #     `t Statistic`   = as.numeric(tst),
  #     `P-Value`       = as.numeric(pvl),
  #     check.names = FALSE
  #   )
  #   
  #   # Append Variance row (no significance test provided)
  #   sigma2 <- tryCatch(as.numeric(fit$sigma2), error = function(e) NA_real_)
  #   if (is.finite(sigma2)) {
  #     out <- rbind(
  #       out,
  #       data.frame(Parameter = "Variance",
  #                  Value = sigma2,
  #                  `Standard Error` = NA_real_,
  #                  `t Statistic` = NA_real_,
  #                  `P-Value` = NA_real_,
  #                  check.names = FALSE)
  #     )
  #   }
  #   
  #   # tidy formatting
  #   num_cols <- c("Value", "Standard Error", "t Statistic", "P-Value")
  #   for (cc in intersect(names(out), num_cols)) {
  #     out[[cc]] <- ifelse(is.na(out[[cc]]), NA, signif(out[[cc]], 6))
  #   }
  #   
  #   out
  # }, rownames = FALSE)
  
  
  # Goodness-of-fit table for Manual SARIMA
  output$manual_gof_table <- renderTable({
    req(manual_fit())
    
    fit <- manual_fit()
    
    # sample size and parameter count
    n <- tryCatch(length(residuals(fit)), error = function(e) NA_integer_)
    k <- tryCatch(length(coef(fit)),      error = function(e) NA_integer_)
    
    # AIC
    AIC_val <- suppressWarnings(tryCatch(stats::AIC(fit), error = function(e) NA_real_))
    
    # BIC (use generic first; if unavailable, compute from logLik)
    BIC_val <- suppressWarnings(tryCatch(stats::BIC(fit), error = function(e) NA_real_))
    if (!is.finite(BIC_val)) {
      ll <- suppressWarnings(tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_))
      if (is.finite(ll) && is.finite(k) && is.finite(n) && n > 0) {
        BIC_val <- (-2 * ll) + k * log(n)
      }
    }
    
    # AICc (use forecast::AICc if available; otherwise use formula)
    AICc_val <- suppressWarnings(tryCatch(forecast::AICc(fit), error = function(e) NA_real_))
    if (!is.finite(AICc_val) && is.finite(AIC_val) && is.finite(k) && is.finite(n) && (n - k - 1) > 0) {
      AICc_val <- AIC_val + (2 * k * (k + 1)) / (n - k - 1)
    }
    
    # Assemble table
    out <- data.frame(
      Metric = c("AIC", "AICc", "BIC"),
      Value  = c(AIC_val, AICc_val, BIC_val),
      check.names = FALSE
    )
    
    # Numeric formatting
    out$Value <- ifelse(is.na(out$Value), NA, signif(out$Value, 6))
    out
  }, rownames = FALSE)
  
  
  
  output$manual_model_equation <- renderUI({
    req(manual_equations())
    eq <- manual_equations()
    
    tagList(
      tags$div(
        style = "text-align:left;",
        
        tags$h4("Manual SARIMA model"),
        tags$p(sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)[%d]", eq$p, eq$d, eq$q, eq$P, eq$D, eq$Q, eq$s)),
        
        tags$h5("Estimated coefficients"),
        tags$ul(lapply(eq$coef_lines, function(x) tags$li(HTML(x)))),
        
        ## parameter significance table
        tags$hr(),
        tags$h4("Table: Estimation Results"),
        tableOutput("manual_param_table"), 
        
        
        tags$hr(),
        tags$h4("Table: Goodness of Fit"),
        tableOutput("manual_gof_table"),
        
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
        HTML(eq$eq_line4),
        
        tags$hr(),
        
  
      ),
      
      # keep MathJax refresh
      tags$script(HTML("if (window.MathJax && MathJax.Hub) MathJax.Hub.Queue(['Typeset', MathJax.Hub]);"))
    )
  })
  
  
  
  
  
  # output$manual_model_equation <- renderUI({
  #   req(manual_equations())
  #   eq <- manual_equations()
  #   
  #   tagList(
  #     # --- MOD: left alignment wrapper for all MathJax blocks ---
  #     tags$div(
  #       style = "text-align:left;",
  #       
  #       tags$h4("Manual SARIMA model"),
  #       tags$p(sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)[%d]", eq$p, eq$d, eq$q, eq$P, eq$D, eq$Q, eq$s)),
  #       
  #       tags$h5("Estimated coefficients"),
  #       tags$ul(lapply(eq$coef_lines, function(x) tags$li(HTML(x)))),
  #       
  #       tags$hr(),
  #       
  #       tags$h4("General SARIMA formulation"),
  #       HTML(eq$eq_general),
  #       
  #       tags$hr(),
  #       
  #       tags$h4("Expanded operator form"),
  #       HTML(eq$eq_expanded),
  #       
  #       tags$hr(),
  #       
  #       tags$h4("Numerical model"),
  #       # --- MOD: show line 3 then line 4 under the same heading ---
  #       HTML(eq$eq_line3),
  #       tags$hr(),
  #       
  #       # HTML(tex_display("\\text{------------}")),
  #       HTML(eq$eq_line4),
  #       
  #       tags$hr(),
  #       
  #       # tags$h4("Equivalent time-domain representation (teaching form)"),
  #       # HTML(eq$eq_time_domain)
  #     ),
  #     
  #     # --- MOD: Force MathJax typesetting for dynamically inserted content ---
  #     tags$script(HTML("if (window.MathJax && MathJax.Hub) MathJax.Hub.Queue(['Typeset', MathJax.Hub]);"))
  #   )
  # })
  
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
  
  
  # ---- Manual SARIMA: Original-scale plot (observed + forecast + CIs) ----
  output$manual_forecast_plot_original <- renderPlot({
    req(manual_fc(), ts_train_test(), prepared())
    
    s <- ts_train_test()
    p <- prepared()
    fc <- manual_fc()$fc
    
    # observed series on ORIGINAL scale
    validate(need("y_filled" %in% names(s$dfm), "Column y_filled not found in ts_train_test()$dfm."))
    obs_df <- s$dfm[, c("x", "y_filled")]
    names(obs_df) <- c("x", "y")
    
    # build forecast df (x alignment) from existing helper, then back-transform values
    fc_df <- plot_forecast_df(obs_df, s$train_n, fc, by = p$by)
    
    # inverse transform helper (must match your global transform choice)
    inv <- function(z) {
      tr <- input$transform %||% "none"
      z <- as.numeric(z)
      
      if (tr == "none") return(z)
      
      if (tr == "log") {
        return(exp(z))
      }
      
      if (tr == "boxcox") {
        y0 <- prepared()$df$y_filled
        validate(need(all(y0 > 0, na.rm = TRUE), "Box-Cox inverse requires strictly positive original values."))
        
        lam <- input$lambda
        if (is.null(lam) || (length(lam) == 1 && is.na(lam))) {
          lam <- forecast::BoxCox.lambda(y0, method = "guerrero")
        } else {
          lam <- as.numeric(lam)
        }
        return(forecast::InvBoxCox(z, lam))
      }
      
      z
    }
    
    # back-transform forecast mean + intervals
    fc_df$mean <- inv(fc_df$mean)
    if ("lo80" %in% names(fc_df)) fc_df$lo80 <- inv(fc_df$lo80)
    if ("hi80" %in% names(fc_df)) fc_df$hi80 <- inv(fc_df$hi80)
    if ("lo95" %in% names(fc_df)) fc_df$lo95 <- inv(fc_df$lo95)
    if ("hi95" %in% names(fc_df)) fc_df$hi95 <- inv(fc_df$hi95)
    
    gg_forecast_plot(
      obs_df, s$train_n, fc_df,
      title = "Manual SARIMA forecast (original scale)"
    )
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
  
  
  
  
  
  
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  
  
  # ---- Manual SARIMA: academic conclusion (cached on Fit click) ----
  
  # =========================
  # Manual SARIMA: FULL academic conclusion object (robust + detailed)
  # =========================
  
  # manual_conclusion_full_obj <- eventReactive(input$fit_manual, {
  #   req(manual_fit(), manual_fc(), manual_equations(), ts_train_test())
  #   
  #   fit <- manual_fit()
  #   fc0 <- manual_fc()
  #   fc  <- fc0$fc
  #   eq  <- manual_equations()
  #   s   <- ts_train_test()
  #   
  #   # ---------- helpers (local + safe)
  #   fmt_num_local <- function(x, d = 3) {
  #     if (length(x) == 0 || all(is.na(x))) return("NA")
  #     x <- suppressWarnings(as.numeric(x[1]))
  #     if (!is.finite(x)) return("NA")
  #     formatC(x, format = "f", digits = d)
  #   }
  #   fmt_p_local <- function(p) {
  #     if (length(p) == 0 || all(is.na(p))) return("NA")
  #     p <- suppressWarnings(as.numeric(p[1]))
  #     if (!is.finite(p)) return("NA")
  #     if (p < .001) "&lt; .001" else sprintf("= %.3f", p)
  #   }
  #   sig_stars <- function(p) {
  #     p <- suppressWarnings(as.numeric(p))
  #     if (!is.finite(p)) return("")
  #     if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else if (p < .10) "†" else ""
  #   }
  #   safe_len <- function(x) if (is.null(x)) 0L else length(x)
  #   
  #   html_table <- function(df) {
  #     if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
  #       return(tags$em("Table unavailable."))
  #     }
  #     tags$table(
  #       class = "table table-striped table-condensed",
  #       tags$thead(
  #         tags$tr(lapply(names(df), function(nm) tags$th(nm)))
  #       ),
  #       tags$tbody(
  #         lapply(seq_len(nrow(df)), function(i) {
  #           tags$tr(lapply(df[i, , drop = FALSE], function(cell) tags$td(HTML(as.character(cell)))))
  #         })
  #       )
  #     )
  #   }
  #   
  #   # ---------- sample sizes (safe)
  #   n_train <- suppressWarnings(as.integer(s$train_n))
  #   if (!is.finite(n_train) || n_train < 1) n_train <- tryCatch(length(residuals(fit)), error = function(e) 0L)
  #   
  #   n_test <- suppressWarnings(as.integer(s$test_n))
  #   if (!is.finite(n_test) || n_test < 0) n_test <- 0L
  #   
  #   N <- n_train + n_test
  #   
  #   # ---------- lag choice (safe)
  #   L_in <- suppressWarnings(as.integer(input$diag_lag))
  #   L <- if (is.finite(L_in) && L_in > 0) L_in else 12L
  #   
  #   # ---------- IC (safe)
  #   AIC_val  <- suppressWarnings(tryCatch(as.numeric(stats::AIC(fit)), error = function(e) NA_real_))
  #   BIC_val  <- suppressWarnings(tryCatch(as.numeric(stats::BIC(fit)), error = function(e) NA_real_))
  #   AICc_val <- suppressWarnings(tryCatch(as.numeric(forecast::AICc(fit)), error = function(e) NA_real_))
  #   
  #   ic_df <- data.frame(
  #     Criterion = c("AIC", "AICc", "BIC"),
  #     Value     = c(fmt_num_local(AIC_val, 2), fmt_num_local(AICc_val, 2), fmt_num_local(BIC_val, 2)),
  #     check.names = FALSE
  #   )
  #   
  #   # ---------- coefficients + significance (robust)
  #   est <- suppressWarnings(tryCatch(stats::coef(fit), error = function(e) NULL))
  #   V   <- suppressWarnings(tryCatch(stats::vcov(fit), error = function(e) NULL))
  #   if (is.null(V)) V <- suppressWarnings(tryCatch(fit$var.coef, error = function(e) NULL))
  #   
  #   coef_df <- NULL
  #   if (!is.null(est) && length(est) > 0) {
  #     est <- as.numeric(est)
  #     nm  <- names(stats::coef(fit))
  #     if (is.null(nm)) nm <- paste0("param_", seq_along(est))
  #     
  #     se <- rep(NA_real_, length(est))
  #     if (!is.null(V)) {
  #       dV <- tryCatch(diag(V), error = function(e) rep(NA_real_, length(est)))
  #       if (length(dV) == length(est)) se <- sqrt(pmax(dV, 0))
  #     }
  #     z  <- est / se
  #     p  <- 2 * stats::pnorm(-abs(z))
  #     
  #     coef_df <- data.frame(
  #       Term     = nm,
  #       Estimate = sprintf("%.6f", est),
  #       SE       = ifelse(is.finite(se), sprintf("%.6f", se), "NA"),
  #       `z/t`    = ifelse(is.finite(z),  sprintf("%.3f",  z),  "NA"),
  #       `p`      = ifelse(is.finite(p),  sprintf("%.3f",  p),  "NA"),
  #       Sig      = vapply(p, sig_stars, character(1)),
  #       check.names = FALSE
  #     )
  #   }
  #   
  #   n_sig <- if (!is.null(coef_df)) sum(suppressWarnings(as.numeric(coef_df$p)) < 0.05, na.rm = TRUE) else 0L
  #   
  #   # ---------- residuals (safe)
  #   res <- suppressWarnings(tryCatch(as.numeric(residuals(fit)), error = function(e) numeric(0)))
  #   res <- res[is.finite(res)]
  #   n_res <- length(res)
  #   
  #   # fitdf for LB: number of estimated parameters (safe)
  #   fitdf <- if (!is.null(est)) length(est) else 0L
  #   
  #   # choose LB lag not exceeding sample
  #   lb_lag <- min(L, max(1L, floor(n_res / 3)))
  #   lb <- if (n_res >= 5) {
  #     tryCatch(stats::Box.test(res, lag = lb_lag, type = "Ljung-Box", fitdf = fitdf), error = function(e) NULL)
  #   } else NULL
  #   
  #   bp <- if (n_res >= 5) {
  #     tryCatch(stats::Box.test(res, lag = lb_lag, type = "Box-Pierce", fitdf = fitdf), error = function(e) NULL)
  #   } else NULL
  #   
  #   jb <- if (requireNamespace("tseries", quietly = TRUE) && n_res >= 5) {
  #     tryCatch(tseries::jarque.bera.test(res), error = function(e) NULL)
  #   } else NULL
  #   
  #   sw <- if (n_res >= 3 && n_res <= 5000) {
  #     tryCatch(stats::shapiro.test(res), error = function(e) NULL)
  #   } else NULL
  #   
  #   arch <- if (requireNamespace("FinTS", quietly = TRUE) && n_res >= 10) {
  #     tryCatch(FinTS::ArchTest(res, lags = min(12L, max(1L, floor(L / 2)))), error = function(e) NULL)
  #   } else NULL
  #   
  #   runs <- if (requireNamespace("tseries", quietly = TRUE) && n_res >= 10) {
  #     tryCatch(tseries::runs.test(res), error = function(e) NULL)
  #   } else NULL
  #   
  #   ad <- if (requireNamespace("nortest", quietly = TRUE) && n_res >= 8) {
  #     tryCatch(nortest::ad.test(res), error = function(e) NULL)
  #   } else NULL
  #   
  #   # ---------- compact residual-test table (academic reporting style)
  #   test_rows <- list()
  #   
  #   add_test <- function(name, stat, pval, note = "") {
  #     test_rows[[length(test_rows) + 1L]] <<- data.frame(
  #       Test = name,
  #       Statistic = if (is.null(stat)) "NA" else fmt_num_local(stat, 3),
  #       `p-value` = if (is.null(pval)) "NA" else fmt_p_local(pval),
  #       Interpretation = note,
  #       check.names = FALSE
  #     )
  #   }
  #   
  #   add_test(
  #     "Ljung–Box (residuals)",
  #     if (!is.null(lb)) unname(lb$statistic) else NULL,
  #     if (!is.null(lb)) unname(lb$p.value) else NULL,
  #     if (!is.null(lb)) {
  #       if (is.finite(lb$p.value) && lb$p.value >= 0.05) "No evidence of remaining autocorrelation (white-noise compatible)."
  #       else "Evidence of remaining autocorrelation (model may be under-specified)."
  #     } else "Not available."
  #   )
  #   
  #   add_test(
  #     "Box–Pierce (residuals)",
  #     if (!is.null(bp)) unname(bp$statistic) else NULL,
  #     if (!is.null(bp)) unname(bp$p.value) else NULL,
  #     if (!is.null(bp)) {
  #       if (is.finite(bp$p.value) && bp$p.value >= 0.05) "Consistent with uncorrelated residuals."
  #       else "Suggests residual autocorrelation."
  #     } else "Not available."
  #   )
  #   
  #   add_test(
  #     "Jarque–Bera (normality)",
  #     if (!is.null(jb)) unname(jb$statistic) else NULL,
  #     if (!is.null(jb)) unname(jb$p.value) else NULL,
  #     if (!is.null(jb)) {
  #       if (is.finite(jb$p.value) && jb$p.value >= 0.05) "No evidence against normality."
  #       else "Residuals deviate from normality (common in real series)."
  #     } else "Package 'tseries' missing or test unavailable."
  #   )
  #   
  #   add_test(
  #     "Shapiro–Wilk (normality)",
  #     if (!is.null(sw)) unname(sw$statistic) else NULL,
  #     if (!is.null(sw)) unname(sw$p.value) else NULL,
  #     if (!is.null(sw)) {
  #       if (is.finite(sw$p.value) && sw$p.value >= 0.05) "No evidence against normality."
  #       else "Evidence against normality."
  #     } else if (n_res > 5000) "Not computed (n > 5000)." else "Not available."
  #   )
  #   
  #   add_test(
  #     "ARCH LM (heteroskedasticity)",
  #     if (!is.null(arch)) unname(arch$statistic) else NULL,
  #     if (!is.null(arch)) unname(arch$p.value) else NULL,
  #     if (!is.null(arch)) {
  #       if (is.finite(arch$p.value) && arch$p.value >= 0.05) "No evidence of remaining ARCH effects."
  #       else "Evidence of ARCH effects → consider GARCH for variance."
  #     } else "Package 'FinTS' missing or test unavailable."
  #   )
  #   
  #   add_test(
  #     "Runs test (randomness)",
  #     if (!is.null(runs)) unname(runs$statistic) else NULL,
  #     if (!is.null(runs)) unname(runs$p.value) else NULL,
  #     if (!is.null(runs)) {
  #       if (is.finite(runs$p.value) && runs$p.value >= 0.05) "No evidence against randomness."
  #       else "Evidence of non-randomness (structure may remain)."
  #     } else "Package 'tseries' missing or test unavailable."
  #   )
  #   
  #   add_test(
  #     "Anderson–Darling (normality)",
  #     if (!is.null(ad)) unname(ad$statistic) else NULL,
  #     if (!is.null(ad)) unname(ad$p.value) else NULL,
  #     if (!is.null(ad)) {
  #       if (is.finite(ad$p.value) && ad$p.value >= 0.05) "No evidence against normality."
  #       else "Evidence against normality (sensitive in tails)."
  #     } else "Package 'nortest' missing or test unavailable."
  #   )
  #   
  #   tests_df <- if (length(test_rows)) do.call(rbind, test_rows) else data.frame()
  #   
  #   # ---------- forecast accuracy (safe)
  #   acc_df <- NULL
  #   acc_sentence <- "No holdout test set was detected; therefore, out-of-sample accuracy was not computed."
  #   
  #   y_test <- s$ts_test
  #   has_test <- !is.null(y_test) && safe_len(y_test) > 0 && n_test > 0
  #   
  #   if (has_test) {
  #     y_test_num <- as.numeric(y_test)
  #     y_hat_num  <- suppressWarnings(tryCatch(as.numeric(fc$mean), error = function(e) rep(NA_real_, length(y_test_num))))
  #     h <- min(length(y_test_num), length(y_hat_num))
  #     if (h >= 1) {
  #       e <- y_test_num[seq_len(h)] - y_hat_num[seq_len(h)]
  #       rmse <- sqrt(mean(e^2, na.rm = TRUE))
  #       mae  <- mean(abs(e), na.rm = TRUE)
  #       mape <- mean(abs(e) / pmax(abs(y_test_num[seq_len(h)]), .Machine$double.eps), na.rm = TRUE)
  #       
  #       acc_df <- data.frame(
  #         Metric = c("RMSE", "MAE", "MAPE"),
  #         Value  = c(fmt_num_local(rmse, 3), fmt_num_local(mae, 3), paste0(fmt_num_local(100*mape, 2), "%")),
  #         check.names = FALSE
  #       )
  #       
  #       acc_sentence <- paste0(
  #         "Over the holdout period (test n = ", h, "), forecast performance was ",
  #         "RMSE = ", fmt_num_local(rmse, 3), ", ",
  #         "MAE = ", fmt_num_local(mae, 3), ", ",
  #         "MAPE = ", fmt_num_local(100*mape, 2), "%."
  #       )
  #     }
  #   }
  #   
  #   
  #   
  #   # =========================
  #   # Manual report plots (used inside manual_conclusion_full_obj)
  #   # =========================
  #   
  #   # helper: apply differencing D (seasonal) then d (non-seasonal)
  #   apply_sarima_diffs <- function(y, d = 0L, D = 0L, s = 1L) {
  #     y <- as.numeric(y)
  #     y <- y[is.finite(y)]
  #     if (length(y) < 3) return(y)
  #     
  #     d <- as.integer(d); if (!is.finite(d) || d < 0) d <- 0L
  #     D <- as.integer(D); if (!is.finite(D) || D < 0) D <- 0L
  #     s <- as.integer(s); if (!is.finite(s) || s < 1) s <- 1L
  #     
  #     yy <- y
  #     
  #     # seasonal differences first (common SARIMA workflow)
  #     if (D > 0 && s > 1) {
  #       for (i in seq_len(D)) {
  #         if (length(yy) <= s + 1) break
  #         yy <- diff(yy, lag = s)
  #       }
  #     }
  #     
  #     # then non-seasonal differences
  #     if (d > 0) {
  #       for (i in seq_len(d)) {
  #         if (length(yy) <= 2) break
  #         yy <- diff(yy, lag = 1)
  #       }
  #     }
  #     
  #     yy
  #   }
  #   
  #   # ---- A) Time series plot (dates + dashed split) ----
  #   output$manual_report_ts_plot <- renderPlot({
  #     req(ts_train_test(), prepared())
  #     s <- ts_train_test()
  #     p <- prepared()
  #     
  #     df <- s$dfm
  #     validate(need(nrow(df) >= 3, "Not enough observations to plot the time series."))
  #     
  #     df$set <- ifelse(seq_len(nrow(df)) <= s$train_n, "Train", "Test/Future")
  #     
  #     # split x-position (draw at boundary between train and test)
  #     has_test <- isTRUE(s$test_n > 0)
  #     x_split  <- if (has_test) df$x[s$train_n] else NA
  #     
  #     g <- ggplot(df, aes(x = x, y = y_trans, color = set)) +
  #       geom_line(linewidth = 0.9) +
  #       theme_minimal() +
  #       labs(
  #         title = "Observed time series (transformed)",
  #         x = p$x_label,
  #         y = "Value",
  #         color = NULL
  #       ) +
  #       theme(legend.position = "bottom")
  #     
  #     # dashed split line
  #     if (has_test && !is.na(x_split)) {
  #       g <- g + geom_vline(xintercept = as.numeric(x_split), linetype = "dashed", linewidth = 0.7, color = "gray40")
  #     }
  #     
  #     # Date-safe axis
  #     if (inherits(df$x, "Date")) {
  #       g <- g + scale_x_date(labels = scales::date_format("%Y-%m"), breaks = scales::pretty_breaks(n = 8))
  #     } else if (inherits(df$x, "POSIXt")) {
  #       g <- g + scale_x_datetime(labels = scales::date_format("%Y-%m"), breaks = scales::pretty_breaks(n = 8))
  #     }
  #     
  #     g
  #   })
  #   
  #   # ---- B) ACF / PACF of training series ----
  #   output$manual_report_acf <- renderPlot({
  #     req(ts_train_test())
  #     x <- ts_train_test()$ts_train
  #     validate(need(length(x) >= 5, "Not enough training observations for ACF."))
  #     forecast::ggAcf(x, lag.max = min(60, length(x) - 1)) +
  #       theme_minimal() +
  #       labs(title = "ACF (training)")
  #   })
  #   
  #   output$manual_report_pacf <- renderPlot({
  #     req(ts_train_test())
  #     x <- ts_train_test()$ts_train
  #     validate(need(length(x) >= 5, "Not enough training observations for PACF."))
  #     forecast::ggPacf(x, lag.max = min(60, length(x) - 1)) +
  #       theme_minimal() +
  #       labs(title = "PACF (training)")
  #   })
  #   
  #   # ---- C) ACF / PACF of modified (differenced) series based on current d, D, s ----
  #   output$manual_report_acf_mod <- renderPlot({
  #     req(ts_train_test(), prepared())
  #     s_obj <- ts_train_test()
  #     p <- prepared()
  #     
  #     # use training series for identification (typical)
  #     y <- as.numeric(s_obj$ts_train)
  #     
  #     s_use <- if (is.na(input$s)) p$freq else as.integer(input$s)
  #     y_mod <- apply_sarima_diffs(y, d = input$d, D = input$D, s = s_use)
  #     
  #     validate(need(length(y_mod) >= 5, "Differencing left too few observations for ACF."))
  #     
  #     ts_mod <- ts(y_mod, frequency = p$freq)
  #     
  #     forecast::ggAcf(ts_mod, lag.max = min(60, length(ts_mod) - 1)) +
  #       theme_minimal() +
  #       labs(title = paste0("ACF (d=", input$d, ", D=", input$D, ", s=", s_use, ")"))
  #   })
  #   
  #   output$manual_report_pacf_mod <- renderPlot({
  #     req(ts_train_test(), prepared())
  #     s_obj <- ts_train_test()
  #     p <- prepared()
  #     
  #     y <- as.numeric(s_obj$ts_train)
  #     
  #     s_use <- if (is.na(input$s)) p$freq else as.integer(input$s)
  #     y_mod <- apply_sarima_diffs(y, d = input$d, D = input$D, s = s_use)
  #     
  #     validate(need(length(y_mod) >= 5, "Differencing left too few observations for PACF."))
  #     
  #     ts_mod <- ts(y_mod, frequency = p$freq)
  #     
  #     forecast::ggPacf(ts_mod, lag.max = min(60, length(ts_mod) - 1)) +
  #       theme_minimal() +
  #       labs(title = paste0("PACF (d=", input$d, ", D=", input$D, ", s=", s_use, ")"))
  #   })
  #   
  #   output$manual_report_ts_trans_and_diff <- renderPlot({
  #     req(ts_train_test(), prepared())
  #     s_obj <- ts_train_test()
  #     p <- prepared()
  #     
  #     df <- s_obj$dfm
  #     validate(need(nrow(df) >= 5, "Not enough observations to plot."))
  #     
  #     # Use ONLY training portion if split exists
  #     train_n <- s_obj$train_n
  #     has_test <- isTRUE(s_obj$test_n > 0)
  #     df_train <- if (has_test) df[seq_len(train_n), , drop = FALSE] else df
  #     
  #     validate(need(nrow(df_train) >= 5, "Not enough training observations to plot."))
  #     
  #     # Base transformed training series
  #     x_train <- df_train$x
  #     y_train <- df_train$y_trans
  #     
  #     # Differenced series (based on current d, D, s)
  #     # Use p$freq if input$s is NA
  #     s_use <- if (is.null(input$s) || is.na(input$s)) p$freq else as.integer(input$s)
  #     
  #     y_mod <- apply_sarima_diffs(y_train, d = input$d, D = input$D, s = s_use)
  #     validate(need(length(y_mod) >= 3, "Differencing left too few observations to plot."))
  #     
  #     # Align differenced series to the LAST dates of training (since diff shortens length)
  #     x_mod <- tail(x_train, length(y_mod))
  #     
  #     plot_df <- rbind(
  #       data.frame(x = x_train, y = y_train, series = "Transformed (train)"),
  #       data.frame(x = x_mod,   y = y_mod,   series = paste0("Differenced (d=", input$d, ", D=", input$D, ", s=", s_use, ")"))
  #     )
  #     
  #     g <- ggplot(plot_df, aes(x = x, y = y, color = series)) +
  #       geom_line(linewidth = 0.9) +
  #       theme_minimal() +
  #       labs(
  #         title = "Training series: transformed vs differenced",
  #         x = p$x_label,
  #         y = "Value",
  #         color = NULL
  #       ) +
  #       theme(legend.position = "bottom")
  #     
  #     # Date-safe axis
  #     if (inherits(plot_df$x, "Date")) {
  #       g <- g + scale_x_date(labels = scales::date_format("%Y-%m"), breaks = scales::pretty_breaks(n = 8))
  #     } else if (inherits(plot_df$x, "POSIXt")) {
  #       g <- g + scale_x_datetime(labels = scales::date_format("%Y-%m"), breaks = scales::pretty_breaks(n = 8))
  #     }
  #     
  #     g
  #   })
  #   
  #   
  #   
  #   
  #   
  #   
  #   # ---------- horizon narrative (as in your PDF style)
  #   horizon_txt <- if (has_test) {
  #     paste0("Validation mode was used: the forecast horizon was forced to match the test length (h = ", n_test, ").")
  #   } else {
  #     paste0("Future mode was used: forecasts were produced beyond the training sample (h = ", fc0$h, ").")
  #   }
  #   
  #   # ---------- model string (safe)
  #   season_txt <- suppressWarnings(as.integer(eq$s))
  #   s_txt <- if (is.finite(season_txt) && season_txt > 0) as.character(season_txt) else "s"
  #   
  #   model_str <- sprintf(
  #     "SARIMA(%d,%d,%d)(%d,%d,%d)[%s]",
  #     eq$p, eq$d, eq$q, eq$P, eq$D, eq$Q, s_txt
  #   )
  #   
  #   # ---------- final decision bullets (lightweight logic)
  #   lb_ok <- !is.null(lb) && is.finite(lb$p.value) && lb$p.value >= 0.05
  #   arch_ok <- is.null(arch) || (is.finite(arch$p.value) && arch$p.value >= 0.05)  # if no test, don't overstate
  #   diag_verdict <- paste0(
  #     if (lb_ok) "Residual autocorrelation was not statistically detected (Ljung–Box p ≥ .05). "
  #     else "Residual autocorrelation may remain (Ljung–Box p < .05). ",
  #     if (arch_ok) "No clear evidence of residual ARCH effects was found (or test unavailable)."
  #     else "Residual ARCH effects were detected → a GARCH extension is recommended."
  #   )
  #   
  #   # ---------- build UI (includes your existing plots/tables)
  #   tagList(
  #     tags$h3("Manual SARIMA: Full academic conclusion (report-ready)"),
  #     
  #     tags$h4("1. Objective and modelling rationale"),
  #     tags$p(
  #       "A manually specified seasonal ARIMA (SARIMA) model was estimated to represent linear temporal dependence,",
  #       " including seasonal structure, and to provide an interpretable baseline for forecasting."
  #     ),
  #     
  #     tags$h4("2. Data design and sample split"),
  #     tags$p(HTML(paste0(
  #       "The analysis used <b>N = ", N, "</b> observations (training <b>n = ", n_train, "</b>",
  #       if (has_test) paste0(", test <b>n = ", n_test, "</b>") else "",
  #       ")."
  #     ))),
  #     tags$p(tags$b("Forecast design. "), horizon_txt),
  #     
  #     
  #     
  #     # --- NEW: Identification visuals for the report ---
  #     tags$h4("3. Identification visuals (time series + ACF/PACF)"),
  #     tags$p(
  #       "The plots below summarize the observed series (with the train/test split, if applicable), ",
  #       "followed by ACF/PACF for the training series and for the differenced series implied by the chosen (d, D, s)."
  #     ),
  #     
  #     tags$h5("Figure A1. Time series with split marker"),
  #     plotOutput("manual_report_ts_plot", height = 360),
  #     
  #     tags$h5("Figure A2. Transformed training series and differenced (d, D, s) series"),
  #     plotOutput("manual_report_ts_trans_and_diff", height = 360),
  #     
  #     tags$h5("Figure B. ACF and PACF (training series)"),
  #     fluidRow(
  #       column(6, plotOutput("manual_report_acf",  height = 280)),
  #       column(6, plotOutput("manual_report_pacf", height = 280))
  #     ),
  #     
  #     tags$h5("Figure C. ACF and PACF (modified / differenced series using current d, D, s)"),
  #     fluidRow(
  #       column(6, plotOutput("manual_report_acf_mod",  height = 280)),
  #       column(6, plotOutput("manual_report_pacf_mod", height = 280))
  #     ),
  #     
  #     
  #     
  #     
  #     
  #     tags$h4("4. Final model specification and fit quality"),
  #     tags$p(HTML(paste0(
  #       "The final manual specification was <b>", model_str, "</b>",
  #       if (isTRUE(input$manual_drift)) " including drift/mean." else " without drift/mean.",
  #       " Model adequacy was assessed using information criteria, coefficient inference, residual diagnostics, and forecast performance."
  #     ))),
  #     
  #     tags$h5("Table A. Goodness-of-fit (information criteria)"),
  #     html_table(ic_df),
  #     
  #     tags$h5("Table B. Parameter estimates and significance (approx. z/t tests)"),
  #     if (!is.null(coef_df)) html_table(coef_df) else tags$em("No coefficients available."),
  #     tags$p(HTML(paste0(
  #       "In total, <b>", n_sig, "</b> parameter(s) were significant at α = .05 (marked by *, **, ***)."
  #     ))),
  #     
  #     tags$h4("5. Model equations (replication-ready)"),
  #     tags$p(
  #       "The fitted model is reported below in operator notation (general form), expanded form, and the numerical equation",
  #       " based on the estimated parameters."
  #     ),
  #     tags$div(
  #       style = "padding:10px;border:1px solid #e5e5e5;border-radius:6px;background:#fcfcfc;",
  #       tags$h5("General SARIMA formulation"),
  #       HTML(eq$eq_general),
  #       tags$hr(),
  #       tags$h5("Expanded operator form"),
  #       HTML(eq$eq_expanded),
  #       tags$hr(),
  #       tags$h5("Numerical model"),
  #       HTML(eq$eq_line3),
  #       tags$hr(),
  #       HTML(eq$eq_line4)
  #     ),
  #     
  #     tags$h4("6. Residual diagnostics (graphical evidence)"),
  #     tags$p(
  #       "Graphical diagnostics evaluate whether residuals resemble white noise (no systematic autocorrelation),",
  #       " approximate normality (Q–Q and histogram), and stable variance."
  #     ),
  #     
  #     fluidRow(
  #       column(6, plotOutput("manual_resid_ts",   height = 220)),
  #       column(6, plotOutput("manual_resid_acf",  height = 220))
  #     ),
  #     fluidRow(
  #       column(6, plotOutput("manual_resid_hist", height = 220)),
  #       column(6, plotOutput("manual_resid_qq",   height = 220))
  #     ),
  #     tags$h5("Ljung–Box p-values by lag"),
  #     plotOutput("manual_resid_lb_pvals", height = 260),
  #     
  #     tags$h4("7. Residual tests (formal inference)"),
  #     tags$p(
  #       "Formal tests complement the plots: Ljung–Box/Box–Pierce assess remaining autocorrelation;",
  #       " Jarque–Bera/Shapiro–Wilk/Anderson–Darling assess normality;",
  #       " ARCH LM tests conditional heteroskedasticity; the runs test checks randomness."
  #     ),
  #     tags$h5("Table C. Residual test summary"),
  #     html_table(tests_df),
  #     tags$p(tags$b("Diagnostic synthesis. "), diag_verdict),
  #     
  #     tags$h4("8. Forecasting results and predictive performance"),
  #     tags$p(acc_sentence),
  #     
  #     if (!is.null(acc_df)) tagList(
  #       tags$h5("Table D. Holdout accuracy (test set)"),
  #       html_table(acc_df)
  #     ) else NULL,
  #     
  #     tags$h5("Forecast plot"),
  #     plotOutput("manual_forecast_plot", height = 420),
  #     
  #     tags$h5("Forecast table"),
  #     tableOutput("manual_forecast_table"),
  #     
  #     tags$h5("Accuracy table (your app output)"),
  #     tableOutput("manual_accuracy_table"),
  #     
  #     tags$h4("9. Final conclusion (academic)"),
  #     tags$p(
  #       "Overall, the manually specified SARIMA model provides a coherent and interpretable representation of seasonal linear dynamics,",
  #       " supported by information criteria, statistically interpretable parameters, and diagnostic checks.",
  #       " When diagnostics indicate remaining autocorrelation, refinement should prioritize revising differencing (d, D) and AR/MA orders guided by ACF/PACF and Ljung–Box.",
  #       " When conditional heteroskedasticity is detected (ARCH LM), a volatility model (e.g., GARCH) should be added to the mean equation to better represent time-varying variance."
  #     ),
  #     tags$p(
  #       "For reporting, the results above provide the full chain of evidence typically expected in academic manuscripts:",
  #       " (i) specification + IC, (ii) parameter inference, (iii) equation reporting, (iv) residual validation with plots and tests, and (v) forecasting with accuracy assessment."
  #     )
  #   )
  # })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # ============================================================
  # CORRECTED STRUCTURE
  #   1) Define ALL report outputs ONCE (outside eventReactive)
  #   2) manual_conclusion_full_obj ONLY computes + returns tagList UI
  # ============================================================
  
  
  # ------------------------------------------------------------
  # (A) REPORT OUTPUTS (define ONCE, outside manual_conclusion_full_obj)
  # ------------------------------------------------------------
  
  # ---- A) Time series plot (dates + dashed split) ----
  output$manual_report_ts_plot <- renderPlot({
    req(manual_conclusion_full_obj())     # ensures this is shown after Fit Manual
    req(ts_train_test(), prepared())
    
    s <- ts_train_test()
    p <- prepared()
    
    df <- s$dfm
    validate(need(nrow(df) >= 3, "Not enough observations to plot the time series."))
    
    df$set <- ifelse(seq_len(nrow(df)) <= s$train_n, "Train", "Test/Future")
    
    has_test <- isTRUE(s$test_n > 0)
    x_split  <- if (has_test) df$x[s$train_n] else NA
    
    g <- ggplot(df, aes(x = x, y = y_trans, color = set)) +
      geom_line(linewidth = 0.9) +
      theme_minimal() +
      labs(
        title = "Observed time series (transformed)",
        x = p$x_label,
        y = "Value",
        color = NULL
      ) +
      theme(legend.position = "bottom")
    
    if (has_test && !is.na(x_split)) {
      g <- g + geom_vline(
        xintercept = as.numeric(x_split),
        linetype = "dashed",
        linewidth = 0.7,
        color = "gray40"
      )
    }
    
    if (inherits(df$x, "Date")) {
      g <- g + scale_x_date(labels = scales::date_format("%Y-%m"),
                            breaks = scales::pretty_breaks(n = 8))
    } else if (inherits(df$x, "POSIXt")) {
      g <- g + scale_x_datetime(labels = scales::date_format("%Y-%m"),
                                breaks = scales::pretty_breaks(n = 8))
    }
    
    g
  })
  
  
  # ---- NEW: Stationarity tests (ADF + PP + KPSS) + conclusion paragraph ----
  output$manual_report_stationarity <- renderUI({
    req(manual_conclusion_full_obj())
    req(ts_train_test())
    
    validate(need(requireNamespace("tseries", quietly = TRUE),
                  "Package 'tseries' is required for stationarity tests (ADF/PP/KPSS)."))
    
    s_obj <- ts_train_test()
    
    # Use training series for stationarity assessment
    x <- as.numeric(s_obj$ts_train)
    x <- x[is.finite(x)]
    validate(need(length(x) >= 10, "Not enough training observations to run stationarity tests (need ≥ 10)."))
    
    fmt_num <- function(z, d = 3) {
      z <- suppressWarnings(as.numeric(z))
      if (length(z) == 0 || !is.finite(z[1])) return("NA")
      formatC(z[1], format = "f", digits = d)
    }
    fmt_p <- function(p) {
      p <- suppressWarnings(as.numeric(p))
      if (length(p) == 0 || !is.finite(p[1])) return("NA")
      if (p[1] < .001) "p < .001" else paste0("p = ", sub("^0\\.", ".", sprintf("%.3f", p[1])))
    }
    
    # safe lag choice for ADF
    k_adf <- max(0, min(12, floor((length(x) - 1)^(1/3))))
    
    adf <- tryCatch(tseries::adf.test(x, k = k_adf), error = function(e) NULL)
    pp  <- tryCatch(tseries::pp.test(x, lshort = TRUE), error = function(e) NULL)
    kpss_level <- tryCatch(tseries::kpss.test(x, null = "Level", lshort = TRUE), error = function(e) NULL)
    kpss_trend <- tryCatch(tseries::kpss.test(x, null = "Trend", lshort = TRUE), error = function(e) NULL)
    
    rows <- list()
    add_row <- function(test, null_h, stat, pval, decision) {
      rows[[length(rows) + 1L]] <<- data.frame(
        Test = test,
        `H0 (null)` = null_h,
        Statistic = stat,
        `p-value` = pval,
        Decision = decision,
        check.names = FALSE
      )
    }
    
    # ADF: H0 = unit root (non-stationary); reject => stationarity evidence
    if (!is.null(adf)) {
      p <- adf$p.value
      dec <- if (is.finite(p) && p < 0.05) "Reject H0 → evidence for stationarity"
      else "Fail to reject H0 → unit root plausible"
      add_row("ADF (Augmented Dickey–Fuller)", "Unit root (non-stationary)",
              fmt_num(unname(adf$statistic)), fmt_p(p), dec)
    } else {
      add_row("ADF (Augmented Dickey–Fuller)", "Unit root (non-stationary)", "NA", "NA", "Not available")
    }
    
    # PP: H0 = unit root (non-stationary); reject => stationarity evidence
    if (!is.null(pp)) {
      p <- pp$p.value
      dec <- if (is.finite(p) && p < 0.05) "Reject H0 → evidence for stationarity"
      else "Fail to reject H0 → unit root plausible"
      add_row("PP (Phillips–Perron)", "Unit root (non-stationary)",
              fmt_num(unname(pp$statistic)), fmt_p(p), dec)
    } else {
      add_row("PP (Phillips–Perron)", "Unit root (non-stationary)", "NA", "NA", "Not available")
    }
    
    # KPSS(Level): H0 = level-stationary; reject => non-stationarity evidence
    if (!is.null(kpss_level)) {
      p <- kpss_level$p.value
      dec <- if (is.finite(p) && p < 0.05) "Reject H0 → evidence against stationarity"
      else "Fail to reject H0 → stationarity plausible"
      add_row("KPSS (Level)", "Level-stationary",
              fmt_num(unname(kpss_level$statistic)), fmt_p(p), dec)
    } else {
      add_row("KPSS (Level)", "Level-stationary", "NA", "NA", "Not available")
    }
    
    # KPSS(Trend): H0 = trend-stationary
    if (!is.null(kpss_trend)) {
      p <- kpss_trend$p.value
      dec <- if (is.finite(p) && p < 0.05) "Reject H0 → evidence against trend-stationarity"
      else "Fail to reject H0 → trend-stationarity plausible"
      add_row("KPSS (Trend)", "Trend-stationary",
              fmt_num(unname(kpss_trend$statistic)), fmt_p(p), dec)
    } else {
      add_row("KPSS (Trend)", "Trend-stationary", "NA", "NA", "Not available")
    }
    
    st_df <- do.call(rbind, rows)
    
    # synthesis
    adf_p <- if (!is.null(adf)) adf$p.value else NA_real_
    pp_p  <- if (!is.null(pp))  pp$p.value  else NA_real_
    kL_p  <- if (!is.null(kpss_level)) kpss_level$p.value else NA_real_
    
    unit_root_rejected <- any(c(adf_p, pp_p) < 0.05, na.rm = TRUE)
    kpss_ok <- is.finite(kL_p) && kL_p >= 0.05
    
    conclusion <- if (unit_root_rejected && kpss_ok) {
      "Across tests, ADF/PP reject the unit-root null (p < .05) while KPSS(Level) does not reject stationarity (p ≥ .05), which is consistent with a stationary series (given the current transformation)."
    } else if (!unit_root_rejected && !kpss_ok) {
      "ADF/PP do not reject the unit-root null (p ≥ .05) while KPSS(Level) rejects stationarity (p < .05), providing convergent evidence of non-stationarity and supporting the need for differencing (d and/or D)."
    } else if (unit_root_rejected && !kpss_ok) {
      "Evidence is mixed: ADF/PP suggest stationarity but KPSS(Level) rejects it. This can occur under breaks, strong seasonality, or test sensitivity; complement these results with differencing checks and ACF/PACF."
    } else {
      "Evidence is inconclusive: ADF/PP do not reject a unit root while KPSS(Level) does not reject stationarity. Because power can be limited, complement these tests with differencing diagnostics and ACF/PACF."
    }
    
    # table renderer (HTML)
    html_tbl <- tags$table(
      class = "table table-striped table-condensed",
      tags$thead(tags$tr(lapply(names(st_df), tags$th))),
      tags$tbody(lapply(seq_len(nrow(st_df)), function(i) {
        tags$tr(lapply(st_df[i, , drop = FALSE], function(cell) tags$td(HTML(as.character(cell)))))
      }))
    )
    
    tagList(
      # tags$h4(tags$strong("3. Stationarity assessment (ADF, KPSS, and Phillips–Perron)")),
      # tags$p("Stationarity tests were applied to the training series to evaluate whether differencing is required before SARIMA identification and estimation."),
      html_tbl,
      tags$p(tags$b("Conclusion. "), conclusion)
    )
  })
  
  
  # ---- helper used by multiple outputs (define ONCE) ----
  apply_sarima_diffs_report <- function(y, d = 0L, D = 0L, s = 1L) {
    y <- as.numeric(y)
    y <- y[is.finite(y)]
    if (length(y) < 3) return(y)
    
    d <- as.integer(d); if (!is.finite(d) || d < 0) d <- 0L
    D <- as.integer(D); if (!is.finite(D) || D < 0) D <- 0L
    s <- as.integer(s); if (!is.finite(s) || s < 1) s <- 1L
    
    yy <- y
    
    if (D > 0 && s > 1) {
      for (i in seq_len(D)) {
        if (length(yy) <= s + 1) break
        yy <- diff(yy, lag = s)
      }
    }
    if (d > 0) {
      for (i in seq_len(d)) {
        if (length(yy) <= 2) break
        yy <- diff(yy, lag = 1)
      }
    }
    yy
  }
  
  
  # ---- B) Transformed training series vs differenced (d,D,s) ----
  output$manual_report_ts_trans_and_diff <- renderPlot({
    req(manual_conclusion_full_obj())
    req(ts_train_test(), prepared())
    
    s_obj <- ts_train_test()
    p <- prepared()
    
    df <- s_obj$dfm
    validate(need(nrow(df) >= 5, "Not enough observations to plot."))
    
    train_n <- s_obj$train_n
    has_test <- isTRUE(s_obj$test_n > 0)
    df_train <- if (has_test) df[seq_len(train_n), , drop = FALSE] else df
    
    validate(need(nrow(df_train) >= 5, "Not enough training observations to plot."))
    
    x_train <- df_train$x
    y_train <- df_train$y_trans
    
    s_use <- if (is.null(input$s) || is.na(input$s)) p$freq else as.integer(input$s)
    y_mod <- apply_sarima_diffs_report(y_train, d = input$d, D = input$D, s = s_use)
    validate(need(length(y_mod) >= 3, "Differencing left too few observations to plot."))
    
    x_mod <- tail(x_train, length(y_mod))
    
    plot_df <- rbind(
      data.frame(x = x_train, y = y_train, series = "Transformed (train)"),
      data.frame(x = x_mod,   y = y_mod,   series = paste0("Differenced (d=", input$d, ", D=", input$D, ", s=", s_use, ")"))
    )
    
    g <- ggplot(plot_df, aes(x = x, y = y, color = series)) +
      geom_line(linewidth = 0.9) +
      theme_minimal() +
      labs(
        title = "Training series: transformed vs differenced",
        x = p$x_label,
        y = "Value",
        color = NULL
      ) +
      theme(legend.position = "bottom")
    
    if (inherits(plot_df$x, "Date")) {
      g <- g + scale_x_date(labels = scales::date_format("%Y-%m"),
                            breaks = scales::pretty_breaks(n = 8))
    } else if (inherits(plot_df$x, "POSIXt")) {
      g <- g + scale_x_datetime(labels = scales::date_format("%Y-%m"),
                                breaks = scales::pretty_breaks(n = 8))
    }
    
    g
  })
  
  
  # ---- C) Seasonal subseries (full observed) ----
  # output$manual_report_subseries <- renderPlot({
  #   req(manual_conclusion_full_obj())
  #   req(ts_train_test())
  #   
  #   s_obj <- ts_train_test()
  #   
  #   x_full <- ts(
  #     c(as.numeric(s_obj$ts_train),
  #       if (!is.null(s_obj$ts_test) && length(s_obj$ts_test) > 0) as.numeric(s_obj$ts_test) else numeric(0)),
  #     start = 1,
  #     frequency = frequency(s_obj$ts_train)
  #   )
  #   
  #   validate(need(frequency(x_full) >= 2, "Seasonal subseries plot requires seasonal frequency (s) >= 2."))
  #   validate(need(length(x_full) >= 2 * frequency(x_full), "Need at least 2 seasonal cycles for a subseries plot."))
  #   
  #   forecast::ggsubseriesplot(x_full) +
  #     theme_minimal() +
  #     labs(title = "Seasonal subseries (observed series)", x = "Seasonal period", y = "Value")
  # })
  
  
  # ---- D) Seasonal box-plot (full observed) ----
  # output$manual_report_seasonal_box <- renderPlot({
  #   req(manual_conclusion_full_obj())
  #   req(ts_train_test())
  #   
  #   s_obj <- ts_train_test()
  #   
  #   x_full <- ts(
  #     c(as.numeric(s_obj$ts_train),
  #       if (!is.null(s_obj$ts_test) && length(s_obj$ts_test) > 0) as.numeric(s_obj$ts_test) else numeric(0)),
  #     start = 1,
  #     frequency = frequency(s_obj$ts_train)
  #   )
  #   
  #   validate(need(frequency(x_full) >= 2, "Seasonal box-plot requires seasonal frequency (s) >= 2."))
  #   validate(need(length(x_full) >= frequency(x_full), "Need at least 1 seasonal cycle for a box-plot."))
  #   
  #   df <- data.frame(
  #     value  = as.numeric(x_full),
  #     season = factor(cycle(x_full), ordered = TRUE)
  #   )
  #   df <- df[is.finite(df$value), , drop = FALSE]
  #   validate(need(nrow(df) >= 5, "Not enough valid observations for seasonal box-plot."))
  #   
  #   ggplot(df, aes(x = season, y = value)) +
  #     geom_boxplot(fill = "#2C7FB8", alpha = 0.45, outlier.alpha = 0.4) +
  #     theme_minimal() +
  #     labs(title = "Seasonal box-plot (observed series)", x = "Seasonal period", y = "Value")
  # })
  
  
  # ---- C) Seasonal subseries (TRAINING ONLY) ----
  output$manual_report_subseries <- renderPlot({
    req(manual_conclusion_full_obj())
    req(ts_train_test())
    
    s_obj <- ts_train_test()
    
    x_train <- s_obj$ts_train
    validate(need(inherits(x_train, "ts"), "Training series is not a 'ts' object."))
    validate(need(stats::frequency(x_train) >= 2, "Seasonal subseries plot requires seasonal frequency (s) >= 2."))
    validate(need(length(x_train) >= 2 * stats::frequency(x_train),
                  "Need at least 2 seasonal cycles in the TRAINING set for a subseries plot."))
    
    forecast::ggsubseriesplot(x_train) +
      theme_minimal() +
      labs(title = "Seasonal subseries (training series)", x = "Seasonal period", y = "Value")
  })
  
  
  # ---- D) Seasonal box-plot (TRAINING ONLY) ----
  output$manual_report_seasonal_box <- renderPlot({
    req(manual_conclusion_full_obj())
    req(ts_train_test())
    
    s_obj <- ts_train_test()
    
    x_train <- s_obj$ts_train
    validate(need(inherits(x_train, "ts"), "Training series is not a 'ts' object."))
    validate(need(stats::frequency(x_train) >= 2, "Seasonal box-plot requires seasonal frequency (s) >= 2."))
    validate(need(length(x_train) >= stats::frequency(x_train),
                  "Need at least 1 seasonal cycle in the TRAINING set for a box-plot."))
    
    df <- data.frame(
      value  = as.numeric(x_train),
      season = factor(stats::cycle(x_train), ordered = TRUE)
    )
    df <- df[is.finite(df$value), , drop = FALSE]
    validate(need(nrow(df) >= 5, "Not enough valid training observations for seasonal box-plot."))
    
    ggplot(df, aes(x = season, y = value)) +
      geom_boxplot(fill = "#2C7FB8", alpha = 0.45, outlier.alpha = 0.4) +
      theme_minimal() +
      labs(title = "Seasonal box-plot (training series)", x = "Seasonal period", y = "Value")
  })
  
  
  # ---- E) ACF / PACF of training series ----
  output$manual_report_acf <- renderPlot({
    req(manual_conclusion_full_obj())
    req(ts_train_test())
    
    x <- ts_train_test()$ts_train
    validate(need(length(x) >= 5, "Not enough training observations for ACF."))
    
    forecast::ggAcf(x, lag.max = min(60, length(x) - 1)) +
      theme_minimal() +
      labs(title = "ACF (training)")
  })
  
  output$manual_report_pacf <- renderPlot({
    req(manual_conclusion_full_obj())
    req(ts_train_test())
    
    x <- ts_train_test()$ts_train
    validate(need(length(x) >= 5, "Not enough training observations for PACF."))
    
    forecast::ggPacf(x, lag.max = min(60, length(x) - 1)) +
      theme_minimal() +
      labs(title = "PACF (training)")
  })
  
  
  # ---- F) ACF / PACF of differenced series ----
  output$manual_report_acf_mod <- renderPlot({
    req(manual_conclusion_full_obj())
    req(ts_train_test(), prepared())
    
    s_obj <- ts_train_test()
    p <- prepared()
    
    y <- as.numeric(s_obj$ts_train)
    s_use <- if (is.null(input$s) || is.na(input$s)) p$freq else as.integer(input$s)
    y_mod <- apply_sarima_diffs_report(y, d = input$d, D = input$D, s = s_use)
    
    validate(need(length(y_mod) >= 5, "Differencing left too few observations for ACF."))
    
    ts_mod <- ts(y_mod, frequency = p$freq)
    
    forecast::ggAcf(ts_mod, lag.max = min(60, length(ts_mod) - 1)) +
      theme_minimal() +
      labs(title = paste0("ACF (d=", input$d, ", D=", input$D, ", s=", s_use, ")"))
  })
  
  output$manual_report_pacf_mod <- renderPlot({
    req(manual_conclusion_full_obj())
    req(ts_train_test(), prepared())
    
    s_obj <- ts_train_test()
    p <- prepared()
    
    y <- as.numeric(s_obj$ts_train)
    s_use <- if (is.null(input$s) || is.na(input$s)) p$freq else as.integer(input$s)
    y_mod <- apply_sarima_diffs_report(y, d = input$d, D = input$D, s = s_use)
    
    validate(need(length(y_mod) >= 5, "Differencing left too few observations for PACF."))
    
    ts_mod <- ts(y_mod, frequency = p$freq)
    
    forecast::ggPacf(ts_mod, lag.max = min(60, length(ts_mod) - 1)) +
      theme_minimal() +
      labs(title = paste0("PACF (d=", input$d, ", D=", input$D, ", s=", s_use, ")"))
  })
  
  
  
  # ---- NEW: Stationarity tests on differenced/transformed series (d, D, s) ----
  output$manual_report_stationarity_mod <- renderUI({
    req(manual_conclusion_full_obj())
    req(ts_train_test(), prepared())
    
    validate(need(requireNamespace("tseries", quietly = TRUE),
                  "Package 'tseries' is required for stationarity tests (ADF/PP/KPSS)."))
    
    s_obj <- ts_train_test()
    p <- prepared()
    
    # training series
    y <- as.numeric(s_obj$ts_train)
    y <- y[is.finite(y)]
    
    # apply differencing implied by current d, D, s
    s_use <- if (is.null(input$s) || is.na(input$s)) p$freq else as.integer(input$s)
    y_mod <- apply_sarima_diffs_report(y, d = input$d, D = input$D, s = s_use)
    y_mod <- as.numeric(y_mod)
    y_mod <- y_mod[is.finite(y_mod)]
    
    validate(need(length(y_mod) >= 10,
                  "Not enough observations after differencing (d, D, s) to run stationarity tests (need ≥ 10)."))
    
    fmt_num <- function(z, d = 3) {
      z <- suppressWarnings(as.numeric(z))
      if (length(z) == 0 || !is.finite(z[1])) return("NA")
      formatC(z[1], format = "f", digits = d)
    }
    fmt_p <- function(pv) {
      pv <- suppressWarnings(as.numeric(pv))
      if (length(pv) == 0 || !is.finite(pv[1])) return("NA")
      if (pv[1] < .001) "p < .001" else paste0("p = ", sub("^0\\.", ".", sprintf("%.3f", pv[1])))
    }
    
    # safe lag choice for ADF (based on effective sample)
    k_adf <- max(0, min(12, floor((length(y_mod) - 1)^(1/3))))
    
    adf <- tryCatch(tseries::adf.test(y_mod, k = k_adf), error = function(e) NULL)
    pp  <- tryCatch(tseries::pp.test(y_mod, lshort = TRUE), error = function(e) NULL)
    kpss_level <- tryCatch(tseries::kpss.test(y_mod, null = "Level", lshort = TRUE), error = function(e) NULL)
    kpss_trend <- tryCatch(tseries::kpss.test(y_mod, null = "Trend", lshort = TRUE), error = function(e) NULL)
    
    rows <- list()
    add_row <- function(test, null_h, stat, pval, decision) {
      rows[[length(rows) + 1L]] <<- data.frame(
        Test = test,
        `H0 (null)` = null_h,
        Statistic = stat,
        `p-value` = pval,
        Decision = decision,
        check.names = FALSE
      )
    }
    
    # ADF
    if (!is.null(adf)) {
      pv <- adf$p.value
      dec <- if (is.finite(pv) && pv < 0.05) "Reject H0 → evidence for stationarity"
      else "Fail to reject H0 → unit root plausible"
      add_row("ADF (Augmented Dickey–Fuller)", "Unit root (non-stationary)",
              fmt_num(unname(adf$statistic)), fmt_p(pv), dec)
    } else {
      add_row("ADF (Augmented Dickey–Fuller)", "Unit root (non-stationary)", "NA", "NA", "Not available")
    }
    
    # PP
    if (!is.null(pp)) {
      pv <- pp$p.value
      dec <- if (is.finite(pv) && pv < 0.05) "Reject H0 → evidence for stationarity"
      else "Fail to reject H0 → unit root plausible"
      add_row("PP (Phillips–Perron)", "Unit root (non-stationary)",
              fmt_num(unname(pp$statistic)), fmt_p(pv), dec)
    } else {
      add_row("PP (Phillips–Perron)", "Unit root (non-stationary)", "NA", "NA", "Not available")
    }
    
    # KPSS Level
    if (!is.null(kpss_level)) {
      pv <- kpss_level$p.value
      dec <- if (is.finite(pv) && pv < 0.05) "Reject H0 → evidence against stationarity"
      else "Fail to reject H0 → stationarity plausible"
      add_row("KPSS (Level)", "Level-stationary",
              fmt_num(unname(kpss_level$statistic)), fmt_p(pv), dec)
    } else {
      add_row("KPSS (Level)", "Level-stationary", "NA", "NA", "Not available")
    }
    
    # KPSS Trend
    if (!is.null(kpss_trend)) {
      pv <- kpss_trend$p.value
      dec <- if (is.finite(pv) && pv < 0.05) "Reject H0 → evidence against trend-stationarity"
      else "Fail to reject H0 → trend-stationarity plausible"
      add_row("KPSS (Trend)", "Trend-stationary",
              fmt_num(unname(kpss_trend$statistic)), fmt_p(pv), dec)
    } else {
      add_row("KPSS (Trend)", "Trend-stationary", "NA", "NA", "Not available")
    }
    
    st_df <- do.call(rbind, rows)
    
    # synthesis
    adf_p <- if (!is.null(adf)) adf$p.value else NA_real_
    pp_p  <- if (!is.null(pp))  pp$p.value  else NA_real_
    kL_p  <- if (!is.null(kpss_level)) kpss_level$p.value else NA_real_
    
    unit_root_rejected <- any(c(adf_p, pp_p) < 0.05, na.rm = TRUE)
    kpss_ok <- is.finite(kL_p) && kL_p >= 0.05
    
    conclusion <- if (unit_root_rejected && kpss_ok) {
      paste0(
        "After applying differencing (d=", input$d, ", D=", input$D, ", s=", s_use, "), ",
        "ADF/PP reject the unit-root null (p < .05) while KPSS(Level) does not reject stationarity (p ≥ .05). ",
        "This pattern is consistent with a stationary series after differencing."
      )
    } else if (!unit_root_rejected && !kpss_ok) {
      paste0(
        "After differencing (d=", input$d, ", D=", input$D, ", s=", s_use, "), ",
        "ADF/PP do not reject a unit root (p ≥ .05) and KPSS(Level) rejects stationarity (p < .05). ",
        "This suggests the series may still be non-stationary (consider revising d/D or checking breaks/seasonality)."
      )
    } else if (unit_root_rejected && !kpss_ok) {
      paste0(
        "After differencing (d=", input$d, ", D=", input$D, ", s=", s_use, "), evidence is mixed: ",
        "ADF/PP suggest stationarity but KPSS(Level) rejects it. This can happen with breaks, strong seasonal effects, ",
        "or finite-sample sensitivity; complement with diagnostics/visual checks."
      )
    } else {
      paste0(
        "After differencing (d=", input$d, ", D=", input$D, ", s=", s_use, "), evidence is inconclusive: ",
        "ADF/PP do not reject a unit root while KPSS(Level) does not reject stationarity. ",
        "Use ACF/PACF of the differenced series and consider alternative lag choices or structural breaks."
      )
    }
    
    html_tbl <- tags$table(
      class = "table table-striped table-condensed",
      tags$thead(tags$tr(lapply(names(st_df), tags$th))),
      tags$tbody(lapply(seq_len(nrow(st_df)), function(i) {
        tags$tr(lapply(st_df[i, , drop = FALSE], function(cell) tags$td(HTML(as.character(cell)))))
      }))
    )
    
    tagList(
      tags$h4(tags$strong(paste0(
        "Stationarity assessment after differencing (d=", input$d,
        ", D=", input$D, ", s=", s_use, ")"
      ))),
      tags$p("The same stationarity tests were re-applied to the training series after applying the current differencing settings to verify that the working series is stationary."),
      tags$br(),
      html_tbl,
      tags$p(tags$b("Conclusion. "), conclusion)
    )
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # # ------------------------------------------------------------
  # # (B) manual_conclusion_full_obj (eventReactive) — UI builder ONLY
  # # ------------------------------------------------------------
  # manual_conclusion_full_obj <- eventReactive(input$fit_manual, {
  #   req(manual_fit(), manual_fc(), manual_equations(), ts_train_test())
  #   
  #   fit <- manual_fit()
  #   fc0 <- manual_fc()
  #   fc  <- fc0$fc
  #   eq  <- manual_equations()
  #   s   <- ts_train_test()
  #   
  #   # ---------- helpers (local + safe)
  #   fmt_num_local <- function(x, d = 3) {
  #     if (length(x) == 0 || all(is.na(x))) return("NA")
  #     x <- suppressWarnings(as.numeric(x[1]))
  #     if (!is.finite(x)) return("NA")
  #     formatC(x, format = "f", digits = d)
  #   }
  #   fmt_p_local <- function(p) {
  #     if (length(p) == 0 || all(is.na(p))) return("NA")
  #     p <- suppressWarnings(as.numeric(p[1]))
  #     if (!is.finite(p)) return("NA")
  #     if (p < .001) "&lt; .001" else sprintf("= %.3f", p)
  #   }
  #   sig_stars <- function(p) {
  #     p <- suppressWarnings(as.numeric(p))
  #     if (!is.finite(p)) return("")
  #     if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else if (p < .10) "†" else ""
  #   }
  #   safe_len <- function(x) if (is.null(x)) 0L else length(x)
  #   
  #   html_table <- function(df) {
  #     if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
  #       return(tags$em("Table unavailable."))
  #     }
  #     tags$table(
  #       class = "table table-striped table-condensed",
  #       tags$thead(tags$tr(lapply(names(df), function(nm) tags$th(nm)))),
  #       tags$tbody(
  #         lapply(seq_len(nrow(df)), function(i) {
  #           tags$tr(lapply(df[i, , drop = FALSE], function(cell) tags$td(HTML(as.character(cell)))))
  #         })
  #       )
  #     )
  #   }
  #   
  #   # ---------- sample sizes (safe)
  #   n_train <- suppressWarnings(as.integer(s$train_n))
  #   if (!is.finite(n_train) || n_train < 1) n_train <- tryCatch(length(residuals(fit)), error = function(e) 0L)
  #   
  #   n_test <- suppressWarnings(as.integer(s$test_n))
  #   if (!is.finite(n_test) || n_test < 0) n_test <- 0L
  #   
  #   N <- n_train + n_test
  #   
  #   # ---------- lag choice (safe)
  #   L_in <- suppressWarnings(as.integer(input$diag_lag))
  #   L <- if (is.finite(L_in) && L_in > 0) L_in else 12L
  #   
  #   # ============================================================
  #   # ---------- IC (safe + robust AICc fallback)
  #   # ============================================================
  #   AIC_val <- suppressWarnings(tryCatch(as.numeric(stats::AIC(fit)), error = function(e) NA_real_))
  #   BIC_val <- suppressWarnings(tryCatch(as.numeric(stats::BIC(fit)), error = function(e) NA_real_))
  #   
  #   AICc_val <- suppressWarnings(tryCatch(as.numeric(forecast::AICc(fit)), error = function(e) NA_real_))
  #   
  #   if (!is.finite(AICc_val) && is.finite(AIC_val)) {
  #     n_fit <- suppressWarnings(tryCatch(stats::nobs(fit), error = function(e) NA_integer_))
  #     if (!is.finite(n_fit) || n_fit <= 0) {
  #       r <- suppressWarnings(tryCatch(as.numeric(residuals(fit)), error = function(e) numeric(0)))
  #       n_fit <- sum(is.finite(r))
  #     }
  #     if (!is.finite(n_fit) || n_fit <= 0) n_fit <- n_train
  #     
  #     k <- suppressWarnings(tryCatch(length(stats::coef(fit)), error = function(e) NA_integer_))
  #     if (!is.finite(k) || k < 0) k <- 0L
  #     k <- k + 1L
  #     
  #     if (is.finite(n_fit) && n_fit > (k + 1)) {
  #       AICc_val <- AIC_val + (2 * k * (k + 1)) / (n_fit - k - 1)
  #     } else {
  #       AICc_val <- NA_real_
  #     }
  #   }
  #   
  #   ic_df <- data.frame(
  #     Criterion = c("AIC", "AICc", "BIC"),
  #     Value     = c(fmt_num_local(AIC_val, 2), fmt_num_local(AICc_val, 2), fmt_num_local(BIC_val, 2)),
  #     check.names = FALSE
  #   )
  #   
  #   # ---------- coefficients + significance (robust)
  #   est <- suppressWarnings(tryCatch(stats::coef(fit), error = function(e) NULL))
  #   V   <- suppressWarnings(tryCatch(stats::vcov(fit), error = function(e) NULL))
  #   if (is.null(V)) V <- suppressWarnings(tryCatch(fit$var.coef, error = function(e) NULL))
  #   
  #   coef_df <- NULL
  #   if (!is.null(est) && length(est) > 0) {
  #     est <- as.numeric(est)
  #     nm  <- names(stats::coef(fit))
  #     if (is.null(nm)) nm <- paste0("param_", seq_along(est))
  #     
  #     se <- rep(NA_real_, length(est))
  #     if (!is.null(V)) {
  #       dV <- tryCatch(diag(V), error = function(e) rep(NA_real_, length(est)))
  #       if (length(dV) == length(est)) se <- sqrt(pmax(dV, 0))
  #     }
  #     z  <- est / se
  #     p  <- 2 * stats::pnorm(-abs(z))
  #     
  #     coef_df <- data.frame(
  #       Term     = nm,
  #       Estimate = sprintf("%.6f", est),
  #       SE       = ifelse(is.finite(se), sprintf("%.6f", se), "NA"),
  #       `z/t`    = ifelse(is.finite(z),  sprintf("%.3f",  z),  "NA"),
  #       `p`      = ifelse(is.finite(p),  sprintf("%.3f",  p),  "NA"),
  #       Sig      = vapply(p, sig_stars, character(1)),
  #       check.names = FALSE
  #     )
  #   }
  #   
  #   n_sig <- if (!is.null(coef_df)) sum(suppressWarnings(as.numeric(coef_df$p)) < 0.05, na.rm = TRUE) else 0L
  #   
  #   # ---------- residuals (safe)
  #   res <- suppressWarnings(tryCatch(as.numeric(residuals(fit)), error = function(e) numeric(0)))
  #   res <- res[is.finite(res)]
  #   n_res <- length(res)
  #   
  #   fitdf <- if (!is.null(est)) length(est) else 0L
  #   
  #   lb_lag <- min(L, max(1L, floor(n_res / 3)))
  #   lb <- if (n_res >= 5) {
  #     tryCatch(stats::Box.test(res, lag = lb_lag, type = "Ljung-Box", fitdf = fitdf), error = function(e) NULL)
  #   } else NULL
  #   
  #   bp <- if (n_res >= 5) {
  #     tryCatch(stats::Box.test(res, lag = lb_lag, type = "Box-Pierce", fitdf = fitdf), error = function(e) NULL)
  #   } else NULL
  #   
  #   jb <- if (requireNamespace("tseries", quietly = TRUE) && n_res >= 5) {
  #     tryCatch(tseries::jarque.bera.test(res), error = function(e) NULL)
  #   } else NULL
  #   
  #   sw <- if (n_res >= 3 && n_res <= 5000) {
  #     tryCatch(stats::shapiro.test(res), error = function(e) NULL)
  #   } else NULL
  #   
  #   arch <- if (requireNamespace("FinTS", quietly = TRUE) && n_res >= 10) {
  #     tryCatch(FinTS::ArchTest(res, lags = min(12L, max(1L, floor(L / 2)))), error = function(e) NULL)
  #   } else NULL
  #   
  #   runs <- if (requireNamespace("tseries", quietly = TRUE) && n_res >= 10) {
  #     tryCatch(tseries::runs.test(res), error = function(e) NULL)
  #   } else NULL
  #   
  #   ad <- if (requireNamespace("nortest", quietly = TRUE) && n_res >= 8) {
  #     tryCatch(nortest::ad.test(res), error = function(e) NULL)
  #   } else NULL
  #   
  #   test_rows <- list()
  #   add_test <- function(name, stat, pval, note = "") {
  #     test_rows[[length(test_rows) + 1L]] <<- data.frame(
  #       Test = name,
  #       Statistic = if (is.null(stat)) "NA" else fmt_num_local(stat, 3),
  #       `p-value` = if (is.null(pval)) "NA" else fmt_p_local(pval),
  #       Interpretation = note,
  #       check.names = FALSE
  #     )
  #   }
  #   
  #   add_test(
  #     "Ljung–Box (residuals)",
  #     if (!is.null(lb)) unname(lb$statistic) else NULL,
  #     if (!is.null(lb)) unname(lb$p.value) else NULL,
  #     if (!is.null(lb)) {
  #       if (is.finite(lb$p.value) && lb$p.value >= 0.05) "No evidence of remaining autocorrelation (white-noise compatible)."
  #       else "Evidence of remaining autocorrelation (model may be under-specified)."
  #     } else "Not available."
  #   )
  #   
  #   add_test(
  #     "Box–Pierce (residuals)",
  #     if (!is.null(bp)) unname(bp$statistic) else NULL,
  #     if (!is.null(bp)) unname(bp$p.value) else NULL,
  #     if (!is.null(bp)) {
  #       if (is.finite(bp$p.value) && bp$p.value >= 0.05) "Consistent with uncorrelated residuals."
  #       else "Suggests residual autocorrelation."
  #     } else "Not available."
  #   )
  #   
  #   add_test(
  #     "Jarque–Bera (normality)",
  #     if (!is.null(jb)) unname(jb$statistic) else NULL,
  #     if (!is.null(jb)) unname(jb$p.value) else NULL,
  #     if (!is.null(jb)) {
  #       if (is.finite(jb$p.value) && jb$p.value >= 0.05) "No evidence against normality."
  #       else "Residuals deviate from normality (common in real series)."
  #     } else "Package 'tseries' missing or test unavailable."
  #   )
  #   
  #   add_test(
  #     "Shapiro–Wilk (normality)",
  #     if (!is.null(sw)) unname(sw$statistic) else NULL,
  #     if (!is.null(sw)) unname(sw$p.value) else NULL,
  #     if (!is.null(sw)) {
  #       if (is.finite(sw$p.value) && sw$p.value >= 0.05) "No evidence against normality."
  #       else "Evidence against normality."
  #     } else if (n_res > 5000) "Not computed (n > 5000)." else "Not available."
  #   )
  #   
  #   add_test(
  #     "ARCH LM (heteroskedasticity)",
  #     if (!is.null(arch)) unname(arch$statistic) else NULL,
  #     if (!is.null(arch)) unname(arch$p.value) else NULL,
  #     if (!is.null(arch)) {
  #       if (is.finite(arch$p.value) && arch$p.value >= 0.05) "No evidence of remaining ARCH effects."
  #       else "Evidence of ARCH effects → consider GARCH for variance."
  #     } else "Package 'FinTS' missing or test unavailable."
  #   )
  #   
  #   add_test(
  #     "Runs test (randomness)",
  #     if (!is.null(runs)) unname(runs$statistic) else NULL,
  #     if (!is.null(runs)) unname(runs$p.value) else NULL,
  #     if (!is.null(runs)) {
  #       if (is.finite(runs$p.value) && runs$p.value >= 0.05) "No evidence against randomness."
  #       else "Evidence of non-randomness (structure may remain)."
  #     } else "Package 'tseries' missing or test unavailable."
  #   )
  #   
  #   add_test(
  #     "Anderson–Darling (normality)",
  #     if (!is.null(ad)) unname(ad$statistic) else NULL,
  #     if (!is.null(ad)) unname(ad$p.value) else NULL,
  #     if (!is.null(ad)) {
  #       if (is.finite(ad$p.value) && ad$p.value >= 0.05) "No evidence against normality."
  #       else "Evidence against normality (sensitive in tails)."
  #     } else "Package 'nortest' missing or test unavailable."
  #   )
  #   
  #   tests_df <- if (length(test_rows)) do.call(rbind, test_rows) else data.frame()
  #   
  #   # ---------- forecast accuracy (safe)
  #   acc_df <- NULL
  #   acc_sentence <- "No holdout test set was detected; therefore, out-of-sample accuracy was not computed."
  #   
  #   y_test <- s$ts_test
  #   has_test <- !is.null(y_test) && safe_len(y_test) > 0 && n_test > 0
  #   
  #   if (has_test) {
  #     y_test_num <- as.numeric(y_test)
  #     y_hat_num  <- suppressWarnings(tryCatch(as.numeric(fc$mean), error = function(e) rep(NA_real_, length(y_test_num))))
  #     h <- min(length(y_test_num), length(y_hat_num))
  #     if (h >= 1) {
  #       e <- y_test_num[seq_len(h)] - y_hat_num[seq_len(h)]
  #       rmse <- sqrt(mean(e^2, na.rm = TRUE))
  #       mae  <- mean(abs(e), na.rm = TRUE)
  #       mape <- mean(abs(e) / pmax(abs(y_test_num[seq_len(h)]), .Machine$double.eps), na.rm = TRUE)
  #       
  #       acc_df <- data.frame(
  #         Metric = c("RMSE", "MAE", "MAPE"),
  #         Value  = c(fmt_num_local(rmse, 3), fmt_num_local(mae, 3), paste0(fmt_num_local(100*mape, 2), "%")),
  #         check.names = FALSE
  #       )
  #       
  #       acc_sentence <- paste0(
  #         "Over the holdout period (test n = ", h, "), forecast performance was ",
  #         "RMSE = ", fmt_num_local(rmse, 3), ", ",
  #         "MAE = ", fmt_num_local(mae, 3), ", ",
  #         "MAPE = ", fmt_num_local(100*mape, 2), "%."
  #       )
  #     }
  #   }
  #   
  #   # ---------- horizon narrative
  #   horizon_txt <- if (has_test) {
  #     paste0("Validation mode was used: the forecast horizon was forced to match the test length (h = ", n_test, ").")
  #   } else {
  #     paste0("Future mode was used: forecasts were produced beyond the training sample (h = ", fc0$h, ").")
  #   }
  #   
  #   # ---------- model string
  #   season_txt <- suppressWarnings(as.integer(eq$s))
  #   s_txt <- if (is.finite(season_txt) && season_txt > 0) as.character(season_txt) else "s"
  #   
  #   model_str <- sprintf(
  #     "SARIMA(%d,%d,%d)(%d,%d,%d)[%s]",
  #     eq$p, eq$d, eq$q, eq$P, eq$D, eq$Q, s_txt
  #   )
  #   
  #   # ---------- diagnostic verdict
  #   lb_ok <- !is.null(lb) && is.finite(lb$p.value) && lb$p.value >= 0.05
  #   arch_ok <- is.null(arch) || (is.finite(arch$p.value) && arch$p.value >= 0.05)
  #   diag_verdict <- paste0(
  #     if (lb_ok) "Residual autocorrelation was not statistically detected (Ljung–Box p ≥ .05). "
  #     else "Residual autocorrelation may remain (Ljung–Box p < .05). ",
  #     if (arch_ok) "No clear evidence of residual ARCH effects was found (or test unavailable)."
  #     else "Residual ARCH effects were detected → a GARCH extension is recommended."
  #   )
  #   
  #   # ---------- build report UI (NO output$ definitions here)
  #   tagList(
  #     tags$h3("Manual SARIMA: Full academic conclusion (report-ready)"),
  #     
  #     
  #     # 1. Objective and modelling rationale
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("1. Objective and modelling rationale")),
  #     tags$p(
  #       "A manually specified seasonal ARIMA (SARIMA) model was estimated to represent linear temporal dependence,",
  #       " including seasonal structure, and to provide an interpretable baseline for forecasting."
  #     ),
  #     
  #     
  #     # 2. Data design and sample split
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("2. Data design and sample split")),
  #     tags$p(HTML(paste0(
  #       "The analysis used <b>N = ", N, "</b> observations (training <b>n = ", n_train, "</b>",
  #       if (has_test) paste0(", test <b>n = ", n_test, "</b>") else "",
  #       ")."
  #     ))),
  #     tags$p(tags$b("Forecast design. "), horizon_txt),
  #     
  #     
  #     # 3. Identification visuals (time series + ACF/PACF)
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("3. Identification visuals (time series + ACF/PACF)")),
  #     tags$p(
  #       "The plots below summarize the observed series (with the train/test split, if applicable), ",
  #       "followed by ACF/PACF for the training series and for the differenced series implied by the chosen (d, D, s)."
  #     ),
  #     
  #     tags$br(), 
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure A. Time series with split marker")),
  #     plotOutput("manual_report_ts_plot", height = 360),
  #     
  #     # 4. Stationarity assessment (ADF, KPSS, and Phillips–Perron)
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("4. Stationarity assessment (ADF, KPSS, and Phillips–Perron)")),
  #     tags$p("Stationarity tests were applied to the training series to evaluate whether differencing is required before SARIMA identification and estimation."),
  #     tags$hr(),
  #     uiOutput("manual_report_stationarity"),
  #     
  #     
  #     # 5. Transformed series: differencing and seasonality
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("5. Transformed series: differencing and seasonality")),
  #     
  # 
  #     tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure B. Seasonal subseries")),
  #     plotOutput("manual_report_subseries", height = 360),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure C. Seasonal box-plot")),
  #     plotOutput("manual_report_seasonal_box", height = 360),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure D. Transformed training series and differenced (d, D, s) series")),
  #     plotOutput("manual_report_ts_trans_and_diff", height = 360),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure E. ACF and PACF (training series)")),
  #     fluidRow(
  #       column(6, plotOutput("manual_report_acf",  height = 280)),
  #       column(6, plotOutput("manual_report_pacf", height = 280))
  #     ),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure F. ACF and PACF (modified / differenced series using current d, D, s)")),
  #     fluidRow(
  #       column(6, plotOutput("manual_report_acf_mod",  height = 280)),
  #       column(6, plotOutput("manual_report_pacf_mod", height = 280))
  #     ),
  #     
  #     
  #     tags$hr(), tags$br(),
  #     uiOutput("manual_report_stationarity_mod"),
  #     
  #     # 6. Final model specification and fit quality
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("6. Final model specification and fit quality")),
  #     tags$p(HTML(paste0(
  #       "The final manual specification was <b>", model_str, "</b>",
  #       if (isTRUE(input$manual_drift)) " including drift/mean." else " without drift/mean.",
  #       " Model adequacy was assessed using information criteria, coefficient inference, residual diagnostics, and forecast performance."
  #     ))),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Table A. Goodness-of-fit (information criteria)")),
  #     html_table(ic_df),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Table B. Parameter estimates and significance (approx. z/t tests)")),
  #     if (!is.null(coef_df)) html_table(coef_df) else tags$em("No coefficients available."),
  #     tags$p(HTML(paste0(
  #       "In total, <b>", n_sig, "</b> parameter(s) were significant at α = .05 (marked by *, **, ***)."
  #     ))),
  #     
  #     
  #     
  #     # 7. Model equations (replication-ready
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("7. Model equations (replication-ready)")),
  #     tags$p(
  #       "The fitted model is reported below in operator notation (general form), expanded form, and the numerical equation",
  #       " based on the estimated parameters."
  #     ),
  #     tags$div(
  #       style = "padding:10px;border:1px solid #e5e5e5;border-radius:6px;background:#fcfcfc;",
  #       
  #       tags$br(),
  #       tags$hr(),
  #       tags$hr(),
  #       
  #       tags$h5("General SARIMA formulation"),
  #       HTML(eq$eq_general),
  #       
  #       tags$br(),
  #       tags$hr(),
  #       tags$hr(),
  #       
  #       tags$h5("Expanded operator form"),
  #       HTML(eq$eq_expanded),
  #       
  #       tags$br(),
  #       tags$hr(),
  #       tags$hr(),
  #       
  #       tags$h5("Numerical model"),
  #       HTML(eq$eq_line3),
  #       
  #       tags$br(),
  #       tags$hr(),
  #       tags$hr(),
  #       
  #       HTML(eq$eq_line4),
  #       
  #       tags$br(),
  #       tags$hr(),
  #       tags$hr(),
  #       
  #     ),
  #     
  #     
  #     
  #     # 8. Residual diagnostics (graphical evidence
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("8. Residual diagnostics (graphical evidence)")),
  #     tags$p(
  #       "Graphical diagnostics evaluate whether residuals resemble white noise (no systematic autocorrelation),",
  #       " approximate normality (Q–Q and histogram), and stable variance."
  #     ),
  #     
  #     fluidRow(
  #       column(6, plotOutput("manual_resid_ts",   height = 220)),
  #       column(6, plotOutput("manual_resid_acf",  height = 220))
  #     ),
  #     fluidRow(
  #       column(6, plotOutput("manual_resid_hist", height = 220)),
  #       column(6, plotOutput("manual_resid_qq",   height = 220))
  #     ),
  #     
  #     
  #     tags$h5("Ljung–Box p-values by lag"),
  #     plotOutput("manual_resid_lb_pvals", height = 260),
  #     
  #     
  #     
  #     # "9. Residual tests (formal inference
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("9. Residual tests (formal inference)")),
  #     tags$p(
  #       "Formal tests complement the plots: Ljung–Box/Box–Pierce assess remaining autocorrelation;",
  #       " Jarque–Bera/Shapiro–Wilk/Anderson–Darling assess normality;",
  #       " ARCH LM tests conditional heteroskedasticity; the runs test checks randomness."
  #     ),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Table C. Residual test summary")),
  #     html_table(tests_df),
  #     tags$p(tags$b("Diagnostic synthesis. "), diag_verdict),
  #     
  #     
  #     
  #     # 10. Forecasting results and predictive performance
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("10. Forecasting results and predictive performance")),
  #     tags$p(acc_sentence),
  #     
  #     if (!is.null(acc_df)) tagList(
  #       tags$h5("Table D. Holdout accuracy (test set)"),
  #       html_table(acc_df)
  #     ) else NULL,
  #     
  #     tags$hr(),  tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Forecast plot")),
  #     plotOutput("manual_forecast_plot", height = 420),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Forecast table")),
  #     tableOutput("manual_forecast_table"),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Accuracy table (your app output)")),
  #     tableOutput("manual_accuracy_table"),
  #     
  #     
  #     # 11. Final conclusion (academic)
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("11. Final conclusion")),
  #     tags$p(
  #       "Overall, the manually specified SARIMA model provides a coherent and interpretable representation of seasonal linear dynamics,",
  #       " supported by information criteria, statistically interpretable parameters, and diagnostic checks.",
  #       " When diagnostics indicate remaining autocorrelation, refinement should prioritize revising differencing (d, D) and AR/MA orders guided by ACF/PACF and Ljung–Box.",
  #       " When conditional heteroskedasticity is detected (ARCH LM), a volatility model (e.g., GARCH) should be added to the mean equation to better represent time-varying variance."
  #     ),
  #     tags$p(
  #       "For reporting, the results above provide the full chain of evidence typically expected in academic manuscripts:",
  #       " (i) specification + IC, (ii) parameter inference, (iii) equation reporting, (iv) residual validation with plots and tests, and (v) forecasting with accuracy assessment."
  #     ),
  #     
  #     tags$br(), tags$hr(), tags$br(), tags$br(),
  #     
  #   )
  # })
  
  
  
  
  
  # # ------------------------------------------------------------
  # # (B) manual_conclusion_full_obj (eventReactive) — UI builder ONLY
  # # ------------------------------------------------------------
  # manual_conclusion_full_obj <- eventReactive(input$fit_manual, {
  #   req(manual_fit(), manual_fc(), manual_equations(), ts_train_test())
  #   
  #   fit <- manual_fit()
  #   fc0 <- manual_fc()
  #   fc  <- fc0$fc
  #   eq  <- manual_equations()
  #   s   <- ts_train_test()
  #   
  #   # ---------- helpers (local + safe)
  #   `%||%` <- function(x, y) {
  #     if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
  #   }
  #   
  #   fmt_num_local <- function(x, d = 3) {
  #     if (length(x) == 0 || all(is.na(x))) return("NA")
  #     x <- suppressWarnings(as.numeric(x[1]))
  #     if (!is.finite(x)) return("NA")
  #     formatC(x, format = "f", digits = d)
  #   }
  #   fmt_p_local <- function(p) {
  #     if (length(p) == 0 || all(is.na(p))) return("NA")
  #     p <- suppressWarnings(as.numeric(p[1]))
  #     if (!is.finite(p)) return("NA")
  #     if (p < .001) "&lt; .001" else sprintf("= %.3f", p)
  #   }
  #   sig_stars <- function(p) {
  #     p <- suppressWarnings(as.numeric(p))
  #     if (!is.finite(p)) return("")
  #     if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else if (p < .10) "†" else ""
  #   }
  #   safe_len <- function(x) if (is.null(x)) 0L else length(x)
  #   
  #   html_table <- function(df) {
  #     if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
  #       return(tags$em("Table unavailable."))
  #     }
  #     tags$table(
  #       class = "table table-striped table-condensed",
  #       tags$thead(tags$tr(lapply(names(df), function(nm) tags$th(nm)))),
  #       tags$tbody(
  #         lapply(seq_len(nrow(df)), function(i) {
  #           tags$tr(lapply(df[i, , drop = FALSE], function(cell) tags$td(HTML(as.character(cell)))))
  #         })
  #       )
  #     )
  #   }
  #   
  #   # ---------- sample sizes (safe)
  #   n_train <- suppressWarnings(as.integer(s$train_n))
  #   if (!is.finite(n_train) || n_train < 1) n_train <- tryCatch(length(residuals(fit)), error = function(e) 0L)
  #   
  #   n_test <- suppressWarnings(as.integer(s$test_n))
  #   if (!is.finite(n_test) || n_test < 0) n_test <- 0L
  #   
  #   N <- n_train + n_test
  #   
  #   # ---------- lag choice (safe)
  #   L_in <- suppressWarnings(as.integer(input$diag_lag))
  #   L <- if (is.finite(L_in) && L_in > 0) L_in else 12L
  #   
  #   # ============================================================
  #   # ---------- IC (safe + robust AICc fallback)
  #   # ============================================================
  #   AIC_val <- suppressWarnings(tryCatch(as.numeric(stats::AIC(fit)), error = function(e) NA_real_))
  #   BIC_val <- suppressWarnings(tryCatch(as.numeric(stats::BIC(fit)), error = function(e) NA_real_))
  #   
  #   AICc_val <- suppressWarnings(tryCatch(as.numeric(forecast::AICc(fit)), error = function(e) NA_real_))
  #   
  #   if (!is.finite(AICc_val) && is.finite(AIC_val)) {
  #     n_fit <- suppressWarnings(tryCatch(stats::nobs(fit), error = function(e) NA_integer_))
  #     if (!is.finite(n_fit) || n_fit <= 0) {
  #       r <- suppressWarnings(tryCatch(as.numeric(residuals(fit)), error = function(e) numeric(0)))
  #       n_fit <- sum(is.finite(r))
  #     }
  #     if (!is.finite(n_fit) || n_fit <= 0) n_fit <- n_train
  #     
  #     k <- suppressWarnings(tryCatch(length(stats::coef(fit)), error = function(e) NA_integer_))
  #     if (!is.finite(k) || k < 0) k <- 0L
  #     k <- k + 1L
  #     
  #     if (is.finite(n_fit) && n_fit > (k + 1)) {
  #       AICc_val <- AIC_val + (2 * k * (k + 1)) / (n_fit - k - 1)
  #     } else {
  #       AICc_val <- NA_real_
  #     }
  #   }
  #   
  #   ic_df <- data.frame(
  #     Criterion = c("AIC", "AICc", "BIC"),
  #     Value     = c(fmt_num_local(AIC_val, 2), fmt_num_local(AICc_val, 2), fmt_num_local(BIC_val, 2)),
  #     check.names = FALSE
  #   )
  #   
  #   # ---------- coefficients + significance (robust)
  #   est <- suppressWarnings(tryCatch(stats::coef(fit), error = function(e) NULL))
  #   V   <- suppressWarnings(tryCatch(stats::vcov(fit), error = function(e) NULL))
  #   if (is.null(V)) V <- suppressWarnings(tryCatch(fit$var.coef, error = function(e) NULL))
  #   
  #   coef_df <- NULL
  #   if (!is.null(est) && length(est) > 0) {
  #     est <- as.numeric(est)
  #     nm  <- names(stats::coef(fit))
  #     if (is.null(nm)) nm <- paste0("param_", seq_along(est))
  #     
  #     se <- rep(NA_real_, length(est))
  #     if (!is.null(V)) {
  #       dV <- tryCatch(diag(V), error = function(e) rep(NA_real_, length(est)))
  #       if (length(dV) == length(est)) se <- sqrt(pmax(dV, 0))
  #     }
  #     z  <- est / se
  #     p  <- 2 * stats::pnorm(-abs(z))
  #     
  #     coef_df <- data.frame(
  #       Term     = nm,
  #       Estimate = sprintf("%.6f", est),
  #       SE       = ifelse(is.finite(se), sprintf("%.6f", se), "NA"),
  #       `z/t`    = ifelse(is.finite(z),  sprintf("%.3f",  z),  "NA"),
  #       `p`      = ifelse(is.finite(p),  sprintf("%.3f",  p),  "NA"),
  #       Sig      = vapply(p, sig_stars, character(1)),
  #       check.names = FALSE
  #     )
  #   }
  #   
  #   n_sig <- if (!is.null(coef_df)) sum(suppressWarnings(as.numeric(coef_df$p)) < 0.05, na.rm = TRUE) else 0L
  #   
  #   # ---------- residuals (safe)
  #   res <- suppressWarnings(tryCatch(as.numeric(residuals(fit)), error = function(e) numeric(0)))
  #   res <- res[is.finite(res)]
  #   n_res <- length(res)
  #   
  #   fitdf <- if (!is.null(est)) length(est) else 0L
  #   
  #   lb_lag <- min(L, max(1L, floor(n_res / 3)))
  #   lb <- if (n_res >= 5) {
  #     tryCatch(stats::Box.test(res, lag = lb_lag, type = "Ljung-Box", fitdf = fitdf), error = function(e) NULL)
  #   } else NULL
  #   
  #   bp <- if (n_res >= 5) {
  #     tryCatch(stats::Box.test(res, lag = lb_lag, type = "Box-Pierce", fitdf = fitdf), error = function(e) NULL)
  #   } else NULL
  #   
  #   jb <- if (requireNamespace("tseries", quietly = TRUE) && n_res >= 5) {
  #     tryCatch(tseries::jarque.bera.test(res), error = function(e) NULL)
  #   } else NULL
  #   
  #   sw <- if (n_res >= 3 && n_res <= 5000) {
  #     tryCatch(stats::shapiro.test(res), error = function(e) NULL)
  #   } else NULL
  #   
  #   arch <- if (requireNamespace("FinTS", quietly = TRUE) && n_res >= 10) {
  #     tryCatch(FinTS::ArchTest(res, lags = min(12L, max(1L, floor(L / 2)))), error = function(e) NULL)
  #   } else NULL
  #   
  #   runs <- if (requireNamespace("tseries", quietly = TRUE) && n_res >= 10) {
  #     tryCatch(tseries::runs.test(res), error = function(e) NULL)
  #   } else NULL
  #   
  #   ad <- if (requireNamespace("nortest", quietly = TRUE) && n_res >= 8) {
  #     tryCatch(nortest::ad.test(res), error = function(e) NULL)
  #   } else NULL
  #   
  #   test_rows <- list()
  #   add_test <- function(name, stat, pval, note = "") {
  #     test_rows[[length(test_rows) + 1L]] <<- data.frame(
  #       Test = name,
  #       Statistic = if (is.null(stat)) "NA" else fmt_num_local(stat, 3),
  #       `p-value` = if (is.null(pval)) "NA" else fmt_p_local(pval),
  #       Interpretation = note,
  #       check.names = FALSE
  #     )
  #   }
  #   
  #   add_test(
  #     "Ljung–Box (residuals)",
  #     if (!is.null(lb)) unname(lb$statistic) else NULL,
  #     if (!is.null(lb)) unname(lb$p.value) else NULL,
  #     if (!is.null(lb)) {
  #       if (is.finite(lb$p.value) && lb$p.value >= 0.05) "No evidence of remaining autocorrelation (white-noise compatible)."
  #       else "Evidence of remaining autocorrelation (model may be under-specified)."
  #     } else "Not available."
  #   )
  #   
  #   add_test(
  #     "Box–Pierce (residuals)",
  #     if (!is.null(bp)) unname(bp$statistic) else NULL,
  #     if (!is.null(bp)) unname(bp$p.value) else NULL,
  #     if (!is.null(bp)) {
  #       if (is.finite(bp$p.value) && bp$p.value >= 0.05) "Consistent with uncorrelated residuals."
  #       else "Suggests residual autocorrelation."
  #     } else "Not available."
  #   )
  #   
  #   add_test(
  #     "Jarque–Bera (normality)",
  #     if (!is.null(jb)) unname(jb$statistic) else NULL,
  #     if (!is.null(jb)) unname(jb$p.value) else NULL,
  #     if (!is.null(jb)) {
  #       if (is.finite(jb$p.value) && jb$p.value >= 0.05) "No evidence against normality."
  #       else "Residuals deviate from normality (common in real series)."
  #     } else "Package 'tseries' missing or test unavailable."
  #   )
  #   
  #   add_test(
  #     "Shapiro–Wilk (normality)",
  #     if (!is.null(sw)) unname(sw$statistic) else NULL,
  #     if (!is.null(sw)) unname(sw$p.value) else NULL,
  #     if (!is.null(sw)) {
  #       if (is.finite(sw$p.value) && sw$p.value >= 0.05) "No evidence against normality."
  #       else "Evidence against normality."
  #     } else if (n_res > 5000) "Not computed (n > 5000)." else "Not available."
  #   )
  #   
  #   add_test(
  #     "ARCH LM (heteroskedasticity)",
  #     if (!is.null(arch)) unname(arch$statistic) else NULL,
  #     if (!is.null(arch)) unname(arch$p.value) else NULL,
  #     if (!is.null(arch)) {
  #       if (is.finite(arch$p.value) && arch$p.value >= 0.05) "No evidence of remaining ARCH effects."
  #       else "Evidence of ARCH effects → consider GARCH for variance."
  #     } else "Package 'FinTS' missing or test unavailable."
  #   )
  #   
  #   add_test(
  #     "Runs test (randomness)",
  #     if (!is.null(runs)) unname(runs$statistic) else NULL,
  #     if (!is.null(runs)) unname(runs$p.value) else NULL,
  #     if (!is.null(runs)) {
  #       if (is.finite(runs$p.value) && runs$p.value >= 0.05) "No evidence against randomness."
  #       else "Evidence of non-randomness (structure may remain)."
  #     } else "Package 'tseries' missing or test unavailable."
  #   )
  #   
  #   add_test(
  #     "Anderson–Darling (normality)",
  #     if (!is.null(ad)) unname(ad$statistic) else NULL,
  #     if (!is.null(ad)) unname(ad$p.value) else NULL,
  #     if (!is.null(ad)) {
  #       if (is.finite(ad$p.value) && ad$p.value >= 0.05) "No evidence against normality."
  #       else "Evidence against normality (sensitive in tails)."
  #     } else "Package 'nortest' missing or test unavailable."
  #   )
  #   
  #   tests_df <- if (length(test_rows)) do.call(rbind, test_rows) else data.frame()
  #   
  #   # ---------- forecast accuracy (safe)  [NOTE: this is on the model scale]
  #   acc_df <- NULL
  #   acc_sentence <- "No holdout test set was detected; therefore, out-of-sample accuracy was not computed."
  #   
  #   y_test <- s$ts_test
  #   has_test <- !is.null(y_test) && safe_len(y_test) > 0 && n_test > 0
  #   
  #   if (has_test) {
  #     y_test_num <- as.numeric(y_test)
  #     y_hat_num  <- suppressWarnings(tryCatch(as.numeric(fc$mean), error = function(e) rep(NA_real_, length(y_test_num))))
  #     h <- min(length(y_test_num), length(y_hat_num))
  #     if (h >= 1) {
  #       e <- y_test_num[seq_len(h)] - y_hat_num[seq_len(h)]
  #       rmse <- sqrt(mean(e^2, na.rm = TRUE))
  #       mae  <- mean(abs(e), na.rm = TRUE)
  #       mape <- mean(abs(e) / pmax(abs(y_test_num[seq_len(h)]), .Machine$double.eps), na.rm = TRUE)
  #       
  #       acc_df <- data.frame(
  #         Metric = c("RMSE", "MAE", "MAPE"),
  #         Value  = c(fmt_num_local(rmse, 3), fmt_num_local(mae, 3), paste0(fmt_num_local(100*mape, 2), "%")),
  #         check.names = FALSE
  #       )
  #       
  #       acc_sentence <- paste0(
  #         "Over the holdout period (test n = ", h, "), forecast performance was ",
  #         "RMSE = ", fmt_num_local(rmse, 3), ", ",
  #         "MAE = ", fmt_num_local(mae, 3), ", ",
  #         "MAPE = ", fmt_num_local(100*mape, 2), "%."
  #       )
  #     }
  #   }
  #   
  #   # ---------- horizon narrative
  #   horizon_txt <- if (has_test) {
  #     paste0("Validation mode was used: the forecast horizon was forced to match the test length (h = ", n_test, ").")
  #   } else {
  #     paste0("Future mode was used: forecasts were produced beyond the training sample (h = ", fc0$h, ").")
  #   }
  #   
  #   # ---------- model string
  #   season_txt <- suppressWarnings(as.integer(eq$s))
  #   s_txt <- if (is.finite(season_txt) && season_txt > 0) as.character(season_txt) else "s"
  #   
  #   model_str <- sprintf(
  #     "SARIMA(%d,%d,%d)(%d,%d,%d)[%s]",
  #     eq$p, eq$d, eq$q, eq$P, eq$D, eq$Q, s_txt
  #   )
  #   
  #   # ---------- diagnostic verdict
  #   lb_ok <- !is.null(lb) && is.finite(lb$p.value) && lb$p.value >= 0.05
  #   arch_ok <- is.null(arch) || (is.finite(arch$p.value) && arch$p.value >= 0.05)
  #   diag_verdict <- paste0(
  #     if (lb_ok) "Residual autocorrelation was not statistically detected (Ljung–Box p ≥ .05). "
  #     else "Residual autocorrelation may remain (Ljung–Box p < .05). ",
  #     if (arch_ok) "No clear evidence of residual ARCH effects was found (or test unavailable)."
  #     else "Residual ARCH effects were detected → a GARCH extension is recommended."
  #   )
  #   
  #   # ============================================================
  #   # ---------- NEW: Inverse transform equation + bias-adjusted back-forecasts
  #   # ============================================================
  #   tr_global <- input$transform %||% "none"
  #   
  #   # We need the same lambda as the transform step uses.
  #   # apply_transform() uses forecast::BoxCox.lambda(y, lower=0) when input$lambda is NA.
  #   p_obj <- NULL
  #   lambda_used <- NA_real_
  #   if (identical(tr_global, "boxcox")) {
  #     p_obj <- tryCatch(prepared(), error = function(e) NULL)
  #     lam_in <- input$lambda
  #     if (is.null(lam_in) || (length(lam_in) == 1 && is.na(lam_in))) {
  #       lambda_used <- tryCatch(
  #         forecast::BoxCox.lambda(p_obj$df$y_filled, lower = 0),
  #         error = function(e) NA_real_
  #       )
  #     } else {
  #       lambda_used <- suppressWarnings(as.numeric(lam_in))
  #     }
  #   }
  #   
  #   # Forecasts are on TRANSFORMED scale:
  #   mu_t <- suppressWarnings(as.numeric(fc$mean))
  #   se_t <- suppressWarnings(as.numeric(fc$se))
  #   sigma2_t <- ifelse(is.finite(se_t), se_t^2, NA_real_)
  #   
  #   lo_t <- tryCatch(fc$lower, error = function(e) NULL)
  #   hi_t <- tryCatch(fc$upper, error = function(e) NULL)
  #   lvl_names <- tryCatch(colnames(fc$lower), error = function(e) NULL)
  #   
  #   inv_equation_html <- NULL
  #   inv_note_html <- NULL
  #   
  #   inv_mean <- mu_t
  #   inv_lo   <- lo_t
  #   inv_hi   <- hi_t
  #   
  #   # helper: nice TeX if tex_display exists, else plain
  #   tex_or_plain <- function(tex, plain) {
  #     if (exists("tex_display", mode = "function")) HTML(tex_display(tex)) else HTML(plain)
  #   }
  #   
  #   if (identical(tr_global, "log")) {
  #     # unbiased mean under lognormal assumption: E[Y] = exp(mu + 0.5*sigma^2)
  #     inv_mean <- exp(mu_t + 0.5 * sigma2_t)
  #     
  #     # interval bounds: back-transform quantiles (standard practice)
  #     if (!is.null(lo_t) && !is.null(hi_t)) {
  #       inv_lo <- exp(lo_t)
  #       inv_hi <- exp(hi_t)
  #     }
  #     
  #     inv_equation_html <- tex_or_plain(
  #       "y = \\exp(z) \\quad \\text{where } z = \\ln(y)",
  #       "y = exp(z), where z = log(y)"
  #     )
  #     inv_note_html <- HTML(
  #       "<span style='font-size:12px;color:#444;'>Point forecasts are bias-adjusted (lognormal): <b>exp(μ + 0.5·σ²)</b>. Interval bounds are back-transformed using <b>exp()</b>.</span>"
  #     )
  #     
  #   } else if (identical(tr_global, "boxcox")) {
  #     lam <- lambda_used
  #     
  #     inv_mean <- tryCatch(
  #       forecast::InvBoxCox(mu_t, lam, biasadj = TRUE, sigma2 = sigma2_t),
  #       error = function(e) rep(NA_real_, length(mu_t))
  #     )
  #     
  #     if (!is.null(lo_t) && !is.null(hi_t)) {
  #       inv_lo <- tryCatch(forecast::InvBoxCox(lo_t, lam), error = function(e) lo_t)
  #       inv_hi <- tryCatch(forecast::InvBoxCox(hi_t, lam), error = function(e) hi_t)
  #     }
  #     
  #     if (is.finite(lam) && abs(lam) < 1e-6) {
  #       inv_equation_html <- tex_or_plain("y = \\exp(z) \\quad (\\lambda \\approx 0)", "y = exp(z)  (lambda ≈ 0)")
  #     } else {
  #       inv_equation_html <- tex_or_plain("y = (\\lambda z + 1)^{1/\\lambda}", "y = (lambda*z + 1)^(1/lambda)")
  #     }
  #     
  #     inv_note_html <- HTML(paste0(
  #       "<span style='font-size:12px;color:#444;'>λ used = <b>",
  #       if (is.finite(lambda_used)) formatC(lambda_used, digits = 4, format = "f") else "NA",
  #       "</b>. Point forecasts use <b>InvBoxCox(..., biasadj=TRUE, sigma2=se²)</b>.</span>"
  #     ))
  #   }
  #   
  #   # Build an ORIGINAL-SCALE forecast table (bias-adjusted mean + back-transformed intervals)
  #   fc_orig_df <- data.frame(
  #     Horizon = seq_along(inv_mean),
  #     Mean_original = as.numeric(inv_mean),
  #     stringsAsFactors = FALSE
  #   )
  #   
  #   if (!is.null(inv_lo) && !is.null(inv_hi) && !is.null(lvl_names)) {
  #     for (j in seq_along(lvl_names)) {
  #       lvl <- lvl_names[j]
  #       fc_orig_df[[paste0("Lo_", lvl)]] <- as.numeric(inv_lo[, j])
  #       fc_orig_df[[paste0("Hi_", lvl)]] <- as.numeric(inv_hi[, j])
  #     }
  #   }
  #   
  #   # format numeric columns for HTML display
  #   fc_orig_df_fmt <- fc_orig_df
  #   for (nm in names(fc_orig_df_fmt)) {
  #     if (is.numeric(fc_orig_df_fmt[[nm]])) {
  #       fc_orig_df_fmt[[nm]] <- ifelse(is.finite(fc_orig_df_fmt[[nm]]),
  #                                      sprintf("%.6f", fc_orig_df_fmt[[nm]]),
  #                                      "NA")
  #     }
  #   }
  #   
  #   # ---------- build report UI (NO output$ definitions here)
  #   tagList(
  #     tags$h3("Manual SARIMA: Full academic conclusion (report-ready)"),
  #     
  #     # 1. Objective and modelling rationale
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("1. Objective and modelling rationale")),
  #     tags$p(
  #       "A manually specified seasonal ARIMA (SARIMA) model was estimated to represent linear temporal dependence,",
  #       " including seasonal structure, and to provide an interpretable baseline for forecasting."
  #     ),
  #     
  #     # 2. Data design and sample split
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("2. Data design and sample split")),
  #     tags$p(HTML(paste0(
  #       "The analysis used <b>N = ", N, "</b> observations (training <b>n = ", n_train, "</b>",
  #       if (has_test) paste0(", test <b>n = ", n_test, "</b>") else "",
  #       ")."
  #     ))),
  #     tags$p(tags$b("Forecast design. "), horizon_txt),
  #     
  #     # 3. Identification visuals (time series + ACF/PACF)
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("3. Identification visuals (time series + ACF/PACF)")),
  #     tags$p(
  #       "The plots below summarize the observed series (with the train/test split, if applicable), ",
  #       "followed by ACF/PACF for the training series and for the differenced series implied by the chosen (d, D, s)."
  #     ),
  #     
  #     tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure A. Time series with split marker")),
  #     plotOutput("manual_report_ts_plot", height = 360),
  #     
  #     # 4. Stationarity assessment (ADF, KPSS, and Phillips–Perron)
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("4. Stationarity assessment (ADF, KPSS, and Phillips–Perron)")),
  #     tags$p("Stationarity tests were applied to the training series to evaluate whether differencing is required before SARIMA identification and estimation."),
  #     tags$hr(),
  #     uiOutput("manual_report_stationarity"),
  #     
  #     # 5. Transformed series: differencing and seasonality
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("5. Transformed series: differencing and seasonality")),
  #     
  #     tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure B. Seasonal subseries")),
  #     plotOutput("manual_report_subseries", height = 360),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure C. Seasonal box-plot")),
  #     plotOutput("manual_report_seasonal_box", height = 360),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure D. Transformed training series and differenced (d, D, s) series")),
  #     plotOutput("manual_report_ts_trans_and_diff", height = 360),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure E. ACF and PACF (training series)")),
  #     fluidRow(
  #       column(6, plotOutput("manual_report_acf",  height = 280)),
  #       column(6, plotOutput("manual_report_pacf", height = 280))
  #     ),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure F. ACF and PACF (modified / differenced series using current d, D, s)")),
  #     fluidRow(
  #       column(6, plotOutput("manual_report_acf_mod",  height = 280)),
  #       column(6, plotOutput("manual_report_pacf_mod", height = 280))
  #     ),
  #     
  #     tags$hr(), tags$br(),
  #     uiOutput("manual_report_stationarity_mod"),
  #     
  #     # 6. Final model specification and fit quality
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("6. Final model specification and fit quality")),
  #     tags$p(HTML(paste0(
  #       "The final manual specification was <b>", model_str, "</b>",
  #       if (isTRUE(input$manual_drift)) " including drift/mean." else " without drift/mean.",
  #       " Model adequacy was assessed using information criteria, coefficient inference, residual diagnostics, and forecast performance."
  #     ))),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Table A. Goodness-of-fit (information criteria)")),
  #     html_table(ic_df),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Table B. Parameter estimates and significance (approx. z/t tests)")),
  #     if (!is.null(coef_df)) html_table(coef_df) else tags$em("No coefficients available."),
  #     tags$p(HTML(paste0(
  #       "In total, <b>", n_sig, "</b> parameter(s) were significant at α = .05 (marked by *, **, ***)."
  #     ))),
  #     
  #     # 7. Model equations (replication-ready)
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("7. Model equations (replication-ready)")),
  #     tags$p(
  #       "The fitted model is reported below in operator notation (general form), expanded form, and the numerical equation",
  #       " based on the estimated parameters."
  #     ),
  #     tags$div(
  #       style = "padding:10px;border:1px solid #e5e5e5;border-radius:6px;background:#fcfcfc;",
  #       
  #       tags$br(),
  #       tags$hr(),
  #       tags$hr(),
  #       
  #       tags$h5("General SARIMA formulation"),
  #       HTML(eq$eq_general),
  #       
  #       tags$br(),
  #       tags$hr(),
  #       tags$hr(),
  #       
  #       tags$h5("Expanded operator form"),
  #       HTML(eq$eq_expanded),
  #       
  #       tags$br(),
  #       tags$hr(),
  #       tags$hr(),
  #       
  #       tags$h5("Numerical model"),
  #       HTML(eq$eq_line3),
  #       
  #       tags$br(),
  #       tags$hr(),
  #       tags$hr(),
  #       
  #       HTML(eq$eq_line4),
  #       
  #       tags$br(),
  #       tags$hr(),
  #       tags$hr()
  #     ),
  #     
  #     # 8. Residual diagnostics (graphical evidence)
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("8. Residual diagnostics (graphical evidence)")),
  #     tags$p(
  #       "Graphical diagnostics evaluate whether residuals resemble white noise (no systematic autocorrelation),",
  #       " approximate normality (Q–Q and histogram), and stable variance."
  #     ),
  #     
  #     fluidRow(
  #       column(6, plotOutput("manual_resid_ts",   height = 220)),
  #       column(6, plotOutput("manual_resid_acf",  height = 220))
  #     ),
  #     fluidRow(
  #       column(6, plotOutput("manual_resid_hist", height = 220)),
  #       column(6, plotOutput("manual_resid_qq",   height = 220))
  #     ),
  #     
  #     tags$h5("Ljung–Box p-values by lag"),
  #     plotOutput("manual_resid_lb_pvals", height = 260),
  #     
  #     # 9. Residual tests (formal inference)
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("9. Residual tests (formal inference)")),
  #     tags$p(
  #       "Formal tests complement the plots: Ljung–Box/Box–Pierce assess remaining autocorrelation;",
  #       " Jarque–Bera/Shapiro–Wilk/Anderson–Darling assess normality;",
  #       " ARCH LM tests conditional heteroskedasticity; the runs test checks randomness."
  #     ),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Table C. Residual test summary")),
  #     html_table(tests_df),
  #     tags$p(tags$b("Diagnostic synthesis. "), diag_verdict),
  #     
  #     # 10. Forecasting results and predictive performance
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("10. Forecasting results and predictive performance")),
  #     tags$p(acc_sentence),
  #     
  #     if (!is.null(acc_df)) tagList(
  #       tags$h5("Table D. Holdout accuracy (test set)"),
  #       html_table(acc_df)
  #     ) else NULL,
  #     
  #     tags$hr(),  tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Forecast plot")),
  #     plotOutput("manual_forecast_plot", height = 420),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Forecast table")),
  #     tableOutput("manual_forecast_table"),
  #     
  #     tags$hr(), tags$br(),
  #     tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Accuracy table (your app output)")),
  #     tableOutput("manual_accuracy_table"),
  #     
  #     # ---------- NEW: Back-transformed forecasts (original scale)
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("10.B Back-transformed forecasts (original scale)")),
  #     if (!identical(tr_global, "none")) tagList(
  #       tags$p("Because the model was estimated on the transformed series, forecasts are reported below on the original measurement scale using an appropriate bias adjustment for the point forecasts."),
  #       tags$div(
  #         style = "padding:10px;border:1px solid #e5e5e5;border-radius:6px;background:#fcfcfc;",
  #         tags$h5("Inverse transform equation (simple form)"),
  #         inv_equation_html,
  #         tags$br(),
  #         inv_note_html
  #       ),
  #       tags$br(),
  #       tags$h5("Back-transformed forecast table (original scale)"),
  #       html_table(fc_orig_df_fmt)
  #     ) else tags$em("No transformation was applied; forecasts are already on the original scale."),
  #     
  #     # 11. Final conclusion (academic)
  #     tags$hr(), tags$br(),
  #     tags$h4(tags$strong("11. Final conclusion")),
  #     tags$p(
  #       "Overall, the manually specified SARIMA model provides a coherent and interpretable representation of seasonal linear dynamics,",
  #       " supported by information criteria, statistically interpretable parameters, and diagnostic checks.",
  #       " When diagnostics indicate remaining autocorrelation, refinement should prioritize revising differencing (d, D) and AR/MA orders guided by ACF/PACF and Ljung–Box.",
  #       " When conditional heteroskedasticity is detected (ARCH LM), a volatility model (e.g., GARCH) should be added to the mean equation to better represent time-varying variance."
  #     ),
  #     tags$p(
  #       "For reporting, the results above provide the full chain of evidence typically expected in academic manuscripts:",
  #       " (i) specification + IC, (ii) parameter inference, (iii) equation reporting, (iv) residual validation with plots and tests, and (v) forecasting with accuracy assessment."
  #     ),
  #     
  #     tags$br(), tags$hr(), tags$br(), tags$br()
  #   )
  # })
  
  
  
  
  # ------------------------------------------------------------
  # (B) manual_conclusion_full_obj (eventReactive) — UI builder ONLY
  # ------------------------------------------------------------
  manual_conclusion_full_obj <- eventReactive(input$fit_manual, {
    req(manual_fit(), manual_fc(), manual_equations(), ts_train_test())
    
    fit <- manual_fit()
    fc0 <- manual_fc()
    fc  <- fc0$fc
    eq  <- manual_equations()
    s   <- ts_train_test()
    
    # ---------- helpers (local + safe)
    `%||%` <- function(x, y) {
      if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
    }
    
    fmt_num_local <- function(x, d = 3) {
      if (length(x) == 0 || all(is.na(x))) return("NA")
      x <- suppressWarnings(as.numeric(x[1]))
      if (!is.finite(x)) return("NA")
      formatC(x, format = "f", digits = d)
    }
    fmt_p_local <- function(p) {
      if (length(p) == 0 || all(is.na(p))) return("NA")
      p <- suppressWarnings(as.numeric(p[1]))
      if (!is.finite(p)) return("NA")
      if (p < .001) "&lt; .001" else sprintf("= %.3f", p)
    }
    sig_stars <- function(p) {
      p <- suppressWarnings(as.numeric(p))
      if (!is.finite(p)) return("")
      if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else if (p < .10) "†" else ""
    }
    safe_len <- function(x) if (is.null(x)) 0L else length(x)
    
    html_table <- function(df) {
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
        return(tags$em("Table unavailable."))
      }
      tags$table(
        class = "table table-striped table-condensed",
        tags$thead(tags$tr(lapply(names(df), function(nm) tags$th(nm)))),
        tags$tbody(
          lapply(seq_len(nrow(df)), function(i) {
            tags$tr(lapply(df[i, , drop = FALSE], function(cell) tags$td(HTML(as.character(cell)))))
          })
        )
      )
    }
    
    # ---------- sample sizes (safe)
    n_train <- suppressWarnings(as.integer(s$train_n))
    if (!is.finite(n_train) || n_train < 1) n_train <- tryCatch(length(residuals(fit)), error = function(e) 0L)
    
    n_test <- suppressWarnings(as.integer(s$test_n))
    if (!is.finite(n_test) || n_test < 0) n_test <- 0L
    
    N <- n_train + n_test
    
    # ---------- lag choice (safe)
    L_in <- suppressWarnings(as.integer(input$diag_lag))
    L <- if (is.finite(L_in) && L_in > 0) L_in else 12L
    
    # ============================================================
    # ---------- IC (safe + robust AICc fallback)
    # ============================================================
    AIC_val <- suppressWarnings(tryCatch(as.numeric(stats::AIC(fit)), error = function(e) NA_real_))
    BIC_val <- suppressWarnings(tryCatch(as.numeric(stats::BIC(fit)), error = function(e) NA_real_))
    
    AICc_val <- suppressWarnings(tryCatch(as.numeric(forecast::AICc(fit)), error = function(e) NA_real_))
    
    if (!is.finite(AICc_val) && is.finite(AIC_val)) {
      n_fit <- suppressWarnings(tryCatch(stats::nobs(fit), error = function(e) NA_integer_))
      if (!is.finite(n_fit) || n_fit <= 0) {
        r <- suppressWarnings(tryCatch(as.numeric(residuals(fit)), error = function(e) numeric(0)))
        n_fit <- sum(is.finite(r))
      }
      if (!is.finite(n_fit) || n_fit <= 0) n_fit <- n_train
      
      k <- suppressWarnings(tryCatch(length(stats::coef(fit)), error = function(e) NA_integer_))
      if (!is.finite(k) || k < 0) k <- 0L
      k <- k + 1L
      
      if (is.finite(n_fit) && n_fit > (k + 1)) {
        AICc_val <- AIC_val + (2 * k * (k + 1)) / (n_fit - k - 1)
      } else {
        AICc_val <- NA_real_
      }
    }
    
    ic_df <- data.frame(
      Criterion = c("AIC", "AICc", "BIC"),
      Value     = c(fmt_num_local(AIC_val, 2), fmt_num_local(AICc_val, 2), fmt_num_local(BIC_val, 2)),
      check.names = FALSE
    )
    
    # ---------- coefficients + significance (robust)
    est <- suppressWarnings(tryCatch(stats::coef(fit), error = function(e) NULL))
    V   <- suppressWarnings(tryCatch(stats::vcov(fit), error = function(e) NULL))
    if (is.null(V)) V <- suppressWarnings(tryCatch(fit$var.coef, error = function(e) NULL))
    
    coef_df <- NULL
    if (!is.null(est) && length(est) > 0) {
      est <- as.numeric(est)
      nm  <- names(stats::coef(fit))
      if (is.null(nm)) nm <- paste0("param_", seq_along(est))
      
      se <- rep(NA_real_, length(est))
      if (!is.null(V)) {
        dV <- tryCatch(diag(V), error = function(e) rep(NA_real_, length(est)))
        if (length(dV) == length(est)) se <- sqrt(pmax(dV, 0))
      }
      z  <- est / se
      p  <- 2 * stats::pnorm(-abs(z))
      
      coef_df <- data.frame(
        Term     = nm,
        Estimate = sprintf("%.6f", est),
        SE       = ifelse(is.finite(se), sprintf("%.6f", se), "NA"),
        `z/t`    = ifelse(is.finite(z),  sprintf("%.3f",  z),  "NA"),
        `p`      = ifelse(is.finite(p),  sprintf("%.3f",  p),  "NA"),
        Sig      = vapply(p, sig_stars, character(1)),
        check.names = FALSE
      )
    }
    
    n_sig <- if (!is.null(coef_df)) sum(suppressWarnings(as.numeric(coef_df$p)) < 0.05, na.rm = TRUE) else 0L
    
    # ---------- residuals (safe)
    res <- suppressWarnings(tryCatch(as.numeric(residuals(fit)), error = function(e) numeric(0)))
    res <- res[is.finite(res)]
    n_res <- length(res)
    
    fitdf <- if (!is.null(est)) length(est) else 0L
    
    lb_lag <- min(L, max(1L, floor(n_res / 3)))
    lb <- if (n_res >= 5) {
      tryCatch(stats::Box.test(res, lag = lb_lag, type = "Ljung-Box", fitdf = fitdf), error = function(e) NULL)
    } else NULL
    
    bp <- if (n_res >= 5) {
      tryCatch(stats::Box.test(res, lag = lb_lag, type = "Box-Pierce", fitdf = fitdf), error = function(e) NULL)
    } else NULL
    
    jb <- if (requireNamespace("tseries", quietly = TRUE) && n_res >= 5) {
      tryCatch(tseries::jarque.bera.test(res), error = function(e) NULL)
    } else NULL
    
    sw <- if (n_res >= 3 && n_res <= 5000) {
      tryCatch(stats::shapiro.test(res), error = function(e) NULL)
    } else NULL
    
    arch <- if (requireNamespace("FinTS", quietly = TRUE) && n_res >= 10) {
      tryCatch(FinTS::ArchTest(res, lags = min(12L, max(1L, floor(L / 2)))), error = function(e) NULL)
    } else NULL
    
    runs <- if (requireNamespace("tseries", quietly = TRUE) && n_res >= 10) {
      tryCatch(tseries::runs.test(res), error = function(e) NULL)
    } else NULL
    
    ad <- if (requireNamespace("nortest", quietly = TRUE) && n_res >= 8) {
      tryCatch(nortest::ad.test(res), error = function(e) NULL)
    } else NULL
    
    test_rows <- list()
    add_test <- function(name, stat, pval, note = "") {
      test_rows[[length(test_rows) + 1L]] <<- data.frame(
        Test = name,
        Statistic = if (is.null(stat)) "NA" else fmt_num_local(stat, 3),
        `p-value` = if (is.null(pval)) "NA" else fmt_p_local(pval),
        Interpretation = note,
        check.names = FALSE
      )
    }
    
    add_test(
      "Ljung–Box (residuals)",
      if (!is.null(lb)) unname(lb$statistic) else NULL,
      if (!is.null(lb)) unname(lb$p.value) else NULL,
      if (!is.null(lb)) {
        if (is.finite(lb$p.value) && lb$p.value >= 0.05) "No evidence of remaining autocorrelation (white-noise compatible)."
        else "Evidence of remaining autocorrelation (model may be under-specified)."
      } else "Not available."
    )
    
    add_test(
      "Box–Pierce (residuals)",
      if (!is.null(bp)) unname(bp$statistic) else NULL,
      if (!is.null(bp)) unname(bp$p.value) else NULL,
      if (!is.null(bp)) {
        if (is.finite(bp$p.value) && bp$p.value >= 0.05) "Consistent with uncorrelated residuals."
        else "Suggests residual autocorrelation."
      } else "Not available."
    )
    
    add_test(
      "Jarque–Bera (normality)",
      if (!is.null(jb)) unname(jb$statistic) else NULL,
      if (!is.null(jb)) unname(jb$p.value) else NULL,
      if (!is.null(jb)) {
        if (is.finite(jb$p.value) && jb$p.value >= 0.05) "No evidence against normality."
        else "Residuals deviate from normality (common in real series)."
      } else "Package 'tseries' missing or test unavailable."
    )
    
    add_test(
      "Shapiro–Wilk (normality)",
      if (!is.null(sw)) unname(sw$statistic) else NULL,
      if (!is.null(sw)) unname(sw$p.value) else NULL,
      if (!is.null(sw)) {
        if (is.finite(sw$p.value) && sw$p.value >= 0.05) "No evidence against normality."
        else "Evidence against normality."
      } else if (n_res > 5000) "Not computed (n > 5000)." else "Not available."
    )
    
    add_test(
      "ARCH LM (heteroskedasticity)",
      if (!is.null(arch)) unname(arch$statistic) else NULL,
      if (!is.null(arch)) unname(arch$p.value) else NULL,
      if (!is.null(arch)) {
        if (is.finite(arch$p.value) && arch$p.value >= 0.05) "No evidence of remaining ARCH effects."
        else "Evidence of ARCH effects → consider GARCH for variance."
      } else "Package 'FinTS' missing or test unavailable."
    )
    
    add_test(
      "Runs test (randomness)",
      if (!is.null(runs)) unname(runs$statistic) else NULL,
      if (!is.null(runs)) unname(runs$p.value) else NULL,
      if (!is.null(runs)) {
        if (is.finite(runs$p.value) && runs$p.value >= 0.05) "No evidence against randomness."
        else "Evidence of non-randomness (structure may remain)."
      } else "Package 'tseries' missing or test unavailable."
    )
    
    add_test(
      "Anderson–Darling (normality)",
      if (!is.null(ad)) unname(ad$statistic) else NULL,
      if (!is.null(ad)) unname(ad$p.value) else NULL,
      if (!is.null(ad)) {
        if (is.finite(ad$p.value) && ad$p.value >= 0.05) "No evidence against normality."
        else "Evidence against normality (sensitive in tails)."
      } else "Package 'nortest' missing or test unavailable."
    )
    
    tests_df <- if (length(test_rows)) do.call(rbind, test_rows) else data.frame()
    
    # ---------- forecast accuracy (safe)  [NOTE: this is on the model scale]
    acc_df <- NULL
    acc_sentence <- "No holdout test set was detected; therefore, out-of-sample accuracy was not computed."
    
    y_test <- s$ts_test
    has_test <- !is.null(y_test) && safe_len(y_test) > 0 && n_test > 0
    
    if (has_test) {
      y_test_num <- as.numeric(y_test)
      y_hat_num  <- suppressWarnings(tryCatch(as.numeric(fc$mean), error = function(e) rep(NA_real_, length(y_test_num))))
      h <- min(length(y_test_num), length(y_hat_num))
      if (h >= 1) {
        e <- y_test_num[seq_len(h)] - y_hat_num[seq_len(h)]
        rmse <- sqrt(mean(e^2, na.rm = TRUE))
        mae  <- mean(abs(e), na.rm = TRUE)
        mape <- mean(abs(e) / pmax(abs(y_test_num[seq_len(h)]), .Machine$double.eps), na.rm = TRUE)
        
        acc_df <- data.frame(
          Metric = c("RMSE", "MAE", "MAPE"),
          Value  = c(fmt_num_local(rmse, 3), fmt_num_local(mae, 3), paste0(fmt_num_local(100*mape, 2), "%")),
          check.names = FALSE
        )
        
        acc_sentence <- paste0(
          "Over the holdout period (test n = ", h, "), forecast performance was ",
          "RMSE = ", fmt_num_local(rmse, 3), ", ",
          "MAE = ", fmt_num_local(mae, 3), ", ",
          "MAPE = ", fmt_num_local(100*mape, 2), "%."
        )
      }
    }
    
    # ---------- horizon narrative
    horizon_txt <- if (has_test) {
      paste0("Validation mode was used: the forecast horizon was forced to match the test length (h = ", n_test, ").")
    } else {
      paste0("Future mode was used: forecasts were produced beyond the training sample (h = ", fc0$h, ").")
    }
    
    # ---------- model string
    season_txt <- suppressWarnings(as.integer(eq$s))
    s_txt <- if (is.finite(season_txt) && season_txt > 0) as.character(season_txt) else "s"
    
    model_str <- sprintf(
      "SARIMA(%d,%d,%d)(%d,%d,%d)[%s]",
      eq$p, eq$d, eq$q, eq$P, eq$D, eq$Q, s_txt
    )
    
    # ---------- diagnostic verdict
    lb_ok <- !is.null(lb) && is.finite(lb$p.value) && lb$p.value >= 0.05
    arch_ok <- is.null(arch) || (is.finite(arch$p.value) && arch$p.value >= 0.05)
    diag_verdict <- paste0(
      if (lb_ok) "Residual autocorrelation was not statistically detected (Ljung–Box p ≥ .05). "
      else "Residual autocorrelation may remain (Ljung–Box p < .05). ",
      if (arch_ok) "No clear evidence of residual ARCH effects was found (or test unavailable)."
      else "Residual ARCH effects were detected → a GARCH extension is recommended."
    )
    
    # ============================================================
    # ---------- Inverse transform equation + bias-adjusted back-forecasts
    # ============================================================
    tr_global <- input$transform %||% "none"
    
    # lambda consistent with your transform step:
    # apply_transform() uses forecast::BoxCox.lambda(y, lower=0) when input$lambda is NA
    lambda_used <- NA_real_
    if (identical(tr_global, "boxcox")) {
      p_obj <- tryCatch(prepared(), error = function(e) NULL)
      lam_in <- input$lambda
      if (is.null(lam_in) || (length(lam_in) == 1 && is.na(lam_in))) {
        lambda_used <- tryCatch(
          forecast::BoxCox.lambda(p_obj$df$y_filled, lower = 0),
          error = function(e) NA_real_
        )
      } else {
        lambda_used <- suppressWarnings(as.numeric(lam_in))
      }
    }
    
    # Forecasts are on TRANSFORMED scale:
    mu_t <- suppressWarnings(as.numeric(fc$mean))
    
    # se can be NULL / missing; handle safely
    se_t <- tryCatch(as.numeric(fc$se), error = function(e) rep(NA_real_, length(mu_t)))
    if (length(se_t) != length(mu_t)) se_t <- rep(NA_real_, length(mu_t))
    var_t <- se_t^2  # forecast variance on transformed scale (when available)
    
    lo_t <- tryCatch(fc$lower, error = function(e) NULL)
    hi_t <- tryCatch(fc$upper, error = function(e) NULL)
    lvl_names <- tryCatch(colnames(fc$lower), error = function(e) NULL)
    
    inv_equation_html <- NULL
    inv_note_html <- NULL
    
    inv_mean <- mu_t
    inv_lo   <- lo_t
    inv_hi   <- hi_t
    
    tex_or_plain <- function(tex, plain) {
      if (exists("tex_display", mode = "function")) HTML(tex_display(tex)) else HTML(plain)
    }
    
    if (identical(tr_global, "log")) {
      
      # Always-available back-transform (median)
      inv_mean <- exp(mu_t)
      
      # Bias-adjusted mean where variance is available
      ok <- is.finite(var_t)
      inv_mean[ok] <- exp(mu_t[ok] + 0.5 * var_t[ok])
      
      # Interval bounds: back-transform quantiles
      if (!is.null(lo_t) && !is.null(hi_t)) {
        inv_lo <- exp(lo_t)
        inv_hi <- exp(hi_t)
      }
      
      inv_equation_html <- tex_or_plain(
        "y = \\exp(z) \\quad \\text{where } z = \\ln(y)",
        "y = exp(z), where z = log(y)"
      )
      inv_note_html <- HTML(
        "<span style='font-size:12px;color:#444;'>
        Point forecasts use <b>exp(μ + 0.5·σ²)</b> when σ² is available; otherwise they fall back to <b>exp(μ)</b>.
        Interval bounds are back-transformed using <b>exp()</b>.
      </span>"
      )
      
    } else if (identical(tr_global, "boxcox")) {
      
      validate(need(is.finite(lambda_used), "Box–Cox lambda is missing/invalid; cannot compute inverse forecasts."))
      lam <- lambda_used
      
      # Always-available back-transform (median)
      inv_mean <- forecast::InvBoxCox(mu_t, lam)
      
      # Bias-adjusted mean where variance is available
      ok <- is.finite(var_t)
      inv_mean[ok] <- forecast::InvBoxCox(mu_t[ok], lam, biasadj = TRUE, var = var_t[ok])
      
      if (!is.null(lo_t) && !is.null(hi_t)) {
        inv_lo <- tryCatch(forecast::InvBoxCox(lo_t, lam), error = function(e) lo_t)
        inv_hi <- tryCatch(forecast::InvBoxCox(hi_t, lam), error = function(e) hi_t)
      }
      
      if (is.finite(lam) && abs(lam) < 1e-6) {
        inv_equation_html <- tex_or_plain(
          "y = \\exp(z) \\quad (\\lambda \\approx 0)",
          "y = exp(z)  (lambda ≈ 0)"
        )
      } else {
        inv_equation_html <- tex_or_plain(
          "y = (\\lambda z + 1)^{1/\\lambda}",
          "y = (lambda*z + 1)^(1/lambda)"
        )
      }
      
      inv_note_html <- HTML(paste0(
        "<span style='font-size:12px;color:#444;'>
        λ used = <b>", if (is.finite(lambda_used)) formatC(lambda_used, digits = 4, format = "f") else "NA", "</b>.
        Point forecasts use <b>InvBoxCox(μ, λ, biasadj=TRUE, var=se²)</b> when se² is available; otherwise they fall back to <b>InvBoxCox(μ, λ)</b>.
        Interval bounds are back-transformed using <b>InvBoxCox()</b>.
      </span>"
      ))
    }
    
    # Build an ORIGINAL-SCALE forecast table (bias-adjusted mean + back-transformed intervals)
    fc_orig_df <- data.frame(
      Horizon = seq_along(inv_mean),
      Mean_original = as.numeric(inv_mean),
      stringsAsFactors = FALSE
    )
    
    if (!is.null(inv_lo) && !is.null(inv_hi) && !is.null(lvl_names)) {
      for (j in seq_along(lvl_names)) {
        lvl <- lvl_names[j]
        fc_orig_df[[paste0("Lo_", lvl)]] <- as.numeric(inv_lo[, j])
        fc_orig_df[[paste0("Hi_", lvl)]] <- as.numeric(inv_hi[, j])
      }
    }
    
    # format numeric columns for HTML display
    fc_orig_df_fmt <- fc_orig_df
    for (nm in names(fc_orig_df_fmt)) {
      if (is.numeric(fc_orig_df_fmt[[nm]])) {
        fc_orig_df_fmt[[nm]] <- ifelse(
          is.finite(fc_orig_df_fmt[[nm]]),
          sprintf("%.6f", fc_orig_df_fmt[[nm]]),
          "NA"
        )
      }
    }
    
    # ---------- build report UI (NO output$ definitions here)
    tagList(
      tags$h3("Manual SARIMA: Full academic conclusion (report-ready)"),
      
      # 1. Objective and modelling rationale
      tags$hr(), tags$br(),
      tags$h4(tags$strong("1. Objective and modelling rationale")),
      tags$p(
        "A manually specified seasonal ARIMA (SARIMA) model was estimated to represent linear temporal dependence,",
        " including seasonal structure, and to provide an interpretable baseline for forecasting."
      ),
      
      # 2. Data design and sample split
      tags$hr(), tags$br(),
      tags$h4(tags$strong("2. Data design and sample split")),
      tags$p(HTML(paste0(
        "The analysis used <b>N = ", N, "</b> observations (training <b>n = ", n_train, "</b>",
        if (has_test) paste0(", test <b>n = ", n_test, "</b>") else "",
        ")."
      ))),
      tags$p(tags$b("Forecast design. "), horizon_txt),
      
      # 3. Identification visuals (time series + ACF/PACF)
      tags$hr(), tags$br(),
      tags$h4(tags$strong("3. Identification visuals (time series + ACF/PACF)")),
      tags$p(
        "The plots below summarize the observed series (with the train/test split, if applicable), ",
        "followed by ACF/PACF for the training series and for the differenced series implied by the chosen (d, D, s)."
      ),
      
      tags$br(),
      tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure A. Time series with split marker")),
      plotOutput("manual_report_ts_plot", height = 360),
      
      # 4. Stationarity assessment (ADF, KPSS, and Phillips–Perron)
      tags$hr(), tags$br(),
      tags$h4(tags$strong("4. Stationarity assessment (ADF, KPSS, and Phillips–Perron)")),
      tags$p("Stationarity tests were applied to the training series to evaluate whether differencing is required before SARIMA identification and estimation."),
      tags$hr(),
      uiOutput("manual_report_stationarity"),
      
      # 5. Transformed series: differencing and seasonality
      tags$hr(), tags$br(),
      tags$h4(tags$strong("5. Transformed series: differencing and seasonality")),
      
      tags$br(),
      tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure B. Seasonal subseries")),
      plotOutput("manual_report_subseries", height = 360),
      
      tags$hr(), tags$br(),
      tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure C. Seasonal box-plot")),
      plotOutput("manual_report_seasonal_box", height = 360),
      
      tags$hr(), tags$br(),
      tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure D. Transformed training series and differenced (d, D, s) series")),
      plotOutput("manual_report_ts_trans_and_diff", height = 360),
      
      tags$hr(), tags$br(),
      tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure E. ACF and PACF (training series)")),
      fluidRow(
        column(6, plotOutput("manual_report_acf",  height = 280)),
        column(6, plotOutput("manual_report_pacf", height = 280))
      ),
      
      tags$hr(), tags$br(),
      tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Figure F. ACF and PACF (modified / differenced series using current d, D, s)")),
      fluidRow(
        column(6, plotOutput("manual_report_acf_mod",  height = 280)),
        column(6, plotOutput("manual_report_pacf_mod", height = 280))
      ),
      
      tags$hr(), tags$br(),
      uiOutput("manual_report_stationarity_mod"),
      
      # 6. Final model specification and fit quality
      tags$hr(), tags$br(),
      tags$h4(tags$strong("6. Final model specification and fit quality")),
      tags$p(HTML(paste0(
        "The final manual specification was <b>", model_str, "</b>",
        if (isTRUE(input$manual_drift)) " including drift/mean." else " without drift/mean.",
        " Model adequacy was assessed using information criteria, coefficient inference, residual diagnostics, and forecast performance."
      ))),
      
      tags$hr(), tags$br(),
      tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Table A. Goodness-of-fit (information criteria)")),
      html_table(ic_df),
      
      tags$hr(), tags$br(),
      tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Table B. Parameter estimates and significance (approx. z/t tests)")),
      if (!is.null(coef_df)) html_table(coef_df) else tags$em("No coefficients available."),
      tags$p(HTML(paste0(
        "In total, <b>", n_sig, "</b> parameter(s) were significant at α = .05 (marked by *, **, ***)."
      ))),
      
      # 7. Model equations (replication-ready)
      tags$hr(), tags$br(),
      tags$h4(tags$strong("7. Model equations (replication-ready)")),
      tags$p(
        "The fitted model is reported below in operator notation (general form), expanded form, and the numerical equation",
        " based on the estimated parameters."
      ),
      tags$div(
        style = "padding:10px;border:1px solid #e5e5e5;border-radius:6px;background:#fcfcfc;",
        
        tags$br(), tags$hr(), tags$hr(),
        tags$h5("General SARIMA formulation"),
        HTML(eq$eq_general),
        
        tags$br(), tags$hr(), tags$hr(),
        tags$h5("Expanded operator form"),
        HTML(eq$eq_expanded),
        
        tags$br(), tags$hr(), tags$hr(),
        tags$h5("Numerical model"),
        HTML(eq$eq_line3),
        
        tags$br(), tags$hr(), tags$hr(),
        HTML(eq$eq_line4),
        
        tags$br(), tags$hr(), tags$hr()
      ),
      
      # 8. Residual diagnostics (graphical evidence)
      tags$hr(), tags$br(),
      tags$h4(tags$strong("8. Residual diagnostics (graphical evidence)")),
      tags$p(
        "Graphical diagnostics evaluate whether residuals resemble white noise (no systematic autocorrelation),",
        " approximate normality (Q–Q and histogram), and stable variance."
      ),
      
      fluidRow(
        column(6, plotOutput("manual_resid_ts",   height = 220)),
        column(6, plotOutput("manual_resid_acf",  height = 220))
      ),
      fluidRow(
        column(6, plotOutput("manual_resid_hist", height = 220)),
        column(6, plotOutput("manual_resid_qq",   height = 220))
      ),
      
      tags$h5("Ljung–Box p-values by lag"),
      plotOutput("manual_resid_lb_pvals", height = 260),
      
      # 9. Residual tests (formal inference)
      tags$hr(), tags$br(),
      tags$h4(tags$strong("9. Residual tests (formal inference)")),
      tags$p(
        "Formal tests complement the plots: Ljung–Box/Box–Pierce assess remaining autocorrelation;",
        " Jarque–Bera/Shapiro–Wilk/Anderson–Darling assess normality;",
        " ARCH LM tests conditional heteroskedasticity; the runs test checks randomness."
      ),
      
      tags$hr(), tags$br(),
      tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Table C. Residual test summary")),
      html_table(tests_df),
      tags$p(tags$b("Diagnostic synthesis. "), diag_verdict),
      
      # 10. Forecasting results and predictive performance
      tags$hr(), tags$br(),
      tags$h4(tags$strong("10. Forecasting results and predictive performance")),
      tags$p(acc_sentence),
      
      if (!is.null(acc_df)) tagList(
        tags$h5("Table D. Holdout accuracy (test set)"),
        html_table(acc_df)
      ) else NULL,
      
      tags$hr(),  tags$br(),
      tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Forecast plot")),
      plotOutput("manual_forecast_plot", height = 420),
      
      tags$hr(), tags$br(),
      tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Forecast table")),
      tableOutput("manual_forecast_table"),
      
      tags$hr(), tags$br(),
      tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Accuracy table (your app output)")),
      tableOutput("manual_accuracy_table"),
      
      # ---------- Back-transformed forecasts (original scale) + plot placeholder
      tags$hr(), tags$br(),
      tags$h4(tags$strong("10.B Back-transformed forecasts (original scale)")),
      if (!identical(tr_global, "none")) tagList(
        tags$p("Because the model was estimated on the transformed series, forecasts are reported below on the original measurement scale using an appropriate bias adjustment for the point forecasts."),
        tags$div(
          style = "padding:10px;border:1px solid #e5e5e5;border-radius:6px;background:#fcfcfc;",
          tags$h5("Inverse transform equation (simple form)"),
          inv_equation_html,
          tags$br(),
          inv_note_html
        ),
        
        tags$hr(), tags$br(),
        tags$h5(strong(" \u00A0\u00A0 \u25A0 \u00A0 Original-scale plot (observed + forecast + CIs)")),
        # NOTE: you must define output$manual_forecast_plot_original <- renderPlot(...) elsewhere in server
        plotOutput("manual_forecast_plot_original", height = 420),
        
        tags$hr(), tags$br(),
        tags$h5("Back-transformed forecast table (original scale)"),
        html_table(fc_orig_df_fmt)
      ) else tags$em("No transformation was applied; forecasts are already on the original scale."),
      
      # 11. Final conclusion (academic)
      tags$hr(), tags$br(),
      tags$h4(tags$strong("11. Final conclusion")),
      tags$p(
        "Overall, the manually specified SARIMA model provides a coherent and interpretable representation of seasonal linear dynamics,",
        " supported by information criteria, statistically interpretable parameters, and diagnostic checks.",
        " When diagnostics indicate remaining autocorrelation, refinement should prioritize revising differencing (d, D) and AR/MA orders guided by ACF/PACF and Ljung–Box.",
        " When conditional heteroskedasticity is detected (ARCH LM), a volatility model (e.g., GARCH) should be added to the mean equation to better represent time-varying variance."
      ),
      tags$p(
        "For reporting, the results above provide the full chain of evidence typically expected in academic manuscripts:",
        " (i) specification + IC, (ii) parameter inference, (iii) equation reporting, (iv) residual validation with plots and tests, and (v) forecasting with accuracy assessment."
      ),
      
      tags$br(), tags$hr(), tags$br(), tags$br()
    )
  })
  
  
  
  
  
  
  
  
  # IMPORTANT: renderUI wrapper + MathJax re-typeset (this is what makes equations show correctly)
  output$manual_conclusion_full <- renderUI({
    ui <- manual_conclusion_full_obj()
    
    # Re-typeset MathJax after the UI is inserted into the DOM
    session$onFlushed(function() {
      session$sendCustomMessage("mathjax-typeset", "manual_conclusion_box")
    }, once = TRUE)
    
    ui
  })
  
  
  
  
  
  # manual_conclusion_full_obj <- eventReactive(input$fit_manual, {
  #   req(manual_fit(), manual_fc(), manual_equations(), ts_train_test())
  #   
  #   fit <- manual_fit()
  #   fc0 <- manual_fc()
  #   fc  <- fc0$fc
  #   eq  <- manual_equations()
  #   s   <- ts_train_test()
  #   
  #   # ---- safe values
  #   n_train <- suppressWarnings(as.integer(s$train_n)); if (!is.finite(n_train)) n_train <- length(residuals(fit))
  #   n_test  <- suppressWarnings(as.integer(s$test_n));  if (!is.finite(n_test))  n_test  <- 0L
  #   N <- n_train + n_test
  #   
  #   L_in <- suppressWarnings(as.integer(input$diag_lag))
  #   L <- if (is.finite(L_in) && L_in > 0) L_in else 12L
  #   fitdf <- length(coef(fit))
  #   
  #   # ---- IC
  #   AIC_val  <- suppressWarnings(as.numeric(fit$aic))
  #   AICc_val <- suppressWarnings(as.numeric(fit$aicc))
  #   BIC_val  <- suppressWarnings(as.numeric(fit$bic))
  #   
  #   # ---- residual test (minimal)
  #   res <- as.numeric(residuals(fit))
  #   res <- res[is.finite(res)]
  #   lb <- tryCatch(Box.test(res, lag = min(L, max(1L, floor(length(res) / 3))), type = "Ljung-Box", fitdf = fitdf),
  #                  error = function(e) NULL)
  #   
  #   # ---- accuracy
  #   acc_line <- NULL
  #   if (n_test > 0) {
  #     acc <- tryCatch(accuracy_df(s$ts_test, fc$mean), error = function(e) NULL)
  #     if (!is.null(acc) && all(c("Metric", "Value") %in% names(acc))) {
  #       rmse <- acc$Value[acc$Metric == "RMSE"][1]
  #       mae  <- acc$Value[acc$Metric == "MAE"][1]
  #       mape <- acc$Value[acc$Metric == "MAPE"][1]
  #       acc_line <- tags$p(
  #         tags$b("Forecast accuracy (test set). "),
  #         HTML(paste0(
  #           "Over the holdout period (n = ", n_test, "), performance was RMSE = ",
  #           fmt_num(rmse, 3), ", MAE = ", fmt_num(mae, 3),
  #           if (is.finite(mape)) paste0(", MAPE = ", fmt_num(100 * mape, 2), "%") else "",
  #           "."
  #         ))
  #       )
  #     }
  #   }
  #   if (is.null(acc_line)) {
  #     acc_line <- tags$p(tags$b("Forecast accuracy. "),
  #                        "No holdout test set was detected; therefore, out-of-sample accuracy was not computed.")
  #   }
  #   
  #   # ---- horizon narrative
  #   horizon_txt <- if (n_test > 0) {
  #     paste0("Validation mode was used: the forecast horizon was forced to match the test length (h = ", n_test, ").")
  #   } else {
  #     paste0("Future mode was used: forecasts were produced beyond the training sample (h = ", fc0$h, ").")
  #   }
  #   
  #   # ---- model text (manual)
  #   season_txt <- if (is.finite(eq$s)) eq$s else NA_integer_
  #   model_str <- sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)[%s]",
  #                        eq$p, eq$d, eq$q, eq$P, eq$D, eq$Q,
  #                        if (is.finite(season_txt)) as.character(season_txt) else "s")
  #   
  #   tagList(
  #     tags$h3("Manual SARIMA: Full academic conclusion"),
  #     
  #     tags$h4("1. Rationale and modelling objective"),
  #     tags$p(
  #       "A manually specified seasonal ARIMA (SARIMA) model was estimated to provide explicit control over non-seasonal and seasonal dynamics. ",
  #       "This approach is appropriate when domain knowledge and diagnostic patterns (ACF/PACF after differencing) motivate targeted structure beyond automated search."
  #     ),
  #     
  #     tags$h4("2. Sample design"),
  #     tags$p(HTML(paste0(
  #       "The analysis used <b>N = ", N, "</b> observations (training <b>n = ", n_train, "</b>",
  #       if (n_test > 0) paste0(", test <b>n = ", n_test, "</b>") else "",
  #       ")."
  #     ))),
  #     tags$p(tags$b("Forecast design. "), horizon_txt),
  #     
  #     tags$h4("3. Final specification and fit"),
  #     tags$p(HTML(paste0(
  #       "The fitted manual specification was <b>", model_str, "</b>",
  #       if (isTRUE(input$manual_drift)) " with drift/mean." else " without drift/mean.",
  #       " The corresponding fitted object was reported as <b>", as.character(fit), "</b>."
  #     ))),
  #     tags$ul(
  #       tags$li(HTML(paste0("AIC = <b>", fmt_num(AIC_val, 2), "</b>"))),
  #       tags$li(HTML(paste0("AICc = <b>", fmt_num(AICc_val, 2), "</b>"))),
  #       tags$li(HTML(paste0("BIC = <b>", fmt_num(BIC_val, 2), "</b>")))
  #     ),
  #     
  #     tags$h4("4. Model equations (reporting-ready)"),
  #     tags$p(
  #       "For academic reporting and replication, the fitted model is expressed in standard operator notation, followed by an expanded form ",
  #       "and a numerical representation using the estimated coefficients."
  #     ),
  #     tags$div(
  #       style = "padding:10px;border:1px solid #e5e5e5;border-radius:6px;background:#fcfcfc;text-align:left;",
  #       tags$h5("General SARIMA formulation"),
  #       HTML(eq$eq_general),
  #       tags$hr(),
  #       tags$h5("Expanded operator form"),
  #       HTML(eq$eq_expanded),
  #       tags$hr(),
  #       tags$h5("Numerical model"),
  #       HTML(eq$eq_line3),
  #       tags$hr(),
  #       HTML(eq$eq_line4)
  #     ),
  #     
  #     tags$h4("5. Residual diagnostics (adequacy of linear dynamics)"),
  #     tags$p(
  #       "Adequacy was evaluated using residual plots and formal tests. A key criterion is that residuals resemble white noise, ",
  #       "indicating that the model has captured the systematic linear dependence."
  #     ),
  #     if (!is.null(lb)) {
  #       tags$p(HTML(paste0(
  #         "<b>Ljung–Box test:</b> Q(", lb$parameter, ") = ", fmt_num(lb$statistic, 3),
  #         ", p ", fmt_p(lb$p.value), "."
  #       )))
  #     } else {
  #       tags$p(tags$b("Ljung–Box test:"), " unavailable (insufficient residuals or test error).")
  #     },
  #     
  #     tags$h4("6. Forecasting and predictive performance"),
  #     acc_line,
  #     
  #     tags$h4("7. Overall conclusion and recommended next steps"),
  #     tags$p(
  #       "In sum, the manual SARIMA specification provides an interpretable representation of linear dependence, contingent on residual whiteness. ",
  #       "If residual autocorrelation persists, revise differencing (d, D) or adjust AR/MA orders guided by diagnostics. ",
  #       "If volatility clustering is evident, consider modelling conditional variance (e.g., GARCH) alongside the SARIMA mean equation."
  #     )
  #   )
  # })
  
  output$manual_conclusion_full <- renderUI({
    validate(need(input$fit_manual > 0, "Click “Fit” in the Manual tab to generate the full academic conclusion."))
    req(manual_conclusion_full_obj())
    
    session$onFlushed(function() {
      session$sendCustomMessage("mathjax-typeset", "manual_conclusion_box")
    }, once = TRUE)
    
    manual_conclusion_full_obj()
  })
  
  
  
  
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  #================================================================================================
  
  
  
  
  
  

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
  # output$tsPlot_Choice <- renderPlot({
  #   req(myData_Choice())
  #   req(input$plot_type_choice)
  #   
  #   ts_obj <- myData_Choice()
  #   p      <- prepared()  # for x-axis label
  #   freq   <- tryCatch(stats::frequency(ts_obj), error = function(e) NA_real_)
  #   
  #   df_ts <- function(z) {
  #     if (inherits(z, "ts")) {
  #       data.frame(t = as.numeric(stats::time(z)), y = as.numeric(z))
  #     } else {
  #       data.frame(t = seq_along(z), y = as.numeric(z))
  #     }
  #   }
  #   
  #   k_ma  <- max(2L, as.integer(input$ma_k %||% 5))
  #   lag_m <- max(1L, as.integer(input$lag_m %||% 12))
  #   
  #   plt <- switch(
  #     input$plot_type_choice,
  #     
  #     "Line" = {
  #       forecast::autoplot(
  #         ts_obj,
  #         size   = 1,
  #         colour = input$ts_line_color %||% "#2C7FB8"
  #       ) +
  #         ggplot2::labs(title = "Transformed series", x = p$x_label, y = "Transformed value")
  #     },
  #     
  #     "Points" = {
  #       d <- df_ts(ts_obj)
  #       ggplot2::ggplot(d, ggplot2::aes(t, y)) +
  #         ggplot2::geom_point(size = 1) +
  #         ggplot2::labs(title = "Points", x = p$x_label, y = "Transformed value")
  #     },
  #     
  #     "Line + Points" = {
  #       d <- df_ts(ts_obj)
  #       ggplot2::ggplot(d, ggplot2::aes(t, y)) +
  #         ggplot2::geom_line() +
  #         ggplot2::geom_point(size = 0.9, alpha = 0.8) +
  #         ggplot2::labs(title = "Line + Points", x = p$x_label, y = "Transformed value")
  #     },
  #     
  #     "Smoothed (LOESS)" = {
  #       d <- df_ts(ts_obj)
  #       ggplot2::ggplot(d, ggplot2::aes(t, y)) +
  #         ggplot2::geom_line(alpha = 0.4) +
  #         ggplot2::geom_smooth(method = "loess", se = FALSE, span = 0.2) +
  #         ggplot2::labs(title = "LOESS smooth", x = p$x_label, y = "Transformed value")
  #     },
  #     
  #     "Moving average" = {
  #       d <- df_ts(ts_obj)
  #       ma <- stats::filter(d$y, rep(1/k_ma, k_ma), sides = 2)
  #       d$ma <- as.numeric(ma)
  #       ggplot2::ggplot(d, ggplot2::aes(t, y)) +
  #         ggplot2::geom_line(alpha = 0.4) +
  #         ggplot2::geom_line(ggplot2::aes(y = ma), size = 1) +
  #         ggplot2::labs(title = sprintf("Moving average (k = %d)", k_ma), x = p$x_label, y = "Transformed value")
  #     },
  #     
  #     "Cumulative sum" = {
  #       d <- df_ts(ts_obj); d$cum <- cumsum(d$y)
  #       ggplot2::ggplot(d, ggplot2::aes(t, cum)) +
  #         ggplot2::geom_line() +
  #         ggplot2::labs(title = "Cumulative sum", x = p$x_label, y = "Cumulative value")
  #     },
  #     
  #     "Seasonal plot" = {
  #       validate(need(is.finite(freq) && freq > 1, "Seasonal plot requires frequency > 1."))
  #       forecast::ggseasonplot(ts_obj, year.labels = TRUE) +
  #         ggplot2::labs(title = "Seasonal plot", x = p$x_label, y = "Value")
  #     },
  #     
  #     "Seasonal subseries" = {
  #       validate(need(is.finite(freq) && freq > 1, "Seasonal subseries requires frequency > 1."))
  #       forecast::ggsubseriesplot(ts_obj) +
  #         ggplot2::labs(title = "Seasonal subseries", x = p$x_label, y = "Value")
  #     },
  #     
  #     "Polar seasonal" = {
  #       validate(need(is.finite(freq) && freq > 1, "Polar seasonal requires frequency > 1."))
  #       forecast::ggseasonplot(ts_obj, polar = TRUE) +
  #         ggplot2::labs(title = "Polar seasonal plot", x = p$x_label, y = "Value")
  #     },
  #     
  #     "Seasonal boxplot" = {
  #       validate(need(is.finite(freq) && freq > 1, "Seasonal boxplot requires frequency > 1."))
  #       d <- if (inherits(ts_obj, "ts")) data.frame(season = stats::cycle(ts_obj), y = as.numeric(ts_obj))
  #       else data.frame(season = factor(1), y = as.numeric(ts_obj))
  #       ggplot2::ggplot(d, ggplot2::aes(x = factor(season), y = y)) +
  #         ggplot2::geom_boxplot() +
  #         ggplot2::labs(title = "Seasonal boxplot", x = "Season", y = "Value")
  #     },
  #     
  #     "Classical decomposition (additive)" = {
  #       validate(need(is.finite(freq) && freq > 1, "Classical decomposition requires frequency > 1."))
  #       dc <- stats::decompose(ts_obj, type = "additive")
  #       forecast::autoplot(dc) + ggplot2::labs(title = "Classical decomposition (additive)")
  #     },
  #     
  #     "Classical decomposition (multiplicative)" = {
  #       validate(need(is.finite(freq) && freq > 1, "Classical decomposition requires frequency > 1."))
  #       dc <- stats::decompose(ts_obj, type = "multiplicative")
  #       forecast::autoplot(dc) + ggplot2::labs(title = "Classical decomposition (multiplicative)")
  #     },
  #     
  #     "STL decomposition" = {
  #       validate(need(is.finite(freq) && freq > 1, "STL decomposition requires frequency > 1."))
  #       decomp <- stats::stl(ts_obj, s.window = "periodic")
  #       forecast::autoplot(decomp) + ggplot2::labs(title = "STL decomposition")
  #     },
  #     
  #     "Histogram" = {
  #       xx <- as.numeric(stats::na.omit(ts_obj))
  #       ggplot2::ggplot(data.frame(x = xx), ggplot2::aes(x)) +
  #         ggplot2::geom_histogram(bins = 30) +
  #         ggplot2::labs(title = "Histogram", x = "Value", y = "Count")
  #     },
  #     
  #     "Density" = {
  #       xx <- as.numeric(stats::na.omit(ts_obj))
  #       ggplot2::ggplot(data.frame(x = xx), ggplot2::aes(x)) +
  #         ggplot2::geom_density() +
  #         ggplot2::labs(title = "Density", x = "Value", y = "Density")
  #     },
  #     
  #     "QQ plot" = {
  #       xx <- as.numeric(stats::na.omit(ts_obj))
  #       ggplot2::ggplot(data.frame(x = xx), ggplot2::aes(sample = x)) +
  #         ggplot2::stat_qq() +
  #         ggplot2::stat_qq_line() +
  #         ggplot2::labs(title = "Normal Q-Q plot", x = "Theoretical quantiles", y = "Sample quantiles")
  #     },
  #     
  #     "Lag-1 scatter" = {
  #       xx <- as.numeric(stats::na.omit(ts_obj))
  #       validate(need(length(xx) >= 2, "Not enough data for lag-1 scatter."))
  #       d <- data.frame(x = xx[-length(xx)], y = xx[-1])
  #       ggplot2::ggplot(d, ggplot2::aes(x, y)) +
  #         ggplot2::geom_point() +
  #         ggplot2::geom_smooth(method = "lm", se = FALSE) +
  #         ggplot2::labs(title = "Lag-1 scatter", x = "y(t-1)", y = "y(t)")
  #     },
  #     
  #     "Lag plot (1..m)" = {
  #       xx <- as.numeric(stats::na.omit(ts_obj))
  #       validate(need(length(xx) > (lag_m + 1), "Increase data or reduce m for lag plots."))
  #       forecast::gglagplot(ts_obj, lags = lag_m) +
  #         ggplot2::labs(title = sprintf("Lag plot (1..%d)", lag_m))
  #     },
  #     
  #     "ACF" = {
  #       forecast::ggAcf(ts_obj) + ggplot2::labs(title = "ACF")
  #     },
  #     
  #     "PACF" = {
  #       forecast::ggPacf(ts_obj) + ggplot2::labs(title = "PACF")
  #     },
  #     
  # 
  #     # inside your switch(input$plot_type_choice, ...)
  #     "ACF+PACF" = {
  #       # ACF (top) & PACF (bottom)
  #       p_acf  <- forecast::ggAcf(ts_obj)  + ggplot2::labs(title = "ACF")
  #       p_pacf <- forecast::ggPacf(ts_obj) + ggplot2::labs(title = "PACF")
  #       
  #       # Apply your theme helper
  #       p_acf  <- add_theme(p_acf)
  #       p_pacf <- add_theme(p_pacf)
  #       
  #       # Vertical layout: ACF on top, PACF below
  #       gridExtra::grid.arrange(p_acf, p_pacf, ncol = 1, heights = c(1, 1))
  #     },
  #     
  #   
  # 
  #     # library(gridExtra)
  #     
  #     "Time + ACF+PACF" = {
  #       # Time plot (top)
  #       p_time <- forecast::autoplot(
  #         ts_obj,
  #         size   = 1,
  #         colour = input$ts_line_color %||% "#2C7FB8"
  #       ) +
  #         ggplot2::labs(
  #           title = "Time plot",
  #           x = p$x_label,
  #           y = "Transformed value"
  #         )
  # 
  #       # ACF (bottom-left) & PACF (bottom-right)
  #       p_acf  <- forecast::ggAcf(ts_obj)  + ggplot2::labs(title = "ACF")
  #       p_pacf <- forecast::ggPacf(ts_obj) + ggplot2::labs(title = "PACF")
  # 
  #       # Apply your selected theme to each plot
  #       p_time <- add_theme(p_time)
  #       p_acf  <- add_theme(p_acf)
  #       p_pacf <- add_theme(p_pacf)
  # 
  #       # Bottom row: ACF (left) | PACF (right)
  #       bottom_row <- gridExtra::arrangeGrob(p_acf, p_pacf, ncol = 2)
  # 
  #       # Final layout: Time on top; ACF+PACF in one row at bottom
  #       gridExtra::grid.arrange(p_time, bottom_row, ncol = 1, heights = c(1.3, 1))
  #     }
  #     
  #     
  # 
  #     "Periodogram" = {
  #       xx <- as.numeric(stats::na.omit(ts_obj))
  #       validate(need(length(xx) > 5, "More data needed for periodogram."))
  #       sp <- stats::spec.pgram(xx, detrend = TRUE, taper = 0.1, plot = FALSE)
  #       d  <- data.frame(freq = sp$freq, spec = sp$spec)
  #       ggplot2::ggplot(d, ggplot2::aes(freq, spec)) +
  #         ggplot2::geom_line() +
  #         ggplot2::labs(title = "Periodogram", x = "Frequency", y = "Spectral density")
  #     }
  #   )
  #   
  #   # apply theme (for ggplot outputs)
  #   if (inherits(plt, "ggplot")) {
  #     plt <- add_theme(plt)
  #   }
  #   plt
  # }, res = 96)
  # 
  
  
  
 
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
          ggplot2::labs(
            title = sprintf("Moving average (k = %d)", k_ma),
            x = p$x_label, y = "Transformed value"
          )
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
      
      "ACF+PACF" = {
        p_acf  <- add_theme(forecast::ggAcf(ts_obj)  + ggplot2::labs(title = "ACF"))
        p_pacf <- add_theme(forecast::ggPacf(ts_obj) + ggplot2::labs(title = "PACF"))
        
        # patchwork vertical
        (p_acf / p_pacf) + patchwork::plot_layout(heights = c(1, 1))
      },
      
      
      
      # "Time + ACF+PACF" = {
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
      #   p_acf  <- forecast::ggAcf(ts_obj)  + ggplot2::labs(title = "ACF")
      #   p_pacf <- forecast::ggPacf(ts_obj) + ggplot2::labs(title = "PACF")
      #   
      #   p_time <- add_theme(p_time)
      #   p_acf  <- add_theme(p_acf)
      #   p_pacf <- add_theme(p_pacf)
      #   
      #   # patchwork: time on top, acf|pacf on bottom
      #   (p_time / (p_acf | p_pacf)) + patchwork::plot_layout(heights = c(1.3, 1))
      # },  
      
      
      # "Time + ACF+PACF" = {
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
      #   p_acf  <- forecast::ggAcf(ts_obj)  + ggplot2::labs(title = "ACF")
      #   p_pacf <- forecast::ggPacf(ts_obj) + ggplot2::labs(title = "PACF")
      #   
      #   p_time <- add_theme(p_time)
      #   p_acf  <- add_theme(p_acf)
      #   p_pacf <- add_theme(p_pacf)
      #   
      #   # turn ggplots into grobs
      #   g_time <- ggplot2::ggplotGrob(p_time)
      #   g_acf  <- ggplot2::ggplotGrob(p_acf)
      #   g_pacf <- ggplot2::ggplotGrob(p_pacf)
      #   
      #   # bottom row: 2 columns
      #   g_bottom <- gridExtra::arrangeGrob(g_acf, g_pacf, ncol = 2)
      #   
      #   # full layout: top then bottom; heights control proportion
      #   g_all <- gridExtra::arrangeGrob(
      #     g_time, g_bottom,
      #     ncol = 1,
      #     heights = grid::unit(c(1.3, 1), "null")
      #   )
      #   
      #   grid::grid.newpage()
      #   grid::grid.draw(g_all)
      # },
      
      
      # "Time + ACF+PACF" = {
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
      #   p_acf  <- forecast::ggAcf(ts_obj)  + ggplot2::labs(title = "ACF")
      #   p_pacf <- forecast::ggPacf(ts_obj) + ggplot2::labs(title = "PACF")
      #   
      #   p_time <- add_theme(p_time)
      #   p_acf  <- add_theme(p_acf)
      #   p_pacf <- add_theme(p_pacf)
      #   
      #   g_time <- ggplot2::ggplotGrob(p_time)
      #   g_acf  <- ggplot2::ggplotGrob(p_acf)
      #   g_pacf <- ggplot2::ggplotGrob(p_pacf)
      #   
      #   g_bottom <- gridExtra::arrangeGrob(g_acf, g_pacf, ncol = 2)
      #   
      #   # empty middle row
      #   g_spacer <- grid::nullGrob()
      #   
      #   g_all <- gridExtra::arrangeGrob(
      #     g_time,
      #     g_spacer,
      #     g_bottom,
      #     ncol = 1,
      #     heights = grid::unit(c(2, 0.1, 1.2), "null")  # adjust 0.15 to control empty space
      #   )
      #   
      #   grid::grid.newpage()
      #   grid::grid.draw(g_all)
      # },
      
      
      
      
      "Time + ACF+PACF" = {
        # Top: time plot
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

        # Bottom: ACF + PACF
        p_acf  <- forecast::ggAcf(ts_obj)  + ggplot2::labs(title = "ACF")
        p_pacf <- forecast::ggPacf(ts_obj) + ggplot2::labs(title = "PACF")

        # Apply theme
        p_time <- add_theme(p_time)
        p_acf  <- add_theme(p_acf)
        p_pacf <- add_theme(p_pacf)

        # Patchwork layout (resizes correctly, no overlap)
        (p_time / (p_acf | p_pacf)) +
          patchwork::plot_layout(heights = c(1.3, 1))
      },
      
      
      
      "Periodogram" = {
        xx <- as.numeric(stats::na.omit(ts_obj))
        validate(need(length(xx) >= 8, "Need at least 8 observations for a periodogram."))
        
        # use your UI taper slider if it exists; otherwise 0.1
        taper <- suppressWarnings(as.numeric(input$stp_spec_taper %||% 0.1))
        if (!is.finite(taper)) taper <- 0.1
        
        sp <- tryCatch(
          stats::spec.pgram(xx, detrend = TRUE, taper = taper, plot = FALSE),
          error = function(e) e
        )
        validate(need(!inherits(sp, "error"), paste("spec.pgram failed:", sp$message)))
        
        d <- data.frame(freq = sp$freq, spec = as.numeric(sp$spec))
        d <- d[is.finite(d$freq) & is.finite(d$spec), , drop = FALSE]
        validate(need(nrow(d) > 1, "Periodogram returned no finite values (check data/transform)."))
        
        ggplot2::ggplot(d, ggplot2::aes(freq, spec)) +
          ggplot2::geom_line(linewidth = 1, colour = input$ts_line_color %||% "#2C7FB8") +
          ggplot2::labs(title = "Periodogram", x = "Frequency", y = "Spectral density") +
          ggplot2::scale_x_continuous()
      }
    )
    
    # apply theme (for plain ggplot outputs only)
    if (inherits(plt, "ggplot")) {
      plt <- add_theme(plt)
    }
    plt
  }, res = 96)
  
  
  
  
  
  
  
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
  
  
  plot_pylike_corr <- function(x, type = c("acf", "pacf"), lag.max = 50, tickSize = 11) {
    type <- match.arg(type)
    x <- as.numeric(na.omit(x))
    n <- length(x)
    stopifnot(n > 3)
    
    # compute acf/pacf without plotting
    obj <- if (type == "acf") stats::acf(x, lag.max = lag.max, plot = FALSE)
    else               stats::pacf(x, lag.max = lag.max, plot = FALSE)
    
    vals <- as.numeric(obj$acf)
    lags <- as.numeric(obj$lag)
    
    # statsmodels-style: usually start at lag 1 for ACF
    if (type == "acf") {
      keep <- lags > 0
      lags <- lags[keep]
      vals <- vals[keep]
    }
    
    df <- data.frame(lag = lags, val = vals)
    
    # common CI band (very similar visually to Python defaults)
    ci <- 1.96 / sqrt(n)
    
    ggplot(df, aes(x = lag, y = val)) +
      geom_hline(yintercept = 0, linewidth = 0.4) +
      geom_ribbon(aes(ymin = -ci, ymax = ci), alpha = 0.2) +
      geom_segment(aes(xend = lag, yend = 0), linewidth = 0.8) +
      geom_point(size = 2) +
      coord_cartesian(ylim = c(-1, 1)) +
      labs(
        title = if (type == "acf") "Autocorrelation" else "Partial Autocorrelation",
        x = "Lag",
        y = toupper(type)
      ) +
      theme_classic(base_size = tickSize) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(fill = NA, linewidth = 0.4)
      )
  }
  
  output$difference2ACFPACF <- renderPlot({
    req(myData_Choice())
    tick <- as.numeric(input$tickSize %||% 11)
    
    p1 <- plot_pylike_corr(myData_Choice(), "acf",  lag.max = 50, tickSize = tick) +
      labs(title = "Autocorrelation of Sales")
    
    p2 <- plot_pylike_corr(myData_Choice(), "pacf", lag.max = 50, tickSize = tick) +
      labs(title = "Partial Autocorrelation of Sales")
    
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
  
  
  output$roadmap_Detailed_Ang_ui <- renderUI({
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
  
  
  
  
  output$roadmap_Detailed_Fr_ui <- renderUI({
    tags$div(
      style = "background:#f7f7f7;padding:14px;border-radius:8px;",
      
      tags$hr(),
      
      tags$p(
        "Ci-dessous, une ",
        tags$b("feuille de route pratique de modélisation SARIMA"),
        "."
      ),
      
      tags$hr(),
      
      tags$h4("[0] - Préparer le terrain : définir le problème de modélisation"),
      
      tags$h5("Ce que les étudiants font"),
      tags$ul(
        tags$li("Définir la ", tags$b("série réponse"), " (y_t) (ce que l’on prévoit)."),
        tags$li("Définir l’", tags$b("indice temporel"), " (quotidien/hebdomadaire/mensuel) et vérifier qu’il est cohérent."),
        tags$li(
          "Définir la ", tags$b("tâche de prévision"), " :",
          tags$ul(
            tags$li("horizon (ex. : 12 mois à l’avance),"),
            tags$li("schéma d’évaluation (origine glissante / rolling-origin ou simple découpage apprentissage/test),"),
            tags$li("métrique de perte (MAE/RMSE/MAPE/sMAPE).")
          )
        ),
        tags$li(
          "Décider si l’on modélise en :",
          tags$ul(
            tags$li(tags$b("niveaux"), " (données brutes),"),
            tags$li(tags$b("log-niveaux"), " (fréquent si la variance augmente avec le niveau),"),
            tags$li("espace transformé ", tags$b("Box–Cox"), " (plus général).")
          )
        )
      ),
      
      tags$h5("Ce qu’ils écrivent (papier)"),
      tags$p(
        tags$b("Méthodes (Données & Objectif). "),
        "« Nous avons modélisé la série temporelle univariée (y_t) observée à une fréquence [mensuelle] de [début] à [fin] (n=...). ",
        "L’objectif était de prévoir à un horizon de (h=...) pas. La performance du modèle a été évaluée à l’aide de [métrique(s)] selon un protocole ",
        "[apprentissage/test ou origine glissante]. »"
      ),
      
      tags$h5("Pièges"),
      tags$ul(
        tags$li("SARIMA suppose un ", tags$b("espacement régulier"), " ; des timestamps irréguliers doivent être corrigés avant toute chose."),
        tags$li("SARIMA modélise ", tags$b("une seule série"), " (sans prédicteurs). Avec des variables explicatives, on parle plutôt de SARIMAX.")
      ),
      
      tags$hr(),
      
      tags$h4("[1] - Décrire les données : taille d’échantillon, valeurs manquantes, statistiques descriptives"),
      
      tags$h5("Ce que les étudiants font"),
      tags$ol(
        tags$li(
          "Rapporter :",
          tags$ul(
            tags$li("taille d’échantillon (n),"),
            tags$li("dates de début/fin,"),
            tags$li("fréquence,"),
            tags$li("nombre/pourcentage de valeurs manquantes.")
          )
        ),
        tags$li(
          "Gérer le manque :",
          tags$ul(
            tags$li("S’il est rare et aléatoire : imputer (interpolation linéaire, interpolation saisonnière)."),
            tags$li("S’il est important : reconsidérer la série, la fréquence ou la source de données.")
          )
        ),
        tags$li(
          "Statistiques descriptives :",
          tags$ul(
            tags$li("moyenne, médiane, écart-type, min/max,"),
            tags$li("éventuellement asymétrie (skewness) / kurtosis,"),
            tags$li("et résumés saisonniers (ex. : moyenne par mois).")
          )
        )
      ),
      
      tags$h5("Ce qu’ils écrivent (style APA)"),
      tags$p(
        tags$b("Résultats (Description des données). "),
        "« La série contient (n=...) observations couvrant [dates] à une fréquence [fréquence]. ",
        "Les valeurs manquantes représentaient (...%) des observations (k=... points). Les observations manquantes ont été traitées par ",
        "[méthode], choisie car [raison]. La distribution de (y_t) présentait une moyenne de (...) (ET=(...)), une médiane (...), ",
        "et un intervalle ([...,...]). »"
      ),
      
      tags$h5("Pièges"),
      tags$ul(
        tags$li("Ne pas imputer « silencieusement » — ", tags$b("toujours justifier"), "."),
        tags$li("Si une transformation logarithmique est appliquée, décrire aussi les statistiques de la série transformée.")
      ),
      
      tags$hr(),
      
      tags$h4("[2] - Explorer visuellement : tendance/saisonnalité/valeurs aberrantes ; rapporter les observations"),
      
      tags$h5("Ce que les étudiants font"),
      tags$p("Produire des graphiques et les commenter :"),
      tags$ul(
        tags$li(tags$b("Courbe"), " de (y_t)."),
        tags$li(tags$b("Graphique saisonnier"), " (ex. : lignes par mois de l’année)."),
        tags$li(tags$b("Boîte à moustaches par saison"), " (mois/trimestre/semaine)."),
        tags$li(
          tags$b("Détection d’outliers"), " :",
          tags$ul(
            tags$li("z-scores, règle IQR, ou méthodes robustes,"),
            tags$li("mais aussi le contexte (fêtes, changements de politique, erreurs de mesure).")
          )
        )
      ),
      
      tags$h5("Ce qu’ils écrivent"),
      tags$p(
        tags$b("Résultats (Analyse exploratoire). "),
        "« L’inspection visuelle a indiqué une tendance [haussière/baissière] et des fluctuations saisonnières récurrentes de période (s=...). ",
        "La variabilité semblait [constante/augmenter avec le niveau], suggérant [aucune transformation / une transformation logarithmique]. ",
        "Plusieurs valeurs potentiellement aberrantes ont été observées autour de [dates], probablement liées à [contexte], et ont été ",
        "[conservées/ajustées] car [raison]. »"
      ),
      
      tags$h5("Pièges"),
      tags$ul(
        tags$li("Les outliers ne sont pas automatiquement « mauvais » : ils peuvent correspondre à des événements réels que la prévision doit respecter."),
        tags$li("Si la variance augmente avec le niveau, SARIMA se comporte souvent mieux après une transformation ", tags$b("log"), " ou ", tags$b("Box–Cox"), ".")
      ),
      
      tags$hr(),
      
      tags$h4("[3] - Décomposer : justifier additif vs multiplicatif ; utiliser STL si robustesse nécessaire"),
      
      tags$h5("Ce que les étudiants font"),
      tags$p("Réaliser une décomposition pour séparer :"),
      tags$ul(
        tags$li("tendance,"),
        tags$li("saisonnalité,"),
        tags$li("reste (bruit).")
      ),
      
      tags$p(tags$b("Choisir la forme du modèle :")),
      tags$ul(
        tags$li(tags$b("Additive :"), " (y_t = T_t + S_t + e_t). À utiliser lorsque l’amplitude saisonnière est à peu près constante."),
        tags$li(tags$b("Multiplicative :"), " (y_t = T_t × S_t × e_t). À utiliser lorsque l’amplitude saisonnière augmente avec le niveau (souvent résolu par log → additif en espace log).")
      ),
      
      tags$p(tags$b("Utiliser la décomposition STL lorsque :")),
      tags$ul(
        tags$li("la saisonnalité évolue lentement au fil du temps,"),
        tags$li("on souhaite une robustesse aux outliers.")
      ),
      
      tags$h5("Ce qu’ils écrivent"),
      tags$p(
        tags$b("Méthodes (Décomposition). "),
        "« Nous avons évalué une structure additive versus multiplicative en examinant si l’amplitude saisonnière évoluait avec le niveau de la série. ",
        "Comme [la variation saisonnière était approximativement constante / augmentait avec le niveau], nous avons utilisé [un modèle additif / une transformation logarithmique] ",
        "et décomposé la série via [décomposition classique / STL]. STL a été retenue pour sa robustesse aux valeurs aberrantes et sa flexibilité ",
        "dans la modélisation d’une saisonnalité évolutive. »"
      ),
      
      tags$h5("Pièges"),
      tags$ul(
        tags$li("« Saison multiplicative » et « transformation log » sont pratiquement meilleurs amis."),
        tags$li("La décomposition STL est descriptive ; l’estimation SARIMA nécessite toujours des vérifications de stationnarité.")
      ),
      
      tags$hr(),
      
      tags$h4("[4] - Vérifier la stationnarité : ADF/KPSS/PP ; justifier la différenciation (d et D)"),
      tags$p("SARIMA requiert la stationnarité ", tags$b("après différenciation"), "."),
      
      tags$h5("Ce que les étudiants font"),
      tags$ol(
        tags$li("Définir la période saisonnière (s) (ex. : 12 pour des données mensuelles, 7 pour des données quotidiennes avec saisonnalité hebdomadaire)."),
        tags$li(
          "Tester la stationnarité sur :",
          tags$ul(
            tags$li("la série originale,"),
            tags$li("après ", tags$b("différenciation ordinaire"), " ((1-B)^d),"),
            tags$li("après ", tags$b("différenciation saisonnière"), " ((1-B^s)^D),"),
            tags$li("et parfois après les deux.")
          )
        )
      ),
      
      tags$h5("Tests de stationnarité : rôle, H0/Ha, et conclusion"),
      tags$ul(
        tags$li(
          tags$b("Test ADF (Augmented Dickey–Fuller)"),
          tags$ul(
            tags$li(tags$b("Rôle :"), " tester si la série se comporte comme si elle avait une ", tags$b("racine unitaire"),
                    " (tendance stochastique), impliquant la non-stationnarité ; le test estime une régression où des différences retardées sont ajoutées pour gérer l’autocorrélation."),
            tags$li(tags$b("H0 :"), " la série a une racine unitaire (non-stationnaire ; les chocs ont des effets permanents)."),
            tags$li(tags$b("Ha :"), " la série n’a pas de racine unitaire (stationnaire autour d’une moyenne ou autour d’une tendance déterministe, selon la spécification ADF)."),
            tags$li(tags$b("Phrase-type de conclusion :"), " « Le test ADF a donné p = [p-value] ; ainsi, au seuil α = [alpha], nous ",
                    tags$b("[rejetons / ne rejetons pas]"), " H0 (racine unitaire). Cela implique que la série est ",
                    tags$b("[stationnaire / non-stationnaire]"), " selon l’ADF ; nous ",
                    tags$b("[n’avons pas appliqué de différenciation supplémentaire / avons appliqué]"), " une différenciation ordinaire [d=…] et/ou saisonnière [D=…] afin d’obtenir une série approximativement stationnaire adaptée à l’estimation SARIMA. »")
          )
        ),
        tags$li(
          tags$b("Test KPSS (Kwiatkowski–Phillips–Schmidt–Shin)"),
          tags$ul(
            tags$li(tags$b("Rôle :"), " tester la stationnarité en examinant si la somme cumulée des résidus (d’une régression de niveau ou de tendance) est trop importante ; c’est un complément à l’ADF en inversant l’hypothèse nulle."),
            tags$li(tags$b("H0 :"), " la série est stationnaire (stationnaire en niveau, ou stationnaire en tendance si une tendance est incluse)."),
            tags$li(tags$b("Ha :"), " la série est non-stationnaire (contient une racine unitaire ou viole la stationnarité)."),
            tags$li(tags$b("Phrase-type de conclusion :"), " « Le test KPSS a donné p = [p-value] ; au seuil α = [alpha], nous ",
                    tags$b("[rejetons / ne rejetons pas]"), " H0 de stationnarité. Cela indique que la série est ",
                    tags$b("[non stationnaire / compatible avec la stationnarité]"), " au sens KPSS ; cela ",
                    tags$b("[soutient l’application / ne nécessite pas]"), " d’une différenciation additionnelle. Nous avons donc retenu [d=…] et [D=…] puis revérifié la stationnarité sur la série transformée. »")
          )
        ),
        tags$li(
          tags$b("Test PP (Phillips–Perron)"),
          tags$ul(
            tags$li(tags$b("Rôle :"), " tester une racine unitaire comme l’ADF, mais en utilisant une correction non paramétrique de l’autocorrélation et de l’hétéroscédasticité (au lieu d’ajouter de nombreux retards)."),
            tags$li(tags$b("H0 :"), " la série a une racine unitaire (non-stationnaire)."),
            tags$li(tags$b("Ha :"), " la série n’a pas de racine unitaire (stationnaire)."),
            tags$li(tags$b("Phrase-type de conclusion :"), " « Le test PP a donné p = [p-value] ; ainsi, au seuil α = [alpha], nous ",
                    tags$b("[rejetons / ne rejetons pas]"), " H0 de racine unitaire. Interprété avec l’ADF et le KPSS, cela suggère que la série est ",
                    tags$b("[stationnaire / non-stationnaire]"), " après application de [d=…] différenciations ordinaires et [D=…] différenciations saisonnières ; cela soutient l’usage d’un SARIMA sur la série différenciée. »")
          )
        )
      ),
      
      tags$p(tags$b("Interpréter ADF/KPSS/PP ensemble (logique à écrire).")),
      tags$ul(
        tags$li(tags$b("Accord idéal :"), " ADF/PP rejettent la racine unitaire (p petit) et KPSS ne rejette pas la stationnarité (p grand) → forte évidence de stationnarité."),
        tags$li(tags$b("Non-stationnarité claire :"), " ADF/PP ne rejettent pas la racine unitaire (p grand) et KPSS rejette la stationnarité (p petit) → forte évidence qu’une différenciation est nécessaire."),
        tags$li(tags$b("Conflits :"), " lorsque les tests divergent, s’appuyer sur l’ensemble des preuves : graphiques + comportement de l’ACF + résultats après différenciation, et indiquer que la conclusion repose sur la convergence des indices plutôt que sur une seule p-value.")
      ),
      
      tags$p(tags$b("Logique de différenciation :")),
      tags$ul(
        tags$li("Choisir ", tags$b("d"), " pour éliminer la tendance / racine unitaire."),
        tags$li("Choisir ", tags$b("D"), " pour éliminer la racine unitaire saisonnière."),
        tags$li("S’arrêter dès que la stationnarité est raisonnable ; éviter la sur-différenciation.")
      ),
      
      tags$h5("Ce qu’ils écrivent"),
      tags$p(
        tags$b("Méthodes (Stationnarité et différenciation). "),
        "« La stationnarité a été évaluée à l’aide des tests ADF, KPSS et PP afin de trianguler l’évidence, ces tests ayant des hypothèses nulles différentes. ",
        "Sur la base des résultats combinés et des diagnostics visuels, nous avons retenu [d=...] différenciations ordinaires et [D=...] différenciations saisonnières avec une période saisonnière (s=...). ",
        "Ce choix visait à supprimer [tendance/racine unitaire saisonnière] tout en évitant la sur-différenciation ; la stationnarité a ensuite été réévaluée sur la série transformée avant d’ajuster les modèles SARIMA. »"
      ),
      
      tags$h5("Pièges (classiques)"),
      tags$ul(
        tags$li(
          tags$b("Sur-différenciation"),
          " provoque :",
          tags$ul(
            tags$li("forte autocorrélation négative au retard 1,"),
            tags$li("variance gonflée,"),
            tags$li("prévisions plus instables.")
          )
        ),
        tags$li("En pratique, D vaut souvent ", tags$b("0 ou 1"), ". Si (D=2) est nécessaire, la série est atypique ou la période saisonnière est mal spécifiée.")
      ),
      
      tags$hr(),
      
      tags$h4("[5] - Ajuster un modèle de référence : Auto-ARIMA pour obtenir un bon point de départ SARIMA"),
      
      tags$h5("Ce que les étudiants font"),
      tags$ul(
        tags$li("Utiliser une procédure ", tags$b("auto-ARIMA"), " (AICc ou similaire) pour proposer : ((p,d,q)(P,D,Q)_s)."),
        tags$li(
          "Documenter :",
          tags$ul(
            tags$li("les transformations utilisées,"),
            tags$li("les contraintes (max p/q, etc.),"),
            tags$li("si une recherche stepwise a été utilisée.")
          )
        ),
        tags$li(tags$b("Important :"), " Auto-ARIMA donne une base solide, pas une vérité gravée dans le marbre.")
      ),
      
      tags$h5("Ce qu’ils écrivent"),
      tags$p(
        tags$b("Méthodes (Modèle de référence). "),
        "« Un modèle SARIMA de référence a été sélectionné via une procédure automatisée basée sur un critère d’information (minimisation de l’AICc) parmi des ordres candidats ((p,q,P,Q)) ",
        "sous contraintes [bornes]. La spécification retenue était SARIMA((p,d,q)(P,D,Q)_s), utilisée comme modèle de référence pour des ajustements ultérieurs guidés par la théorie. »"
      ),
      
      tags$h5("Pièges"),
      tags$ul(
        tags$li("Auto-ARIMA peut sélectionner des modèles statistiquement corrects mais ", tags$b("difficiles à interpréter"), " ou légèrement instables."),
        tags$li("Si l’objectif est la prévision, il est acceptable de préférer des modèles ", tags$b("plus simples"), " à précision comparable.")
      ),
      
      tags$hr(),
      
      tags$h4("[6] - Ajuster un modèle guidé par la théorie : SARIMA manuel via ACF/PACF + tests"),
      
      tags$h5("Ce que les étudiants font"),
      tags$p("À partir de la série différenciée (après choix de (d, D)) :"),
      tags$ol(
        tags$li("Tracer ", tags$b("ACF/PACF"), "."),
        tags$li(
          "Proposer des structures candidates :",
          tags$ul(
            tags$li(
              "Non saisonnier :",
              tags$ul(
                tags$li("AR(p) : la PACF se coupe autour de p ; l’ACF décroît."),
                tags$li("MA(q) : l’ACF se coupe autour de q ; la PACF décroît.")
              )
            ),
            tags$li(
              "Saisonnier :",
              tags$ul(
                tags$li("AR saisonnier (P) : pics PACF aux retards (s, 2s, ...)."),
                tags$li("MA saisonnier (Q) : pics ACF aux retards (s, 2s, ...).")
              )
            )
          )
        ),
        tags$li("Ajuster un ", tags$b("petit ensemble"), " de candidats plausibles (ex. : 3–8 modèles)."),
        tags$li(
          "Comparer via :",
          tags$ul(
            tags$li("AICc/BIC,"),
            tags$li("significativité des paramètres (avec prudence),"),
            tags$li("stabilité / inversibilité.")
          )
        )
      ),
      
      tags$h5("Ce qu’ils écrivent"),
      tags$p(
        tags$b("Méthodes (Construction guidée par la théorie). "),
        "« Les structures SARIMA candidates ont été proposées d’après le comportement ACF/PACF de la série différenciée. Des pics aux retards [lags] suggéraient des composantes non saisonnières [AR/MA], ",
        "tandis que des autocorrélations aux multiples de (s) indiquaient des termes saisonniers [AR/MA]. Plusieurs modèles ont été ajustés et comparés via [AICc/BIC], la sélection finale tenant compte de la parcimonie et de l’adéquation diagnostique. »"
      ),
      
      tags$h5("Pièges"),
      tags$ul(
        tags$li("Les heuristiques ACF/PACF sont des ", tags$b("guides"), ", pas des commandements."),
        tags$li("Éviter de brute-forcer 200 modèles puis d’appeler ça « théorie ». Préférer un ", tags$b("petit ensemble raisonné"), ".")
      ),
      
      tags$hr(),
      
      tags$h4("[7] - Diagnostiquer & comparer : tests des résidus + précision de prévision ; choisir le modèle final"),
      
      tags$h5("Ce que les étudiants font"),
      tags$p(tags$b("Diagnostics des résidus (indispensable) :")),
      tags$ul(
        tags$li("Courbe des résidus (doit ressembler à du bruit)."),
        tags$li("ACF des résidus (pas de gros pics)."),
        tags$li(tags$b("Test de Ljung–Box"), " pour l’autocorrélation des résidus."),
        tags$li("Vérification de normalité (QQ-plot ; Shapiro-Wilk est trop sensible pour grands n)."),
        tags$li("Vérifier l’hétéroscédasticité (variance changeante).")
      ),
      
      tags$p(tags$b("Évaluation de la prévision (indispensable) :")),
      tags$ul(
        tags$li("Échantillon test ou validation croisée rolling."),
        tags$li("Métriques : MAE/RMSE ; MAPE seulement si la série n’est jamais proche de zéro."),
        tags$li(
          "Comparer :",
          tags$ul(
            tags$li("baseline auto-ARIMA,"),
            tags$li("candidats SARIMA manuels,"),
            tags$li("éventuellement un benchmark simple (naïf saisonnier).")
          )
        )
      ),
      
      tags$p(tags$b("Règle de choix (saine) :")),
      tags$ul(
        tags$li("Doit passer les diagnostics de façon raisonnable."),
        tags$li("Doit battre le benchmark naïf."),
        tags$li("Préférer le modèle le plus simple si la précision est quasi identique.")
      ),
      
      tags$h5("Ce qu’ils écrivent"),
      tags$p(
        tags$b("Résultats (Diagnostics et performance). "),
        "« Les diagnostics des résidus indiquaient un comportement proche du bruit blanc : les autocorrélations résiduelles étaient faibles et le test de Ljung–Box était [non significatif/significatif] au seuil (α=...). ",
        "La performance de prévision sur la fenêtre d’évaluation donnait MAE=(...) et RMSE=(...), surpassant les modèles de référence et de benchmark. Sur la base de l’adéquation diagnostique et de la performance prédictive, le modèle final retenu était SARIMA((p,d,q)(P,D,Q)_s). »"
      ),
      
      tags$h5("Pièges"),
      tags$ul(
        tags$li("Un modèle avec un excellent AIC mais des résidus autocorrélés te ", tags$i("raconte une belle histoire"), " — mais fausse."),
        tags$li("Des résidus non normaux peuvent encore donner de bonnes prévisions ; le problème majeur est l’", tags$b("autocorrélation"), " résiduelle.")
      ),
      
      tags$hr(),
      
      tags$h4("[8] - Rédiger le rapport : paragraphes APA à chaque étape ; assembler Méthodes/Résultats"),
      
      tags$h5("Ce que les étudiants font (checklist d’assemblage)"),
      tags$p(tags$b("Section Méthodes")),
      tags$ul(
        tags$li("Données (source, fréquence, traitement du manque, transformation)."),
        tags$li("Approche exploratoire (graphiques, décomposition)."),
        tags$li("Tests de stationnarité et choix de différenciation."),
        tags$li("Baseline (paramètres auto-ARIMA)."),
        tags$li("Justification de la sélection manuelle (ACF/PACF + candidats)."),
        tags$li("Diagnostics et protocole d’évaluation.")
      ),
      
      tags$p(tags$b("Section Résultats")),
      tags$ul(
        tags$li("Résumé des données + observations visuelles clés."),
        tags$li("Résultats de décomposition (tendance/saisonnalité)."),
        tags$li("Résultats des tests de stationnarité et choix (d, D)."),
        tags$li("Paramètres du modèle final."),
        tags$li("Diagnostics et métriques de précision."),
        tags$li("Graphique de prévision + tableau d’erreurs.")
      ),
      
      tags$h5("Guidance (structure APA)"),
      tags$ul(
        tags$li("Pour chaque sous-section : ", tags$b("Ce qu’on a fait → Pourquoi → Ce qu’on a observé → Conclusion"), "."),
        tags$li("Passé pour Méthodes ; passé orienté résultats pour Résultats."),
        tags$li("Écrire les maths avec parcimonie ; donner la spécification complète du modèle une seule fois, clairement.")
      ),
      
      tags$hr(),
      
      tags$h4("Un « pack livrable » propre que les étudiants doivent rendre"),
      tags$ul(
        tags$li(
          "Un notebook/script qui :",
          tags$ul(
            tags$li("charge les données,"),
            tags$li("traite les valeurs manquantes,"),
            tags$li("réalise l’EDA (graphiques),"),
            tags$li("décomposition,"),
            tags$li("tests de stationnarité,"),
            tags$li("baseline auto-ARIMA,"),
            tags$li("candidats manuels,"),
            tags$li("diagnostics,"),
            tags$li("évaluation,"),
            tags$li("prévision finale.")
          )
        ),
        tags$li(
          "Un court rapport avec :",
          tags$ul(
            tags$li("Méthodes + Résultats alignés aux étapes 1–7,"),
            tags$li("figures : courbe temporelle, décomposition, ACF/PACF, ACF des résidus, graphique de prévision,"),
            tags$li("tableau comparatif des modèles candidats (AICc + métriques).")
          )
        )
      ),
      
      tags$hr()
    )
  })
  
  
  
  
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  

  
  
  
  
  # =========================
  # 8 — GARCH (NEW SERVER LOGIC)
  # Requires: rugarch
  # =========================
  
  # Teaching notes for the tab
  output$garch_notes <- renderUI({
    note_box(list(
      "Use GARCH when the residual variance changes over time (volatility clustering).",
      "Model the mean (ARMA) and the conditional variance (GARCH family).",
      "If you model log-returns, interpret forecasts as return forecasts and sigma as volatility."
    ))
  })
  
  output$garch_what_to_do <- renderUI({
    tags$div(
      style = "background:#eef5ff;padding:10px;border-radius:6px;margin-bottom:10px;",
      tags$b("What to do:"),
      tags$ul(
        tags$li("Choose whether to model levels or log-returns (returns usually preferred)."),
        tags$li("Pick a variance model (sGARCH is the baseline; gjrGARCH/eGARCH capture asymmetry)."),
        tags$li("Check diagnostics: standardized residuals should look like white noise; squared residuals should also be clean."),
        tags$li("If residuals are heavy-tailed, switch Normal → Student-t or skewed t."),
        tags$li("Compare models by AIC/BIC and out-of-sample accuracy (if you have a test set).")
      )
    )
  })
  
  # ---- Helper: safe pkg check (you already have has_pkg() globally)
  need_pkg <- function(pkg) {
    validate(need(has_pkg(pkg), paste0("Please install package '", pkg, "' (install.packages('", pkg, "')) to use this tab.")))
  }
  
  # ---- Helper: accuracy metrics (same spirit as your existing accuracy_table)
  garch_accuracy <- function(actual, forecast) {
    a <- as.numeric(actual)
    f <- as.numeric(forecast)
    e <- a - f
    rmse <- sqrt(mean(e^2, na.rm = TRUE))
    mae  <- mean(abs(e), na.rm = TRUE)
    mape <- mean(abs(e / a), na.rm = TRUE)
    smape <- mean(2 * abs(e) / (abs(a) + abs(f)), na.rm = TRUE)
    
    data.frame(
      Metric = c("RMSE", "MAE", "MAPE", "sMAPE"),
      Value = c(rmse, mae, mape, smape),
      stringsAsFactors = FALSE
    )
  }
  
  # ---- Helper: make the modeling series from prepared()/ts split
  garch_series <- reactive({
    req(prepared())
    df <- prepared()$df
    
    # try to use your existing train/test logic if present
    train_n <- NULL
    test_n  <- 0L
    
    if (exists("ts_train_test", mode = "function")) {
      s <- tryCatch(ts_train_test(), error = function(e) NULL)
      if (!is.null(s) && !is.null(s$train_n)) {
        train_n <- s$train_n
        # preferred modeling target (in your app you use y_trans for modeling)
        y_full <- df$y_trans
        if (isTRUE(input$garch_series_type == "logret")) {
          y_full <- as.numeric(y_full)
          y_full <- y_full[is.finite(y_full)]
          y_mod <- diff(log(y_full))
          if (isTRUE(input$garch_scale_100)) y_mod <- 100 * y_mod
          # returns length reduced by 1
          # align train_n accordingly (roughly)
          train_n_mod <- max(5L, min(length(y_mod), train_n - 1L))
          list(y = y_mod, x = df$x[-1], train_n = train_n_mod, test_n = max(0L, length(y_mod) - train_n_mod))
        } else {
          y_mod <- as.numeric(df$y_trans)
          ok <- is.finite(y_mod)
          y_mod <- y_mod[ok]
          x_mod <- df$x[ok]
          train_n_mod <- max(5L, min(length(y_mod), train_n))
          list(y = y_mod, x = x_mod, train_n = train_n_mod, test_n = max(0L, length(y_mod) - train_n_mod))
        }
      }
    }
    
    # fallback: use all data
    if (isTRUE(input$garch_series_type == "logret")) {
      y_full <- as.numeric(df$y_trans)
      y_full <- y_full[is.finite(y_full)]
      y_mod <- diff(log(y_full))
      if (isTRUE(input$garch_scale_100)) y_mod <- 100 * y_mod
      list(y = y_mod, x = df$x[-1], train_n = length(y_mod), test_n = 0L)
    } else {
      y_mod <- as.numeric(df$y_trans)
      ok <- is.finite(y_mod)
      list(y = y_mod[ok], x = df$x[ok], train_n = sum(ok), test_n = 0L)
    }
  })
  
  # ---- Fit event
  garch_fit <- eventReactive(input$fit_garch, {
    need_pkg("rugarch")
    s <- garch_series()
    y <- s$y
    validate(need(length(y) >= 30, "Need at least ~30 observations for a stable GARCH fit."))
    
    # training window
    y_train <- y[seq_len(s$train_n)]
    
    spec <- rugarch::ugarchspec(
      variance.model = list(
        model = input$garch_vmodel,
        garchOrder = c(as.integer(input$garch_p), as.integer(input$garch_q))
      ),
      mean.model = list(
        armaOrder = c(as.integer(input$garch_ar), as.integer(input$garch_ma)),
        include.mean = isTRUE(input$garch_include_mean)
      ),
      distribution.model = input$garch_dist
    )
    
    fit <- tryCatch(
      rugarch::ugarchfit(spec = spec, data = y_train, solver = "hybrid"),
      error = function(e) {
        validate(paste("GARCH fit failed:", e$message))
        NULL
      }
    )
    
    list(spec = spec, fit = fit, series = s, y_train = y_train)
  })
  
  # ---- Model spec text
  output$garch_model_spec <- renderPrint({
    req(garch_fit())
    gf <- garch_fit()
    s <- gf$series
    
    cat("GARCH model specification\n")
    cat("------------------------------------------------------------\n")
    cat("Series modeled        :", if (input$garch_series_type == "logret") "Log-returns" else "Level", "\n")
    cat("Train N               :", s$train_n, "\n")
    cat("Test N                :", s$test_n, "\n\n")
    
    cat("Mean model (ARMA)\n")
    cat("  ARMA(p,q)           : (", input$garch_ar, ",", input$garch_ma, ")\n", sep = "")
    cat("  Include mean (mu)   :", if (isTRUE(input$garch_include_mean)) "Yes" else "No", "\n\n")
    
    cat("Variance model\n")
    cat("  Variant             :", input$garch_vmodel, "\n")
    cat("  Order (p,q)         : (", input$garch_p, ",", input$garch_q, ")\n\n", sep = "")
    
    cat("Innovations\n")
    cat("  Distribution        :", input$garch_dist, "\n\n")
    
    show_methods <- tryCatch(rugarch::infocriteria(gf$fit), error = function(e) NULL)
    if (!is.null(show_methods)) {
      cat("Information criteria\n")
      print(round(show_methods, 4))
    }
  })
  
  # ---- Coef table
  output$garch_coef_table <- renderTable({
    req(garch_fit())
    gf <- garch_fit()
    
    m <- tryCatch(gf$fit@fit$matcoef, error = function(e) NULL)
    validate(need(!is.null(m), "Could not extract coefficient table."))
    
    out <- as.data.frame(m)
    out$term <- rownames(out)
    rownames(out) <- NULL
    out <- out[, c("term", colnames(m)), drop = FALSE]
    names(out) <- c("term", "Estimate", "Std. Error", "t value", "Pr(>|t|)")
    out
  }, rownames = FALSE)
  
  # ---- Model equation (MathJax)
  output$garch_model_equation <- renderUI({
    req(garch_fit())
    
    p  <- as.integer(input$garch_ar)
    q  <- as.integer(input$garch_ma)
    vp <- as.integer(input$garch_p)
    vq <- as.integer(input$garch_q)
    
    mean_eq <- if (p == 0 && q == 0) {
      if (isTRUE(input$garch_include_mean)) "\\mu_t = \\mu" else "\\mu_t = 0"
    } else {
      paste0(
        "y_t = \\mu + \\sum_{i=1}^{", p, "} \\phi_i y_{t-i} + ",
        "\\sum_{j=1}^{", q, "} \\theta_j \\varepsilon_{t-j} + \\varepsilon_t"
      )
    }
    
    var_eq <- switch(
      input$garch_vmodel,
      "sGARCH"   = paste0("\\sigma_t^2 = \\omega + \\sum_{i=1}^{", vq, "} \\alpha_i \\varepsilon_{t-i}^2 + \\sum_{j=1}^{", vp, "} \\beta_j \\sigma_{t-j}^2"),
      "gjrGARCH" = paste0("\\sigma_t^2 = \\omega + \\sum_{i=1}^{", vq, "} (\\alpha_i \\varepsilon_{t-i}^2 + \\gamma_i \\varepsilon_{t-i}^2 \\mathbb{I}(\\varepsilon_{t-i}<0)) + \\sum_{j=1}^{", vp, "} \\beta_j \\sigma_{t-j}^2"),
      "eGARCH"   = paste0("\\log(\\sigma_t^2) = \\omega + \\sum_{i=1}^{", vq, "} (\\alpha_i \\left|\\frac{\\varepsilon_{t-i}}{\\sigma_{t-i}}\\right| + \\gamma_i \\frac{\\varepsilon_{t-i}}{\\sigma_{t-i}}\\right) + \\sum_{j=1}^{", vp, "} \\beta_j \\log(\\sigma_{t-j}^2)"),
      "\\sigma_t^2 = \\omega + \\sum \\alpha \\varepsilon^2 + \\sum \\beta \\sigma^2"
    )
    
    html <- paste0(
      "<p><b>Mean equation:</b></p>",
      "<div>$$", mean_eq, "$$</div>",
      "<p><b>Variance equation:</b></p>",
      "<div>$$", var_eq, "$$</div>",
      "<p>Where $$\\varepsilon_t = \\sigma_t z_t$$ and $$z_t$$ follows the chosen innovation distribution.</p>"
    )
    
    # force MathJax to typeset the updated UI
    session$onFlushed(function() {
      session$sendCustomMessage("mathjax-typeset", "garch_eq_box")
    }, once = TRUE)
    
    HTML(html)
  })
  
  
  
  
  
  # output$garch_model_equation <- renderUI({
  #   req(garch_fit())
  #   
  #   p <- as.integer(input$garch_ar)
  #   q <- as.integer(input$garch_ma)
  #   vp <- as.integer(input$garch_p)
  #   vq <- as.integer(input$garch_q)
  #   
  #   mean_eq <- if (p == 0 && q == 0) {
  #     if (isTRUE(input$garch_include_mean)) {
  #       "\\mu_t = \\mu"
  #     } else {
  #       "\\mu_t = 0"
  #     }
  #   } else {
  #     # compact ARMA notation
  #     paste0(
  #       "y_t = \\mu + \\sum_{i=1}^{", p, "} \\phi_i y_{t-i} + \\sum_{j=1}^{", q, "} \\theta_j \\varepsilon_{t-j} + \\varepsilon_t"
  #     )
  #   }
  #   
  #   var_eq <- switch(
  #     input$garch_vmodel,
  #     "sGARCH"   = paste0("\\sigma_t^2 = \\omega + \\sum_{i=1}^{", vq, "} \\alpha_i \\varepsilon_{t-i}^2 + \\sum_{j=1}^{", vp, "} \\beta_j \\sigma_{t-j}^2"),
  #     "gjrGARCH" = paste0("\\sigma_t^2 = \\omega + \\sum_{i=1}^{", vq, "} \\left(\\alpha_i \\varepsilon_{t-i}^2 + \\gamma_i \\varepsilon_{t-i}^2 \\mathbb{I}(\\varepsilon_{t-i}<0)\\right) + \\sum_{j=1}^{", vp, "} \\beta_j \\sigma_{t-j}^2"),
  #     "eGARCH"   = paste0("\\log(\\sigma_t^2) = \\omega + \\sum_{i=1}^{", vq, "} \\left(\\alpha_i \\left|\\frac{\\varepsilon_{t-i}}{\\sigma_{t-i}}\\right| + \\gamma_i \\frac{\\varepsilon_{t-i}}{\\sigma_{t-i}}\\right) + \\sum_{j=1}^{", vp, "} \\beta_j \\log(\\sigma_{t-j}^2)"),
  #     paste0("\\sigma_t^2 = \\omega + \\sum \\alpha \\varepsilon^2 + \\sum \\beta \\sigma^2")
  #   )
  #   
  #   tags$div(
  #     tags$p(tags$b("Mean equation:")),
  #     tags$div(HTML(paste0("$$", mean_eq, "$$"))),
  #     tags$p(tags$b("Variance equation:")),
  #     tags$div(HTML(paste0("$$", var_eq, "$$"))),
  #     tags$p(
  #       HTML("Where $$\\varepsilon_t = \\sigma_t z_t$$ and $$z_t$$ follows the chosen innovation distribution.")
  #     )
  #   )
  # })
  
  # ---- Diagnostics plots (standardized residuals etc.)
  output$garch_resid_ts <- renderPlot({
    req(garch_fit())
    gf <- garch_fit()
    z <- tryCatch(rugarch::residuals(gf$fit, standardize = TRUE), error = function(e) NULL)
    validate(need(!is.null(z), "No residuals available."))
    z <- as.numeric(z)
    
    plot(z, type = "l", col = "#2C7FB8", main = "Standardized residuals", xlab = "t", ylab = "z_t")
    abline(h = 0, lty = 2, col = "gray50")
  })
  
  output$garch_resid_acf <- renderPlot({
    req(garch_fit())
    gf <- garch_fit()
    z <- tryCatch(rugarch::residuals(gf$fit, standardize = TRUE), error = function(e) NULL)
    validate(need(!is.null(z), "No residuals available."))
    forecast::ggAcf(as.numeric(z)) + ggplot2::labs(title = "ACF of standardized residuals")
  })
  
  output$garch_resid_hist <- renderPlot({
    req(garch_fit())
    gf <- garch_fit()
    z <- tryCatch(rugarch::residuals(gf$fit, standardize = TRUE), error = function(e) NULL)
    validate(need(!is.null(z), "No residuals available."))
    z <- as.numeric(z)
    
    ggplot2::ggplot(data.frame(z = z), ggplot2::aes(x = z)) +
      ggplot2::geom_histogram(bins = 30, fill = "#74a9cf", color = "white") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Histogram (standardized residuals)", x = "z_t", y = "Count")
  })
  
  output$garch_resid_qq <- renderPlot({
    req(garch_fit())
    gf <- garch_fit()
    z <- tryCatch(rugarch::residuals(gf$fit, standardize = TRUE), error = function(e) NULL)
    validate(need(!is.null(z), "No residuals available."))
    z <- as.numeric(z)
    
    qqnorm(z, main = "QQ plot (standardized residuals)", col = "#2C7FB8")
    qqline(z, col = "gray40", lwd = 2)
  })
  
  # Conditional sigma plot
  output$garch_sigma_plot <- renderPlot({
    req(garch_fit())
    gf <- garch_fit()
    sig <- tryCatch(rugarch::sigma(gf$fit), error = function(e) NULL)
    validate(need(!is.null(sig), "No sigma available."))
    sig <- as.numeric(sig)
    
    plot(sig, type = "l", col = "#d95f0e", main = "Conditional volatility (sigma_t)", xlab = "t", ylab = expression(sigma[t]))
  })
  
  # ---- Residual tests (text)
  output$garch_diag_tests <- renderPrint({
    req(garch_fit())
    gf <- garch_fit()
    
    z <- as.numeric(rugarch::residuals(gf$fit, standardize = TRUE))
    z2 <- z^2
    
    L <- as.integer(input$diag_lag %||% 12)
    L <- max(1L, min(L, floor(length(z) / 2)))
    
    cat("GARCH residual diagnostics (standardized residuals)\n")
    cat("------------------------------------------------------------\n")
    
    # LB on z
    lb1 <- tryCatch(Box.test(z, lag = L, type = "Ljung-Box"), error = function(e) NULL)
    if (!is.null(lb1)) cat(sprintf("- Ljung-Box on z_t: Q(%d)=%.4f, p=%.4g\n", L, lb1$statistic, lb1$p.value))
    
    # LB on z^2
    lb2 <- tryCatch(Box.test(z2, lag = L, type = "Ljung-Box"), error = function(e) NULL)
    if (!is.null(lb2)) cat(sprintf("- Ljung-Box on z_t^2: Q(%d)=%.4f, p=%.4g\n", L, lb2$statistic, lb2$p.value))
    
    # Jarque-Bera
    jb <- tryCatch(tseries::jarque.bera.test(z), error = function(e) NULL)
    if (!is.null(jb)) cat(sprintf("- Jarque-Bera: JB=%.4f, p=%.4g\n", jb$statistic, jb$p.value))
    
    # ARCH LM (on standardized residuals)
    if (has_pkg("FinTS")) {
      arch <- tryCatch(FinTS::ArchTest(z, lags = L), error = function(e) NULL)
      if (!is.null(arch)) cat(sprintf("- ARCH LM: LM=%.4f, p=%.4g (lags=%d)\n", arch$statistic, arch$p.value, L))
    } else {
      cat("- ARCH LM: install.packages('FinTS') to enable this test.\n")
    }
    
    cat("\nRule of thumb:\n")
    cat("• We want Ljung-Box p-values on z_t and z_t^2 to be > 0.05 (no leftover autocorrelation / ARCH).\n")
    cat("• Heavy tails: JB often rejects; use Student-t / skewed t innovations.\n")
  })
  
  # ---- Forecast & accuracy
  garch_forecast <- reactive({
    req(garch_fit())
    need_pkg("rugarch")
    
    gf <- garch_fit()
    s <- gf$series
    y <- s$y
    
    # horizon: if test exists => align to test length; else input$garch_h
    h <- if (s$test_n > 0) s$test_n else {
      hh <- suppressWarnings(as.integer(input$garch_h))
      if (!is.finite(hh) || is.na(hh) || hh < 1) 12L else hh
    }
    
    fc <- tryCatch(
      rugarch::ugarchforecast(gf$fit, n.ahead = h),
      error = function(e) {
        validate(paste("Forecast failed:", e$message))
        NULL
      }
    )
    
    # extract mean and sigma forecasts
    mu_hat <- as.numeric(rugarch::fitted(fc))
    sig_hat <- as.numeric(rugarch::sigma(fc))
    
    list(h = h, mu = mu_hat, sigma = sig_hat, series = s, y = y)
  })
  
  output$garch_horizon_note <- renderPrint({
    req(garch_forecast())
    fc <- garch_forecast()
    s <- fc$series
    if (s$test_n > 0) {
      cat("Forecast horizon equals the test set length (h = ", fc$h, ").\n", sep = "")
    } else {
      cat("No test set detected; using future horizon h = ", fc$h, ".\n", sep = "")
    }
  })
  
  output$garch_forecast_table <- renderTable({
    req(garch_forecast())
    fc <- garch_forecast()
    out <- data.frame(
      step = seq_len(fc$h),
      mean_forecast = fc$mu,
      sigma_forecast = fc$sigma
    )
    if (!isTRUE(input$garch_forecast_sigma)) out$sigma_forecast <- NULL
    out
  }, rownames = FALSE)
  
  output$garch_accuracy_table <- renderTable({
    req(garch_forecast())
    fc <- garch_forecast()
    s <- fc$series
    
    if (s$test_n <= 0) {
      return(data.frame(Metric = "Note", Value = "No test set available; accuracy not computed.", stringsAsFactors = FALSE))
    }
    
    y_test <- fc$y[(s$train_n + 1):(s$train_n + s$test_n)]
    y_hat  <- fc$mu[seq_len(length(y_test))]
    garch_accuracy(y_test, y_hat)
  }, rownames = FALSE)
  
  output$garch_forecast_plot <- renderPlot({
    req(garch_forecast())
    fc <- garch_forecast()
    s <- fc$series
    
    y <- fc$y
    train_n <- s$train_n
    h <- fc$h
    
    # Build x for plotting
    x <- s$x
    x_train <- x[seq_len(train_n)]
    x_fc <- x[(train_n + 1):min(length(x), train_n + h)]
    if (length(x_fc) < h) {
      # fallback: simple index extension
      x_fc <- (length(x_train) + 1):(length(x_train) + h)
    }
    
    df_train <- data.frame(x = x_train, y = y[seq_len(train_n)], set = "Train")
    df_fc <- data.frame(x = x_fc, mean = fc$mu[seq_len(h)], sigma = fc$sigma[seq_len(h)])
    
    p <- ggplot() +
      geom_line(data = df_train, aes(x = x, y = y), color = "#2C7FB8") +
      geom_line(data = df_fc, aes(x = x, y = mean), color = "#d95f0e", linewidth = 1) +
      theme_minimal() +
      labs(
        title = "GARCH forecast (mean) and optional volatility band",
        x = "Time",
        y = if (input$garch_series_type == "logret") "Return" else "Level"
      )
    
    if (isTRUE(input$garch_forecast_sigma)) {
      p <- p + geom_ribbon(
        data = df_fc,
        aes(x = x, ymin = mean - 1.96 * sigma, ymax = mean + 1.96 * sigma),
        alpha = 0.15, fill = "#d95f0e"
      )
    }
    
    # if test exists, overlay actual test
    if (s$test_n > 0) {
      df_test <- data.frame(
        x = x[(train_n + 1):(train_n + s$test_n)],
        y = y[(train_n + 1):(train_n + s$test_n)]
      )
      p <- p + geom_line(data = df_test, aes(x = x, y = y), color = "gray40")
    }
    
    p
  })
  
  
  # ---- APA paragraph (GARCH)
  output$garch_apa_paragraph <- renderPrint({
    req(garch_fit())
    gf <- garch_fit()
    
    # pull series/test info
    s <- gf$series
    n_train <- s$train_n
    n_test  <- s$test_n
    
    # info criteria
    ic <- tryCatch(rugarch::infocriteria(gf$fit), error = function(e) NULL)
    aic <- if (!is.null(ic) && "Akaike" %in% names(ic)) as.numeric(ic["Akaike"]) else NA_real_
    bic <- if (!is.null(ic) && "Bayes"  %in% names(ic)) as.numeric(ic["Bayes"])  else NA_real_
    
    # coefs
    mat <- tryCatch(gf$fit@fit$matcoef, error = function(e) NULL)
    
    get_coef <- function(name) {
      if (is.null(mat)) return(list(est = NA_real_, p = NA_real_))
      if (!name %in% rownames(mat)) return(list(est = NA_real_, p = NA_real_))
      list(est = as.numeric(mat[name, 1]), p = as.numeric(mat[name, 4]))
    }
    
    # common variance params
    omega <- get_coef("omega")
    alpha1 <- get_coef("alpha1")
    beta1  <- get_coef("beta1")
    gamma1 <- get_coef("gamma1")   # gjrGARCH
    shape  <- get_coef("shape")    # t / ged
    skew   <- get_coef("skew")     # sstd
    
    # residual diagnostics
    z <- tryCatch(as.numeric(rugarch::residuals(gf$fit, standardize = TRUE)), error = function(e) NULL)
    validate(need(!is.null(z), "Could not compute standardized residuals for APA paragraph."))
    
    L <- as.integer(input$diag_lag %||% 12)
    L <- max(1L, min(L, floor(length(z) / 2)))
    
    lb_z  <- tryCatch(Box.test(z,  lag = L, type = "Ljung-Box"), error = function(e) NULL)
    lb_z2 <- tryCatch(Box.test(z^2, lag = L, type = "Ljung-Box"), error = function(e) NULL)
    
    jb <- tryCatch(tseries::jarque.bera.test(z), error = function(e) NULL)
    
    arch <- NULL
    if (has_pkg("FinTS")) {
      arch <- tryCatch(FinTS::ArchTest(z, lags = L), error = function(e) NULL)
    }
    
    # Forecast accuracy (if test exists)
    acc_text <- ""
    if (n_test > 0) {
      fc <- garch_forecast()
      y_test <- fc$y[(n_train + 1):(n_train + n_test)]
      y_hat  <- fc$mu[seq_len(length(y_test))]
      
      e <- y_test - y_hat
      rmse <- sqrt(mean(e^2, na.rm = TRUE))
      mae  <- mean(abs(e), na.rm = TRUE)
      mape <- mean(abs(e / y_test), na.rm = TRUE)
      
      acc_text <- paste0(
        " Out-of-sample forecast accuracy over the test set (n = ", n_test,
        ") was RMSE = ", fmt_num(rmse), ", MAE = ", fmt_num(mae),
        ", and MAPE = ", fmt_pct(mape), "."
      )
    } else {
      acc_text <- " No holdout test set was detected, so out-of-sample accuracy was not computed."
    }
    
    # build spec text
    mean_part <- paste0(
      "An ARMA(", input$garch_ar, ", ", input$garch_ma, ") mean model",
      if (isTRUE(input$garch_include_mean)) " with an intercept" else " without an intercept",
      " was specified."
    )
    
    var_part <- paste0(
      " Conditional variance was modeled using ", input$garch_vmodel,
      " with order (", input$garch_p, ", ", input$garch_q, ")."
    )
    
    dist_part <- paste0(" Innovations were assumed to follow a ", input$garch_dist, " distribution.")
    
    # parameter highlights (only mention if available)
    parm_bits <- c()
    
    if (is.finite(omega$est)) {
      parm_bits <- c(parm_bits, paste0("ω = ", fmt_num(omega$est), " (p ", fmt_p(omega$p), ")"))
    }
    if (is.finite(alpha1$est)) {
      parm_bits <- c(parm_bits, paste0("α₁ = ", fmt_num(alpha1$est), " (p ", fmt_p(alpha1$p), ")"))
    }
    if (is.finite(beta1$est)) {
      parm_bits <- c(parm_bits, paste0("β₁ = ", fmt_num(beta1$est), " (p ", fmt_p(beta1$p), ")"))
    }
    if (is.finite(gamma1$est)) {
      parm_bits <- c(parm_bits, paste0("γ₁ = ", fmt_num(gamma1$est), " (p ", fmt_p(gamma1$p), ")"))
    }
    if (is.finite(shape$est)) {
      parm_bits <- c(parm_bits, paste0("shape = ", fmt_num(shape$est), " (p ", fmt_p(shape$p), ")"))
    }
    if (is.finite(skew$est)) {
      parm_bits <- c(parm_bits, paste0("skew = ", fmt_num(skew$est), " (p ", fmt_p(skew$p), ")"))
    }
    
    parms_text <- if (length(parm_bits) > 0) {
      paste0(" Key parameter estimates included ", paste(parm_bits, collapse = "; "), ".")
    } else {
      ""
    }
    
    # diagnostics sentence
    diag_bits <- c()
    if (!is.null(lb_z))  diag_bits <- c(diag_bits, paste0("Ljung–Box on zₜ: Q(", L, ") = ", fmt_num(lb_z$statistic), ", p ", fmt_p(lb_z$p.value)))
    if (!is.null(lb_z2)) diag_bits <- c(diag_bits, paste0("Ljung–Box on zₜ²: Q(", L, ") = ", fmt_num(lb_z2$statistic), ", p ", fmt_p(lb_z2$p.value)))
    if (!is.null(arch))  diag_bits <- c(diag_bits, paste0("ARCH LM: LM = ", fmt_num(arch$statistic), ", p ", fmt_p(arch$p.value)))
    if (!is.null(jb))    diag_bits <- c(diag_bits, paste0("Jarque–Bera: JB = ", fmt_num(jb$statistic), ", p ", fmt_p(jb$p.value)))
    
    diag_text <- if (length(diag_bits) > 0) {
      paste0(" Residual diagnostics indicated ", paste(diag_bits, collapse = "; "), ".")
    } else {
      " Residual diagnostics were not available for reporting."
    }
    
    # AIC/BIC text
    ic_text <- if (is.finite(aic) || is.finite(bic)) {
      paste0(
        " Model fit was summarized by information criteria (AIC = ",
        if (is.finite(aic)) fmt_num(aic) else "NA",
        ", BIC = ",
        if (is.finite(bic)) fmt_num(bic) else "NA",
        ")."
      )
    } else {
      ""
    }
    
    series_name <- if (input$garch_series_type == "logret") "log-returns" else "the level series"
    series_scale <- if (input$garch_series_type == "logret" && isTRUE(input$garch_scale_100)) " (scaled by 100)" else ""
    
    paragraph <- paste0(
      "A GARCH model was fitted to ", series_name, series_scale, " using the training sample (n = ", n_train, "). ",
      mean_part, var_part, dist_part,
      ic_text,
      parms_text,
      diag_text,
      acc_text
    )
    
    cat(paragraph)
  })
  
  
  output$garch_conclusion <- renderUI({
    req(garch_fit())
    
    # ---- helpers (use yours if present)
    fmt_num_local <- function(x, d = 3) {
      if (!is.finite(x) || is.na(x)) return("NA")
      formatC(x, format = "f", digits = d)
    }
    fmt_p_local <- function(p) {
      if (!is.finite(p) || is.na(p)) return("= NA")
      if (p < 0.001) return("< .001")
      paste0("= ", sub("^0", "", fmt_num_local(p, 3)))
    }
    fmt_pct_local <- function(x) {
      if (!is.finite(x) || is.na(x)) return("NA")
      paste0(fmt_num_local(100 * x, 2), "%")
    }
    
    gf <- garch_fit()
    s <- gf$series
    n_train <- s$train_n
    n_test  <- s$test_n
    
    # ---- equations
    p  <- as.integer(input$garch_ar)
    q  <- as.integer(input$garch_ma)
    vp <- as.integer(input$garch_p)
    vq <- as.integer(input$garch_q)
    
    mean_eq <- if (p == 0 && q == 0) {
      if (isTRUE(input$garch_include_mean)) "y_t = \\mu + \\varepsilon_t" else "y_t = \\varepsilon_t"
    } else {
      paste0(
        "y_t = \\mu + \\sum_{i=1}^{", p, "} \\phi_i y_{t-i} + ",
        "\\sum_{j=1}^{", q, "} \\theta_j \\varepsilon_{t-j} + \\varepsilon_t"
      )
    }
    
    var_eq <- switch(
      input$garch_vmodel,
      "sGARCH"   = paste0("\\sigma_t^2 = \\omega + \\sum_{i=1}^{", vq, "} \\alpha_i \\varepsilon_{t-i}^2 + \\sum_{j=1}^{", vp, "} \\beta_j \\sigma_{t-j}^2"),
      "gjrGARCH" = paste0("\\sigma_t^2 = \\omega + \\sum_{i=1}^{", vq, "} (\\alpha_i \\varepsilon_{t-i}^2 + \\gamma_i \\varepsilon_{t-i}^2\\,\\mathbb{I}(\\varepsilon_{t-i}<0)) + \\sum_{j=1}^{", vp, "} \\beta_j \\sigma_{t-j}^2"),
      "eGARCH"   = paste0("\\log(\\sigma_t^2) = \\omega + \\sum_{i=1}^{", vq, "} (\\alpha_i |z_{t-i}| + \\gamma_i z_{t-i}) + \\sum_{j=1}^{", vp, "} \\beta_j \\log(\\sigma_{t-j}^2)"),
      "\\sigma_t^2 = \\omega + \\sum \\alpha \\varepsilon^2 + \\sum \\beta \\sigma^2"
    )
    
    # ---- info criteria
    ic <- tryCatch(rugarch::infocriteria(gf$fit), error = function(e) NULL)
    aic <- if (!is.null(ic) && "Akaike" %in% names(ic)) as.numeric(ic["Akaike"]) else NA_real_
    bic <- if (!is.null(ic) && "Bayes"  %in% names(ic)) as.numeric(ic["Bayes"])  else NA_real_
    
    # ---- coefficients and p-values
    mat <- tryCatch(gf$fit@fit$matcoef, error = function(e) NULL)
    
    get_coef <- function(name) {
      if (is.null(mat) || !name %in% rownames(mat)) return(list(est = NA_real_, p = NA_real_))
      list(est = as.numeric(mat[name, 1]), p = as.numeric(mat[name, 4]))
    }
    
    omega <- get_coef("omega")
    alpha1 <- get_coef("alpha1")
    beta1  <- get_coef("beta1")
    gamma1 <- get_coef("gamma1")
    shape  <- get_coef("shape")
    skew   <- get_coef("skew")
    mu     <- get_coef("mu")
    
    # ---- persistence + half-life (only meaningful for sGARCH/gjrGARCH)
    persistence <- NA_real_
    halflife <- NA_real_
    if (input$garch_vmodel %in% c("sGARCH", "gjrGARCH") && is.finite(alpha1$est) && is.finite(beta1$est)) {
      # common approximation: alpha + beta (gjr has extra terms; keep conservative summary)
      persistence <- alpha1$est + beta1$est
      if (is.finite(persistence) && persistence > 0 && persistence < 1) {
        halflife <- log(0.5) / log(persistence)
      }
    }
    
    # ---- residual diagnostics
    z <- tryCatch(as.numeric(rugarch::residuals(gf$fit, standardize = TRUE)), error = function(e) NULL)
    validate(need(!is.null(z), "Could not compute standardized residuals."))
    
    L <- as.integer(input$diag_lag %||% 12)
    L <- max(1L, min(L, floor(length(z) / 2)))
    
    lb_z  <- tryCatch(Box.test(z, lag = L, type = "Ljung-Box"), error = function(e) NULL)
    lb_z2 <- tryCatch(Box.test(z^2, lag = L, type = "Ljung-Box"), error = function(e) NULL)
    jb    <- tryCatch(tseries::jarque.bera.test(z), error = function(e) NULL)
    arch  <- if (has_pkg("FinTS")) tryCatch(FinTS::ArchTest(z, lags = L), error = function(e) NULL) else NULL
    
    # ---- forecast & accuracy (if test exists)
    acc_line <- "No holdout test set was detected; therefore, out-of-sample accuracy statistics were not computed."
    fc_tbl <- NULL
    
    if (n_test > 0) {
      fc <- garch_forecast()   # uses your earlier reactive
      y_test <- fc$y[(n_train + 1):(n_train + n_test)]
      y_hat  <- fc$mu[seq_len(length(y_test))]
      
      e <- y_test - y_hat
      rmse <- sqrt(mean(e^2, na.rm = TRUE))
      mae  <- mean(abs(e), na.rm = TRUE)
      mape <- mean(abs(e / y_test), na.rm = TRUE)
      
      acc_line <- paste0(
        "Forecast accuracy on the test set (n = ", n_test, ") was RMSE = ",
        fmt_num_local(rmse), ", MAE = ", fmt_num_local(mae),
        ", and MAPE = ", fmt_pct_local(mape), "."
      )
      
      # small forecast snippet for conclusion
      fc_tbl <- data.frame(
        step = seq_len(min(5, fc$h)),
        mean_forecast = fc$mu[seq_len(min(5, fc$h))],
        sigma_forecast = fc$sigma[seq_len(min(5, fc$h))]
      )
    }
    
    # ---- narrative bits
    series_text <- if (input$garch_series_type == "logret") {
      if (isTRUE(input$garch_scale_100)) "log-returns (scaled by 100)" else "log-returns"
    } else {
      "the level series"
    }
    
    dist_long <- switch(
      input$garch_dist,
      "norm" = "Gaussian",
      "std"  = "Student’s t",
      "sstd" = "skewed Student’s t",
      "ged"  = "generalized error (GED)",
      input$garch_dist
    )
    
    # ---- build HTML (academic format)
    html <- paste0(
      "<h4 style='margin-top:0;'>Conclusion (GARCH modelling)</h4>",
      
      "<p><b>Model overview.</b> A ", input$garch_vmodel, " model was estimated for ", series_text,
      " using a training sample of <i>n</i> = ", n_train, " observations. The conditional mean was specified as ARMA(",
      p, ", ", q, ") ", if (isTRUE(input$garch_include_mean)) "with an intercept" else "without an intercept",
      ", and innovations were assumed to follow a ", dist_long, " distribution.</p>",
      
      "<p><b>Model equations.</b></p>",
      "<div style='margin-left:6px;'>$$", mean_eq, "$$</div>",
      "<div style='margin-left:6px;'>$$", var_eq, "$$</div>",
      "<p>with $$\\varepsilon_t = \\sigma_t z_t$$.</p>",
      
      "<p><b>Key estimation results.</b> ",
      if (is.finite(aic) || is.finite(bic)) paste0("Information criteria indicated AIC = ", fmt_num_local(aic), " and BIC = ", fmt_num_local(bic), ". ") else "",
      if (is.finite(mu$est)) paste0("The estimated mean term was \\(\\mu\\) = ", fmt_num_local(mu$est), " (p ", fmt_p_local(mu$p), "). ") else "",
      if (is.finite(omega$est)) paste0("The variance intercept was \\(\\omega\\) = ", fmt_num_local(omega$est), " (p ", fmt_p_local(omega$p), "). ") else "",
      if (is.finite(alpha1$est)) paste0("The ARCH effect \\(\\alpha_1\\) = ", fmt_num_local(alpha1$est), " (p ", fmt_p_local(alpha1$p), "), ") else "",
      if (is.finite(beta1$est))  paste0("and the GARCH effect \\(\\beta_1\\) = ", fmt_num_local(beta1$est), " (p ", fmt_p_local(beta1$p), "). ") else "",
      if (is.finite(gamma1$est)) paste0("An asymmetric (leverage) component \\(\\gamma_1\\) = ", fmt_num_local(gamma1$est), " (p ", fmt_p_local(gamma1$p), ") was also estimated. ") else "",
      if (is.finite(shape$est))  paste0("Tail thickness (shape) was estimated as ", fmt_num_local(shape$est), " (p ", fmt_p_local(shape$p), "). ") else "",
      if (is.finite(skew$est))   paste0("Skewness (skew) was estimated as ", fmt_num_local(skew$est), " (p ", fmt_p_local(skew$p), "). ") else "",
      if (is.finite(persistence)) paste0("Volatility persistence (approx. \\(\\alpha_1 + \\beta_1\\)) was ", fmt_num_local(persistence), 
                                         if (is.finite(halflife)) paste0(", corresponding to an approximate half-life of ", fmt_num_local(halflife, 2), " periods. ") else ". ")
      else "",
      "</p>",
      
      "<p><b>Residual diagnostics.</b> ",
      if (!is.null(lb_z))  paste0("Ljung–Box tests suggested ", if (lb_z$p.value > 0.05) "no" else "remaining", " autocorrelation in standardized residuals (Q(", L, ") = ", fmt_num_local(lb_z$statistic), ", p ", fmt_p_local(lb_z$p.value), "). ") else "",
      if (!is.null(lb_z2)) paste0("For squared residuals, the Ljung–Box test ", if (lb_z2$p.value > 0.05) "did not indicate" else "indicated", " remaining ARCH structure (Q(", L, ") = ", fmt_num_local(lb_z2$statistic), ", p ", fmt_p_local(lb_z2$p.value), "). ") else "",
      if (!is.null(arch))  paste0("The ARCH LM test ", if (arch$p.value > 0.05) "did not provide evidence" else "provided evidence", " of additional conditional heteroskedasticity (LM = ", fmt_num_local(arch$statistic), ", p ", fmt_p_local(arch$p.value), "). ") else "",
      if (!is.null(jb))    paste0("Normality was assessed using Jarque–Bera (JB = ", fmt_num_local(jb$statistic), ", p ", fmt_p_local(jb$p.value), "), which commonly rejects under heavy tails—consistent with adopting non-Gaussian innovations when appropriate. ") else "",
      "</p>",
      
      "<p><b>Forecast performance.</b> ", acc_line, "</p>",
      
      if (!is.null(fc_tbl)) {
        paste0(
          "<p><b>Forecast excerpt (first 5 steps).</b></p>",
          "<div style='margin-left:6px;'>",
          paste0(
            "<table class='table table-condensed' style='width:100%;max-width:520px;'>",
            "<thead><tr><th>Step</th><th>Mean forecast</th><th>Sigma forecast</th></tr></thead><tbody>",
            paste(
              apply(fc_tbl, 1, function(r) {
                paste0(
                  "<tr><td>", r[[1]], "</td><td>", fmt_num_local(as.numeric(r[[2]])),
                  "</td><td>", fmt_num_local(as.numeric(r[[3]])), "</td></tr>"
                )
              }),
              collapse = ""
            ),
            "</tbody></table>"
          ),
          "</div>"
        )
      } else "",
      
      "<p><b>Overall conclusion.</b> Collectively, the estimated parameters and residual diagnostics ",
      "support the use of a conditional heteroskedasticity framework for capturing time-varying volatility in the series. ",
      "When diagnostics indicate limited remaining autocorrelation in \\(z_t\\) and \\(z_t^2\\), the fitted GARCH specification ",
      "may be considered adequate for inference and forecasting. Where residual tests suggest remaining dependence or heavy tails, ",
      "improvements may include refining the ARMA mean orders, increasing the GARCH order, or adopting heavier-tailed/skewed innovation distributions.</p>"
    )
    
    # trigger MathJax typesetting for this conclusion box
    session$onFlushed(function() {
      session$sendCustomMessage("mathjax-typeset", "garch_conclusion_box")
    }, once = TRUE)
    
    HTML(html)
  })
  
  
  
  
  # =========================
  # 9 — Auto-GARCH (SERVER)
  # =========================
  
  output$autogarch_notes <- renderUI({
    if (isTRUE(input$show_teaching_notes)) {
      tags$div(
        class = "alert alert-info",
        tags$b("Auto-GARCH:"),
        " searches across GARCH variants, orders, and innovation distributions, and selects the best model using AIC/BIC.",
        tags$br(),
        "Tip: Start with (1,1) + {sGARCH, gjrGARCH} and {std, norm} to keep search fast."
      )
    }
  })
  
  output$autogarch_what_to_do <- renderUI({
    tags$div(
      tags$ol(
        tags$li("Upload data and set frequency/transform in the left sidebar."),
        tags$li("Choose series type (level vs log-returns)."),
        tags$li("Choose search space (mean ARMA grid, variance models, distributions)."),
        tags$li("Click “Run Auto-GARCH Search”."),
        tags$li("Inspect top candidates, diagnostics, residual tests, then forecast performance.")
      )
    )
  })
  
  # ---- Helper: build series for Auto-GARCH
  autogarch_series <- reactive({
    req(ts_train_test())
    s <- ts_train_test()
    y_full <- as.numeric(s$ts_full)
    
    if (input$autogarch_series_type == "logret") {
      # log-returns: diff(log(y))
      y_full <- diff(log(y_full))
      y_full <- y_full[is.finite(y_full)]
      if (isTRUE(input$autogarch_scale_100)) y_full <- 100 * y_full
    } else {
      y_full <- y_full[is.finite(y_full)]
    }
    
    # train/test sizes based on original split proportion
    # If train_prop == 1 => test size 0
    n_full <- length(y_full)
    train_prop <- suppressWarnings(as.numeric(input$train_prop))
    if (!is.finite(train_prop) || train_prop <= 0) train_prop <- 1
    n_train <- if (train_prop >= 0.999) n_full else max(10L, floor(train_prop * n_full))
    n_test  <- max(0L, n_full - n_train)
    
    list(
      y = y_full,
      train_n = n_train,
      test_n  = n_test
    )
  })
  
  # ---- Auto search (eventReactive on button)
  autogarch_search <- eventReactive(input$fit_autogarch, {
    validate(need(has_pkg("rugarch"), "Package 'rugarch' is required for Auto-GARCH. Please install it."))
    s <- autogarch_series()
    y <- s$y
    validate(need(length(y) >= 80, "Need at least ~80 observations for a stable Auto-GARCH search."))
    
    y_train <- y[seq_len(s$train_n)]
    
    # grids
    vmodels <- input$autogarch_vmodels
    dists   <- input$autogarch_dists
    validate(need(length(vmodels) > 0, "Select at least one GARCH variant."))
    validate(need(length(dists) > 0, "Select at least one distribution."))
    
    gorders <- switch(
      input$autogarch_orders,
      "11"    = list(c(1L, 1L)),
      "small" = list(c(1L, 1L), c(1L, 2L), c(2L, 1L)),
      "22"    = list(c(1L, 1L), c(1L, 2L), c(2L, 1L), c(2L, 2L)),
      list(c(1L, 1L))
    )
    
    if (isTRUE(input$autogarch_search_mean)) {
      pmax <- as.integer(input$autogarch_pmax); if (!is.finite(pmax) || pmax < 0) pmax <- 0L
      qmax <- as.integer(input$autogarch_qmax); if (!is.finite(qmax) || qmax < 0) qmax <- 0L
      arma_grid <- expand.grid(ar = 0:pmax, ma = 0:qmax)
    } else {
      arma_grid <- data.frame(ar = 0L, ma = 0L)
    }
    
    include_mean <- isTRUE(input$autogarch_include_mean)
    
    # safe fit one
    fit_one <- function(vmodel, dist, go, ar, ma) {
      spec <- rugarch::ugarchspec(
        variance.model = list(model = vmodel, garchOrder = c(go[1], go[2])),
        mean.model = list(armaOrder = c(ar, ma), include.mean = include_mean),
        distribution.model = dist
      )
      fit <- tryCatch(
        rugarch::ugarchfit(spec = spec, data = y_train, solver = "hybrid"),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      ic <- tryCatch(rugarch::infocriteria(fit), error = function(e) NULL)
      if (is.null(ic)) return(NULL)
      
      data.frame(
        vmodel = vmodel, dist = dist,
        garch_p = go[1], garch_q = go[2],
        ar = ar, ma = ma,
        AIC = as.numeric(ic["Akaike"]),
        BIC = as.numeric(ic["Bayes"]),
        stringsAsFactors = FALSE,
        fit = I(list(fit))
      )
    }
    
    # run search with progress
    total <- length(vmodels) * length(dists) * length(gorders) * nrow(arma_grid)
    res_list <- vector("list", total)
    k <- 1L
    
    withProgress(message = "Auto-GARCH search", value = 0, {
      for (vm in vmodels) {
        for (di in dists) {
          for (go in gorders) {
            for (i in seq_len(nrow(arma_grid))) {
              incProgress(1 / total, detail = paste(vm, di, paste0("(", go[1], ",", go[2], ")"), "ARMA", arma_grid$ar[i], arma_grid$ma[i]))
              out <- fit_one(vm, di, go, arma_grid$ar[i], arma_grid$ma[i])
              if (!is.null(out)) {
                res_list[[k]] <- out
                k <- k + 1L
              }
            }
          }
        }
      }
    })
    
    res_list <- res_list[seq_len(k - 1L)]
    validate(need(length(res_list) > 0, "All candidate models failed. Try smaller ARMA grid, fewer distributions, or only sGARCH(1,1)."))
    
    res <- do.call(rbind, res_list)
    
    # rank by selection criterion
    ic_col <- input$autogarch_select_ic
    res <- res[order(res[[ic_col]]), , drop = FALSE]
    res
  })
  
  # ---- best fit + series bundle
  autogarch_best <- reactive({
    req(autogarch_search())
    res <- autogarch_search()
    list(
      fit = res$fit[[1]],
      meta = res[1, c("vmodel","dist","garch_p","garch_q","ar","ma","AIC","BIC"), drop = FALSE],
      series = autogarch_series()
    )
  })
  
  # ---- search results table
  output$autogarch_rank_table <- renderTable({
    validate(need(input$fit_autogarch > 0, "Click “Run Auto-GARCH Search” to see results."))
    req(autogarch_search())
    topk <- as.integer(input$autogarch_topk)
    if (!is.finite(topk) || topk < 5) topk <- 10L
    
    res <- autogarch_search()
    out <- head(res, topk)
    out$fit <- NULL
    out
  }, rownames = FALSE)
  
  output$autogarch_best_spec <- renderPrint({
    validate(need(input$fit_autogarch > 0, "Click “Run Auto-GARCH Search” first."))
    req(autogarch_best())
    b <- autogarch_best()
    cat("Best Auto-GARCH model selected.\n\n")
    print(b$meta)
    cat("\nInfo criteria (rugarch::infocriteria):\n")
    print(round(rugarch::infocriteria(b$fit), 4))
  })
  
  # ---- model spec + coefficient table
  output$autogarch_model_spec <- renderPrint({
    validate(need(input$fit_autogarch > 0, "Click “Run Auto-GARCH Search” first."))
    req(autogarch_best())
    b <- autogarch_best()
    
    cat("Auto-GARCH best specification:\n")
    print(b$meta)
    cat("\n\nRUGARCH fit summary:\n")
    show(b$fit)
  })
  
  output$autogarch_coef_table <- renderTable({
    validate(need(input$fit_autogarch > 0, "Click “Run Auto-GARCH Search” first."))
    req(autogarch_best())
    fit <- autogarch_best()$fit
    mat <- tryCatch(fit@fit$matcoef, error = function(e) NULL)
    validate(need(!is.null(mat), "Coefficient table unavailable."))
    
    data.frame(
      Term = rownames(mat),
      Estimate = as.numeric(mat[, 1]),
      `Std. Error` = as.numeric(mat[, 2]),
      `t value` = as.numeric(mat[, 3]),
      `Pr(>|t|)` = as.numeric(mat[, 4]),
      row.names = NULL,
      check.names = FALSE
    )
  }, digits = 5)
  
  # ---- equation (same style as your GARCH equation)
  output$autogarch_model_equation <- renderUI({
    validate(need(input$fit_autogarch > 0, "Run Auto-GARCH Search to generate the equation."))
    req(autogarch_best())
    b <- autogarch_best()
    m <- b$meta[1, ]
    
    p  <- as.integer(m$ar); if (!is.finite(p)) p <- 0L
    q  <- as.integer(m$ma); if (!is.finite(q)) q <- 0L
    vp <- as.integer(m$garch_p); if (!is.finite(vp)) vp <- 1L
    vq <- as.integer(m$garch_q); if (!is.finite(vq)) vq <- 1L
    vmodel <- as.character(m$vmodel)
    
    mean_eq <- if (p == 0 && q == 0) {
      if (isTRUE(input$autogarch_include_mean)) "y_t = \\mu + \\varepsilon_t" else "y_t = \\varepsilon_t"
    } else {
      paste0("y_t = \\mu + \\sum_{i=1}^{", p, "} \\phi_i y_{t-i} + \\sum_{j=1}^{", q, "} \\theta_j \\varepsilon_{t-j} + \\varepsilon_t")
    }
    
    var_eq <- switch(
      vmodel,
      "sGARCH"   = paste0("\\sigma_t^2 = \\omega + \\sum_{i=1}^{", vq, "} \\alpha_i \\varepsilon_{t-i}^2 + \\sum_{j=1}^{", vp, "} \\beta_j \\sigma_{t-j}^2"),
      "gjrGARCH" = paste0("\\sigma_t^2 = \\omega + \\sum_{i=1}^{", vq, "} (\\alpha_i \\varepsilon_{t-i}^2 + \\gamma_i \\varepsilon_{t-i}^2\\,\\mathbb{I}(\\varepsilon_{t-i}<0)) + \\sum_{j=1}^{", vp, "} \\beta_j \\sigma_{t-j}^2"),
      "eGARCH"   = paste0("\\log(\\sigma_t^2) = \\omega + \\sum_{i=1}^{", vq, "} (\\alpha_i |z_{t-i}| + \\gamma_i z_{t-i}) + \\sum_{j=1}^{", vp, "} \\beta_j \\log(\\sigma_{t-j}^2)"),
      paste0("\\sigma_t^2 = \\omega + \\sum \\alpha \\varepsilon^2 + \\sum \\beta \\sigma^2")
    )
    
    html <- paste0(
      "<p><b>Mean equation:</b></p><div>$$", mean_eq, "$$</div>",
      "<p><b>Variance equation:</b></p><div>$$", var_eq, "$$</div>",
      "<p>Where $$\\varepsilon_t = \\sigma_t z_t$$.</p>"
    )
    
    session$onFlushed(function() {
      session$sendCustomMessage("mathjax-typeset", "autogarch_eq_box")
    }, once = TRUE)
    
    HTML(html)
  })
  
  # ---- residuals + sigma
  autogarch_resids <- reactive({
    req(autogarch_best())
    fit <- autogarch_best()$fit
    z <- tryCatch(as.numeric(rugarch::residuals(fit, standardize = TRUE)), error = function(e) NULL)
    validate(need(!is.null(z), "Cannot compute standardized residuals."))
    z
  })
  
  autogarch_sigma <- reactive({
    req(autogarch_best())
    fit <- autogarch_best()$fit
    sig <- tryCatch(as.numeric(rugarch::sigma(fit)), error = function(e) NULL)
    validate(need(!is.null(sig), "Cannot compute sigma."))
    sig
  })
  
  # ---- diagnostics plots
  output$autogarch_resid_ts <- renderPlot({
    validate(need(input$fit_autogarch > 0, "Run Auto-GARCH Search first."))
    z <- autogarch_resids()
    plot(z, type = "l", main = "Standardized residuals", ylab = "z_t", xlab = "t")
    abline(h = 0, col = "red", lty = 2)
  })
  
  output$autogarch_resid_acf <- renderPlot({
    validate(need(input$fit_autogarch > 0, "Run Auto-GARCH Search first."))
    z <- autogarch_resids()
    stats::acf(z, main = "ACF of standardized residuals")
  })
  
  output$autogarch_resid_hist <- renderPlot({
    validate(need(input$fit_autogarch > 0, "Run Auto-GARCH Search first."))
    z <- autogarch_resids()
    hist(z, breaks = 30, main = "Histogram of z_t", xlab = "z_t")
  })
  
  output$autogarch_resid_qq <- renderPlot({
    validate(need(input$fit_autogarch > 0, "Run Auto-GARCH Search first."))
    z <- autogarch_resids()
    qqnorm(z, main = "Q–Q plot of z_t"); qqline(z, col = "red")
  })
  
  output$autogarch_sigma_plot <- renderPlot({
    validate(need(input$fit_autogarch > 0, "Run Auto-GARCH Search first."))
    sig <- autogarch_sigma()
    plot(sig, type = "l", main = "Conditional volatility (sigma)", ylab = "sigma_t", xlab = "t")
  })
  
  # ---- residual tests text
  output$autogarch_diag_tests <- renderPrint({
    validate(need(input$fit_autogarch > 0, "Run Auto-GARCH Search first."))
    z <- autogarch_resids()
    
    L <- as.integer(input$diag_lag)
    if (!is.finite(L) || L < 1) L <- 12L
    L <- max(1L, min(L, floor(length(z) / 2)))
    
    cat("Residual tests (standardized residuals z_t)\n")
    cat("=======================================\n\n")
    
    lb1 <- Box.test(z,  lag = L, type = "Ljung-Box")
    lb2 <- Box.test(z^2, lag = L, type = "Ljung-Box")
    
    cat("Ljung–Box on z_t:\n"); print(lb1); cat("\n")
    cat("Ljung–Box on z_t^2:\n"); print(lb2); cat("\n")
    
    if (has_pkg("FinTS")) {
      cat("ARCH LM test (FinTS::ArchTest):\n")
      print(FinTS::ArchTest(z, lags = L)); cat("\n")
    } else {
      cat("ARCH LM test: FinTS not installed.\n\n")
    }
    
    if (has_pkg("tseries")) {
      cat("Jarque–Bera normality test (tseries):\n")
      print(tseries::jarque.bera.test(z)); cat("\n")
    } else {
      cat("Jarque–Bera: tseries not installed.\n\n")
    }
  })
  
  # ---- forecast + accuracy
  autogarch_forecast <- reactive({
    req(autogarch_best())
    b <- autogarch_best()
    fit <- b$fit
    s <- b$series
    
    n_test <- as.integer(s$test_n); if (!is.finite(n_test)) n_test <- 0L
    n_train <- as.integer(s$train_n)
    
    if (n_test > 0) {
      h <- n_test
    } else {
      h_in <- suppressWarnings(as.integer(input$autogarch_h))
      if (!is.finite(h_in) || h_in < 1) h_in <- 20L
      h <- h_in
    }
    
    fc <- rugarch::ugarchforecast(fit, n.ahead = h)
    
    mu <- tryCatch(as.numeric(rugarch::fitted(fc)), error = function(e) rep(NA_real_, h))
    sig <- tryCatch(as.numeric(rugarch::sigma(fc)), error = function(e) rep(NA_real_, h))
    
    list(h = h, fc = fc, mu = mu, sigma = sig)
  })
  
  output$autogarch_horizon_note <- renderPrint({
    validate(need(input$fit_autogarch > 0, "Run Auto-GARCH Search first."))
    s <- autogarch_series()
    if (s$test_n > 0) {
      cat("Holdout test set detected. Forecast horizon h was set to the test length (h =", s$test_n, ").\n")
    } else {
      h <- autogarch_forecast()$h
      cat("No test set detected. Forecast horizon h was set to", h, "using the Auto-GARCH horizon input.\n")
    }
  })
  
  output$autogarch_accuracy_table <- renderTable({
    validate(need(input$fit_autogarch > 0, "Run Auto-GARCH Search first."))
    b <- autogarch_best()
    s <- b$series
    if (s$test_n <= 0) return(NULL)
    
    y <- b$series$y
    y_test <- y[(s$train_n + 1):(s$train_n + s$test_n)]
    
    mu <- autogarch_forecast()$mu[seq_along(y_test)]
    validate(need(length(mu) == length(y_test), "Forecast length mismatch."))
    
    accuracy_df(y_test, mu)
  }, rownames = FALSE)
  
  output$autogarch_forecast_table <- renderTable({
    validate(need(input$fit_autogarch > 0, "Run Auto-GARCH Search first."))
    f <- autogarch_forecast()
    h <- f$h
    out <- data.frame(
      step = 1:h,
      mean_forecast = f$mu,
      sigma_forecast = if (isTRUE(input$autogarch_forecast_sigma)) f$sigma else NA_real_
    )
    out
  }, digits = 6)
  
  output$autogarch_forecast_plot <- renderPlot({
    validate(need(input$fit_autogarch > 0, "Run Auto-GARCH Search first."))
    b <- autogarch_best()
    s <- b$series
    y <- b$series$y
    
    f <- autogarch_forecast()
    h <- f$h
    
    # plot last window + forecast
    n <- length(y)
    window <- min(200, n)
    y_tail <- y[(n - window + 1):n]
    
    plot(y_tail, type = "l", main = "Auto-GARCH: mean forecast", ylab = "Series", xlab = "t (tail)")
    lines((window + 1):(window + h), f$mu, col = "blue", lwd = 2)
    
    if (isTRUE(input$autogarch_forecast_sigma)) {
      # simple +/- 2*sigma band for visualization
      up <- f$mu + 2 * f$sigma
      lo <- f$mu - 2 * f$sigma
      lines((window + 1):(window + h), up, col = "gray40", lty = 2)
      lines((window + 1):(window + h), lo, col = "gray40", lty = 2)
      legend("topleft", legend = c("history", "mean forecast", "±2 sigma"), col = c("black","blue","gray40"),
             lty = c(1,1,2), bty = "n")
    } else {
      legend("topleft", legend = c("history", "mean forecast"), col = c("black","blue"),
             lty = c(1,1), bty = "n")
    }
  })
  
  # ---- APA paragraph (Auto-GARCH)
  output$autogarch_apa_paragraph <- renderPrint({
    validate(need(input$fit_autogarch > 0, "Run Auto-GARCH Search first."))
    req(autogarch_best())
    b <- autogarch_best()
    m <- b$meta[1, ]
    
    ic <- rugarch::infocriteria(b$fit)
    aic <- as.numeric(ic["Akaike"])
    bic <- as.numeric(ic["Bayes"])
    
    z <- autogarch_resids()
    L <- as.integer(input$diag_lag); if (!is.finite(L) || L < 1) L <- 12L
    L <- max(1L, min(L, floor(length(z) / 2)))
    
    lbz  <- Box.test(z, lag = L, type = "Ljung-Box")
    lbz2 <- Box.test(z^2, lag = L, type = "Ljung-Box")
    
    cat(
      "An Auto-GARCH search selected a ", m$vmodel, " model with GARCH order (",
      m$garch_p, ", ", m$garch_q, "), ARMA(", m$ar, ", ", m$ma, ") mean, and ",
      m$dist, " innovations. Model fit was summarized by AIC = ", fmt_num(aic, 3),
      " and BIC = ", fmt_num(bic, 3), ". Residual diagnostics indicated Ljung–Box ",
      "tests on standardized residuals z_t (Q(", L, ") = ", fmt_num(lbz$statistic, 3),
      ", p ", fmt_p(lbz$p.value), ") and on squared residuals z_t^2 (Q(", L, ") = ",
      fmt_num(lbz2$statistic, 3), ", p ", fmt_p(lbz2$p.value), ").",
      sep = ""
    )
  })
  
  # ---- Conclusion (academic) (Auto-GARCH)
  output$autogarch_conclusion <- renderUI({
    validate(need(input$fit_autogarch > 0, "Run Auto-GARCH Search first."))
    req(autogarch_best())
    
    b <- autogarch_best()
    m <- b$meta[1, ]
    s <- b$series
    n_train <- s$train_n
    n_test  <- s$test_n
    
    # equation rendered elsewhere; conclusion references it + tests + accuracy
    z <- autogarch_resids()
    L <- as.integer(input$diag_lag); if (!is.finite(L) || L < 1) L <- 12L
    L <- max(1L, min(L, floor(length(z) / 2)))
    
    lbz  <- Box.test(z, lag = L, type = "Ljung-Box")
    lbz2 <- Box.test(z^2, lag = L, type = "Ljung-Box")
    
    ic <- rugarch::infocriteria(b$fit)
    aic <- as.numeric(ic["Akaike"])
    bic <- as.numeric(ic["Bayes"])
    
    # accuracy summary (if test exists)
    acc_txt <- "No holdout test set was detected; therefore, out-of-sample accuracy was not computed."
    if (n_test > 0) {
      acc <- tryCatch(output$autogarch_accuracy_table(), error = function(e) NULL)
      # better: recompute quickly
      y <- b$series$y
      y_test <- y[(n_train + 1):(n_train + n_test)]
      mu <- autogarch_forecast()$mu[seq_along(y_test)]
      e <- y_test - mu
      rmse <- sqrt(mean(e^2, na.rm = TRUE))
      mae  <- mean(abs(e), na.rm = TRUE)
      acc_txt <- paste0("Forecast accuracy on the test set (n = ", n_test, ") was RMSE = ", fmt_num(rmse, 3),
                        " and MAE = ", fmt_num(mae, 3), ".")
    }
    
    html <- paste0(
      "<h4 style='margin-top:0;'>Auto-GARCH conclusion (academic)</h4>",
      "<p><b>Selected specification.</b> The Auto-GARCH search selected <b>", m$vmodel,
      "</b>(", m$garch_p, ",", m$garch_q, ") with ARMA(", m$ar, ",", m$ma,
      ") mean and <b>", m$dist, "</b> innovations. Fit indices were AIC = ",
      fmt_num(aic, 3), " and BIC = ", fmt_num(bic, 3), ".</p>",
      "<p><b>Model adequacy.</b> Ljung–Box tests suggested ",
      "Q(", L, ") = ", fmt_num(lbz$statistic, 3), " (p ", fmt_p(lbz$p.value),
      ") for z<sub>t</sub> and Q(", L, ") = ", fmt_num(lbz2$statistic, 3),
      " (p ", fmt_p(lbz2$p.value), ") for z<sub>t</sub><sup>2</sup>. ",
      "When these are non-significant, the fitted variance dynamics are typically considered adequate.</p>",
      "<p><b>Forecasting.</b> ", acc_txt, "</p>",
      "<p><b>Overall.</b> The automated selection provides a defensible volatility model candidate. ",
      "If residual tests indicate remaining dependence, refine the mean ARMA orders, broaden the order search, or consider alternative innovation distributions.</p>"
    )
    
    session$onFlushed(function() {
      session$sendCustomMessage("mathjax-typeset", "autogarch_conclusion_box")
    }, once = TRUE)
    
    HTML(html)
  })
  
  
  
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  
  
  
  # # Manual SARIMAX selector UI
  # output$sarimax_xreg_ui <- renderUI({
  #   req(prepared())
  #   df <- prepared()$df
  #   
  #   # exclude date/value columns if you have them; adjust names as needed
  #   candidates <- setdiff(names(df), c(prepared()$date_col, prepared()$value_col))
  #   candidates <- candidates[candidates != ""]
  #   selectizeInput(
  #     "sarimax_x_cols",
  #     "Select X variables",
  #     choices = candidates,
  #     multiple = TRUE,
  #     options = list(placeholder = "Choose one or more regressors…")
  #   )
  # })
  # 
  # 
  # # Auto-SARIMAX selector UI
  # output$autox_xreg_ui <- renderUI({
  #   req(prepared())
  #   df <- prepared()$df
  #   candidates <- setdiff(names(df), c(prepared()$date_col, prepared()$value_col))
  #   candidates <- candidates[candidates != ""]
  #   selectizeInput(
  #     "autox_x_cols",
  #     "Select X variables",
  #     choices = candidates,
  #     multiple = TRUE,
  #     options = list(placeholder = "Choose one or more regressors…")
  #   )
  # })
  # 
  # 
  #   #   3) SARIMAX (Manual): fit, equation, diagnostics, forecast
  #   
  # # a) Fit
  # sarimax_fit <- eventReactive(input$fit_sarimax, {
  #   req(ts_train_test(), prepared())
  #   
  #   s <- ts_train_test()
  #   y_train <- as.numeric(s$ts_train)
  #   y_test  <- as.numeric(s$ts_test)
  #   
  #   train_n <- length(y_train)
  #   test_n  <- length(y_test)
  #   
  #   df <- prepared()$df
  #   xcols <- input$sarimax_x_cols %||% character(0)
  #   
  #   xs <- build_xreg_split(df, xcols, train_n, test_n, scale_x = isTRUE(input$sarimax_scale_x))
  #   
  #   # seasonal period
  #   s_in <- suppressWarnings(as.integer(input$sx_s))
  #   if (!is.finite(s_in)) s_in <- as.integer(prepared()$freq)
  #   if (!is.finite(s_in) || s_in < 1) s_in <- 1L
  #   
  #   fit <- forecast::Arima(
  #     y_train,
  #     order = c(as.integer(input$sx_p), as.integer(input$sx_d), as.integer(input$sx_q)),
  #     seasonal = list(order = c(as.integer(input$sx_P), as.integer(input$sx_D), as.integer(input$sx_Q)), period = s_in),
  #     xreg = xs$x_train,
  #     include.mean = isTRUE(input$sarimax_drift),
  #     include.drift = isTRUE(input$sarimax_drift),
  #     method = "ML"
  #   )
  #   
  #   list(
  #     fit = fit,
  #     xcols = xcols,
  #     x_train = xs$x_train,
  #     x_test = xs$x_test,
  #     y_train = y_train,
  #     y_test = y_test,
  #     period = s_in
  #   )
  # })
  # 
  # 
  # # b) Model spec + coefs
  # output$sarimax_model_spec <- renderPrint({
  #   validate(need(input$fit_sarimax > 0, "Click Fit to estimate SARIMAX."))
  #   req(sarimax_fit())
  #   print(sarimax_fit()$fit)
  # })
  # 
  # output$sarimax_coef_table <- renderTable({
  #   req(sarimax_fit())
  #   sm <- summary(sarimax_fit()$fit)
  #   co <- as.data.frame(sm$coef)
  #   co$term <- rownames(co)
  #   rownames(co) <- NULL
  #   co <- co[, c("term", names(co)[1:(ncol(co)-1)]), drop = FALSE]
  #   co
  # }, digits = 5)
  # 
  # 
  # 
  # # c) Equation (MathJax)
  # output$sarimax_model_equation <- renderUI({
  #   validate(need(input$fit_sarimax > 0, "Fit SARIMAX to show the equation."))
  #   req(sarimax_fit())
  #   obj <- sarimax_fit()
  #   xcols <- obj$xcols
  #   
  #   # build LaTeX for regression part
  #   x_part <- if (length(xcols) == 0) {
  #     "0"
  #   } else {
  #     paste0("\\sum_{k=1}^{", length(xcols), "} \\beta_k x_{k,t}")
  #   }
  #   
  #   # mean equation (generic SARIMAX)
  #   mean_eq <- paste0(
  #     "y_t = c + ", x_part, " + \\varepsilon_t"
  #   )
  #   
  #   # full SARIMA operator form (generic)
  #   op_eq <- paste0(
  #     "\\Phi(B^s)\\phi(B)(1-B)^d(1-B^s)^D y_t = c + ",
  #     x_part, " + \\Theta(B^s)\\theta(B)\\varepsilon_t"
  #   )
  #   
  #   html <- paste0(
  #     "<p><b>Regression mean equation:</b></p><div>$$", mean_eq, "$$</div>",
  #     "<p><b>SARIMAX operator form:</b></p><div>$$", op_eq, "$$</div>"
  #   )
  #   
  #   session$onFlushed(function() {
  #     session$sendCustomMessage("mathjax-typeset", "sarimax_eq_box")
  #   }, once = TRUE)
  #   
  #   HTML(html)
  # })
  # 
  # 
  # 
  # # d) Diagnostics plots + residual tests
  # sarimax_resid <- reactive({
  #   req(sarimax_fit())
  #   as.numeric(residuals(sarimax_fit()$fit))
  # })
  # 
  # output$sarimax_resid_ts <- renderPlot({
  #   validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
  #   r <- sarimax_resid()
  #   plot(r, type = "l", main = "Residuals", ylab = "e_t", xlab = "t")
  #   abline(h = 0, col = "red", lty = 2)
  # })
  # 
  # output$sarimax_resid_acf <- renderPlot({
  #   validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
  #   stats::acf(sarimax_resid(), main = "ACF of residuals")
  # })
  # 
  # output$sarimax_resid_hist <- renderPlot({
  #   validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
  #   hist(sarimax_resid(), breaks = 30, main = "Residual histogram", xlab = "e_t")
  # })
  # 
  # output$sarimax_resid_qq <- renderPlot({
  #   validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
  #   r <- sarimax_resid()
  #   qqnorm(r, main = "Q–Q plot"); qqline(r, col = "red")
  # })
  # 
  # output$sarimax_resid_lb_pvals <- renderPlot({
  #   validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
  #   r <- sarimax_resid()
  #   maxL <- max(5, as.integer(input$diag_lag))
  #   lags <- 1:maxL
  #   pvals <- sapply(lags, function(L) {
  #     out <- tryCatch(Box.test(r, lag = L, type = "Ljung-Box"), error = function(e) NULL)
  #     if (is.null(out)) NA_real_ else out$p.value
  #   })
  #   plot(lags, pvals, type = "b", pch = 16, ylim = c(0,1),
  #        main = "Ljung–Box p-values by lag", xlab = "Lag", ylab = "p-value")
  #   abline(h = 0.05, col = "red", lty = 2)
  # })
  # 
  # output$sarimax_diag_tests <- renderPrint({
  #   validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
  #   r <- sarimax_resid()
  #   L <- as.integer(input$diag_lag); if (!is.finite(L) || L < 1) L <- 12L
  #   cat("Residual tests (SARIMAX)\n=======================\n\n")
  #   print(Box.test(r, lag = L, type = "Ljung-Box"))
  #   cat("\n")
  #   if (requireNamespace("tseries", quietly = TRUE)) {
  #     print(tseries::jarque.bera.test(r))
  #     cat("\n")
  #   }
  # })
  # 
  # # e) Forecast & accuracy
  # sarimax_forecast <- reactive({
  #   req(sarimax_fit())
  #   obj <- sarimax_fit()
  #   
  #   n_test <- length(obj$y_test)
  #   if (n_test > 0) {
  #     h <- n_test
  #   } else {
  #     h_in <- suppressWarnings(as.integer(input$sarimax_h))
  #     if (!is.finite(h_in) || h_in < 1) h_in <- 20L
  #     h <- h_in
  #   }
  #   
  #   xfuture <- if (is.null(obj$x_test) || n_test == 0) {
  #     # future mode: require user-provided future X? If absent, set zeros.
  #     if (is.null(obj$x_train)) NULL else matrix(0, nrow = h, ncol = ncol(obj$x_train),
  #                                                dimnames = list(NULL, colnames(obj$x_train)))
  #   } else {
  #     obj$x_test
  #   }
  #   
  #   fc <- forecast::forecast(obj$fit, h = h, xreg = xfuture)
  #   list(fc = fc, h = h, xfuture = xfuture)
  # })
  # 
  # output$sarimax_horizon_note <- renderPrint({
  #   validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
  #   obj <- sarimax_fit()
  #   if (length(obj$y_test) > 0) {
  #     cat("Holdout test detected. Forecast horizon h was set to test length (h =", length(obj$y_test), ").\n")
  #   } else {
  #     cat("No test set detected. Forecast horizon uses sarimax_h (or default).\n")
  #     cat("Note: future X values are set to 0 unless you implement a future-X input.\n")
  #   }
  # })
  # 
  # output$sarimax_forecast_plot <- renderPlot({
  #   validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
  #   plot(sarimax_forecast()$fc, main = "SARIMAX forecast")
  # })
  # 
  # output$sarimax_forecast_table <- renderTable({
  #   validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
  #   fc <- sarimax_forecast()$fc
  #   out <- data.frame(
  #     t = seq_along(fc$mean),
  #     mean = as.numeric(fc$mean),
  #     lo80 = as.numeric(fc$lower[,1]),
  #     hi80 = as.numeric(fc$upper[,1]),
  #     lo95 = as.numeric(fc$lower[,2]),
  #     hi95 = as.numeric(fc$upper[,2])
  #   )
  #   out
  # }, digits = 6)
  # 
  # output$sarimax_accuracy_table <- renderTable({
  #   validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
  #   obj <- sarimax_fit()
  #   if (length(obj$y_test) == 0) return(NULL)
  #   
  #   y_test <- obj$y_test
  #   y_hat <- as.numeric(sarimax_forecast()$fc$mean)
  #   accuracy_df(y_test, y_hat)
  # }, rownames = FALSE)
  # 
  # output$apa_sarimax_paragraph <- renderPrint({
  #   validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
  #   obj <- sarimax_fit()
  #   fit <- obj$fit
  #   sm <- summary(fit)
  #   cat("A SARIMAX model was estimated with exogenous regressors (xreg) to model the mean structure while capturing seasonal and non-seasonal dynamics. ",
  #       "The fitted specification was ", as.character(fit),
  #       ". Model adequacy was assessed via residual diagnostics (ACF and Ljung–Box), and forecasts were generated using the available exogenous information.",
  #       sep = "")
  # })
  # 
  # 
  # 
  # 
  # #==============================================================
  # # 4) Auto-SARIMAX: fit + outputs (parallel to Auto-ARIMA)
  # #==============================================================
  # 
  # # a) Fit
  # autox_fit <- eventReactive(input$fit_autox, {
  #   req(ts_train_test(), prepared())
  #   
  #   s <- ts_train_test()
  #   y_train <- as.numeric(s$ts_train)
  #   y_test  <- as.numeric(s$ts_test)
  #   
  #   train_n <- length(y_train)
  #   test_n  <- length(y_test)
  #   
  #   df <- prepared()$df
  #   xcols <- input$autox_x_cols %||% character(0)
  #   xs <- build_xreg_split(df, xcols, train_n, test_n, scale_x = isTRUE(input$autox_scale_x))
  #   
  #   fit <- forecast::auto.arima(
  #     y_train,
  #     xreg = xs$x_train,
  #     seasonal = isTRUE(input$autox_seasonal),
  #     stepwise = isTRUE(input$autox_stepwise),
  #     approximation = isTRUE(input$autox_approx),
  #     allowmean = isTRUE(input$autox_allow_mean),
  #     allowdrift = isTRUE(input$autox_allow_mean),
  #     max.order = as.integer(input$autox_max_order)
  #   )
  #   
  #   list(
  #     fit = fit,
  #     xcols = xcols,
  #     x_train = xs$x_train,
  #     x_test = xs$x_test,
  #     y_train = y_train,
  #     y_test = y_test
  #   )
  # })
  # 
  # 
  # # b) Outputs
  # output$autox_model_spec <- renderPrint({
  #   validate(need(input$fit_autox > 0, "Click Fit Auto-SARIMAX first."))
  #   req(autox_fit())
  #   print(autox_fit()$fit)
  # })
  # 
  # output$autox_coef_table <- renderTable({
  #   req(autox_fit())
  #   sm <- summary(autox_fit()$fit)
  #   co <- as.data.frame(sm$coef)
  #   co$term <- rownames(co)
  #   rownames(co) <- NULL
  #   co
  # }, digits = 5)
  # 
  # output$autox_model_equation <- renderUI({
  #   validate(need(input$fit_autox > 0, "Fit Auto-SARIMAX to show the equation."))
  #   req(autox_fit())
  #   obj <- autox_fit()
  #   xcols <- obj$xcols
  #   
  #   x_part <- if (length(xcols) == 0) "0" else paste0("\\sum_{k=1}^{", length(xcols), "} \\beta_k x_{k,t}")
  #   mean_eq <- paste0("y_t = c + ", x_part, " + \\varepsilon_t")
  #   
  #   html <- paste0(
  #     "<p><b>Regression mean equation:</b></p><div>$$", mean_eq, "$$</div>",
  #     "<p><b>Selected ARIMA structure (from auto.arima):</b> ", as.character(obj$fit), "</p>"
  #   )
  #   
  #   session$onFlushed(function() {
  #     session$sendCustomMessage("mathjax-typeset", "autox_eq_box")
  #   }, once = TRUE)
  #   
  #   HTML(html)
  # })
  
  
  

  
  
  # =============================
  # Auto-SARIMAX (auto.arima + xreg)
  # =============================
  
  output$step5b_notes <- renderUI({
    tags$div(
      tags$ol(
        tags$li("Select one or more exogenous regressors (X)."),
        tags$li("Click “Fit Auto-SARIMAX”."),
        tags$li("Inspect model specification and coefficients."),
        tags$li("Check residual diagnostics and formal tests."),
        tags$li("Review forecasts and (if available) test-set accuracy."),
        tags$li(tags$b("Note:"), " If train_prop = 1 (no test set), future X values are assumed 0 unless you add a future-X input.")
      )
    )
  })
  
  # ---- X selector UI (uses prepared()$df and excludes date/value columns if available)
  output$autox_xreg_ui <- renderUI({
    req(prepared())
    df <- prepared()$df
    
    date_col  <- prepared()$date_col %||% character(0)
    value_col <- prepared()$value_col %||% character(0)
    
    candidates <- setdiff(names(df), c(date_col, value_col))
    candidates <- candidates[nzchar(candidates)]
    
    selectizeInput(
      "autox_x_cols",
      "Select X variables",
      choices = candidates,
      multiple = TRUE,
      options = list(placeholder = "Choose one or more regressors…")
    )
  })
  
  # ---- Fit model (cached on Fit button)
  autox_fit <- eventReactive(input$fit_autox, {
    req(ts_train_test(), prepared())
    validate(need(requireNamespace("forecast", quietly = TRUE), "Package 'forecast' is required for Auto-SARIMAX."))
    
    s <- ts_train_test()
    y_train <- as.numeric(s$ts_train)
    y_test  <- as.numeric(s$ts_test)
    
    train_n <- length(y_train)
    test_n  <- length(y_test)
    
    df <- prepared()$df
    xcols <- input$autox_x_cols %||% character(0)
    
    xs <- build_xreg_split(
      df = df,
      cols = xcols,
      train_n = train_n,
      test_n  = test_n,
      scale_x = isTRUE(input$autox_scale_x)
    )
    
    fit <- forecast::auto.arima(
      y_train,
      xreg = xs$x_train,
      seasonal = isTRUE(input$autox_seasonal),
      stepwise = isTRUE(input$autox_stepwise),
      approximation = isTRUE(input$autox_approx),
      allowmean = isTRUE(input$autox_allow_mean),
      allowdrift = isTRUE(input$autox_allow_mean),
      max.order = as.integer(input$autox_max_order)
    )
    
    list(
      fit = fit,
      xcols = xcols,
      x_train = xs$x_train,
      x_test  = xs$x_test,
      y_train = y_train,
      y_test  = y_test
    )
  })
  
  # ---- Spec
  output$autox_model_spec <- renderPrint({
    validate(need(input$fit_autox > 0, "Click “Fit Auto-SARIMAX” first."))
    req(autox_fit())
    print(autox_fit()$fit)
  })
  
  # ---- Coefs
  output$autox_coef_table <- renderTable({
    req(autox_fit())
    sm <- summary(autox_fit()$fit)
    co <- as.data.frame(sm$coef)
    co$Term <- rownames(co)
    rownames(co) <- NULL
    co <- co[, c("Term", setdiff(names(co), "Term")), drop = FALSE]
    co
  }, digits = 6)
  
  # ---- Equation (MathJax, rendered like GARCH style)
  output$autox_model_equation <- renderUI({
    validate(need(input$fit_autox > 0, "Fit Auto-SARIMAX to show the equation."))
    req(autox_fit())
    obj <- autox_fit()
    
    xcols <- obj$xcols
    x_part <- if (length(xcols) == 0) {
      "0"
    } else {
      paste0("\\sum_{k=1}^{", length(xcols), "} \\beta_k x_{k,t}")
    }
    
    mean_eq <- paste0("y_t = c + ", x_part, " + \\varepsilon_t")
    
    html <- paste0(
      "<p><b>Mean (regression) equation:</b></p><div>$$", mean_eq, "$$</div>",
      "<p><b>Selected ARIMA structure (auto.arima):</b> ", as.character(obj$fit), "</p>"
    )
    
    session$onFlushed(function() {
      session$sendCustomMessage("mathjax-typeset", "autox_eq_box")
    }, once = TRUE)
    
    HTML(html)
  })
  
  # ---- Residuals
  autox_resid <- reactive({
    req(autox_fit())
    as.numeric(residuals(autox_fit()$fit))
  })
  
  # ---- Diagnostics plots
  output$autox_resid_ts <- renderPlot({
    validate(need(input$fit_autox > 0, "Fit Auto-SARIMAX first."))
    r <- autox_resid()
    plot(r, type = "l", main = "Residuals", ylab = "e_t", xlab = "t")
    abline(h = 0, col = "red", lty = 2)
  })
  
  output$autox_resid_acf <- renderPlot({
    validate(need(input$fit_autox > 0, "Fit Auto-SARIMAX first."))
    stats::acf(autox_resid(), main = "ACF of residuals")
  })
  
  output$autox_resid_hist <- renderPlot({
    validate(need(input$fit_autox > 0, "Fit Auto-SARIMAX first."))
    hist(autox_resid(), breaks = 30, main = "Residual histogram", xlab = "e_t")
  })
  
  output$autox_resid_qq <- renderPlot({
    validate(need(input$fit_autox > 0, "Fit Auto-SARIMAX first."))
    r <- autox_resid()
    qqnorm(r, main = "Q–Q plot")
    qqline(r, col = "red")
  })
  
  output$autox_resid_lb_pvals <- renderPlot({
    validate(need(input$fit_autox > 0, "Fit Auto-SARIMAX first."))
    r <- autox_resid()
    maxL <- suppressWarnings(as.integer(input$diag_lag))
    if (!is.finite(maxL) || maxL < 5) maxL <- 12L
    lags <- 1:maxL
    pvals <- sapply(lags, function(L) {
      out <- tryCatch(Box.test(r, lag = L, type = "Ljung-Box"), error = function(e) NULL)
      if (is.null(out)) NA_real_ else out$p.value
    })
    plot(lags, pvals, type = "b", pch = 16, ylim = c(0, 1),
         main = "Ljung–Box p-values by lag", xlab = "Lag", ylab = "p-value")
    abline(h = 0.05, col = "red", lty = 2)
  })
  
  # ---- Residual tests
  output$autox_diag_tests <- renderPrint({
    validate(need(input$fit_autox > 0, "Fit Auto-SARIMAX first."))
    r <- autox_resid()
    
    L <- suppressWarnings(as.integer(input$diag_lag))
    if (!is.finite(L) || L < 1) L <- 12L
    
    cat("Residual tests (Auto-SARIMAX)\n========================\n\n")
    cat("Ljung–Box:\n")
    print(Box.test(r, lag = L, type = "Ljung-Box"))
    cat("\n")
    
    if (requireNamespace("tseries", quietly = TRUE)) {
      cat("Jarque–Bera normality test:\n")
      print(tseries::jarque.bera.test(r))
      cat("\n")
    } else {
      cat("Jarque–Bera: package 'tseries' not installed.\n\n")
    }
  })
  
  # ---- Forecast
  autox_forecast <- reactive({
    req(autox_fit())
    validate(need(requireNamespace("forecast", quietly = TRUE), "Package 'forecast' is required."))
    
    obj <- autox_fit()
    test_n <- length(obj$y_test)
    
    if (test_n > 0) {
      h <- test_n
      xfuture <- obj$x_test
    } else {
      h_in <- suppressWarnings(as.integer(input$autox_h))
      if (!is.finite(h_in) || h_in < 1) h_in <- 20L
      h <- h_in
      
      # No future X UI -> default zeros matrix with same columns as training X
      xfuture <- if (is.null(obj$x_train)) NULL else
        matrix(0, nrow = h, ncol = ncol(obj$x_train), dimnames = list(NULL, colnames(obj$x_train)))
    }
    
    fc <- forecast::forecast(obj$fit, h = h, xreg = xfuture)
    list(fc = fc, h = h)
  })
  
  output$autox_horizon_note <- renderPrint({
    validate(need(input$fit_autox > 0, "Fit Auto-SARIMAX first."))
    obj <- autox_fit()
    if (length(obj$y_test) > 0) {
      cat("Holdout test detected. Forecast horizon set to test length (h =", length(obj$y_test), ").\n")
    } else {
      cat("No test detected. Forecast horizon uses autox_h (or default). Future X is assumed 0 unless you add a future-X input.\n")
    }
  })
  
  output$autox_forecast_plot <- renderPlot({
    validate(need(input$fit_autox > 0, "Fit Auto-SARIMAX first."))
    plot(autox_forecast()$fc, main = "Auto-SARIMAX forecast")
  })
  
  output$autox_forecast_table <- renderTable({
    validate(need(input$fit_autox > 0, "Fit Auto-SARIMAX first."))
    fc <- autox_forecast()$fc
    data.frame(
      step = seq_along(fc$mean),
      mean = as.numeric(fc$mean),
      lo80 = as.numeric(fc$lower[, 1]),
      hi80 = as.numeric(fc$upper[, 1]),
      lo95 = as.numeric(fc$lower[, 2]),
      hi95 = as.numeric(fc$upper[, 2])
    )
  }, digits = 6)
  
  output$autox_accuracy_table <- renderTable({
    validate(need(input$fit_autox > 0, "Fit Auto-SARIMAX first."))
    obj <- autox_fit()
    if (length(obj$y_test) == 0) return(NULL)
    
    y_test <- obj$y_test
    y_hat  <- as.numeric(autox_forecast()$fc$mean)
    
    if (exists("accuracy_df", mode = "function")) {
      accuracy_df(y_test, y_hat)
    } else {
      e <- y_test - y_hat
      data.frame(
        Metric = c("RMSE", "MAE"),
        Value  = c(sqrt(mean(e^2, na.rm = TRUE)), mean(abs(e), na.rm = TRUE))
      )
    }
  }, rownames = FALSE)
  
  output$apa_autox_paragraph <- renderPrint({
    validate(need(input$fit_autox > 0, "Fit Auto-SARIMAX first."))
    obj <- autox_fit()
    cat(
      "An Auto-SARIMAX model was estimated using forecast::auto.arima with exogenous regressors (xreg). ",
      "The selected specification was ", as.character(obj$fit),
      ". Residual diagnostics (plots and Ljung–Box testing) were used to assess adequacy, and forecasts were generated conditional on the available xreg information.",
      sep = ""
    )
  })
  
  
  
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  
  # =============================
  # Manual SARIMAX (Arima + xreg)
  # =============================
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b  # only if you don't already have it
  
  output$step6b_notes <- renderUI({
    tags$div(
      tags$ol(
        tags$li("Select one or more exogenous regressors (X)."),
        tags$li("Choose the SARIMA orders (p,d,q)(P,D,Q)[s]."),
        tags$li("Click Fit to estimate SARIMAX."),
        tags$li("Check diagnostics and residual tests."),
        tags$li("Evaluate forecasts and (if available) test-set accuracy."),
        tags$li(tags$b("Note:"), " If train_prop = 1 (no test set), future X values are assumed 0 unless you add a future-X input.")
      )
    )
  })
  
  # ---- X selector UI
  output$sarimax_xreg_ui <- renderUI({
    req(prepared())
    df <- prepared()$df
    
    date_col  <- prepared()$date_col %||% character(0)
    value_col <- prepared()$value_col %||% character(0)
    
    candidates <- setdiff(names(df), c(date_col, value_col))
    candidates <- candidates[nzchar(candidates)]
    
    selectizeInput(
      "sarimax_x_cols",
      "Select X variables",
      choices = candidates,
      multiple = TRUE,
      options = list(placeholder = "Choose one or more regressors…")
    )
  })
  
  # ---- Split text/plot (simple; adjust if you already have a nicer split plot helper)
  output$sarimax_split_text <- renderPrint({
    req(ts_train_test())
    s <- ts_train_test()
    cat("Training length:", length(s$ts_train), "\n")
    cat("Test length:", length(s$ts_test), "\n")
  })
  
  output$sarimax_split_plot <- renderPlot({
    req(ts_train_test())
    s <- ts_train_test()
    y <- as.numeric(s$ts_full)
    n_train <- length(s$ts_train)
    
    plot(y, type = "l", main = "Train/Test split", xlab = "t", ylab = "y")
    abline(v = n_train, col = "red", lty = 2)
    legend("topleft", legend = c("Series", "Train/Test boundary"), col = c("black", "red"), lty = c(1,2), bty = "n")
  })
  
  # ---- Fit SARIMAX (cached)
  sarimax_fit <- eventReactive(input$fit_sarimax, {
    req(ts_train_test(), prepared())
    validate(need(requireNamespace("forecast", quietly = TRUE), "Package 'forecast' is required for SARIMAX."))
    
    s <- ts_train_test()
    y_train <- as.numeric(s$ts_train)
    y_test  <- as.numeric(s$ts_test)
    
    train_n <- length(y_train)
    test_n  <- length(y_test)
    
    df <- prepared()$df
    xcols <- input$sarimax_x_cols %||% character(0)
    
    xs <- build_xreg_split(
      df = df,
      cols = xcols,
      train_n = train_n,
      test_n  = test_n,
      scale_x = isTRUE(input$sarimax_scale_x)
    )
    
    # seasonal period: sx_s overrides sidebar frequency
    s_in <- suppressWarnings(as.integer(input$sx_s))
    if (!is.finite(s_in)) s_in <- suppressWarnings(as.integer(prepared()$freq))
    if (!is.finite(s_in) || s_in < 1) s_in <- 1L
    
    p <- as.integer(input$sx_p); if (!is.finite(p) || p < 0) p <- 0L
    d <- as.integer(input$sx_d); if (!is.finite(d) || d < 0) d <- 0L
    q <- as.integer(input$sx_q); if (!is.finite(q) || q < 0) q <- 0L
    
    P <- as.integer(input$sx_P); if (!is.finite(P) || P < 0) P <- 0L
    D <- as.integer(input$sx_D); if (!is.finite(D) || D < 0) D <- 0L
    Q <- as.integer(input$sx_Q); if (!is.finite(Q) || Q < 0) Q <- 0L
    
    drift <- isTRUE(input$sarimax_drift)
    
    fit <- forecast::Arima(
      y_train,
      order = c(p, d, q),
      seasonal = list(order = c(P, D, Q), period = s_in),
      xreg = xs$x_train,
      include.mean = drift,
      include.drift = drift,
      method = "ML"
    )
    
    list(
      fit = fit,
      xcols = xcols,
      x_train = xs$x_train,
      x_test  = xs$x_test,
      y_train = y_train,
      y_test  = y_test,
      period  = s_in,
      orders  = list(p=p,d=d,q=q,P=P,D=D,Q=Q)
    )
  })
  
  # ---- Model spec
  output$sarimax_model_spec <- renderPrint({
    validate(need(input$fit_sarimax > 0, "Click Fit SARIMAX first."))
    req(sarimax_fit())
    print(sarimax_fit()$fit)
  })
  
  # ---- Coef table
  output$sarimax_coef_table <- renderTable({
    req(sarimax_fit())
    sm <- summary(sarimax_fit()$fit)
    co <- as.data.frame(sm$coef)
    co$Term <- rownames(co)
    rownames(co) <- NULL
    co <- co[, c("Term", setdiff(names(co), "Term")), drop = FALSE]
    co
  }, digits = 6)
  
  # ---- Equation (MathJax)
  output$sarimax_model_equation <- renderUI({
    validate(need(input$fit_sarimax > 0, "Fit SARIMAX to show the equation."))
    req(sarimax_fit())
    obj <- sarimax_fit()
    
    k <- length(obj$xcols)
    
    x_part <- if (k == 0) {
      "0"
    } else {
      paste0("\\sum_{k=1}^{", k, "} \\beta_k x_{k,t}")
    }
    
    mean_eq <- paste0("y_t = c + ", x_part, " + \\varepsilon_t")
    
    # Operator form for SARIMAX (generic)
    op_eq <- paste0(
      "\\Phi(B^s)\\phi(B)(1-B)^d(1-B^s)^D y_t = c + ",
      x_part,
      " + \\Theta(B^s)\\theta(B)\\varepsilon_t"
    )
    
    html <- paste0(
      "<p><b>Regression mean equation:</b></p><div>$$", mean_eq, "$$</div>",
      "<p><b>SARIMAX operator form:</b></p><div>$$", op_eq, "$$</div>"
    )
    
    session$onFlushed(function() {
      session$sendCustomMessage("mathjax-typeset", "sarimax_eq_box")
    }, once = TRUE)
    
    HTML(html)
  })
  
  # ---- Residuals
  sarimax_resid <- reactive({
    req(sarimax_fit())
    as.numeric(residuals(sarimax_fit()$fit))
  })
  
  # ---- Diagnostics plots
  output$sarimax_resid_ts <- renderPlot({
    validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
    r <- sarimax_resid()
    plot(r, type = "l", main = "Residuals", ylab = "e_t", xlab = "t")
    abline(h = 0, col = "red", lty = 2)
  })
  
  output$sarimax_resid_acf <- renderPlot({
    validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
    stats::acf(sarimax_resid(), main = "ACF of residuals")
  })
  
  output$sarimax_resid_hist <- renderPlot({
    validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
    hist(sarimax_resid(), breaks = 30, main = "Residual histogram", xlab = "e_t")
  })
  
  output$sarimax_resid_qq <- renderPlot({
    validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
    r <- sarimax_resid()
    qqnorm(r, main = "Q–Q plot"); qqline(r, col = "red")
  })
  
  output$sarimax_resid_lb_pvals <- renderPlot({
    validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
    r <- sarimax_resid()
    maxL <- suppressWarnings(as.integer(input$diag_lag))
    if (!is.finite(maxL) || maxL < 5) maxL <- 12L
    lags <- 1:maxL
    pvals <- sapply(lags, function(L) {
      out <- tryCatch(Box.test(r, lag = L, type = "Ljung-Box"), error = function(e) NULL)
      if (is.null(out)) NA_real_ else out$p.value
    })
    plot(lags, pvals, type = "b", pch = 16, ylim = c(0, 1),
         main = "Ljung–Box p-values by lag", xlab = "Lag", ylab = "p-value")
    abline(h = 0.05, col = "red", lty = 2)
  })
  
  # ---- Residual tests
  output$sarimax_diag_tests <- renderPrint({
    validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
    r <- sarimax_resid()
    
    L <- suppressWarnings(as.integer(input$diag_lag))
    if (!is.finite(L) || L < 1) L <- 12L
    
    cat("Residual tests (Manual SARIMAX)\n==========================\n\n")
    
    cat("Ljung–Box:\n")
    print(Box.test(r, lag = L, type = "Ljung-Box"))
    cat("\n")
    
    if (requireNamespace("tseries", quietly = TRUE)) {
      cat("Jarque–Bera normality test:\n")
      print(tseries::jarque.bera.test(r))
      cat("\n")
    } else {
      cat("Jarque–Bera: package 'tseries' not installed.\n\n")
    }
  })
  
  # ---- Forecast
  sarimax_forecast <- reactive({
    req(sarimax_fit())
    validate(need(requireNamespace("forecast", quietly = TRUE), "Package 'forecast' is required."))
    
    obj <- sarimax_fit()
    test_n <- length(obj$y_test)
    
    if (test_n > 0) {
      h <- test_n
      xfuture <- obj$x_test
    } else {
      h_in <- suppressWarnings(as.integer(input$sarimax_h))
      if (!is.finite(h_in) || h_in < 1) h_in <- 20L
      h <- h_in
      
      # No future X UI -> default zeros
      xfuture <- if (is.null(obj$x_train)) NULL else
        matrix(0, nrow = h, ncol = ncol(obj$x_train), dimnames = list(NULL, colnames(obj$x_train)))
    }
    
    fc <- forecast::forecast(obj$fit, h = h, xreg = xfuture)
    list(fc = fc, h = h)
  })
  
  output$sarimax_horizon_note <- renderPrint({
    validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
    obj <- sarimax_fit()
    if (length(obj$y_test) > 0) {
      cat("Holdout test detected. Forecast horizon set to test length (h =", length(obj$y_test), ").\n")
    } else {
      cat("No test detected. Forecast horizon uses sarimax_h (or default). Future X is assumed 0 unless you add a future-X input.\n")
    }
  })
  
  output$sarimax_forecast_plot <- renderPlot({
    validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
    plot(sarimax_forecast()$fc, main = "SARIMAX forecast")
  })
  
  output$sarimax_forecast_table <- renderTable({
    validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
    fc <- sarimax_forecast()$fc
    data.frame(
      step = seq_along(fc$mean),
      mean = as.numeric(fc$mean),
      lo80 = as.numeric(fc$lower[, 1]),
      hi80 = as.numeric(fc$upper[, 1]),
      lo95 = as.numeric(fc$lower[, 2]),
      hi95 = as.numeric(fc$upper[, 2])
    )
  }, digits = 6)
  
  output$sarimax_accuracy_table <- renderTable({
    validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
    obj <- sarimax_fit()
    if (length(obj$y_test) == 0) return(NULL)
    
    y_test <- obj$y_test
    y_hat  <- as.numeric(sarimax_forecast()$fc$mean)
    
    if (exists("accuracy_df", mode = "function")) {
      accuracy_df(y_test, y_hat)
    } else {
      e <- y_test - y_hat
      data.frame(
        Metric = c("RMSE", "MAE"),
        Value  = c(sqrt(mean(e^2, na.rm = TRUE)), mean(abs(e), na.rm = TRUE))
      )
    }
  }, rownames = FALSE)
  
  # ---- APA paragraph
  output$apa_sarimax_paragraph <- renderPrint({
    validate(need(input$fit_sarimax > 0, "Fit SARIMAX first."))
    obj <- sarimax_fit()
    fit <- obj$fit
    
    cat(
      "A SARIMAX model (seasonal ARIMA with exogenous regressors) was estimated using forecast::Arima with xreg predictors. ",
      "The specified model was ", as.character(fit),
      if (isTRUE(input$sarimax_drift)) " including a drift/mean term." else " without a drift/mean term.",
      " Model adequacy was assessed using residual plots and Ljung–Box testing, and forecasts were generated conditional on the available exogenous information.",
      sep = ""
    )
  })
  
  
  
  
  
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  #=====================================================================================================
  
  
  
  # ============================
  # ChatGPT-SARIMA (Manual) TAB
  # ============================
  
  # ---- small helpers (safe if you already have them)
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  cg_fmt_num <- function(x, d = 3) {
    if (length(x) == 0 || !is.finite(x)) return(NA_character_)
    formatC(x, format = "f", digits = d)
  }
  cg_fmt_p <- function(p) {
    if (!is.finite(p)) return("= NA")
    if (p < 0.001) "< .001" else paste0("= ", cg_fmt_num(p, 3))
  }
  cg_sig_stars <- function(p) {
    if (!is.finite(p)) return("")
    if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else if (p < 0.1) "†" else ""
  }
  
  output$cg_notes <- renderUI({
    if (!isTRUE(input$cg_show_teaching)) return(NULL)
    tags$div(
      class = "alert alert-info",
      tags$b("ChatGPT-SARIMA workflow:"),
      tags$ol(
        tags$li("Confirm stationarity/differencing decisions (d, D) on the Stationarity tab."),
        tags$li("Use ACF/PACF patterns to propose p/q and P/Q (seasonal at lag s)."),
        tags$li("Fit the model; then check residuals: they should resemble white noise."),
        tags$li("If residual autocorrelation remains, adjust orders; if volatility clustering exists, consider GARCH.")
      )
    )
  })
  
  # ---- Split summary + plot
  output$cg_split_text <- renderPrint({
    req(ts_train_test())
    s <- ts_train_test()
    cat("Training length:", length(s$ts_train), "\n")
    cat("Test length:", length(s$ts_test), "\n")
  })
  
  output$cg_split_plot <- renderPlot({
    req(ts_train_test())
    s <- ts_train_test()
    y <- as.numeric(s$ts_full)
    n_train <- length(s$ts_train)
    
    plot(y, type = "l", main = "Train/Test split", xlab = "t", ylab = "y")
    abline(v = n_train, col = "red", lty = 2)
    legend("topleft", legend = c("Series", "Train/Test boundary"),
           col = c("black", "red"), lty = c(1, 2), bty = "n")
  })
  
  # ---- Fit SARIMA (cached on button)
  cg_fit <- eventReactive(input$fit_cg_sarima, {
    req(ts_train_test())
    validate(need(requireNamespace("forecast", quietly = TRUE), "Package 'forecast' is required."))
    
    s <- ts_train_test()
    y_train <- as.numeric(s$ts_train)
    y_test  <- as.numeric(s$ts_test)
    
    p <- as.integer(input$cg_p); if (!is.finite(p) || p < 0) p <- 0L
    d <- as.integer(input$cg_d); if (!is.finite(d) || d < 0) d <- 0L
    q <- as.integer(input$cg_q); if (!is.finite(q) || q < 0) q <- 0L
    
    P <- as.integer(input$cg_P); if (!is.finite(P) || P < 0) P <- 0L
    D <- as.integer(input$cg_D); if (!is.finite(D) || D < 0) D <- 0L
    Q <- as.integer(input$cg_Q); if (!is.finite(Q) || Q < 0) Q <- 0L
    
    s_in <- suppressWarnings(as.integer(input$cg_s))
    if (!is.finite(s_in)) {
      # try prepared()$freq if available, else default 1
      s_in <- tryCatch(as.integer(prepared()$freq), error = function(e) NA_integer_)
    }
    if (!is.finite(s_in) || s_in < 1) s_in <- 1L
    
    include_mean <- isTRUE(input$cg_include_mean)
    method <- as.character(input$cg_method) %||% "ML"
    
    fit <- forecast::Arima(
      y_train,
      order = c(p, d, q),
      seasonal = list(order = c(P, D, Q), period = s_in),
      include.mean = include_mean,
      include.drift = include_mean,
      method = method,
      biasadj = isTRUE(input$cg_biasadj)
    )
    
    list(
      fit = fit,
      y_train = y_train,
      y_test  = y_test,
      p=p, d=d, q=q, P=P, D=D, Q=Q, s=s_in,
      include_mean = include_mean
    )
  })
  
  # ---- Model spec + IC table
  output$cg_model_spec <- renderPrint({
    validate(need(input$fit_cg_sarima > 0, "Click Fit to estimate ChatGPT-SARIMA."))
    req(cg_fit())
    print(cg_fit()$fit)
  })
  
  output$cg_ic_table <- renderTable({
    req(cg_fit())
    fit <- cg_fit()$fit
    data.frame(
      Metric = c("AIC", "AICc", "BIC", "logLik"),
      Value  = c(
        as.numeric(fit$aic),
        as.numeric(fit$aicc),
        as.numeric(fit$bic),
        as.numeric(stats::logLik(fit))
      )
    )
  }, digits = 4)
  
  # ---- Coef table with stars
  output$cg_coef_table <- renderTable({
    req(cg_fit())
    sm <- summary(cg_fit()$fit)
    co <- as.data.frame(sm$coef)
    co$Term <- rownames(co); rownames(co) <- NULL
    names(co) <- c("Estimate", "Std.Error", "t.value", "Pr(>|t|)", "Term")
    
    co <- co[, c("Term", "Estimate", "Std.Error", "t.value", "Pr(>|t|)"), drop = FALSE]
    co$Estimate <- as.numeric(co$Estimate)
    co$Std.Error <- as.numeric(co$Std.Error)
    co$`t.value` <- as.numeric(co$`t.value`)
    co$`Pr(>|t|)` <- as.numeric(co$`Pr(>|t|)`)
    
    if (isTRUE(input$cg_show_stars)) {
      co$Sig <- vapply(co$`Pr(>|t|)`, cg_sig_stars, character(1))
    }
    
    co
  }, digits = 6)
  
  # ---- Equation (MathJax) — rendered as HTML + forced typeset
  output$cg_model_equation <- renderUI({
    validate(need(input$fit_cg_sarima > 0, "Fit the model to generate the equation."))
    req(cg_fit())
    obj <- cg_fit()
    
    # Generic academic operator form
    # (This avoids fragile “expanded numeric” equations and renders reliably.)
    eq <- paste0(
      "\\Phi(B^{", obj$s, "})\\,\\phi(B)\\,(1-B)^{", obj$d, "}(1-B^{", obj$s, "})^{", obj$D, "}\\,y_t = ",
      if (obj$include_mean) "c + " else "",
      "\\Theta(B^{", obj$s, "})\\,\\theta(B)\\,\\varepsilon_t"
    )
    
    html <- paste0(
      "<p><b>SARIMA operator form:</b></p>",
      "<div>$$", eq, "$$</div>",
      "<p><b>Selected orders:</b> SARIMA(",
      obj$p, ",", obj$d, ",", obj$q, ")(",
      obj$P, ",", obj$D, ",", obj$Q, ")[", obj$s, "]</p>"
    )
    
    session$onFlushed(function() {
      session$sendCustomMessage("mathjax-typeset", "cg_eq_box")
    }, once = TRUE)
    
    HTML(html)
  })
  
  # ---- Residuals
  cg_resid <- reactive({
    req(cg_fit())
    as.numeric(residuals(cg_fit()$fit))
  })
  
  # ---- Diagnostics plots
  output$cg_resid_ts <- renderPlot({
    validate(need(input$fit_cg_sarima > 0, "Fit ChatGPT-SARIMA first."))
    r <- cg_resid()
    plot(r, type = "l", main = "Residuals", ylab = "e_t", xlab = "t")
    abline(h = 0, col = "red", lty = 2)
  })
  
  output$cg_resid_acf <- renderPlot({
    validate(need(input$fit_cg_sarima > 0, "Fit ChatGPT-SARIMA first."))
    stats::acf(cg_resid(), main = "ACF of residuals")
  })
  
  output$cg_resid_hist <- renderPlot({
    validate(need(input$fit_cg_sarima > 0, "Fit ChatGPT-SARIMA first."))
    hist(cg_resid(), breaks = 30, main = "Residual histogram", xlab = "e_t")
  })
  
  output$cg_resid_qq <- renderPlot({
    validate(need(input$fit_cg_sarima > 0, "Fit ChatGPT-SARIMA first."))
    r <- cg_resid()
    qqnorm(r, main = "Q–Q plot"); qqline(r, col = "red")
  })
  
  output$cg_lb_pvals <- renderPlot({
    validate(need(input$fit_cg_sarima > 0, "Fit ChatGPT-SARIMA first."))
    r <- cg_resid()
    maxL <- suppressWarnings(as.integer(input$diag_lag))
    if (!is.finite(maxL) || maxL < 5) maxL <- 12L
    
    lags <- 1:maxL
    pvals <- sapply(lags, function(L) {
      out <- tryCatch(Box.test(r, lag = L, type = "Ljung-Box"), error = function(e) NULL)
      if (is.null(out)) NA_real_ else out$p.value
    })
    
    plot(lags, pvals, type = "b", pch = 16, ylim = c(0, 1),
         main = "Ljung–Box p-values by lag", xlab = "Lag", ylab = "p-value")
    abline(h = 0.05, col = "red", lty = 2)
  })
  
  # ---- Residual tests (text)
  output$cg_resid_tests <- renderPrint({
    validate(need(input$fit_cg_sarima > 0, "Fit ChatGPT-SARIMA first."))
    r <- cg_resid()
    
    L <- suppressWarnings(as.integer(input$diag_lag))
    if (!is.finite(L) || L < 1) L <- 12L
    
    cat("Residual tests (ChatGPT-SARIMA)\n")
    cat("================================\n\n")
    
    cat("Ljung–Box test on residuals:\n")
    print(Box.test(r, lag = L, type = "Ljung-Box"))
    cat("\n")
    
    if (requireNamespace("tseries", quietly = TRUE)) {
      cat("Jarque–Bera normality test:\n")
      print(tseries::jarque.bera.test(r))
      cat("\n")
    } else {
      cat("Jarque–Bera: package 'tseries' not installed.\n\n")
    }
    
    if (requireNamespace("FinTS", quietly = TRUE)) {
      cat("ARCH LM test (FinTS::ArchTest):\n")
      print(FinTS::ArchTest(r, lags = L))
      cat("\n")
    } else {
      cat("ARCH LM: package 'FinTS' not installed.\n\n")
    }
  })
  
  # ---- Forecast + accuracy
  cg_forecast <- reactive({
    req(cg_fit())
    validate(need(requireNamespace("forecast", quietly = TRUE), "Package 'forecast' is required."))
    
    obj <- cg_fit()
    test_n <- length(obj$y_test)
    
    if (test_n > 0) {
      h <- test_n
    } else {
      h_in <- suppressWarnings(as.integer(input$cg_h))
      if (!is.finite(h_in) || h_in < 1) h_in <- 20L
      h <- h_in
    }
    
    fc <- forecast::forecast(obj$fit, h = h)
    list(fc = fc, h = h)
  })
  
  output$cg_horizon_note <- renderPrint({
    validate(need(input$fit_cg_sarima > 0, "Fit ChatGPT-SARIMA first."))
    obj <- cg_fit()
    if (length(obj$y_test) > 0) {
      cat("Holdout test detected. Horizon set to test length (h =", length(obj$y_test), ").\n")
    } else {
      cat("No test detected. Horizon uses cg_h (or default).\n")
    }
  })
  
  output$cg_forecast_plot <- renderPlot({
    validate(need(input$fit_cg_sarima > 0, "Fit ChatGPT-SARIMA first."))
    plot(cg_forecast()$fc, main = "ChatGPT-SARIMA forecast")
  })
  
  output$cg_forecast_table <- renderTable({
    validate(need(input$fit_cg_sarima > 0, "Fit ChatGPT-SARIMA first."))
    fc <- cg_forecast()$fc
    
    out <- data.frame(
      step = seq_along(fc$mean),
      mean = as.numeric(fc$mean)
    )
    
    # intervals if present
    if (!is.null(fc$lower) && !is.null(fc$upper) && isTRUE(input$cg_show_pi)) {
      out$lo80 <- as.numeric(fc$lower[, 1])
      out$hi80 <- as.numeric(fc$upper[, 1])
      out$lo95 <- as.numeric(fc$lower[, 2])
      out$hi95 <- as.numeric(fc$upper[, 2])
    }
    out
  }, digits = 6)
  
  output$cg_accuracy_table <- renderTable({
    validate(need(input$fit_cg_sarima > 0, "Fit ChatGPT-SARIMA first."))
    obj <- cg_fit()
    if (length(obj$y_test) == 0) return(NULL)
    
    y_test <- obj$y_test
    y_hat <- as.numeric(cg_forecast()$fc$mean)
    
    if (exists("accuracy_df", mode = "function")) {
      accuracy_df(y_test, y_hat)
    } else {
      e <- y_test - y_hat
      data.frame(
        Metric = c("RMSE", "MAE", "MAPE"),
        Value = c(
          sqrt(mean(e^2, na.rm = TRUE)),
          mean(abs(e), na.rm = TRUE),
          mean(abs(e / y_test), na.rm = TRUE)
        )
      )
    }
  }, rownames = FALSE)
  
  # ---- Academic conclusion: one UI that includes everything (tables, plots, tests, equation, narrative)
  output$cg_conclusion_ui <- renderUI({
    validate(
      need(input$fit_cg_sarima > 0, "Click Fit to generate the academic conclusion.")
    )
    req(cg_fit())
    
    obj <- cg_fit()
    fit <- obj$fit
    r   <- cg_resid()
    
    # ---- core IC (safe)
    aic  <- suppressWarnings(as.numeric(fit$aic))
    aicc <- suppressWarnings(as.numeric(fit$aicc))
    bic  <- suppressWarnings(as.numeric(fit$bic))
    
    aic_txt  <- if (is.finite(aic))  cg_fmt_num(aic,  2) else "NA"
    aicc_txt <- if (is.finite(aicc)) cg_fmt_num(aicc, 2) else "NA"
    bic_txt  <- if (is.finite(bic))  cg_fmt_num(bic,  2) else "NA"
    
    # ---- coefficient summary for narrative (ROBUST)
    sm <- tryCatch(summary(fit), error = function(e) NULL)
    
    co_raw <- NULL
    if (!is.null(sm)) {
      if (!is.null(sm$coef)) co_raw <- sm$coef
      if (is.null(co_raw) && !is.null(sm$coefficients)) co_raw <- sm$coefficients
    }
    
    co <- if (!is.null(co_raw)) {
      df <- as.data.frame(co_raw)
      df$term <- rownames(df)
      rownames(df) <- NULL
      
      # normalize column names to detect estimate/se/stat/p
      nm  <- names(df)
      nm0 <- tolower(gsub("[^a-z]+", "", nm))
      
      pick <- function(keys) {
        idx <- which(nm0 %in% keys)
        if (length(idx) == 0) NA_integer_ else idx[1]
      }
      
      i_est <- pick(c("estimate", "est", "coef", "value"))
      i_se  <- pick(c("se", "stderror", "stderr"))
      i_st  <- pick(c("tvalue", "tstat", "zvalue", "zstat", "statistic"))
      i_p   <- pick(c("prtt", "prgt", "prgtz", "pvalue", "p"))
      
      out <- data.frame(
        term  = df$term,
        est   = if (!is.na(i_est)) suppressWarnings(as.numeric(df[[i_est]])) else NA_real_,
        se    = if (!is.na(i_se))  suppressWarnings(as.numeric(df[[i_se]]))  else NA_real_,
        stat  = if (!is.na(i_st))  suppressWarnings(as.numeric(df[[i_st]]))  else NA_real_,
        p     = if (!is.na(i_p))   suppressWarnings(as.numeric(df[[i_p]]))   else NA_real_,
        stringsAsFactors = FALSE
      )
      
      out$stars <- vapply(out$p, cg_sig_stars, character(1))
      out
    } else {
      data.frame(term = character(0), est = numeric(0), se = numeric(0), stat = numeric(0), p = numeric(0), stars = character(0))
    }
    
    n_sig <- sum(is.finite(co$p) & co$p < 0.05)
    
    # ---- residual tests quick stats for narrative
    L <- suppressWarnings(as.integer(input$diag_lag))
    if (!is.finite(L) || L < 1) L <- 12L
    
    lb <- tryCatch(Box.test(r, lag = L, type = "Ljung-Box"), error = function(e) NULL)
    
    jb <- if (requireNamespace("tseries", quietly = TRUE)) {
      tryCatch(tseries::jarque.bera.test(r), error = function(e) NULL)
    } else NULL
    
    arch <- if (requireNamespace("FinTS", quietly = TRUE)) {
      tryCatch(FinTS::ArchTest(r, lags = L), error = function(e) NULL)
    } else NULL
    
    # ---- accuracy quick stats (if test)
    acc_txt <- "No holdout test set was available; out-of-sample accuracy was not computed."
    if (length(obj$y_test) > 0) {
      y_test <- obj$y_test
      fc_obj <- cg_forecast()
      y_hat  <- if (!is.null(fc_obj) && !is.null(fc_obj$fc)) as.numeric(fc_obj$fc$mean) else rep(NA_real_, length(y_test))
      
      if (length(y_hat) == length(y_test) && any(is.finite(y_hat))) {
        e <- y_test - y_hat
        rmse <- sqrt(mean(e^2, na.rm = TRUE))
        mae  <- mean(abs(e), na.rm = TRUE)
        acc_txt <- paste0(
          "Out-of-sample performance (test n = ", length(y_test),
          ") was RMSE = ", cg_fmt_num(rmse, 3),
          " and MAE = ", cg_fmt_num(mae, 3), "."
        )
      } else {
        acc_txt <- paste0(
          "A test set was detected (n = ", length(y_test),
          "), but forecast values were not available or not finite; accuracy could not be computed."
        )
      }
    }
    
    # ---- typeset equation in conclusion box
    session$onFlushed(function() {
      session$sendCustomMessage("mathjax-typeset", "cg_conclusion_box")
    }, once = TRUE)
    
    tagList(
      tags$h3("ChatGPT-SARIMA: Academic conclusion"),
      
      tags$div(
        class = "cg-h",
        tags$h4("1. Model objective and specification"),
        tags$p(
          "A manually specified seasonal ARIMA (SARIMA) model was estimated to capture both non-seasonal and seasonal dependence in the series.",
          " The final specification was ",
          tags$b(sprintf(
            "SARIMA(%d,%d,%d)(%d,%d,%d)[%d]%s",
            obj$p, obj$d, obj$q, obj$P, obj$D, obj$Q, obj$s,
            if (isTRUE(obj$include_mean)) " with mean/drift." else " without mean/drift."
          )),
          "."
        ),
        tags$p(HTML(paste0(
          "<b>Information criteria:</b> AIC = ", aic_txt,
          ", AICc = ", aicc_txt,
          ", BIC = ", bic_txt, "."
        )))
      ),
      
      tags$div(
        class = "cg-h",
        tags$h4("2. Coefficient inference and practical interpretation"),
        tags$p(HTML(paste0(
          "Coefficient inference was evaluated using approximate t-tests. ",
          "A total of <b>", n_sig, "</b> parameters were statistically significant at α = .05 (see coefficient table)."
        ))),
        tags$p("Interpretation should prioritize model adequacy (white-noise residuals) and forecast performance rather than isolated p-values.")
      ),
      
      tags$div(
        class = "cg-h",
        tags$h4("3. Model equation"),
        tags$div(
          style = "border:1px solid #e5e5e5; border-radius:6px; background:#fcfcfc; padding:10px;",
          uiOutput("cg_model_equation")
        )
      ),
      
      tags$div(
        class = "cg-h",
        tags$h4("4. Diagnostic evidence"),
        tags$p("Residual diagnostics were examined using time-domain plots, autocorrelation diagnostics, and distributional checks."),
        fluidRow(
          column(6, plotOutput("cg_resid_ts", height = 220)),
          column(6, plotOutput("cg_resid_acf", height = 220))
        ),
        fluidRow(
          column(6, plotOutput("cg_resid_hist", height = 220)),
          column(6, plotOutput("cg_resid_qq", height = 220))
        ),
        plotOutput("cg_lb_pvals", height = 260)
      ),
      
      tags$div(
        class = "cg-h",
        tags$h4("5. Residual tests (formal)"),
        tags$p(
          "The following tests provide formal evidence regarding remaining autocorrelation (Ljung–Box), ",
          "departures from normality (Jarque–Bera), and conditional heteroskedasticity (ARCH LM)."
        ),
        verbatimTextOutput("cg_resid_tests"),
        if (!is.null(lb)) tags$p(HTML(paste0(
          "<b>Ljung–Box:</b> Q(", lb$parameter, ") = ", cg_fmt_num(lb$statistic, 3),
          ", p ", cg_fmt_p(lb$p.value), "."
        ))) else NULL,
        if (!is.null(jb)) tags$p(HTML(paste0(
          "<b>Jarque–Bera:</b> JB = ", cg_fmt_num(jb$statistic, 3),
          ", p ", cg_fmt_p(jb$p.value), "."
        ))) else NULL,
        if (!is.null(arch)) tags$p(HTML(paste0(
          "<b>ARCH LM:</b> TR^2 = ", cg_fmt_num(arch$statistic, 3),
          ", p ", cg_fmt_p(arch$p.value), "."
        ))) else NULL
      ),
      
      tags$div(
        class = "cg-h",
        tags$h4("6. Forecasting and predictive performance"),
        tags$p(acc_txt),
        plotOutput("cg_forecast_plot", height = 380),
        tags$h5("Forecast table"),
        tableOutput("cg_forecast_table"),
        tags$h5("Accuracy table (if test set available)"),
        tableOutput("cg_accuracy_table")
      ),
      
      tags$div(
        class = "cg-h",
        tags$h4("7. Final conclusion"),
        tags$p(
          "Overall, the fitted SARIMA model constitutes an interpretable baseline for linear seasonal dependence.",
          " If residual tests indicate non-white-noise behavior (e.g., significant Ljung–Box), the model should be refined by adjusting AR/MA orders or differencing.",
          " If ARCH effects remain, a conditional variance model (e.g., GARCH) is recommended as an extension for volatility clustering."
        )
      )
    )
  })
  
  
  
  
  # output$cg_conclusion_ui <- renderUI({
  #   validate(
  #     need(input$fit_cg_sarima > 0, "Click Fit to generate the academic conclusion.")
  #   )
  #   req(cg_fit())
  # 
  # obj <- cg_fit()
  # fit <- obj$fit
  # r <- cg_resid()
  # 
  # # core IC
  # aic  <- as.numeric(fit$aic)
  # aicc <- as.numeric(fit$aicc)
  # bic  <- as.numeric(fit$bic)
  # 
  # 
  # # coefficient summary for narrative
  # sm <- summary(fit)
  # co <- as.data.frame(sm$coef)
  # co$term <- rownames(co); rownames(co) <- NULL
  # names(co) <- c("est","se","t","p","term")
  # co$stars <- vapply(co$p, cg_sig_stars, character(1))
  # 
  # n_sig <- sum(is.finite(co$p) & co$p < 0.05)
  # 
  # # residual tests quick stats for narrative
  # L <- suppressWarnings(as.integer(input$diag_lag))
  # if (!is.finite(L) || L < 1) L <- 12L
  # lb <- tryCatch(Box.test(r, lag = L, type = "Ljung-Box"), error = function(e) NULL)
  # 
  # jb <- if (requireNamespace("tseries", quietly = TRUE)) {
  #   tryCatch(tseries::jarque.bera.test(r), error = function(e) NULL)
  # } else NULL
  # 
  # arch <- if (requireNamespace("FinTS", quietly = TRUE)) {
  #   tryCatch(FinTS::ArchTest(r, lags = L), error = function(e) NULL)
  # } else NULL
  # 
  # # accuracy quick stats (if test)
  # acc_txt <- "No holdout test set was available; out-of-sample accuracy was not computed."
  # if (length(obj$y_test) > 0) {
  #   y_test <- obj$y_test
  #   y_hat <- as.numeric(cg_forecast()$fc$mean)
  #   e <- y_test - y_hat
  #   rmse <- sqrt(mean(e^2, na.rm = TRUE))
  #   mae  <- mean(abs(e), na.rm = TRUE)
  #   acc_txt <- paste0("Out-of-sample performance (test n = ", length(y_test),
  #                     ") was RMSE = ", cg_fmt_num(rmse, 3),
  #                     " and MAE = ", cg_fmt_num(mae, 3), ".")
  # }
  # 
  # # typeset equation in conclusion box as well (safe)
  # session$onFlushed(function() {
  #   session$sendCustomMessage("mathjax-typeset", "cg_conclusion_box")
  # }, once = TRUE)
  # 
  # tagList(
  #   tags$h3("ChatGPT-SARIMA: Academic conclusion"),
  #   
  #   tags$div(class="cg-h",
  #            tags$h4("1. Model objective and specification"),
  #            tags$p(
  #              "A manually specified seasonal ARIMA (SARIMA) model was estimated to capture both non-seasonal and seasonal dependence in the series.",
  #              " The final specification was ",
  #              tags$b(sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)[%d]%s",
  #                             obj$p, obj$d, obj$q, obj$P, obj$D, obj$Q, obj$s,
  #                             if (obj$include_mean) " with mean/drift." else " without mean/drift.")),
  #              "."
  #            ),
  #            tags$p(HTML(paste0(
  #              "<b>Information criteria:</b> AIC = ", cg_fmt_num(aic, 2),
  #              ", AICc = ", cg_fmt_num(aicc, 2),
  #              ", BIC = ", cg_fmt_num(bic, 2), "."
  #            )))
  #   ),
  #   
  #   tags$div(class="cg-h",
  #            tags$h4("2. Coefficient inference and practical interpretation"),
  #            tags$p(HTML(paste0(
  #              "Coefficient inference was evaluated using approximate t-tests. ",
  #              "A total of <b>", n_sig, "</b> parameters were statistically significant at α = .05 (see coefficient table)."
  #            ))),
  #            tags$p("Interpretation should prioritize model adequacy (white-noise residuals) and forecast performance rather than isolated p-values.")
  #   ),
  #   
  #   tags$div(class="cg-h",
  #            tags$h4("3. Model equation"),
  #            tags$div(
  #              style="border:1px solid #e5e5e5; border-radius:6px; background:#fcfcfc; padding:10px;",
  #              uiOutput("cg_model_equation")
  #            )
  #   ),
  #   
  #   tags$div(class="cg-h",
  #            tags$h4("4. Diagnostic evidence"),
  #            tags$p("Residual diagnostics were examined using time-domain plots, autocorrelation diagnostics, and distributional checks."),
  #            fluidRow(
  #              column(6, plotOutput("cg_resid_ts", height = 220)),
  #              column(6, plotOutput("cg_resid_acf", height = 220))
  #            ),
  #            fluidRow(
  #              column(6, plotOutput("cg_resid_hist", height = 220)),
  #              column(6, plotOutput("cg_resid_qq", height = 220))
  #            ),
  #            plotOutput("cg_lb_pvals", height = 260)
  #   ),
  #   
  #   tags$div(class="cg-h",
  #            tags$h4("5. Residual tests (formal)"),
  #            tags$p(
  #              "The following tests provide formal evidence regarding remaining autocorrelation (Ljung–Box), ",
  #              "departures from normality (Jarque–Bera), and conditional heteroskedasticity (ARCH LM)."
  #            ),
  #            verbatimTextOutput("cg_resid_tests"),
  #            if (!is.null(lb)) tags$p(HTML(paste0(
  #              "<b>Ljung–Box:</b> Q(", lb$parameter, ") = ", cg_fmt_num(lb$statistic, 3),
  #              ", p ", cg_fmt_p(lb$p.value), "."
  #            ))) else NULL,
  #            if (!is.null(jb)) tags$p(HTML(paste0(
  #              "<b>Jarque–Bera:</b> JB = ", cg_fmt_num(jb$statistic, 3),
  #              ", p ", cg_fmt_p(jb$p.value), "."
  #            ))) else NULL,
  #            if (!is.null(arch)) tags$p(HTML(paste0(
  #              "<b>ARCH LM:</b> TR^2 = ", cg_fmt_num(arch$statistic, 3),
  #              ", p ", cg_fmt_p(arch$p.value), "."
  #            ))) else NULL
  #   ),
  #   
  #   tags$div(class="cg-h",
  #            tags$h4("6. Forecasting and predictive performance"),
  #            tags$p(acc_txt),
  #            plotOutput("cg_forecast_plot", height = 380),
  #            tags$h5("Forecast table"),
  #            tableOutput("cg_forecast_table"),
  #            tags$h5("Accuracy table (if test set available)"),
  #            tableOutput("cg_accuracy_table")
  #   ),
  #   
  #   tags$div(class="cg-h",
  #            tags$h4("7. Final conclusion"),
  #            tags$p(
  #              "Overall, the fitted SARIMA model constitutes an interpretable baseline for linear seasonal dependence.",
  #              " If residual tests indicate non-white-noise behavior (e.g., significant Ljung–Box), the model should be refined by adjusting AR/MA orders or differencing.",
  #              " If ARCH effects remain, a conditional variance model (e.g., GARCH) is recommended as an extension for volatility clustering."
  #            )
  #   )
  # )
  # })

  
  
  
  
  
  
  
}
