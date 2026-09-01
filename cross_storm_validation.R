# ============================================================
# CROSS-STORM VALIDATION FOR SYM-H ARDL MODELS
# Leave-One-Storm-Out (LOSO) cross-validation across 10 storms
# ============================================================
#
# This script keeps the same basic setup as the single-storm script:
#   - 5-minute OMNI data
#   - linear interpolation of Bz and solar-wind speed
#   - direct future SYM-H forecast
#   - current SYM-H as an autoregressive predictor
#   - B-spline-smoothed lag effects for Bz and speed
#
# IMPORTANT:
# Linear interpolation here is retrospective. It can use observations on both
# sides of a missing value, so it should not be described as operational
# real-time preprocessing.
#
# TO CHANGE THE LAG WINDOW:
# Edit ONLY max_lag below.
#   max_lag <- 12  gives lag0 ... lag12  = 60 minutes of history
#   max_lag <- 24  gives lag0 ... lag24  = 120 minutes of history
#   max_lag <- 36  gives lag0 ... lag36  = 180 minutes of history
# Note: max_lag = 24 means 25 lagged Bz/speed values because lag0 is included.
# ============================================================

library(dplyr)
library(splines)
library(zoo)

# ============================================================
# 1. USER SETTINGS
# ============================================================

# Put the 10 CSV files in this folder.
# If they are in the same folder as this R script, leave this as ".".
# If they are in a data folder, use: data_dir <- "data"
data_dir <- "."

# Maximum predictor lag in 5-minute steps.
# This is the main value to change if you want 12, 24, 36, etc.
max_lag <- 24

# B-spline degrees of freedom for the lag-response curve.
bs_df <- 5

# 120-minute-ahead forecast = 24 five-minute steps.
forecast_horizon <- 24

# Evaluate a forecast every 120 minutes = every 24 eligible rows.
refresh_steps <- 24

# Data resolution.
data_interval_minutes <- 5

# Output folders.
results_dir <- "results"
figures_dir <- "figures"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 2. STORM FILES
# ============================================================
# The names on the left become storm IDs in the output tables.
# Change a filename here if your local filename differs.

storm_files <- c(
  "1998-05-04" = "omni_5min_1998-05-02_07_bz_speed_symh_cleaned.csv",
  "2001-03-31" = "omni_5min_2001-03-29_04-03_bz_speed_symh_cleaned.csv",
  "2001-04-12" = "omni_5min_2001-04-10_15_bz_speed_symh_cleaned.csv",
  "2003-11-20" = "omni_5min_2003-11-18_23_bz_speed_symh_cleaned.csv",
  "2004-11-08" = "omni_5min_2004-11-07_12_bz_speed_symh_cleaned.csv",
  "2005-05-15" = "omni_5min_2005-05-13_18_bz_speed_symh_cleaned.csv",
  "2012-07-15" = "omni_5min_2012-07-13_18_bz_speed_symh_cleaned.csv",
  "2015-06-23" = "omni_5min_2015-06-20_25_bz_speed_symh_cleaned.csv",
  "2017-05-28" = "omni_5min_2017-05-26_31_bz_speed_symh_cleaned.csv",
  "2023-04-24" = "omni_5min_2023-04-22_27_bz_speed_symh_cleaned.csv"
)

storm_paths <- file.path(data_dir, storm_files)
missing_files <- storm_paths[!file.exists(storm_paths)]

if (length(missing_files) > 0) {
  stop(
    "These storm files were not found:\n",
    paste(missing_files, collapse = "\n"),
    "\n\nChange data_dir or the filenames in storm_files."
  )
}

# ============================================================
# 3. GLOBAL LAG/BASIS SETTINGS
# ============================================================

if (max_lag < 1) {
  stop("max_lag must be at least 1.")
}

if (length(0:max_lag) <= bs_df) {
  stop("max_lag is too small for the chosen bs_df. Reduce bs_df or increase max_lag.")
}

lag_values <- 0:max_lag
lag_names <- paste0("lag", lag_values)
speed_lag_names <- paste0("speed_lag", lag_values)

# One common lag basis is used for every storm and every LOSO fold.
# Bz and speed use the same basis shape, but receive different fitted coefficients.
bs_basis <- bs(lag_values, df = bs_df)

bz_bs_names <- paste0("Bz_BS", seq_len(ncol(bs_basis)))
speed_bs_names <- paste0("Speed_BS", seq_len(ncol(bs_basis)))

cat("============================================================\n")
cat("Cross-storm LOSO validation\n")
cat("Maximum lag:", max_lag, "steps =",
    max_lag * data_interval_minutes, "minutes\n")
cat("Lagged predictor values per variable:", length(lag_values), "\n")
cat("Forecast horizon:", forecast_horizon * data_interval_minutes, "minutes\n")
cat("Evaluation refresh:", refresh_steps * data_interval_minutes, "minutes\n")
cat("============================================================\n\n")

# ============================================================
# 4. FUNCTION: PREPARE ONE STORM
# ============================================================
# Critical rule: interpolation, lags, and future targets are constructed
# WITHIN each storm. No lag is ever allowed to cross a storm boundary.

prepare_storm <- function(file_path, storm_id) {
  
  raw <- read.csv(file_path, stringsAsFactors = FALSE)
  
  required_columns <- c(
    "datetime_utc",
    "bz_gsm_nt",
    "flow_speed_kms",
    "sym_h_nt"
  )
  
  if (!all(required_columns %in% names(raw))) {
    stop(
      "File for storm ", storm_id,
      " does not contain all required columns: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  data <- raw %>%
    transmute(
      storm_id = storm_id,
      Time = as.POSIXct(datetime_utc, tz = "UTC"),
      Bz = as.numeric(bz_gsm_nt),
      Speed = as.numeric(flow_speed_kms),
      SYMH = as.numeric(sym_h_nt)
    ) %>%
    arrange(Time)
  
  # Protect against OMNI fill values in case any raw fill code remains.
  data$Bz[data$Bz > 9000] <- NA
  data$Speed[data$Speed > 90000] <- NA
  data$SYMH[abs(data$SYMH) > 9000] <- NA
  
  raw_missing_bz <- sum(is.na(data$Bz))
  raw_missing_speed <- sum(is.na(data$Speed))
  raw_missing_symh <- sum(is.na(data$SYMH))
  
  # Retrospective linear interpolation, matching the current project script.
  data$Bz_filled <- zoo::na.approx(
    data$Bz,
    x = as.numeric(data$Time),
    na.rm = FALSE
  )
  
  data$Speed_filled <- zoo::na.approx(
    data$Speed,
    x = as.numeric(data$Time),
    na.rm = FALSE
  )
  
  # Build Bz and speed lags separately inside this storm.
  for (i in lag_values) {
    data[[paste0("lag", i)]] <- dplyr::lag(data$Bz_filled, i)
    data[[paste0("speed_lag", i)]] <- dplyr::lag(data$Speed_filled, i)
  }
  
  # Direct future target.
  data$SYMH_future <- dplyr::lead(data$SYMH, forecast_horizon)
  data$Forecast_Time <- dplyr::lead(data$Time, forecast_horizon)
  
  # Keep only rows that can actually be used by all models.
  required_for_model <- c(
    "SYMH", "SYMH_future", "Forecast_Time",
    lag_names, speed_lag_names
  )
  
  data <- data[complete.cases(data[, required_for_model]), ]
  
  # B-spline compression of the lag histories.
  bz_lag_matrix <- as.matrix(data[, lag_names, drop = FALSE])
  speed_lag_matrix <- as.matrix(data[, speed_lag_names, drop = FALSE])
  
  bz_bs_predictors <- bz_lag_matrix %*% bs_basis
  speed_bs_predictors <- speed_lag_matrix %*% bs_basis
  
  colnames(bz_bs_predictors) <- bz_bs_names
  colnames(speed_bs_predictors) <- speed_bs_names
  
  data <- cbind(
    data,
    as.data.frame(bz_bs_predictors),
    as.data.frame(speed_bs_predictors)
  )
  
  cat(
    storm_id, ": rows =", nrow(raw),
    "| missing Bz =", raw_missing_bz,
    "| missing Speed =", raw_missing_speed,
    "| missing SYM-H =", raw_missing_symh,
    "| usable model rows =", nrow(data), "\n"
  )
  
  return(data)
}

# ============================================================
# 5. PREPARE ALL 10 STORMS
# ============================================================

storm_data <- vector("list", length(storm_files))
names(storm_data) <- names(storm_files)

for (storm_id in names(storm_files)) {
  storm_data[[storm_id]] <- prepare_storm(
    file_path = file.path(data_dir, storm_files[[storm_id]]),
    storm_id = storm_id
  )
}

cat("\nAll storms prepared separately.\n\n")

# ============================================================
# 6. METRIC FUNCTIONS
# ============================================================

rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2, na.rm = TRUE))
}

mae <- function(actual, predicted) {
  mean(abs(actual - predicted), na.rm = TRUE)
}

skill_vs_persistence <- function(model_rmse, persistence_rmse) {
  1 - model_rmse / persistence_rmse
}

# ============================================================
# 7. LEAVE-ONE-STORM-OUT CROSS-VALIDATION
# ============================================================
# For each fold:
#   - hold out exactly one entire storm
#   - train on the other nine storms
#   - predict only the held-out storm
#   - evaluate every refresh_steps rows on the held-out storm
#
# Models compared:
#   1. Persistence: future SYM-H = current SYM-H
#   2. SYM-H-only autoregression
#   3. Bz B-spline ARDL
#   4. Bz + Speed B-spline ARDL

all_fold_predictions <- list()
per_storm_metrics <- list()

for (test_id in names(storm_data)) {
  
  cat("------------------------------------------------------------\n")
  cat("Holding out storm:", test_id, "\n")
  
  train_ids <- setdiff(names(storm_data), test_id)
  
  train_data <- bind_rows(storm_data[train_ids])
  test_data <- storm_data[[test_id]]
  
  # Explicit safety check: held-out storm cannot appear in training data.
  if (test_id %in% unique(train_data$storm_id)) {
    stop("LOSO error: held-out storm appears in the training data.")
  }
  
  # ----------------------------------------------------------
  # MODEL A: SYM-H ONLY
  # ----------------------------------------------------------
  symh_only_model <- lm(
    SYMH_future ~ SYMH,
    data = train_data
  )
  
  # ----------------------------------------------------------
  # MODEL B: Bz-ONLY B-SPLINE ARDL
  # ----------------------------------------------------------
  bz_bs_formula <- reformulate(
    c("SYMH", bz_bs_names),
    response = "SYMH_future"
  )
  
  bz_bs_model <- lm(
    bz_bs_formula,
    data = train_data
  )
  
  # ----------------------------------------------------------
  # MODEL C: Bz + SPEED B-SPLINE ARDL
  # ----------------------------------------------------------
  bz_speed_bs_formula <- reformulate(
    c("SYMH", bz_bs_names, speed_bs_names),
    response = "SYMH_future"
  )
  
  bz_speed_bs_model <- lm(
    bz_speed_bs_formula,
    data = train_data
  )
  
  # ----------------------------------------------------------
  # PREDICT THE HELD-OUT STORM
  # ----------------------------------------------------------
  test_data$Predicted_Persistence <- test_data$SYMH
  
  test_data$Predicted_SYMH_Only <- predict(
    symh_only_model,
    newdata = test_data
  )
  
  test_data$Predicted_Bz_BS <- predict(
    bz_bs_model,
    newdata = test_data
  )
  
  test_data$Predicted_Bz_Speed_BS <- predict(
    bz_speed_bs_model,
    newdata = test_data
  )
  
  # Mimic the project's refresh design: evaluate every Nth eligible row.
  refresh_rows <- seq(1, nrow(test_data), by = refresh_steps)
  test_eval <- test_data[refresh_rows, ]
  
  # Keep only columns needed in final prediction output.
  fold_predictions <- test_eval %>%
    transmute(
      storm_id,
      Prediction_Made_At = Time,
      Forecast_Time,
      Current_SYMH = SYMH,
      Actual_Future_SYMH = SYMH_future,
      Predicted_Persistence,
      Predicted_SYMH_Only,
      Predicted_Bz_BS,
      Predicted_Bz_Speed_BS
    )
  
  all_fold_predictions[[test_id]] <- fold_predictions
  
  # ----------------------------------------------------------
  # METRICS FOR THIS HELD-OUT STORM
  # ----------------------------------------------------------
  persistence_rmse <- rmse(
    fold_predictions$Actual_Future_SYMH,
    fold_predictions$Predicted_Persistence
  )
  
  persistence_mae <- mae(
    fold_predictions$Actual_Future_SYMH,
    fold_predictions$Predicted_Persistence
  )
  
  symh_only_rmse <- rmse(
    fold_predictions$Actual_Future_SYMH,
    fold_predictions$Predicted_SYMH_Only
  )
  
  symh_only_mae <- mae(
    fold_predictions$Actual_Future_SYMH,
    fold_predictions$Predicted_SYMH_Only
  )
  
  bz_bs_rmse <- rmse(
    fold_predictions$Actual_Future_SYMH,
    fold_predictions$Predicted_Bz_BS
  )
  
  bz_bs_mae <- mae(
    fold_predictions$Actual_Future_SYMH,
    fold_predictions$Predicted_Bz_BS
  )
  
  bz_speed_bs_rmse <- rmse(
    fold_predictions$Actual_Future_SYMH,
    fold_predictions$Predicted_Bz_Speed_BS
  )
  
  bz_speed_bs_mae <- mae(
    fold_predictions$Actual_Future_SYMH,
    fold_predictions$Predicted_Bz_Speed_BS
  )
  
  per_storm_metrics[[test_id]] <- data.frame(
    storm_id = test_id,
    n_forecasts = nrow(fold_predictions),
    
    Persistence_RMSE = persistence_rmse,
    SYMH_Only_RMSE = symh_only_rmse,
    Bz_BS_RMSE = bz_bs_rmse,
    Bz_Speed_BS_RMSE = bz_speed_bs_rmse,
    
    Persistence_MAE = persistence_mae,
    SYMH_Only_MAE = symh_only_mae,
    Bz_BS_MAE = bz_bs_mae,
    Bz_Speed_BS_MAE = bz_speed_bs_mae,
    
    SYMH_Only_Skill = skill_vs_persistence(symh_only_rmse, persistence_rmse),
    Bz_BS_Skill = skill_vs_persistence(bz_bs_rmse, persistence_rmse),
    Bz_Speed_BS_Skill = skill_vs_persistence(bz_speed_bs_rmse, persistence_rmse),
    
    stringsAsFactors = FALSE
  )
  
  cat("Held-out forecasts:", nrow(fold_predictions), "\n")
  cat("Persistence RMSE:", round(persistence_rmse, 3), "\n")
  cat("SYM-H-only RMSE:", round(symh_only_rmse, 3), "\n")
  cat("Bz BS RMSE:", round(bz_bs_rmse, 3), "\n")
  cat("Bz + Speed BS RMSE:", round(bz_speed_bs_rmse, 3), "\n")
  cat("Bz + Speed skill vs persistence:",
      round(skill_vs_persistence(bz_speed_bs_rmse, persistence_rmse), 3), "\n")
}

# Combine all held-out predictions and per-storm metrics.
loso_predictions <- bind_rows(all_fold_predictions)
loso_metrics <- bind_rows(per_storm_metrics)

# ============================================================
# 8. OVERALL POOLED LOSO METRICS
# ============================================================
# These metrics pool every held-out forecast from all 10 folds.

pooled_persistence_rmse <- rmse(
  loso_predictions$Actual_Future_SYMH,
  loso_predictions$Predicted_Persistence
)

pooled_symh_only_rmse <- rmse(
  loso_predictions$Actual_Future_SYMH,
  loso_predictions$Predicted_SYMH_Only
)

pooled_bz_bs_rmse <- rmse(
  loso_predictions$Actual_Future_SYMH,
  loso_predictions$Predicted_Bz_BS
)

pooled_bz_speed_bs_rmse <- rmse(
  loso_predictions$Actual_Future_SYMH,
  loso_predictions$Predicted_Bz_Speed_BS
)

pooled_persistence_mae <- mae(
  loso_predictions$Actual_Future_SYMH,
  loso_predictions$Predicted_Persistence
)

pooled_symh_only_mae <- mae(
  loso_predictions$Actual_Future_SYMH,
  loso_predictions$Predicted_SYMH_Only
)

pooled_bz_bs_mae <- mae(
  loso_predictions$Actual_Future_SYMH,
  loso_predictions$Predicted_Bz_BS
)

pooled_bz_speed_bs_mae <- mae(
  loso_predictions$Actual_Future_SYMH,
  loso_predictions$Predicted_Bz_Speed_BS
)

pooled_metrics <- data.frame(
  Model = c(
    "Persistence",
    "SYM-H only",
    "Bz B-spline ARDL",
    "Bz + Speed B-spline ARDL"
  ),
  RMSE = c(
    pooled_persistence_rmse,
    pooled_symh_only_rmse,
    pooled_bz_bs_rmse,
    pooled_bz_speed_bs_rmse
  ),
  MAE = c(
    pooled_persistence_mae,
    pooled_symh_only_mae,
    pooled_bz_bs_mae,
    pooled_bz_speed_bs_mae
  ),
  Skill_vs_Persistence = c(
    0,
    skill_vs_persistence(pooled_symh_only_rmse, pooled_persistence_rmse),
    skill_vs_persistence(pooled_bz_bs_rmse, pooled_persistence_rmse),
    skill_vs_persistence(pooled_bz_speed_bs_rmse, pooled_persistence_rmse)
  ),
  stringsAsFactors = FALSE
)

# Also report the simple mean of the 10 storm-level RMSE/MAE values.
mean_storm_metrics <- data.frame(
  Model = c(
    "Persistence",
    "SYM-H only",
    "Bz B-spline ARDL",
    "Bz + Speed B-spline ARDL"
  ),
  Mean_Storm_RMSE = c(
    mean(loso_metrics$Persistence_RMSE),
    mean(loso_metrics$SYMH_Only_RMSE),
    mean(loso_metrics$Bz_BS_RMSE),
    mean(loso_metrics$Bz_Speed_BS_RMSE)
  ),
  Mean_Storm_MAE = c(
    mean(loso_metrics$Persistence_MAE),
    mean(loso_metrics$SYMH_Only_MAE),
    mean(loso_metrics$Bz_BS_MAE),
    mean(loso_metrics$Bz_Speed_BS_MAE)
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# 9. SAVE TABLES
# ============================================================

write.csv(
  loso_predictions,
  file.path(results_dir, "loso_all_predictions.csv"),
  row.names = FALSE
)

write.csv(
  loso_metrics,
  file.path(results_dir, "loso_per_storm_metrics.csv"),
  row.names = FALSE
)

write.csv(
  pooled_metrics,
  file.path(results_dir, "loso_pooled_metrics.csv"),
  row.names = FALSE
)

write.csv(
  mean_storm_metrics,
  file.path(results_dir, "loso_mean_storm_metrics.csv"),
  row.names = FALSE
)

# ============================================================
# 10. PLOT: ACTUAL VS PREDICTED THROUGH TIME
# ============================================================
# One PDF page per held-out storm.
# Every prediction shown here is out-of-sample for that storm.

pdf(
  file.path(figures_dir, "loso_actual_vs_predicted_by_storm.pdf"),
  width = 10,
  height = 6
)

for (storm_id in names(storm_files)) {
  
  plot_data <- loso_predictions %>%
    filter(.data$storm_id == storm_id) %>%
    arrange(Forecast_Time)
  
  y_limits <- range(
    plot_data$Actual_Future_SYMH,
    plot_data$Predicted_Persistence,
    plot_data$Predicted_Bz_Speed_BS,
    na.rm = TRUE
  )
  
  plot(
    plot_data$Forecast_Time,
    plot_data$Actual_Future_SYMH,
    type = "l",
    lwd = 2,
    col = "black",
    ylim = y_limits,
    xlab = "Forecast time",
    ylab = "SYM-H (nT)",
    main = paste("LOSO actual vs predicted - held-out storm", storm_id)
  )
  
  lines(
    plot_data$Forecast_Time,
    plot_data$Predicted_Bz_Speed_BS,
    type = "b",
    pch = 16,
    lwd = 1.5,
    col = "purple"
  )
  
  lines(
    plot_data$Forecast_Time,
    plot_data$Predicted_Persistence,
    lwd = 1.2,
    lty = 2,
    col = "gray40"
  )
  
  legend(
    "bottomright",
    legend = c(
      "Actual future SYM-H",
      "Predicted Bz + Speed BS",
      "Persistence"
    ),
    col = c("black", "purple", "gray40"),
    lty = c(1, 1, 2),
    pch = c(NA, 16, NA),
    lwd = c(2, 1.5, 1.2),
    cex = 0.8,
    bg = "white"
  )
}

dev.off()

# ============================================================
# 11. PLOT: PREDICTED VS ACTUAL SCATTER
# ============================================================
# A perfect forecast would lie on the 1:1 line.

png(
  file.path(figures_dir, "loso_predicted_vs_actual_scatter.png"),
  width = 1000,
  height = 850,
  res = 120
)

scatter_limits <- range(
  loso_predictions$Actual_Future_SYMH,
  loso_predictions$Predicted_Bz_Speed_BS,
  na.rm = TRUE
)

plot(
  loso_predictions$Actual_Future_SYMH,
  loso_predictions$Predicted_Bz_Speed_BS,
  pch = 16,
  cex = 0.7,
  col = "purple",
  xlim = scatter_limits,
  ylim = scatter_limits,
  xlab = "Actual future SYM-H (nT)",
  ylab = "Predicted future SYM-H (nT)",
  main = "LOSO predicted vs actual: Bz + Speed B-spline ARDL"
)

abline(0, 1, lty = 2, lwd = 2, col = "black")

dev.off()

# ============================================================
# 12. PRINT FINAL RESULTS
# ============================================================

cat("\n============================================================\n")
cat("PER-STORM LOSO METRICS\n")
cat("============================================================\n")
print(loso_metrics)

cat("\n============================================================\n")
cat("POOLED HELD-OUT METRICS\n")
cat("============================================================\n")
print(pooled_metrics)

cat("\n============================================================\n")
cat("MEAN OF THE 10 STORM-LEVEL METRICS\n")
cat("============================================================\n")
print(mean_storm_metrics)

cat("\nFiles written:\n")
cat(" -", file.path(results_dir, "loso_all_predictions.csv"), "\n")
cat(" -", file.path(results_dir, "loso_per_storm_metrics.csv"), "\n")
cat(" -", file.path(results_dir, "loso_pooled_metrics.csv"), "\n")
cat(" -", file.path(results_dir, "loso_mean_storm_metrics.csv"), "\n")
cat(" -", file.path(figures_dir, "loso_actual_vs_predicted_by_storm.pdf"), "\n")
cat(" -", file.path(figures_dir, "loso_predicted_vs_actual_scatter.png"), "\n")
cat("============================================================\n")
