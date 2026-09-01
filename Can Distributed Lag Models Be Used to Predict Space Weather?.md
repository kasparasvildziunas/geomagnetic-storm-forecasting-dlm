# geomagnetic-storm-forecasting-dlm
SYM-H Distributed Lag Modelling

This repository contains the data, R code, and report for a study investigating the use of distributed lag models for forecasting the SYM-H geomagnetic index from solar-wind conditions.

Files

cross_storm_validation.R — main analysis and cross-validation script

SYM-H_DLM_report.pdf — full research report

data/ — cleaned 5-minute OMNIWeb storm datasets

Data

The data were obtained from NASA's OMNIWeb High Resolution OMNI database.

Ten geomagnetic storm periods are included. Each CSV contains:

SYM-H

Bz in GSM coordinates

solar-wind flow speed

UTC timestamp

The files were cleaned from the original OMNIWeb downloads by converting missing-data/fill values into missing observations and retaining the relevant variables in a consistent format.

Missing Bz and flow-speed observations are not pre-interpolated in the CSV files. Linear interpolation is carried out inside the R analysis script.

R code

cross_storm_validation.R performs leave-one-storm-out (LOSO) cross-validation across the ten storms.

For each validation fold, one complete storm is held out while the remaining nine storms are used for model fitting. The held-out storm is then used only for evaluation. Predictor histories and future SYM-H targets are constructed separately within each storm.

The script:

constructs lagged histories of Bz and solar-wind speed;

represents the lag effects using a B-spline basis;

fits a SYM-H-only model, Bz ARDL, and Bz + solar-wind-speed ARDL;

compares these models with a persistence forecast;

calculates RMSE, MAE, and RMSE skill relative to persistence;

reports both storm-specific and pooled LOSO results;

generates actual-versus-predicted and predicted-versus-actual figures.

The main modelling parameters can be changed near the beginning of the script:

max_lag <- 24
forecast_horizon <- 24
refresh_steps <- 24

The observations are recorded at 5-minute resolution, so a value of 24 corresponds to 120 minutes.

Running the analysis

Install the required R packages:

install.packages(c("dplyr", "splines", "zoo"))

Place the storm CSV files in the data/ directory and set:

data_dir <- "data"

Then run:

source("cross_storm_validation.R")

The script automatically creates results/ and figures/ directories containing the validation outputs.

Note

Linear interpolation is retrospective and can use observations on both sides of a missing value. The analysis should therefore be interpreted as retrospective model validation rather than a complete real-time forecasting system.
