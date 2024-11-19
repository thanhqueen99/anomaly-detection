# READ ME
Detecting Anomalies in Satellite Orbit Data - Data Science Research Project

## Overview

This research project uses machine learning and statistical models for anomaly detection, to predict satellite manoeuvres based on given orbital data. The methods used are XGBoost and ARIMA. XGBoost is trained to forecast future orbital elements and detect anomalous events based on the residuals between the predicted and observed values. ARIMA is used for time series analysis of orbital features. Precision-Recall metrics and F1-Scores are calculated to evaluate the models' performance across multiple satellites.

## File Overview

- **man.readme**: A file with notes about manoeuvre files format.

- **manoeuvres/**: Contains manoeuvre data for satellites, used to label the observed manoeuvres.
  
- **orbital_elements/**: 
  - **propagated_elements/**: Contains orbital elements that have been propagated.
  - **unpropagated_elements/**: Contains orbital elements data.

- **cryosat2.R**: R Script for running the anomaly detection on CryoSat 2 satellite.

- **haiyang2a.R**: R Script for running the anomaly detection on Haiyang 2A satellite.

- **saral.R**: R Script for running the anomaly detection on SARAL satellite.

- **sentinel3a.R**: R Script for running the anomaly detection on Sentinel 3A satellite.

- **topex.R**: R Script for running the anomaly detection on TOPEX satellite.

## Improvements from Progress Report
- ARIMA method implementation for anomaly detection
- Precision-Recall analysis with Precision-Recall curves
- Clearer plots with increased text size and easier-to-read colours

## How to Run the Code
Ensure you set the correct working directory as all file paths within the R Scripts are relative to this directory.

### Software

Have the following installed:

- **R** 
- **RStudio** 

### All Required R Packages

The following R packages are required. You can install them by running:

```R
install.packages(c("caret", "readr", "tidyr", "dplyr", "lubridate", "janitor", "ggplot2", "xgboost", "lightgbm", "corrplot", "catboost", "forecast", "purrr"))


