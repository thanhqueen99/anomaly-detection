# READ ME
Detecting Anomalies in Satellite Orbit Data - Data Science Research Project

## Overview

This research project uses machine learning to predict satellite manoeuvres based on given orbital data. The method used is XGBoost, trained to forecast future orbital elements and detect anomalous events based on the residuals between the predicted and observed values.

## File Overview

- **manoeuvres/**: Contains manoeuvre data for satellites, used to label the observed manoeuvres.
  
- **orbital_elements/**: 
  - **propagated_elements/**: Contains orbital elements that have been propagated.
  - **unpropagated_elements/**: Contains orbital elements data.

- **cryosat2.R**: R Script for running the anomaly detection on CryoSat 2 satellite.

- **haiyang2a.R**: R Script for running the anomaly detection on Haiyang 2A satellite.

- **man.readme**: A file with notes about manoeuvre files format.

- **saral.R**: R Script for running the anomaly detection on SARAL satellite.

- **sentinel3a.R**: R Script for running the anomaly detection on Sentinel 3A satellite.

- **topex.R**: R Script for running the anomaly detection on TOPEX satellite.

## How to Run the Code
Ensure you set the correct working directory as all file paths within the R Scripts are relative to this directory.

### Software

Have the following installed:

- **R** 
- **RStudio** 

### All Required R Packages

The following R packages are required. You can install them by running:

```R
install.packages(c("caret", "readr", "tidyr", "dplyr", "lubridate", "janitor", "ggplot2", "xgboost", "lightgbm", "corrplot", "catboost", "forecast"))


