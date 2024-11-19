# Load packages
pacman::p_load(caret, readr, tidyr, dplyr, lubridate, janitor, ggplot2, xgboost, 
               corrplot, forecast, purrr)

# Load unpropagated dataset into R
cryosat2 <- read_csv(
  '/Users/vuminh/Documents/satellite_data/orbital_elements/unpropagated_elements/unpropagated_elements_CryoSat-2.csv')

# Convert first column to timestamp 
cryosat2 <- cryosat2 %>%
  rename(timestamp = `...1`) %>%
  mutate(timestamp = ymd_hms(timestamp))
head(cryosat2)

# Summary statistics of features
sumstat <- cryosat2 %>%
  summarise(
    `Mean Eccentricity` = mean(eccentricity),
    `Std Dev Eccentricity` = sd(eccentricity),
    `Mean Argument of Perigee` = mean(`argument of perigee`),
    `Std Dev Argument of Perigee` = sd(`argument of perigee`),
    `Mean Inclination` = mean(inclination),
    `Std Dev Inclination` = sd(inclination),
    `Mean Mean Anomaly` = mean(`mean anomaly`),
    `Std Dev Mean Anomaly` = sd(`mean anomaly`),
    `Mean Brouwer Mean Motion` = mean(`Brouwer mean motion`),
    `Std Dev Brouwer Mean Motion` = sd(`Brouwer mean motion`),
    `Mean Right Ascension` = mean(`right ascension`),
    `Std Dev Right Ascension` = sd(`right ascension`)
  )
sumstat

# Function for Exploratory data analysis plots
eda <- function(df, ft) {
  # Line graphs
  lin <- ggplot(df, aes(x = timestamp, y = .data[[ft]])) +
    geom_line(col = 'blue') +
    xlab('Years') +
    ylab(ft) +
    ggtitle(paste('Line graph of', ft, 'over years for CryoSat 2')) +
    theme_bw() + 
    theme(plot.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 12))
  
  # Histograms
  hst <- ggplot(df, aes(x = .data[[ft]])) +
    geom_histogram(aes(y = after_stat(density)), fill = 'pink') +
    geom_density(col = 'blue') +
    xlab(ft) +
    ylab('Distribution frequency') +
    ggtitle(paste('Histogram of', ft, 'distribution for CryoSat 2')) +
    theme_bw() + 
    theme(plot.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 12))
  
  # Box plots
  bxp <- ggplot(df, aes(y = .data[[ft]])) +
    geom_boxplot(fill = 'pink') +
    ylab(ft) +
    ggtitle(paste('Box plot of', ft, 'for CryoSat 2')) +
    theme_bw() + 
    theme(plot.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 12))
  
  print(hst)
  print(lin)
  print(bxp)
}

# EDA plots for features
fts <- colnames(cryosat2)[2:7]
for (ft in fts) {
  eda(cryosat2, ft)
}

# Centre and scale Brouwer mean motion
cryosat2 <- cryosat2 %>%
  mutate(`Brouwer mean motion` = scale(`Brouwer mean motion`))
head(cryosat2$`Brouwer mean motion`)

# Initialise satellite data
satdt <- list()
mandt <- function(fle) {
  man <- readLines(fle)
  man <- strsplit(man, "\\s+")
  man <- do.call(rbind, man)
  man <- as.data.frame(man, stringsAsFactors = FALSE)
  return(man)
}

# Convert manoeuvre timestamps
mandate <- function(man) {
  man <- man %>%
    mutate(
      begindate = make_datetime(
        as.integer(V2), as.integer(V3), as.integer(V4), as.integer(V5)),
      enddate = make_datetime(
        as.integer(V6), as.integer(V7), as.integer(V8), as.integer(V9))
    )
  return(man)
}

# Observed manoeuvres
obsman <- function(dt, man) {
  man$beginman <- as.POSIXct(
    paste(man$V2, man$V3, man$V4, man$V5, sep = "-"), 
    format = "%Y-%j-%H-%M")
  
  man$endman <- as.POSIXct(
    paste(man$V6, man$V7, man$V8, man$V9, sep = "-"), 
    format = "%Y-%j-%H-%M")
  dt$mn <- "n"
  
  for (i in seq_len(nrow(man))) {
    st <- man$beginman[i]
    ed <- man$endman[i]
    flobs <- which(dt$timestamp > ed)[1]
    
    if (!is.na(flobs)) {
      dt$mn[flobs] <- "y"
    } else {
      dt$mn[nrow(dt)] <- "y"
    }
  }
  return(dt)
}

# Update function
upd <- function(dt, falpos) {
  man <- mandt(falpos)
  man <- mandate(man)
  dt <- obsman(dt, man)
  return(dt)
}

# POSIX timestamps
cryosat2$timestamp <- as.POSIXct(
  cryosat2$timestamp, 
  format="%Y-%m-%d %H:%M:%S")

# Update satellite with Manoeuvre data 
fp <- '/Users/vuminh/Documents/satellite_data/manoeuvres/cs2man.txt'
cryosat2 <- upd(cryosat2, fp)

satdt <- list(CryoSat_2 = list(unpr = cryosat2))

# Function to get Precision, Recall, F-score
prerecft <- function(basedt, resis, obsman) {
  # Sort unique residuals in ascending order
  thrhol <- sort(unique(resis))
  
  # Apply to each element
  prerec <- thrhol %>%
    purrr::map_dfr(~{
      # Assign to variable
      thrd <- .x
      # Get predicted manoeuvres
      basedt <- basedt %>%
        mutate(predman = ifelse(resis > thrd, "y", "n"))
      
      # Get True Positive, False Positive, False Negative
      results <- basedt %>%
        mutate(
          # TP (predicted y == observed y)
          trupos = ifelse(predman == "y" & 
                            (obsman == "y" | 
                               lag(obsman, 1, default = "n") == "y" | 
                               lag(obsman, 2, default = "n") == "y" | 
                               lag(obsman, 3, default = "n") == "y"), 1, 0),
          
          # FP (predicted y != observed)
          falpos = ifelse(predman == "y" & trupos == 0, 1, 0),
          
          # FN (predicted n != observed)
          falneg = ifelse(predman == "n" & obsman == "y", 1, 0)
        ) %>%
        
        # Total counts
        summarise(
          trupos = sum(trupos),
          falpos = sum(falpos),
          falneg = sum(falneg)
        )
      
      # Precision
      prec <- ifelse(
        results$trupos + results$falpos == 0, 0, 
        results$trupos / (results$trupos + results$falpos))
      
      # Recall
      rec <- ifelse(
        results$trupos + results$falneg == 0, 0, 
        results$trupos / (results$trupos + results$falneg))
      
      # F-Score
      fs <- ifelse(
        prec + rec == 0, 0, 2 * (prec * rec) / (prec + rec))
      
      tibble(threshold = thrd, precision = prec, recall = rec, fscore = fs)
    })
  
  return(prerec)
}

# Precision Recall plot
plotprerec <- function(prerec, model, satname, feat) {
  ggplot(prerec, aes(x = recall, y = precision)) +
    geom_line(color = "blue") +
    labs(
      title = "Precision-Recall of ARIMA Predicting Brouwer mean motion Residuals on CryoSat 2",
      x = "Recall",
      y = "Precision"
    ) +
    theme_classic() + 
    theme(plot.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 12))
}

# ARIMA function
arimaft <- function(cryosat2, feat, satdt, thrd = 0.0128) {
  # Initialise variables
  model <- "ARIMA"  
  satname <- cryosat2  
  feat <- feat
  
  # Sort data by timestamps
  df <- satdt[[cryosat2]]$unpr
  df <- df %>%
    dplyr::arrange(timestamp) %>%
    dplyr::rename(Feature = !!sym(feat)) 
  
  # Split into 80% training set, 20% testing set
  dfsplit <- split(df, seq_len(nrow(df)) <= round(0.8 * nrow(df)))
  # Training set
  trainset <- dfsplit$`TRUE`
  # Testing set
  testset <- dfsplit$`FALSE`
  
  # Box-Cox transformation, making sure values are positive
  minfeat <- min(trainset$Feature)
  shift <- ifelse(minfeat <= 0, abs(minfeat) + 1, 0)
  trainset$shiftfeat <- trainset$Feature + shift
  
  # Box-Cox lambda
  lbda <- BoxCox.lambda(trainset$shiftfeat)
  # Box-Cox transformation on training set
  trainset <- trainset %>%
    mutate(featfit = BoxCox(
      shiftfeat, lbda))
  
  # Fit ARIMA model
  arimamod <- trainset %>%
    pull(featfit) %>%
    auto.arima(seasonal = TRUE)
  
  # Inverse Box-Cox transformation
  frmod <- forecast(arimamod, 
                    h = nrow(testset))
  bxcx <- frmod$mean %>% 
    as.numeric() %>% 
    InvBoxCox(lbda)
  
  # Compute residuals
  resis <- diff(testset$Feature - bxcx)
  bxcx <- bxcx[-1]
  testset <- testset[-1, ]
  
  # Anomaly detection
  testset <- testset %>%
    mutate(predman = ifelse(
      abs(resis) > thrd, "y", "n"))
  
  # Observed manoeuvre timestamps
  mandt <- testset %>%
    filter(mn == "y") %>%
    pull(timestamp)
  
  # Predicted manoeuvre timestamps
  predmandt <- testset %>%
    filter(predman == "y") %>%
    pull(timestamp)
  
  # Put in y or n format
  testset <- testset %>%
    mutate(obsman = ifelse(
      timestamp %in% mandt, "y", "n"))
  
  # Compute Precision, Recall
  prerec <- prerecft(
    testset,
    abs(resis),
    obsman
  )
  
  # Plot residuals
  arimaplot <- ggplot(testset, 
                      aes(x = timestamp, y = abs(resis))) +
    # Residuals shown in red colour
    geom_line(color = "red") + 
    # Solid green line = observed manoeuvres
    geom_vline(aes(xintercept = timestamp), 
               data = tibble(timestamp = mandt),
               color = "green",
               linetype = "solid") +
    # Dashed blue line = predicted manoeuvres
    geom_vline(aes(xintercept = timestamp),
               data = tibble(timestamp = predmandt),
               color = "darkblue",
               linetype = "dashed") + labs(
                 subtitle = "Observed Manoeuvres: Green solid lines, Predicted Manoeuvres: Blue dashed lines",
                 title = glue::glue("ARIMA: {cryosat2} Manoeuvres Observed versus Predicted on {feat}"),
                 x = "Time",
                 y = "Residuals"
               ) +
    theme_classic() + 
    theme(plot.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 12))
  
  # Plot Precision Recall curve
  prerecplot <- plotprerec(prerec, model, satname, feat)
  
  print(arimaplot)
  print(prerecplot)
  
  # Return results
  return(list(forecast = bxcx, model = arimamod, prerec = prerec))
}

# Anomaly detection with ARIMA 
aricall <- arimaft(
  cryosat2 = "CryoSat_2",
  feat = "Brouwer mean motion",
  satdt = satdt,
  thrd = 0.0128
)

# ARIMA Precision, Recall, F-score
arimametrics <- function(satdt, feat, thrd = 0.0128) {
  # Anomaly detection with ARIMA 
  aricall <- arimaft(
    cryosat2 = "CryoSat_2",
    feat = feat,
    satdt = satdt,
    thrd = thrd
  )
  
  # Precision, Recall, F score 
  finalprerec <- aricall$prerec
  
  # Select best threshold with F score
  bestthr <- finalprerec %>%
    filter(fscore == max(fscore)) %>%
    head(1)
  
  # Return best threshold
  return(list(
    `ARIMA Precision Recall Metrics` = finalprerec,
    `ARIMA Best threshold` = bestthr
  ))
}

# ARIMA Metrics: Precision, Recall, F-Score
arimametrics(satdt, "Brouwer mean motion", thrd = 0.0128)

# XGBoost
# XGBoost Function for Anomaly Detection 
predsat <- function(df, feat, upth, lwth) {
  # Number of lag features
  lagft <- 3
  
  # Sort data by timestamps
  df <- df %>%
    arrange(timestamp) %>%
    rename(Feature = all_of(feat))
  
  # Reference date
  ref <- as.Date("1999-06-02")
  # Convert timestamps to numeric days
  dat <- as.numeric(as.Date(df$timestamp) - ref)
  # Get rid of NA values
  dat <- dat[!is.na(dat)]
  # Get median day
  medday <- median(dat, na.rm = TRUE)
  # Convert median day to date
  meddat <- ref + medday
  
  # Normalise features
  df$Feature <- scale(df$Feature, center = TRUE, scale = TRUE)
  
  # Select timestamp and feature
  res <- df %>% select(timestamp, Feature)
  lgres <- res
  
  # Get lag features
  for (lag in 1:lagft) {
    lgres <- lgres %>%
      mutate(!!paste0("lag_", lag) := dplyr::lag(Feature, n = lag))
  }
  
  # Filter out rows
  lgres <- lgres %>%
    filter(row_number() > lagft)
  res <- res %>%
    filter(row_number() > lagft)
  
  # Match timestamps to data
  lgres <- lgres %>%
    mutate(timestamp = df$timestamp[lagft + 1:n()])
  res <- res %>%
    mutate(timestamp = df$timestamp[lagft + 1:n()])
  
  meddat <- as.Date(meddat)
  
  # Split data
  lgrestr <- lgres[lgres$timestamp < meddat, -which(
    names(lgres) == "timestamp")]
  lgrests <- lgres[lgres$timestamp >= meddat, -which(
    names(lgres) == "timestamp")]
  
  restr <- res[res$timestamp < meddat, -which(
    names(res) == "timestamp")]
  rests <- res[res$timestamp >= meddat, -which(
    names(res) == "timestamp")]
  
  # Convert into DMatrix format
  mattr <- xgb.DMatrix(data = as.matrix(
    lgrestr[, grep("^lag", names(lgrestr))]), label = restr$Feature)
  matts <- xgb.DMatrix(data = as.matrix(
    lgrests[, grep("^lag", names(lgrests))]), label = rests$Feature)
  
  # Define hyperparameters
  params <- list(
    objective = "reg:squarederror", # Regression
    eval_metric = "mae",  # Evaluation 
    max_depth = 2,      # Max trees depth  
    eta = 0.1           # Learning rate
  )
  
  # Train model
  xgbmod <- xgb.train(
    params = params,
    data = mattr, # Input Training data
    nrounds = 100,   # Number of rounds        
    watchlist = list(train = mattr, test = matts),
    early_stopping_rounds = 10,  # Early stopping
    verbose = TRUE    # Print
  )
  
  # Check feature importance
  impmat <- xgb.importance(model = xgbmod)
  print(xgb.plot.importance(impmat))
  
  # Predict for testing set
  preds <- predict(xgbmod, newdata = matts)
  # Get residuals
  resi <- preds - rests$Feature
  
  respl <- data.frame(
    timestamp = df$timestamp[df$timestamp >= meddat][1:length(preds)],
    resi = resi
  )
  
  # Anomaly detection
  respl <- respl %>%
    mutate(
      # Residuals outside of thresholds = anomalous events
      manpred = ifelse(resi > upth | resi < lwth, "y", "n"),
      # Convert into catagorical variable
      manpred = factor(manpred)
    )
  
  # Get manoeuvre timestamps
  obsmanti <- df %>%
    filter(mn == "y", timestamp >= meddat) %>%
    select(timestamp) %>%
    pull()
  
  # Get predicted timestamps
  predmanti <- respl %>%
    filter(manpred == "y") %>%
    select(timestamp) %>%
    pull()
  
  # Manoeuvres predicted vs observed
  obsmanti_df <- data.frame(timestamp = obsmanti)
  predmanti_df <- data.frame(timestamp = predmanti)
  
  # Plot residuals + observed and predicted manoeuvres
  pred_obs <- ggplot(respl, aes(x = timestamp, y = abs(resi))) +
    # Residuals shown in colour red
    geom_line(color = "red") +
    # Observed manoeuvres = solid green line
    geom_vline(aes(xintercept = timestamp), 
               data = obsmanti_df, 
               color = "green", 
               linetype = "solid") +
    # Predicted manoeuvres = dashed blue line
    geom_vline(aes(xintercept = timestamp), 
               data = predmanti_df, 
               color = "darkblue", 
               linetype = "dashed") +
    labs(subtitle = "Observed Manoeuvres: Green solid lines, Predicted Manoeuvres: Blue dashed lines",
         title = "XGBoost: CryoSat 2 Manoeuvres Observed versus Predicted on Brouwer mean motion",
         x = "Time",
         y = "Residuals"
    ) +
    theme_classic() + 
    theme(plot.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 12))
  print(pred_obs)
  
  return(respl)
}

# Precision, Recall, F-score function
prerecft <- function(basedt, resis, obsman) {
  # Sort residuals in ascending order
  thrhol <- sort(unique(resis))
  
  # Precision, Recall, F-Score for each threshold
  prerec <- thrhol %>%
    purrr::map_dfr(~{
      # Assign to variable
      thrd <- .x
      # Get predicted manoeuvres
      basedt <- basedt %>%
        mutate(predman = ifelse(resis > thrd, "y", "n"))
      
      # Get True Positive, False Positive, False Negative
      results <- basedt %>%
        mutate(
          # True Positive (predicted y == observed y)
          trupos = ifelse(predman == "y" & 
                            (obsman == "y" | 
                               lag(obsman, 1, default = "n") == "y" | 
                               lag(obsman, 2, default = "n") == "y" | 
                               lag(obsman, 3, default = "n") == "y"), 1, 0),
          # False Positive (predicted y != observed)
          falpos = ifelse(predman == "y" & trupos == 0, 1, 0),
          # False Negative (predicted n != observed y)
          falneg = ifelse(predman == "n" & obsman == "y", 1, 0)
        ) %>%
        # TP, FP, FN counts 
        summarise(
          trupos = sum(trupos),
          falpos = sum(falpos),
          falneg = sum(falneg)
        )
      
      # Precision
      prec <- ifelse(
        results$trupos + results$falpos == 0, 0, 
        results$trupos / (results$trupos + results$falpos))
      
      # Recall
      rec <- ifelse(
        results$trupos + results$falneg == 0, 0, 
        results$trupos / (results$trupos + results$falneg))
      
      # F-Score
      fs <- ifelse(
        prec + rec == 0, 0, 2 * (prec * rec) / (prec + rec))
      
      # Precision Recall Metrics
      tibble(threshold = thrd, precision = prec, recall = rec, fscore = fs)
    })
  
  return(prerec)
}

# Plot Precision Recall Curve
plotprerec <- function(prerec) {
  ggplot(prerec, 
         aes(x = recall, y = precision)) +
    geom_line(color = "blue") +
    labs(
      title = "Precision-Recall of XGBoost Predicting Brouwer mean motion Residuals on CryoSat 2", 
      x = "Recall", y = "Precision") +
    theme_classic() + 
    theme(plot.title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 12))
}

resplot <- function(satdt, feat, upth, lwth) {
  # Get data for satellite
  dt <- satdt[["CryoSat_2"]]$unpr
  
  # Run XGBoost function
  respl <- predsat(dt, feat, upth, lwth)
  
  # Mark y or n
  obsman <- ifelse(respl$manpred == "y", "y", "n")
  
  # Get Precision, Recall, F-score 
  finalprerec <- prerecft(respl, resis = respl$resi, obsman = obsman)
  
  # Plot Precision Recall curve
  print(plotprerec(finalprerec))
  
  # Best threshold 
  xgbestthr <- finalprerec %>%
    filter(fscore == max(fscore)) %>%
    head(1)  
  
  print("XGBoost Best Threshold:")
  print(xgbestthr)
}

satdt <- list()  
satdt[["CryoSat_2"]] <- list(unpr = cryosat2)

# Run function
resplot(satdt, "Brouwer mean motion", upth = 1, lwth = -1)
