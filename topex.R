pacman::p_load(caret, readr, tidyr, dplyr, lubridate, janitor, ggplot2, xgboost, lightgbm, corrplot, catboost, forecast)

# Load unpropagated TOPEX dataset into R
topex <- read_csv('/Users/vuminh/Documents/satellite_data/orbital_elements/unpropagated_elements/unpropagated_elements_TOPEX.csv')

# Convert first column to timestamp 
topex <- topex %>%
  rename(timestamp = `...1`) %>%
  mutate(timestamp = ymd_hms(timestamp))
topex

# Centre and scale Brouwer mean motion
topex <- topex %>%
  mutate(`Brouwer mean motion` = scale(`Brouwer mean motion`))
topex

# Summary statistics of features
sumstat <- topex %>%
  summarise(
    eccmean = mean(eccentricity),
    eccsd = sd(eccentricity),
    aopmean = mean(`argument of perigee`),
    aopsd = sd(`argument of perigee`),
    incmean = mean(inclination),
    incsd = sd(inclination),
    mamean = mean(`mean anomaly`),
    masd = sd(`mean anomaly`),
    bmmmean = mean(`Brouwer mean motion`),
    bmmsd = sd(`Brouwer mean motion`),
    right_ascension_mean = mean(`right ascension`),
    right_ascension_sd = sd(`right ascension`)
  )
sumstat

fts <- colnames(topex)[2:7]

# Function for Exploratory data analysis plots
eda <- function(df, ft) {
  lin <- ggplot(df, aes(x = timestamp, y = .data[[ft]])) +
    geom_line(col = 'blue') +
    xlab('Years') +
    ylab(ft) +
    ggtitle(paste('Line graph of', ft, 'over years for TOPEX')) +
    theme_bw()
  
  hst <- ggplot(df, aes(x = .data[[ft]])) +
    geom_histogram(aes(y = after_stat(density)), fill = 'pink') +
    geom_density(col = 'blue') +
    xlab(ft) +
    ylab('Distribution frequency') +
    ggtitle(paste('Histogram of', ft, 'distribution for TOPEX')) +
    theme_bw()
  
  bxp <- ggplot(df, aes(y = .data[[ft]])) +
    geom_boxplot(fill = 'pink') +
    ylab(ft) +
    ggtitle(paste('Box plot of', ft, 'for TOPEX')) +
    theme_bw()
  
  print(hst)
  print(lin)
  print(bxp)
}

# EDA plots for TOPEX features
for (ft in fts) {
  eda(topex, ft)
}

# Initialise satellite data
satdt <- list()
mandt <- function(fle) {
  man <- readLines(fle)
  man <- strsplit(man, "\\s+")
  man <- do.call(rbind, man)
  man <- as.data.frame(man, stringsAsFactors = FALSE)
  return(man)
}

# Convert dates of manoeuvres
mandate <- function(man) {
  man <- man %>%
    mutate(
      begindate = make_datetime(as.integer(V2), as.integer(V3), as.integer(V4), as.integer(V5)),
      enddate = make_datetime(as.integer(V6), as.integer(V7), as.integer(V8), as.integer(V9))
    )
  return(man)
}

# Observed manoeuvres
obsman <- function(dt, man) {
  man$beginman <- as.POSIXct(paste(man$V2, man$V3, man$V4, man$V5, sep = "-"), format = "%Y-%j-%H-%M")
  man$endman <- as.POSIXct(paste(man$V6, man$V7, man$V8, man$V9, sep = "-"), format = "%Y-%j-%H-%M")
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
upd <- function(dt, fp) {
  man <- mandt(fp)
  man <- mandate(man)
  dt <- obsman(dt, man)
  return(dt)
}

# POSIX timestamps
topex$timestamp <- as.POSIXct(topex$timestamp, format="%Y-%m-%d %H:%M:%S")

# Manoeuvres for TOPEX
fp <- '/Users/vuminh/Documents/satellite_data/manoeuvres/topman.txt'
topex <- upd(topex, fp)

predsat <- function(df, feat, upth, lwth) {
  lagft <- 3
  df <- df %>%
    arrange(timestamp) %>%
    rename(Feature = all_of(feat))
  
  ref <- as.Date("1999-06-02")
  dat <- as.numeric(as.Date(df$timestamp) - ref)
  dat <- dat[!is.na(dat)]  
  medday <- median(dat, na.rm = TRUE)
  meddat <- ref + medday
  
  df$Feature <- scale(df$Feature, center = TRUE, scale = TRUE)
  
  # Create data frame with timestamps and normalised feature
  res <- df %>% select(timestamp, Feature)
  lgres <- res
  for (lag in 1:lagft) {
    lgres <- lgres %>%
      mutate(!!paste0("lag_", lag) := dplyr::lag(Feature, n = lag))
  }
  
  lgres <- lgres %>%
    filter(row_number() > lagft)
  res <- res %>%
    filter(row_number() > lagft)
  
  lgres <- lgres %>%
    mutate(timestamp = df$timestamp[lagft + 1:n()])
  res <- res %>%
    mutate(timestamp = df$timestamp[lagft + 1:n()])
  
  meddat <- as.Date(meddat)
  
  lgrestr <- lgres[lgres$timestamp < meddat, -which(names(lgres) == "timestamp")]
  lgrests <- lgres[lgres$timestamp >= meddat, -which(names(lgres) == "timestamp")]
  
  restr <- res[res$timestamp < meddat, -which(names(res) == "timestamp")]
  rests <- res[res$timestamp >= meddat, -which(names(res) == "timestamp")]
  
  mattr <- xgb.DMatrix(data = as.matrix(lgrestr[, grep("^lag", names(lgrestr))]), label = restr$Feature)
  matts <- xgb.DMatrix(data = as.matrix(lgrests[, grep("^lag", names(lgrests))]), label = rests$Feature)
  
  params <- list(
    objective = "reg:squarederror",
    eval_metric = "mae",   
    max_depth = 2,          
    eta = 0.1
  )
  
  xgbmod <- xgb.train(
    params = params,
    data = mattr,
    nrounds = 100,           
    watchlist = list(train = mattr, test = matts),
    early_stopping_rounds = 10,  
    verbose = TRUE
  )
  
  # Check feature importance
  impmat <- xgb.importance(model = xgbmod)
  print(xgb.plot.importance(impmat))
  
  preds <- predict(xgbmod, newdata = matts)
  resi <- preds - rests$Feature
  
  respl <- data.frame(
    timestamp = df$timestamp[df$timestamp >= meddat][1:length(preds)],
    resi = resi
  )
  
  respl <- respl %>%
    mutate(
      manpred = ifelse(resi > upth | resi < lwth, "y", "n"),
      manpred = factor(manpred)
    )
  
  obsmanti <- df %>%
    filter(mn == "y", timestamp >= meddat) %>%
    select(timestamp) %>%
    pull()
  
  predmanti <- respl %>%
    filter(manpred == "y") %>%
    select(timestamp) %>%
    pull()
  
  # Manoeuvres predicted versus observed
  # Get observed and predicted manoeuvre timestamps
  obsmanti_df <- data.frame(timestamp = obsmanti)
  predmanti_df <- data.frame(timestamp = predmanti)
  
  # Plot
  pred_obs <- ggplot(respl, aes(x = timestamp, y = resi)) +
    geom_line(color = "sienna4") +
    geom_vline(aes(xintercept = timestamp), data = obsmanti_df, color = "darkorange", linetype = "dashed") +
    geom_vline(aes(xintercept = timestamp), data = predmanti_df, color = "blue", linetype = "dotted") +
    labs(
      title = "TOPEX Manoeuvres predicted versus observed",
      x = "Years",
      y = "Residuals"
    ) +
    theme_minimal()
  print(pred_obs)
}

satdt[["TOPEX"]] <- list(unpr = topex)

resplot <- function(satdt, feat, upth, lwth) {
  dt <- satdt[["TOPEX"]]$unpr
  predsat(dt, feat, upth, lwth)
}
resplot(satdt, "Brouwer mean motion", upth = 1, lwth = -1)