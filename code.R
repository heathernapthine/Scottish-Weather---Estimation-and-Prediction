# Heather Napthine (s2065896)

#' test_stat
#'
#'  Function to to return test statistic, 
#'  T = |winter average − summer average|.
#'
#' @param values values to be permuted
#' @param group groups to assign to data
#' @param perm permutation to apply to data
#' @param statistic_fn function used for test statistic
#' 
#' @return A test statistic, that gives the absolute value of the difference
#'         in means between summer and winter data.  

test_stat <- function(values, group, perm = NULL, statistic_fn = mean){
  
  if(!is.null(perm)){
    values <- c(values[perm], values[-perm])
  }
  
  summeraverage = statistic_fn(values[group == "summer"]) 
  winteraverage = statistic_fn(values[group == "winter"])
  
  return(abs(summeraverage - winteraverage))
}

#' p_value_CI
#'
#'  Function to Construct a Monte Carlo permutation test
#'  and compute a p value and approximate confidence
#'  95% interval for the p value, for each station.
#'
#' @param StationID ID of station to compute scores for
#' @param alpha specified significance level
#' @param  N number of permutations for test
#' 
#' @return A vector containing the p value, computed using the 'test_stat'
#'         function, and its associated confidence interval.

p_value_CI <- function(StationID, alpha = 0.05, N = 10000){
  
  stationdata <- ghcnd %>%
    filter(ID == StationID)
  
  stationdata <- stationdata %>%
    filter(Element %in% c("PRCP"))
  
  nsummer <- length((stationdata %>%
                       filter(Summer == TRUE))$Summer)
  
  nwinter <- length((stationdata %>%
                       filter(Summer == FALSE))$Summer)
  
  null_data <- data.frame(
    values = stationdata$Value,  
    group = stationdata$Season)
  
  t_null <- numeric(N)
  for(i in seq_len(N)){
    if(i == 1){
      
      xindex <- 1:nsummer
      
    }else{
      
      xindex <- sample(seq_len(nsummer+nwinter), size = nsummer,
                       replace = FALSE)
      
    }
    
    t_null[i] <- test_stat(values = null_data$values, 
                           group = null_data$group, 
                           perm = xindex)
  }
  
  t_stat <- t_null[1]
  
  x = -1
  for (t in t_null){
    if(t >= t_stat){
      x = x + 1
    }
  }
  
  p <- x/N
  
  if (x > 0){
    
    CI = x/N - (sqrt((x*N - x**2)/(N**3)))*qnorm(c(1 - alpha / 2, alpha / 2))
    
  } else {
    
    CI = c(0, 1 - 0.025**(1/N))
  }
  
  monte <- sqrt(p*(1 - p)/N)
  
  Pvalue_CI_monte <- c(p, CI , monte)
  Pvalue_CI_monte
  
}

#' Model
#'
#'  Function to estimate models for the square root of the monthly
#'  averaged precipitation values in Scotland.
#'
#' @param data data for the model estimation
#' @param k specified frequency
#' 
#' @return An object of class "lm", containing the fitted model defined in 
#'         the Project2 document.

Model <- function(data, k){
  
  if (k == 0){
    fit <- lm(Value_sqrt_avg ~ 1 + Longitude + Latitude
              + Elevation + DecYear, data = data)
  } 
  
  if (k == 1){
    fit <- lm(Value_sqrt_avg ~ 1 + Longitude + Latitude
              + Elevation + DecYear + sin(2*pi*DecYear) + cos(2*pi*DecYear),
              data = data)
  }
  if (k == 2){
    fit <- lm(Value_sqrt_avg ~ 1 + Longitude + Latitude
              + Elevation + DecYear + sin(2*pi*DecYear) + cos(2*pi*DecYear) + sin(2*2*pi*DecYear)
              + cos(2*2*pi*DecYear),
              data = data)
  }
  if (k == 3){
    fit <- lm(Value_sqrt_avg ~ 1 + Longitude + Latitude
              + Elevation + DecYear + sin(2*pi*DecYear) + cos(2*pi*DecYear) + sin(2*2*pi*DecYear)
              + cos(2*2*pi*DecYear) + sin(3*2*pi*DecYear) + cos(3*2*pi*DecYear),
              data = data)
  }
  if (k == 4){
    fit <- lm(Value_sqrt_avg ~ 1 + Longitude + Latitude
              + Elevation + DecYear + sin(2*pi*DecYear) + cos(2*pi*DecYear) + sin(2*2*pi*DecYear)
              + cos(2*2*pi*DecYear) + sin(3*2*pi*DecYear) + cos(3*2*pi*DecYear)
              + sin(4*2*pi*DecYear) + cos(4*2*pi*DecYear),
              data = data)
  }
  
  fit
  
}

#' Scores
#'
#'  Function to estimate models for the square root of the monthly
#'  averaged precipitation values in Scotland.
#'
#' @param data
#' @param kvalue k parameter for model
#' 
#' @return  A `data.frame` with columns StationID, Name, Month,
#'          Value, and Element, where the values give the prediction scores for
#'          each station, as well as the overall cross-validated average scores,
#'          for both Squared Error and Dawid-Sebastiani scores, for each of the
#'          12 months of the year.

Scores <- function(data, kvalue){
  
  ScoreDataFrame <- data.frame()
  
  average_score_se <- c()
  average_score_ds <- c()
  
  for (m in (1:12)){
    
    monthdata <- data %>%
      filter(Month == m)
    
    scores_se <- c()
    scores_ds <- c()
    
    grp = 0
    
    for (id in StationIDs) {
      
      grp = grp + 1
      
      bigdata <- data %>% filter(ID != id)
      
      fit <- Model(bigdata,kvalue)
      
      smalldata <- monthdata %>% filter(ID == id)
      
      colnames(smalldata) <- colnames(bigdata)
      
      pred <- predict(fit, newdata = smalldata, se.fit = TRUE)
      sd <- sqrt((pred$se.fit)**2 + (pred$residual.scale)**2)
      
      scores_se[grp] <- mean(proper_score(
        "se",
        obs = smalldata %>% pull("Value_sqrt_avg"),
        mean = (pred$fit)))
      
      scores_ds[grp] <- mean(proper_score(
        "ds",
        smalldata %>% pull("Value_sqrt_avg"),
        mean = (pred$fit), sd = sd))
      
      row1 <- cbind(id,Names[grp],m,scores_se[grp],"StationScoreSE")
      row2 <- cbind(id,Names[grp],m,scores_ds[grp],"StationScoreDS")
      colnames(row2) <- colnames(row1)
      ScoreDataFrame  <- rbind(ScoreDataFrame,row1,row2)
    }
    average_score_ds[m] <- mean(scores_ds)
    average_score_se[m] <- mean(scores_se)
    
  }
  
  grp = 0
  
  for (id in StationIDs) {
    
    grp = grp + 1
    
    for (m in (1:12)){
      
      row3 <- cbind(id,Names[grp],m,average_score_se[m],"AverageScoreSE")
      row4 <- cbind(id,Names[grp],m,average_score_ds[m],"AverageScoreDS")
      
      colnames(row3) <- colnames(row1)
      colnames(row4) <- colnames(row1)
      ScoreDataFrame  <- rbind(ScoreDataFrame,row3,row4)
      
    }
  }
  
  colnames(ScoreDataFrame) <- c("StationID", "Name", "Month", "Value",
                                "Element")
  
  ScoreDataFrame
  
}
