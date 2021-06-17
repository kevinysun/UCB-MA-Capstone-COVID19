#Fitting

require(zoo) #Data wrangling
require(tidyverse) #Data wrangling
require(lubridate) #For date variables
require(naniar) # analysis of missing values
require(broom) #Tidy regression output
require(kableExtra) #For tables
library(reshape2)
require(purrr)
require(msm)
require(boot)


#' Frequentist fitting of COVID-19 data to estimate beta
#'
#' @param measure_data dataframe of COVID case counts
#' @param start_date starting date for fitting
#' @param end_date ending date for fitting
#' @param MA_window how many days to include in moving average calculation
#'
#' @return parameters by age
#' @export
#'
#' @examples
frequentist_betas <- function(measure_data, # dataframe
                              start_date,   # start date for frequentist fitting
                              end_date,     # end date for frequentist fitting
                              MA_window     # how many days for rolling average
){  
  #Calculate daily lagged data, then right-aligned rolling averages
  measure_data_final <- measure_data %>%
    group_by(age_group) %>% 
    mutate(lag = lag(date), 
           deaths_daily = ifelse(lag == lag(date), deaths - lag(deaths), NA),
           cases_daily = ifelse(lag == lag(date), cases - lag(cases), NA)) %>%
    arrange(date) %>%
    ungroup() %>%
    group_by(age_group) %>%
    mutate(cases_daily_ma = rollmeanr(cases_daily,MA_window, fill = NA),
           deaths_daily_ma = rollmeanr(deaths_daily,MA_window, fill = NA)) %>%
    arrange(date)
  #Obtain models by age group
  models_by_age <- measure_data_final %>%
    filter(date >= start_date & date <= end_date) %>% 
    mutate(time = row_number()) %>%
    group_by(age_group) %>%
    do(model = lm(log(cases_daily_ma) ~ time, data = .))
  
  # calculate the variance
  summary_by_age <- lapply(models_by_age$model, summary)
  cov_by_age <- lapply(models_by_age$model, vcov)
  
  
  #Tidy this output to obtain coefficients
  coefficients <- unlist(lapply(1:nrow(models_by_age), function(x) {
    val <- tidy(models_by_age$model[[x]]) %>%
      filter(term == "time") %>%
      pull(estimate)
    
    return(val)
  }
  ))
  
  
  coeff_sesq <- cov_by_age %>% map(~.x[2, 2]) %>% unlist()
  
  
  #delta method for se of R_0, beta
  R0_sesq <- rep(0,nrow(models_by_age))
  for(i in 1:nrow(models_by_age)){
    R0_sesq[i] <- deltamethod(~(1+10*x1)*(1+6*x1), coefficients[i], 
                              coeff_sesq[i])
  }
  beta0_sesq <- R0_sesq/((10+6)^2)
  #Calculate R0 and beta by age
  params <- tibble(age_group = models_by_age$age_group,
                   a = coefficients,
                   R0 = (1+10*a)*(1+6*a),
                   beta = R0/(10 + 6),
                   R0_se = sqrt(R0_sesq),
                   beta_se = sqrt(beta0_sesq))
  
  params <- params %>%
    rename(Age = "age_group")
  return(params)
}


##################

