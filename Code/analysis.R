#Analysis

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


#' Perform sensitivity analysis on length of the moving average window
#'
#' @param start_date 
#' @param end_date 
#' @param MA_windows 
#'
#' @return
#' @export
#'
#' @examples
MA_sens_analysis <- function(start_date, end_date, 
                             MA_windows = seq(from = 2, to = 28, by = 2)){
  
  betas <- vector("list", length(MA_windows))
  for (i in 1:length(MA_windows)){
    param <- frequentist_betas(measure_data = measure_data, 
                               start_date = start_date,
                               end_date = end_date,
                               MA_window = MA_windows[i])
    betas[[i]] <- param
  }
  
  names(betas) <- MA_windows
  sensitivity <- list("0-17"=NULL, "18-49"=NULL, "50-64"=NULL, "65+"=NULL)
  for (i in 1:4){
    
    for (j in 1:length(MA_windows)){
      sensitivity[[i]][j] <- betas[[j]]$R0[i]
    }
  }
  
  sensitivity_se <- list("0-17"=NULL, "18-49"=NULL, "50-64"=NULL, "65+"=NULL)
  for (i in 1:4){
    
    for (j in 1:length(MA_windows)){
      sensitivity_se[[i]][j] <- betas[[j]]$R0_se[i]
    }
  }
  
  
  sensitivities_MA <- data.frame(sensitivity)
  sensitivities_MA_se <- data.frame(sensitivity_se)
  names(sensitivities_MA) <- c("0-17", "18-49", "50-64", "65+")
  names(sensitivities_MA_se) <- c("0-17", "18-49", "50-64", "65+")
  
  sensitivities_MA_se <- sensitivities_MA_se %>% mutate("window" = MA_windows) %>% 
    pivot_longer(!window, names_to = "Age", values_to = "R0_se") %>% select(R0_se)
  
  sensitivities_MA <- sensitivities_MA %>% mutate("window" = MA_windows) %>% 
    pivot_longer(!window, names_to = "Age", values_to = "R0") %>% mutate(sensitivities_MA_se)
  
  return(sensitivities_MA)
}


#################
#### Bootstrap for nonparametric measure of uncertainty
#' Perform bootstrap analysis on beta estimates to get nonparametric CIs
#'
#' @param start_date 
#' @param end_date 
#'
#' @return a dataframe of beta estimates with 
#' @export
#'
#' @examples
bootstrapBeta <- function(start_date, end_date){
  
  start_window <- c(as.Date(start_date)-7, as.Date(start_date)+7)
  start_range <- seq(start_window[1], start_window[2], "days")
  end_window <- c(as.Date(end_date)-7, as.Date(end_date)+7)
  end_range <-  seq(end_window[1], end_window[2], "days")
  #ma_window_range <-  c(3, 7, 14)
  
  nboot <- 1000
  start_boot <-  sample(start_range, size = nboot, replace = TRUE)
  end_boot <-  sample(end_range, size = nboot, replace = TRUE)
  #ma_boot <- sample(ma_window_range, size = nboot, replace = TRUE) 
  
  betas <- vector("list", nboot)
  for (i in 1:nboot){
    
    param <- frequentist_betas(measure_data = measure_data, 
                               start_date = start_boot[i],
                               end_date = end_boot[i],
                               MA_window = 7)
    betas[[i]] <- param
    
  }
  
  
  beta_boot <- list("0-17"=NULL, "18-49"=NULL, "50-64"=NULL, "65+"=NULL)
  R0_boot <- list("0-17"=NULL, "18-49"=NULL, "50-64"=NULL, "65+"=NULL)
  for (i in 1:4){
    
    for (j in 1:nboot){
      beta_boot[[i]][j] <- betas[[j]]$beta[i]
      R0_boot[[i]][j] <- betas[[j]]$R0[i]
    }
  }
  
  
  beta_boot <- data.frame(beta_boot)
  R0_boot <- data.frame(R0_boot)
  names(beta_boot) <- c("0-17", "18-49", "50-64", "65+")
  names(R0_boot) <- c("0-17", "18-49", "50-64", "65+")
  
  beta <- sapply(beta_boot, mean)
  R0 <- sapply(R0_boot, mean)
  beta_boot_se <- sapply(beta_boot, sd)
  R0_boot_se <- sapply(R0_boot, sd)
  
  beta_quant <- sapply(beta_boot, quantile, probs = c(0.025, 0.975))
  R0_quant <- sapply(R0_boot, quantile, probs = c(0.025, 0.975))
  
  beta_df <- data.frame(age_group = names(beta), beta = beta, 
                        beta_lower = beta_quant[1,], beta_upper = beta_quant[2,],
                        row.names = NULL)
  
  return(beta_df)
}