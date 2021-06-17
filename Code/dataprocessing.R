## Data Processing

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

#If covidAgeData is not installed, run this step 
# remotes::install_github("eshom/covid-age-data")
require(covidAgeData) #Obtaining Covid-19 case and death data

################################
#' Download COVID-19 data from COVerAGE-DB
#'
#'
#'
#' @return a processed dataframe, filtered on california
#' 
#'
#' @examples
getData <- function(){
  coverage <- download_covid("inputDB") #Full data
  ca_data <- subset_covid(coverage, Country = "USA", #Data for CA only
                          Region = "California")
  measure_data <- ca_data %>%
    filter(Measure %in% c("Deaths", "Cases") & Sex == "b" &
             Age %in% c("0", "18", "50", "65")) %>% #These are lower bounds
    mutate(age_group = if_else(Age == "0", "0-17", #Desired age intervals
                               if_else(Age == "18", "18-49",
                                       if_else(Age == "50", "50-64",
                                               if_else(Age == "65", "65+", NA_character_)))),
           date = as_date(Date, #Converting string date to object
                          format = "%d.%m.%Y")) %>% 
    filter(age_group %ni% NA_character_) %>% #Removing any potential problem obs
    select(date, Measure, Value, age_group) %>% #Removing irrelevant identifiers
    pivot_wider(names_from = "Measure", values_from = "Value")
  
  #Change variable names to lowercase for greater convenience
  names(measure_data) <- tolower(names(measure_data))
  
  
  return(measure_data)
}


#######################
