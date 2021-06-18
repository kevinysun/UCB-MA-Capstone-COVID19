# UC Berkeley MA Capstone Project: California COVID-19 Epidemiology
 
## Description 
This repository contains the code and data for a semester-long final project in UC Berkeley's MA Capstone class, done with Mallika Snyder and Ian Shen. For this project we were tasked with finding data and formulating and answering a question related to the epidmiological behavior of COVID-19. We worked with the [COVerAGE-DB dataset](https://osf.io/mpwjq/) from the Max Planck Institute and investigated how different vaccination prioritization schemes would affect the outcome of the virus based on SEIR compartmental modeling. The project had an interesting characteristic in that it developed at the same time that vaccinations were being released and distribution was occuring in California. This project involved literature review in order to gain familiarity with epidemiological practices, data processing and filtering, and statistical analysis in terms of fitting parameters for the SEIR model and nonparametric methods for estimation of the variance.  

For the complete methodology and analysis, refer to [the individual final report](Final_Writeup.pdf) submitted by Kevin Sun. 

## Structure
### Code
+ **dataprocessing.R** contains functions used to import and clean the full data. Note that it uses a package that downloads the COVerAGE data
+ **fitting.R** contains functions for fitting regressions for each age compartment
+ **analysis.R** contains functions for performing sensitivity analysis and bootstrap sampling on fitted models
+ **forecasting.R** contains functions for predicting SEIR trajectories based on various vaccination distribution schemes (represented by functions) as well as extracting the actual SEIR curve from more recent data
+ **execution.R** is a script used to generate outputs using the functions of the above R scripts

