# Forecasting

# Arguments that go into the seir() function:
# population: Total number of people, living and dead.
# init_infectious, init_exposed, etc: People infectious, exposed, etc. at the start of the model.
# beta_f: Average people infected by a given infectious person per day. Can act as a function of time.
# incubation_days: Time from when a person is exposed to when they're infectious.
# recovery_days: Time it takes someone to either recover or die from the disease.
# ifr: Infection fatality rate- percentage of infections that are fatal. This will vary by age group.
# birthrate, deathrate: Additional births and deaths outside of COVID, measured per person per day. Default values are based on US figures.
# vacc_f: Rate of vaccination per day- is a function of time to reflect changes in production rates.
# days: Amount of days we want the model to simulate.

#library(dplyr)
seir= function(population,init_infectious,
               init_exposed=init_infectious/10,init_recovered=0,init_dead=0,
               beta_f,incubation_days=6,recovery_days=10,
               ifr,vacc_f,days=100){
  # Initial values for S, E, I, R, and D populations
  susc = population-init_infectious-init_exposed-init_recovered-init_dead
  exp = init_exposed
  inf = init_infectious
  rec = init_recovered
  dead = init_dead
  living = population-init_dead
  # Initialize our vectors that will be storing SEIRD populations over time.
  Svec = susc
  Evec = exp
  Ivec = inf
  Rvec = rec
  Dvec = dead
  Lvec = living
  for(i in 1:days-1){
    # Calculate the number of infections- people moving from S to E.
    infections = round(beta_f(i)*susc*inf/population)
    # Calculate how many exposed people become infectious
    incubations = round(exp/incubation_days)
    # Calculate how many infectious people recover or die- this will vary by IFR, which varies by age.
    deaths = round(inf/recovery_days*ifr)
    recoveries = round((inf/recovery_days)*(1-ifr))
    # Calculate how many vaccinations are done upon the susceptible.
    # We assume vaccinations are done equally upon all living groups and have no effect on the non-susceptible.
    vaccs = round(vacc_f(i)*susc/living)
    # Update all groups:
    susc = susc-infections-vaccs
    exp = exp+infections-incubations
    inf = inf+incubations-recoveries-deaths
    rec = rec+recoveries+vaccs
    dead = dead+deaths
    living = living-deaths
    # Update all vectors:
    Svec = c(Svec,susc)
    Evec = c(Evec,exp)
    Ivec = c(Ivec,inf)
    Rvec = c(Rvec,rec)
    Dvec = c(Dvec,dead)
    Lvec = c(Lvec,living)
  }
  # Output a data frame with all our vectors
  return(data.frame(Day=0:days,Susceptible=Svec,Exposed=Evec,Infectious=Ivec,
             Recovered=Rvec,Dead=Dvec,Living=Lvec))
}


