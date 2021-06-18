# executable

source("dataprocessing.R")
source("fitting.R")
source("analysis.R")
source("forecasting.R")


# load and process data
measure_data <- getData()


# sensitivity analysis
sensitivities_MA <- MA_sens_analysis(start_date = "2020-11-01",
                                     end_date = "2020-12-15")


point_est <- MA_sens_analysis(start_date = "2020-11-01",
                              end_date = "2020-12-15", MA_windows = 7)


sensitivities_MA %>% mutate(beta = R0/16, beta_se = R0_se/16) %>% 
  ggplot() + aes(x=window, y=R0/16, color=Age)+
  geom_line()+geom_point()+  
  geom_ribbon(aes(ymin=beta-beta_se, ymax=beta+beta_se), linetype=2, alpha=0.1)+
  labs(x = "Days for Moving Average Calculation", 
       y = "Frequentist Beta estimate", color = "Age Group",
       title = "Frequentist Beta sensitivity to MA window length (11-1-20 to 12-15-20)")

# bootstrapping

beta_df <- bootstrapBeta(start_date = "2020-11-01",
                         end_date = "2020-12-15")

beta_df %>% mutate(beta = round(beta, digits = 3)) %>% 
  ggplot() + aes(x = age_group, y = beta) + 
  geom_point(size = 5, shape = 18) + 
  geom_errorbar(aes(ymin=beta_lower, ymax=beta_upper), width = 0.5) + ylim(0.10, 0.125) +
  labs(title = 'Beta estimates by age group, with bootstrap 95% percentile interval', x = 'Age Group', y = 'Beta') +
  geom_text(mapping = aes(label = beta, y = 0.107), size = 7)


# simulated vaccination outcomes
capop = 39500000
f_beta = function(t){0.117}
# Define the vaccination function with coefficients from a prior linear regression of CDC
# vaccine data, flattening in a change point on Day 105(April 15th)
f_vacc_linear = function(t){
  if(t<=105){34089+3625*t}
  else{414714}}

sixtyfiveplus = 5850000
# Estimate init_exp as number of new infections from Jan 1 to Jan 7
init_exp_old = 267362 - 243143
ifr_old = 0.09
init_rec_old = 243143*2
# Estimate init_inf as number of new infections from Dec 21 to Jan 1
init_inf_old = 243143 - 201338
init_dead_old = 19545
f_vacc_line = function(t){if(t<15){0}
  else{f_vacc_linear(t-14)/2}
}
seir_old = seir(sixtyfiveplus,init_inf_old,init_exp_old,
                init_rec_old, init_dead_old,f_beta,ifr=ifr_old,
                vacc_f=f_vacc_line,days=200)

midagepop = 6900000
init_exp_mid = 487547 - 444733
ifr_mid = 0.006
init_rec_mid = 444733*2
init_inf_mid = 444733 - 362733
init_dead_mid = 4874
f_vacc_linear = function(t){34089+3625*t}
# Complete prioritization means seniors are vaccinated by Day 63; younger age groups receive no vaccine until then.
f_vacc_lineM = function(t){if(t<(15+48)){0}
  else{f_vacc_linear(t-14)/2}
}
seir_mid = seir(midagepop,init_inf_mid,init_exp_mid,
                init_rec_mid, init_dead_mid,f_beta,ifr=ifr_mid,
                vacc_f=f_vacc_lineM,days=200)

youngadultpop = 17750000
init_exp_ya = 1497244 - 1372058
ifr_ya = 0.0005
init_rec_ya = 1372058*2
init_inf_ya = 1372058 - 1134542
init_dead_ya = 1807
# Complete prioritization means those 50+ are vaccinated by Day 89; younger age groups receive no vaccine until then.
f_vacc_lineY = function(t){if(t<(15+74)){0}
  else{f_vacc_linear(t-14)/2}
}
seir_ya = seir(youngadultpop,init_inf_ya,init_exp_ya,
               init_rec_ya, init_dead_ya,f_beta,ifr=ifr_ya,
               vacc_f=f_vacc_lineY,days=200)

childpop = 9000000
init_exp_c = 314694 - 284404
ifr_c = 0.00002
init_rec_c = 284404*2
init_inf_c = 284404 - 229284
init_dead_c = 6
# Children (outside of a few 16-17 for Pfizer) are not approved to get vaccinated, thus our functional vaccination rate for them is 0.
f_vacc_lineC = function(t){0}
seir_c = seir(childpop,init_inf_c,init_exp_c,
              init_rec_c, init_dead_c,f_beta,ifr=ifr_c,
              vacc_f=f_vacc_lineC,days=200)

# Our total prioritized model adds together the values from all age brackets.
seir_priority = seir_old+seir_mid+seir_ya+seir_c
seir_priority$Day = 1:201
seir_priority1 = seir_priority[,c(1,3,4,6)]
seir_long = seir_priority1 %>% gather(key="group",value="value",-Day)
priorplot = ggplot(seir_long,aes(x=Day,y=value,color=group)) + geom_line()+
  labs(x="Day",y="People In Group",title="Population, Vaccines Prioritized")

init_exp_t = sum(init_exp_old,init_exp_mid,init_exp_ya,init_exp_c)
ifr_t = (sixtyfiveplus*ifr_old+midagepop*ifr_mid+youngadultpop*ifr_ya+childpop*ifr_c)/capop
init_rec_t = sum(init_rec_old,init_rec_mid,init_rec_ya,init_exp_c)
init_inf_t = sum(init_inf_old,init_inf_mid,init_inf_ya,init_inf_c)
init_dead_t = sum(init_dead_old,init_dead_mid,init_dead_ya,init_dead_c)

seir_equal = seir(capop,init_inf_t,init_exp_t,
                  init_rec_t, init_dead_t,f_beta,ifr=ifr_t,
                  vacc_f=f_vacc_line,days=200)
seir_equal1 = seir_equal[,c(1,3,4,6)]
seir_long = seir_equal1 %>% gather(key="group",value="value",-Day)
eqplot = ggplot(seir_long,aes(x=Day,y=value,color=group)) + geom_line()+
  labs(x="Day",y="People In Group",title="Population, Vaccines Not Prioritized")

## Upper Bounds

f_betahigh = function(t){0.121}
f_betahigh2 = function(t){0.119}

useir_old = seir(sixtyfiveplus,init_inf_old,init_exp_old,
                 init_rec_old, init_dead_old,f_betahigh,ifr=ifr_old,
                 vacc_f=f_vacc_line,days=200)
useir_mid = seir(midagepop,init_inf_mid,init_exp_mid,
                 init_rec_mid, init_dead_mid,f_betahigh,ifr=ifr_mid,
                 vacc_f=f_vacc_lineM,days=200)
useir_ya = seir(youngadultpop,init_inf_ya,init_exp_ya,
                init_rec_ya, init_dead_ya,f_betahigh2,ifr=ifr_ya,
                vacc_f=f_vacc_lineY,days=200)
useir_c = seir(childpop,init_inf_c,init_exp_c,
               init_rec_c, init_dead_c,f_betahigh,ifr=ifr_c,
               vacc_f=f_vacc_lineC,days=200)
useir_priority = useir_old+useir_mid+useir_ya+seir_c
useir_priority$Day = 1:201
useir_priority1 = useir_priority[,c(1,3,4,6)]
useir_long = useir_priority1 %>% gather(key="group",value="value",-Day)

## Low Bounds

f_betalow = function(t){0.11}

lseir_old = seir(sixtyfiveplus,init_inf_old,init_exp_old,
                 init_rec_old, init_dead_old,f_betalow,ifr=ifr_old,
                 vacc_f=f_vacc_line,days=200)
lseir_mid = seir(midagepop,init_inf_mid,init_exp_mid,
                 init_rec_mid, init_dead_mid,f_betalow,ifr=ifr_mid,
                 vacc_f=f_vacc_lineM,days=200)
lseir_ya = seir(youngadultpop,init_inf_ya,init_exp_ya,
                init_rec_ya, init_dead_ya,f_betalow,ifr=ifr_ya,
                vacc_f=f_vacc_lineY,days=200)
lseir_c = seir(childpop,init_inf_c,init_exp_c,
               init_rec_c, init_dead_c,f_betalow,ifr=ifr_c,
               vacc_f=f_vacc_lineC,days=200)
lseir_priority = lseir_old+lseir_mid+lseir_ya+lseir_c
lseir_priority$Day = 1:201
lseir_priority1 = lseir_priority[,c(1,3,4,6)]
lseir_long = lseir_priority1 %>% gather(key="group",value="value",-Day)
priorplot + geom_line(data=lseir_long,linetype="dashed") + geom_line(data=useir_long,linetype="dashed")

lseir_equal = seir(capop,init_inf_t,init_exp_t,
                   init_rec_t, init_dead_t,f_betalow,ifr=ifr_t,
                   vacc_f=f_vacc_line,days=200)
lseir_equal1 = lseir_equal[,c(1,3,4,6)]
lseir_long = lseir_equal1 %>% gather(key="group",value="value",-Day)
useir_equal = seir(capop,init_inf_t,init_exp_t,
                   init_rec_t, init_dead_t,f_betahigh,ifr=ifr_t,
                   vacc_f=f_vacc_line,days=200)


useir_equal1 = useir_equal[,c(1,3,4,6)]
useir_long = useir_equal1 %>% gather(key="group",value="value",-Day)
eqplot + geom_line(data=lseir_long,linetype="dashed") + geom_line(data=useir_long,linetype="dashed")



# actual vaccination outcomes
seir_measured <- measuredVacc(measure_data)

seir_measured %>% ggplot() + geom_line(aes(x = date, y = infected, color = 'Infected')) + 
  geom_line(aes(x = date, y = deathsum, color = 'Dead')) + 
  labs(title = 'Population, Measured since 01-01-2021', 
       x = 'Date', y = 'Number of People in Group')


