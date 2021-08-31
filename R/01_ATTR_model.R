
library(heemod)
library(tidyverse)
library(ggplot2)
library(dplyr)

#### ATTR_CM Markov model ####
#-------------------------------------------------------------------------------
#                              DEFINE PARAMETERS
#...............................................................................
param <- define_parameters(
  
  fudge_up = 0,
  fudge_down = 0,
  p_N1N2_1 = 0.002536125,
  p_N2N3_1 = 0.003078164,
  
  p_N1N2_2 = 0.0905894,
  p_N2N3_2 = 0.2091134,
  
  p_N1N2_3 = 0.09210348 ,
  p_N2N3_3 = 0.06838006,
  
  p_N1N2_4 = p_N1N2_3,
  p_N2N3_4 = p_N1N2_3,
  
  p_N1D_1 = 0.002755 + (fudge_up/15), # month 0 -6
  p_N2D_1 = 0.005372,
  p_N3D_1 = 0.012067, 
  
  p_N1D_2 = 0.00241 + fudge_up, # month 6 - 12
  p_N2D_2 = 0.0590,
  p_N3D_2 = 0.1096,
  
  #p_N1D_3 = 0.103,
  #p_N2D_3 = 0.252,
  #p_N3D_3 = 0.441,
  
  p_N1D_3 = 0.05284029638 + fudge_up,# month 12 - 18
  p_N2D_3 = 0.135174477,
  p_N3D_3 = 0.2525989344,
  
  p_N1D_4 = p_N1D_3, # month 18 -24
  p_N2D_4 = p_N2D_3,
  p_N3D_4 = p_N3D_3,
  #p_N3D_4 = 0.23,
  #FROM N1
  p_N1N2 = ifelse(markov_cycle <= 1 ,p_N1N2_1,
           ifelse(markov_cycle >= 2 & markov_cycle < 3, p_N1N2_2,
           ifelse(markov_cycle >= 3 & markov_cycle < 4, p_N1N2_3,
            p_N1N2_4))),  #If greater than 24, then constant at p_N1N2_3. since NAC data only goes to 24.
  
  p_N1D = ifelse(markov_cycle <= 1 ,p_N1D_1,
          ifelse(markov_cycle >= 2 & markov_cycle < 3, p_N1D_2,
          ifelse(markov_cycle >= 3 & markov_cycle < 4, p_N1D_3,
          p_N1D_4))),
                               
  #FROM N2
  p_N2N3 = ifelse(markov_cycle <= 1 ,p_N2N3_1,
           ifelse(markov_cycle >= 2 & markov_cycle < 3, p_N2N3_2,
           ifelse(markov_cycle >= 3 & markov_cycle < 4, p_N2N3_3,
            p_N2N3_4))),
                                
  
  p_N2D = ifelse(markov_cycle <= 1 ,p_N2D_1,
          ifelse(markov_cycle >= 2 & markov_cycle < 3, p_N2D_2,
          ifelse(markov_cycle >= 3 & markov_cycle < 4, p_N2D_3,
          p_N2D_4))),
                              
  #FROM N3
  p_N3D = ifelse(markov_cycle <= 1 ,p_N3D_1,
          ifelse(markov_cycle >= 2 & markov_cycle < 3, p_N3D_2,
          ifelse(markov_cycle >= 3 & markov_cycle < 4, p_N3D_3,
          p_N3D_4))),
  
  HR_6 = 0.64,  #HR after month 30, included for OWSA 
 
  #### TREATMENT EFFECT ####
  HR = ifelse(markov_cycle <= 5 ,0.7,
              HR_6),
  
  #### HOSPITALISATION ####
  #bsc cv hosp
  r_hosp = 0.70, # CV hosp rate per person per year
  p_hosp = 1- exp(-r_hosp*1/2), #CV hosp 6 month probability 
  
  #tafamidis cv hosp
  HR_hosp_tafa = 0.68,
  r_hosp_tafa = r_hosp * HR_hosp_tafa, # CV hosp rate per person per year
  p_hosp_tafa = 1- exp(-r_hosp_tafa*1/2), # CV hosp 6 month probability
  
  disutility_hosp = 0.10, # https://onlinelibrary.wiley.com/doi/full/10.1002/ehf2.12844
  
  #### COSTS AND UTILITIES ####  
  
  # EXPERT ELLICITATION: weighting for NYHA to NAC
  # NAC: 
  #divde by two due to 6 month cycle
  u_NAC_I = 0.893/2, # Source: Rozenbaum NYHA utilities. Assume NYHA I = NAC I,..., AVG(NYHA_III, NYHA_IV) = NAC_III
  u_NAC_II = 0.802/2,
  u_NAC_III = (0.7*0.706+0.3*0.406)/2, 
  u_NAC_I_tafa = 0.874/2,
  u_NAC_II_tafa = 0.832/2,
  u_NAC_III_tafa = (0.7*0.707 + 0.3*0.558)/2,
  c_tafa = 10840.82*6, # cost per 6 months; GBP; source: NICE tafamidis STA pp. 137
  c_bsc = 7031.92/2, # GBP; source nuffield trust estimate for health care sepnding on 75 year olds (UK) inflation adjusted to 2020 ((find actual source)) https://www.theguardian.com/society/2016/feb/01/ageing-britain-two-fifths-nhs-budget-spent-over-65s 
  #c_bsc = 18462/2, # background hc cost 6 months; dollars US Kazi
  c_hosp = 2536.88, # GBP source: NICE tafamidis STA pp.140
  c_eol = 9287.86, # # GBP source: NICE tafamidis STA pp.140, end of life costs
  R = 0.03  #discount rate
)
#-------------------------------------------------------------------------------
#                       DEFINE TRANSITION MATRIX
#...............................................................................
mat_bsc <- define_transition(state_names = c(
  "NAC1",
  "NAC2",
  "NAC3",
  "Death"
),

  1-(p_N1N2 + p_N1D),                 p_N1N2,           0,  p_N1D,
                   0,     1- (p_N2N3+ p_N2D),      p_N2N3,  p_N2D,
                   0,                      0,   1-(p_N3D),  p_N3D,
                   0,                      0,           0,   1.00
)

mat_tafa <- define_transition(state_names = c(
  "NAC1",
  "NAC2",
  "NAC3",
  "Death"
),
 
  1-(p_N1N2 + p_N1D*HR),                 p_N1N2,           0,  p_N1D*HR,
                   0,     1- (p_N2N3+ p_N2D*HR),      p_N2N3,  p_N2D*HR,
                   0,                      0,   1-(p_N3D*HR),  p_N3D*HR,
                   0,                      0,              0,   1.00
)

plot(mat_bsc) #model diagram
plot(mat_tafa) #model diagram

#-------------------------------------------------------------------------------
#           DEFINE strategy & state dependent costs and utilities
#...............................................................................
R <- (1+0.03)^(1/3) - 1
NAC1 <- define_state(
  cost_background_health = discount(c_bsc, R),
  cost_drugs = discount(dispatch_strategy(
      bsc = 0,
      tafa = c_tafa, 
  ), R),
  cost_hosp = discount(dispatch_strategy(
      bsc = p_hosp * c_hosp,
      tafa = p_hosp_tafa* c_hosp), R),
  cost_total = discount(cost_background_health + cost_drugs + cost_hosp, 0),
  
  qaly = dispatch_strategy(
      bsc = u_NAC_I,
      tafa = u_NAC_I_tafa
      ),
  cv_hosp_disutility = dispatch_strategy(
      bsc = p_hosp * disutility_hosp,
      tafa = p_hosp_tafa* disutility_hosp),
  qaly_total =  discount(qaly-cv_hosp_disutility, R),
  life_year = discount(0.5, R)
)
NAC2 <- define_state(
  cost_background_health = discount(c_bsc, R),
  cost_drugs = discount(dispatch_strategy(
        bsc = 0,
        tafa = c_tafa, 
      ), R),
  cost_hosp = discount(dispatch_strategy(
      bsc = p_hosp * c_hosp,
      tafa = p_hosp_tafa* c_hosp), R),
  cost_total = discount(cost_background_health + cost_drugs + cost_hosp, 0),
  
  qaly = dispatch_strategy(
      bsc = u_NAC_II,
      tafa = u_NAC_II_tafa
  ),
  cv_hosp_disutility = dispatch_strategy(
      bsc = p_hosp * disutility_hosp,
      tafa = p_hosp_tafa* disutility_hosp),
  qaly_total =  discount(qaly- cv_hosp_disutility, R),
  life_year = discount(0.5, R)
)
NAC3 <- define_state(
  cost_background_health = discount(c_bsc, R),
  cost_drugs = discount(dispatch_strategy(
    bsc = 0,
    tafa = c_tafa, 
  ), R),
  cost_hosp = discount(dispatch_strategy(
    bsc = p_hosp * c_hosp,
    tafa = p_hosp_tafa* c_hosp), R),
  cost_total = discount(cost_background_health + cost_drugs + cost_hosp, 0),

  qaly = dispatch_strategy(
    bsc = u_NAC_III,
    tafa = u_NAC_III_tafa
  ),
  cv_hosp_disutility = dispatch_strategy(
    bsc = p_hosp * disutility_hosp,
    tafa = p_hosp_tafa* disutility_hosp),
  qaly_total =  discount(qaly - cv_hosp_disutility, R),
  life_year = discount(0.5, R)
)
Death <- define_state(
  cost_background_health = 0,
  cost_drugs = 0,
  cost_hosp =0,
  cost_total = discount(cost_background_health + cost_drugs + cost_hosp, .03),
  qaly = 0,
  cv_hosp_disutility = 0,
  qaly_total =  discount(qaly - cv_hosp_disutility, R),
  life_year = discount(0, R)
)
strat_bsc <- define_strategy(
  transition = mat_bsc,
  NAC1 = NAC1,
  NAC2 = NAC2,
  NAC3 = NAC3,
  Death = Death
)      
strat_tafa <- define_strategy(
  transition = mat_tafa,
  NAC1 = NAC1,
  NAC2 = NAC2,
  NAC3 = NAC3,
  Death = Death
)  
res_mod <- run_model(
  bsc = strat_bsc,
  tafa = strat_tafa,
  parameters = param,
  cycles = 40,
  cost = cost_total,
  #effect = qaly_total,
  effect = life_year,
  method = "end",
  init = c(436, 350, 159, 0)
)
summary(res_mod)
plot(res_mod) # count per state     
plot(res_mod, type = "ce") # CE frontier

#-------------------------------------------------------------------------------
#                                  PSA
#...............................................................................
beta_params <- function (mean, sigma) #function to calculate beta params from the mean and se
{
  alpha <- ((1 - mean)/sigma^2 - 1/mean) * mean^2
  beta <- alpha * (1/mean - 1)
  params<- list(alpha = alpha, beta = beta)
  return(params)
}  
      
gamma_params <- function (mu, sigma, scale = TRUE)  #function to calculate gamma params from the mean and se
{
  if (scale) {
    shape <- (mu^2)/(sigma^2)
    scale <- (sigma^2)/mu
    params <- list(shape = shape, scale = scale)
  }
  else {
    shape <- (mu^2)/(sigma^2)
    rate <- mu/(sigma^2)
    params <- list(shape = shape, rate = rate)
  }
  return(params)
}

beta_u_NAC_I <-beta_params(mean = .893/2, sigma = .02)
beta_u_NAC_II <-beta_params(mean = .802/2, sigma = .01)
beta_u_NAC_III <-beta_params(mean = ((0.706+0.406)/2)/2, sigma = .04)
beta_u_NAC_I_tafa <-beta_params(mean = .874/2, sigma = .01)
beta_u_NAC_II_tafa <-beta_params(mean = .832/2, sigma = .01)
beta_u_NAC_III_tafa <-beta_params(mean = ((0.707+0.558)/2)/2, sigma = .04) 

gamma_c_hosp <- gamma_params(mu = 2546, sigma = 507.38)

rsp <- define_psa(
  u_NAC_I ~ beta(beta_u_NAC_I$alpha, beta_u_NAC_I$beta),
  u_NAC_II ~ beta(beta_u_NAC_II$alpha, beta_u_NAC_II$beta),
  u_NAC_III ~ beta(beta_u_NAC_III$alpha, beta_u_NAC_III$beta),
  u_NAC_I_tafa ~ beta(beta_u_NAC_I_tafa$alpha, beta_u_NAC_I_tafa$beta),
  u_NAC_II_tafa ~ beta(beta_u_NAC_II_tafa$alpha, beta_u_NAC_II_tafa$beta),
  u_NAC_III_tafa ~ beta(beta_u_NAC_III_tafa$alpha, beta_u_NAC_III_tafa$beta),
  
  disutility_hosp ~ beta(10, 90), # source: https://onlinelibrary.wiley.com/doi/full/10.1002/ehf2.12844
  
  p_N1N2_1 ~ beta(param$p_N1N2_1, param$p_N1N2_1*0.15),
  p_N2N3_1 ~ beta(param$p_N2N3_1, param$p_N2N3_1*0.15),
  p_N1N2_2 ~ beta(param$p_N1N2_2, param$p_N1N2_2*0.15),
  p_N2N3_2 ~ beta(param$p_N2N3_2, param$p_N2N3_2*0.15),
  p_N1N2_3 ~ beta(param$p_N1N2_3, param$p_N1N2_3*0.15),
  p_N2N3_3 ~ beta(param$p_N2N3_3, param$p_N2N3_3*0.15),
  
  p_N1D_1 ~ beta(param$p_N1D_1, param$p_N1D_1*0.15),   # alpha = count of events; beta = number at risk - alpha 
  p_N2D_1 ~ beta(param$p_N2D_1, param$p_N2D_1*0.15),
  p_N3D_1 ~ beta(param$p_N3D_1, param$p_N3D_1*0.15),
  p_N1D_2 ~ beta(param$p_N1D_2, param$p_N1D_2*0.15),   
  p_N2D_2 ~ beta(param$p_N2D_2, param$p_N2D_2*0.15),
  p_N3D_2 ~ beta(param$p_N3D_2, param$p_N3D_2*0.15),
  p_N1D_3 ~ beta(param$p_N1D_3, param$p_N1D_3*0.15),   
  p_N2D_3 ~ beta(param$p_N2D_3, param$p_N2D_3*0.15),
  p_N3D_3 ~ beta(param$p_N3D_3, param$p_N3D_3*0.15),
  p_N1D_4 ~ beta(param$p_N1D_4, param$p_N1D_4*0.15),   
  p_N2D_4 ~ beta(param$p_N2D_4, param$p_N2D_4*0.15),
  p_N3D_4 ~ beta(param$p_N3D_4, param$p_N3D_4*0.15),
 
  c_bsc ~ gamma(mean = 7031.92/2, sd = sqrt((7031.92/2))), # "gamma" seems to not require the above gamma function. ?
  c_hosp ~ gamma(mean = 2546, sd= 507.38),
  
  HR ~ beta(12.02719, 4.639787), #from script "estimate beta parameters"
  HR_hosp_tafa ~ beta (36.16638, 16.63898)  #from script "estimate beta parameters"
  
  #TO add:
  #p_N1N2
  #p_N2N3
  #HR_6
  #c_bsc
)

#"run_psa" is masked from dampack, must specify heemod
pm <- run_psa( 
  model = res_mod,
  psa = rsp,
  N = 500
)   

summary(pm)
plot(pm, type = "ac", max_wtp = 1000000, log_scale = FALSE)        
plot(pm, type = "cov")
plot(pm, type = "cov", diff = TRUE, threshold = 10000)
plot(pm, type = "ce")
#library(BCEA)
#bcea <- run_bcea(pm, plot = TRUE, Kmax = 1500000)

#-------------------------------------------------------------------------------
#                                   OWSA
#...............................................................................
## best case; HR decreases to 0.64 after 30 months (from the ATTR-ACT extension)
      
## base case; HR remains constant over time 
se <- define_dsa(
  #u_NAC_I, 0.2, 0.5,
  #u_NAC_II, 0.2, 0.5,
  #u_NAC_III, 0.1, 0.5,
  #u_NAC_I_tafa, 0.2, 0.5,
  #u_NAC_II_tafa, 0.2, 0.5, 
  #u_NAC_III_tafa, 0.2, 0.5, 
  
  fudge_up, 0.01, 0.1,
  
  disutility_hosp, 0.025, 0.225, 

  p_N1D_1, 0, 0.5,     
  p_N2D_1, 0, 0.5, 
  p_N3D_1, 0, 0.5, 
  
  p_N1D_2, 0, 0.5, 
  p_N2D_2, 0, 0.5, 
  p_N3D_2, 0, 0.5, 
  
  p_N1D_3, 0, 0.5, 
  p_N2D_3, 0, 0.5, 
  p_N3D_3, 0, 0.5, 
  
  p_N1D_4, 0, 0.5, 
  p_N2D_4, 0, 0.5, 
  p_N3D_4, 0, 0.5, 

  p_N2N3_3, 0, 0.1,
  p_N2N3_4, 0, 0.11,
  #p_N1D_1, param$p_N1D_1*0.33, param$p_N1D_1*3,     
  #p_N2D_1, param$p_N2D_1*0.33, param$p_N2D_1*3,
  #p_N3D_1, param$p_N3D_1*0.33, param$p_N3D_1*3,
  #
  #p_N1D_2, param$p_N1D_2*0.33, param$p_N1D_2*3,
  #p_N2D_2, param$p_N2D_2*0.33, param$p_N2D_2*3,
  #p_N3D_2, param$p_N3D_2*0.33, param$p_N3D_2*3,
  #
  #p_N1D_3, param$p_N1D_3*0.33, param$p_N1D_3*3,
  #p_N2D_3, param$p_N2D_3*0.33, param$p_N2D_3*3,
  #p_N3D_3, param$p_N3D_3*0.33, param$p_N3D_3*3,
  #
  #p_N1D_4, param$p_N1D_4*0.33, param$p_N1D_4*3,
  #p_N2D_4, param$p_N2D_4*0.33, param$p_N2D_4*3,
  #p_N3D_4, param$p_N3D_4*0.33, param$p_N3D_4*3,

  c_bsc, (7031.92/2)*0.33, (7031.92/2)*3,
  c_hosp,  2546*0.33,  2546*3,
  c_tafa, 100, 100000
  
  #HR, 0.5, 0.9, #from script "estimate beta parameters"
  #HR_hosp_tafa, 0.45, 0.85  #from script "estimate beta parameters"
)

res_dsa <- run_dsa(
  model = res_mod,
  dsa = se
)

res_dsa

plot(res_dsa,
     strategy = "tafa",
     result = "effect",
     type = "difference")  

plot(res_dsa, 
     strategy = "bsc",
     result = "effect",
     type = "simple")

plot(res_dsa, 
     strategy = "tafa",
     result = "effect",
     type = "simple")
## worse case; HR is zero after month 30      
    
#-------------------------------------------------------------------------------
#                       Extra tables and graphs
#...............................................................................
c_n <- sum(c(436, 350, 159, 0))   # cohort size    

# Tidy data
#get state counts
c_state <- as.data.frame(get_counts(res_mod)) 
c_state_wide <- pivot_wider(c_state, names_from = state_names, values_from = count) #convert long table into wide table

#### Table 1: state distribution ####
#BSC
nac_dist <-c_state_wide %>%
  rowwise() %>%
  mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
  mutate(proportion_alive = total_alive/rowSums(across(Death:total_alive)))%>%
  mutate(proportion_N1= NAC1/total_alive)%>%
  mutate(proportion_N2= NAC2/total_alive)%>%
  mutate(proportion_N3= NAC3/total_alive)%>%
  filter(.strategy_names =="bsc")%>%
  select(.strategy_names,markov_cycle, proportion_alive, proportion_N1,proportion_N2, proportion_N3)

#Tafmidis
nac_dist_tafa <-c_state_wide %>%
  rowwise() %>%
  mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
  mutate(proportion_alive = total_alive/rowSums(across(Death:total_alive)))%>%
  mutate(proportion_N1= NAC1/total_alive)%>%
  mutate(proportion_N2= NAC2/total_alive)%>%
  mutate(proportion_N3= NAC3/total_alive)%>%
  filter(.strategy_names =="tafa")%>%
  select(.strategy_names,markov_cycle, proportion_alive, proportion_N1,proportion_N2, proportion_N3)

#View table 1:
nac_dist 
nac_dist_tafa

#### Table 2: mean LY ####    EDIT: better approach in table 3 
#BSC
#life_years_bsc <-c_state_wide %>%
#  rowwise() %>%
#  mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
#  filter(.strategy_names =="bsc")%>%
#  select(.strategy_names,markov_cycle, total_alive)
#
#avg_life_years_bsc <- sum(life_years_bsc$total_alive/2)/c_n    #c_n cohort size
#avg_life_years_bsc
#
##BSC
#life_years_tafa <-c_state_wide %>%
#  rowwise() %>%
#  mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
#  filter(.strategy_names =="tafa")%>%
#  select(.strategy_names,markov_cycle, total_alive)
#
#avg_life_years_tafa <- sum(life_years_tafa$total_alive/2)/c_n    #c_n cohort size
#avg_life_years_tafa

#### Table 3: mean outcomes ##### 
#### Prepare data: value total by cycle 
df_values <- as.data.frame(get_values(res_mod))
df_values_bsc <- df_values %>% # BSC
  filter(.strategy_names =="bsc") 
df_values_tafa <- df_values %>% # Tafa
  filter(.strategy_names =="tafa") 

values_bsc_wide <- pivot_wider(df_values_bsc, names_from = value_names, values_from = value)
values_tafa_wide <- pivot_wider(df_values_tafa, names_from = value_names, values_from = value)

## n prop
#nac_dist <-c_state_wide %>%
#  rowwise() %>%
#  mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
#  mutate(proportion_alive = total_alive/rowSums(across(Death:total_alive)))%>%
#  mutate(proportion_N1= NAC1/total_alive)%>%
#  mutate(proportion_N2= NAC2/total_alive)%>%
#  mutate(proportion_N3= NAC3/total_alive)%>%
#  filter(.strategy_names =="bsc")%>%
#  select(total_alive, proportion_N1,proportion_N2, proportion_N3)
#
#
#values_bsc_1 <- cbind(values_bsc_wide, nac_dist)
#
#
#cal_tar <- values_bsc_1  %>%
#  rowwise() %>%
#  mutate(targ = qaly_total/life_year) %>%
#  filter(markov_cycle %in% c(2, 20)) %>%
#  select(targ)
#
#z <- values_bsc_1  %>%
#  rowwise() %>%
#  mutate(targ = qaly_total/life_year) %>%
#  
#  select(targ)  
#
#
#qaly_check <- values_bsc_1 %>%
#  rowwise() %>%
#  mutate(qaly_1 = qaly/total_alive) %>%
#  select(markov_cycle, qaly_1)
#

#### create column: cohort size per cycle
c_alive_bsc <-c_state_wide %>%
  rowwise() %>%
  mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
  filter(.strategy_names =="bsc")%>%
  select(total_alive)

c_alive_tafa <-c_state_wide %>%
  rowwise() %>%
  mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
  filter(.strategy_names =="tafa")%>%
  select(total_alive)

values_bsc <- cbind(values_bsc_wide, c_alive_bsc)
values_tafa <- cbind(values_tafa_wide, c_alive_tafa)        
        
#### Calculate mean values ####        
avg_qaly_bsc <-  sum(values_bsc_wide$qaly_total)/c_n        
avg_ly_bsc <-  sum(values_bsc_wide$life_year)/c_n        
avg_cost_bsc <-  sum(values_bsc_wide$cost_total)/c_n         
avg_cost_hosp_bsc <- sum(values_bsc_wide$cost_hosp)/c_n

avg_qaly_tafa <-  sum(values_tafa_wide$qaly_total)/c_n        
avg_ly_tafa <-  sum(values_tafa_wide$life_year)/c_n        
avg_cost_tafa <-  sum(values_tafa_wide$cost_total)/c_n         
avg_cost_hosp_tafa <- sum(values_tafa_wide$cost_hosp)/c_n

avg_drug_cost_tafa <- sum(values_tafa_wide$cost_drugs)/c_n

incremental_cost <- avg_cost_tafa - avg_cost_bsc
incremental_ly <- avg_ly_tafa - avg_ly_bsc
incremental_qaly <- avg_qaly_tafa - avg_qaly_bsc

result_summary <- data.frame(QALY_tafa =  avg_qaly_tafa,
                             QALY_bsc =  avg_qaly_bsc,
                             Life_years_tafa = avg_ly_tafa,
                             Life_years_bsc = avg_ly_bsc,
                             cost_tafa =  avg_cost_tafa,
                             cost_bsc = avg_cost_bsc,
                             drug_cost_tafa = avg_drug_cost_tafa,
                             hosp_cost_tafa = avg_cost_hosp_tafa,
                             hosp_cost_bsc = avg_cost_hosp_bsc,
                             incremental_cost = incremental_cost,
                             incremental_ly = incremental_ly,
                             incremental_qaly = incremental_qaly)       
result_summary        
        
#----------------------------------------------------------
#plot survival stuff
#-------------------------------------------------------------
#plot death state count
plot(res_mod, type = "counts", states = "Death", panel = "by_state", free_y = TRUE) +
  theme_bw() +
  scale_color_brewer(
    name = "Strategy",
    palette = "Set1"
  )

#get state counts
c_state <- as.data.frame(get_counts(res_mod))

c_state_bsc <- filter(c_state, .strategy_names == "bsc")
c_state_tafamidis <- filter(c_state, .strategy_names == "tafamidis")


c_state_wide_2 <- pivot_wider(c_state_bsc, names_from = state_names, values_from = count)

#### BSC survival ####
#convert long table into wide table
c_state_wide <- pivot_wider(c_state, names_from = c(state_names,.strategy_names), values_from = count) 

#create new columns: total_alive & proportion_alive
surv_prop <-c_state_wide %>%
  rowwise() %>%
  mutate(total_alive_tafa = sum(across(c(NAC1_tafa, NAC2_tafa, NAC3_tafa)), na.rm = T)) %>%
  mutate(total_dead_tafa = Death_tafa) %>%
  mutate(proportion_alive = total_alive_tafa/rowSums(across(Death_tafa:total_alive_tafa)))%>%
  select(markov_cycle, proportion_alive)
surv_prop

surv_prop_bsc <- c_state_wide %>%
  rowwise() %>%
  mutate(total_alive_bsc = sum(across(c(NAC1_bsc, NAC2_bsc, NAC3_bsc)), na.rm = T)) %>%
  mutate(proportion_alive = total_alive_bsc/rowSums(across(Death_bsc:total_alive_bsc))) %>%
  mutate(total_dead_bsc = Death_bsc) %>%  
  select(markov_cycle, proportion_alive, total_alive_bsc, total_dead_bsc)
surv_prop_bsc 


# Kazi et al
kazi_surv_tafa<- read.csv("data_raw/kazi_tafa.csv", header = FALSE)
kazi_x_tafa<-kazi_surv_tafa[,1]/6
kazi_y_tafa<-kazi_surv_tafa[,2]
kazi_tafa <- data.frame(kazi_x = kazi_x_tafa, kazi_y = kazi_y_tafa)

kazi_surv_bsc<- read.csv("data_raw/kazi_control.csv", header = FALSE)
kazi_x_bsc<-kazi_surv_bsc[,1]/6
kazi_y_bsc<-kazi_surv_bsc[,2]
kazi_bsc <- data.frame(kazi_x = kazi_x_bsc, kazi_y = kazi_y_bsc)



#plot surv_prop curve
ggplot() + 
  geom_line(data=surv_prop, aes(x=markov_cycle, y=proportion_alive), color='green') + 
  geom_line(data=surv_prop_bsc, aes(x=markov_cycle, y=proportion_alive), color='red') +
  geom_line(data=kazi_tafa, aes(x=kazi_x, y=kazi_y), color='blue') +
  geom_line(data=kazi_bsc, aes(x=kazi_x, y=kazi_y), color='orange') 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
         
        
        
   