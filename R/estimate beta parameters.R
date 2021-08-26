
#################################################################
# estimate beta parameters given 95% confidence interval and mean
#################################################################

###################################
#### HAZARD RATIO# mortality: Tafa vs BSC ####
### Mean = 0.7; CI: 0.51 - 0.96 ###
###################################
objective.function <- function(params) {
  alpha <- params[1]
  beta <- params[2]
  
  intended.quantiles <- c(0.51, 0.96)
  calculated.quantiles <- qbeta(p=c(0.025, 0.975), shape1=alpha, shape2=beta)
  squared.error.quantiles <- sum((intended.quantiles - calculated.quantiles)^2)
  
  intended.mean<- 0.70
  calculated.mean <- calculate.mean(alpha, beta)
  squared.error.mean <- (intended.mean - calculated.mean)^2
  
  return(squared.error.quantiles + squared.error.mean)
}

calculate.mean <- function(alpha, beta) {
  return((alpha-1) / (alpha+beta-2))
}


starting.params <- c(20, 200)
#starting params?

nlm.result <- nlm(f = objective.function, p = starting.params)
optimal.alpha <- nlm.result$estimate[1]
optimal.beta <- nlm.result$estimate[2]

optimal.alpha
optimal.beta 


qbeta(p=c(0.025, 0.975), shape1=optimal.alpha, shape2=optimal.beta)

calculate.mean(optimal.alpha, optimal.beta)

##########################################################
#### HAZARD RATIO# cv hospitalizations : Tafa vs BSC ####
### Mean = 0.68; CI, 0.56 to 0.81 ###
#########################################################
objective.function <- function(params) {
  alpha <- params[1]
  beta <- params[2]
  
  intended.quantiles <- c(0.56, 0.81)
  calculated.quantiles <- qbeta(p=c(0.025, 0.975), shape1=alpha, shape2=beta)
  squared.error.quantiles <- sum((intended.quantiles - calculated.quantiles)^2)
  
  intended.mean<- 0.68
  calculated.mean <- calculate.mean(alpha, beta)
  squared.error.mean <- (intended.mean - calculated.mean)^2
  
  return(squared.error.quantiles + squared.error.mean)
}

calculate.mean <- function(alpha, beta) {
  return((alpha-1) / (alpha+beta-2))
}


starting.params_hosp <- c(20, 200)
#starting params?

nlm.result_hosp <- nlm(f = objective.function, p = starting.params_hosp)
optimal.alpha_hosp <- nlm.result_hosp$estimate[1]
optimal.beta_hosp <- nlm.result_hosp$estimate[2]

optimal.alpha_hosp
optimal.beta_hosp


qbeta(p=c(0.025, 0.975), shape1=optimal.alpha_hosp, shape2=optimal.beta_hosp)

calculate.mean(optimal.alpha_hosp, optimal.beta_hosp)



