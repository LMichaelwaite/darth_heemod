###########################################################################
# This code was created by the DARTH workgroup (www.darthworkgroup.com). 
# When using or modifying this code, please do so with attribution and 
# cite our publications:

# - Alarid-Escudero F, Maclehose RF, Peralta Y, Kuntz KM, Enns EA. 
#   Non-identifiability in model calibration and implications for 
#   medical decision making. Med Decis Making. 2018; 38(7):810-821.

# - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E, 
#   Hunink MG. An Overview of R in Health Decision Sciences. 
#   Med Decis Making. 2017; 37(3): 735-746. 

# A walkthrough of the code could be found in the follwing link:
# - https://darth-git.github.io/calibSMDM2018-materials/
###########################################################################

###################  Calibration Specifications  ###################

# Model: 3-State Cancer Relative Survival (CRS) Markov Model
# Inputs to be calibrated: p_Mets, p_DieMets
# Targets: Surv

# Calibration method: Incremental mixture importance sampling (IMIS)
# Goodness-of-fit measure: Sum of Log-Likelihood

####################################################################


####################################################################
######  Load packages and function files  ####
####################################################################
# calibration functionality
library(lhs)
library(IMIS)
library(matrixStats) # package used for sumamry statistics

# visualization
library(plotrix)
library(psych)


####################################################################
######  Load target data  ######
####################################################################
load("data/ATTR_CalibTargets.RData")
lst_targets <- ATTR_targets
lst_targets


# Plot the targets

# TARGET 1: NAC 1 proportion
plotrix::plotCI(x = lst_targets$dist$N1$time, y = lst_targets$dist$N1$value, 
                ui = lst_targets$dist$N1$ub,
                li = lst_targets$dist$N1$lb,
                ylim = c(0, 1),
                xlim = c(0, 30),
                xlab = "Time", ylab = "Prop. N1")

# TARGET 2: NAC 2 proportion
plotrix::plotCI(x = lst_targets$dist$N2$time, y = lst_targets$dist$N2$value, 
                ui = lst_targets$dist$N2$ub,
                li = lst_targets$dist$N2$lb,
                ylim = c(0, 1),
                xlim = c(0, 30),
                xlab = "Time", ylab = "Prop. N2")

# TARGET 3: NAC 3 proportion
plotrix::plotCI(x = lst_targets$dist$N3$time, y = lst_targets$dist$N3$value, 
                ui = lst_targets$dist$N3$ub,
                li = lst_targets$dist$N3$lb,
                ylim = c(0, 1), 
                xlim = c(0, 30),
                xlab = "Time", ylab = "Prop. N3")


# plotrix::plotCI(x = lst_targets$Target2$time, y = lst_targets$Target2$value, 
#                 ui = lst_targets$Target2$ub,
#                 li = lst_targets$Target2$lb,
#                 ylim = c(0, 1), 
#                 xlab = "Time", ylab = "Target 2")


####################################################################
######  Load model as a function  ######
####################################################################
# - inputs are parameters to be estimated through calibration
# - outputs correspond to the target data

source("R/01_ATTR_model_new.R") # creates the function run_crs_markov()

library(heemod)
library(tidyverse)
library(ggplot2)
library(dplyr)
# Check that it works
v_params_test <- c(p_N1N2 = 0.2,
                   p_N1D = 0.05,
                   p_N2N3 = 0.15,
                   p_N2D = 0.05
                  
)


test_results <- run_sick_sicker_markov2(v_params_test) # It works!

str(test_results)

plot(test_results$prop_N1)

##### compare test output to targets ######
# TARGET 1: NAC 1 proportion
plotrix::plotCI(x = lst_targets$dist$N1$time, y = lst_targets$dist$N1$value, 
                ui = lst_targets$dist$N1$ub,
                li = lst_targets$dist$N1$lb,
                ylim = c(0, 1), 
                xlim = c(0, 30),
                xlab = "Time", ylab = "Prop. N1")
points(test_results$prop_N1,
       col = "green", pch = 20)
legend("bottomright", legend = c("Targets", "Outputs"),
       col = c("black", "green"),
       pch = c(1, 20))

# TARGET 2: NAC 2 proportion
plotrix::plotCI(x = lst_targets$dist$N2$time, y = lst_targets$dist$N2$value, 
                ui = lst_targets$dist$N2$ub,
                li = lst_targets$dist$N2$lb,
                ylim = c(0, 1), 
                xlim = c(0, 30),
                xlab = "Time", ylab = "Prop. N2")
points(test_results$prop_N2,
       col = "green", pch = 20)
legend("bottomright", legend = c("Targets", "Outputs"),
       col = c("black", "green"),
       pch = c(1, 20))

# TARGET 3: NAC 3 proportion
plotrix::plotCI(x = lst_targets$dist$N3$time, y = lst_targets$dist$N3$value, 
                ui = lst_targets$dist$N3$ub,
                li = lst_targets$dist$N3$lb,
                ylim = c(0, 1), 
                xlim = c(0, 30),
                xlab = "Time", ylab = "Prop. N3")
points(test_results$prop_N3,
       col = "green", pch = 20)
legend("bottomright", legend = c("Targets", "Outputs"),
       col = c("black", "green"),
       pch = c(1, 20))


####################################################################
######  Specify calibration parameters  ######
####################################################################
# Specify seed (for reproducible sequence of random numbers)
set.seed(072218)

# number of rand samples
n_resamp <- 1000

# names and number of input parameters to be calibrated
# names and number of input parameters to be calibrated

v_param_names <- c("p_N1N2",
                   
                   
                   "p_N1D", 
                   
                   
                   "p_N2N3",
                  
                   
                   "p_N2D"
                  
)

n_param <- length(v_param_names)

# range on input search space
#lb <- c(p_Mets = 0.04, p_DieMets = 0.04) # lower bound
#ub <- c(p_Mets = 0.16, p_DieMets = 0.16) # upper bound

lb <- c(p_N1N2 = 0.04,
       
        p_N1D = 0.04,
        
        p_N2N3 = 0.04,
        
        p_N2D = 0.04
       
)

ub <- c(p_N1N2 = 0.2,
        
        p_N1D = 0.2,
        
        p_N2N3 = 0.2,
       
        p_N2D = 0.2
        )



# number of calibration targets
v_target_names <- c("prop_N1", "prop_N2", "prop_N1" )
n_target <- length(v_target_names)


### Calibration functions
#  Write function to sample from prior
sample_prior <- function(n_samp){
  m_lhs_unit   <- randomLHS(n = n_samp, k = n_param)
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  colnames(m_param_samp) <- v_param_names
  for (i in 1:n_param){
    m_param_samp[, i] <- qunif(m_lhs_unit[,i],
                               min = lb[i],
                               max = ub[i])
    # ALTERNATIVE prior using beta (or other) distributions
    # m_param_samp[, i] <- qbeta(m_lhs_unit[,i],
    #                            shape1 = 1,
    #                            shape2 = 1)
  }
  return(m_param_samp)
}

# view resulting parameter set samples
pairs.panels(sample_prior(1000))


###  PRIOR  ### 
# Write functions to evaluate log-prior and prior

# function that calculates the log-prior
calc_log_prior <- function(v_params){
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  n_samp <- nrow(v_params)
  colnames(v_params) <- v_param_names
  lprior <- rep(0, n_samp)
  for (i in 1:n_param){
    lprior <- lprior + dunif(v_params[, i],
                             min = lb[i],
                             max = ub[i], 
                             log = T)
    # ALTERNATIVE prior using beta distributions
    # lprior <- lprior + dbeta(v_params[, i],
    #                          shape1 = 1,
    #                          shape2 = 1, 
    #                          log = T)
  }
  return(lprior)
}
calc_log_prior(v_params = v_params_test)
calc_log_prior(v_params = sample_prior(10))




