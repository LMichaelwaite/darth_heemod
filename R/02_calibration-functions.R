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
                xlab = "Time", ylab = "Prop. N1")



# TARGET 2: (if you had more...)
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

source("R/01_ATTR_model.R") # creates the function run_crs_markov()

library(heemod)
library(tidyverse)
library(ggplot2)
library(dplyr)
# Check that it works
v_params_test <- c(p_N1N1_6 = 0.6, p_N1N3_6 = 0.15)
test_results <- run_sick_sicker_markov(v_params_test) # It works!

str(test_results)

plot(test_results$prop_N1)

# compare test output to targets 
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



