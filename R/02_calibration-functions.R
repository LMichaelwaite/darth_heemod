### takes a while to run, make faster by decreasing
#n_samp <- 500



###########################################################################
# adapted DARTH random search calibration 

###########################################################################

###################  Calibration Specifications  ###################

# Model: ATTR Markov Model
# Inputs to be calibrated: N1N2_i, N1D_i, N2N3_i, N2D_i, where i = (6, 12, 18) 
# Note: the first NAC landmark KM surves as the baseline distribution. THe second as month 6, and so on
# Targets: prop_N1, prop_N2, prop_N3

# Search method: Random search using Latin-Hypercube Sampling
# Goodness-of-fit measure: Sum of log-likelihoods

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

source("R/01_ATTR_model.R") # creates the function run_crs_markov()

library(heemod)
library(tidyverse)
library(ggplot2)
library(dplyr)
# Check that it works
v_params_test <- c(p_N1N2_6 = 0.2,
                   p_N1N2_12 = 0.01,
                   p_N1N2_18 = 0.2,
                   p_N1N2_24 = 0.2,
                   p_N1N2_30 = 0.2,
                   p_N1N2_36 = 0.2,
                   
                   p_N1D_6 = 0.05,
                   p_N1D_12 = 0.97,
                   p_N1D_18 = 0.05,
                   p_N1D_24 = 0.05,
                   p_N1D_30 = 0.05,
                   p_N1D_36 = 0.05,
                   
                   p_N2N3_6 = 0.15,
                   p_N2N3_12 = 0.15,
                   p_N2N3_18 = 0.15,
                   p_N2N3_24 = 0.15,
                   p_N2N3_30 = 0.15,
                   p_N2N3_36 = 0.15,
                   
                   p_N2D_6 = 0.05,
                   p_N2D_12 = 0.05,
                   p_N2D_18 = 0.05,
                   p_N2D_24 = 0.05,
                   p_N2D_30 = 0.05,
                   p_N2D_36 = 0.05
                   )


test_results <- run_ATTR_markov(v_params_test) # It works!

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

# number of random samples
n_samp <- 1000

# names and number of input parameters to be calibrated
v_param_names <- c("p_N1N2_6",
                   "p_N1N2_12",
                   "p_N1N2_18",
                   
                   "p_N1D_6", 
                   "p_N1D_12",
                   "p_N1D_18",
                   
                   "p_N2N3_6",
                   "p_N2N3_12",
                   "p_N2N3_18",
                  
                   "p_N2D_6", 
                   "p_N2D_12",
                   "p_N2D_18"
                   
)
  
n_param <- length(v_param_names)

# range on input search space
#lb <- c(p_Mets = 0.04, p_DieMets = 0.04) # lower bound
#ub <- c(p_Mets = 0.16, p_DieMets = 0.16) # upper bound

lb <- c(p_N1N2_6 = 0.04,
        p_N1N2_12 = 0.04,
        p_N1N2_18 = 0.04,
        
        
        p_N1D_6 = 0.04,
        p_N1D_12 = 0.04,
        p_N1D_18 = 0.04,
        
        
        p_N2N3_6 = 0.04,
        p_N2N3_12 = 0.04,
        p_N2N3_18 = 0.04,
        
        
        p_N2D_6 = 0.04,
        p_N2D_12 = 0.04,
        p_N2D_18 = 0.04
       
)

ub <- c(p_N1N2_6 = 0.2,
        p_N1N2_12 = 0.2,
        p_N1N2_18 = 0.2,
       
        
        p_N1D_6 = 0.2,
        p_N1D_12 = 0.2,
        p_N1D_18 = 0.2,
        
        
        p_N2N3_6 = 0.2,
        p_N2N3_12 = 0.2,
        p_N2N3_18 = 0.2,
       
        
        p_N2D_6 = 0.2,
        p_N2D_12 = 0.2,
        p_N2D_18 = 0.2
       
)



################
# number of calibration targets
v_target_names <- c("prop_N1", "prop_N2", "prop_N1" )
n_target <- length(v_target_names)


####################################################################
######  Calibrate!  ######
####################################################################
# record start time of calibration
t_init <- Sys.time()

###  Generate a random sample of input values  ###

# Sample unit Latin Hypercube
m_lhs_unit <- randomLHS(n_samp, n_param)


# Rescale to min/max of each parameter
m_param_samp <- matrix(nrow=n_samp,ncol=n_param)
for (i in 1:n_param){
  m_param_samp[,i] <- qunif(m_lhs_unit[,i],
                            min = lb[i],
                            max = ub[i])
}
colnames(m_param_samp) <- v_param_names

# view resulting parameter set samples
#pairs.panels(m_param_samp)


##################
# multiple targets
##################

# initialize goodness-of-fit vector
m_GOF <- matrix(nrow = n_samp, ncol = n_target)
colnames(m_GOF) <- paste0(v_target_names, "_fit")

# loop through sampled sets of input values
for (j in 1:n_samp){
  
  ###  Run model for a given parameter set  ###
  model_res <- run_ATTR_markov(v_params = m_param_samp[j, ])
  
  
  ###  Calculate goodness-of-fit of model outputs to targets  ###
  
  # TARGET 1: Survival ("Surv")
  # log likelihood  
  m_GOF[j,1] <- sum(dnorm(x = lst_targets$dist$N1$value,
                          mean = model_res$prop_N1_use,
                          sd = lst_targets$dist$N1$se,
                          log = T))
  
  m_GOF[j,2] <- sum(dnorm(x = lst_targets$dist$N2$value,
                          mean = model_res$prop_N2_use,
                          sd = lst_targets$dist$N2$se,
                          log = T))
  
  m_GOF[j,3] <- sum(dnorm(x = lst_targets$dist$N3$value,
                          mean = model_res$prop_N3_use,
                          sd = lst_targets$dist$N3$se,
                          log = T))
  
} # End loop over sampled parameter sets


###  Combine fits to the different targets into single GOF  ###
# can give different targets different weights
v_weights <- matrix(1, nrow = n_target, ncol = 1)
# matrix multiplication to calculate weight sum of each GOF matrix row
v_GOF_overall <- c(m_GOF%*%v_weights)
# Store in GOF matrix with column name "Overall"
m_GOF <- cbind(m_GOF,Overall_fit=v_GOF_overall)

# Calculate computation time
comp_time <- Sys.time() - t_init

####################################################################
######  Exploring best-fitting input sets  ######
####################################################################

# Arrange parameter sets in order of fit
m_calib_res <- cbind(m_param_samp,m_GOF)
m_calib_res <- m_calib_res[order(-m_calib_res[,"Overall_fit"]),]

# Examine the top 10 best-fitting sets
m_calib_res[1:10,]
m_calib_res[1:100,1]

# Plot the top 100 (top 10%)
plot(m_calib_res[1:100,1],m_calib_res[1:100,2],
     xlim=c(lb[1],ub[1]),ylim=c(lb[2],ub[2]),
     xlab = colnames(m_calib_res)[1],ylab = colnames(m_calib_res)[2])

# Pairwise comparison of top 100 sets
#pairs.panels(m_calib_res[1:100,v_param_names])

### Plot model-predicted output at best set vs targets ###
v_out_best <- run_ATTR_markov(m_calib_res[1,])


# TARGET 1: NAC 1 proportion
plotrix::plotCI(x = lst_targets$dist$N1$time, y = lst_targets$dist$N1$value, 
                ui = lst_targets$dist$N1$ub,
                li = lst_targets$dist$N1$lb,
                ylim = c(0, 1),
                xlim = c(0, 30),
                xlab = "Time", ylab = "Prop. N1")
points(x = lst_targets$dist$N1$time, 
       y = v_out_best$prop_N1_use, 
       pch = 8, col = "red")
legend("topright", 
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))




# TARGET 2: NAC 2 proportion
plotrix::plotCI(x = lst_targets$dist$N2$time, y = lst_targets$dist$N2$value, 
                ui = lst_targets$dist$N2$ub,
                li = lst_targets$dist$N2$lb,
                ylim = c(0, 1),
                xlim = c(0, 30),
                xlab = "Time", ylab = "Prop. N2")
points(x = lst_targets$dist$N2$time, 
       y = v_out_best$prop_N2_use, 
       pch = 8, col = "red")
legend("topright", 
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))


# TARGET 3: NAC 3 proportion
plotrix::plotCI(x = lst_targets$dist$N3$time, y = lst_targets$dist$N3$value, 
                ui = lst_targets$dist$N3$ub,
                li = lst_targets$dist$N3$lb,
                ylim = c(0, 1),
                xlim = c(0, 30),
                xlab = "Time", ylab = "Prop. N3")
points(x = lst_targets$dist$N3$time, 
       y = v_out_best$prop_N3_use, 
       pch = 8, col = "red")
legend("topright", 
       legend = c("Target", "Model-predicted output"),
       col = c("black", "red"), pch = c(1, 8))







































