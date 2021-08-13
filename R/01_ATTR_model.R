
library(heemod)
library(tidyverse)
library(ggplot2)
library(dplyr)

#### ATTR Markov model in a function ####
run_ATTR_markov <- function(v_params) {
  with(as.list(v_params), {
   
        #-------------------
        #PLACEHOLDERS
        #....................
        # HASHED transition intensities are for calibration
        #FROM N1
        p_N1N1_6 <- 0.6 
        p_N1N1_12 <- 0.01
        p_N1N1_18 <- 0.6
        p_N1N1_24 <- 0.6
        p_N1N1_30 <- 0.6
        p_N1N1_36 <- 0.6
        
        #p_N1N2_6 <- 0.2
        #p_N1N2_12 <- 0.01
        #p_N1N2_18 <- 0.2
        #p_N1N2_24 <- 0.2
        #p_N1N2_30 <- 0.2
        #p_N1N2_36 <- 0.2
        
        p_N1N3_6 <- 0.15
        p_N1N3_12 <- 0.01
        p_N1N3_18 <- 0.15
        p_N1N3_24 <- 0.15
        p_N1N3_30 <- 0.15
        p_N1N3_36 <- 0.15
        
        #p_N1D_6 <- 0.05
        #p_N1D_12 <- 0.97
        #p_N1D_18 <- 0.05
        #p_N1D_24 <- 0.05
        #p_N1D_30 <- 0.05
        #p_N1D_36 <- 0.05
        
        #FROM N2
        p_N2N1_6 <- 0.5
        p_N2N1_12 <- 0.5
        p_N2N1_18 <- 0.5
        p_N2N1_24 <- 0.5
        p_N2N1_30 <- 0.5
        p_N2N1_36 <- 0.5
        
        p_N2N2_6 <- 0.3
        p_N2N2_12 <- 0.3
        p_N2N2_18 <- 0.3
        p_N2N2_24 <- 0.3
        p_N2N2_30 <- 0.3
        p_N2N2_36 <- 0.3
        
        #p_N2N3_6 <- 0.15
        #p_N2N3_12 <- 0.15
        #p_N2N3_18 <- 0.15
        #p_N2N3_24 <- 0.15
        #p_N2N3_30 <- 0.15
        #p_N2N3_36 <- 0.15
        
        #p_N2D_6 <- 0.05
        #p_N2D_12 <- 0.05
        #p_N2D_18 <- 0.05
        #p_N2D_24 <- 0.05
        #p_N2D_30 <- 0.05
        #p_N2D_36 <- 0.05
        
        #FROM N3
        p_N3N1_6 <- 0.5
        p_N3N1_12 <- 0.5
        p_N3N1_18 <- 0.5
        p_N3N1_24 <- 0.5
        p_N3N1_30 <- 0.5
        p_N3N1_36 <- 0.5
        
        p_N3N2_6 <- 0.3
        p_N3N2_12 <- 0.3
        p_N3N2_18 <- 0.3
        p_N3N2_24 <- 0.3
        p_N3N2_30 <- 0.3
        p_N3N2_36 <- 0.3
        
        p_N3N3_6 <- 0.15
        p_N3N3_12 <- 0.15
        p_N3N3_18 <- 0.15
        p_N3N3_24 <- 0.15
        p_N3N3_30 <- 0.15
        p_N3N3_36 <- 0.15
        
        p_N3D_6 <- 0.15
        p_N3D_12 <- 0.15
        p_N3D_18 <- 0.15
        p_N3D_24 <- 0.15
        p_N3D_30 <- 0.15
        p_N3D_36 <- 0.15
        
        #-----------------------------------------
        #DEFINE PARAMETERS
        #.........................................
        
        param <- define_parameters(
          rr = ifelse(markov_cycle <= 2, .509, 1),
          #FROM N1
          p_N1N1 = ifelse(markov_cycle <= 6 ,p_N1N1_6,
                          ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N1N1_12,
                                 ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N1N1_18,
                                        ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N1N1_24,
                                               ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N1N1_30,
                                                      ifelse(markov_cycle >= 31 & markov_cycle <= 36, p_N1N1_24,
                                                             p_N1N1_36
                                                      )))))),
          
          p_N1N2 = ifelse(markov_cycle <= 6 ,p_N1N2_6,
                          ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N1N2_12,
                                 ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N1N2_18,
                                        p_N1N2_18))),                                         #If greater than 18, then constant at p_N1N2_18. since NAC data only goes to 18.
                                        
          
          p_N1N3 = ifelse(markov_cycle <= 6 ,p_N1N3_6,
                          ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N1N3_12,
                                 ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N1N3_18,
                                        ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N1N3_24,
                                               ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N1N3_30,
                                                      ifelse(markov_cycle >= 31 & markov_cycle <= 36, p_N1N3_24,
                                                             p_N1N3_36
                                                      )))))),
          
          p_N1D = ifelse(markov_cycle <= 6 ,p_N1D_6,
                         ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N1D_12,
                                ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N1D_18,
                                       p_N1D_18))),
                                       
          
          
          #FROM N2
          p_N2N1 = ifelse(markov_cycle <= 6 ,p_N2N1_6,
                          ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N2N1_12,
                                 ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N2N1_18,
                                        ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N2N1_24,
                                               ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N2N1_30,
                                                      ifelse(markov_cycle >= 31 & markov_cycle <= 36, p_N2N1_24,
                                                             p_N2N1_36
                                                      )))))),
          
          p_N2N2 = ifelse(markov_cycle <= 6 ,p_N2N2_6,
                          ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N2N2_12,
                                 ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N2N2_18,
                                        ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N2N2_24,
                                               ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N2N2_30,
                                                      ifelse(markov_cycle >= 31 & markov_cycle <= 36, p_N2N2_24,
                                                             p_N2N2_36
                                                      )))))),
          
          p_N2N3 = ifelse(markov_cycle <= 6 ,p_N2N3_6,
                          ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N2N3_12,
                                 ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N2N3_18,
                                        p_N2N3_18))),
                                        
          
          p_N2D = ifelse(markov_cycle <= 6 ,p_N2D_6,
                         ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N2D_12,
                                ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N2D_18,
                                       p_N2D_18))),
                                      
          
          #FROM N3
          p_N3N1 = ifelse(markov_cycle <= 6 ,p_N3N1_6,
                          ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N3N1_12,
                                 ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N3N1_18,
                                        ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N3N1_24,
                                               ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N3N1_30,
                                                      ifelse(markov_cycle >= 31 & markov_cycle <= 36, p_N3N1_24,
                                                             p_N3N1_36
                                                      )))))),
          
          p_N3N2 = ifelse(markov_cycle <= 6 ,p_N3N2_6,
                          ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N3N2_12,
                                 ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N3N2_18,
                                        ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N3N2_24,
                                               ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N3N2_30,
                                                      ifelse(markov_cycle >= 31 & markov_cycle <= 36, p_N3N2_24,
                                                             p_N3N2_36
                                                      )))))),
          
          p_N3N3 = ifelse(markov_cycle <= 6 ,p_N3N3_6,
                          ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N3N3_12,
                                 ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N3N3_18,
                                        ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N3N3_24,
                                               ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N3N3_30,
                                                      ifelse(markov_cycle >= 31 & markov_cycle <= 36, p_N3N3_24,
                                                             p_N3N3_36
                                                      )))))),
          
          p_N3D = ifelse(markov_cycle <= 6 ,p_N3D_6,
                         ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N3D_12,
                                ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N3D_18,
                                       ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N3D_24,
                                              ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N3D_30,
                                                     ifelse(markov_cycle >= 31 & markov_cycle <= 36, p_N3D_24,
                                                            p_N3D_36
                                                     )))))),  
          
          cost_lami = ifelse(markov_cycle <= 2, 2086.5, 0),
          cost_zido = 2278,
          u_successP = .85,
          u_revisionTHR = .30,
          u_successR = .75,
          c_revisionTHR = 5294
        )
        
        #-----------------------------------------
        #DEFINE TRANSITION
        #.........................................
        
     
        
        
        mat_bsc <- define_transition(
          state_names = c(
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
        
        mat_tafa <- define_transition(
          state_names = c(
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
        
        #-------------------------------------------------------
        #DEFINE strategy & state dependent costs and utilities
        #.......................................................
        mod_bsc <- define_strategy(
          transition = mat_bsc,
          NAC1 = define_state(
            utility = 0,
            cost = 394
          ),
          NAC2= define_state(
            utility = discount(u_successP, .015),
            cost = 0
          ),
          NAC3 = define_state(
            utility = discount(u_revisionTHR, .015),
            cost = discount(c_revisionTHR, .06)
          ),
          Death = define_state(
            utility = 0,
            cost = 0
          )
        )
        
        mod_tafamidis <- define_strategy(
          transition = mat_tafa,
          NAC1 = define_state(
            utility = 0,
            cost = 579
          ),
          NAC2= define_state(
            utility = discount(u_successP, .015),
            cost = 0
          ),
          NAC3 = define_state(
            utility = discount(u_revisionTHR, .015),
            cost = discount(c_revisionTHR, .06)
          ),
          Death = define_state(
            utility = 0,
            cost = 0
          )
        )
        
        res_mod <- run_model(
          bsc = mod_bsc,
          tafamidis = mod_tafamidis,
          parameters = param,
          cycles = 60,
          cost = cost,
          effect = utility,
          method = "beginning",
          init = c(436, 350, 159, 0)
        )
        summary(res_mod)
    out <-    plot(res_mod)
       
    
    
    
    ####### OUTPUT FOR CALIBRATION ###########################################
    #### NAC stage distribution ####
    
    #get state counts
    c_state <- as.data.frame(get_counts(res_mod))
    
    #convert long table into wide table
    c_state_wide <- pivot_wider(c_state, names_from = state_names, values_from = count) 
    
    #Filter for i) only month 6, 12, 18, 24; ii) control arm
    
    nac_dist <-c_state_wide %>%
      rowwise() %>%
      mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
      mutate(proportion_alive = total_alive/rowSums(across(Death:total_alive)))%>%
      mutate(proportion_N1= NAC1/total_alive)%>%
      mutate(proportion_N2= NAC2/total_alive)%>%
      mutate(proportion_N3= NAC3/total_alive)%>%
      filter(.strategy_names =="bsc")%>%
      select(.strategy_names,markov_cycle, proportion_alive, proportion_N1,proportion_N2, proportion_N3)
    
    nac_dist
    #create NAC dist. outs
    v_prop_N1 <- nac_dist %>% 
      select(markov_cycle, proportion_N1) %>%
      filter(markov_cycle %in% c(6, 12, 18))
    
    v_prop_N2 <- nac_dist %>% 
      select(markov_cycle, proportion_N2) %>%
      filter(markov_cycle %in% c(6, 12, 18))
    
    v_prop_N3 <- nac_dist %>% 
      select(markov_cycle, proportion_N3) %>%
      filter(markov_cycle %in% c(6, 12, 18))
    
  
    v_prop_N1_use <- as.numeric(unlist(v_prop_N1['proportion_N1']))
    v_prop_N2_use <- as.numeric(unlist(v_prop_N2['proportion_N2']))
    v_prop_N3_use <- as.numeric(unlist(v_prop_N3['proportion_N3']))
    

    
    ####### RETURN OUTPUT  ###########################################
    out <- list(prop_N1 = v_prop_N1, 
                prop_N2 = v_prop_N2, 
                prop_N3 = v_prop_N3,
                prop_N1_use = v_prop_N1_use, 
                prop_N2_use = v_prop_N2_use, 
                prop_N3_use = v_prop_N3_use
                 )

    
    return(out)
  }
  )
}