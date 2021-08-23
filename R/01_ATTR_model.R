
library(heemod)
library(tidyverse)
library(ggplot2)
library(dplyr)

#### ATTR Markov model in a function ####
run_ATTR_markov <- function(v_params) {
  with(as.list(v_params), {
   # 0 = diagnosis time
    # time stop 1: 0 - 6
    # time stop 2: 6 - 12
    # time stop 3: 12 - 18
    # time stop 4: 18 - 24
    
    #### input parameters to calculate y1 and y2 ####
    prop_NAC_0 <- data.frame(N1 = 0.4614, N2 = 0.3704, N3 = 0.16825) #prportion in N stages at time 0
    prop_NAC_1 <- data.frame(N1 = 0.4306, N2 = 0.3403, N3 = 0.2292)
    prop_NAC_2 <- data.frame(N1 = 0.3843, N2= 0.3754, N3 = 0.2402)
    prop_NAC_3 <- data.frame(N1 = 0.3583, N2 =  0.3760, N3 = 0.2657)
    prop_NAC_4 <- data.frame(N1 = 0.3323, N2 =  0.3766, N3 = 0.2911)
    
    #prop_NAC_3 <- data.frame(N1 = 0.4193, N2 =  0.3740, N3 = 0.2068)
    #prop_NAC_4 <- data.frame(N1 = 0.3323, N2 =  0.3766, N3 = 0.2911)
    
    p_N1D_1 <- 0.002755 # month 0 -6
    p_N2D_1 <- 0.005372
    p_N3D_1 <- 0.012067 
    
    p_N1D_2 <- 0.00241 # month 6 - 12
    p_N2D_2 <- 0.0590
    p_N3D_2 <- 0.1096
    
    #p_N1D_3 <- 0.103
    #p_N2D_3 <- 0.252
    #p_N3D_3 <- 0.441
    
    p_N1D_3 <- 0.05284029638 # month 12 - 18
    p_N2D_3 <- 0.135174477
    p_N3D_3 <- 0.2525989344
    
    p_N1D_4 <- 0.04965 # month 18 -24
    p_N2D_4 <- 0.12163
    p_N3D_4 <- 0.21298
    
    v_D1 <- 0.005291 # proportion dead after timestop 1
    v_D2 <- 0.055555
    v_D3 <- 0.131737
    v_D4 <- 0.124935
    
    x0 <- c(prop_NAC_0$N1, prop_NAC_0$N2, prop_NAC_0$N3, 0)    # known values
    
    x1 <- c((1-v_D1)*prop_NAC_0$N1, 
            (1-v_D1)*prop_NAC_0$N2,
            (1-v_D1)*prop_NAC_0$N3,
            v_D1
    )
    x2 <- c((1-v_D2)*prop_NAC_1$N1, 
            (1-v_D2)*prop_NAC_1$N2,
            (1-v_D2)*prop_NAC_1$N3,
            v_D2
    )
    x3 <- c((1-v_D3)*prop_NAC_2$N1, 
            (1-v_D3)*prop_NAC_2$N2,
            (1-v_D3)*prop_NAC_2$N3,
            v_D3
    )
    x4 <- c((1-v_D4)*prop_NAC_3$N1, 
            (1-v_D4)*prop_NAC_3$N2,
            (1-v_D4)*prop_NAC_3$N3,
            v_D4
    )
    
    
    #### calculate y1 and y2 ####
    #0 -6 
    ofunc_1 <- function(y){
      y1 <- y[1]
      y2 <- y[2]
      # Manually construct your matrix here using y1, y2 where necessary)
      Q <- matrix(                      
        c(1- p_N1D_1 - y1, y1,          0, p_N1D_1, 
          0, 1 - p_N2D_1 - y2,         y2,  p_N2D_1, 
          0,                0, 1- p_N3D_1,  p_N3D_1,
          0,                0,          0, 1.00),
        ncol =  4,
        byrow = TRUE     
      )
      x1_hat <- t(t(x0) %*% Q)
      d <- x1_hat - x1
      
      return(sum(d * d))
      
    }
    
    solution_1 <- optim(par = c(0.01, 0.01), fn = ofunc_1, lower = c(0, 0), upper = c( 1, 1),
                        method = "L-BFGS-B" )
    
    p_N1N2_1 <- solution_1$par[1]
    p_N2N3_1 <- solution_1$par[2]
    
    # 6 - 12
    ofunc_2 <- function(y){
      y1 <- y[1]
      y2 <- y[2]
      # Manually construct your matrix here using y1, y2 where necessary)
      Q <- matrix(                      
        c(1- p_N1D_2 - y1, y1,          0, p_N1D_2, 
          0, 1 - p_N2D_1 - y2,         y2,  p_N2D_2, 
          0,                0, 1- p_N3D_2,  p_N3D_2,
          0,                0,          0, 1.00),
        ncol =  4,
        byrow = TRUE     
      )
      x2_hat <- t(t(x1) %*% Q)
      d <- x2_hat - x2
      
      return(sum(d * d))
      
    }
    
    solution_2 <- optim(par = c(0.01, 0.01), fn = ofunc_2, lower = c(0.0, 0.0), upper = c( 1, 1),
                        method = "L-BFGS-B" )
    
    p_N1N2_2 <- solution_2$par[1]
    p_N2N3_2 <- solution_2$par[2]
    
    
    # 12 -18 ACTUAL    
    
    ofunc_3 <- function(y){
      y1 <- y[1]
      y2 <- y[2]
      # Manually construct your matrix here using y1, y2 where necessary)
      Q <- matrix(                      
        c(1- p_N1D_3 - y1, y1,          0, p_N1D_3, 
          0, 1 - p_N2D_3 - y2,         y2,  p_N2D_3, 
          0,                0, 1- p_N3D_3,  p_N3D_3,
          0,                0,          0, 1.00),
        ncol =  4,
        byrow = TRUE     
      )
      x3_hat <- t(t(x2) %*% Q)
      d <- x3_hat - x3
      
      return(sum(d * d))
      
    }
    
    solution_3 <- optim(par = c(0.01, 0.01), fn = ofunc_3, lower = c(0.0, 0.0), upper = c( 1, 1),
                        method = "L-BFGS-B" )
    
    p_N1N2_3 <- solution_3$par[1]
    p_N2N3_3 <- solution_3$par[2]
    
    # 18 - 24
    
    ofunc_4 <- function(y){
      y1 <- y[1]
      y2 <- y[2]
      # Manually construct your matrix here using y1, y2 where necessary)
      Q <- matrix(                      
        c(1- p_N1D_4 - y1, y1,          0, p_N1D_4, 
          0, 1 - p_N2D_4 - y2,         y2,  p_N2D_4, 
          0,                0, 1- p_N3D_4,  p_N3D_4,
          0,                0,          0, 1.00),
        ncol =  4,
        byrow = TRUE     
      )
      x4_hat <- t(t(x3) %*% Q)
      d <- x4_hat - x4
      
      return(sum(d * d))
      
    }
    
    solution_4 <- optim(par = c(0.01, 0.01), fn = ofunc_4, lower = c(0.0, 0.0), upper = c( 1, 1),
                        method = "L-BFGS-B" )
    
    p_N1N2_4 <- solution_4$par[1]
    p_N2N3_4 <- solution_4$par[2]
    
        #-----------------------------------------
        #DEFINE PARAMETERS
        #.........................................
        
        param <- define_parameters(
          rr = ifelse(markov_cycle <= 2, .509, 1),
          
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
          
          
          #### TREATMENT EFFECT ####
          HR = ifelse(markov_cycle <= 1 , 1,
               ifelse(markov_cycle >= 2 & markov_cycle < 3, 1,
               ifelse(markov_cycle >= 3 & markov_cycle < 4, 0.7,
                0.7))),
          
          
          ##### COSTS AND UTILITIES ####                                         
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
          1-(p_N1N2 + p_N1D*HR),                 p_N1N2,           0,  p_N1D*HR,
                           0,     1- (p_N2N3+ p_N2D*HR),      p_N2N3,  p_N2D*HR,
                           0,                      0,   1-(p_N3D*HR),  p_N3D*HR,
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
          cycles = 21,
          cost = cost,
          effect = utility,
          method = "beginning",
          init = c(436, 350, 159, 0)
        )
        summary(res_mod)

# Run1: NAT baseline: init = c(436, 350, 159, 0)       
# Run2: NAT baseline: init = c(945, 0, 0, 0)       
    
    
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
    
    ##### for tafa
    
    nac_dist_tafa <-c_state_wide %>%
      rowwise() %>%
      mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
      mutate(proportion_alive = total_alive/rowSums(across(Death:total_alive)))%>%
      mutate(proportion_N1= NAC1/total_alive)%>%
      mutate(proportion_N2= NAC2/total_alive)%>%
      mutate(proportion_N3= NAC3/total_alive)%>%
      filter(.strategy_names =="tafamidis")%>%
      select(.strategy_names,markov_cycle, proportion_alive, proportion_N1,proportion_N2, proportion_N3)
    
    nac_dist_tafa
    
    
    
    
    
    
    
    #create NAC dist. outs
    v_prop_N1 <- nac_dist %>% 
      select(markov_cycle, proportion_N1) %>%
      filter(markov_cycle %in% c(2, 3, 4, 5))
    
    v_prop_N2 <- nac_dist %>% 
      select(markov_cycle, proportion_N2) %>%
      filter(markov_cycle %in% c(2, 3, 4, 5))
    
    v_prop_N3 <- nac_dist %>% 
      select(markov_cycle, proportion_N3) %>%
      filter(markov_cycle %in% c(2, 3, 4, 5))
    
  
    v_prop_N1_use <- as.numeric(unlist(v_prop_N1['proportion_N1']))
    v_prop_N2_use <- as.numeric(unlist(v_prop_N2['proportion_N2']))
    v_prop_N3_use <- as.numeric(unlist(v_prop_N3['proportion_N3']))
    
### SURV target
 
    v_surv <- nac_dist %>% 
      select(markov_cycle, proportion_alive)
  
    v_surv_use <- as.numeric(unlist(v_surv['proportion_alive']))

#### plot ####   
    plot(res_mod, type = "counts", panel = "by_state", free_y = TRUE) +
      theme_bw() +
      scale_color_brewer(
        name = "Strategy",
        palette = "Set1"
      )
    
    plot(res_mod, type = "counts", panel = "by_strategy") +
      xlab("Time") +
      theme_bw() +
      scale_color_brewer(
        name = "State",
        palette = "Set1"
      )
    
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
      mutate(total_alive_tafamidis = sum(across(c(NAC1_tafamidis, NAC2_tafamidis, NAC3_tafamidis)), na.rm = T)) %>%
      mutate(total_dead_tafamidis = Death_tafamidis) %>%
      mutate(proportion_alive = total_alive_tafamidis/rowSums(across(Death_tafamidis:total_alive_tafamidis)))%>%
      select(markov_cycle, proportion_alive)
      surv_prop
   
    surv_prop_bsc <- c_state_wide %>%
      rowwise() %>%
      mutate(total_alive_bsc = sum(across(c(NAC1_bsc, NAC2_bsc, NAC3_bsc)), na.rm = T)) %>%
      mutate(proportion_alive = total_alive_bsc/rowSums(across(Death_bsc:total_alive_bsc))) %>%
      mutate(total_dead_bsc = Death_bsc) %>%  
      select(markov_cycle, proportion_alive, total_alive_bsc, total_dead_bsc)
      surv_prop_bsc 
    
    surv_prop
    #plot surv_prop curve
    ggplot() + 
      geom_line(data=surv_prop, aes(x=markov_cycle, y=proportion_alive), color='green') + 
      geom_line(data=surv_prop_bsc, aes(x=markov_cycle, y=proportion_alive), color='red')
    
  
#### tafamidis survival ####
    #convert long table into wide table
    c_state_wide_tafamidis <- pivot_wider(c_state_tafamidis, names_from = state_names, values_from = count) 
    
    #create new columns: total_alive & proportion_alive
    surv_prop_tafamidis <-c_state_wide_tafamidis %>%
      rowwise() %>%
      mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
      mutate(proportion_alive = total_alive/rowSums(across(Death:total_alive)))%>%
      mutate(proportion_alive = prop_alive_bsc$surv_prop) %>%
      select(markov_cycle, proportion_alive)
    
    #plot surv_prop curve
    ggplot(data=surv_prop_tafamidis, aes(x=markov_cycle, proportion_alive, group=1)) +
      geom_line(color="red")+
      geom_point()
    
    
    
    
    
    ############################################################
    
   
    
    
    ####### RETURN OUTPUT  ###########################################
    #prop_N1 is a DF with cols: cycle number, prop. N1, prop. N2, prop. N3. used to generate graphs
    #prop_N3_use is data type "double" used as the model output to assess GOF in calibration
    
    out <- list(prop_N1 = v_prop_N1, 
                prop_N2 = v_prop_N2, 
                prop_N3 = v_prop_N3,
                prop_N1_use = v_prop_N1_use, 
                prop_N2_use = v_prop_N2_use, 
                prop_N3_use = v_prop_N3_use,
                surv = v_surv,
                surv_use = v_surv_use
                 )

    
    return(out)
  }
  )
}