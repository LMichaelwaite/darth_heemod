
library(heemod)
library(tidyverse)
library(ggplot2)
library(dplyr)

#### ATTR Markov model in a function ####
run_ATTR_markov <- function(v_params) {
  with(as.list(v_params), {
   
        #-----------------------------------------
        #DEFINE PARAMETERS
        #.........................................
        
        param <- define_parameters(
          p_N1N2_1 = 0.002536125,
          p_N2N3_1 = 0.003078164,
          
          p_N1N2_2 = 0.0905894,
          p_N2N3_2 = 0.2091134,
          
          p_N1N2_3 = 0.1549904,
          p_N2N3_3 = 0.08997529,
          
          p_N1N2_4 = 0.1181116,
          p_N2N3_4 = 0.09979286,
          
          p_N1D_1 = 0.002755, # month 0 -6
          p_N2D_1 = 0.005372,
          p_N3D_1 = 0.012067, 
          
          p_N1D_2 = 0.00241 ,# month 6 - 12
          p_N2D_2 = 0.0590,
          p_N3D_2 = 0.1096,
          
          #p_N1D_3 = 0.103,
          #p_N2D_3 = 0.252,
          #p_N3D_3 = 0.441,
          
          p_N1D_3 = 0.05144 ,# month 12 - 18
          p_N2D_3 = 0.12604,
          p_N3D_3 = 0.22070,
          
          p_N1D_4 = 0.04965 ,# month 18 -24
          p_N2D_4 = 0.12163,
          p_N3D_4 = 0.21298,
          
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
          HR = 0.7,
          
         # HR = ifelse(markov_cycle <= 1 , 1,
           #           ifelse(markov_cycle >= 2 & markov_cycle < 3, 1,
           #                  ifelse(markov_cycle >= 3 & markov_cycle < 4, 0.7,
            #                        0.7))),
          
         ##### HOSPITALISATION ####
         # hospitalisation 
         
         #bsc cv hosp
         r_hosp = 0.70, # CV hosp rate per person per year
         p_hosp = 1- exp(-r_hosp*1/2), #CV hosp 6 month probability 
         
         #tafamidis cv hosp
         HR_hosp_tafa = 0.68,
         r_hosp_tafa = r_hosp * HR_hosp_tafa, # CV hosp rate per person per year
         p_hosp_tafa = 1- exp(-r_hosp_tafa*1/2), # CV hosp 6 month probability
         
         disutility_hosp = 0.10, # https://onlinelibrary.wiley.com/doi/full/10.1002/ehf2.12844
         
         ##### COSTS AND UTILITIES ####                                         
          #cost_lami = ifelse(markov_cycle <= 2, 2086.5, 0),
         # cost_zido = 2278,
         #divde by two due to 6 month cycle
          u_NAC_I = 0.893/2,
          u_NAC_II = 0.802/2,
          u_NAC_III = ((0.706 + 0.406)/2)/2, 
          u_NAC_I_tafa = 0.874/2,
          u_NAC_II_tafa = 0.832/2,
          u_NAC_III_tafa = ((0.707 + 0.558)/2)/2,
          c_tafa = 10840.82*6, # cost per 6 months; GBP; source: NICE tafamidis STA pp. 137
          c_bsc = 7031.92/2, # GBP; source nuffield trust estimate for health care sepnding on 75 year olds (UK) inflation adjusted to 2020 ((find actual source)) https://www.theguardian.com/society/2016/feb/01/ageing-britain-two-fifths-nhs-budget-spent-over-65s 
         #c_bsc = 18462/2, # background hc cost 6 months; dollars US Kazi
          c_hosp = 2536.88, # GBP source: NICE tafamidis STA pp.140
          c_eol = 9287.86 # # GBP source: NICE tafamidis STA pp.140, end of life costs
        )
        #-----------------------------------------
        #DEFINE TRANSITION
        #.........................................
        
        mat_bsc <- define_transition(
      
          1-(p_N1N2 + p_N1D),                 p_N1N2,           0,  p_N1D,
                           0,     1- (p_N2N3+ p_N2D),      p_N2N3,  p_N2D,
                           0,                      0,   1-(p_N3D),  p_N3D,
                           0,                      0,           0,   1.00
        )
        
        mat_tafa <- define_transition(
         
          1-(p_N1N2 + p_N1D*HR),                 p_N1N2,           0,  p_N1D*HR,
                           0,     1- (p_N2N3+ p_N2D*HR),      p_N2N3,  p_N2D*HR,
                           0,                      0,   1-(p_N3D*HR),  p_N3D*HR,
                           0,                      0,           0,   1.00
        )
        
        
        plot(mat_bsc)
        plot(mat_tafa)
        
        #-------------------------------------------------------
        #DEFINE strategy & state dependent costs and utilities
        #.......................................................
      NAC1 <- define_state(
          cost_health = 0,
          cost_drugs = dispatch_strategy(
            bsc = c_bsc,
            tafa = c_tafa + c_bsc
          ),
          cost_hosp = dispatch_strategy(
            bsc = p_hosp * c_hosp,
            tafa = p_hosp_tafa* c_hosp),
          
          cost_total = discount(cost_health + cost_drugs + cost_hosp, .03),
          
          qaly = dispatch_strategy(
            bsc = u_NAC_I,
            tafa = u_NAC_I_tafa
            ),
          cv_hosp_disutility = dispatch_strategy(
            bsc = p_hosp * disutility_hosp,
            tafa = p_hosp_tafa* disutility_hosp),
          
          qaly_total =  discount(qaly-cv_hosp_disutility, 0.0),
          life_year = 0.5
        )
      NAC2 <- define_state(
          cost_health = 0,
          cost_drugs = dispatch_strategy(
            bsc = c_bsc,
            tafa = c_tafa + c_bsc
          ),
          cost_hosp = dispatch_strategy(
            bsc = p_hosp * c_hosp,
            tafa = p_hosp_tafa* c_hosp),
          cost_total = discount(cost_health + cost_drugs + cost_hosp, .03),
          
          qaly = dispatch_strategy(
            bsc = u_NAC_II,
            tafa = u_NAC_II_tafa
          ),
          cv_hosp_disutility = dispatch_strategy(
            bsc = p_hosp * disutility_hosp,
            tafa = p_hosp_tafa* disutility_hosp),
          qaly_total =  discount(qaly- cv_hosp_disutility, 0.0),
          life_year = 0.5
        )
      NAC3 <- define_state(
          cost_health = 0,
          cost_drugs = dispatch_strategy(
            bsc = c_bsc,
            tafa = c_tafa + c_bsc
          ),
          cost_hosp = dispatch_strategy(
            bsc = p_hosp * c_hosp,
            tafa = p_hosp_tafa* c_hosp),
          cost_total = discount(cost_health + cost_drugs + cost_hosp, .03),
          
          qaly = dispatch_strategy(
            bsc = u_NAC_III,
            tafa = u_NAC_III_tafa
          ),
          cv_hosp_disutility = dispatch_strategy(
            bsc = p_hosp * disutility_hosp,
            tafa = p_hosp_tafa* disutility_hosp),
          qaly_total =  discount(qaly - cv_hosp_disutility, 0.0),
          life_year = 0.5
        )
      Death <- define_state(
        cost_health = 0,
        cost_drugs = 0,
        cost_hosp =0,
        cost_total = discount(cost_health + cost_drugs + cost_hosp, .03),
        qaly = 0,
        cv_hosp_disutility = 0,
        qaly_total =  discount(qaly - cv_hosp_disutility, 0.0),
        life_year = 0
      )
        
      
      strat_bsc <- define_strategy(
        transition = mat_bsc,
        NAC1,
        NAC2,
        NAC3,
        Death
      )      
      
      strat_tafa <- define_strategy(
        transition = mat_tafa,
        NAC1,
        NAC2,
        NAC3,
        Death
      )  
      
      res_mod <- run_model(
        bsc = strat_bsc,
        tafa = strat_tafa,
        parameters = param,
        cycles = 21,
        cost = cost_total,
        effect = qaly_total,
        method = "end",
        init = c(436, 350, 159, 0)
      )
      
      
      summary(res_mod)
      plot(res_mod)     
      #plot(res_mod, type = "ce") 
     
      #################################
      ##### PSA ####
     #################################
      
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
        
        p_N1D_1 ~ beta(1.201171875, 434.7988281),   # alpha = count of events; beta = number at risk - alpha 
        p_N2D_1 ~ beta(1.880274781, 348.1197252),
        p_N3D_1 ~ beta(1.918624122, 157.0813759),
        
        p_N1D_2 ~ beta(4.4792, 181.5208),
        p_N2D_2 ~ beta(8.6731, 138.3269),
        p_N3D_2 ~ beta(10.8477, 88.1523),
        
        p_N1D_3 ~ beta(11.112, 204.888),
        p_N2D_3 ~ beta(26.594, 184.406),
        p_N3D_3 ~ beta(29.794, 105.206),
        
        p_N1D_4 ~ beta(9.997, 191.373),
        p_N2D_4 ~ beta(25.703, 185.617),
        p_N3D_4 ~ beta(31.800, 117.510),
      
        c_bsc ~ gamma(mean = 7031.92/2, sd = sqrt((7031.92/2))),
        c_hosp ~ gamma(mean = 2546, sd= 507.38),
        
        HR ~ beta(12.02719, 4.639787), #from script "estimate beta parameters"
        HR_hosp_tafa ~ beta (36.16638, 16.63898)  #from script "estimate beta parameters"
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
      library(BCEA)
      bcea <- run_bcea(pm, plot = TRUE, Kmax = 1500000)
      
      
      
      
      
        
            
      
      
      
      
      
      
      
      
        
        
         mod_bsc <- define_strategy(
          transition = mat_bsc,
          NAC1 = define_state(
            qaly = u_NAC_I,
            cost = c_bsc
          ),
          NAC2= define_state(
            qaly = u_NAC_II,
            cost = c_bsc
          ),
          NAC3 = define_state(
            qaly = u_NAC_III,
            cost = c_bsc
          ),
          Death = define_state(
            qaly = 0,
            cost = 0
          )
        )
      
        mod_tafamidis <- define_strategy(
          transition = mat_tafa,
          NAC1 = define_state(
            qaly = u_NAC_I_tafa,
            cost = c_tafa
          ),
          NAC2= define_state(
            qaly = u_NAC_II_tafa,
            cost = c_tafa
          ),
          NAC3 = define_state(
            qaly = u_NAC_III_tafa,
            cost = c_tafa
          ),
          Death = define_state(
            qaly = 0,
            cost = 0
          )
        )
        
        res_mod <- run_model(
          bsc = mod_bsc,
          tafamidis = mod_tafamidis,
          parameters = param,
          cycles = 21,
          cost = cost,
          effect = qaly,
          method = "end",
          init = c(436, 350, 159, 0)
        )
        c_n <- sum(c(436, 350, 159, 0))
        
        summary(res_mod)
        head(get_values(res_mod ))
        
        
        plot(res_mod, type = "values", panel = "by_value",
             free_y = TRUE) +
          xlab("Time") +
          theme_bw() +
          scale_color_brewer(
            name = "Strategy",
            palette = "Set1"
          )

        
        ##### PSA ####
        
beta_params <- function (mean, sigma) 
{
  alpha <- ((1 - mean)/sigma^2 - 1/mean) * mean^2
  beta <- alpha * (1/mean - 1)
  params<- list(alpha = alpha, beta = beta)
  return(params)
}    
beta_u_NAC_I <-beta_params(mean = .893/2, sigma = .02)
beta_u_NAC_II <-beta_params(mean = .802/2, sigma = .01)
beta_u_NAC_III <-beta_params(mean = ((0.706+0.406)/2)/2, sigma = .04)

beta_u_NAC_I_tafa <-beta_params(mean = .874/2, sigma = .01)
beta_u_NAC_II_tafa <-beta_params(mean = .832/2, sigma = .01)
beta_u_NAC_III_tafa <-beta_params(mean = ((0.707+0.558)/2)/2, sigma = .04) 
      
        rsp <- define_psa(
           u_NAC_I ~ beta(beta_u_NAC_I$alpha, beta_u_NAC_I$beta),
           u_NAC_II ~ beta(beta_u_NAC_II$alpha, beta_u_NAC_II$beta),
           u_NAC_III ~ beta(beta_u_NAC_III$alpha, beta_u_NAC_III$beta),
           u_NAC_I_tafa ~ beta(beta_u_NAC_I_tafa$alpha, beta_u_NAC_I_tafa$beta),
           u_NAC_II_tafa ~ beta(beta_u_NAC_II_tafa$alpha, beta_u_NAC_II_tafa$beta),
           u_NAC_III_tafa ~ beta(beta_u_NAC_III_tafa$alpha, beta_u_NAC_III_tafa$beta),
           c_tafa ~ gamma(mean = 225000/2, sd = sqrt((225000/2))),
           c_bsc ~ gamma(mean = 18462/2, sd = sqrt((18462/2))))
           
        #"run_psa" is masked from dampack, must specify heemod
        pm <- heemod::run_psa( 
          model = res_mod,
          psa = rsp,
          N = 100
        )   
        
        summary(pm)
        plot(pm, type = "ce")
        plot(pm, type = "ac", max_wtp = 100000, log_scale = FALSE)        
        
        
        
        
        
        
        
        
        
# Run1: NAT baseline: init = c(436, 350, 159, 0)       
# Run2: NAT baseline: init = c(945, 0, 0, 0)       
    
    #### COMPUTE AVG QALY AND COST ####
        #get state counts
        c_state <- as.data.frame(get_counts(res_mod))
        
        #convert long table into wide table
        c_state_wide <- pivot_wider(c_state, names_from = state_names, values_from = count) 
        
        
        
        #BSC
        c_QALYS <- as.data.frame(get_values(res_mod))
        
        QALYS_bsc <-   c_QALYS %>%
          filter(.strategy_names =="bsc")  
        
        QALYS_bsc_wide <- pivot_wider(QALYS_bsc, names_from = value_names, values_from = value)   
        
        c_alive_bsc <-c_state_wide %>%
          rowwise() %>%
          mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
          filter(.strategy_names =="bsc")%>%
          select(total_alive)
        
        result_summary_bsc <- cbind(QALYS_bsc_wide, c_alive_bsc)
        
        adjust_qaly_bsc <- result_summary_bsc %>%
          rowwise() %>% 
          mutate(count_cv_hosp = total_alive*p_hosp)%>%
          mutate(disutility_cv_hosp = count_cv_hosp*disutility_hosp)%>%
          mutate(adj_qaly_bsc = qaly - disutility_cv_hosp) %>%
          mutate(adj_cost_bsc = cost + (count_cv_hosp*c_hosp ))
        
        avg_QALY_bsc <- sum(adjust_qaly_bsc$qaly)/c_n     
        avg_cost_bsc <- sum(adjust_qaly_bsc$cost)/c_n
        avg_QALY_adj_bsc <- sum(adjust_qaly_bsc$adj_qaly_bsc)/c_n
        avg_cost_adj_bsc <- sum(adjust_qaly_bsc$adj_cost_bsc)/c_n
        
        #TAFA 
        c_QALYS <- as.data.frame(get_values(res_mod))
        
        QALYS_tafa <-   c_QALYS %>%
          filter(.strategy_names =="tafamidis")  
      
        QALYS_tafa_wide <-  pivot_wider(QALYS_tafa, names_from = value_names, values_from = value)   
      
        c_alive_tafamidis <-c_state_wide %>%
          rowwise() %>%
          mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
          filter(.strategy_names =="tafamidis")%>%
          select(total_alive)
        
        result_summary_tafamidis <- cbind(QALYS_tafa_wide, c_alive_tafamidis)
        
        adjust_qaly_tafa <- result_summary_tafamidis %>%
          rowwise() %>% 
          mutate(count_cv_hosp = total_alive*p_hosp)%>%
          mutate(disutility_cv_hosp = count_cv_hosp*disutility_hosp)%>%
          mutate(adj_qaly_tafamidis = qaly - disutility_cv_hosp)%>%
          mutate(adj_cost_tafamidis = cost + (count_cv_hosp*c_hosp ))
        
        avg_QALY_tafa <- sum(adjust_qaly_tafa$qaly)/c_n     
        avg_cost_tafa <- sum(adjust_qaly_tafa$cost)/c_n 
        avg_QALY_adj_tafa <- sum(adjust_qaly_tafa$adj_qaly_tafamidis)/c_n
        avg_cost_adj_tafa <- sum(adjust_qaly_tafa$adj_cost_tafamidis)/c_n
      
        result_summary <- data.frame(QALY_tafa =  avg_QALY_tafa, adj_QALY_tafa = avg_QALY_adj_tafa,
                                     QALY_bsc =  avg_QALY_bsc, adj_QALY_bsc = avg_QALY_adj_bsc,
                                     cost_tafa =  avg_cost_tafa, adj_cost_tafa = avg_cost_adj_tafa,
                                     cost_bsc = avg_cost_bsc, adj_cost_bsc = avg_cost_adj_bsc)
        
        avg_QALY_tafa -  avg_QALY_bsc
        avg_QALY_adj_tafa - avg_QALY_adj_bsc
        
        ICER <- (avg_cost_tafa- avg_cost_bsc)/ ( avg_QALY_tafa - avg_QALY_bsc )
        ICER_2 <- (avg_cost_tafa- avg_cost_bsc)/ ( avg_QALY_adj_tafa - avg_QALY_adj_bsc )
        ICER_3 <- (avg_cost_adj_tafa- avg_cost_adj_bsc)/ ( avg_QALY_adj_tafa - avg_QALY_adj_bsc )
        
       
       
  #### adjust QALYs for CV hosp ####
      #bsc
      c_alive_bsc <-c_state_wide %>%
          rowwise() %>%
          mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
          filter(.strategy_names =="bsc")%>%
          select(total_alive)
        
        
       result_summary_bsc <- cbind(QALYS_bsc_wide, c_alive_bsc)
        
        adjust_qaly_bsc <- result_summary_bsc %>%
        rowwise() %>% 
        mutate(count_cv_hosp = total_alive*p_hosp)%>%
        mutate(disutility_cv_hosp = count_cv_hosp*disutility_hosp)%>%
        mutate(adj_qaly_bsc = qaly - disutility_cv_hosp)
        
      #tafamidis
        c_alive_tafamidis <-c_state_wide %>%
          rowwise() %>%
          mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
          filter(.strategy_names =="tafamidis")%>%
          select(total_alive)
        
        result_summary_tafamidis <- cbind(QALYS_tafa_wide, c_alive_tafamidis)
        
        adjust_qaly_tafa <- result_summary_tafamidis %>%
          rowwise() %>% 
          mutate(count_cv_hosp = total_alive*p_hosp)%>%
          mutate(disutility_cv_hosp = count_cv_hosp*disutility_hosp)%>%
          mutate(adj_qaly_tafamidis = qaly - disutility_cv_hosp)
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
         
        
        
        
        ####### OUTPUT FOR CALIBRATION ###########################################
    #### NAC stage distribution ####
    

         #get state counts
    c_state <- as.data.frame(get_counts(res_mod))
    
    #convert long table into wide table
    c_state_wide <- pivot_wider(c_state, names_from = state_names, values_from = count) 
    
    #Filter for i) only month 6, 12, 18, 24; ii) control arm
    # life years 
   
     #bsc
     life_years_bsc <-c_state_wide %>%
        rowwise() %>%
        mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
        mutate(proportion_alive = total_alive/rowSums(across(Death:total_alive)))%>%
        mutate(proportion_N1= NAC1/total_alive)%>%
        mutate(proportion_N2= NAC2/total_alive)%>%
        mutate(proportion_N3= NAC3/total_alive)%>%
        filter(.strategy_names =="bsc")%>%
        select(.strategy_names,markov_cycle, total_alive)
      
     avg_life_years_bsc <- sum(life_years_bsc$total_alive/2)/c_n
     avg_life_years_bsc
    
     # tafa
     life_years_tafamidis <-c_state_wide %>%
       rowwise() %>%
       mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
       mutate(proportion_alive = total_alive/rowSums(across(Death:total_alive)))%>%
       mutate(proportion_N1= NAC1/total_alive)%>%
       mutate(proportion_N2= NAC2/total_alive)%>%
       mutate(proportion_N3= NAC3/total_alive)%>%
       filter(.strategy_names =="tafamidis")%>%
       select(.strategy_names,markov_cycle, total_alive)
     
     avg_life_years_tafamidis <- sum(life_years_tafamidis$total_alive/2)/c_n
     avg_life_years_tafamidis
     
     
     
       
    
    # state dist
    
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