

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








##################################################################










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