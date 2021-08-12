rm(list = ls()) # to clean the workspace


library(heemod)
library(tidyverse)
library(ggplot2)

#-------------------
#PLACEHOLDERS
#....................

#FROM N1
p_N1N1_6 <- 0.6
p_N1N1_12 <- 0.01
p_N1N1_18 <- 0.6
p_N1N1_24 <- 0.6
p_N1N1_30 <- 0.6
p_N1N1_36 <- 0.6

p_N1N2_6 <- 0.2
p_N1N2_12 <- 0.01
p_N1N2_18 <- 0.2
p_N1N2_24 <- 0.2
p_N1N2_30 <- 0.2
p_N1N2_36 <- 0.2

p_N1N3_6 <- 0.15
p_N1N3_12 <- 0.01
p_N1N3_18 <- 0.15
p_N1N3_24 <- 0.15
p_N1N3_30 <- 0.15
p_N1N3_36 <- 0.15

p_N1D_6 <- 0.05
p_N1D_12 <- 0.97
p_N1D_18 <- 0.05
p_N1D_24 <- 0.05
p_N1D_30 <- 0.05
p_N1D_36 <- 0.05

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

p_N2N3_6 <- 0.15
p_N2N3_12 <- 0.15
p_N2N3_18 <- 0.15
p_N2N3_24 <- 0.15
p_N2N3_30 <- 0.15
p_N2N3_36 <- 0.15

p_N2D_6 <- 0.05
p_N2D_12 <- 0.05
p_N2D_18 <- 0.05
p_N2D_24 <- 0.05
p_N2D_30 <- 0.05
p_N2D_36 <- 0.05

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

p_N3D_6 <- 0.05
p_N3D_12 <- 0.05
p_N3D_18 <- 0.05
p_N3D_24 <- 0.05
p_N3D_30 <- 0.05
p_N3D_36 <- 0.05

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
          ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N1N1_24,
          p_N1N1_36
          )))))),

p_N1N2 = ifelse(markov_cycle <= 6 ,p_N1N2_6,
                ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N1N2_12,
                       ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N1N2_18,
                              ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N1N2_24,
                                     ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N1N2_30,
                                            ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N1N2_24,
                                                   p_N1N2_36
                                            )))))),

p_N1N3 = ifelse(markov_cycle <= 6 ,p_N1N3_6,
                ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N1N3_12,
                       ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N1N3_18,
                              ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N1N3_24,
                                     ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N1N3_30,
                                            ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N1N3_24,
                                                   p_N1N3_36
                                            )))))),

p_N1D = ifelse(markov_cycle <= 6 ,p_N1D_6,
               ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N1D_12,
                      ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N1D_18,
                             ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N1D_24,
                                    ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N1D_30,
                                           ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N1D_24,
                                                  p_N1D_36
                                           )))))),


#FROM N2
p_N2N1 = ifelse(markov_cycle <= 6 ,p_N2N1_6,
                ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N2N1_12,
                       ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N2N1_18,
                              ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N2N1_24,
                                     ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N2N1_30,
                                            ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N2N1_24,
                                                   p_N2N1_36
                                            )))))),

p_N2N2 = ifelse(markov_cycle <= 6 ,p_N2N2_6,
                ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N2N2_12,
                       ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N2N2_18,
                              ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N2N2_24,
                                     ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N2N2_30,
                                            ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N2N2_24,
                                                   p_N2N2_36
                                            )))))),

p_N2N3 = ifelse(markov_cycle <= 6 ,p_N2N3_6,
                ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N2N3_12,
                       ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N2N3_18,
                              ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N2N3_24,
                                     ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N2N3_30,
                                            ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N2N3_24,
                                                   p_N2N3_36
                                            )))))),

p_N2D = ifelse(markov_cycle <= 6 ,p_N2D_6,
               ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N2D_12,
                      ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N2D_18,
                             ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N2D_24,
                                    ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N2D_30,
                                           ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N2D_24,
                                                  p_N2D_36
                                           )))))),

#FROM N3
p_N3N1 = ifelse(markov_cycle <= 6 ,p_N3N1_6,
                ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N3N1_12,
                       ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N3N1_18,
                              ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N3N1_24,
                                     ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N3N1_30,
                                            ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N3N1_24,
                                                   p_N3N1_36
                                            )))))),

p_N3N2 = ifelse(markov_cycle <= 6 ,p_N3N2_6,
                ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N3N2_12,
                       ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N3N2_18,
                              ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N3N2_24,
                                     ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N3N2_30,
                                            ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N3N2_24,
                                                   p_N3N2_36
                                            )))))),

p_N3N3 = ifelse(markov_cycle <= 6 ,p_N3N3_6,
                ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N3N3_12,
                       ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N3N3_18,
                              ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N3N3_24,
                                     ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N3N3_30,
                                            ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N3N3_24,
                                                   p_N3N3_36
                                            )))))),

p_N3D = ifelse(markov_cycle <= 6 ,p_N3D_6,
               ifelse(markov_cycle >= 7 & markov_cycle <= 12, p_N3D_12,
                      ifelse(markov_cycle >= 13 & markov_cycle <= 18, p_N3D_18,
                             ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N3D_24,
                                    ifelse(markov_cycle >= 25 & markov_cycle <= 30, p_N3D_30,
                                           ifelse(markov_cycle >= 19 & markov_cycle <= 24, p_N3D_24,
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

#   C: The complement of other probabilities in a given row
mat_bsc <- define_transition(
  state_names = c(
    "NAC1",
    "NAC2",
    "NAC3",
    "Death"
  ),
  C, 0, 0,  p_N1D,
  C, 0, 0,  p_N2D,
  C, 0, 0, p_N3D,
  0,  0,      0,   1.00
)

#   C: The complement of other probabilities in a given row
mat_tafa <- define_transition(
  state_names = c(
    "NAC1",
    "NAC2",
    "NAC3",
    "Death"
  ),
  C, 0, 0,  p_N1D,
  C, 0, 0,  p_N2D,
  C, 0, 0, p_N3D,
  0,  0,      0,   1.00
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
  cycles = 20,
  cost = cost,
  effect = utility,
  method = "beginning",
  init = c(1000, 0, 0, 0)
)
summary(res_mod)
plot(res_mod)

plot(res_mod, type = "counts", panel = "by_state", free_y = TRUE) +
  theme_bw() +
  scale_color_brewer(
    name = "Strategy",
    palette = "Set1"
  )



#--------------------------------------------------------------------------------
#plot survival stuff
#................................................................................
#plot death state count
plot(res_mod, type = "counts", states = "Death", panel = "by_state", free_y = TRUE) +
  theme_bw() +
  scale_color_brewer(
    name = "Strategy",
    palette = "Set1"
  )

#get state counts
c_state <- as.data.frame(get_counts(res_mod))

#convert long table into wide table
c_state_wide <- pivot_wider(c_state, names_from = state_names, values_from = count) 


#create new columns: total_alive & proportion_alive
surv_prop <-c_state_wide %>%
  rowwise() %>%
  mutate(total_alive = sum(across(starts_with("NAC")), na.rm = T)) %>%
  mutate(proportion_alive = total_alive/rowSums(across(Death:total_alive)))%>% 
  select(markov_cycle, proportion_alive)

#plot surv curve
ggplot(data=surv_prop, aes(x=markov_cycle, proportion_alive, group=1)) +
  geom_line(color="red")+
  geom_point()


#---------------------------------------------------------------
#Calibration 
#...............................................................

get_counts(res_mod) %>% 
  dplyr::filter(markov_cycle == 12 & state_names == "Death")

#define a function to extract the values we want to change from the model and return them as a numeric vector 
extract_values <- function(x) {
  dplyr::filter(
    get_counts(x),
    markov_cycle == 12 & state_names == "Death"
  )$count
}
extract_values(res_mod)


calib_fn <- define_calibration_fn(
  type = "count",
  strategy_names = c("bsc", "tafamidis"),
  element_names = c("Death", "Death"),
  cycles = c(12, 12)
)
calib_fn(res_mod)


#now call calibrate_model(), and give the values we want to reach as target_values.
res_cal <- calibrate_model(
  res_mod,
  parameter_names = "p_N1D",
  fn_values = extract_values,
  target_values = 700,
  lower = 0, upper = 1
)


res_cal
############################################################
