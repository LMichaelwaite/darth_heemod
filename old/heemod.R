library(heemod)

#-------------------
#PLACEHOLDERS
#....................

#FROM N1
p_N1N1_6 <- 0.5
p_N1N1_12 <- 0.01
p_N1N1_18 <- 0.5
p_N1N1_24 <- 0.5
p_N1N1_30 <- 0.5
p_N1N1_36 <- 0.5

p_N1N2_6 <- 0.3
p_N1N2_12 <- 0.01
p_N1N2_18 <- 0.3
p_N1N2_24 <- 0.3
p_N1N2_30 <- 0.3
p_N1N2_36 <- 0.3

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

par_mod <- define_parameters(
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
  cost_zido = 2278
)

#-----------------------------------------
#DEFINE TRANSITION
#.........................................



mat_bsc <- define_transition(
  p_N1N1, p_N1N2, p_N1N3,  p_N1D,
  p_N2N1, p_N2N2, p_N2N3,  p_N2D,
  p_N3N1, p_N3N2, p_N3N3, p_N3D,
       0,      0,      0,   1.00
)

mat_tafa <- define_transition(
  p_N1N1, p_N1N2, p_N1N3,  p_N1D,
  p_N2N1, p_N2N2, p_N2N3,  p_N2D,
  p_N3N1, p_N3N2, p_N3N3, p_N3D,
  0,      0,      0,   1.00
)

#-----------------------------------------
#DEFINE STATES
#.........................................


A_bsc <- define_state(
  cost_health = 2756,
  cost_drugs = cost_zido,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
  life_year = 1
)
B_bsc <- define_state(
  cost_health = 3052,
  cost_drugs = cost_zido,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
  life_year = 1
)
C_bsc <- define_state(
  cost_health = 9007,
  cost_drugs = cost_zido,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
  life_year = 1
)
D_bsc <- define_state(
  cost_health = 0,
  cost_drugs = 0,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
  life_year = 0
)

A_tafa <- define_state(
  cost_health = 2756,
  cost_drugs = cost_zido + cost_lami,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
  life_year = 1
)
B_tafa <- define_state(
  cost_health = 3052,
  cost_drugs = cost_zido + cost_lami,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
  life_year = 1
)
C_tafa <- define_state(
  cost_health = 9007,
  cost_drugs = cost_zido + cost_lami,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
  life_year = 1
)
D_tafa <- define_state(
  cost_health = 0,
  cost_drugs = 0,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
  life_year = 0
)

mod_bsc <- define_strategy(
  transition = mat_bsc,
  A_bsc,
  B_bsc,
  C_bsc,
  D_bsc
)
mod_tafa <- define_strategy(
  transition = mat_tafa,
  A_tafa,
  B_tafa,
  C_tafa,
  D_tafa
)

res_mod <- run_model(
  bsc = mod_bsc,
  tafa = mod_tafa,
  parameters = par_mod,
  cycles = 40,
  cost = cost_total,
  effect = life_year,
  method = "end",
  init = c(600, 200, 100, 0)
)


plot(res_mod)
summary(res_mod)

