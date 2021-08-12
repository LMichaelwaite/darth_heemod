#.strategy_names markov_cycle proportion_N3


target_dist_N1 <- data.frame(time = c(6, 12, 18, 24), 
                             value = c(0.46, 0.43, 0.38, 0.34),
                             ub = c(0.46, 0.43, 0.38, 0.34),
                             lb =  c(0.46, 0.43, 0.38, 0.34),
                             se = c(0.01, 0.01, 0.01, 0.01)
                             )

target_dist_N2 <- data.frame(time = c(6, 12, 18, 24), 
                             value = c(0.37, 0.34, 0.38, 0.37),
                             ub = c(0.37, 0.34, 0.38, 0.37),
                             lb = c(0.37, 0.34, 0.38, 0.37),
                             se = c(0.01, 0.01, 0.01, 0.01)
                             )

target_dist_N3 <- data.frame(time = c(6, 12, 18, 24), 
                             value = c(0.17, 0.23, 0.24, 0.29),
                             ub= c(0.17, 0.23, 0.24, 0.29),
                             lb = c(0.17, 0.23, 0.24, 0.29),
                             se = c(0.01, 0.01, 0.01, 0.01)
                             )

target_dist <- list(N1 = target_dist_N1,
                    N2 = target_dist_N2,
                    N3 = target_dist_N3)


#lst_targets$Surv$time, y = lst_targets$Surv$value


ATTR_targets <- list(dist = target_dist)
ATTR_targets 

save(ATTR_targets , file = "data/ATTR_CalibTargets.RData")