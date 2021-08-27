# estimates for p_N1N2 ; p_N2N3 
# 0 = diagnosis time
# time stop 1: 0 - 6
# time stop 2: 6 - 12
# time stop 3: 12 - 24


#### input parameters to calculate y1 and y2 ####
prop_NAC_0 <- data.frame(N1 = 0.4614, N2 = 0.3704, N3 = 0.16825) #prportion in N stages at time 0
prop_NAC_1 <- data.frame(N1 = 0.4306, N2 = 0.3403, N3 = 0.2292)
prop_NAC_2 <- data.frame(N1 = 0.3843, N2= 0.3754, N3 = 0.2402)


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

p_N1D_3 <- 0.05284029638 # month 12 - 24
p_N2D_3 <- 0.135174477
p_N3D_3 <- 0.2525989344

v_D1 <- 0.005291 # proportion dead after timestop 1
v_D2 <- 0.055555
v_D3 <- 0.240214


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


# 12 - 24 ACTUAL    

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


##### convert 12 month probability to 6 month probability ####

rescale_prob <- function (p, to = 1, from = 1) 
{
  stopifnot(p >= 0, p <= 1, to > 0, from > 0)
  r <- -log(1 - p)/from
  rate_to_prob(r, to = to)
}

p_N1N2_3 <- rescale_prob(p_N1N2_3, to = 0.5, from = 1)
p_N2N3_3 <- rescale_prob(p_N2N3_3, to = 0.5, from = 1)

summary_NiD <- data.frame(p_N1N2_1, p_N2N3_1,
                          p_N1N2_2, p_N2N3_2,
                          p_N1N2_3, p_N2N3_3)



