
prop_NAC_0 <- data.frame(N1 = 0.46, N2 = 0.37, N3 = 0.17) #prportion in N stages at time 0
prop_NAC_1 <- data.frame(N1 = 0.43, N2 = 0.34, N3 = 0.23)
prop_NAC_2 <- data.frame(N1 = 0.38, N2= 0.38, N3 = 0.24)
prop_NAC_3 <- data.frame(N1 = 0.34, N2 =  0.37, N3 = 0.29)

p_N1D_0 <- 0.002755
p_N2D_0 <- 0.005372
p_N3D_0 <- 0.012067 

p_N1D_1 <- 0.00241
p_N2D_1 <- 0.0590
p_N3D_1 <- 0.1096

p_N1D_2 <- 0.103
p_N2D_2 <- 0.252
p_N3D_2 <- 0.441

v_D1 <- 0.005291 # proportion dead after timestop 1
v_D2 <- 0.055555
v_D3 <- 0.240213

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

##### 0 -6 #####
ofunc_1 <- function(y){
  y1 <- y[1]
  y2 <- y[2]
  # Manually construct your matrix here using y1, y2 where necessary)
  Q <- matrix(                      
    c(1- p_N1D_0 - y1, y1,          0, p_N1D_0, 
      0, 1 - p_N2D_0 - y2,         y2,  p_N2D_0, 
      0,                0, 1- p_N3D_0,  p_N3D_0,
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

solution_1

##### 6 - 12 #####
ofunc_2 <- function(y){
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
    x2_hat <- t(t(x1) %*% Q)
    d <- x2_hat - x2
   
    return(sum(d * d))
   
}

solution_2 <- optim(par = c(0.01, 0.01), fn = ofunc_2, lower = c(0.0, 0.0), upper = c( 1, 1),
                  method = "L-BFGS-B" )

solution_2

##### 12 - 18 #####

ofunc_3 <- function(y){
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
  x3_hat <- t(t(x2) %*% Q)
  d <- x3_hat - x3
  
  return(sum(d * d))
  
}

solution_3 <- optim(par = c(0.01, 0.01), fn = ofunc_2, lower = c(0.0, 0.0), upper = c( 1, 1),
                    method = "L-BFGS-B" )

solution_3

