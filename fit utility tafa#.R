a <- array(c(0.461, 0.37, 0.168, 0.487, 0.326, 0.237, 0.5919, 0.32, 0.0875), dim=c(3, 3))
u <- c(0.788, 0.744, 0.57)
a
b



#### interval (0,1) yrs

u_0 <- c(0.788)

utility_0 <- function(y){
  y1 <- y[1]
  y2 <- y[2]
  y3 <- y[3]
  # Manually construct your matrix here using y1, y2 where necessary)
  a <- array(c(0.461, 0.37, 0.168), dim=c(1, 3))
  Q <- c(y1, y2, y3)
  
  u_hat <- a %*% Q 
  
  d <- u_hat - u_0
  
  return(sum(d * d))
  
}

solution_u_0 <- optim(par = c(0.001, 0.001, 0.001), fn = utility_0, lower = c(0, 0, 0), upper = c( 0.99, 0.9, 0.9),
                      method = "L-BFGS-B" )
solution_u_0


####
#### year 1 

u_1 <- c(0.744)

utility_1 <- function(y){
  y1 <- y[1]
  y2 <- y[2]
  y3 <- y[3]
  # Manually construct your matrix here using y1, y2 where necessary)
  a <- array(c(0.437, 0.326, 0.237), dim=c(1, 3))
  Q <- c(y1, y2, y3)
  
  u_hat <- a %*% Q 
  
  d <- u_hat - u_1
  
  return(sum(d * d))
  
}

solution_u_1 <- optim(par = c(0.01, 0.01, 0.01), fn = utility_1, lower = c(0, 0, 0), upper = c( 1, 0.7720690, 0.3511070),
                      method = "L-BFGS-B" )
solution_u_1



###############
#5
##########

u_5 <- c(0.57)

utility_5 <- function(y){
  y1 <- y[1]
  y2 <- y[2]
  y3 <- y[3]
  # Manually construct your matrix here using y1, y2 where necessary)
  a <- array(c(0.273, 0.385, 0.329), dim=c(1, 3))
  Q <- c(y1, y2, y3)
  
  u_hat <- a %*% Q 
  
  d <- u_hat - u_5
  
  return(sum(d * d))
  
}

solution_u_5 <- optim(par = c(0.01, 0.01, 0.01), fn = utility_5, lower = c(0.6, 0, 0), upper = c( 0.9617104, 1, 0.3511070),
                      method = "L-BFGS-B" )
solution_u_5



###############
#10
##########

u_10 <- c(0.439)

utility_10 <- function(y){
  y1 <- y[1]
  y2 <- y[2]
  y3 <- y[3]
  # Manually construct your matrix here using y1, y2 where necessary)
  a <- array(c(0.1984834, 0.3559003,0.4456163), dim=c(1, 3))
  Q <- c(y1, y2, y3)
  
  u_hat <- a %*% Q 
  
  d <- u_hat - u_10
  
  return(sum(d * d))
  
}

solution_u_10 <- optim(par = c(0.01, 0.01, 0.01), fn = utility_10, lower = c(0, 0, 0), upper = c( 0.9617104, 0.4985425, 0.3511070),
                       method = "L-BFGS-B" )
solution_u_10




results <- list(solution_u_0$par, solution_u_1$par, solution_u_5$par, solution_u_10$par )
results