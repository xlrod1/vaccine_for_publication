# M monotonicity check function
compute_EY_am_no_S <- function(data, a, m){
  EY_am <- 0

  # Step 2: compute E[Y^{a,m}]
 
    for(u in 0:1){
      EY_amu <- mean(data$Y[data$A == a & data$B == m & data$U == u])
      P_u <- mean(data$U == u)
      EY_am <- EY_am + EY_amu * P_u 
    }
  
  
  return(EY_am)
}
