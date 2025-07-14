# M monotonicity check function
compute_EY_am <- function(data, a, m){
  EY_am <- 0
  
  # Step 1: compute P*(S = s | A = a)
  P_star_s_given_a <- c()  # vector of length 2: P*(S=0|A=a), P*(S=1|A=a)
  for(s in 0:1){
    sum_prob <- 0
    for(u in 0:1){
      P_s_given_au <- mean(data$S[data$A == a & data$U == u] == s)
      P_u <- mean(data$U == u)
      sum_prob <- sum_prob + P_s_given_au * P_u
    }
    P_star_s_given_a[s + 1] <- sum_prob  # store in index 1 for s=0, index 2 for s=1
  }
  
  # Step 2: compute E[Y^{a,m}]
  for(s in 0:1){
    for(u in 0:1){
      EY_asmu <- mean(data$Y[data$A == a & data$B == m & data$S == s & data$U == u])
      P_u <- mean(data$U == u)
      EY_am <- EY_am + EY_asmu * P_u * P_star_s_given_a[s + 1]
    }
  }
  
  return(EY_am)
}
