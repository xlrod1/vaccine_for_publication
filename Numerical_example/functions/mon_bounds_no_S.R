compute_monotonicity_bounds_no_S <- function(df) {
  # Define basic probabilities
  p_110 <- mean(df$Y == 1 & df$A == 1 & df$B == 0) / mean(df$A == 1 & df$B == 0)
  p_111 <- mean(df$Y == 1 & df$A == 1 & df$B == 1) / mean(df$A == 1 & df$B == 1)
  p_011 <- mean(df$Y == 1 & df$A == 0 & df$B == 1) / mean(df$A == 0 & df$B == 1)
  p_001 <- mean(df$Y == 1 & df$A == 0 & df$B == 0) / mean(df$A == 0 & df$B == 0)
  
  # P(Y = 1 | A = a)
  p_11 <- mean(df$Y == 1 & df$A == 1) / mean(df$A == 1)
  p_10 <- mean(df$Y == 1 & df$A == 0) / mean(df$A == 0)
  
  # VE(0)
  ve0_lower <- 1 - p_11 / max(p_001, 1e-6)
  ve0_upper <- 1 - p_110 / max(p_10, 1e-6)
  
  # VE(1)
  ve1_lower <- 1 - p_111 / max(p_10, 1e-6)
  ve1_upper <- 1 - p_11 / max(p_011, 1e-6)
  
  # VEt
  vet_lower <- 1 - p_111 / max(p_001, 1e-6)
  vet_upper <- 1 - p_11 / max(p_10, 1e-6)
  
  bounds_df <- data.frame(
    VE0_lower = ve0_lower,
    VE0_upper = ve0_upper,
    VE1_lower = ve1_lower,
    VE1_upper = ve1_upper,
    VEt_lower = vet_lower,
    VEt_upper = vet_upper
  )
  
  return(bounds_df)
}
