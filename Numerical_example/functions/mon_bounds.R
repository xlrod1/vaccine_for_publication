compute_monotonicity_bounds <- function(df,monotonicity_constraint ="NO") {
  bounds_list <- list()
  
  if(monotonicity_constraint=="positive_M"){
  for (s in 0:1) {
    df_s <- subset(df, S == s)
    
    # P(Y = 1 | A = a, B = b, S = s)
    p_110s <- mean(df_s$Y == 1 & df_s$A == 1 & df_s$B == 0) / mean(df_s$A == 1 & df_s$B == 0)
    p_111s <- mean(df_s$Y == 1 & df_s$A == 1 & df_s$B == 1) / mean(df_s$A == 1 & df_s$B == 1)
    p_011s <- mean(df_s$Y == 1 & df_s$A == 0 & df_s$B == 1) / mean(df_s$A == 0 & df_s$B == 1)
    p_001s <- mean(df_s$Y == 1 & df_s$A == 0 & df_s$B == 0) / mean(df_s$A == 0 & df_s$B == 0)
    
    # P(Y = 1 | A = a)
    p_11 <- mean(df$Y == 1 & df$A == 1) / mean(df$A == 1)
    p_10 <- mean(df$Y == 1 & df$A == 0) / mean(df$A == 0)
    
    # VE(0)
    ve0_lower <- 1 - p_11 / max(p_001s, 1e-6)
    ve0_upper <- 1 - p_110s / max(p_10, 1e-6)
    
    # VE(1)
    ve1_lower <- 1 - p_111s / max(p_10, 1e-6)
    ve1_upper <- 1 - p_11 / max(p_011s, 1e-6)
    
    # VEt
    vet_lower <- 1 - p_111s / max(p_001s, 1e-6)
    vet_upper <- 1 - p_11 / p_10
    
    bounds_list[[paste0("s", s)]] <- c(
      VE0_lower = ve0_lower,
      VE0_upper = ve0_upper,
      VE1_lower = ve1_lower,
      VE1_upper = ve1_upper,
      VEt_lower = vet_lower,
      VEt_upper = vet_upper
    )
  }
  }else{
    
    for (s in 0:1) {
      df_s <- subset(df, S == s)
      
      # P(Y = 1 | A = a, B = b, S = s)
      p_110s <- mean(df_s$Y == 1 & df_s$A == 1 & df_s$B == 0) / mean(df_s$A == 1 & df_s$B == 0)
      p_111s <- mean(df_s$Y == 1 & df_s$A == 1 & df_s$B == 1) / mean(df_s$A == 1 & df_s$B == 1)
      p_011s <- mean(df_s$Y == 1 & df_s$A == 0 & df_s$B == 1) / mean(df_s$A == 0 & df_s$B == 1)
      p_001s <- mean(df_s$Y == 1 & df_s$A == 0 & df_s$B == 0) / mean(df_s$A == 0 & df_s$B == 0)
      
      
      
      # VE(0)
    
      ve0_lower <- 1 - 1 / max(p_001s, 1e-6)
      ve0_upper <- 1 - p_110s 
      # VE(1)
      ve1_lower <- -10^6
      ve1_upper <- 1 
      
      # VEt
      vet_lower <- 1 - p_111s / max(p_001s, 1e-6)
      vet_upper <- 1 
      
      bounds_list[[paste0("s", s)]] <- c(
        VE0_lower = ve0_lower,
        VE0_upper = ve0_upper,
        VE1_lower = ve1_lower,
        VE1_upper = ve1_upper,
        VEt_lower = vet_lower,
        VEt_upper = vet_upper
      )
    }
  }
  
  bounds_df <- do.call(rbind, bounds_list)
  
  final_bounds <- data.frame(
    VE0_lower = max(bounds_df[, "VE0_lower"], na.rm = TRUE),
    VE0_upper = min(bounds_df[, "VE0_upper"], na.rm = TRUE),
    VE1_lower = max(bounds_df[, "VE1_lower"], na.rm = TRUE),
    VE1_upper = min(bounds_df[, "VE1_upper"], na.rm = TRUE),
    VEt_lower = max(bounds_df[, "VEt_lower"], na.rm = TRUE),
    VEt_upper = min(bounds_df[, "VEt_upper"], na.rm = TRUE)
  )
  
  return(final_bounds)
}
