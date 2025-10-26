compute_LP_formula_bounds_S <- function(df,monotonicity_constraint = NULL) {
  bounds_list <- list()
  
  if(monotonicity_constraint=="NO"){
    
    for (s1 in 0:1) {
      for (s2 in 0:1) {
        df_s1 <- subset(df, S == s1)
        df_s2 <- subset(df, S == s2)
        
        # Denominators for conditional probabilities
        denom_1_s1 <- sum(df_s1$A == 1)
        denom_0_s2 <- sum(df_s2$A == 0)
        
        # Probabilities from S = s1 (numerators)
        p_001_s1 <- sum(df_s1$Y == 0 & df_s1$B == 0 & df_s1$A == 1) / denom_1_s1
        p_101_s1 <- sum(df_s1$Y == 1 & df_s1$B == 0 & df_s1$A == 1) / denom_1_s1
        p_011_s1 <- sum(df_s1$Y == 0 & df_s1$B == 1 & df_s1$A == 1) / denom_1_s1
        p_111_s1 <- sum(df_s1$Y == 1 & df_s1$B == 1 & df_s1$A == 1) / denom_1_s1
        
        # Probabilities from S = s2 (denominators)
        p_100_s2 <- sum(df_s2$Y == 1 & df_s2$B == 0 & df_s2$A == 0) / denom_0_s2
        p_000_s2 <- sum(df_s2$Y == 0 & df_s2$B == 0 & df_s2$A == 0) / denom_0_s2
        p_110_s2 <- sum(df_s2$Y == 1 & df_s2$B == 1 & df_s2$A == 0) / denom_0_s2
        p_010_s2 <- sum(df_s2$Y == 0 & df_s2$B == 1 & df_s2$A == 0) / denom_0_s2
        
        # VE0
        ve0_lower <- 1 - (1 - p_001_s1) / max(p_100_s2, 1e-6)
        ve0_upper <- 1 - p_101_s1 / max(1 - p_000_s2, 1e-6)
        
        ve0_upper <- 1 - (probs$p_01 * exps$EY_10+probs$p_11 * 0)/ (probs$p_00 * exps$EY_00+probs$p_10*1)
        ve0_lower <- 1 - (probs$p_01 * exps$EY_10 + probs$p_11 * 1)/( probs$p_00 * exps$EY_00 + probs$p_10 * 0)
        
        # VE1
        ve1_lower <- 1 - (1 - p_011_s1) / max(p_110_s2, 1e-6)
        ve1_upper <- 1 - p_111_s1 / max(1 - p_010_s2, 1e-6)
        
        # VE1
        ve1_upper<- 1 -(probs$p_01 * 0 + probs$p_11 * exps$EY_11)/(probs$p_00 *1 + probs$p_10 * exps$EY_01)
        ve1_lower<- 1 - probs$p_01 * 1 + probs$p_11 * exps$EY_11/(probs$p_00 *0 + probs$p_10 * exps$EY_01)
        
        # VEt
        vet_lower <- 1 - (1 - p_011_s1) / max(p_100_s2, 1e-6)
        vet_upper <- 1 - p_111_s1 / max(1 - p_000_s2, 1e-6)
        
        
        # VEt
        vet_upper <- 1 - (probs$p_01 *0 + probs$p_11 * exps$EY_11) / (probs$p_00 * exps$EY_00 + probs$p_10 * 1)
        vet_lower <- 1 - (probs$p_01 * 1 + probs$p_11 * exps$EY_11) / (probs$p_00 * exps$EY_00 + probs$p_10 * 0)
        
        
        bounds_list[[paste0("s", s1, s2)]] <- c(
          VE0_lower = ve0_lower,
          VE0_upper = ve0_upper,
          VE1_lower = ve1_lower,
          VE1_upper = ve1_upper,
          VEt_lower = vet_lower,
          VEt_upper = vet_upper
        )
      }
    }
    
    
    
    
    }else{
  
  
  for (s1 in 0:1) {
    for (s2 in 0:1) {
      df_s1 <- subset(df, S == s1)
      df_s2 <- subset(df, S == s2)
      
      # Denominators for conditional probabilities
      denom_1_s1 <- sum(df_s1$A == 1)
      denom_0_s2 <- sum(df_s2$A == 0)
      
      # Probabilities from S = s1 (numerators)
      p_001_s1 <- sum(df_s1$Y == 0 & df_s1$B == 0 & df_s1$A == 1) / denom_1_s1
      p_101_s1 <- sum(df_s1$Y == 1 & df_s1$B == 0 & df_s1$A == 1) / denom_1_s1
      p_011_s1 <- sum(df_s1$Y == 0 & df_s1$B == 1 & df_s1$A == 1) / denom_1_s1
      p_111_s1 <- sum(df_s1$Y == 1 & df_s1$B == 1 & df_s1$A == 1) / denom_1_s1
      
      # Probabilities from S = s2 (denominators)
      p_100_s2 <- sum(df_s2$Y == 1 & df_s2$B == 0 & df_s2$A == 0) / denom_0_s2
      p_000_s2 <- sum(df_s2$Y == 0 & df_s2$B == 0 & df_s2$A == 0) / denom_0_s2
      p_110_s2 <- sum(df_s2$Y == 1 & df_s2$B == 1 & df_s2$A == 0) / denom_0_s2
      p_010_s2 <- sum(df_s2$Y == 0 & df_s2$B == 1 & df_s2$A == 0) / denom_0_s2
      
      # VE0
      ve0_lower <- 1 - (1 - p_001_s1) / max(p_100_s2, 1e-6)
      ve0_upper <- 1 - p_101_s1 / max(1 - p_000_s2, 1e-6)
      
      # VE1
      ve1_lower <- 1 - (1 - p_011_s1) / max(p_110_s2, 1e-6)
      ve1_upper <- 1 - p_111_s1 / max(1 - p_010_s2, 1e-6)
      
      # VEt
      vet_lower <- 1 - (1 - p_011_s1) / max(p_100_s2, 1e-6)
      vet_upper <- 1 - p_111_s1 / max(1 - p_000_s2, 1e-6)
      
      bounds_list[[paste0("s", s1, s2)]] <- c(
        VE0_lower = ve0_lower,
        VE0_upper = ve0_upper,
        VE1_lower = ve1_lower,
        VE1_upper = ve1_upper,
        VEt_lower = vet_lower,
        VEt_upper = vet_upper
      )
    }
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
