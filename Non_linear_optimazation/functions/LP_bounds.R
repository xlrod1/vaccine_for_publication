library(nloptr)

compute_probabilities_S <- function(df, s_value) {
  df_s <- subset(df, S == s_value)
  with(df_s, list(
    p_01 = sum(A == 1 & B == 0) / sum(A == 1),  # P(B=0 | A=1, S=s)
    p_11 = sum(A == 1 & B == 1) / sum(A == 1),  # P(B=1 | A=1, S=s)
    p_00 = sum(A == 0 & B == 0) / sum(A == 0),  # P(B=0 | A=0, S=s)
    p_10 = sum(A == 0 & B == 1) / sum(A == 0)   # P(B=1 | A=0, S=s)
  ))
}

compute_expectations_S <- function(df, s_value) {
  df_s <- subset(df, S == s_value)
  with(df_s, list(
    EY_00 = mean(Y[A == 0 & B == 0], na.rm = TRUE),
    EY_10 = mean(Y[A == 1 & B == 0], na.rm = TRUE),
    EY_11 = mean(Y[A == 1 & B == 1], na.rm = TRUE),
    EY_01 = mean(Y[A == 0 & B == 1], na.rm = TRUE)
  ))
}

bounds_optimized_VE_with_S <- function(df, VE_type = "VE0", monotonicity_constraint = NULL) {
  results_all <- list()
  
  for (s1 in 0:1) {
    for (s2 in 0:1) {
      probs_s1 <- compute_probabilities_S(df, s1)
      exps_s1 <- compute_expectations_S(df, s1)
      probs_s2 <- compute_probabilities_S(df, s2)
      exps_s2 <- compute_expectations_S(df, s2)
      
      p_01_1 <- probs_s1$p_01
      p_11_1 <- probs_s1$p_11
      p_00_2 <- probs_s2$p_00
      p_10_2 <- probs_s2$p_10
      
      EY_00_2 <- exps_s2$EY_00
      EY_10_1 <- exps_s1$EY_10
      EY_11_1 <- exps_s1$EY_11
      EY_01_2 <- exps_s2$EY_01
      
      objective <- function(params) {
        x <- params[1]
        y <- params[2]
        z <- params[3]
        q <- params[4]
        
        if (tolower(VE_type) == "ve0") {
          N <- p_01_1 * EY_10_1 + p_11_1 * x
          D <- p_00_2 * EY_00_2 + p_10_2 * y
        } else if (tolower(VE_type) == "ve1") {
          N <- p_01_1 * z + p_11_1 * EY_11_1
          D <- p_00_2 * q + p_10_2 * EY_01_2
        } else if (tolower(VE_type) == "vet") {
          N <- p_01_1 * z + p_11_1 * EY_11_1
          D <- p_00_2 * EY_00_2 + p_10_2 * y
        } else stop("Invalid VE_type")
        
        if (D == 0) return(Inf)
        return(1 - (N / D))
      }
      
      constraint <- function(params) {
        x <- params[1]; y <- params[2]; z <- params[3]; q <- params[4]
        constraints <- numeric(0)
        
        if (!is.null(monotonicity_constraint)) {
          if (monotonicity_constraint == "positive_M") {
            constraints <- c(
              EY_10_1 - z,  # Z >= EY_10_1 => EY_10_1 - z <= 0
              EY_00_2 - q,  # Q >= EY_00_2
              x - EY_11_1,  # X <= EY_11_1 => x - EY_11_1 <= 0
              y - EY_01_2   # Y <= EY_01_2
            )
          } else if (monotonicity_constraint == "negative_M") {
            constraints <- c(
              z - EY_10_1,  # Z <= EY_10_1 => z - EY_10_1 <= 0
              q - EY_00_2,  # Q <= EY_00_2
              EY_11_1 - x,  # X >= EY_11_1
              EY_01_2 - y   # Y >= EY_01_2
            )
          }
        }
        
        return(as.numeric(constraints))
      }
      
      opt_upper <- nloptr(x0 = rep(0.5, 4),
                          eval_f = function(p) -objective(p),
                          lb = rep(0, 4), ub = rep(1, 4),
                          eval_g_ineq = constraint,
                          opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-7, maxeval = 1000))
      
      opt_lower <- nloptr(x0 = rep(0.5, 4),
                          eval_f = objective,
                          lb = rep(0, 4), ub = rep(1, 4),
                          eval_g_ineq = constraint,
                          opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-7, maxeval = 1000))
      
      results_all[[paste0("s", s1, s2)]] <- list(
        lower = opt_lower$objective,
        upper = -opt_upper$objective
      )
    }
  }
  
  lower_bound <- max(sapply(results_all, function(res) res$lower), na.rm = TRUE)
  upper_bound <- min(sapply(results_all, function(res) res$upper), na.rm = TRUE)
  
  return(list(
    VE_type = VE_type,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    all_results = results_all
  ))
}
get_all_LP_bounds_summary <- function(df, monotonicity_constraint = NULL) {
  res_ve0 <- bounds_optimized_VE_with_S(df, VE_type = "VE0", monotonicity_constraint = monotonicity_constraint)
  res_ve1 <- bounds_optimized_VE_with_S(df, VE_type = "VE1", monotonicity_constraint = monotonicity_constraint)
  res_vet <- bounds_optimized_VE_with_S(df, VE_type = "VEt", monotonicity_constraint = monotonicity_constraint)
  
  return(data.frame(
    VE0_lower = res_ve0$lower_bound,
    VE0_upper = res_ve0$upper_bound,
    VE1_lower = res_ve1$lower_bound,
    VE1_upper = res_ve1$upper_bound,
    VEt_lower = res_vet$lower_bound,
    VEt_upper = res_vet$upper_bound
  ))
}
