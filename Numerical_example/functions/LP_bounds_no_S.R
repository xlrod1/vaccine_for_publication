library(nloptr)

compute_probabilities_no_s <- function(df) {
  with(df, list(
    p_01 = sum(A == 1 & B == 0) / sum(A == 1),  # P(B=0 | A=1)
    p_11 = sum(A == 1 & B == 1) / sum(A == 1),  # P(B=1 | A=1)
    p_00 = sum(A == 0 & B == 0) / sum(A == 0),  # P(B=0 | A=0)
    p_10 = sum(A == 0 & B == 1) / sum(A == 0)   # P(B=1 | A=0)
  ))
}

compute_expectations_no_s <- function(df) {
  with(df, list(
    EY_00 = mean(Y[A == 0 & B == 0], na.rm = TRUE),
    EY_10 = mean(Y[A == 1 & B == 0], na.rm = TRUE),
    EY_11 = mean(Y[A == 1 & B == 1], na.rm = TRUE),
    EY_01 = mean(Y[A == 0 & B == 1], na.rm = TRUE)
  ))
}

bounds_optimized_VE_no_s <- function(df, VE_type = "VE0", monotonicity_constraint = NULL) {
  probs <- compute_probabilities_no_s(df)
  exps <- compute_expectations_no_s(df)
  
  objective <- function(params) {
    x <- params[1]
    y <- params[2]
    z <- params[3]
    q <- params[4]
    
    if (tolower(VE_type) == "ve0") {
      N <- probs$p_01 * exps$EY_10 + probs$p_11 * x
      D <- probs$p_00 * exps$EY_00 + probs$p_10 * y
    } else if (tolower(VE_type) == "ve1") {
      N <- probs$p_01 * z + probs$p_11 * exps$EY_11
      D <- probs$p_00 * q + probs$p_10 * exps$EY_01
    } else if (tolower(VE_type) == "vet") {
      N <- probs$p_01 * z + probs$p_11 * exps$EY_11
      D <- probs$p_00 * exps$EY_00 + probs$p_10 * y
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
          exps$EY_10 - z,  # Z >= EY_10
          exps$EY_00 - q,  # Q >= EY_00
          x - exps$EY_11,  # X <= EY_11
          y - exps$EY_01   # Y <= EY_01
        )
      } else if (monotonicity_constraint == "negative_M") {
        constraints <- c(
          z - exps$EY_10,  # Z <= EY_10
          q - exps$EY_00,  # Q <= EY_00
          exps$EY_11 - x,  # X >= EY_11
          exps$EY_01 - y   # Y >= EY_01
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
  
  return(list(
    VE_type = VE_type,
    lower_bound = opt_lower$objective,
    upper_bound = -opt_upper$objective
  ))
}

get_all_LP_bounds_summary_no_s <- function(df, monotonicity_constraint = NULL) {
  res_ve0 <- bounds_optimized_VE_no_s(df, VE_type = "VE0", monotonicity_constraint = monotonicity_constraint)
  res_ve1 <- bounds_optimized_VE_no_s(df, VE_type = "VE1", monotonicity_constraint = monotonicity_constraint)
  res_vet <- bounds_optimized_VE_no_s(df, VE_type = "VEt", monotonicity_constraint = monotonicity_constraint)
  
  return(data.frame(
    VE0_lower = res_ve0$lower_bound,
    VE0_upper = res_ve0$upper_bound,
    VE1_lower = res_ve1$lower_bound,
    VE1_upper = res_ve1$upper_bound,
    VEt_lower = res_vet$lower_bound,
    VEt_upper = res_vet$upper_bound
  ))
}
