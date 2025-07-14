compute_LP_formula_bounds_with_S<- function(df, monotonicity_constraint = NULL) {
  
  get_pe <- function(df_s) {
    p_01 <- sum(df_s$A == 1 & df_s$B == 0) / sum(df_s$A == 1)
    p_11 <- sum(df_s$A == 1 & df_s$B == 1) / sum(df_s$A == 1)
    p_00 <- sum(df_s$A == 0 & df_s$B == 0) / sum(df_s$A == 0)
    p_10 <- sum(df_s$A == 0 & df_s$B == 1) / sum(df_s$A == 0)
    EY_00 <- mean(df_s$Y[df_s$A == 0 & df_s$B == 0], na.rm=TRUE)
    EY_10 <- mean(df_s$Y[df_s$A == 1 & df_s$B == 0], na.rm=TRUE)
    EY_11 <- mean(df_s$Y[df_s$A == 1 & df_s$B == 1], na.rm=TRUE)
    EY_01 <- mean(df_s$Y[df_s$A == 0 & df_s$B == 1], na.rm=TRUE)
    list(p_01=p_01, p_11=p_11, p_00=p_00, p_10=p_10,
         EY_00=EY_00, EY_10=EY_10, EY_11=EY_11, EY_01=EY_01)
  }
  
  pe_list <- lapply(0:1, function(s) get_pe(subset(df, S == s)))
  
  nums_denoms <- lapply(pe_list, function(pe) {
    with(pe, {
      if (is.null(monotonicity_constraint)) {
        
      
        # --- NO monotonicity ---
        n0_u <- p_01*EY_10 + p_11*0
        d0_u <- p_00*EY_00 + p_10*1
        n0_l <- p_01*EY_10 + p_11*1
        d0_l <- p_00*EY_00 + p_10*0
        
        n1_u <- p_01*0      + p_11*EY_11
        d1_u <- p_00*1      + p_10*EY_01
        n1_l <- p_01*1      + p_11*EY_11
        d1_l <- p_00*0      + p_10*EY_01
        
        nT_u <- n1_u
        dT_u <- d0_u
        nT_l <- n1_l
        dT_l <- d0_l
        
      } else if (monotonicity_constraint == "positive_M") {
      
        
        # --- POS monotonicity ---
        n0_u <- p_01*EY_10 + p_11*0
        d0_u <- p_00*EY_00 + p_10*EY_01
        n0_l <- p_01*EY_10 + p_11*EY_11
        d0_l <- p_00*EY_00 + p_10*0
        
        n1_u <- p_01*EY_10 + p_11*EY_11
        d1_u <- p_00*1      + p_10*EY_01
        n1_l <- p_01*1      + p_11*EY_11
        d1_l <- p_00*EY_00 + p_10*EY_01
        
        nT_u <- n1_u
        dT_u <- d0_u
        nT_l <- n1_l
        dT_l <- d0_l
        
      } else {
        stop("Unknown monotonicity_constraint")
      }
      list(
        num0_u=n0_u, den0_u=d0_u, num0_l=n0_l, den0_l=d0_l,
        num1_u=n1_u, den1_u=d1_u, num1_l=n1_l, den1_l=d1_l,
        numT_u=nT_u, denT_u=dT_u, numT_l=nT_l, denT_l=dT_l
      )
    })
  })
  
  # ——— HERE is the fix ———
  # turn that list‐of‐lists into a true numeric data.frame:
  vecs <- as.data.frame(t(sapply(nums_denoms, unlist)))
  
  # now all these are numeric vectors and outer() will work
  ratio0_u <- outer(vecs$num0_u, vecs$den0_u, `/`)
  ratio0_l <- outer(vecs$num0_l, vecs$den0_l, `/`)
  ratio1_u <- outer(vecs$num1_u, vecs$den1_u, `/`)
  ratio1_l <- outer(vecs$num1_l, vecs$den1_l, `/`)
  ratioT_u <- outer(vecs$numT_u, vecs$denT_u, `/`)
  ratioT_l <- outer(vecs$numT_l, vecs$denT_l, `/`)
  
  out <- data.frame(
    VE0_lower = max(1 - ratio0_l), VE0_upper = min(1 - ratio0_u),
    VE1_lower = max(1-ratio1_l), VE1_upper = min(1 - ratio1_u),
    VEt_lower = max(1 - ratioT_l), VEt_upper =  min(1 -ratioT_u)
  )
  
  return(out)
}