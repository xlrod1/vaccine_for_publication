check_Y_monotonicity_U_nondecreasing <- function(df) {
  results <- list()
  for (a_val in unique(df$A)) {
    for (s_val in unique(df$S)) {
      subdf <- df[df$A == a_val & df$S == s_val & df$B == 1, ]
      if (nrow(subdf) == 0) next
      pr_Y_given_U <- tapply(subdf$Y, subdf$U, mean)
      nondecreasing <- all(diff(pr_Y_given_U) >= 0)
      results[[paste0("A=", a_val, "_S=", s_val)]] <- list(
        pr_Y_given_U = pr_Y_given_U,
        nondecreasing = nondecreasing
      )
    }
  }
  return(results)
}


check_Y_monotonicity_B_nondecreasing <- function(df) {
  results <- list()
  for (a_val in unique(df$A)) {
    for (s_val in unique(df$S)) {
      for (u_val in unique(df$U)) {
        subdf <- df[df$A == a_val & df$S == s_val & df$U == u_val, ]
        if (nrow(subdf) == 0) next
        pr_Y_given_B <- tapply(subdf$Y, subdf$B, mean)
        nondecreasing <- all(diff(pr_Y_given_B) >= 0)
        results[[paste0("A=", a_val, "_S=", s_val, "_U=", u_val)]] <- list(
          pr_Y_given_B = pr_Y_given_B,
          nondecreasing = nondecreasing
        )
      }
    }
  }
  return(results)
}


check_M_monotonicity <- function(df) {
  results <- list()
  for (a_prime in unique(df$A)) {
    for (b in unique(df$B)) {
      for (s in unique(df$S)) {
        idx <- which(df$A == a_prime & df$B == b & df$S == s)
        if (length(idx) == 0) next
        # Empirical P(U=u)
        p_u <- table(df$U) / nrow(df)
        # For m = 0 and m = 1
        EY_am <- c()
        for (m in 0:1) {
          EY_um <- c()
          for (u in 0:1) {
            sub_idx <- which(df$A == a_prime & df$B == m & df$S == s & df$U == u)
            if (length(sub_idx) == 0) {
              EY_um <- c(EY_um, 0)
            } else {
              EY_um <- c(EY_um, mean(df$Y[sub_idx]))
            }
          }
          # weighted sum over U
          EY_am <- c(EY_am, sum(EY_um * as.numeric(p_u)))
        }
        # Observed mean
        EY_obs <- mean(df$Y[idx])
        # Check monotonicity
        monotone <- (EY_am[1] <= EY_obs) && (EY_obs <= EY_am[2])
        results[[paste0("A=", a_prime, "_B=", b, "_S=", s)]] <- list(
          EY_am0 = EY_am[1],
          EY_obs = EY_obs,
          EY_am1 = EY_am[2],
          monotone = monotone
        )
      }
    }
  }
  return(results)
}
