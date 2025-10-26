# Function to compute bootstrap estimates
compute_bootstrap <- function(data, compute_func, monotonicity_constraint) {
  n <- nrow(data)
  boot_results <- replicate(n_boot, {
    boot_sample <- data[sample(1:n, replace = TRUE), ]
    bounds <- compute_func(
      boot_sample,
      monotonicity_constraint = monotonicity_constraint
    )
    unlist(bounds)
  })
  boot_results <- t(boot_results)
  colnames(boot_results) <- colnames(compute_func(data, monotonicity_constraint))
  boot_results <- as.data.frame(boot_results)
  return(boot_results)
}

# Helper function to get percentile intervals
get_percentile_CI <- function(boot_df, level = 0.95) {
  lower_bound <- (1 - level) / 2
  upper_bound <- 1 - lower_bound
  
  ci <- apply(boot_df, 2, quantile, probs = c(lower_bound, upper_bound), na.rm = TRUE)
  ci <- as.data.frame(t(ci))
  names(ci) <- c("CI_lower", "CI_upper")
  
  return(ci)
}

# Function to combine point estimate and bootstrap CI
summarize_bounds_with_CI <- function(original_bounds, boot_bounds, level = 0.95) {
  lower <- (1 - level) / 2
  upper <- 1 - lower
  
  CI_mat <- apply(boot_bounds, 2, quantile, probs = c(lower, upper), na.rm = TRUE)
  
  summary_df <- data.frame(
    Bound = names(original_bounds),
    Estimate = unlist(original_bounds),
    CI_lower = CI_mat[1, ],
    CI_upper = CI_mat[2, ]
  )
  
  return(summary_df)
}
