bootstrap_LP_bounds_CI <- function(data, n_boot = 500, conf_level = 0.95, seed = NULL, monotonicity_constraint = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Column names expected from get_all_LP_bounds_summary()
  bound_names <- c("VE0_lower", "VE0_upper", "VE1_lower", "VE1_upper", "VEt_lower", "VEt_upper")
  boot_matrix <- matrix(NA, nrow = n_boot, ncol = length(bound_names))
  colnames(boot_matrix) <- bound_names
  
  for (i in 1:n_boot) {
    # Resample with replacement
    boot_sample <- data[sample(nrow(data), replace = TRUE), ]
    
    # Compute bounds, skip iteration on failure
    try({
      bounds_df <- get_all_LP_bounds_summary(boot_sample, monotonicity_constraint = monotonicity_constraint)
      if (all(bound_names %in% names(bounds_df))) {
        boot_matrix[i, ] <- as.numeric(bounds_df[1, bound_names])
      }
    }, silent = TRUE)
  }
  
  # Remove failed iterations
  boot_matrix <- na.omit(boot_matrix)
  
  # Compute percentile CIs
  alpha <- (1 - conf_level) / 2
  CI_bounds <- apply(boot_matrix, 2, quantile, probs = c(alpha, 1 - alpha), na.rm = TRUE)
  
  # Return as a tidy data frame
  data.frame(
    estimand = bound_names,
    CI_lower = CI_bounds[1, ],
    CI_upper = CI_bounds[2, ]
  )
}
