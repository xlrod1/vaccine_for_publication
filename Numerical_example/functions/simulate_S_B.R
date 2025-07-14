simulate_U_B <- function(data, scenario = c(1, 2, 3, 4),
                         mu_0 = -0.5, mu_Y = 1.0, mu_S = 1.0,
                         beta_0 = -0.5, beta_S = 0.8, beta_Y = 1.0, beta_U = 1.5) {
  scenario <- match.arg(as.character(scenario), choices = c("1", "2", "3", "4"))
  n <- nrow(data)
  
  # Simulate U
  if (scenario %in% c("1", "2")) {
    data$U <- rnorm(n, mean = mu_0 + mu_Y * data$Y, sd = 1)
  } else if (scenario %in% c("3", "4")) {
    data$U <- rnorm(n, mean = mu_0 + mu_S * data$S + mu_Y * data$Y, sd = 1)
  }
  
  # Simulate B
  logit_B <- beta_0 + beta_S * data$S + beta_Y * data$Y + beta_U * data$U
  prob_B <- 1 / (1 + exp(-logit_B))
  data[[paste0("B_scenario", scenario)]] <- rbinom(n, 1, prob_B)
  
  return(data)
}
