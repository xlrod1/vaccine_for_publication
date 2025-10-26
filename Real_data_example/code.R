library(tidyverse)
library(haven)
library(dplyr)
library(tidyr)
library(stringr)
library(gt)

# -------------------------------------------------------------------
# Mock data generator (added so the script runs without real data)
# -------------------------------------------------------------------
set.seed(137)
n <- 1000
df <- tibble(
  A = rbinom(n, 1, 0.5),
  S = rbinom(n, 1, plogis(-0.3 + 1.0 * A)),
  Y = rbinom(n, 1, plogis(-2 + 0.8 * A + 0.4 * S))
)
# -------------------------------------------------------------------

# Function to compute VE(0), VE(1), VE_T
compute_VEs <- function(data) {
  VE0 <- 1 - mean(data$Y[data$A == 1 & data$B == 0]) /
    mean(data$Y[data$A == 0 & data$B == 0])
  VE1 <- 1 - mean(data$Y[data$A == 1 & data$B == 1]) /
    mean(data$Y[data$A == 0 & data$B == 1])
  VEt <- 1 - mean(data$Y[data$A == 1 & data$B == 1]) /
    mean(data$Y[data$A == 0 & data$B == 1])
  return(c(VE0 = VE0, VE1 = VE1, VEt = VEt))
}

# upload the functions
folder_path <- "functions path"
r_files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
sapply(r_files, source)


#-------------------------------------------------------------------
#Simulate B
#-------------------------------------------------------------------
set.seed(137)
df$B[df$S == 1 & df$Y == 1] <- rbinom(sum(df$S == 1 & df$Y == 1), 1, 0.7)
df$B[df$S == 0 & df$Y == 1] <- rbinom(sum(df$S == 0 & df$Y == 1), 1, 0.2)
df$B[df$S == 1 & df$Y == 0] <- rbinom(sum(df$S == 1 & df$Y == 0), 1, 0.5)
df$B[df$S == 0 & df$Y == 0] <- rbinom(sum(df$S == 0 & df$Y == 0), 1, 0.1)

#-------------------------------------------------------------------
#Estimate the point-estimates of the VEs, using Assumptoion 1-3
#-------------------------------------------------------------------
estimates <- compute_VEs(df)
estimates

# -------------------------------------------------------------------
# Add Bootstrap CI for the point-estimates
# -------------------------------------------------------------------
B <- 1000  # number of bootstrap samples
boot_results <- replicate(B, {
  sample_idx <- sample(seq_len(nrow(df)), replace = TRUE)
  compute_VEs(df[sample_idx, ])
})

boot_results <- t(boot_results)
colnames(boot_results) <- c("VE0", "VE1", "VEt")

# 95% percentile CI
boot_ci <- apply(boot_results, 2, quantile, probs = c(0.025, 0.975))

list(
  estimates = round(estimates * 100, 1),
  CI = round(boot_ci * 100, 1)
)

# -------------------------------------------------------------------
# Estiamte non-parametric bounds
# -------------------------------------------------------------------
n_boot <- 1000
#set.seed(123)
#----without monotonicity assumption, under figure 3d
LP_no_S  <- compute_LP_formula_bounds_no_S(df, monotonicity_constraint = "NO")
mon_no_S <- compute_monotonicity_bounds_no_S(df, monotonicity_constraint = "NO")
#----without monotonicity assumption, under figure 3a
LP_S  <- compute_LP_formula_bounds_with_S(df, monotonicity_constraint = "NO")
mon_S <- compute_monotonicity_bounds(df, monotonicity_constraint = "NO")

# ----with monotonicity assumption, under figure 3d
LP_no_S_M  <- compute_LP_formula_bounds_no_S(df, monotonicity_constraint = "positive_M")
mon_no_S_M <- compute_monotonicity_bounds_no_S(df, monotonicity_constraint = "positive_M")
#----with monotonicity assumption, under figure 3a
LP_S_M  <- compute_LP_formula_bounds_with_S(df, monotonicity_constraint = "positive_M")
mon_S_M <- compute_monotonicity_bounds(df, monotonicity_constraint = "positive_M")

# -------------------------------------------------------------------
# Ading Bootsrap CIs to the bounds
# -------------------------------------------------------------------
#----without monotonicity assumption, under figure 3d
boot_LP_no_S  <- compute_bootstrap(df, compute_LP_formula_bounds_no_S, "NO")
boot_mon_no_S <- compute_bootstrap(df, compute_monotonicity_bounds_no_S, "NO")
#----without monotonicity assumption, under figure 3a
boot_LP_S  <- compute_bootstrap(df, compute_LP_formula_bounds_with_S, "NO")
boot_mon_S <- compute_bootstrap(df, compute_monotonicity_bounds, "NO")

#----with monotonicity assumption, under figure 3d
boot_LP_no_S_M  <- compute_bootstrap(df, compute_LP_formula_bounds_no_S, "positive_M")
boot_mon_no_S_M <- compute_bootstrap(df, compute_monotonicity_bounds_no_S, "positive_M")
# ----with monotonicity assumption, under figure 3d
boot_LP_S_M  <- compute_bootstrap(df, compute_LP_formula_bounds_with_S, "positive_M")
boot_mon_S_M <- compute_bootstrap(df, compute_monotonicity_bounds, "positive_M")

# -------------------------------------------------------------------
# Example: Get 95% percentile CI
# -------------------------------------------------------------------
summary_LP_no_S  <- summarize_bounds_with_CI(LP_no_S, boot_LP_no_S)
summary_mon_no_S <- summarize_bounds_with_CI(mon_no_S, boot_mon_no_S)
summary_LP_no_S_M  <- summarize_bounds_with_CI(LP_no_S_M, boot_LP_no_S_M)
summary_mon_no_S_M <- summarize_bounds_with_CI(mon_no_S_M, boot_mon_no_S_M)

summary_LP_S  <- summarize_bounds_with_CI(LP_S, boot_LP_S)
summary_mon_S <- summarize_bounds_with_CI(mon_S, boot_mon_S)
summary_LP_S_M  <- summarize_bounds_with_CI(LP_S_M, boot_LP_S_M)
summary_mon_S_M <- summarize_bounds_with_CI(mon_S_M, boot_mon_S_M)

# -------------------------------------------------------------------
# Figure 4 – display numeric summaries
# -------------------------------------------------------------------
summary_LP_no_S  %>% mutate(across(where(is.numeric), ~round(. * 100, 1)))
summary_mon_no_S %>% mutate(across(where(is.numeric), ~round(. * 100, 1)))
summary_LP_no_S_M  %>% mutate(across(where(is.numeric), ~round(. * 100, 1)))
summary_mon_no_S_M %>% mutate(across(where(is.numeric), ~round(. * 100, 1)))

# -------------------------------------------------------------------
# Figure 3 – display numeric summaries
# -------------------------------------------------------------------
summary_LP_S  %>% mutate(across(where(is.numeric), ~round(. * 100, 1)))
summary_mon_S %>% mutate(across(where(is.numeric), ~round(. * 100, 1)))
summary_LP_S_M  %>% mutate(across(where(is.numeric), ~round(. * 100, 1)))
summary_mon_S_M %>% mutate(across(where(is.numeric), ~round(. * 100, 1)))
