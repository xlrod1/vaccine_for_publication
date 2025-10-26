library(dplyr)
library(ggplot2)

set.seed(12345)

folder_path <- "yours function path"
# List all R files in the folder
r_files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
# Source each file
sapply(r_files, source)

simulate_data_PO <- function(n,beta_S) {
  # Randomized A and U
  A <- rbinom(n, 1, 0.5)
  U <- rbinom(n, 1, 0.5)
  
  # Parameters (to be tuned for constraints)
  # S depends on A and U
  p_S_0 <- plogis(-1 + 1.5*0 )
  p_S_1 <- plogis(-1 + 1.5*1)
  S_0 <- rbinom(n, 1, p_S_0)
  S_1 <- rbinom(n, 1, p_S_1)
  S<-ifelse(A==1,S_1,S_0)
  # B depends on S and U (increasing in U)
  p_B_0 <- plogis(-1 + beta_S*S_0 + 1.5*U)
  p_B_1 <- plogis(-1 + beta_S*S_1 + 1.5*U)
  
  B_0<-rbinom(n, 1, p_B_0)
  B_1<-rbinom(n, 1, p_B_1)
  
  B <- ifelse(A==1,B_1,B_0)
  
  
  #I wnat to create all teh ppotentoal outcome. under the assumption that 
  #Y^{a,m=-1}|B=m=Y^{a,m}|B=m
  p_Y_00<-plogis(-2 - 1.7*0 +  0.5*0 + 1*U)
  p_Y_01<-plogis(-2 - 1.7*0 +  0.5*1 + 1*U)
  p_Y_10<-plogis(-2 - 1.7*1 +  0.5*0 + 1*U)
  p_Y_11<-plogis(-2 - 1.7*1 +  0.5*1 + 1*U)
  
  Y_00 <- rbinom(n, 1, p_Y_00)
  Y_01 <- rbinom(n, 1, p_Y_01)
  Y_10 <- rbinom(n, 1, p_Y_10)
  Y_11 <- rbinom(n, 1, p_Y_11)
  
  
  df<-data.frame(A, U, S, B, Y_00,Y_01,Y_10,Y_11)
  # And now to observed data
  df <- df %>%
    mutate(Y = case_when(
      A == 0 & B == 0 ~ Y_00,
      A == 0 & B == 1 ~ Y_01,
      A == 1 & B == 0 ~ Y_10,
      A == 1 & B == 1 ~ Y_11
    ))
  
  
  # Y depends on A, S, B, U (monotonically increasing in U and S)
  # Ensure monotonicity in B as well (M-monotonicity)
  # p_Y <- plogis(-2 - 1.7*A + 1.5*S + 0.5*B + 1*U)
  # Y <- rbinom(n, 1, p_Y)
  
  return(df)
}


# ——————————————————————————————————————————————————————————————————————
# 1. Compute “true” bounds on a huge cohort
# ——————————————————————————————————————————————————————————————————————
beta_S =log(2)
set.seed(2025)
N_true  <- 1e6
df_true <- simulate_data_PO(N_true, beta_S =beta_S)

true_LP  <-compute_LP_formula_bounds_with_S(
  df_true,
  monotonicity_constraint = "positive_M"
)
true_mon <- compute_monotonicity_bounds(df_true)
true_LP
true_mon
# Pack into a single named vector:
truth_LP_vec  <- setNames(as.numeric(true_LP[1,]),
                          paste0("LP_",  names(true_LP)))
truth_mon_vec <- setNames(as.numeric(true_mon[1,]),
                          paste0("MON_", names(true_mon)))
truth_vec     <- c(truth_LP_vec, truth_mon_vec)


# ——————————————————————————————————————————————————————————————————————
# 2. Monte Carlo + bootstrap coverage check
# ——————————————————————————————————————————————————————————————————————

set.seed(123)
n <- 1000    # sample size per replicate
B <- 500    # bootstrap resamples
R <- 200     # Monte Carlo replicates

cover_counts <- setNames(rep(0, length(truth_vec)), names(truth_vec))

for (i in seq_len(R)) {
  # 2.1 draw a fresh sample
  df_i <- simulate_data_PO(n, beta_S = beta_S)
  
  # 2.2 bootstrap distribution of bounds
  boot_mat <- replicate(B, {
    idx   <- sample.int(n, n, replace = TRUE)
    dfb   <- df_i[idx, ]
    lp_b  <- compute_LP_formula_bounds_with_S(dfb,
                                            monotonicity_constraint = "positive_M")
    mon_b <- compute_monotonicity_bounds(dfb)
    
    c(setNames(as.numeric(lp_b[1,]),  paste0("LP_",  names(lp_b))),
      setNames(as.numeric(mon_b[1,]), paste0("MON_", names(mon_b))))
  }, simplify = "matrix")
  
  # 2.3 percentile CIs
  cis_i <- apply(boot_mat, 1, quantile, probs = c(0.025, 0.975))
  
  # 2.4 check coverage
  covered      <- (cis_i["2.5%", ] <= truth_vec) &
    (truth_vec <= cis_i["97.5%", ])
  cover_counts <- cover_counts + covered
}

# ——————————————————————————————————————————————————————————————————————
# 3. Summarize empirical coverage
# ——————————————————————————————————————————————————————————————————————

coverage_prob <- cover_counts / R
print(coverage_prob)
# Expect values close to 0.95 if your bootstrap-percentile CIs are well calibrated.