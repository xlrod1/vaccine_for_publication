library(dplyr)
library(ggplot2)

set.seed(12345)

folder_path <- "yours function path"
# List all R files in the folder
r_files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
# Source each file
sapply(r_files, source)

simulate_data_PO <- function(n,gamma_M_0,gamma_M_1) {
  # Randomized A and U
  A <- rbinom(n, 1, 0.5)
  U <- rbinom(n, 1, 0.5)
  
  # Parameters (to be tuned for constraints)
  # S depends on A and U
  p_S_0 <- plogis(-1 + 1.5*0 + 1*U)
  p_S_1 <- plogis(-1 + 1.5*1 + 1*U)
  S_0 <- rbinom(n, 1, p_S_0)
  S_1 <- rbinom(n, 1, p_S_1)
  S<-ifelse(A==1,S_1,S_0)
  # B depends on S and U (increasing in U)
  p_B_0 <- plogis(-1 + 1.5*S_0 + 1.5*U)
  p_B_1 <- plogis(-1 + 1.5*S_1 + 1.5*U)
  
  B_0<-rbinom(n, 1, p_B_0)
  B_1<-rbinom(n, 1, p_B_1)
  
  B <- ifelse(A==1,B_1,B_0)
  
  
  #I wnat to create all teh ppotentoal outcome. under the assumption that 
  #Y^{a,m=-1}|B=m=Y^{a,m}|B=m
  p_Y_00<-plogis(-2 - 1.7*0 + 1.5*S_0 + gamma_M_0*0 + 1*U)
  p_Y_01<-plogis(-2 - 1.7*0 + 1.5*S_0 + gamma_M_0*1 + 1*U)
  p_Y_10<-plogis(-2 - 1.7*1 + 1.5*S_1 + gamma_M_1*0 + 1*U)
  p_Y_11<-plogis(-2 - 1.7*1 + 1.5*S_1 + gamma_M_1*1 + 1*U)
  
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

#1. investigate hoe beta_S effcts the width of the boounds 

all_bounds <- list()

# Define beta_S values
# Fixed parameters
gamma_M0_fixed <- log(1.5)
gamma_M1_vals <- c(log(1.5), 0, -log(1.5),-log(3),-log(5))  # specific values only
# Loop over beta_S values
for (gamma_M1 in gamma_M1_vals) {
  df <- simulate_data_PO(1000000,gamma_M_0=gamma_M0_fixed ,gamma_M_1=gamma_M1)
  
  # True VE values
  VE_0 <- 1 - mean(df$Y_10) / mean(df$Y_00)
  VE_1 <- 1 - mean(df$Y_11) / mean(df$Y_01)
  VE_T <- 1 - mean(df$Y_11) / mean(df$Y_00)
  
  # Bounds from both methods
  mon_no_s <- compute_monotonicity_bounds_no_S(df,monotonicity_constraint = "positive_M")
  LP_M_no_S <- get_all_LP_bounds_summary_no_s(df, monotonicity_constraint = "positive_M")
  
  # Combine bounds into long-format data
  bounds_df <- bind_rows(
    data.frame(
      gamma_M1 = gamma_M1,
      VE_type = c("VE0", "VE1", "VEt"),
      Method = "Monotonicity",
      Lower = c(mon_no_s$VE0_lower, mon_no_s$VE1_lower, mon_no_s$VEt_lower),
      Upper = c(mon_no_s$VE0_upper, mon_no_s$VE1_upper, mon_no_s$VEt_upper),
      True = c(VE_0, VE_1, VE_T)
    ),
    data.frame(
      gamma_M1 = gamma_M1,
      VE_type = c("VE0", "VE1", "VEt"),
      Method = "LP",
      Lower = c(LP_M_no_S$VE0_lower, LP_M_no_S$VE1_lower, LP_M_no_S$VEt_lower),
      Upper = c(LP_M_no_S$VE0_upper, LP_M_no_S$VE1_upper, LP_M_no_S$VEt_upper),
      True = c(VE_0, VE_1, VE_T)
    )
  )
  
  all_bounds[[as.character(gamma_M1)]] <- bounds_df
}

# Combine all into one dataframe
plot_df <- bind_rows(all_bounds)
plot_df<-plot_df %>% filter(Lower<Upper)

# Treat beta_S as a categorical variable with the correct order
pos    <- position_dodge(width = 0.6)
cap_hw <- 0.12
tol    <- 1e-8


# helper to map numeric gamma_M1 to math-label strings
make_gamma_label <- function(v, tol = 1e-8) {
  
  if (abs(v - (-log(5)))   < tol) return("-log(5)")
  if (abs(v - (-log(3)))   < tol) return("-log(3)")
  if (abs(v - (-log(1.5))) < tol) return("-log(1.5)")
  if (abs(v) < tol) return("0")
  if (abs(v -  log(1.5)) < tol) return("log(1.5)")
  

  # fallback: plain number
  format(v, digits = 3, trim = TRUE, scientific = FALSE)
}


# ----- Prep data -----
plot_df2 <- plot_df %>%
  mutate(
    Lower100 = Lower * 100,
    Upper100 = Upper * 100,
    True100  = True  * 100,
    VE_facet = factor(
      VE_type,
      levels = c("VE0","VE1","VE_T","VEt"),
      labels = c("VE(0)","VE(1)","VE[T]","VE[T]")
    ),
    gamma_lab = sapply(gamma_M1, make_gamma_label)
  )
desired_order <- c("-log(5)", "-log(3)", "-log(1.5)", "0", "log(1.5)")
plot_df2 <- plot_df2 %>% mutate(gamma_lab = factor(gamma_lab, levels = desired_order))

df_ok   <- plot_df2 %>% filter(Lower100 >= 0)
df_clip <- plot_df2 %>% filter(Lower100 < 0)

# arrows at 0 (per method)
arrows_below <- df_clip %>% mutate(y_start = 2, y_end = 0)
#   mutate(y_start = 2, y_end = 0)

df_true <- plot_df2 %>% distinct(VE_facet, gamma_lab, True100)

cap_w <- 0.28                       # cap width

# ----- Plot -----
ggplot(plot_df2, aes(x = gamma_lab, group = VE_type,type=Method)) +
  # regular intervals (both caps)
  geom_errorbar(
    data = df_ok,
    aes(ymin = pmax(Lower100, 0), ymax = pmin(Upper100, 100)),
    width = 0.24, position = pos
  ) +
  # vertical line for Lower < 0: from 0 up to upper
  geom_segment(
    data = df_clip,
    aes(x = gamma_lab, xend = gamma_lab, y = 0, yend = pmin(Upper100, 100)),
    position = pos
  ) +
  # top cap at the upper end (ymin==ymax)
  geom_errorbar(
    data = df_clip,
    aes(ymin = pmin(Upper100, 100), ymax = pmin(Upper100, 100)),
    width = 0.24, position = pos
  ) +
  # down arrow at 0
  geom_segment(
    data = arrows_below,
    aes(x = gamma_lab, xend = gamma_lab, y = y_start, yend = y_end),
    position = pos,
    arrow = arrow(length = unit(3, "pt"), ends = "last", type = "closed")
  ) +
  # true value
  geom_point(aes(y = True100), size = 1.3, position = pos) +
  facet_grid(VE_facet ~ Method, labeller = label_parsed) +
  labs(x = expression(gamma[B^1]), y = "True VE and estimated bounds (%)") +
  scale_x_discrete(labels = function(x) parse(text = x), drop = FALSE) +
  coord_cartesian(ylim = c(0, 100), expand = FALSE) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_text(hjust = 0.5))















