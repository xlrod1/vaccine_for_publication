library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

set.seed(12345)

folder_path <- "yours fnctio path"
# List all R files in the folder
r_files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
# Source each file
sapply(r_files, source)

simulate_data_PO_not_mon <- function(n,beta_U,gamma_U) {
  # Randomized A and U
  A <- rbinom(n, 1, 0.5)
  U <- rnorm(n,0,1)
  U_2<-U*U
  # Parameters (to be tuned for constraints)
  # S depends on A and U
  p_S_0 <- plogis(-1 + 1.5*0 + 1*U)
  p_S_1 <- plogis(-1 + 1.5*1 + 1*U)
  S_0 <- rbinom(n, 1, p_S_0)
  S_1 <- rbinom(n, 1, p_S_1)
  S<-ifelse(A==1,S_1,S_0)
  # B depends on S and U (increasing in U)
  p_B_0 <- plogis(-1 + 1.5*S_0 + beta_U*U_2)
  p_B_1 <- plogis(-1 + 1.5*S_1 + beta_U*U_2)
  
  B_0<-rbinom(n, 1, p_B_0)
  B_1<-rbinom(n, 1, p_B_1)
  
  B <- ifelse(A==1,B_1,B_0)
  
  
  #I wnat to create all teh ppotentoal outcome. under the assumption that 
  #Y^{a,m=-1}|B=m=Y^{a,m}|B=m
  p_Y_00<-plogis(-2 - 1.7*0 + 1.5*S_0 + 0.5*0 +gamma_U*U_2)
  p_Y_01<-plogis(-2 - 1.7*0 + 1.5*S_0 + 0.5*1 +gamma_U*U_2)
  p_Y_10<-plogis(-2 - 1.7*1 + 1.5*S_1 + 0.5*0 +gamma_U*U_2)
  p_Y_11<-plogis(-2 - 1.7*1 + 1.5*S_1 + 0.5*1 +gamma_U*U_2)
  
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
}
  
  # Define beta_U values
  beta_vals <- c(-log(5),-log(2), 0, log(2),log(5))
  
  # Define gamma_U values
  gamma_vals <- c(-log(2), log(2))
  
  all_bounds <- list()
  
  # Outer loop over gamma_U
  for (g in gamma_vals) {
    gamma_label <- paste0("gamma[U]==", g)
    
    for (b in beta_vals) {
      df <- simulate_data_PO_not_mon(1000000, beta_U = b, gamma_U = g)
      
      VE_0 <- 1 - mean(df$Y_10) / mean(df$Y_00)
      VE_1 <- 1 - mean(df$Y_11) / mean(df$Y_01)
      VE_T <- 1 - mean(df$Y_11) / mean(df$Y_00)
      
      mon_no_s <- compute_monotonicity_bounds_no_S(df,monotonicity_constraint = "positive_M")
      
      bounds_df <- data.frame(
        beta_U = b,
        gamma_U = gamma_label,  # now a string, works with label_parsed
        VE_type = c("VE0", "VE1", "VEt"),
        Method = "Monotonicity",
        Lower = c(mon_no_s$VE0_lower, mon_no_s$VE1_lower, mon_no_s$VEt_lower),
        Upper = c(mon_no_s$VE0_upper, mon_no_s$VE1_upper, mon_no_s$VEt_upper),
        True = c(VE_0, VE_1, VE_T)
      )
      
      all_bounds[[paste(g, b, sep = "_")]] <- bounds_df
    }
  }
  
  
  # Combine all results
  plot_df <- bind_rows(all_bounds) %>% 
    filter(Lower<Upper)
  
  # Ensure beta_U is a factor with correct ordering
  plot_df <- plot_df %>%
    mutate(beta_U = factor(beta_U, levels = c(-log(5),-log(2), 0, log(2),log(5)))) %>% 
    mutate(Lower=100*Lower,
           Upper=100*Upper,
           True=100*True) %>% 
    mutate(VE_type = factor(
      VE_type,
      levels = c("VE0", "VE1", "VEt"),
      labels = c("paste('VE(', 0, ')')", "paste('VE(', 1, ')')", "VE[T]")
    )) %>% 
    # 1) turn gamma_U into a factor with the exact labels we want
    mutate(
      gamma_U = factor(
        gamma_U,
        levels = c("gamma[U]==-0.693147180559945", "gamma[U]==0.693147180559945"),
        labels = c("gamma[U]==-log(2)", "gamma[U]==log(2)")
      )
    ) 
 
  
  plot_df2<-plot_df
  df_ok   <- plot_df2 %>% filter(Lower >= 0)
  df_clip <- plot_df2 %>% filter(Lower < 0)
  pos    <- position_dodge(width = 0.6)
  
  
  
  
  # arrows at 0 (per method)
  arrows_below <- df_clip %>% mutate(y_start = 2, y_end = 0)
  #   mutate(y_start = 2, y_end = 0)
  
  df_true <- plot_df2 %>% distinct(VE_type, beta_U, True)
  
  cap_w <- 0.28     
  
  
  
  ggplot(plot_df2, aes(x = beta_U, group = VE_type)) +
    # regular intervals (both caps)
    geom_errorbar(
      data = df_ok,
      aes(ymin = pmax(Lower, 0), ymax = pmin(Upper, 100)),
      width = 0.24, position = pos
    ) +
    # vertical line for Lower < 0: from 0 up to upper
    geom_segment(
      data = df_clip,
      aes(x = beta_U, xend = beta_U, y = 0, yend = pmin(Upper, 100)),
      position = pos
    ) +
    # top cap at the upper end (ymin==ymax)
    geom_errorbar(
      data = df_clip,
      aes(ymin = pmin(Upper, 100), ymax = pmin(Upper, 100)),
      width = 0.24, position = pos
    ) +
    # down arrow at 0
    geom_segment(
      data = arrows_below,
      aes(x =beta_U, xend = beta_U, y = y_start, yend = y_end),
      position = pos,
      arrow = arrow(length = unit(3, "pt"), ends = "last", type = "closed")
    ) +
    # true value
    geom_point(aes(y = True), size = 1.3, position = pos) +
    facet_grid(gamma_U ~ VE_type, labeller = label_parsed) +
    labs(x = expression(gamma[U]), y = "True VE and estimated bounds (%)") +
    scale_x_discrete(
      labels = c("0" = "0", 
                 "-1.6094379124341"="-log(5)",
                 "1.6094379124341"="log(5)",
                 "-0.693147180559945" = "-log(2)", 
                 "0.693147180559945" = "log(2)")
    )+
    coord_cartesian(ylim = c(0, 100), expand = FALSE) +
    theme_minimal() +
    theme(strip.text = element_text(size = 12),
          axis.text.x = element_text(hjust = 0.5),
          panel.spacing.x = unit(0.8, "cm"),   # space between columns
          panel.spacing.y = unit(0.8, "cm")    # space between rows
          )
  

#LP bounds for Appendix
all_bounds_lp <- list()

# Outer loop over gamma_U
for (g in gamma_vals) {
  gamma_label <- paste0("gamma[U]==", g)  # for math expression in facet labels
  
  # Inner loop over beta_U
  for (b in beta_vals) {
    df <- simulate_data_PO_not_mon(1000000, beta_U = b, gamma_U = g)
    
    # True VE values
    VE_0 <- 1 - mean(df$Y_10) / mean(df$Y_00)
    VE_1 <- 1 - mean(df$Y_11) / mean(df$Y_01)
    VE_T <- 1 - mean(df$Y_11) / mean(df$Y_00)
    
    # LP bounds with monotonicity
    LP_M_no_S <- get_all_LP_bounds_summary_no_s(df, monotonicity_constraint = "positive_M")
    
    bounds_df <- data.frame(
      beta_U = b,
      gamma_U = gamma_label,  # character string for label_parsed
      VE_type = c("VE0", "VE1", "VEt"),
      Method = "LP",
      Lower = c(LP_M_no_S$VE0_lower, LP_M_no_S$VE1_lower, LP_M_no_S$VEt_lower),
      Upper = c(LP_M_no_S$VE0_upper, LP_M_no_S$VE1_upper, LP_M_no_S$VEt_upper),
      True = c(VE_0, VE_1, VE_T)
    )
    
    all_bounds_lp[[paste(g, b, sep = "_")]] <- bounds_df
  }
}

# Combine results
plot_df_lp <- bind_rows(all_bounds_lp)

# Ensure beta_U is a factor with correct ordering
plot_df_lp <- plot_df_lp %>%
  mutate(beta_U = factor(beta_U, levels = c(-5, -2, 0, 2, 5)))

# Plot
ggplot(plot_df_lp, aes(x = beta_U, group = VE_type)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = True),
             size = 1, position = position_dodge(width = 0.5)) +
  facet_grid(gamma_U ~ VE_type, labeller = label_parsed, scales = "free_y") +
  labs(
    x = expression(beta[U]),
    y = "Vaccine Efficacy (VE)",
#    title = "LP Bounds under Varying gamma[U] and beta[U]"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )


