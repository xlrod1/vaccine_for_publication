library(dplyr)
library(ggplot2)

set.seed(12345)

folder_path <- "G:/My Drive/PhD3/data_example/code/functions"
# List all R files in the folder
r_files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
# Source each file
sapply(r_files, source)

simulate_data_PO <- function(n,beta_S,gamma_S) {
  # Randomized A and U
  A <- rbinom(n, 1, 0.5)
  U <- rbinom(n, 1, 0.5)
  
  # Parameters (to be tuned for constraints)
  # S depends on A and U
  p_S_0 <- plogis(-1 + 1.5*0 + beta_S*U)
  p_S_1 <- plogis(-1 + 1.5*1 + beta_S*U)
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
  p_Y_00<-plogis(-2 - 1.7*0 + gamma_S*S_0 + 0.5*0 + 1*U)
  p_Y_01<-plogis(-2 - 1.7*0 + gamma_S*S_0 + 0.5*1 + 1*U)
  p_Y_10<-plogis(-2 - 1.7*1 + gamma_S*S_1 + 0.5*0 + 1*U)
  p_Y_11<-plogis(-2 - 1.7*1 + gamma_S*S_1 + 0.5*1 + 1*U)
  
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


# 1. Store results for all beta_S
all_bounds <- list()
beta_vals <- c(-2,0, 2, 5, 10)

for (b in beta_vals) {
  df <- simulate_data_PO(1000000, beta_S = b,gamma_S=0)
  
  # True VE values
  VE_0 <- 1 - mean(df$Y_10) / mean(df$Y_00)
  VE_1 <- 1 - mean(df$Y_11) / mean(df$Y_01)
  VE_T <- 1 - mean(df$Y_11) / mean(df$Y_00)
  
  # Bounds with and without S
  mon       <- compute_monotonicity_bounds(df)
  LP_M      <- get_all_LP_bounds_summary(df, monotonicity_constraint = "positive_M")
  mon_no_s  <- compute_monotonicity_bounds_no_S(df)
  LP_M_no_S <- get_all_LP_bounds_summary_no_s(df, monotonicity_constraint = "positive_M")
  
  # Combine all four methods
  bounds_df <- bind_rows(
    data.frame(
      beta_S = b,
      VE_type = c("VE0", "VE1", "VEt"),
      Method = "Monotonicity",
      Use_S = "With S",
      Lower = c(mon$VE0_lower, mon$VE1_lower, mon$VEt_lower),
      Upper = c(mon$VE0_upper, mon$VE1_upper, mon$VEt_upper),
      True  = c(VE_0, VE_1, VE_T)
    ),
    data.frame(
      beta_S = b,
      VE_type = c("VE0", "VE1", "VEt"),
      Method = "LP",
      Use_S = "With S",
      Lower = c(LP_M$VE0_lower, LP_M$VE1_lower, LP_M$VEt_lower),
      Upper = c(LP_M$VE0_upper, LP_M$VE1_upper, LP_M$VEt_upper),
      True  = c(VE_0, VE_1, VE_T)
    ),
    data.frame(
      beta_S = b,
      VE_type = c("VE0", "VE1", "VEt"),
      Method = "Monotonicity",
      Use_S = "No S",
      Lower = c(mon_no_s$VE0_lower, mon_no_s$VE1_lower, mon_no_s$VEt_lower),
      Upper = c(mon_no_s$VE0_upper, mon_no_s$VE1_upper, mon_no_s$VEt_upper),
      True  = c(VE_0, VE_1, VE_T)
    ),
    data.frame(
      beta_S = b,
      VE_type = c("VE0", "VE1", "VEt"),
      Method = "LP",
      Use_S = "No S",
      Lower = c(LP_M_no_S$VE0_lower, LP_M_no_S$VE1_lower, LP_M_no_S$VEt_lower),
      Upper = c(LP_M_no_S$VE0_upper, LP_M_no_S$VE1_upper, LP_M_no_S$VEt_upper),
      True  = c(VE_0, VE_1, VE_T)
    )
  )
  
  all_bounds[[as.character(b)]] <- bounds_df
}

# 2. Combine and format for plotting
plot_df <- bind_rows(all_bounds) %>%
  mutate(
    beta_S = factor(beta_S, levels = c(0, 2, 5, 10)),
    Method = factor(Method, levels = c("LP", "Monotonicity")),
    Use_S = factor(Use_S, levels = c("With S", "No S"))
  )

# 3. Plot
plot_df_mon<-plot_df %>% filter(beta_S%in%c("0","2","5","10"),Use_S!="No S",Method!="LP") %>% 
  filter(Lower<Upper)
ggplot(plot_df_mon, aes(x = beta_S)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Use_S),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = True, color = Use_S),
             size = 1, position = position_dodge(width = 0.5)) +
  facet_grid(VE_type ~ Method, scales = "free_y") +
  labs(
    #title = "VE Bounds by Method, Use of S, and β_S",
       x = expression(beta[S]),
       y = "Vaccine Efficacy (VE)") +
  scale_color_manual(values = c("With S" = "#1f78b4", "No S" = "#33a02c")) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("G:\\My Drive\\PhD3\\data_example\\Figures\\app_ve_bounds_plot_SU_assump_mon.pdf", width = 6.5, height = 3, units = "in")



plot_df_LP<-plot_df %>% filter(beta_S%in%c("0","2","5","10"),Use_S!="No S",Method=="LP") %>% 
  filter(Lower<Upper)
ggplot(plot_df_LP, aes(x = beta_S)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Use_S),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = True, color = Use_S),
             size = 1, position = position_dodge(width = 0.5)) +
  facet_grid(VE_type ~ Method, scales = "free_y") +
  labs(
    #title = "VE Bounds by Method, Use of S, and β_S",
       x = expression(beta[S]),
       y = "Vaccine Efficacy (VE)") +
  scale_color_manual(values = c("With S" = "#1f78b4", "No S" = "#33a02c")) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank())

ggsave("G:\\My Drive\\PhD3\\data_example\\Figures\\app_ve_bounds_plot_SU_assump_LP.pdf", width = 6.5, height = 3, units = "in")






###################################When gamma_S>0
# 1. Store results for all beta_S
all_bounds <- list()
beta_vals <- c(-2,0, 2, 5, 10)

for (b in beta_vals) {
  df <- simulate_data_PO(1000000, beta_S = b,gamma_S=2)
  
  # True VE values
  VE_0 <- 1 - mean(df$Y_10) / mean(df$Y_00)
  VE_1 <- 1 - mean(df$Y_11) / mean(df$Y_01)
  VE_T <- 1 - mean(df$Y_11) / mean(df$Y_00)
  
  # Bounds with and without S
  mon       <- compute_monotonicity_bounds(df)
  LP_M      <- get_all_LP_bounds_summary(df, monotonicity_constraint = "positive_M")
  mon_no_s  <- compute_monotonicity_bounds_no_S(df)
  LP_M_no_S <- get_all_LP_bounds_summary_no_s(df, monotonicity_constraint = "positive_M")
  
  # Combine all four methods
  bounds_df <- bind_rows(
    data.frame(
      beta_S = b,
      VE_type = c("VE0", "VE1", "VEt"),
      Method = "Monotonicity",
      Use_S = "With S",
      Lower = c(mon$VE0_lower, mon$VE1_lower, mon$VEt_lower),
      Upper = c(mon$VE0_upper, mon$VE1_upper, mon$VEt_upper),
      True  = c(VE_0, VE_1, VE_T)
    ),
    data.frame(
      beta_S = b,
      VE_type = c("VE0", "VE1", "VEt"),
      Method = "LP",
      Use_S = "With S",
      Lower = c(LP_M$VE0_lower, LP_M$VE1_lower, LP_M$VEt_lower),
      Upper = c(LP_M$VE0_upper, LP_M$VE1_upper, LP_M$VEt_upper),
      True  = c(VE_0, VE_1, VE_T)
    ),
    data.frame(
      beta_S = b,
      VE_type = c("VE0", "VE1", "VEt"),
      Method = "Monotonicity",
      Use_S = "No S",
      Lower = c(mon_no_s$VE0_lower, mon_no_s$VE1_lower, mon_no_s$VEt_lower),
      Upper = c(mon_no_s$VE0_upper, mon_no_s$VE1_upper, mon_no_s$VEt_upper),
      True  = c(VE_0, VE_1, VE_T)
    ),
    data.frame(
      beta_S = b,
      VE_type = c("VE0", "VE1", "VEt"),
      Method = "LP",
      Use_S = "No S",
      Lower = c(LP_M_no_S$VE0_lower, LP_M_no_S$VE1_lower, LP_M_no_S$VEt_lower),
      Upper = c(LP_M_no_S$VE0_upper, LP_M_no_S$VE1_upper, LP_M_no_S$VEt_upper),
      True  = c(VE_0, VE_1, VE_T)
    )
  )
  
  all_bounds[[as.character(b)]] <- bounds_df
}

# 2. Combine and format for plotting
plot_df <- bind_rows(all_bounds) %>%
  mutate(
    beta_S = factor(beta_S, levels = c(0, 2, 5, 10)),
    Method = factor(Method, levels = c("LP", "Monotonicity")),
    Use_S = factor(Use_S, levels = c("With S", "No S"))
  )

# 3. Plot
plot_df_mon<-plot_df %>% filter(beta_S%in%c("0","2","5","10"),Use_S!="No S",Method!="LP") %>% 
  filter(Lower<Upper)
ggplot(plot_df_mon, aes(x = beta_S)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Use_S),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = True, color = Use_S),
             size = 1, position = position_dodge(width = 0.5)) +
  facet_grid(VE_type ~ Method, scales = "free_y") +
  labs(
    #title = "VE Bounds by Method, Use of S, and β_S",
    x = expression(beta[S]),
    y = "Vaccine Efficacy (VE)") +
  scale_color_manual(values = c("With S" = "#1f78b4", "No S" = "#33a02c")) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank())


ggsave("G:\\My Drive\\PhD3\\data_example\\Figures\\app_ve_bounds_plot_SU_SY_assump_mon.pdf", width = 6.5, height = 3, units = "in")



plot_df_LP<-plot_df %>% filter(beta_S%in%c("0","2","5","10"),Use_S!="No S",Method=="LP") %>% 
  filter(Lower<Upper)
ggplot(plot_df_LP, aes(x = beta_S)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Use_S),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = True, color = Use_S),
             size = 1, position = position_dodge(width = 0.5)) +
  facet_grid(VE_type ~ Method, scales = "free_y") +
  labs(
    #title = "VE Bounds by Method, Use of S, and β_S",
    x = expression(beta[S]),
    y = "Vaccine Efficacy (VE)") +
  scale_color_manual(values = c("With S" = "#1f78b4", "No S" = "#33a02c")) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank())

ggsave("G:\\My Drive\\PhD3\\data_example\\Figures\\app_ve_bounds_plot_SU_SY_assump_LP.pdf", width = 6.5, height = 3, units = "in")



