library(dplyr)
library(ggplot2)

set.seed(12345)

folder_path <- "G:/My Drive/PhD3/data_example/code/functions"
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
  p_S_0 <- plogis(-1 + 1.5*0 + 1*U)
  p_S_1 <- plogis(-1 + 1.5*1 + 1*U)
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
  p_Y_00<-plogis(-2 - 1.7*0 + 1.5*S_0 + 0.5*0 + 1*U)
  p_Y_01<-plogis(-2 - 1.7*0 + 1.5*S_0 + 0.5*1 + 1*U)
  p_Y_10<-plogis(-2 - 1.7*1 + 1.5*S_1 + 0.5*0 + 1*U)
  p_Y_11<-plogis(-2 - 1.7*1 + 1.5*S_1 + 0.5*1 + 1*U)
  
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
beta_vals <- c(0, log(1.5),log(2))

# Loop over beta_S values
for (b in beta_vals) {
  df <- simulate_data_PO(1000000, beta_S = b)
  
  # True VE values
  VE_0 <- 1 - mean(df$Y_10) / mean(df$Y_00)
  VE_1 <- 1 - mean(df$Y_11) / mean(df$Y_01)
  VE_T <- 1 - mean(df$Y_11) / mean(df$Y_00)
  
  # Bounds from both methods
  mon_no_s <- compute_monotonicity_bounds_no_S(df)
  LP_M_no_S <- get_all_LP_bounds_summary_no_s(df, monotonicity_constraint = "positive_M")
  
  # Combine bounds into long-format data
  bounds_df <- bind_rows(
    data.frame(
      beta_S = b,
      VE_type = c("VE0", "VE1", "VEt"),
      Method = "Monotonicity",
      Lower = c(mon_no_s$VE0_lower, mon_no_s$VE1_lower, mon_no_s$VEt_lower),
      Upper = c(mon_no_s$VE0_upper, mon_no_s$VE1_upper, mon_no_s$VEt_upper),
      True = c(VE_0, VE_1, VE_T)
    ),
    data.frame(
      beta_S = b,
      VE_type = c("VE0", "VE1", "VEt"),
      Method = "LP",
      Lower = c(LP_M_no_S$VE0_lower, LP_M_no_S$VE1_lower, LP_M_no_S$VEt_lower),
      Upper = c(LP_M_no_S$VE0_upper, LP_M_no_S$VE1_upper, LP_M_no_S$VEt_upper),
      True = c(VE_0, VE_1, VE_T)
    )
  )
  
  all_bounds[[as.character(b)]] <- bounds_df
}

# Combine all into one dataframe
plot_df <- bind_rows(all_bounds)
#in precentaeges
plot_df_p<-plot_df
plot_df_p$Lower<-100*plot_df_p$Lower
plot_df_p$Upper<-100*plot_df_p$Upper
plot_df_p$True<-100*plot_df_p$True
# Plot
#library(ggplot2)
#library(dplyr)

# Treat beta_S as a categorical variable with the correct order
plot_df_p <- plot_df_p %>%
  mutate(beta_S = factor(beta_S, levels = c(0, log(1.5),log(2)))) %>% 
  mutate(VE_type = factor(
    VE_type,
    levels = c("VE0", "VE1", "VEt"),
    labels = c("paste('VE(', 0, ')')", "paste('VE(', 1, ')')", "VE[T]")
  ))

plot_df<- plot_df%>%
  mutate(beta_S = factor(beta_S, levels = c(0, log(1.5),log(2)))) %>% 
  mutate(VE_type = factor(
    VE_type,
    levels = c("VE0", "VE1", "VEt"),
    labels = c("paste('VE(', 0, ')')", "paste('VE(', 1, ')')", "VE[T]")
  ))
#%>% 
  #only for summary
#filter(VE_type=="VE0")
# Plot: one panel per Method, grouped by VE_type and beta_S (categorical x-axis)
plot_df %>% 
 # filter(VE_type%in%c("paste('VE(', 0, ')')", "paste('VE(', 1, ')')")) %>% 
  ggplot( aes(x = beta_S, group = VE_type)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = True),
             size = 1, position = position_dodge(width = 0.5)) +
  # geom_hline(yintercept = c(-50, 50), linetype = "dashed", color = "red", size = 0.3) +
  facet_grid(Method ~ VE_type, labeller = label_parsed) +
  scale_x_discrete(
    labels = c("0" = "0", 
               "0.405465108108164" = "log(1.5)", 
               "0.693147180559945" = "log(2)")
  ) +
  # scale_y_continuous(
  #  limits = c(0, 50),
  # oob = scales::squish,
  #  name = "VE (%)"
  # ) +
  labs(x = expression(beta[S]),y = "True VE and Estimated Bounds") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.5))



plot_df_p %>%
  filter(VE_type %in% c("paste('VE(', 0, ')')", "paste('VE(', 1, ')')")) %>%
  ggplot(aes(x = beta_S, group = VE_type)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(y = True),
             size = 2.5,  # ⬅️ Increased point size
             shape = 19,  # solid circle
             position = position_dodge(width = 0.5)) +
  facet_grid(Method ~ VE_type, labeller = label_parsed) +
  scale_x_discrete(
    labels = c("0" = "0",
               "0.405465108108164" = "log(1.5)",
               "0.693147180559945" = "log(2)")
  ) +
  labs(
    x = expression(beta[S]),
    y = "True VE and Estimated Bounds (%)"
  ) +
  theme_minimal(base_size = 14) +  # ⬅️ Increase base font size for all text
  theme(
    strip.text = element_text(size = 16, face = "bold"),  # Facet titles
    axis.text.x = element_text(size = 13, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold")
  )


ggsave("G:\\My Drive\\PhD3\\data_example\\Figures\\ve_bounds_plot_blinding_effect_log_per.pdf", width = 12, height = 6, dpi = 300)


# ggsave("G:\\My Drive\\PhD3\\data_example\\Figures\\ve_bounds_plot_blinding_effect.pdf", width = 6.5, height = 3, units = "in")
# ggsave("G:\\My Drive\\PhD3\\data_example\\Figures\\ve_bounds_plot_blinding_effect_summary.pdf", width = 6.5, height = 3, units = "in")
# 
# 
# 
# 
# ggsave("G:\\My Drive\\PhD3\\data_example\\Figures\\ve_bounds_plot_blinding_effect_log.pdf", width = 6.5, height = 3, units = "in")

ggsave("G:\\My Drive\\PhD3\\data_example\\Figures\\ve_bounds_plot_blinding_effect_log_per_full.pdf", width = 12, height = 6, dpi = 300)





