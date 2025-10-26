# --- Packages ---
library(dplyr)
library(ggplot2)
library(grid)  # for unit() used in arrow()

set.seed(12345)

# --- Load your helper functions from the folder ---
#folder_path <- "your path of functions"
r_files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
sapply(r_files, source)

# --- Simulator (PO = potential outcomes) ---
simulate_data_PO <- function(n, beta_S) {
  # Randomized A and U
  A <- rbinom(n, 1, 0.5)
  U <- rbinom(n, 1, 0.5)
  
  # S depends on A and U
  p_S_0 <- plogis(-1 + 1.5*0 + 1*U)
  p_S_1 <- plogis(-1 + 1.5*1 + 1*U)
  S_0 <- rbinom(n, 1, p_S_0)
  S_1 <- rbinom(n, 1, p_S_1)
  S <- ifelse(A == 1, S_1, S_0)
  
  # B depends on S and U (increasing in U); beta_S controls S→B strength
  p_B_0 <- plogis(-1 + beta_S * S_0 + 1.5*U)
  p_B_1 <- plogis(-1 + beta_S * S_1 + 1.5*U)
  B_0 <- rbinom(n, 1, p_B_0)
  B_1 <- rbinom(n, 1, p_B_1)
  B   <- ifelse(A == 1, B_1, B_0)
  
  # Create all potential outcomes Y^{a,m}
  # Assumption used: Y^{a,m=-1} | B=m = Y^{a,m} | B=m
  p_Y_00 <- plogis(-2 - 1.7*0 + 1.5*S_0 + 0.5*0 + 1*U)
  p_Y_01 <- plogis(-2 - 1.7*0 + 1.5*S_0 + 0.5*1 + 1*U)
  p_Y_10 <- plogis(-2 - 1.7*1 + 1.5*S_1 + 0.5*0 + 1*U)
  p_Y_11 <- plogis(-2 - 1.7*1 + 1.5*S_1 + 0.5*1 + 1*U)
  
  Y_00 <- rbinom(n, 1, p_Y_00)
  Y_01 <- rbinom(n, 1, p_Y_01)
  Y_10 <- rbinom(n, 1, p_Y_10)
  Y_11 <- rbinom(n, 1, p_Y_11)
  
  # Observed outcome by (A,B)
  df <- data.frame(A, U, S, B, Y_00, Y_01, Y_10, Y_11) %>%
    mutate(Y = dplyr::case_when(
      A == 0 & B == 0 ~ Y_00,
      A == 0 & B == 1 ~ Y_01,
      A == 1 & B == 0 ~ Y_10,
      A == 1 & B == 1 ~ Y_11
    ))
  
  return(df)
}

# 1) Investigate how beta_S affects bounds’ width
all_bounds <- list()
beta_vals  <- c(0, log(1.5), log(2))

for (b in beta_vals) {
  df <- simulate_data_PO(1000000, beta_S = b)
  
  # True VE values
  VE_0 <- 1 - mean(df$Y_10) / mean(df$Y_00)
  VE_1 <- 1 - mean(df$Y_11) / mean(df$Y_01)
  VE_T <- 1 - mean(df$Y_11) / mean(df$Y_00)
  
  # Bounds from both methods
  mon_no_s <- compute_monotonicity_bounds_no_S(df, monotonicity_constraint = "positive_M")
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

# Combine results from all beta_S values
plot_df <- dplyr::bind_rows(all_bounds)

# Percentage scale for plotting
plot_df_p <- plot_df %>%
  mutate(
    Lower = 100 * Lower,
    Upper = 100 * Upper,
    True  = 100 * True
  )

# Order factors for facets/axes
plot_df_p <- plot_df_p %>%
  mutate(
    beta_S  = factor(beta_S, levels = c(0, log(1.5), log(2))),
    VE_type = factor(
      VE_type,
      levels = c("VE0", "VE1", "VEt"),
      labels = c("paste('VE(', 0, ')')", "paste('VE(', 1, ')')", "VE[T]")
    )
  )

# (Optional) keep a non-percentage copy with the same factor setup
plot_df <- plot_df %>%
  mutate(
    beta_S  = factor(beta_S, levels = c(0, log(1.5), log(2))),
    VE_type = factor(
      VE_type,
      levels = c("VE0", "VE1", "VEt"),
      labels = c("paste('VE(', 0, ')')", "paste('VE(', 1, ')')", "VE[T]")
    )
  )

# Split rows where lower < 0 to draw arrows from 0 downwards
plot_df2 <- plot_df_p
df_ok    <- plot_df2 %>% filter(Lower >= 0)
df_clip  <- plot_df2 %>% filter(Lower < 0)
pos      <- position_dodge(width = 0.6)

# Arrow segments for the “below 0” cases
arrows_below <- df_clip %>% mutate(y_start = 2, y_end = 0)

# (For over-plotting a dot at the true VE)
df_true <- plot_df2 %>% distinct(VE_type, beta_S, True)

cap_w <- 0.28

p <- ggplot(plot_df2, aes(x = beta_S, group = VE_type, type = Method)) +
  # regular intervals (both caps) when Lower >= 0
  geom_errorbar(
    data = df_ok,
    aes(ymin = pmax(Lower, 0), ymax = pmin(Upper, 100)),
    width = 0.24, position = pos
  ) +
  # vertical line when Lower < 0: from 0 up to upper
  geom_segment(
    data = df_clip,
    aes(x = beta_S, xend = beta_S, y = 0, yend = pmin(Upper, 100)),
    position = pos
  ) +
  # top cap at the upper end for the clipped rows
  geom_errorbar(
    data = df_clip,
    aes(ymin = pmin(Upper, 100), ymax = pmin(Upper, 100)),
    width = 0.24, position = pos
  ) +
  # down arrow indicating “extends below 0”
  geom_segment(
    data = arrows_below,
    aes(x = beta_S, xend = beta_S, y = y_start, yend = y_end),
    position = pos,
    arrow = arrow(length = unit(3, "pt"), ends = "last", type = "closed")
  ) +
  # True VE (percent)
  geom_point(aes(y = True), size = 1.3, position = pos) +
  facet_grid(VE_type ~ Method, labeller = label_parsed) +
  labs(x = expression(beta[S]), y = "True VE and estimated bounds (%)") +
  scale_x_discrete(
    labels = c("0" = "0",
               "0.405465108108164" = "log(1.5)",
               "0.693147180559945" = "log(2)")
  ) +
  coord_cartesian(ylim = c(0, 100), expand = FALSE) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(hjust = 0.5)
  )

print(p)

# (Optional) save
# ggsave("bounds_vs_betaS.png", p, width = 8, height = 6, dpi = 300)
