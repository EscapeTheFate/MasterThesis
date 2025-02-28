library(ggplot2)
library(dplyr)
library(tidyverse)

data = read.csv(paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"))

filtered_data = data %>% filter(beta_free == 1, omega_var == 0.01, rho == 0, Tfull == 2,
                                n_reps == 500) %>% filter(!(N == 20))

no_omega_filter = data %>% filter(beta_free == 1, rho == 0, Tfull == 2,
                                n_reps == 500) %>% filter(!(N == 20))


p <- ggplot(filtered_data, aes(x = N)) + 
  list(
  geom_line(aes(y = beta_error, color = "Beta Error")),
  geom_line(aes(y = omega_error, color = "Omega Error")),
  geom_line(aes(y = rho_error, color = "Rho Error")),
  geom_point(aes(y = beta_error, color = "Beta Error"), size = 1),
  geom_point(aes(y = omega_error, color = "Omega Error"), size = 1),
  geom_point(aes(y = rho_error, color = "Rho Error"), size = 1)) +
  labs(title = "Whatever",
       x = "Sample Size (N)",
       y = "Value",
       fill = "Confidence Interval") +
  ylim(min(c(no_omega_filter$beta_error,
             no_omega_filter$omega_error,
             no_omega_filter$rho_error)),
       max(c(no_omega_filter$beta_error,
             no_omega_filter$omega_error,
             no_omega_filter$rho_error))) +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.85)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  scale_x_continuous(breaks = seq(0, max(filtered_data$N, na.rm = TRUE), by = 50))

print(p)
