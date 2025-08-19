library(mvtnorm)


rho_vals <- seq(0, 1, by = 0.001)

# Function to compute choice probability for car
compute_car_prob <- function(rho) {
  # Define mean vector for utility differences
  mean_vec <- c(0, 0)  
  
  # Covariance matrix for utility differences:
  # Define V1 = ε_RedBus - ε_Car
  #       V2 = ε_BlueBus - ε_Car
  # So,
  # Var(V1) = Var(ε_RedBus) + Var(ε_Car) = 1 + 1 = 2
  # Var(V2) = Var(ε_BlueBus) + Var(ε_Car) = 1 + 1 = 2
  # Cov(V1, V2) = Cov(ε_RedBus, ε_BlueBus) + Cov(ε_Car, ε_Car) = rho + 1 = rho + 1
  
  sigma_mat <- matrix(c(2, rho + 1,
                        rho + 1, 2), nrow = 2)
  
  # We compute: P(ε_Car > ε_RedBus and ε_Car > ε_BlueBus)
  #             P(ε_RedBus - ε_Car < 0 and ε_BlueBus - ε_Car < 0)
  #             P(V1 < 0, V2 < 0) = CDF of bivariate normal at (0,0)
  prob <- pmvnorm(upper = c(0, 0), mean = mean_vec, sigma = sigma_mat)
  
  return(as.numeric(prob))
}

# apply to each rho
results <- data.frame(
  rho = rho_vals,
  car_prob = sapply(rho_vals, compute_car_prob)
)



# graphics ----------------------------------------------------------------

library(ggplot2)

ggplot(results, aes(x = rho, y = car_prob)) +
  geom_line(color = "black", size = 1.2) +
  scale_x_continuous(
    name = expression("Correlation" ~ rho ~ "between Red Bus and Blue Bus"),
    breaks = seq(0, 1, by = 0.25)
  ) +
  scale_y_continuous(
    name = expression("Prob(Car)"),
    limits = c(0, 0.5),
    breaks = seq(0, 1, by = 0.1)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
  ) +
  geom_segment(x = 0, y = 0, xend = 0, yend = 0.52, 
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
             color = "black") +
  geom_segment(aes(x = 0, y = 0, xend = 1.05, yend = 0), 
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               color = "black")

