library(ggplot2)
library(dplyr)
library(tidyverse)

data = read.csv(paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"))

plot_omega = c(0.01, 0.04, 1, 16, 25)
for (i in 1:length(plot_omega)){
filtered_data = data %>% filter(beta_free == 1, omega_var == plot_omega[i], rho == 0, Tfull == 2,
                                n_reps == 500) %>% filter(!(N == 20))

no_omega_filter = data %>% filter(beta_free == 1, rho == 0, Tfull == 2,
                                n_reps == 500) %>% filter(!(N == 20)) %>% filter(omega_error < 5)


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
  coord_cartesian(ylim = c(min(c(no_omega_filter$beta_error,
             no_omega_filter$omega_error,
             no_omega_filter$rho_error)),
       max(c(no_omega_filter$beta_error,
             no_omega_filter$omega_error,
             no_omega_filter$rho_error)))) +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.85)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  scale_x_continuous(breaks = seq(0, max(filtered_data$N, na.rm = TRUE), by = 50))

print(p)
}


# -------------------------------------------------------------------------
library(purrr)
library(future.apply)

get_estimate_results <- function(n_reps, N){
  
  get_estimates <- function(N){
    
    create_df <- function(N){
      c = sample(x = c(2, 4, 6), size = N, replace = T)
      x = c(rep.int(NA, N))
      for (i in 1:length(x)){
        x[i] = rnorm(N, mean = c[i], sd = 4)
      }
      y = 3*x + rnorm(N, mean = 0, sd = 3)
      return(df = data.frame(y, x, c))
    }
    
    df = create_df(N)
    lm = lm(y ~ x, data = df)
    beta = c(beta_0 = coef(lm)[1], beta_1 = coef(lm)[2])
    names(beta) = c("beta_0", "beta_1")
    sd = sqrt(diag(vcov(lm))) # structure may not be good enough here
    names(sd) = c("sd_0", "sd_1")
    res = list(beta = beta, sd = sd)
    
    return(res)
  }
  
  print("Starting true var calcs")
  
  true_var = foo$sd 
  
  print("Finished true var calcs")
  
  plan(multisession)
  
  sim_replications <- future_replicate(n_reps, get_estimates(N = N), # check if N = N makes a difference
    simplify = FALSE) 
  
  # Get estimates
  beta_0_estimates <- sapply(sim_replications, function(x) x[[1]]["beta_0"])
  beta_1_estimates <- sapply(sim_replications, function(x) x[[1]]["beta_1"])
  
  true_parameters <- c(beta_0 = 0, beta_1 = 3)
  
  # Get CI-Overlap info by first adjusting the true var to the appropriate sd 
  # for current N and not use the one based on large N from quasi-true variance
  adjust_var_to_N <- function(unadjusted_true_var, N_new, N_old){
    return((N_old*unadjusted_true_var)/N_new)
  }
  
  adjusted_true_var <- adjust_var_to_N(true_var, N_new = N, N_old = 20000)
  alpha = 0.05
  
  if(alpha > 0 & alpha < 1){
    z_value = qnorm(1-alpha/2) # For different CIs
  } else {
    print("Error: Alpha value invalid")
  }
  
  ci_lower <- true_parameters - z_value * sqrt(adjusted_true_var)
  ci_lower <- unlist(ci_lower)
  
  ci_upper <- true_parameters + z_value * sqrt(adjusted_true_var)
  ci_upper <- unlist(ci_upper)
  
  # Check cuccess or failure for each parameter
  success_beta_0 <- sum(beta_0_estimates >= ci_lower[1] & beta_0_estimates <= ci_upper[1])
  failure_beta_0 <- length(beta_0_estimates) - success_beta_0
  
  success_beta_1 <- sum(beta_1_estimates >= ci_lower[2] & beta_1_estimates <= ci_upper[2])
  failure_beta_1 <- length(beta_1_estimates) - success_beta_1
  
  # Obtain uncertainty (upper and lower) of CI-overlap result by normal approximation
  eps = 1e-10
  p_0 = success_beta_0/length(beta_0_estimates)
  norm_uncertainty_upper_b0 = p_0 + z_value * sqrt( (p_0*(1-p_0) + eps)/N )
  norm_uncertainty_lower_b0 = p_0 - z_value * sqrt( (p_0*(1-p_0) + eps)/N )
  
  p_1 = success_beta_1/length(beta_1_estimates)
  norm_uncertainty_upper_b1 = p_1 + z_value * sqrt( (p_1*(1-p_1) + eps)/N )
  norm_uncertainty_lower_b1 = p_1 - z_value * sqrt( (p_1*(1-p_1) + eps)/N )
  
  # Now, do the same, but to obtain exact ones (by clopper and pearson)
  p_0_ub = 1 - qbeta(p = alpha/2, shape1 = n_reps - success_beta_0, shape2 = success_beta_0 + 1)
  p_0_lb = 1 - qbeta(p = 1 - alpha/2, shape1 = n_reps - success_beta_0 + 1, shape2 = success_beta_0)
  
  p_1_ub = 1 - qbeta(p = alpha/2, shape1 = n_reps - success_beta_1, shape2 = success_beta_1 + 1)
  p_1_lb = 1 - qbeta(p = 1 - alpha/2, shape1 = n_reps - success_beta_1 + 1, shape2 = success_beta_1)
  
  # Compute Wilson Score Interval
  wilson_ci <- function(p, n, z) {
    denominator <- 1 + (z^2 / n)
    center <- p + (z^2 / (2 * n))
    margin <- z * sqrt((p * (1 - p) / n) + (z^2 / (4 * n^2)))
    
    lower <- (center - margin) / denominator
    upper <- (center + margin) / denominator
    
    return(c(lower, upper))
  }
  
  # Compute Wilson CI for beta_0 and beta_1
  wilson_beta_0 <- wilson_ci(p_0, n_reps, z_value)
  wilson_beta_1 <- wilson_ci(p_1, n_reps, z_value)
  
  # Add to results table
  results_table <- data.frame(
    Param = c("beta_0", "beta_1"),
    n_reps = n_reps,
    Success_Count = c(success_beta_0, success_beta_1),
    Success_rate = c(p_0, p_1),
    CI_Lower = ci_lower,
    True_param = true_parameters,
    CI_Upper = ci_upper,
    alpha = alpha,
    N = N,
    CI_norm_Uncertainty_upper = c(norm_uncertainty_upper_b0, norm_uncertainty_upper_b1),
    CI_norm_uncertainty_lower = c(norm_uncertainty_lower_b0, norm_uncertainty_lower_b1),
    CI_ex_Uncertainty_upper = c(p_0_ub, p_1_ub),
    CI_ex_uncertainty_lower = c(p_0_lb, p_1_lb),
    CI_wilson_upper = c(wilson_beta_0[2], wilson_beta_1[2]),
    CI_wilson_lower = c(wilson_beta_0[1], wilson_beta_1[1]),
    row.names = NULL
  )
  
  return(results_table)
}

N = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 400, 500, 1000)
n_reps = c(1000)
set.seed(1) 

create_grid <- function(...) {
  vars <- list(...)
  var_names <- sapply(substitute(list(...))[-1], deparse)
  setNames(expand.grid(vars), var_names)
}

(parameters = create_grid(N, n_reps))
(foo = get_estimates(N = 20000))
results <- pmap(parameters, get_estimate_results) # see of plan(multisession) outside is more efficient + produces same results as before under seed
combined_df = do.call(rbind, results)

p <- ggplot(combined_df, aes(x = N, color = Param)) + 
  geom_line(aes(y = Success_rate)) +
  geom_point(aes(y = Success_rate), size = 1) +
  labs(title = "Success Rate of estimate lying within 95%-CI",
       x = "Sample Size (N)",
       y = "Within CI Success Rate") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(legend.position = c(0.15, 0.15)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black") +
  scale_x_continuous(breaks = seq(from = 0, to = 500, by = 50))

print(p)


# -------------------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(ggplot2)

df = read.csv(paste0(getwd(), "/GitHub/MasterThesis/DGP1_3000reps_Results.csv"))
df <- df %>% filter(rho == 0,
         Tfull == 2,
         n_reps == 500,
         omega_var == 1)
df <- df %>%
  mutate(beta_lower = beta_error - 1.96 * beta_estimate_sd,
         beta_upper = beta_error + 1.96 * beta_estimate_sd,
         omega_lower = omega_error - 1.96 * omega_estimate_sd,
         omega_upper = omega_error + 1.96 * omega_estimate_sd,
         rho_lower = rho_error - 1.96 * rho_estimate_sd,
         rho_upper = rho_error + 1.96 * rho_estimate_sd)
foo = df
foo <- foo %>% filter(beta_error < 10, beta_error > -10, omega_error < 10,
                      omega_error > -10, rho_error < 10, rho_error > -10)
min = min( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
max = max( c(foo$beta_error, foo$omega_error, foo$rho_error), na.rm = TRUE )
p <- ggplot(df, aes(x = N)) + list(
  geom_ribbon(aes(ymin = beta_lower, ymax = beta_upper, fill = "Beta CI"), alpha = 0.3),
  geom_ribbon(aes(ymin = omega_lower, ymax = omega_upper, fill = "Omega CI"), alpha = 0.3),
  geom_ribbon(aes(ymin = rho_lower, ymax = rho_upper, fill = "Rho CI"), alpha = 0.3),
  geom_line(aes(y = beta_error, color = "Beta Error")),
  geom_line(aes(y = omega_error, color = "Omega Error")),
  geom_line(aes(y = rho_error, color = "Rho Error")),
  geom_point(aes(y = beta_error, color = "Beta Error"), size = 1),
  geom_point(aes(y = omega_error, color = "Omega Error"), size = 1),
  geom_point(aes(y = rho_error, color = "Rho Error"), size = 1)) +
  labs(title = "title_text",
       x = "Sample Size (N)",
       y = "Value",,
       fill = "Confidence Interval") +
  coord_cartesian(ylim = c(min, max)) +
  theme_minimal() +
  theme(legend.position = c(0.85, 0.85)) +
  scale_fill_manual(values = c("Beta CI" = "#E69F00",   # Orange
                               "Omega CI" = "#009900",  # Purple
                               "Rho CI" = "#0066ff"))   # Teal
print(p)  
  
  
  
  
  
  
  
  