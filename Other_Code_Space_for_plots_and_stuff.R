
# The trash/space inbetween different code environments

# reworked simulated dataframes are slower than original one
draw_simulation_data_reworked <- function(Tfull, N, quest, rho, system){
  
  time <- rep(1:Tfull, each = quest)
  Tp <- length(time)
  latent_means <- c(2, 4, 7)
  
  # Setup dataframe by preallocating list structure
  data_list <- vector("list", N)
  
  for (j in 1:N) {
    latent_class <- sample(1:3, 1, prob = c(1/3, 1/3, 1/3))
    
    X = matrix(c(rnorm(n = Tp, mean = latent_means[latent_class], sd = 1),
                 rnorm(n = Tp, mean = 0, sd = 1)),
               nrow = Tp, ncol = 2)
    
    LG <- t(chol(system$GammaT))
    LO <- t(chol(system$Omega))
    lRE = nrow(LO)
    
    gam = LO %*% rnorm(lRE)
    
    U = matrix(0, nrow = Tp)
    U = X %*% system$beta + X[, c(1:lRE)] %*% gam + LG %*% rnorm(dim(LG)[2])
    
    y = U*0
    y <- rowSums(outer(U, system$tauk, ">")) + 1
    
    data_list[[j]] <- data.table(ID = j, X1 = X[,1], X2 = X[,2], y = y, time = time, quest = rep(1:quest, Tfull))
  }
  
  return(rbindlist(data_list))
}


# -------------------------------------------------------------------------
# First out-of-code comparison of two different dataframe creations


library(data.table)
library(Rprobit)
library(tictoc)
library(purrr)

for (z in 1:2){
  set.seed(1)
  Tfull = 3
  N = 200
  quest = 1
  beta_free = 1
  rho = 0.5
  omega_var = 1
  return_val = ""
  timer = F
  DontSkipFit = T
  
  mod <- mod_AR_cl$new(alt  = 7,
                       Hb   = diag(2)[,-2,drop=FALSE],
                       fb   = matrix(0,2,1),
                       HO   = matrix(1,1,1),
                       fO   = matrix(0,1,1),
                       HL   = matrix(0,1,0),
                       fL   = matrix(1,1,1),
                       ordered = TRUE
  )
  
  mod$fb[2] = beta_fixed = 1 # <- Since b = [b_1 (free), 1 (fixed)], we need to change fb entry
  mod$lag_length = 1 # The p of the error AR(p)-process
  mod$stationary <- TRUE # indicates whether the stationary distribution shall be used to start the state process
  
  # Setup for theta_0 true parameter value
  tau = c(0,rep.int(log(2),5))
  theta_0 = c(beta_free, chol(omega_var), tau, rho) # (beta (-coef1), Omega, Sigma (-all), tau, rho)
  
  # Advance
  time = c(1:Tfull) # choice occasion ID
  Tp = length(time)*quest # number of questions per person 
  
  # Setup dataframe and system for DGP:
  system <- build_system_from_model_AR_R(theta_0, mod, time)
  
  if(z == 1){
    df_old = draw_simulation_data(Tfull, N, quest, rho, system)
    print(head(df_old))
  } else {
    df_new = draw_new_simulation_data(Tfull, N, quest, rho, system)
    print(head(df_new))
  }
}



# -------------------------------------------------------------------------
# Functions

draw_simulation_data <- function(Tfull, N, quest, rho, system){
  
  time = c(1:Tfull)
  Tp = length(time)*quest # number of questions per person 
  
  # Setup dataframe by first individual
  latent_class = sample(x = c(1,2,3), size = 1, prob = c(rep.int(1/3,3)))
  latent_means = c(2,4,7)
  X = matrix(0,Tp,2)
  X[,1] = rnorm(n = Tp, mean = latent_means[latent_class], sd = 1) # b_1 is free
  X[,2] = rnorm(n = Tp, mean = 0, sd = 1) # b_2 is fixed 
  tau = system$tauk
  b = system$beta
  #beta = rnorm(n = Tp, mean = 0, sd = 0.5) # mixed part
  #epsilon = rnorm(n = 1, sd = system$Sigma)
  LG = t(chol(system$GammaT))
  LO = t(chol(system$Omega))
  lRE = dim(LO)[1]
  # omega_sd = system$Omega
  # gam = LO %*% rnorm(lRE, sd = omega_sd) # mixed part here
  gam = LO %*% rnorm(lRE) # mixed part here
  
  # Obtain utilities for each choice occasion for individual 1
  U = matrix(0, nrow = Tp)
  U = X %*% system$beta + X[,c(1:lRE)] %*% gam + LG %*% rnorm(dim(LG)[2])
  
  # Choose appropriate alternative by utility value
  y = U*0
  for (tj in 1:Tp){
    y[tj] = sum(U[tj]>tau)+1
  }
  
  # Now, initialise the df
  time = c(1:Tfull) # choice occasion ID
  quest_df = rep(c(1:quest),length(time))
  iota_q = rep(1,quest)
  time_df = matrix(0,Tp,1)
  for (j in 1:length(time)){
    time_df[(j-1)*quest+c(1:quest),1] <- iota_q * time[j]
  }
  df = data.frame(ID = 1, X1 = X[,1], X2 = X[,2], y = y, time = time_df, quest = quest_df)
  
  # cycle over deciders 
  for (j in 2:N){
    latent_class = sample(x = c(1,2,3), size = 1, prob = c(rep.int(1/3,3)))
    X = matrix(0,Tp,2)
    (X[,1] = rnorm(n = Tp, mean = latent_means[latent_class], sd = 1)) # b_1 is free
    X[,2] = rnorm(n = Tp, mean = 0, sd = 1) # b_2 is fixed 
    LG = t(chol(system$GammaT))
    LO = t(chol(system$Omega))
    lRE = dim(LO)[1]
    gam = LO %*% rnorm(lRE)
    
    # Obtain utilities for each choice occasion for individual 1
    U = matrix(0, nrow = Tp)
    U = X %*% system$beta + X[,c(1:lRE)] %*% gam + LG %*% rnorm(dim(LG)[2])
    
    # Choose appropriate alternative by utility value
    y = U*0
    for (tj in 1:Tp){
      y[tj] = sum(U[tj]>system$tauk)+1
    }
    df = rbind(df, data.frame(ID = j, X1 = X[,1], X2 = X[,2], y = y, time = time_df, quest = quest_df))
  } 
  
  return(df)
}

draw_new_simulation_data <- function(Tfull, N, quest, rho, system) {
  
  time <- rep(1:Tfull, each = quest)
  Tp <- length(time)
  latent_means <- c(2, 4, 7)
  
  # ??? Preallocate list to avoid slow `rbind()`
  data_list <- vector("list", N)
  for (j in 1:N) {
    latent_class <- sample(1:3, 1, prob = c(1/3, 1/3, 1/3))
    
    X = matrix(c(rnorm(n = Tp, mean = latent_means[latent_class], sd = 1),
                 rnorm(n = Tp, mean = 0, sd = 1)),
               nrow = Tp, ncol = 2)
    
    LG <- t(chol(system$GammaT))
    LO <- t(chol(system$Omega))
    lRE = nrow(LO)
    
    gam = LO %*% rnorm(lRE)
    
    U = matrix(0, nrow = Tp)
    U = X %*% system$beta + X[, c(1:lRE)] %*% gam + LG %*% rnorm(dim(LG)[2])
    
    y = U*0
    y <- rowSums(outer(U, system$tauk, ">")) + 1
    
    data_list[[j]] <- data.table(ID = j, X1 = X[,1], X2 = X[,2], y = y, time = time, quest = rep(1:quest, Tfull))
  }
  return(rbindlist(data_list))
}

get_simulation_estimates <- function(Tfull, # number of choice occasions (time)
                                     N, # individuals
                                     quest = 1, beta_free, rho, omega_var, return_val = "", timer = F, DontSkipFit = T){
  mod <- mod_AR_cl$new(alt  = 7,
                       Hb   = diag(2)[,-2,drop=FALSE],
                       fb   = matrix(0,2,1),
                       HO   = matrix(1,1,1),
                       fO   = matrix(0,1,1),
                       HL   = matrix(0,1,0),
                       fL   = matrix(1,1,1),
                       ordered = TRUE
  )
  
  mod$fb[2] = beta_fixed = 1 # <- Since b = [b_1 (free), 1 (fixed)], we need to change fb entry
  mod$lag_length = 1 # The p of the error AR(p)-process
  mod$stationary <- TRUE # indicates whether the stationary distribution shall be used to start the state process
  
  # Setup for theta_0 true parameter value
  tau = c(0,rep.int(log(2),5))
  theta_0 = c(beta_free, chol(omega_var), tau, rho) # (beta (-coef1), Omega, Sigma (-all), tau, rho)
  
  # Advance
  time = c(1:Tfull) # choice occasion ID
  Tp = length(time)*quest # number of questions per person 
  
  # Setup dataframe and system for DGP:
  system <- build_system_from_model_AR_R(theta_0, mod, time)
  df = draw_simulation_data(Tfull = Tfull, N = N, quest = quest, rho = rho, system = system)
  
  # Convert dataframe to data_raw_StSp_cl for estimation:
  data_raw <- data_raw_StSp_cl$new(df = df, alt_names = c("1","2","3","4","5","6","7"), id = "ID", choice = "y", ordered = TRUE, varying = "", dec_char = c("X1","X2"))
  data_raw$set_time_col("time")
  data_raw$set_quest_col("quest")
  data_raw$alt_names <- alt_names <- sprintf('%d',1:7)
  
  # set up Rprobit_obj: only type 2 regressors as usual in ordered probit situations. 
  form = y ~ 0 | X1 + X2 | 0
  
  control_simul <- list(Tp = Tp)
  control <- list(control_simulation = control_simul)
  re <- c("X1")
  
  # generate Rprobit_obj 
  Rprobit_obj <- setup_Rprobit(form = form, data_raw = data_raw,
                               mod = mod,
                               re = re,
                               seed = 17,
                               theta_0 = theta_0,
                               control = control_simul
  )
  
  if(timer == T){
    tic("Fitting of Rprobit Object") # Only of interest in case of close to true variance estimation, which is later used as baseline for CI-calcs
  }
  # fit the system with cml_pair_type = 2 (only adjacent pairs + first and last decision)
  if(DontSkipFit){
    result <- tryCatch(
      { 
        Rprob_mod <- suppressMessages(fit_Rprobit(Rprobit_obj, init_method = "theta", cml_pair_type = 0))
        
        if(nzchar(return_val)){
          switch(return_val,
                 "summary" = summary(Rprob_mod),
                 "theta" = replace(Rprob_mod$theta, 2, Rprob_mod$theta[2]^2), # replace second entry (omega) by delta method
                 "var" = replace(diag(Rprob_mod$vv)^2, 2, 4*Rprob_mod$theta[2]^2*diag(Rprob_mod$vv)[2]), # replace second entry (omega) by delta method
                 "sd" = replace(sqrt(diag(Rprob_mod$vv)), 2, sqrt(4*Rprob_mod$theta[2]^2*diag(Rprob_mod$vv)[2])))
        } else {
          # Omega is Var-Cov, d.h. für Vergleich entweder sd quadrieren oder var wurzeln 
          b = c(b_coef = Rprob_mod$theta[1], b_sd = sqrt(diag(Rprob_mod$vv))[1])
          omega = c(omega_coef = Rprob_mod$theta[2]^2, omega_sd = sqrt(4*Rprob_mod$theta[2]^2*diag(Rprob_mod$vv)[2])) # coef_omega = theta_omega^2 (omega matrix entry), 
          rho = c(rho_coef = Rprob_mod$theta[9], rho_sd = sqrt(diag(Rprob_mod$vv))[9])                                # sd = delta method (sd of omega mat entry)
          conv = c(nlm_code = Rprob_mod$info$nlm_info$code)
          tau_raw = c(tau_coef = Rprob_mod$theta[3:8], tau_sd = sqrt(diag(Rprob_mod$vv)[3:8]))
          #Rprob_conv = c(H_conv = T, J_conv = T)
          list(b, omega, rho, conv, tau_raw)
        }
      },
      error = function(e) {
        if(nzchar(return_val)){
          switch(return_val,
                 "summary" = print("Model fitting failed: nlm crashed"),
                 "theta" = c(rep.int(NA, length(theta_0))),
                 "var" = c(rep.int(NA, length(theta_0))),
                 "sd" = c(rep.int(NA, length(theta_0))))
          
          message("Model fitting failed: ", conditionMessage(e))
        } else {
          # If an error occurs, return NA values with the same structure
          b <- c(b_coef = NA, b_sd = NA)
          omega <- c(omega_coef = NA, omega_sd = NA)
          rho <- c(rho_coef = NA, rho_sd = NA)
          conv <- c(nlm_code = 99)
          tau_raw = c(tau_coef = rep.int(NA, 6), tau_sd = rep.int(NA, 6))
          
          message("Model fitting failed: ", conditionMessage(e))
          
          # Return the same structure with NA values
          list(b = b, omega = omega, rho = rho, conv = conv, tau_raw = tau_raw)
        }
      }
    )
  } else { # if skipped, return NA same as if error occured
    b <- c(b_coef = NA, b_sd = NA)
    omega <- c(omega_coef = NA, omega_sd = NA)
    rho <- c(rho_coef = NA, rho_sd = NA)
    conv <- c(nlm_code = 99)
    tau_raw = c(tau_coef = rep.int(NA, 6), tau_sd = rep.int(NA, 6))
    result = list(b = b, omega = omega, rho = rho, conv = conv, tau_raw = tau_raw)
  }
  if(timer == T){
    toc()
  }
  # Return results 
  return(result)
}




# Some graphic related things ---------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyverse)
library(viridis)
library(ggnewscale)
library(grid)
# par(mfrow=c(1, 2))
original = read.csv(paste0(getwd(), "/GitHub/MasterThesis/DGP1_3000reps_Results_allB_v3.csv")) %>% filter(N != 500)
original2 = read.csv(paste0(getwd(), "/GitHub/MasterThesis/DGP1_3000reps_Results_allB_v5.csv"))
original = rbind(original, original2)

# Beta Bias Pics ----------------------------------------------------------
## Varying mean beta coefficient increases bias ----
data = original %>% filter(rho == -0.9, Tfull == 2, n_reps == 3000, omega_var == 1)
alpha = 0.05
z_value = qnorm(1-alpha/2)
sqrtN_to_divide_by = sqrt(as.numeric(data$success_beta + data$failure_beta))
data = data %>%
  mutate(beta_lower = beta_bias - z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         beta_upper = beta_bias + z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         omega_lower = omega_bias - z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         omega_upper = omega_bias + z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         rho_lower = rho_bias - z_value * sd_of_rho_estimates/sqrtN_to_divide_by,
         rho_upper = rho_bias + z_value * sd_of_rho_estimates/sqrtN_to_divide_by)

ggplot(data, aes(x = N, y = beta_bias, group = beta_free, color = as.factor(beta_free))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = beta_lower, ymax = beta_upper, fill = as.factor(beta_free)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Bias of Mean β1 Estimate vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Bias of Mean β1",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = c(0, 1)) + 
  scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 100)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")

## Varying omega_var increases mean beta bias ----
data = original %>% filter(rho == -0.9, Tfull == 2, n_reps == 3000, beta_free == 2)
alpha = 0.05
z_value = qnorm(1-alpha/2)
sqrtN_to_divide_by = sqrt(as.numeric(data$success_beta + data$failure_beta))
data = data %>%
  mutate(beta_lower = beta_bias - z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         beta_upper = beta_bias + z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         omega_lower = omega_bias - z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         omega_upper = omega_bias + z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         rho_lower = rho_bias - z_value * sd_of_rho_estimates/sqrtN_to_divide_by,
         rho_upper = rho_bias + z_value * sd_of_rho_estimates/sqrtN_to_divide_by)

ggplot(data, aes(x = N, y = beta_bias, group = omega_var, color = as.factor(omega_var))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = beta_lower, ymax = beta_upper, fill = as.factor(omega_var)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nβ = 2,\nρ = -0.9\n\n# of Reps:\n~3000\n\nVar ω²:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nβ = 2,\nρ = -0.9\n\n# of Reps:\n~3000\n\nVar ω²:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Bias of Mean β1 Estimate vs. Sample Size (N)",
    x = "Sample Size (N)",
    y = "Bias of Mean β1",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  coord_cartesian(ylim = c(0, 1)) + 
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")


## Compare avg. model-based SE with empirical SE ----
data = original %>% filter(rho == -0.9, Tfull == 2, n_reps == 3000, beta_free == 1) %>% 
  filter(beta_free %in% c(0,1,2,3,4), omega_var %in% c(0.25^2, 1, 4^2)) 

ggplot(data, aes(x = N, color = as.factor(omega_var))) + 
  # Model-based SE (solid line)
  geom_line(aes(y = avg_of_beta_se, linetype = "Avg. Model-\nbased")) + 
  geom_point(aes(y = avg_of_beta_se), size = 1) +
  
  # Empirical SE (dashed line)
  geom_line(aes(y = sd_of_beta_estimates, linetype = "Empirical")) +
  geom_point(aes(y = sd_of_beta_estimates), shape = 17, size = 1) +
  
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Mean β1") +
  
  # Add linetype legend
  scale_linetype_manual(
    name = "SE Type",
    values = c("Avg. Model-\nbased" = "solid", "Empirical" = "dashed")
  ) +
  
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    )
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  
  labs(
    title = expression("SE Comparison for Mean β1 Estimates vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Standard Error"
  ) +
  
  coord_cartesian(ylim = c(0, 4)) + 
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")



## Rho vs. Beta Empirical SE -----
data = original %>% filter(Tfull == 2, N == 300, n_reps == 3000, omega_var == 1)

ggplot(data, aes(x = rho, y = sd_of_beta_estimates, group = beta_free, color = as.factor(beta_free))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nN = 300,\nω² = 1\n\n# of Reps:\n~3000\n\nMean β1:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nN = 300,\nω² = 1\n\n# of Reps:\n~3000\n\nMean β1:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Empirical SE of beta estimate vs. AR(1) coefficient ρ",
    x = "Rho (ρ)",
    y = "Empirical standard error",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) + 
  coord_cartesian(ylim = c(0, 1))


## Rho vs. Beta Bias -----
data = original %>% filter(Tfull == 2, N == 100, n_reps == 3000, omega_var == 0.25^2) %>%
  filter(beta_free %in% c(0, 1, 2, 3, 4, 0.4, 0.8))
ggplot(data, aes(x = rho, y = beta_bias, group = beta_free, color = as.factor(beta_free))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 3,\nN = 300,\nω² = 0.25^2\n\n# of Reps:\n~3000\n\nMean β1:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 3,\nN = 300,\nω² = 0.25^2\n\n# of Reps:\n~3000\n\nMean β1:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Rho vs. Beta Bias",
    x = "Rho",
    y = "Beta Bias",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0, 0.2))


## TRASH: Varying rho and the effect on bias over N -------------------------------
data = original %>% filter(Tfull == 2, n_reps == 3000, beta_free == 2, omega_var == 16)
alpha = 0.05
z_value = qnorm(1-alpha/2)
sqrtN_to_divide_by = sqrt(as.numeric(data$success_beta + data$failure_beta))
data = data %>%
  mutate(beta_lower = beta_bias - z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         beta_upper = beta_bias + z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         omega_lower = omega_bias - z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         omega_upper = omega_bias + z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         rho_lower = rho_bias - z_value * sd_of_rho_estimates/sqrtN_to_divide_by,
         rho_upper = rho_bias + z_value * sd_of_rho_estimates/sqrtN_to_divide_by)

ggplot(data, aes(x = N, y = beta_bias, group = rho, color = rho)) + 
  geom_line(size = 1) + 
  geom_point(size = 1) +
  scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 0) +
  
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Bias of Mean β1 Estimate vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Bias of Mean β1",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = c(0, 0.5)) + 
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")

## Bias plots as T increases  -------------------------------
data = original %>% filter(rho == 0.3, n_reps == 3000, omega_var == 1, beta_free == 2)
alpha = 0.05
z_value = qnorm(1-alpha/2)
sqrtN_to_divide_by = sqrt(as.numeric(data$success_beta + data$failure_beta))
data = data %>%
  mutate(beta_lower = beta_bias - z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         beta_upper = beta_bias + z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         omega_lower = omega_bias - z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         omega_upper = omega_bias + z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         rho_lower = rho_bias - z_value * sd_of_rho_estimates/sqrtN_to_divide_by,
         rho_upper = rho_bias + z_value * sd_of_rho_estimates/sqrtN_to_divide_by)

ggplot(data, aes(x = N, y = beta_bias, group = Tfull, color = as.factor(Tfull))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = beta_lower, ymax = beta_upper, fill = as.factor(Tfull)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nβ = 2,\nω² = 1,\nρ = 0.3\n\n# of Reps:\n~3000\n\nTime T:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nβ = 2,\nω² = 1,\nρ = 0.3\n\n# of Reps:\n~3000\n\nTime T:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Bias of Mean β1 Estimate vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Bias of Mean β1",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = c(0, 0.3)) + 
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")

## Rho vs. Coverage Probs with beta increasing  -----
data = original %>% filter(Tfull == 2, N == 200, n_reps == 3000, omega_var == 1)
ggplot(data, aes(x = rho, y = p_beta, group = beta_free, color = as.factor(beta_free))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.4, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.4, end = 0.9) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nN = 200,\nω² = 1\n\n# of Reps:\n~3000\n\nMean β1:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nN = 200,\nω² = 1\n\n# of Reps:\n~3000\n\nMean β1:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Coverage Probability for Mean β vs. AR(1) coefficient ρ",
    x = "Rho",
    y = "Coverage Probability",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0.85, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")

data$beta_CI_width <- data$bino_ci_exact_upper_beta - data$bino_ci_exact_lower_beta
data <- data %>% dplyr::relocate(beta_CI_width, .before = 1)
data <- data %>% dplyr::relocate(rho, .before = 1)
data <- data %>% dplyr::relocate(beta_free, .before = 1)

data <- data %>%
  mutate(beta_CI_width_rank = rank(-beta_CI_width)) 

## Varying Omega in coverage plots -----------------------------------------
data = original %>% filter(rho == -0.9, Tfull == 2, n_reps == 3000, beta_free == 0)

ggplot(data, aes(x = N, y = p_beta, group = omega_var, color = as.factor(omega_var))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_beta, ymax = bino_ci_exact_upper_beta, fill = as.factor(omega_var)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nβ = 4,\nρ = -0.9\n\n# of Reps:\n~3000\n\nVar ω²:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nβ = 4,\nρ = -0.9\n\n# of Reps:\n~3000\n\nVar ω²:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Coverage Probability for Mean β vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = c(0.7, 1)) + 
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black") 

## Other param vs. coverage Probs ------
data = original %>% filter(N == 200, n_reps == 3000, beta_free == 2, omega_var == 1)
ggplot(data, aes(x = rho, y = p_beta, group = Tfull, color = as.factor(Tfull))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.4, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.4, end = 0.9) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nβ = 2,\nN = 200,\nω² = 1\n\n# of Reps:\n~3000\n\nTime T:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nβ = 2,\nN = 200,\nω² = 1\n\n# of Reps:\n~3000\n\nTime T:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Coverage Probability for Mean β vs. AR(1) coefficient ρ",
    x = "Mean of Beta Coefficient",
    y = "Coverage Probability",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0.85, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")
## Varying T in coverage plots -----------------------------------------
data = original %>% filter(rho == -0.9, beta_free == 3, n_reps == 3000, omega_var == 0.25^2)

ggplot(data, aes(x = N, y = p_beta, group = Tfull, color = as.factor(Tfull))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_beta,
                  ymax = bino_ci_exact_upper_beta,
                  fill = as.factor(Tfull)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nβ = 3,\nρ = -0.9\nω² = 0.25^2\n\n# of Reps:\n~3000\n\nTime T:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nβ = 3,\nρ = -0.9\nω² = 0.25^2\n\n# of Reps:\n~3000\n\nTime T:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Coverage Probability for Mean β vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = c(0.825, 1)) + 
  scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 100)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")


## Compare coverage rates between SE types while varying Beta -------------------------------
data = original %>% filter(rho == 0.6, Tfull == 2, n_reps == 3000, omega_var == 1, beta_free %in% c(1,3))

ggplot(data, aes(x = N, group = beta_free)) + 
  # First set: original p_beta
  geom_line(aes(y = p_beta, color = as.factor(beta_free))) + 
  geom_point(aes(y = p_beta, color = as.factor(beta_free)), size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_beta, ymax = bino_ci_exact_upper_beta, fill = as.factor(beta_free)), 
              alpha = 0.3, color = NA) +
  
  # First color & fill scale: blue-green
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Mean β and\nModel-based") +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Mean β and\nModel-based") +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_line(aes(y = p_beta_rescaled, color = as.factor(beta_free))) + 
  geom_point(aes(y = p_beta_rescaled, color = as.factor(beta_free)), shape = 1, size = 1.2) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_beta_rescaled, ymax = bino_ci_exact_upper_beta_rescaled, fill = as.factor(beta_free)), 
              alpha = 0.3, color = NA) +
  
  scale_color_manual(
    values = c("red", "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  scale_fill_manual(
    values = c("red", "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = 0.6\n\n# of Reps:\n~3000\n\nMean β and\nRescaled:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = 0.6\n\n# of Reps:\n~3000\n\nMean β and\nRescaled:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = expression("CI Comparison for Mean β1 Estimates vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability"
  ) +
  coord_cartesian(ylim = c(0.75, 1)) + 
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")
## Compare coverage rates between SE types while varying Omega -------------------------------
data = original %>% filter(rho == 0.6, Tfull == 2, n_reps == 3000, beta_free == 2)

ggplot(data, aes(x = N, group = omega_var)) + 
  # First set: original p_beta
  geom_line(aes(y = p_beta, color = as.factor(omega_var))) + 
  geom_point(aes(y = p_beta, color = as.factor(omega_var)), size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_beta, ymax = bino_ci_exact_upper_beta, fill = as.factor(omega_var)), 
              alpha = 0.3, color = NA) +
  
  # First color & fill scale: blue-green
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Var ω² and\nModel-based") +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Var ω² and\nModel-based") +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_line(aes(y = p_beta_rescaled, color = as.factor(omega_var))) + 
  geom_point(aes(y = p_beta_rescaled, color = as.factor(omega_var)), shape = 1, size = 1.2) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_beta_rescaled, ymax = bino_ci_exact_upper_beta_rescaled, fill = as.factor(omega_var)), 
              alpha = 0.3, color = NA) +
  
  scale_color_manual(
    values = c("purple", "red",  "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  scale_fill_manual(
    values = c("purple", "red",  "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nβ = 2,\nρ = 0.6\n\n# of Reps:\n~3000\n\nVar ω² and\nRescaled:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nβ = 2,\nρ = 0.6\n\n# of Reps:\n~3000\n\nVar ω² and\nRescaled:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = expression("CI Comparison for Mean β1 Estimates vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability"
  ) +
  coord_cartesian(ylim = c(0.8, 1)) + 
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")

# Omega Plots -------------------------------------------------------------
## Omega variance bias plot for varying beta ----
data = original %>% filter(rho == -0.9, Tfull == 2, n_reps == 3000, omega_var == 1) %>%
  filter(beta_free %in% c(0, 2, 4))
alpha = 0.05
z_value = qnorm(1-alpha/2)
sqrtN_to_divide_by = sqrt(as.numeric(data$success_beta + data$failure_beta))
data = data %>%
  mutate(beta_lower = beta_bias - z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         beta_upper = beta_bias + z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         omega_lower = omega_bias - z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         omega_upper = omega_bias + z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         rho_lower = rho_bias - z_value * sd_of_rho_estimates/sqrtN_to_divide_by,
         rho_upper = rho_bias + z_value * sd_of_rho_estimates/sqrtN_to_divide_by)

ggplot(data, aes(x = N, y = omega_bias, group = beta_free, color = as.factor(beta_free))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = omega_lower, ymax = omega_upper, fill = as.factor(beta_free)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Bias of Mean β1 Estimate vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Bias of Mean β1",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = c(0, 1)) + 
  scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 100)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")

## Omega variance bias plot for varying omega ----
data = original %>% filter(rho == -0.9, Tfull == 2, n_reps == 3000, beta_free == 2) %>%
  filter(beta_free %in% c(0, 2, 4), omega_var %in% c(0.25^2, 1, 4^2))
alpha = 0.05
z_value = qnorm(1-alpha/2)
sqrtN_to_divide_by = sqrt(as.numeric(data$success_beta + data$failure_beta))
data = data %>%
  mutate(beta_lower = beta_bias - z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         beta_upper = beta_bias + z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         omega_lower = omega_bias - z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         omega_upper = omega_bias + z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         rho_lower = rho_bias - z_value * sd_of_rho_estimates/sqrtN_to_divide_by,
         rho_upper = rho_bias + z_value * sd_of_rho_estimates/sqrtN_to_divide_by)

ggplot(data, aes(x = N, y = omega_bias, group = omega_var, color = as.factor(omega_var))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = omega_lower, ymax = omega_upper, fill = as.factor(omega_var)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Bias of Mean β1 Estimate vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Bias of Mean β1",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = c(0, 5)) + 
  scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 100)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")

## Rho vs. Omega Bias (param in color) -----
data = original %>% filter(Tfull == 2, N == 500, n_reps == 3000, beta_free == 1)
ggplot(data, aes(x = rho, y = omega_bias, group = omega_var, color = as.factor(omega_var))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  guides(
    color = guide_legend(
      title = paste0("Fixed\nParams:\nT = ",unique(data$Tfull),
                     "\nN = ",unique(data$N),"\nω² = ",
                     unique(data$omega_var),
                     "\n\n# of Reps:\n~3000\n\nMean β1:")
    ),
    fill = guide_legend(
      title = paste0("Fixed\nParams:\nT = ",
                     unique(data$Tfull),"\nN = ",
                     unique(data$N),"\nω² = ",
                     unique(data$omega_var),
                     "\n\n# of Reps:\n~3000\n\nMean β1:")
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Rho vs. Omega Bias",
    x = "Rho",
    y = "Omega Bias",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0, 0.04))

# ---

what_to_plot = list(omega_var = 0.25^2,
                    N = 90,
                    Tfull = 2)

data_dense  <- original %>%
  filter(Tfull == what_to_plot$T, N == what_to_plot$N, n_reps == 3000, omega_var == what_to_plot$omega_var) %>%
  filter(beta_free <= 1) %>%
  filter(!(beta_free %in% c(0.2, 0.6)))

data_coarse <- original %>%
  filter(Tfull == what_to_plot$T, N == what_to_plot$N, n_reps == 3000, omega_var == what_to_plot$omega_var) %>% filter(beta_free %in% c(1,2,3,4))

ggplot() +
  geom_line(data = data_coarse,
            aes(x = rho, y = omega_bias, group = beta_free, color = as.factor(beta_free)),
            size = 0.9) +
  geom_point(data = data_coarse,
             aes(x = rho, y = omega_bias, color = as.factor(beta_free)),
             size = 1.8) +
  scale_color_manual(
    name = NULL,
    values = c(`2` = "#FDE725", `3` = "#FDB366", `4` = "#D73027"),
    guide = guide_legend(order = 2)   # put this legend BELOW
  ) +
  ggnewscale::new_scale_color() +
  geom_line(data = data_dense,
            aes(x = rho, y = omega_bias, group = beta_free, color = as.factor(beta_free)),
            size = 0.75) +
  geom_point(data = data_dense,
             aes(x = rho, y = omega_bias, color = as.factor(beta_free)),
             size = 1.5) +
  scale_color_viridis_d(
    option = "D", begin = 0.4, end = 0.9,
    name = paste0(
      "Fixed\nParams:\nT = ", unique(data_dense$Tfull),
      "\nN = ", unique(data_dense$N),
      "\nω² = ", unique(data_dense$omega_var),
      "\n\n# of Reps:\n~3000\n\nMean β1:"
    ),
    guide = guide_legend(order = 1) 
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0, "cm"),
    plot.title = element_text(hjust = 0.1)
  ) +
  labs(
    title = "Rho vs. beta Bias",
    x = "Rho",
    y = "beta Bias"
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0, 0.08)) # 0.1/0.04 for var=1 & T=2/3, 4 for Var=16 & both T

## Rho vs. Omega Empirical SE -----
data = original %>% filter(Tfull == 2, N == 300, n_reps == 3000, omega_var == 1)
ggplot(data, aes(x = rho, y = sd_of_omega_estimates, group = beta_free, color = as.factor(beta_free))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  guides(
    color = guide_legend(
      title = paste0("Fixed\nParams:\nT = ",unique(data$Tfull),"\nN = ",unique(data$N),"\nω² = ",unique(data$omega_var),"\n\n# of Reps:\n~3000\n\nMean β1:")
    ),
    fill = guide_legend(
      title = paste0("Fixed\nParams:\nT = ",unique(data$Tfull),"\nN = ",unique(data$N),"\nω² = ",unique(data$omega_var),"\n\n# of Reps:\n~3000\n\nMean β1:")
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Rho vs. Omega Bias",
    x = "Rho",
    y = "Omega Bias",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0, 0.2))

# ---

what_to_plot = list(omega_var = 2^2,
                 N = 500,
                 Tfull = 2)

data_dense  <- original %>%
  filter(Tfull == what_to_plot$T, N == what_to_plot$N, n_reps == 3000, omega_var == what_to_plot$omega_var) %>%
  filter(beta_free <= 1) %>%
  filter(!(beta_free %in% c(0.2, 0.6)))

data_coarse <- original %>%
  filter(Tfull == what_to_plot$T, N == what_to_plot$N, n_reps == 3000, omega_var == what_to_plot$omega_var) %>% filter(beta_free %in% c(1,2,3,4))

ggplot() +
  geom_line(data = data_coarse,
            aes(x = rho, y = sd_of_omega_estimates, group = beta_free, color = as.factor(beta_free)),
            size = 0.9) +
  geom_point(data = data_coarse,
             aes(x = rho, y = sd_of_omega_estimates, color = as.factor(beta_free)),
             size = 1.8) +
  scale_color_manual(
    name = NULL,
    values = c(`2` = "#FDE725", `3` = "#FDB366", `4` = "#D73027"),
    guide = guide_legend(order = 2)   # put this legend BELOW
  ) +
  ggnewscale::new_scale_color() +
  geom_line(data = data_dense,
            aes(x = rho, y = sd_of_omega_estimates, group = beta_free, color = as.factor(beta_free)),
            size = 0.75) +
  geom_point(data = data_dense,
             aes(x = rho, y = sd_of_omega_estimates, color = as.factor(beta_free)),
             size = 1.5) +
  scale_color_viridis_d(
    option = "D", begin = 0.2, end = 0.7,
    name = paste0(
      "Fixed\nParams:\nT = ", unique(data_dense$Tfull),
      "\nN = ", unique(data_dense$N),
      "\nω² = ", unique(data_dense$omega_var),
      "\n\n# of Reps:\n~3000\n\nMean β1:"
    ),
    guide = guide_legend(order = 1)   # put this legend ABOVE
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0, "cm"),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Rho vs. Omega Empirical SE",
    x = "Rho",
    y = "Empirical Standard Error"
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0, 3)) # 0.1/0.04 for var=1 & T=2/3, 4 for Var=16 & both T


## Omega Model-based SE vs. Empirical SE (Vary beta) -----------------------------------
data = original %>% filter(rho == 0.9, Tfull == 2, n_reps == 3000, omega_var == 0.25^2) %>% filter(beta_free %in% c(0,1,2,3,4))

ggplot(data, aes(x = N, color = as.factor(beta_free))) + 
  # Model-based SE (solid line)
  geom_line(aes(y = avg_of_omega_se, linetype = "Avg. Model-\nbased")) + 
  geom_point(aes(y = avg_of_omega_se), size = 1) +
  
  # Empirical SE (dashed line)
  geom_line(aes(y = sd_of_omega_estimates, linetype = "Empirical")) +
  geom_point(aes(y = sd_of_omega_estimates), shape = 17, size = 1) +
  
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Mean β1") +
  
  # Add linetype legend
  scale_linetype_manual(
    name = "SE Type",
    values = c("Avg. Model-\nbased" = "solid", "Empirical" = "dashed")
  ) +
  
  guides(
    color = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nω² = ", unique(data$omega_var),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nMean β1:"
      )
    ),
    fill = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nω² = ", unique(data$omega_var),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nMean β1:"
      )
    )
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  
  labs(
    title = expression("SE Comparison for Omega Variance Estimates vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Standard Error"
  ) +
  
  coord_cartesian(ylim = c(0, 3)) + 
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")

## Omega Model-based SE vs. Empirical SE (Vary Omega) -----------------------------------
data = original %>% filter(rho == 0.9, Tfull == 2, n_reps == 3000, beta_free == 2) %>% filter(omega_var %in% c(0.25^2, 1, 4^2))

ggplot(data, aes(x = N, color = as.factor(omega_var))) + 
  # Model-based SE (solid line)
  geom_line(aes(y = avg_of_omega_se, linetype = "Avg. Model-\nbased")) + 
  geom_point(aes(y = avg_of_omega_se), size = 1) +
  
  # Empirical SE (dashed line)
  geom_line(aes(y = sd_of_omega_estimates, linetype = "Empirical")) +
  geom_point(aes(y = sd_of_omega_estimates), shape = 17, size = 1) +
  
  scale_color_viridis_d(option = "D", begin = 0.4, end = 0.9, name = "Mean β1") +
  
  # Add linetype legend
  scale_linetype_manual(
    name = "SE Type",
    values = c("Avg. Model-\nbased" = "solid", "Empirical" = "dashed")
  ) +
  
  guides(
    color = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nβ = ", unique(data$beta_free),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nVar ω²:"
      )
    ),
    fill = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nβ = ", unique(data$beta_free),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nVar ω²:"
      )
    )
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  
  labs(
    title = expression("SE Comparison for Omega Variance Estimates vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Standard Error"
  ) +
  
  coord_cartesian(ylim = c(0, 7.5)) + 
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")


## Omega Bias figure ----
data = original %>% filter(rho == 0, Tfull == 2, n_reps == 3000, omega_var == 1) %>% filter(beta_free %in% c(0,2,4))
alpha = 0.05
z_value = qnorm(1-alpha/2)
sqrtN_to_divide_by = sqrt(as.numeric(data$success_beta + data$failure_beta))
data = data %>%
  mutate(beta_lower = beta_bias - z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         beta_upper = beta_bias + z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         omega_lower = omega_bias - z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         omega_upper = omega_bias + z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         rho_lower = rho_bias - z_value * sd_of_rho_estimates/sqrtN_to_divide_by,
         rho_upper = rho_bias + z_value * sd_of_rho_estimates/sqrtN_to_divide_by)

ggplot(data, aes(x = N, y = omega_bias, group = beta_free, color = as.factor(beta_free))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = omega_lower, ymax = omega_upper, fill = as.factor(beta_free)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Bias of Mean β1 Estimate vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Bias of Mean β1",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = c(0, 0.5)) + 
  scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 100)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")



## Rho vs. Omega Bias (+ Sample size in Color) -----------------------------------
data = original %>% filter(Tfull == 2, n_reps == 3000, beta_free == 2, omega_var == 1^2) %>% filter(N %in% c(50, 80, 100, 150, 200))
ggplot(data, aes(x = rho, y = omega_bias, group = N, color = as.factor(N))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  guides(
    color = guide_legend(
      title = paste0("Fixed\nParams:\nT = ",unique(data$Tfull),"\nN = ",unique(data$N),"\nω² = ",unique(data$omega_var),"\n\n# of Reps:\n~3000\n\nMean β1:")
    ),
    fill = guide_legend(
      title = paste0("Fixed\nParams:\nT = ",unique(data$Tfull),"\nN = ",unique(data$N),"\nω² = ",unique(data$omega_var),"\n\n# of Reps:\n~3000\n\nMean β1:")
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Rho vs. Omega Bias",
    x = "Rho",
    y = "Omega Bias",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0, 1))

## Rho vs. omega Coverage Probs with beta increasing  -----
data = original %>% filter(Tfull == 2, N == 500, n_reps == 3000, omega_var == 1^2) %>%
  filter(omega_var %in% c(0.25^2, 1, 4^2), beta_free %in% c(0, 1, 2, 3, 4))
ggplot(data, aes(x = rho, y = p_omega_rescaled, group = beta_free, color = as.factor(beta_free))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.4, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.4, end = 0.9) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nN = 200,\nω² = 1\n\n# of Reps:\n~3000\n\nMean β1:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nN = 200,\nω² = 1\n\n# of Reps:\n~3000\n\nMean β1:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Coverage Probability for Mean β vs. AR(1) coefficient ρ",
    x = "Rho",
    y = "Coverage Probability",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0.85, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")

data$beta_CI_width <- data$bino_ci_exact_upper_beta - data$bino_ci_exact_lower_beta
data <- data %>% dplyr::relocate(beta_CI_width, .before = 1)
data <- data %>% dplyr::relocate(rho, .before = 1)
data <- data %>% dplyr::relocate(beta_free, .before = 1)

data <- data %>%
  mutate(beta_CI_width_rank = rank(-beta_CI_width)) 

## Rho vs. omega Coverage Probs with Omega increasing  -----
data = original %>% filter(Tfull == 3, N == 100, n_reps == 3000, beta_free == 0) %>%
  filter(omega_var %in% c(0.25^2, 1, 4^2), beta_free %in% c(0, 1, 2, 3, 4))
ggplot(data, aes(x = rho, y = p_omega_rescaled, group = omega_var, color = as.factor(omega_var))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.4, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.4, end = 0.9) +
  guides(
    color = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nN = ", unique(data$N),
        "\nβ = ", unique(data$beta_free),
        "\n\n# of Reps:\n~3000\n\nVar ω²:"
      )
    ),
    fill = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nN = ", unique(data$N),
        "\nβ = ", unique(data$beta_free),
        "\n\n# of Reps:\n~3000\n\nVar ω²:"
      )
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Coverage Probability for Mean β vs. AR(1) coefficient ρ",
    x = "Rho",
    y = "Coverage Probability",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0.85, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")

## Rho vs. Omega Coverage Probs + Fix Omega and vary N  -----
data = original %>% filter(Tfull == 2, N %in% c(50, 70, 100, 150, 200, 300),
                           n_reps == 3000, beta_free == 4, omega_var == 1^2) %>%
  filter(omega_var %in% c(0.25^2, 1, 4^2), beta_free %in% c(0, 1, 2, 3, 4))
ggplot(data, aes(x = rho, y = p_omega_rescaled, group = N, color = as.factor(N))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.4, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.4, end = 0.9) +
  guides(
    color = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nω² = ", unique(data$omega_var),
        "\nβ = ", unique(data$beta_free),
        "\n\n# of Reps:\n~3000\n\nSample\nSize N:"
      )
    ),
    fill = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nω² = ", unique(data$omega_var),
        "\nβ = ", unique(data$beta_free),
        "\n\n# of Reps:\n~3000\n\nSample\nSize N:"
      )
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Comparison of Omega Coverage Probability vs. AR(1) coefficient ρ",
    x = "Rho",
    y = "Coverage Probability",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0.75, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")


## Rho vs. Omega Coverage Probs + Fix Omega and vary N and T  -----
data = original %>% filter(Tfull %in% c(3, 4) , N %in% c(50, 100, 150, 200), n_reps == 3000, beta_free == 4, omega_var == 4^2) %>%
  filter(omega_var %in% c(0.25^2, 1, 4^2), beta_free %in% c(0, 1, 2, 3, 4))
ggplot(data, aes(
  x = rho,
  y = p_omega_rescaled,
  group = interaction(N, Tfull),  # group by both
  color = as.factor(N),           # N as color
  linetype = as.factor(Tfull)     # T as line type
)) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.4, end = 0.9) +
  scale_linetype_manual(
    values = setNames(
      c("dotted", "solid"),  # first is for lower T, second for higher T
      as.character(c(min(unique(data$Tfull)), max(unique(data$Tfull))))
    )
  ) +
  scale_fill_viridis_d(option = "D", begin = 0.4, end = 0.9) +
  guides(
    color = guide_legend(order = 1,
      title = paste0(
        "Fixed\nParams:",
        "\nω² = ", unique(data$omega_var),
        "\nβ = ", unique(data$beta_free),
        "\n\n# of Reps:\n~3000\n\nSample\nSize N:"
      )
    ),
    linetype = guide_legend(
      title = "Time T:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Comparison of Coverage Success over AR(1) coefficient ρ for different N and T",
    x = "Rho",
    y = "Coverage Probability",
    color = "N",
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0.80, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")

## N vs. Omega coverage plots + Vary Omega -----------------------------------------
data = original %>% filter(rho == -0.9, Tfull == 2, n_reps == 3000, beta_free == 4) %>%
  filter(omega_var %in% c(0.25^2, 1, 4^2))

ggplot(data, aes(x = N, y = p_omega_rescaled, group = omega_var, color = as.factor(omega_var))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_omega_rescaled, ymax = bino_ci_exact_upper_omega_rescaled, fill = as.factor(omega_var)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nβ = ", unique(data$beta_free),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nVar ω²:"
      )
    ),
    fill = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nβ = ", unique(data$beta_free),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nVar ω²:"
      )
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Coverage Probability for Mean β vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = c(0.5, 1)) + 
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black") 


## N vs. Omega Coverage Probability + Varying beta ------------------------------------------
data = original %>% filter(rho == -0.9, Tfull == 2, n_reps == 3000, omega_var == 1^2) %>%
  filter(omega_var %in% c(0.25^2, 1, 4^2), beta_free %in% c(0, 2, 4))

ggplot(data, aes(x = N, y = p_omega_rescaled, group = beta_free, color = as.factor(beta_free))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_omega_rescaled, ymax = bino_ci_exact_upper_omega_rescaled, fill = as.factor(beta_free)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nω² = ", unique(data$omega_var),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nMean β:"
      )
    ),
    fill = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nω² = ", unique(data$omega_var),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nMean β:"
      )
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Coverage Probability for Mean β vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = c(0.75, 1)) + 
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black") 


## Other param vs. fixed N coverage rate ------
data = original %>% filter(N == 200, n_reps == 3000, beta_free == 2, omega_var == 1)
ggplot(data, aes(x = rho, y = p_beta, group = Tfull, color = as.factor(Tfull))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.4, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.4, end = 0.9) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nβ = 2,\nN = 200,\nω² = 1\n\n# of Reps:\n~3000\n\nTime T:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nβ = 2,\nN = 200,\nω² = 1\n\n# of Reps:\n~3000\n\nTime T:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Coverage Probability for Mean β vs. AR(1) coefficient ρ",
    x = "Mean of Beta Coefficient",
    y = "Coverage Probability",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0.85, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")
## Varying T in coverage plots -----------------------------------------
data = original %>% filter(rho == -0.9, beta_free == 3, n_reps == 3000, omega_var == 0.25^2)

ggplot(data, aes(x = N, y = p_beta, group = Tfull, color = as.factor(Tfull))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_beta,
                  ymax = bino_ci_exact_upper_beta,
                  fill = as.factor(Tfull)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nβ = 3,\nρ = -0.9\nω² = 0.25^2\n\n# of Reps:\n~3000\n\nTime T:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nβ = 3,\nρ = -0.9\nω² = 0.25^2\n\n# of Reps:\n~3000\n\nTime T:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Coverage Probability for Mean β vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = c(0.825, 1)) + 
  scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 100)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")


## Compare Omega vs. Beta coverage rate + Vary Omega -------------------------------
data = original %>% filter(rho == -0.9, Tfull == 2, n_reps == 3000, omega_var == 4^2, beta_free %in% c(1, 3))

ggplot(data, aes(x = N, group = beta_free)) + 
  # First set: original p_beta
  geom_line(aes(y = p_omega_rescaled, color = as.factor(beta_free))) + 
  geom_point(aes(y = p_omega_rescaled, color = as.factor(beta_free)), size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_omega_rescaled, ymax = bino_ci_exact_upper_omega_rescaled, fill = as.factor(beta_free)), 
              alpha = 0.3, color = NA) +
  
  # First color & fill scale: blue-green
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Omega Parameter") +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Omega Parameter") +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_line(aes(y = p_beta_rescaled, color = as.factor(beta_free))) + 
  geom_point(aes(y = p_beta_rescaled, color = as.factor(beta_free)), shape = 1, size = 1.2) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_beta_rescaled, ymax = bino_ci_exact_upper_beta_rescaled, fill = as.factor(beta_free)), 
              alpha = 0.3, color = NA) +
  
  scale_color_manual(
    values = c("red", "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  scale_fill_manual(
    values = c("red", "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  guides(
    color = guide_legend( 
      title = paste0("Fixed\nParams:\nT = ", unique(data$Tfull),
                                 "\nω² = ", unique(data$omega_var),
                                 "\nρ = ", unique(data$rho),
                                 "\n\n# of Reps:\n~3000\n\nMean β and\nRescaled:"
      )
    ),
    fill = guide_legend(
      title = paste0("Fixed\nParams:\nT = ", unique(data$Tfull),
                     "\nω² = ", unique(data$omega_var),
                     "\nρ = ", unique(data$rho),
                     "\n\n# of Reps:\n~3000\n\nMean β and\nRescaled:"
      )
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = expression("CI Comparison for Mean β1 Estimates vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability"
  ) +
  coord_cartesian(ylim = c(0.75, 1)) + 
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")
data = original %>% filter(rho == -0.9, Tfull == 4, n_reps == 3000, beta_free == 4) %>% filter(omega_var %in% c(0.25^2, 1, 4^2))

ggplot(data, aes(x = N, group = omega_var)) + 
  # First set: original p_beta
  geom_line(aes(y = p_omega_rescaled, color = as.factor(omega_var))) + 
  geom_point(aes(y = p_omega_rescaled, color = as.factor(omega_var)), size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_omega_rescaled, ymax = bino_ci_exact_upper_omega_rescaled, fill = as.factor(omega_var)), 
              alpha = 0.3, color = NA) +
  
  # First color & fill scale: blue-green
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Omega Parameter") +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Omega Parameter") +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_line(aes(y = p_beta_rescaled, color = as.factor(omega_var))) + 
  geom_point(aes(y = p_beta_rescaled, color = as.factor(omega_var)), shape = 1, size = 1.2) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_beta_rescaled, ymax = bino_ci_exact_upper_beta_rescaled, fill = as.factor(omega_var)), 
              alpha = 0.3, color = NA) +
  
  scale_color_manual(
    values = c("red", "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  scale_fill_manual(
    values = c("red", "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  guides(
    color = guide_legend( 
      title = paste0("Fixed\nParams:\nT = ", unique(data$Tfull),
                     "\nβ = ", unique(data$beta_free),
                     "\nρ = ", unique(data$rho),
                     "\n\n# of Reps:\n~3000\n\nMean β and\nRescaled:"
      )
    ),
    fill = guide_legend(
      title = paste0("Fixed\nParams:\nT = ", unique(data$Tfull),
                     "\nβ = ", unique(data$beta_free),
                     "\nρ = ", unique(data$rho),
                     "\n\n# of Reps:\n~3000\n\nMean β and\nRescaled:"
      )
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = expression("CI Comparison for Mean β1 Estimates vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability"
  ) +
  coord_cartesian(ylim = c(0.75, 1)) + 
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")
## Compare coverage rates between SE types while varying Omega -------------------------------
data = original %>% filter(rho == 0.9, Tfull == 2, n_reps == 3000, beta_free == 4) %>%
  filter(beta_free %in% c(0,1,2,3,4), omega_var %in% c(1^2))

ggplot(data, aes(x = N, group = omega_var)) + 
  # First set: original p_beta
  geom_line(aes(y = p_omega, color = as.factor(omega_var))) + 
  geom_point(aes(y = p_omega, color = as.factor(omega_var)), size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_omega, ymax = bino_ci_exact_upper_omega, fill = as.factor(omega_var)), 
              alpha = 0.3, color = NA) +
  
  # First color & fill scale: blue-green
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Var ω² and\nModel-based") +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Var ω² and\nModel-based") +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_line(aes(y = p_omega_rescaled, color = as.factor(omega_var))) + 
  geom_point(aes(y = p_omega_rescaled, color = as.factor(omega_var)), shape = 1, size = 1.2) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_omega_rescaled, ymax = bino_ci_exact_upper_omega_rescaled, fill = as.factor(omega_var)), 
              alpha = 0.3, color = NA) +
  
  scale_color_manual(
    values = c("purple", "red",  "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  scale_fill_manual(
    values = c("purple", "red",  "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  guides(
    color = guide_legend(
      title = paste0("Fixed\nParams:\nT = ", unique(data$Tfull),
                     "\nβ = ", unique(data$beta_free),
                     "\nρ = ", unique(data$rho),
                     "\n\n# of Reps:\n~3000\n\nVar ω² and\nRescaled:"
        )
    ),
    fill = guide_legend(
      title = paste0("Fixed\nParams:\nT = ", unique(data$Tfull),
                     "\nβ = ", unique(data$beta_free),
                     "\nρ = ", unique(data$rho),
                     "\n\n# of Reps:\n~3000\n\nVar ω² and\nRescaled:"
      )
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = expression("CI Comparison for Mean β1 Estimates vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability"
  ) +
  coord_cartesian(ylim = c(0.5, 1)) + 
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")


# Rho Parameter section ---------------------------------------------------
## Rho bias plot for varying beta ----
data = original %>% filter(rho == 0.6, Tfull == 2, n_reps == 3000, omega_var == 1) %>%
  filter(beta_free %in% c(0,2,4))
alpha = 0.05
z_value = qnorm(1-alpha/2)
sqrtN_to_divide_by = sqrt(as.numeric(data$success_beta + data$failure_beta))
data = data %>%
  mutate(beta_lower = beta_bias - z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         beta_upper = beta_bias + z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         omega_lower = omega_bias - z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         omega_upper = omega_bias + z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         rho_lower = rho_bias - z_value * sd_of_rho_estimates/sqrtN_to_divide_by,
         rho_upper = rho_bias + z_value * sd_of_rho_estimates/sqrtN_to_divide_by)

if(data$rho_bias[which(data$N == 300)[1]] > 0){
  yaxis = c(0, 0.25)
} else {
  yaxis = c(-0.25, 0)
}

ggplot(data, aes(x = N, y = rho_bias, group = beta_free, color = as.factor(beta_free))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = rho_lower, ymax = rho_upper, fill = as.factor(beta_free)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Bias of Mean β1 Estimate vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Bias of Mean β1",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = yaxis) + 
  scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 100)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")

## Rho bias plot for varying omega ----
data = original %>% filter(rho == 0.6, Tfull == 2, n_reps == 3000, beta_free == 2) %>%
  filter(beta_free %in% c(0, 1, 2, 3, 4), omega_var %in% c(0.25^2, 1, 4^2))
alpha = 0.05
z_value = qnorm(1-alpha/2)
sqrtN_to_divide_by = sqrt(as.numeric(data$success_beta + data$failure_beta))
data = data %>%
  mutate(beta_lower = beta_bias - z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         beta_upper = beta_bias + z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         omega_lower = omega_bias - z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         omega_upper = omega_bias + z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         rho_lower = rho_bias - z_value * sd_of_rho_estimates/sqrtN_to_divide_by,
         rho_upper = rho_bias + z_value * sd_of_rho_estimates/sqrtN_to_divide_by)
if(data$rho_bias[which(data$N == 300)[1]] > 0){
  yaxis = c(0, 0.25)
} else {
  yaxis = c(-0.25, 0)
}

ggplot(data, aes(x = N, y = rho_bias, group = omega_var, color = as.factor(omega_var))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = rho_lower, ymax = rho_upper, fill = as.factor(omega_var)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Bias of Mean β1 Estimate vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Bias of Mean β1",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = yaxis) + 
  scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 100)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")


## N vs. Rho Bias + Vary rho -----------------------------------------------
data <- original %>%
  filter(omega_var == 4^2, Tfull == 3, n_reps == 3000, beta_free == 4)

alpha <- 0.05
z_value <- qnorm(1 - alpha/2)
sqrtN_to_divide_by <- sqrt(data$success_beta + data$failure_beta)
data <- data %>%
  mutate(
    beta_lower  = beta_bias  - z_value * sd_of_beta_estimates  / sqrtN_to_divide_by,
    beta_upper  = beta_bias  + z_value * sd_of_beta_estimates  / sqrtN_to_divide_by,
    omega_lower = omega_bias - z_value * sd_of_omega_estimates / sqrtN_to_divide_by,
    omega_upper = omega_bias + z_value * sd_of_omega_estimates / sqrtN_to_divide_by,
    rho_lower   = rho_bias   - z_value * sd_of_rho_estimates   / sqrtN_to_divide_by,
    rho_upper   = rho_bias   + z_value * sd_of_rho_estimates   / sqrtN_to_divide_by
  )

yaxis <- c(-0.2, 0.2)
brks  <- round(seq(-0.9, 0.9, by = 0.3), 1)

ggplot(
  data,
  aes(x = N, y = rho_bias, group = rho, color = rho)  # rho as numeric for gradient
) +
  geom_ribbon(
    aes(ymin = rho_lower, ymax = rho_upper, fill = rho, group = rho),  # group ribbons by rho
    alpha = 0.25, color = NA, show.legend = FALSE
  ) +
  geom_line() +
  geom_point(size = 1) +
  scale_color_gradientn(
    colors = c("#08306B", "#4EA3D9", "#2FBF71", "#F9D057", "#D7191C"),
    values = scales::rescale(c(-0.9, -0.3, 0, 0.3, 0.9)),
    limits = c(-0.9, 0.9),
    breaks = brks,
    labels = scales::label_number(accuracy = 0.1, trim = TRUE),
    name = expression(rho)
  ) +
  scale_fill_gradientn(  # match ribbon fill to line color
    colors = c("#08306B", "#4EA3D9", "#2FBF71", "#F9D057", "#D7191C"),
    values = scales::rescale(c(-0.9, -0.3, 0, 0.3, 0.9)),
    limits = c(-0.9, 0.9),
    guide = "none"
  ) +
  guides(
    color = guide_colorbar(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nω² = ", unique(data$omega_var),
        "\nβ = ", unique(data$beta_free),
        "\n\n# of Reps:\n~3000\n\nAR(1) ρ:\n"
      ),
      barheight = grid::unit(140, "pt"),
      barwidth  = grid::unit(10,  "pt")
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = expression("Bias of AR(1) ρ Estimate vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Rho (ρ) Bias"
  ) +
  coord_cartesian(ylim = yaxis) +
  scale_x_continuous(breaks = seq(0, 1000, by = 100)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")



## Rho  vs. Bias Rho + any param in color -----
data = original %>% filter(Tfull == 2, N == 500, n_reps == 3000, omega_var == 1)
ggplot(data, aes(x = rho, y = rho_bias, group = beta_free, color = as.factor(beta_free))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  guides(
    color = guide_legend(
      title = paste0("Fixed\nParams:\nT = ",unique(data$Tfull),
                     "\nN = ",unique(data$N),"\nω² = ",
                     unique(data$omega_var),
                     "\n\n# of Reps:\n~3000\n\nMean β1:")
    ),
    fill = guide_legend(
      title = paste0("Fixed\nParams:\nT = ",
                     unique(data$Tfull),"\nN = ",
                     unique(data$N),"\nω² = ",
                     unique(data$omega_var),
                     "\n\n# of Reps:\n~3000\n\nMean β1:")
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Rho Bias vs. Rho and Omega Variance",
    x = "Rho",
    y = "Rho Bias",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(-0.1, 0.1))

# ---

what_to_plot = list(omega_var = 0.5^2,
                    N = 500,
                    Tfull = 2)

data_dense  <- original %>%
  filter(Tfull == what_to_plot$T, N == what_to_plot$N, n_reps == 3000, omega_var == what_to_plot$omega_var) %>%
  filter(beta_free <= 1) %>%
  filter(!(beta_free %in% c(0.2, 0.6)))

data_coarse <- original %>%
  filter(Tfull == what_to_plot$T, N == what_to_plot$N, n_reps == 3000, omega_var == what_to_plot$omega_var) %>% filter(beta_free %in% c(1,2,3,4))

ggplot() +
  geom_line(data = data_coarse,
            aes(x = rho, y = rho_bias, group = beta_free, color = as.factor(beta_free)),
            size = 0.9) +
  geom_point(data = data_coarse,
             aes(x = rho, y = rho_bias, color = as.factor(beta_free)),
             size = 1.8) +
  scale_color_manual(
    name = NULL,
    values = c(`2` = "#FDE725", `3` = "#FDB366", `4` = "#D73027"),
    guide = guide_legend(order = 2)   # put this legend BELOW
  ) +
  ggnewscale::new_scale_color() +
  geom_line(data = data_dense,
            aes(x = rho, y = rho_bias, group = beta_free, color = as.factor(beta_free)),
            size = 0.75) +
  geom_point(data = data_dense,
             aes(x = rho, y = rho_bias, color = as.factor(beta_free)),
             size = 1.5) +
  scale_color_viridis_d(
    option = "D", begin = 0.4, end = 0.9,
    name = paste0(
      "Fixed\nParams:\nT = ", unique(data_dense$Tfull),
      "\nN = ", unique(data_dense$N),
      "\nω² = ", unique(data_dense$omega_var),
      "\n\n# of Reps:\n~3000\n\nMean β:"
    ),
    guide = guide_legend(order = 1) 
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0, "cm"),
    plot.title = element_text(hjust = 0.1)
  ) +
  labs(
    title = "Rho vs. Rho Bias",
    x = "Rho",
    y = "Rho Bias"
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(-0.025, 0.025)) # 0.1/0.04 for var=1 & T=2/3, 4 for Var=16 & both T

## Rho vs. Rho Empirical SE -----
data = original %>% filter(Tfull == 2, N == 500, n_reps == 3000, omega_var == 16)
ggplot(data, aes(x = rho, y = sd_of_rho_estimates, group = beta_free, color = as.factor(beta_free))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  guides(
    color = guide_legend(
      title = paste0("Fixed\nParams:\nT = ",unique(data$Tfull),"\nN = ",unique(data$N),"\nω² = ",unique(data$omega_var),"\n\n# of Reps:\n~3000\n\nMean β1:")
    ),
    fill = guide_legend(
      title = paste0("Fixed\nParams:\nT = ",unique(data$Tfull),"\nN = ",unique(data$N),"\nω² = ",unique(data$omega_var),"\n\n# of Reps:\n~3000\n\nMean β1:")
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Rho vs. Omega Bias",
    x = "Rho",
    y = "Omega Bias",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0, 0.4))

# ---

what_to_plot = list(omega_var = 2^2,
                    N = 500,
                    Tfull = 2)

data_dense  <- original %>%
  filter(Tfull == what_to_plot$T, N == what_to_plot$N, n_reps == 3000, omega_var == what_to_plot$omega_var) %>%
  filter(beta_free <= 1) %>%
  filter(!(beta_free %in% c(0.2, 0.6)))

data_coarse <- original %>%
  filter(Tfull == what_to_plot$T, N == what_to_plot$N, n_reps == 3000, omega_var == what_to_plot$omega_var) %>% filter(beta_free %in% c(1,2,3,4))

ggplot() +
  geom_line(data = data_coarse,
            aes(x = rho, y = sd_of_rho_estimates, group = beta_free, color = as.factor(beta_free)),
            size = 0.9) +
  geom_point(data = data_coarse,
             aes(x = rho, y = sd_of_rho_estimates, color = as.factor(beta_free)),
             size = 1.8) +
  scale_color_manual(
    name = NULL,
    values = c(`2` = "#FDE725", `3` = "#FDB366", `4` = "#D73027"),
    guide = guide_legend(order = 2)   # put this legend BELOW
  ) +
  ggnewscale::new_scale_color() +
  geom_line(data = data_dense,
            aes(x = rho, y = sd_of_rho_estimates, group = beta_free, color = as.factor(beta_free)),
            size = 0.75) +
  geom_point(data = data_dense,
             aes(x = rho, y = sd_of_rho_estimates, color = as.factor(beta_free)),
             size = 1.5) +
  scale_color_viridis_d(
    option = "D", begin = 0.2, end = 0.7,
    name = paste0(
      "Fixed\nParams:\nT = ", unique(data_dense$Tfull),
      "\nN = ", unique(data_dense$N),
      "\nω² = ", unique(data_dense$omega_var),
      "\n\n# of Reps:\n~3000\n\nMean β1:"
    ),
    guide = guide_legend(order = 1)   # put this legend ABOVE
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0, "cm"),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Rho vs. Rho Empirical SE",
    x = "Rho",
    y = "Empirical Standard Error"
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0, 0.3)) # 0.1/0.04 for var=1 & T=2/3, 4 for Var=16 & both T

# --- lets try a percentage comparison (very messy)
library(scales)

# Preliminary settings here
what <- list(omega_var = 1^2, N = 500)
metric_col <- "rho_bias"   # or "sd_of_rho_estimates"
baseline_T <- 2
compare_T  <- c(3, 4)
beta_range <- c(0, 0.4, 0.8, 1.2)

df <- original %>%
  filter(N == what$N, n_reps == 3000, omega_var == what$omega_var, beta_free %in% beta_range,
         Tfull %in% c(baseline_T, compare_T))

base <- df %>%
  filter(Tfull == baseline_T, beta_free %in% beta_range,) %>%
  select(rho, beta_free, base_metric = !!sym(metric_col))

comp <- df %>%
  filter(Tfull %in% compare_T) %>%
  left_join(base, by = c("rho", "beta_free")) %>%
  mutate(
    pct_reduction = if_else(
      abs(base_metric) > 0,
      100 * (abs(base_metric) - abs(!!sym(metric_col))) / abs(base_metric),
      NA_real_
    )
  )


# plot % reduction vs rho, colored by beta, facets for T 
ggplot(
  comp,
  aes(x = rho, y = pct_reduction, color = as.factor(beta_free), group = beta_free)
) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(size = 0.9) +
  geom_point(size = 1.8) +
  facet_wrap(~ Tfull, labeller = as_labeller(function(x)
    paste0("T = ", x, " vs T = ", baseline_T))) +
  scale_color_viridis_d(
    option = "D", begin = 0.2, end = 0.7,
    name = paste0(
      "Fixed\nParams:\nT = ", paste(compare_T, collapse = ", "),
      "\nN = ", unique(df$N),
      "\nω² = ", unique(df$omega_var),
      "\n\n# of Reps:\n~3000\n\nMean β1:"
    )
  ) +
  scale_x_continuous(breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9), limits = c(-0.9, 0.9)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0, "cm"),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Percentage reduction in absolute ρ-bias relative to T = 2",
    x = "Rho",
    y = "Percent reduction"
  )

# settings for empirical SE comparison
what <- list(omega_var = 1^2, N = 500)
baseline_T <- 2
compare_T  <- c(3, 4)              # which T to compare against baseline
metric_col <- "sd_of_rho_estimates"
beta_range <- c(0, 1, 2, 3, 4)

df <- original %>%
  filter(N == what$N, n_reps == 3000, omega_var == what$omega_var, beta_free %in% beta_range,
         Tfull %in% c(baseline_T, compare_T))

base <- df %>%
  filter(Tfull == baseline_T, beta_free %in% beta_range) %>%
  select(rho, beta_free, base_metric = !!sym(metric_col))

comp <- df %>%
  filter(Tfull %in% compare_T) %>%
  left_join(base, by = c("rho", "beta_free")) %>%
  mutate(
    pct_reduction = if_else(
      base_metric > 0,
      100 * (base_metric - !!sym(metric_col)) / base_metric,  # positive = improvement
      NA_real_
    )
  )

ggplot(
  comp,
  aes(x = rho, y = pct_reduction, color = as.factor(beta_free), group = beta_free)
) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(size = 0.9) +
  geom_point(size = 1.8) +
  facet_wrap(~ Tfull, labeller = as_labeller(function(x)
    paste0("T = ", x, " vs T = ", baseline_T))) +
  scale_color_viridis_d(
    option = "D", begin = 0.2, end = 0.7,
    name = paste0(
      "Fixed\nParams:\nT = ", paste(compare_T, collapse = ", "),
      "\nN = ", unique(df$N),
      "\nω² = ", unique(df$omega_var),
      "\n\n# of Reps:\n~3000\n\nMean β1:"
    )
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0, "cm"),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Percent reduction in empirical SE of ρ relative to T = 2",
    x = "Rho",
    y = "Percent reduction"
  )

what <- list(omega_var = 1^2, N = 500)
pairs <- list(c(2, 3), c(3, 4))   # baseline and comparison T pairs
metric_col <- "sd_of_rho_estimates"

plots <- list()

for (p in seq_along(pairs)) {
  baseline_T <- pairs[[p]][1]
  compare_T  <- pairs[[p]][2]
  
  df <- original %>%
    filter(N == what$N, n_reps == 3000, omega_var == what$omega_var,
           Tfull %in% c(baseline_T, compare_T))
  
  # baseline
  base <- df %>%
    filter(Tfull == baseline_T) %>%
    select(rho, beta_free, base_metric = !!sym(metric_col))
  
  # compare
  comp <- df %>%
    filter(Tfull == compare_T) %>%
    left_join(base, by = c("rho", "beta_free")) %>%
    mutate(
      pct_reduction = if_else(
        base_metric > 0,
        100 * (base_metric - !!sym(metric_col)) / base_metric,
        NA_real_
      )
    )
  
  # plot
  plots[[p]] <- ggplot(
    comp,
    aes(x = rho, y = pct_reduction, color = as.factor(beta_free), group = beta_free)
  ) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(size = 0.9) +
    geom_point(size = 1.8) +
    scale_color_viridis_d(
      option = "D", begin = 0.2, end = 0.7,
      name = paste0("Fixed\nParams:\nT = ", compare_T,
                    "\nN = ", unique(df$N),
                    "\nω² = ", unique(df$omega_var),
                    "\n\n# of Reps:\n~3000\n\nMean β1:")
    ) +
    scale_x_continuous(
      breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
      limits = c(-0.9, 0.9)
    ) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(0, "cm"),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = paste0("Percent reduction in empirical SE of ρ: T = ",
                     baseline_T, " → T = ", compare_T),
      x = "Rho",
      y = "Percent reduction"
    )
}

# Display them
plots[[1]]  # T=2 → T=3
plots[[2]]  # T=3 → T=4

## Rho Model-based SE vs. Empirical SE (Vary beta) -----------------------------------
data = original %>% filter(rho == 0.9, Tfull == 2, n_reps == 3000, omega_var == 1^2) %>% filter(beta_free %in% c(0,2,4))

ggplot(data, aes(x = N, color = as.factor(beta_free))) + 
  # Model-based SE (solid line)
  geom_line(aes(y = avg_of_rho_se, linetype = "Avg. Model-\nbased")) + 
  geom_point(aes(y = avg_of_rho_se), size = 1) +
  
  # Empirical SE (dashed line)
  geom_line(aes(y = sd_of_rho_estimates, linetype = "Empirical")) +
  geom_point(aes(y = sd_of_rho_estimates), shape = 17, size = 1) +
  
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Mean β1") +
  
  # Add linetype legend
  scale_linetype_manual(
    name = "SE Type",
    values = c("Avg. Model-\nbased" = "solid", "Empirical" = "dashed")
  ) +
  
  guides(
    color = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nω² = ", unique(data$omega_var),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nMean β1:"
      )
    ),
    fill = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nω² = ", unique(data$omega_var),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nMean β1:"
      )
    )
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  
  labs(
    title = expression("SE Comparison for Omega Variance Estimates vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Standard Error"
  ) +
  
  coord_cartesian(ylim = c(0, 0.3)) + 
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")


## Are the situations where model-based SE < empirical SE?  -----------------
df <- original %>%
  filter(Tfull == 2, n_reps == 3000, N == 500) %>%
  mutate(
    se_emp   = sd_of_rho_estimates,
    se_mod   = avg_of_rho_se,
    ratio    = se_emp / se_mod,
    emp_gt   = se_emp > se_mod,
    abs_gap  = se_emp - se_mod,
    rel_gap  = (se_emp - se_mod) / se_mod
  )

# combos where empirical > model-based
cases_emp_gt <- df %>%
  filter(emp_gt) %>%
  arrange(desc(ratio))

# argest gaps
head(cases_emp_gt %>% select(rho, omega_var, beta_free, se_emp, se_mod, ratio, abs_gap, rel_gap), 20)



# by rho, omega_var, and beta_free:
sum_rho_omega_beta <- df %>%
  group_by(rho, omega_var, beta_free) %>%
  summarise(
    n = n(),
    emp_gt = any(emp_gt),
    mean_ratio = mean(ratio, na.rm = TRUE),
    max_ratio  = max(ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(emp_gt) %>%
  arrange(desc(max_ratio))

## Rho Model-based SE vs. Empirical SE (Vary Omega) -----------------------------------
data = original %>% filter(rho == 0.6, Tfull == 2, n_reps == 3000, beta_free == 2) %>% filter(omega_var %in% c(0.25^2, 1, 4^2))

ggplot(data, aes(x = N, color = as.factor(omega_var))) + 
  # Model-based SE (solid line)
  geom_line(aes(y = avg_of_rho_se, linetype = "Avg. Model-\nbased")) + 
  geom_point(aes(y = avg_of_rho_se), size = 1) +
  
  # Empirical SE (dashed line)
  geom_line(aes(y = sd_of_rho_estimates, linetype = "Empirical")) +
  geom_point(aes(y = sd_of_rho_estimates), shape = 17, size = 1) +
  
  scale_color_viridis_d(option = "D", begin = 0.4, end = 0.9, name = "Mean β1") +
  
  # Add linetype legend
  scale_linetype_manual(
    name = "SE Type",
    values = c("Avg. Model-\nbased" = "solid", "Empirical" = "dashed")
  ) +
  
  guides(
    color = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nβ = ", unique(data$beta_free),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nVar ω²:"
      )
    ),
    fill = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nβ = ", unique(data$beta_free),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nVar ω²:"
      )
    )
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  
  labs(
    title = expression("SE Comparison for Omega Variance Estimates vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Standard Error"
  ) +
  
  coord_cartesian(ylim = c(0, 1)) + 
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")


## Rho Bias figure ----
data = original %>% filter(rho == 0, Tfull == 2, n_reps == 3000, omega_var == 1) %>% filter(beta_free %in% c(0,2,4))
alpha = 0.05
z_value = qnorm(1-alpha/2)
sqrtN_to_divide_by = sqrt(as.numeric(data$success_beta + data$failure_beta))
data = data %>%
  mutate(beta_lower = beta_bias - z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         beta_upper = beta_bias + z_value * sd_of_beta_estimates/sqrtN_to_divide_by,
         omega_lower = omega_bias - z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         omega_upper = omega_bias + z_value * sd_of_omega_estimates/sqrtN_to_divide_by,
         rho_lower = rho_bias - z_value * sd_of_rho_estimates/sqrtN_to_divide_by,
         rho_upper = rho_bias + z_value * sd_of_rho_estimates/sqrtN_to_divide_by)

ggplot(data, aes(x = N, y = rho_bias, group = beta_free, color = as.factor(beta_free))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = rho_lower, ymax = rho_upper, fill = as.factor(beta_free)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nT = 2,\nω² = 1,\nρ = -0.9\n\n# of Reps:\n~3000\n\nMean β1:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Bias of Mean β1 Estimate vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Bias of Mean β1",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = c(0, 0.5)) + 
  scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 100)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")



## Rho vs. Rho Bias (+ Sample size in Color) -----------------------------------
data = original %>% filter(Tfull == 2, n_reps == 3000, beta_free == 2, omega_var == 1^2) %>% filter(N %in% c(50, 80, 100, 150, 200))
ggplot(data, aes(x = rho, y = rho_bias, group = N, color = as.factor(N))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  guides(
    color = guide_legend(
      title = paste0("Fixed\nParams:\nT = ",unique(data$Tfull),"\nN = ",unique(data$N),"\nω² = ",unique(data$omega_var),"\n\n# of Reps:\n~3000\n\nMean β1:")
    ),
    fill = guide_legend(
      title = paste0("Fixed\nParams:\nT = ",unique(data$Tfull),"\nN = ",unique(data$N),"\nω² = ",unique(data$omega_var),"\n\n# of Reps:\n~3000\n\nMean β1:")
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Rho vs. Omega Bias",
    x = "Rho",
    y = "Omega Bias",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(-0.25, 0.25))

## Rho vs. Rho Coverage Probs with beta increasing  -----
data = original %>% filter(Tfull == 2, N == 300, n_reps == 3000, omega_var == 0.25^2) %>%
  filter(omega_var %in% c(0.25^2, 1, 4^2), beta_free %in% c(0, 1, 2, 3, 4))
ggplot(data, aes(x = rho, y = p_rho_rescaled, group = beta_free, color = as.factor(beta_free))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  guides(
    color = guide_legend(
      title = paste0("Fixed\nParams:\nT = ", unique(data$Tfull),
                     "\nN = ", unique(data$N),
                     "\nω² = ", unique(data$omega_var),
                     "\n\n# of Reps:\n~3000\n\nMean β:"
      )
                     
    ),
    fill = guide_legend(
      title = paste0("Fixed\nParams:\nT =", unique(data$Tfull),
                     "\nN = ", unique(data$N),
                     "\nω² = ", unique(data$omega_var),
                     "\n\n# of Reps:\n~3000\n\nMean β:"
      )
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Coverage Probability for AR(1) coefficient ρ vs. Rho",
    x = "Rho",
    y = "Coverage Probability",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0.85, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")


## Rho vs. Rho Coverage Probs with Omega increasing  -----
data = original %>% filter(Tfull == 2, N == 100, n_reps == 3000, beta_free == 0) %>%
  filter(omega_var %in% c(0.25^2, 1, 4^2), beta_free %in% c(0, 1, 2, 3, 4))
ggplot(data, aes(x = rho, y = p_rho_rescaled, group = omega_var, color = as.factor(omega_var))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.4, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.4, end = 0.9) +
  guides(
    color = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nN = ", unique(data$N),
        "\nβ = ", unique(data$beta_free),
        "\n\n# of Reps:\n~3000\n\nVar ω²:"
      )
    ),
    fill = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nN = ", unique(data$N),
        "\nβ = ", unique(data$beta_free),
        "\n\n# of Reps:\n~3000\n\nVar ω²:"
      )
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Coverage Probability for AR(1) coefficient ρ vs. Rho",
    x = "Rho",
    y = "Coverage Probability",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0.85, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")

## Rho vs. Rho Coverage Probs + Fix Omega and vary N  -----
data = original %>% filter(Tfull == 2, N %in% c(150, 200, 300, 500, 1000),
                           n_reps == 3000, beta_free == 1, omega_var == 4^2) %>%
  filter(omega_var %in% c(0.25^2, 1, 4^2), beta_free %in% c(0, 1, 2, 3, 4))
ggplot(data, aes(x = rho, y = p_rho_rescaled, group = N, color = as.factor(N))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  guides(
    color = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nω² = ", unique(data$omega_var),
        "\nβ = ", unique(data$beta_free),
        "\n\n# of Reps:\n~3000\n\nSample\nSize N:"
      )
    ),
    fill = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nω² = ", unique(data$omega_var),
        "\nβ = ", unique(data$beta_free),
        "\n\n# of Reps:\n~3000\n\nSample\nSize N:"
      )
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Comparison of Omega Coverage Probability vs. AR(1) coefficient ρ",
    x = "Rho",
    y = "Coverage Probability",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0.75, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")


## Rho vs. Rho Coverage Probs + Fix Omega and vary N and T  -----
data = original %>% filter(Tfull %in% c(3, 4) , N %in% c(150, 200, 300, 500, 1000),
                           n_reps == 3000, beta_free == 3, omega_var == 4^2) %>%
  filter(omega_var %in% c(0.25^2, 1, 4^2), beta_free %in% c(0, 1, 2, 3, 4))
ggplot(data, aes(
  x = rho,
  y = p_rho_rescaled,
  group = interaction(N, Tfull),  # group by both
  color = as.factor(N),           # N as color
  linetype = as.factor(Tfull)     # T as line type
)) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  scale_linetype_manual(
    values = setNames(
      c("dashed", "solid"),  # first is for lower T, second for higher T
      as.character(c(min(unique(data$Tfull)), max(unique(data$Tfull))))
    )
  ) +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  guides(
    color = guide_legend(order = 1,
                         title = paste0(
                           "Fixed\nParams:",
                           "\nω² = ", unique(data$omega_var),
                           "\nβ = ", unique(data$beta_free),
                           "\n\n# of Reps:\n~3000\n\nSample\nSize N:"
                         )
    ),
    linetype = guide_legend(
      title = "Time T:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Comparison of Coverage Success over AR(1) coefficient ρ for different N and T",
    x = "Rho",
    y = "Coverage Probability",
    color = "N",
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0.75, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")

grid.arrange(p1, p2, ncol = 2) 
## N vs. Rho coverage plots + Vary Omega -----------------------------------------
data = original %>% filter(rho == 0.3, Tfull == 2, n_reps == 3000, beta_free == 0) %>%
  filter(omega_var %in% c(0.25^2, 1, 4^2))

ggplot(data, aes(x = N, y = p_rho_rescaled, group = omega_var, color = as.factor(omega_var))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_rho_rescaled, ymax = bino_ci_exact_upper_rho_rescaled, fill = as.factor(omega_var)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nβ = ", unique(data$beta_free),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nVar ω²:"
      )
    ),
    fill = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nβ = ", unique(data$beta_free),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nVar ω²:"
      )
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Coverage Probability for Mean β vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = c(0.75, 1)) + 
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black") 


## N vs. Rho Coverage Probability + Varying beta ------------------------------------------
data = original %>% filter(rho == 0, Tfull == 2, n_reps == 3000, omega_var == 0.25^2) %>%
  filter(omega_var %in% c(0.25^2, 1, 4^2), beta_free %in% c(0, 2, 4))

ggplot(data, aes(x = N, y = p_rho_rescaled, group = beta_free, color = as.factor(beta_free))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_rho_rescaled, ymax = bino_ci_exact_upper_rho_rescaled, fill = as.factor(beta_free)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nω² = ", unique(data$omega_var),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nMean β:"
      )
    ),
    fill = guide_legend(
      title = paste0(
        "Fixed\nParams:\nT = ", unique(data$Tfull),
        "\nω² = ", unique(data$omega_var),
        "\nρ = ", unique(data$rho),
        "\n\n# of Reps:\n~3000\n\nMean β:"
      )
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Coverage Probability for Mean β vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = c(0.75, 1)) + 
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black") 


## Other param vs. fixed N coverage rate (not rescaled) ------
data = original %>% filter(N == 200, n_reps == 3000, beta_free == 2, omega_var == 1)
ggplot(data, aes(x = rho, y = p_rho, group = Tfull, color = as.factor(Tfull))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) + 
  scale_color_viridis_d(option = "D", begin = 0.4, end = 0.9) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.4, end = 0.9) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nβ = 2,\nN = 200,\nω² = 1\n\n# of Reps:\n~3000\n\nTime T:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nβ = 2,\nN = 200,\nω² = 1\n\n# of Reps:\n~3000\n\nTime T:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = "Coverage Probability for Mean β vs. AR(1) coefficient ρ",
    x = "Mean of Beta Coefficient",
    y = "Coverage Probability",
    color = expression(omega^2),
    fill = expression(omega^2)
  ) +
  scale_x_continuous(
    breaks = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9),
    limits = c(-0.9, 0.9)
  ) +
  coord_cartesian(ylim = c(0.85, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")
## Varying T in coverage plots (not rescaled) -----------------------------------------
data = original %>% filter(rho == -0.9, beta_free == 3, n_reps == 3000, omega_var == 0.25^2)

ggplot(data, aes(x = N, y = p_rho, group = Tfull, color = as.factor(Tfull))) + 
  geom_line() + 
  geom_point(size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_rho,
                  ymax = bino_ci_exact_upper_rho,
                  fill = as.factor(Tfull)), alpha = 0.3, color = NA) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7) +  # nicer gradient
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7) +
  guides(
    color = guide_legend(
      title = "Fixed\nParams:\nβ = 3,\nρ = -0.9\nω² = 0.25^2\n\n# of Reps:\n~3000\n\nTime T:"
    ),
    fill = guide_legend(
      title = "Fixed\nParams:\nβ = 3,\nρ = -0.9\nω² = 0.25^2\n\n# of Reps:\n~3000\n\nTime T:"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) + 
  labs(
    title = expression("Coverage Probability for Mean β vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability",
    color = "Mean β1",
    fill = "Mean β1"
  ) +
  coord_cartesian(ylim = c(0.825, 1)) + 
  scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 100)) +
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")


## Compare Rho vs. Beta/Omega coverage rate + Vary Omega -------------------------------
rho_range = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
for (i in 1:length(rho_range)){
data = original %>% filter(rho == rho_range[i], Tfull == 2, n_reps == 3000,
                           omega_var == 4^2, beta_free %in% c(1, 3))

p <- ggplot(data, aes(x = N, group = beta_free)) + 
  # First set: original p_beta
  geom_line(aes(y = p_rho_rescaled, color = as.factor(beta_free))) + 
  geom_point(aes(y = p_rho_rescaled, color = as.factor(beta_free)), size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_rho_rescaled, ymax = bino_ci_exact_upper_rho_rescaled, fill = as.factor(beta_free)), 
              alpha = 0.3, color = NA) +
  
  # First color & fill scale: blue-green
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Rho Parameter") +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Rho Parameter") +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_line(aes(y = p_omega_rescaled, color = as.factor(beta_free))) + 
  geom_point(aes(y = p_omega_rescaled, color = as.factor(beta_free)), shape = 1, size = 1.2) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_omega_rescaled, ymax = bino_ci_exact_upper_omega_rescaled, fill = as.factor(beta_free)), 
              alpha = 0.3, color = NA) +
  
  scale_color_manual(
    values = c("red", "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  scale_fill_manual(
    values = c("red", "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  guides(
    color = guide_legend( 
      title = paste0("Fixed\nParams:\nT = ", unique(data$Tfull),
                     "\nω² = ", unique(data$omega_var),
                     "\nρ = ", unique(data$rho),
                     "\n\n# of Reps:\n~3000\n\nMean β/Omega Var:"
      )
    ),
    fill = guide_legend(
      title = paste0("Fixed\nParams:\nT = ", unique(data$Tfull),
                     "\nω² = ", unique(data$omega_var),
                     "\nρ = ", unique(data$rho),
                     "\n\n# of Reps:\n~3000\n\nMean β/Omega Var:"
      )
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = expression("CI Comparison for Mean β1 Estimates vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability"
  ) +
  coord_cartesian(ylim = c(0.75, 1)) + 
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")
print(p)
}
data = original %>% filter(rho == 0.3, Tfull == 2, n_reps == 3000, beta_free == 2) %>% 
  filter(omega_var %in% c(0.25^2, 1, 4^2))

ggplot(data, aes(x = N, group = omega_var)) + 
  # First set: original p_beta
  geom_line(aes(y = p_rho_rescaled, color = as.factor(omega_var))) + 
  geom_point(aes(y = p_rho_rescaled, color = as.factor(omega_var)), size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_rho_rescaled, ymax = bino_ci_exact_upper_rho_rescaled, fill = as.factor(omega_var)), 
              alpha = 0.3, color = NA) +
  
  # First color & fill scale: blue-green
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Omega Parameter") +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Omega Parameter") +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_line(aes(y = p_beta_rescaled, color = as.factor(omega_var))) + 
  geom_point(aes(y = p_beta_rescaled, color = as.factor(omega_var)), shape = 1, size = 1.2) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_beta_rescaled, ymax = bino_ci_exact_upper_beta_rescaled, fill = as.factor(omega_var)), 
              alpha = 0.3, color = NA) +
  
  scale_color_manual(
    values = c("red", "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  scale_fill_manual(
    values = c("red", "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  guides(
    color = guide_legend( 
      title = paste0("Fixed\nParams:\nT = ", unique(data$Tfull),
                     "\nβ = ", unique(data$beta_free),
                     "\nρ = ", unique(data$rho),
                     "\n\n# of Reps:\n~3000\n\nMean β and\nRescaled:"
      )
    ),
    fill = guide_legend(
      title = paste0("Fixed\nParams:\nT = ", unique(data$Tfull),
                     "\nβ = ", unique(data$beta_free),
                     "\nρ = ", unique(data$rho),
                     "\n\n# of Reps:\n~3000\n\nMean β and\nRescaled:"
      )
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = expression("CI Comparison for Mean β1 Estimates vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability"
  ) +
  coord_cartesian(ylim = c(0.75, 1)) + 
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")
## Compare coverage rates between SE types while varying Omega (not rescaled) -------------------------------
data = original %>% filter(rho == -0.9, Tfull == 4, n_reps == 3000, beta_free == 4) %>%
  filter(beta_free %in% c(0,1,2,3,4), omega_var %in% c(4^2))

ggplot(data, aes(x = N, group = omega_var)) + 
  # First set: original p_beta
  geom_line(aes(y = p_rho, color = as.factor(omega_var))) + 
  geom_point(aes(y = p_rho, color = as.factor(omega_var)), size = 1) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_rho, ymax = bino_ci_exact_upper_rho, fill = as.factor(omega_var)), 
              alpha = 0.3, color = NA) +
  
  # First color & fill scale: blue-green
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Var ω² and\nModel-based") +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.7, name = "Var ω² and\nModel-based") +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  geom_line(aes(y = p_omega_rescaled, color = as.factor(omega_var))) + 
  geom_point(aes(y = p_omega_rescaled, color = as.factor(omega_var)), shape = 1, size = 1.2) + 
  geom_ribbon(aes(ymin = bino_ci_exact_lower_omega_rescaled, ymax = bino_ci_exact_upper_omega_rescaled, fill = as.factor(omega_var)), 
              alpha = 0.3, color = NA) +
  
  scale_color_manual(
    values = c("purple", "red",  "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  scale_fill_manual(
    values = c("purple", "red",  "orange", "gold", "darkred", "tomato", "darkorange"),
    name = "ω² and\nRescaled"
  ) +
  guides(
    color = guide_legend(
      title = paste0("Fixed\nParams:\nT = ", unique(data$Tfull),
                     "\nβ = ", unique(data$beta_free),
                     "\nρ = ", unique(data$rho),
                     "\n\n# of Reps:\n~3000\n\nVar ω² and\nRescaled:"
      )
    ),
    fill = guide_legend(
      title = paste0("Fixed\nParams:\nT = ", unique(data$Tfull),
                     "\nβ = ", unique(data$beta_free),
                     "\nρ = ", unique(data$rho),
                     "\n\n# of Reps:\n~3000\n\nVar ω² and\nRescaled:"
      )
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = expression("CI Comparison for Mean β1 Estimates vs. Sample Size (N)"),
    x = "Sample Size (N)",
    y = "Coverage Probability"
  ) +
  coord_cartesian(ylim = c(0.5, 1)) + 
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "black")






