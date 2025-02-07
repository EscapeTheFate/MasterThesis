#
# Master thesis - Code Space for Data Generating Process (DGP) of
# Ordered MNP-based Discrete Choice Models in Rprobit environment
#

# Install package from source ---------------------------------------------
#file.exists("~/Rprobit_0.3.2update.tar.gz")
#install.packages("~/Rprobit_0.3.2update.tar.gz", repos = NULL, type = "source")

# Load package
library(Rprobit)
library(tictoc)


# Data-Generating-Process Structure

check_if_stationary <- function(rho){
  if(!is.vector(rho)){
    rho = as.vector(rho)
  }
  dim = length(rho)
  if(dim == 1){
    if(abs(rho[1]) < 1){
      return(T)
    }
    else{
      return(F)
    }
  }
  if(dim == 2){
    if(abs(rho[2]) < 1 & (rho[1] + rho[2]) < 1 & (rho[2] - rho[1]) < 1){
      return(T)
    }
    else{
      return(F)
    }
  }
  else{
    return("Error: Only up to AR(2) is specified")
  }
}

create_grid <- function(...) {
  vars <- list(...)
  var_names <- sapply(substitute(list(...))[-1], deparse)
  setNames(expand.grid(vars), var_names)
}

draw_simulation_data <- function(Tfull, N, quest, rho, system){
  
  time = c(1:Tfull)
  Tp = length(time)*quest # number of questions per person 
  
  # Setup dataframe by first individual
  latent_class = sample(x = c(1,2,3), size = 1, prob = c(rep.int(1/3,3)))
  latent_means = c(2,4,7)
  X = matrix(0,Tp,2)
  X[,1] = rnorm(n = Tp, mean = latent_means[latent_class], sd = 1)
  X[,2] = rnorm(n = Tp, mean = 0, sd = 1) # b=1 is fixed 
  tau = system$tauk
  b = system$beta
  #beta = rnorm(n = Tp, mean = 0, sd = 0.5) # mixed part
  #epsilon = rnorm(n = 1, sd = system$Sigma)
  LG = t(chol(system$GammaT))
  LO = t(chol(system$Omega))
  lRE = dim(LO)[1]
  omega_sd = system$Omega
  gam = LO %*% rnorm(lRE, sd = omega_sd) # mixed part here
  
  # get_AR_epsilon <- function(time, rho, previous_error, Sigma_stationary = c(1), Sigma = c(1)){
  #   
  #   get_stationary_epsilon <- function(Sigma = c(1)){
  #     epsilon = rnorm(n=1, mean = 0, sd = Sigma)
  #     return(epsilon)
  #   }
  #   
  #   if(time == 1){
  #     epsilon_1 = get_stationary_epsilon(Sigma_stationary)
  #     return(epsilon_1)
  #   }
  #   else{
  #     u = get_stationary_epsilon(Sigma)
  #     epsilon_new = rho * previous_error + u
  #     return(epsilon_new)
  #   }
  # }
  
  # Obtain utilities for each choice occasion for individual 1
  U = matrix(0, nrow = Tp)
  # for (t in 1:Tp){
  #   if(t == 1){
  #     epsilon_old = get_AR_epsilon(time = t)
  #     U[t] = X[t,] %*% b + X[t,1] %*% beta[t] + epsilon_old
  #   }
  #   else{
  #     epsilon_new = get_AR_epsilon(time = t, rho = system$factor, previous_error = epsilon_old)
  #     U[t] = X[t,] %*% b + X[t,1] %*% beta[t] + epsilon_new
  #   }
  #   
  # }
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
    latent_means = c(3,5,8)
    X = matrix(0,Tp,2)
    (X[,1] = rnorm(n = Tp, mean = latent_means[latent_class], sd = 1)) # b=0.5 is free
    X[,2] = rnorm(n = Tp, mean = 0, sd = 1) # b=1 is fixed 
    tau = system$tauk
    b = system$beta
    #beta = rnorm(n = Tp, mean = 0, sd = 0.5) # mixed part
    #epsilon = rnorm(n = 1, sd = system$Sigma)
    LG = t(chol(system$GammaT))
    LO = t(chol(system$Omega))
    lRE = dim(LO)[1]
    gam = LO %*% rnorm(lRE, sd = omega_sd)
    
    # Obtain utilities for each choice occasion for individual 1
    U = matrix(0, nrow = Tp)
    # for (t in 1:Tp){
    #   if(t == 1){
    #     epsilon_old = get_AR_epsilon(time = t)
    #     U[t] = X[t,] %*% b + X[t,1] %*% beta[t] + epsilon_old
    #   }
    #   else{
    #     epsilon_new = get_AR_epsilon(time = t, rho = system$factor, previous_error = epsilon_old)
    #     U[t] = X[t,] %*% b + X[t,1] %*% beta[t] + epsilon_new
    #   }
    #   
    # }
    U = X %*% system$beta + X[,c(1:lRE)] %*% gam + LG %*% rnorm(dim(LG)[2])
    
    # Choose appropriate alternative by utility value
    y = U*0
    for (tj in 1:Tp){
      y[tj] = sum(U[tj]>tau)+1
    }
    df = rbind(df, data.frame(ID = j, X1 = X[,1], X2 = X[,2], y = y, time = time_df, quest = quest_df))
  } 
  
  return(df)
}

get_simulation_estimates <- function(Tfull = 5, # number of choice occasions (time)
                                     N = 100, # individuals
                                     quest = 1, beta_coefs = c(1, 1), rho = 0.5, return_val = "", timer = F, DontSkipFit = T){
  mod <- mod_AR_cl$new(alt  = 7,
                       Hb   = diag(2)[,-2,drop=FALSE],
                       fb   = matrix(0,2,1),
                       HO   = matrix(1,1,1),
                       fO   = matrix(0,1,1),
                       HL   = matrix(0,1,0),
                       fL   = matrix(1,1,1),
                       ordered = TRUE
  )
  
  mod$fb[2] = beta_coefs[2] # <- Since b = [b_1 (free), 2 (fixed)], we need to change fb entry
  mod$lag_length = 1 # The p of the error AR(p)-process
  mod$stationary <- TRUE # indicates whether the stationary distribution shall be used to start the state process
  
  # Setup for theta_0 true parameter value
  tau = c(0,rep.int(log(2),5))
  omega = 1
  theta_0 = c(beta_coefs[1], omega, tau, rho) # (beta (-coef1), Omega, Sigma (-all), tau, rho)
  
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
          #Rprob_conv = c(H_conv = T, J_conv = T)
          list(b, omega, rho, conv)
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
      
          message("Model fitting failed: ", conditionMessage(e))
      
          # Return the same structure with NA values
          list(b = b, omega = omega, rho = rho, conv = conv)
      }
    }
  )
  } else { # if skipped, return NA same as if error occured
    b <- c(b_coef = NA, b_sd = NA)
    omega <- c(omega_coef = NA, omega_sd = NA)
    rho <- c(rho_coef = NA, rho_sd = NA)
    conv <- c(nlm_code = 99)
    result = list(b = b, omega = omega, rho = rho, conv = conv)
  }
  if(timer == T){
    toc()
  }
  # Return results 
  return(result)
}

set.seed(1)
get_simulation_estimates(Tfull = 5, # number of choice occasions (time)
                         N = 100, # individuals
                         quest = 1, # num of questions (master thesis is for now limited to 1)
                         beta_coefs = c(1, 1), # beta coefs; here: (free, fixed)
                         rho = 0.5, # rho in error AR(1)-process
                         return_val = "summary", # if return should be something else
                         DontSkipFit = T) # If simulation 1 to K takes too long, 
                                          # one may stop at k and may want to continue 
                                          # estimating at k+1, then the DGP of iteration
                                          # 1 to k will be reproduced without estimation (time-consuming),
                                          # to then continue estimating afterwards (don't 
                                          # forget to pass on DontSkipFit = F for 1 to k)

# Obtain bias stats ------------------------------------------------------
# To-Do:
# - Add H/J-Mat conv code (nlm convergence code done)
# - Improve get_simultion_estimates to only depend on draw_simulation_data
# - Slot in stability check infront of parameter grid

get_bias_and_averages <- function(n_reps = 200, Tfull = 5, N = 100, quest = 1, beta_coefs = c(1,1), rho = 0.5,
                                  return_val = "", timer_bias = F, timer_data = F, DontSkipFit = T, saveEstimates = T){

  if(timer_bias){ # timer = c(inner_timer for get_bias)
    parameterMappedTo <- paste0("N=", N, ",rho=", rho, ",Tfull=", Tfull, ",n_reps=", n_reps)
    msg = paste0("Fitting of Rprobit Models: ", parameterMappedTo)
    tic(msg = msg)
  }
  sim_replications <- replicate(n_reps, get_simulation_estimates(Tfull, # number of choice occasions (time)
                                                                 N, # individuals
                                                                 quest, # num of questions
                                                                 beta_coefs, # beta coefs; here: (free, fixed)
                                                                 rho, # rho in error AR(1)-process
                                                                 return_val = return_val, # return coefs & sds
                                                                 timer = timer_data,
                                                                 DontSkipFit = DontSkipFit)) 
  if(timer_bias){
    toc()
  }
  # Overview of sim_replications:
  # - Rows = Parameter list in chronological order: (b, omega, rho, nlm_conv_code)
  # - Columns = n_reps
  
  # Calculate convergence success rate and filter out non-successful coefficient estimates
  conv_tab <- table(unlist(sim_replications[4,]))
  conv_success_index <- which(sim_replications[4,] == 1 | sim_replications[4,] == 2 | sim_replications[4,] == 3) # for whatever reason %in% does not work
  conv_failed_index <- which(sim_replications[4,] == 4 | sim_replications[4,] == 5 | sim_replications[4,] == 99) # code 99 = nlm crashed (no return value)
  if (is.integer(conv_failed_index) && length(conv_failed_index) == 0) {
    # All models successfully estimated
    msg <- paste0("All nlm models were successfully estimated (N=", N, ", Tfull=", Tfull, ", rho=", rho, ")")
  } else {
    # Some (or all) models failed
    if (!DontSkipFit) {
      msg <- paste0("Skipped nlm model iterations (", length(conv_failed_index), "/", n_reps, 
                    ") (N=", N, ", Tfull=", Tfull, ", rho=", rho, ")")
    } else {
      msg <- paste0("Some nlm models (", length(conv_failed_index), "/", n_reps, 
                    ") failed to converge (N=", N, ", Tfull=", Tfull, ", rho=", rho, ")")
    }
    # Adjust sim_replications to exclude failed indices
    sim_replications <- sim_replications[, -conv_failed_index]
  }
  print(msg)
  cat("\n")
  
  
  # Calculate bias:
  if(length(conv_success_index) == 0){
    beta_bias = beta_sd_avg = omega_bias = omega_sd_avg = rho_bias = rho_sd_avg = NA
    
    return(list(beta_bias = beta_bias, 
                beta_sd_avg = beta_sd_avg,
                omega_bias = omega_bias,
                omega_sd_avg = omega_sd_avg,
                rho_bias = rho_bias,
                rho_sd_avg = rho_sd_avg))
  }
  
  true_parameters = c(beta = beta_coefs[1], omega = 1, rho = rho)
  
  # Get beta bias
  beta_estimates <- sapply(sim_replications[1,], function(x) x["b_coef"])
  beta_bias <- mean(beta_estimates - true_parameters["beta"])
  
  # Repeat for omega and rho
  omega_bias <- mean((omega_estimates <- sapply(sim_replications[2,], function(x) x["omega_coef"])) - true_parameters["omega"])
  rho_bias <- mean((rho_estimates <- sapply(sim_replications[3,], function(x) x["rho_coef"])) - true_parameters["rho"])
  
  # Now, calculate average standard deviation for each estimate
  # Get average sds:
  beta_sds = sapply(sim_replications[1,], function(x) x["b_sd"])
  beta_sd_avg = mean(beta_sds)
  
  # Repeat for omega and rho
  omega_sd_avg = mean(omega_sds <- sapply(sim_replications[2,], function(x) x["omega_sd"]))
  rho_sd_avg = mean(rho_sds <- sapply(sim_replications[3,], function(x) x["rho_sd"]))
  
  if(saveEstimates == T){
    filename <- paste0("BiasAvgSD_", "N=", N, "rho=", rho, "Tfull=", Tfull, "nreps=", n_reps, ".csv")
    df = data.frame(beta = beta_estimates, 
                    beta_sd = beta_sds,
                    omega = omega_estimates,
                    omega_sd = omega_sds,
                    rho = rho_estimates,
                    rho_sd = rho_sds)
    if(file.exists(paste0(getwd(),"/GitHub/MasterThesis/BiasAvgSD_Folder"))){
      old_path = getwd()
      setwd(dir = paste0(getwd(),"/GitHub/MasterThesis/BiasAvgSD_Folder"))
      write.csv(df, file = filename, row.names = F)
      msg = paste0("Estimates saved: ", filename, " @ ", getwd()) 
      print(msg)
      setwd(dir = old_path)
    } else {
      write.csv(df, file = filename, row.names = F)
      msg = paste0("Estimates saved: ", filename, " @: ", getwd()) 
      print(msg)
    }
    cat("\n")
  }
  
  return(list(beta_bias = beta_bias, 
              beta_sd_avg = beta_sd_avg,
              omega_bias = omega_bias,
              omega_sd_avg = omega_sd_avg,
              rho_bias = rho_bias,
              rho_sd_avg = rho_sd_avg))
}

#foo = get_bias_and_averages(n_reps = 20, Tfull = 5, N = 100, quest = 1, beta_coefs = c(1,1), rho = 0.9, timer_bias = F, timer_data = F, DontSkipFit = F, saveEstimates = F)
#foo$beta_bias

## Extract bias and averages for different models --------------------------

# Define grid level to map get_bias on:
#rho = c(0.9, 0.5, 0.25, 0, -0.25, -0.5, -0.9)
rho = 0.5
Tfull = c(5)
N = c(5, 10, 20, 30, 40, 50, 100, 200, 500)
n_reps = c(20)
timer_bias = T

parameters = create_grid(N, n_reps, Tfull, rho, timer_bias)
(parameters)

# Fit models to grid
bias <- Map(get_bias_and_averages, N = parameters$N, n_reps = parameters$n_reps, Tfull = parameters$Tfull, rho = parameters$rho, timer_bias = parameters$timer_bias)

# Formating results
bias <- do.call(rbind, bias)
sim_results <- cbind(parameters, bias)
str(sim_results)
sim_results$beta_bias = unlist(sim_results$beta_bias)
sim_results$omega_bias = unlist(sim_results$omega_bias)
sim_results$rho_bias = unlist(sim_results$rho_bias)
sim_results

# Save calculated dataframe for later
write.csv(sim_results, file = "BiasResultsDifferentRho.csv", row.names = F)
#sim_results = read.csv(file = "~/GitHub/MasterThesis/InitialBIASresults.csv")


library(ggplot2)
library(tidyr)

# Reshape the dataframe using tidyr::pivot_longer
sim_results_long <- sim_results %>%
  pivot_longer(
    cols = c(beta_bias, omega_bias, rho_bias), # Columns to reshape
    names_to = "Coefficient",                 # New column name for variable names
    values_to = "Bias"                        # New column name for values
  )

# simple plot
i = 6
tx = paste0("Beta Bias vs. Sample Size (N) for rho =", sim_results[(i*9+1):((i+1)*9),]$rho[1])
# Beta plot
plot(x = sim_results[(i*9+1):((i+1)*9),]$N, y = sim_results[(i*9+1):((i+1)*9),]$beta_bias, type = "b", 
     col = "blue",                    
     pch = 19,                        
     lwd = 2,                         
     xlab = "Number of Individuals (N)", 
     ylab = "Beta Bias",              
     main = tx, 
     cex.axis = 0.9,                  
     cex.lab = 1,                   
     cex.main = 1.2)
grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")
abline(h = 0, col = "black", lty = "dotted", lwd = 2)
# Omega plot
tx = paste0("Omega Bias vs. Sample Size (N) for rho =", sim_results[(i*9+1):((i+1)*9),]$rho[1])
plot(x = sim_results[(i*9+1):((i+1)*9),]$N, y = sim_results[(i*9+1):((i+1)*9),]$omega_bias, type = "b", 
     col = "red",                    
     pch = 19,                       
     lwd = 2,                         
     xlab = "Number of Individuals (N)", 
     ylab = "Omega Bias",             
     main = tx, 
     cex.axis = 0.9,                 
     cex.lab = 1,                     
     cex.main = 1.2)
grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")
abline(h = 0, col = "black", lty = "dotted", lwd = 2)
# Rho plot
tx = paste0("Rho Bias vs. Sample Size (N) for rho =", sim_results[(i*9+1):((i+1)*9),]$rho[1])
plot(x = sim_results[(i*9+1):((i+1)*9),]$N, y = sim_results[(i*9+1):((i+1)*9),]$rho_bias, type = "b", 
     col = "darkgreen",                    # Line color
     pch = 19,                        # Point symbol
     lwd = 2,                         # Line width
     xlab = "Number of Individuals (N)", # X-axis label
     ylab = "Rho Bias",              # Y-axis label
     main = tx, # Title
     cex.axis = 0.9,                  # Axis text size
     cex.lab = 1,                     # Label text size
     cex.main = 1.2)
grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")
abline(h = 0, col = "black", lty = "dotted", lwd = 2)

#
i = 0
# Beta plot
plot(x = sim_results[(i*9+1):((i+1)*9),]$N, y = sim_results[(i*9+1):((i+1)*9),]$beta_bias, type = "b", 
     col = "blue",                    
     pch = 19,                        
     lwd = 2,                         
     xlab = "Number of Individuals (N)", 
     ylab = "Beta Bias",
     ylim = c(-0.1, 0.5),
     main = "Beta Bias vs. Sample Size (N) for varying rho", 
     cex.axis = 0.9,                  
     cex.lab = 1,                   
     cex.main = 1.2)
grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")
abline(h = 0, col = "black", lty = "dotted", lwd = 2)
colors = c("lightblue", "grey", "green", "yellow", "orange", "red")
for (i in 1:length(colors)){
  points(x = sim_results[(i*9+1):((i+1)*9),]$N, y = sim_results[(i*9+1):((i+1)*9),]$beta_bias,
         type = "b", pch = 19, col = colors[i], lwd = 2)
}



# Create the ggplot
ggplot(sim_results_long, aes(x = N, y = Bias, color = Coefficient, group = Coefficient)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~ Coefficient, scales = "free_y") +  # Separate plot for each coefficient
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Bias of Coefficient Estimates Across N",
    x = "Number of Individuals (N)",
    y = "Bias",
    color = "Coefficient") #+
  # annotate("text", x = max(sim_results$N) * 0.7, y = max(sim_results_long$Bias) * 0.8, 
  #          label = paste0("Tfull = ", unique(sim_results$Tfull), 
  #                         "\nrho = ", unique(sim_results$rho), 
  #                         "\nn_reps = ", unique(sim_results$n_reps)),
  #          hjust = 0, color = "black", size = 5, alpha = 0.7)
foo <- sim_results[,1:5] %>% group_by(rho)
ggplot(sim_results_long %>% filter(Coefficient == "beta_bias"), aes(x = N, y = Bias, color = rho)) +
  geom_line(size = 1.2)

# Print the plot
print(bias_plot)

# Create the ggplot
#ggplot(sim_results, aes(x = N, y = beta_bias)) + geom_line(size = 1)
ggplot(CI_results_long, aes(x = N, y = rate*100, color = parameter)) +  # rate * 100 for %
  geom_line(size = 1) + geom_point(size = 2) + # line plot with points
  geom_hline(yintercept = 95, linetype = "dotted", color = "black", size = 0.8) + 
  labs(title = "Percentage of Estimates within 95%-Confidence Interval",
       x = "Number of Individuals (N)",
       y = "Percentage within CI (%)",
       color = "Parameter") +
  scale_color_manual(values = c("beta_pc" = "blue", "omega_pc" = "red", "rho_pc" = "green")) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "top") + 
  annotate("text",
           x = max(CI_results$N) * 0.7, y = 0.3,  # Adjust position
           label = paste0(
             "rho: ", CI_results$rho[1], "\n",
             "Tfull: ", CI_results$Tfull[1], "\n",
             "n_reps: ", CI_results$n_reps[1]
           ),
           hjust = 0, vjust = 0, size = 10, color = "gray50", fontface = "italic"
  )




# Obtain Variance resembling True Variance (close enough for large N)------
rho = seq(from = -0.9, to = 0.9, by = 0.1)
rho = 0.9
Tfull = 5
N = 20000 # maybe 20000?
n_reps = 1
timer_bias = T
saveEstimates = F

parameters = create_grid(N, n_reps, Tfull, rho, timer_bias)
(parameters)

set.seed(1)
BiasAvgSD <- Map(get_bias_and_averages, N = parameters$N, n_reps = parameters$n_reps, Tfull = parameters$Tfull, rho = parameters$rho, timer_bias = parameters$timer_bias)

# Extract the SD entries from the list
beta_sd_avg <- sapply(BiasAvgSD, function(x) x$beta_sd_avg)
omega_sd_avg <- sapply(BiasAvgSD, function(x) x$omega_sd_avg)
rho_sd_avg <- sapply(BiasAvgSD, function(x) x$rho_sd_avg)

# Combine them into data frame
sd_dataframe <- data.frame(
  rho = parameters$rho,
  beta_sd_avg = beta_sd_avg,
  omega_sd_avg = omega_sd_avg,
  rho_sd_avg = rho_sd_avg
)

sd_dataframe
old_path = getwd()
setwd(dir = paste0(getwd(), "/GitHub/MasterThesis/TrueVarianceInfoCollection"))
write.csv(sd_dataframe, file = "CIresults.csv", row.names = F)
setwd(dir = old_path)


# Obtain Confidence Interval stats ----------------------------------------
get_CI_results <- function(n_reps = 200, Tfull = 5, N = 100, quest = 1, beta_coefs = c(1,1), rho = 0.5, alpha = 0.05, return_val = ""){
  
  sim_replications <- replicate(n_reps, get_simulation_estimates(Tfull, # number of choice occasions (time)
                                                                 N, # individuals
                                                                 quest, # num of questions
                                                                 beta_coefs, # beta coefs; here: (free, fixed)
                                                                 rho, # rho in error AR(1)-process
                                                                 return_val = "")) # return coefs & sds
                          
  # Overview of sim_replications:
  # - Rows = Parameter list in chronological order: (b, omega, rho, nlm_conv_code)
  # - Columns = n_reps
  
  # Calculate convergence success rate and filter out non-successful coefficient estimates
  conv_tab <- table(unlist(sim_replications[4,]))
  conv_success_index <- which(sim_replications[4,] == 1 | sim_replications[4,] == 2 | sim_replications[4,] == 3) # for whatever reason %in% does not work
  conv_failed_index <- which(sim_replications[4,] == 4 | sim_replications[4,] == 5)
  if(is.integer(conv_failed_index) && length(conv_failed_index) == 0){
    msg = paste0("All nlm models were successfully estimated (N=", N, ", Tfull=", Tfull, ", rho=", rho, ")")
    print(msg)
  } else {
    msg = paste0("Some nlm models (",length(conv_failed_index), "out of", n_reps,") failed to converge (N=", N, ", Tfull=", Tfull, ", rho=", rho, ")")
    print(msg)
    sim_replications = sim_replications[,-conv_failed_index]
  }
  
  # Get estimates
  beta_estimates <- sapply(sim_replications[1,], function(x) x["b_coef"])
  omega_estimates <- sapply(sim_replications[2,], function(x) x["omega_coef"])
  rho_estimates <- sapply(sim_replications[3,], function(x) x["rho_coef"])
  
  
  # Calculate the proportion of model estimates that lie within (1-alpha)%-CI
  
  # First, load in correct true variance file that was calculated on N = 10000
  N_to_load <- 10000
  rho_to_load <- rho
  Tfull_to_load <- Tfull
  
  # construct appropriate filename and routing
  github_folder <- paste0(getwd(), "/GitHub/MasterThesis/TrueVarianceInfoCollection")
  filename_to_load <- paste0("trueVarN=", N_to_load, "rho=", rho_to_load, "Tfull=", Tfull_to_load, ".csv")
  full_path_to_load <- file.path(github_folder, filename_to_load)
  
  # Load the file
  if (file.exists(full_path_to_load)) {
    true_var <- read.csv(file = full_path_to_load)
    names(true_var) = "true_var"
    true_var <- as.vector(unlist(true_var))
  } else {
    print(paste("File not found:", full_path_to_load))
  }
  
  
  # if(file.exists("~/trueVarRho0.5Tfull5.csv")){
  #   true_var <- load("trueVarRho0.5Tfull5.csv")
  # } else {
  #   true_var <- read.csv("~/GitHub/MasterThesis/trueVarRho0.5Tfull5.csv")
  #   names(true_var) = "true_var"
  #   true_var <- as.vector(unlist(true_var))
  # }
  
  adjust_var_to_N <- function(unadjusted_true_var, N_new, N_old){
    return((N_old*unadjusted_true_var)/N_new)
  }
  
  true_var <- true_var[c(1,2,9)] # beta(1), omega(2), tau(3:8), rho(9)
  adjusted_variances <- adjust_var_to_N(true_var, N_new = N, N_old = N_to_load)
  true_parameters = c(beta = beta_coefs[1], omega = 1, rho = rho)
  
  # Calculate (1-alpha)%-CI boundaries
  if(alpha > 0 & alpha < 1){
    z_value = qnorm(1-alpha/2) # For different CIs
  } else {
    print("Error: Alpha value invalid")
  }
  
  ci_lower <- true_parameters - z_value * sqrt(adjusted_variances)
  ci_upper <- true_parameters + z_value * sqrt(adjusted_variances)
  
  # Check cuccess or failure for each parameter
  success_beta <- sum(beta_estimates >= ci_lower[1] & beta_estimates <= ci_upper[1])
  failure_beta <- length(beta_estimates) - success_beta
  
  success_omega <- sum(omega_estimates >= ci_lower[2] & omega_estimates <= ci_upper[2])
  failure_omega <- length(omega_estimates) - success_omega
  
  success_rho <- sum(rho_estimates >= ci_lower[3] & rho_estimates <= ci_upper[3])
  failure_rho <- length(rho_estimates) - success_rho
  
  
  if(nzchar(return_val)){
    return(Success_rate = c(success_beta, success_omega, success_rho)/length(beta_estimates))
  } else {
  # Display results in table
  results_table <- data.frame(
    Parameter = c("beta", "omega", "rho"),
    Total_Estimates = c(rep.int(length(beta_estimates), 3)),
    Success_Count = c(success_beta, success_omega, success_rho),
    Success_rate = c(success_beta, success_omega, success_rho)/length(beta_estimates),
    CI_Lower = ci_lower,
    True_param = true_parameters,
    CI_Upper = ci_upper
  )
  
  # Output
  print(paste0("Confidence Level: ", (1-alpha)*100, "%"))
  print(results_table)
  }
 
}

get_CI_results(n_reps = 100,
               Tfull = 5,
               N = 50,
               quest = 1,
               beta_coefs = c(1,1),
               rho = 0.5,
               alpha = 0.05,
               return_val = "success_rate")

create_grid <- function(...) {
  vars <- list(...)
  var_names <- sapply(substitute(list(...))[-1], deparse)
  setNames(expand.grid(vars), var_names)
}

# ----- Calculate first set (hold everything but N)
N = c(5, 10, 20, 30, 40, 50, 100, 200, 500, 750, 1000)
Tfull = c(5)
n_reps = c(200)
rho = c(0.5)

parameters = create_grid(N, n_reps, Tfull, rho)
suppressWarnings({
success_rates <- Map(get_CI_results, N = parameters$N, n_reps = parameters$n_reps, Tfull = parameters$Tfull, rho = parameters$rho, return_val = "success_rate")
})

# Formating results
success_rates <- do.call(rbind, success_rates)

CI_results <- cbind(parameters, success_rates)
names(CI_results)[5:7] = c("beta_pc", "omega_pc", "rho_pc")
CI_results

library(tidyr)
library(ggplot2)

# Reshape to long format
CI_results_long <- CI_results %>%
  pivot_longer(
    cols = c(beta_pc, omega_pc, rho_pc), # Specify columns to reshape
    names_to = "parameter",              # New column for parameter names
    values_to = "rate"                   # New column for percentage rates
  )

# Save calculated dataframe for later
#write.csv(CI_results, file = "InitialCIresults.csv", row.names = F)
CI_results = read.csv(file = "~/GitHub/MasterThesis/InitialCIresults.csv")

# Plot using ggplot2
ggplot(CI_results_long, aes(x = N, y = rate*100, color = parameter)) +  # rate * 100 for %
  geom_line(size = 1) + geom_point(size = 2) + # line plot with points
  geom_hline(yintercept = 95, linetype = "dotted", color = "black", size = 0.8) + 
  labs(title = "Percentage of Estimates within 95%-Confidence Interval",
       x = "Number of Individuals (N)",
       y = "Percentage within CI (%)",
       color = "Parameter") +
  scale_color_manual(values = c("beta_pc" = "blue", "omega_pc" = "red", "rho_pc" = "green")) +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "top") + 
  annotate("text",
           x = max(CI_results$N) * 0.7, y = 0.3,  # Adjust position
           label = paste0(
             "rho: ", CI_results$rho[1], "\n",
             "Tfull: ", CI_results$Tfull[1], "\n",
             "n_reps: ", CI_results$n_reps[1]
           ),
           hjust = 0, vjust = 0, size = 10, color = "gray50", fontface = "italic"
  )

