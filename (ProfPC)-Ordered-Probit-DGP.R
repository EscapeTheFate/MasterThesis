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
library(purrr)
library(future.apply)

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
  X[,1] = rnorm(n = Tp, mean = latent_means[latent_class], sd = 1) # b_1 is free
  X[,2] = rnorm(n = Tp, mean = 0, sd = 1) # b_2 is fixed 
  tau = system$tauk
  b = system$beta
  
  LG = t(chol(system$GammaT))
  LO = t(chol(system$Omega))
  lRE = dim(LO)[1]
  
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
    tau = system$tauk
    b = system$beta
    
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
      y[tj] = sum(U[tj]>tau)+1
    }
    df = rbind(df, data.frame(ID = j, X1 = X[,1], X2 = X[,2], y = y, time = time_df, quest = quest_df))
  } 
  
  return(df)
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

# set.seed(1) ----
# get_simulation_estimates(Tfull = 5, # number of choice occasions (time)
#                          N = 100, # individuals
#                          quest = 1, # num of questions (master thesis is for now limited to 1)
#                          beta_free = 1, # beta coefs; here: (free, fixed)
#                          rho = 0.5, # rho in error AR(1)-process
#                          omega_var = 0.5,
#                          return_val = "summary", # if return should be something else
#                          DontSkipFit = T, # If simulation 1 to K takes too long, 
#                          timer = F)       # one may stop at k and may want to continue 
#                                           # estimating at k+1, then the DGP of iteration
#                                           # 1 to k will be reproduced without estimation (time-consuming),
#                                           # to then continue estimating afterwards (don't 
#                                           # forget to pass on DontSkipFit = F for 1 to k)

# Obtain error stats ------------------------------------------------------
# To-Do:
# - Slot in stability check infront of parameter grid

get_error_and_averages <- function(n_reps, Tfull, N, quest = 1, beta_free, rho, omega_var,
                                  return_val = "", timer_error = T, timer_data = F, DontSkipFit = T, saveEstimates = T){

  if(timer_error){ 
    parameterMappedTo <- paste0("N=", N, ", Tfull=", Tfull, ", beta=", beta_free, ", rho=", rho, ", omega=", omega_var)
    msg = paste0("Fitting of Rprobit Models: ", parameterMappedTo)
    tic(msg = msg)
  }
  
  plan(multisession)
  sim_replications <- future_replicate(n_reps, get_simulation_estimates(
    Tfull, N, quest, beta_free, rho, omega_var, return_val, timer_data, DontSkipFit),
    simplify = FALSE) 
  
  if(timer_error){
    # toc()
    foo = toc(quiet = TRUE)
    elapsed_time = foo$toc - foo$tic
    elapsed_minutes = round(elapsed_time/60, 2)
    cat(paste0("Elapsed fitting time (n_reps=", n_reps, "): ", elapsed_minutes, " min, Parameters: ", parameterMappedTo, " , Finished: [", Sys.time(), "]\n"))
  }
  # Overview of sim_replications:
  # - Rows = Parameter list in chronological order: (b, omega, rho, nlm_conv_code)
  # - Columns = n_reps
  
  # Calculate convergence success rate and filter out non-successful coefficient estimates
  conv_status <- sapply(sim_replications, function(x) x[[4]]["nlm_code"])
  conv_success_index <- conv_status %in% c(1, 2, 3)
  conv_failed_index <- conv_status %in% c(4, 5, 99) # code 99 for nlm crash
  
  if(sum(conv_failed_index) == 0 & DontSkipFit == T){
    msg <- paste0("All nlm models (n_reps=", n_reps, ") were successfully estimated (N=", N, 
                  ", Tfull=", Tfull, ", rho=", rho, ", beta=", beta_free, ", omega=", omega_var, ")")
  } else {
    if (!DontSkipFit) {
      msg <- paste0("Skipped nlm model iterations (", n_reps, ") (N=", N, ", Tfull=", Tfull, 
                    ", rho=", rho, ", beta=", beta_free, ", omega=", omega_var, ")")
    } else {
      msg <- paste0("Some nlm models (", sum(conv_failed_index), "/", n_reps, 
                    ") failed to converge (N=", N, ", Tfull=", Tfull, ", rho=", rho, ", beta=", beta_free, 
                    ", omega=", omega_var, ")")
    }
    
    # # Remove failed iterations
    # sim_replications <- sim_replications[conv_success_index] # ---> New version: Use na.rm in mean calc, so that conv info is in estimate df collection
  }
  print(msg)
  cat("\n")
  
  if(sum(conv_success_index) == 0){ # only if estimating is skipped, or if everything failed
    return(list(beta_error = NA, beta_sd_avg = NA, omega_error = NA, omega_sd_avg = NA,
                rho_error = NA, rho_sd_avg = NA))
  }
  
  ## Calculate error
  
  # Get estimates
  beta_estimates <- sapply(sim_replications, function(x) x[[1]]["b_coef"])
  omega_estimates <- sapply(sim_replications, function(x) x[[2]]["omega_coef"]) 
  rho_estimates <- sapply(sim_replications, function(x) x[[3]]["rho_coef"])
  
  true_parameters <- c(beta = beta_free, omega = omega_var, rho = rho)
  
  # Get mean error
  beta_error <- mean(beta_estimates - true_parameters["beta"], na.rm = T)
  omega_error <- mean(omega_estimates - true_parameters["omega"], na.rm = T) # we are comparing omega_var entry with estimated omega_var parameter
  rho_error <- mean(rho_estimates - true_parameters["rho"], na.rm = T)
  
  # Now, calculate average standard deviation for each estimate
  beta_sd_avg <- mean(beta_sds <- sapply(sim_replications, function(x) x[[1]]["b_sd"]), na.rm = T)
  omega_sd_avg <- mean(omega_sds <- sapply(sim_replications, function(x) x[[2]]["omega_sd"]), na.rm = T) # these are sd's of omega_var
  rho_sd_avg <- mean(rho_sds <- sapply(sim_replications, function(x) x[[3]]["rho_sd"]), na.rm = T)
  
  if(saveEstimates){
    filename <- paste0("ErrAvgSD_N=", N, "_Tfull=", Tfull, "_nreps=", n_reps, 
                       "_rho=", rho, "_beta=", beta_free, "_omega=", omega_var, ".csv")
    
    # Create df of model estimates and save for later CI calcs/if something crashed
    df = data.frame(beta = beta_estimates, beta_sd = beta_sds,
                    omega = omega_estimates, omega_sd = omega_sds,
                    rho = rho_estimates, rho_sd = rho_sds,
                    conv_code = conv_status)
    
    for (i in 1:6){
      # Extract (raw) tau_i + it's (raw) sd in iteration i
      correct_tau = paste0("tau_coef", i)
      correct_tau_sd = paste0("tau_sd", i)
      
      tau_i_estimates <- sapply(sim_replications, function(x) x[[5]][correct_tau])
      tau_i_sd <- sapply(sim_replications, function(x) x[[5]][correct_tau_sd])
      
      old_row_count = dim(df)[2]
      df = cbind(df, tau_i_estimates, tau_i_sd)
      names(df)[(old_row_count+1):(old_row_count+2)] = c(paste0("tau", i), paste0("tau", i, "_sd"))
    }
    
    if(N == 10000 | N == 20000){
      if(file.exists(paste0(getwd(),"/GitHub/MasterThesis/True_Variance_Estimate_Collection"))){
        old_path = getwd()
        setwd(dir = paste0(getwd(),"/GitHub/MasterThesis/True_Variance_Estimate_Collection"))
        write.csv(df, file = filename, row.names = F)
        msg = paste0("Estimates saved: ", filename, " @ ", getwd()) 
        print(msg)
        setwd(dir = old_path)
      } else {
        write.csv(df, file = filename, row.names = F)
        msg = paste0("Estimates saved: ", filename, " @: ", getwd()) 
        print(msg)
      }
      
    } else {
      if(file.exists(paste0(getwd(),"/GitHub/MasterThesis/Estimate_Collection"))){
        old_path = getwd()
        setwd(dir = paste0(getwd(),"/GitHub/MasterThesis/Estimate_Collection"))
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
  }
  
  return(list(beta_error = beta_error, 
              beta_sd_avg = beta_sd_avg,
              omega_error = omega_error,
              omega_sd_avg = omega_sd_avg,
              rho_error = rho_error,
              rho_sd_avg = rho_sd_avg))
}

#set.seed(1) # ---------
#foo = get_error_and_averages(n_reps = 1, Tfull = 2, N = 20000, quest = 1, beta_free = 1, rho = 0.9, omega_var = 0.2^2, timer_error = T, timer_data = F, DontSkipFit = T, saveEstimates = F)

## Extract error and average sd for different models --------------------------

# Define grid level to map get_error on (before on all pcs):

## On Home PC --------------------------------------------------------------
rho = seq(from = -0.9, to = 0.9, by = 0.1)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = 1
omega_var = 1
Tfull = c(4) 
set.seed(5) 

rho = seq(from = -0.9, to = 0.9, by = 0.1)
N = 10000
n_reps = 5
beta_free = 1
omega_var = 1
Tfull = c(2:4)
set.seed(6) # set.seed(9) for row 45:57

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = 1
omega_var = c(0.1^2, 0.2^2)
Tfull = c(2:4) 
set.seed(8) # set.seed(9) for N = 125, 150, 200, 300, 500, 1000

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = 1
omega_var = c(16, 25)
Tfull = c(2:4) 
set.seed(12) # set.seed(14) starting from ro 34: parameters = parameters[34:nrow(parameters),]
# set.seed(15) starting from row 63, set.seed(16) starting from row 92

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(0.5, 0.6)
omega_var = 1
Tfull = c(2:4) # c(2:5)
set.seed(17) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(0, 1.1)
omega_var = 1
Tfull = c(2:4) # c(2:5)
set.seed(24) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(1.2, 1.3)
omega_var = 1
Tfull = c(2:4)
set.seed(30) 

## On work PC --------------------------------------------------------------
rho = seq(from = -0.9, to = 0.9, by = 0.1)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(1000)
beta_free = 1
omega_var = 1
Tfull = c(2) # c(2:5)
set.seed(1) 

rho = seq(from = -0.9, to = 0.9, by = 0.1)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = 1
omega_var = 1
Tfull = c(2) # c(2:5)
set.seed(2) 

rho = seq(from = -0.9, to = 0.9, by = 0.1)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = 1
omega_var = 1
Tfull = c(3) # c(2:5)
set.seed(3) 

rho = seq(from = -0.9, to = 0.9, by = 0.1)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(1000)
beta_free = 1
omega_var = 1
Tfull = c(3:4) # c(2:5)
set.seed(4) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = 0.9
omega_var = 1
Tfull = c(2:4) # c(2:5)
set.seed(11) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = 0.8
omega_var = 1
Tfull = c(2:4) # c(2:5)
set.seed(18) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = 0.3
omega_var = 1
Tfull = c(2:4) # c(2:5)
set.seed(21) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = 0.1
omega_var = 1
Tfull = c(2:4) # c(2:5)
set.seed(23) 


## On Prof PC --------------------------------------------------------------
rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) 
n_reps = c(500)
beta_free = 0.7
omega_var = 1
Tfull = c(2:4) 
set.seed(19)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) 
n_reps = c(500)
beta_free = 0.4
omega_var = 1
Tfull = c(2:4) 
set.seed(20)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) 
n_reps = c(500)
beta_free = 0.2
omega_var = 1
Tfull = c(2:4) 
set.seed(22)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(2.1, 2.2, 2.3)
omega_var = 1
Tfull = c(2:4) # c(2:5)
set.seed(39) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(0.5, 0.6, 0.7, 0.8, 0.9)
omega_var = 0.25^2
Tfull = c(2:4) # c(2:5)
set.seed(44) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(1.7, 1.8, 1.9, 2)
omega_var = 0.25^2
Tfull = c(2:4) # c(2:5)
set.seed(48) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) # 20:60 missing due to crashes, adjust later
n_reps = c(500)
beta_free = c(1.2, 1.3, 1.4, 1.5, 1.6)
omega_var = 4^2
Tfull = c(2:4)
set.seed(57) 
(parameters = create_grid(rho, n_reps, beta_free, omega_var, Tfull, N))
error <- pmap(parameters, get_error_and_averages)

# rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
# N = c(50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) # 20:40 missing due to crashes, adjust later
# n_reps = c(500)
# beta_free = c(2.3, 2.4, 2.5)
# omega_var = 4^2
# Tfull = c(2:4)
# set.seed(62) # set.seed(63) starting from row Y (=ind+1), adjust later
# (parameters = create_grid(rho, n_reps, beta_free, omega_var, Tfull, N))
# (ind = which(parameters$N == 90 & parameters$Tfull == 3  & parameters$rho == 0.6 & parameters$beta == 2.5 & parameters$omega == 16))
# parameters = parameters[(ind+1):nrow(parameters),]
# set.seed(63)
# error <- pmap(parameters, get_error_and_averages)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) # 20:40 missing due to crashes, adjust later
n_reps = c(500)
beta_free = c(2.3, 2.4, 2.5)
omega_var = 4^2
Tfull = c(2:4)
set.seed(63) # <---- adjusted seed due to error

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) # 20:40 missing due to crashes, adjust later
n_reps = c(500)
beta_free = c(2.8, 2.9, 3)
omega_var = 4^2
Tfull = c(2:4)
set.seed(72)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) # 20:40 missing due to crashes, adjust later
n_reps = c(2500)
beta_free = c(1)
omega_var = c(0.25^2)
Tfull = c(2:4)
set.seed(78)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(500) # 20:300 missing due to calcs on other PC, adjust later
n_reps = c(2500)
beta_free = c(1)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(83)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) 
n_reps = c(500)
beta_free = c(1)
omega_var = c(0.25^2)
Tfull = c(2:4)
set.seed(184) # set.seed(186) starting from row 5
parameters = parameters[-c(1:4),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(300, 500, 1000)
n_reps = c(500)
beta_free = c(3)
omega_var = c(0.25^2)
Tfull = c(2:4)
set.seed(197)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(100, 125, 150, 200, 300, 500, 1000) 
n_reps = c(500)
beta_free = c(3.1)
omega_var = c(0.25^2, 1, 4^2)
Tfull = c(2:4)
set.seed(200)

# Fit models to grid
(parameters = create_grid(rho, n_reps, beta_free, omega_var, Tfull, N))

error <- pmap(parameters, get_error_and_averages)
(error_results = map_dfr(error, as.data.frame))
(error_results = cbind(error_results, parameters))

# Save calculated dataframe for later
file.exists(paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_laptop.csv"))
previous_error_results = read.csv(file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"))
error_results = dplyr::bind_rows(previous_error_results, error_results)
write.csv(error_results, file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_large_Nreps_omega_var16.csv"), row.names = F)


# Calculate true variance on large sample ---------------------------------
rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20000) 
n_reps = c(1)
beta_free = c(1, 1.1, 1.2, 1.3, 1.4, 1.5)
omega_var = 0.0625
Tfull = c(2:4)
set.seed(32) # set.seed(34) starting from row 61

(parameters = create_grid(rho, n_reps, beta_free, omega_var, Tfull, N))

# Fit models to grid
error <- pmap(parameters, get_error_and_averages)
(error_results = map_dfr(error, as.data.frame))
(error_results = cbind(error_results, parameters))

# -------------------------------------------------------------------------


library(ggplot2)
library(tidyr)

# Reshape the dataframe using tidyr::pivot_longer
sim_results_long <- sim_results %>%
  pivot_longer(
    cols = c(beta_error, omega_error, rho_error), # Columns to reshape
    names_to = "Coefficient",                 # New column name for variable names
    values_to = "Error"                        # New column name for values
  )

# simple plot
i = 6
tx = paste0("Beta Error vs. Sample Size (N) for rho =", sim_results[(i*9+1):((i+1)*9),]$rho[1])
# Beta plot
plot(x = sim_results[(i*9+1):((i+1)*9),]$N, y = sim_results[(i*9+1):((i+1)*9),]$beta_error, type = "b", 
     col = "blue",                    
     pch = 19,                        
     lwd = 2,                         
     xlab = "Number of Individuals (N)", 
     ylab = "Beta Error",              
     main = tx, 
     cex.axis = 0.9,                  
     cex.lab = 1,                   
     cex.main = 1.2)
grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")
abline(h = 0, col = "black", lty = "dotted", lwd = 2)
# Omega plot
tx = paste0("Omega Error vs. Sample Size (N) for rho =", sim_results[(i*9+1):((i+1)*9),]$rho[1])
plot(x = sim_results[(i*9+1):((i+1)*9),]$N, y = sim_results[(i*9+1):((i+1)*9),]$omega_error, type = "b", 
     col = "red",                    
     pch = 19,                       
     lwd = 2,                         
     xlab = "Number of Individuals (N)", 
     ylab = "Omega Error",             
     main = tx, 
     cex.axis = 0.9,                 
     cex.lab = 1,                     
     cex.main = 1.2)
grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")
abline(h = 0, col = "black", lty = "dotted", lwd = 2)
# Rho plot
tx = paste0("Rho Error vs. Sample Size (N) for rho =", sim_results[(i*9+1):((i+1)*9),]$rho[1])
plot(x = sim_results[(i*9+1):((i+1)*9),]$N, y = sim_results[(i*9+1):((i+1)*9),]$rho_error, type = "b", 
     col = "darkgreen",                    # Line color
     pch = 19,                        # Point symbol
     lwd = 2,                         # Line width
     xlab = "Number of Individuals (N)", # X-axis label
     ylab = "Rho Error",              # Y-axis label
     main = tx, # Title
     cex.axis = 0.9,                  # Axis text size
     cex.lab = 1,                     # Label text size
     cex.main = 1.2)
grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")
abline(h = 0, col = "black", lty = "dotted", lwd = 2)

#
i = 0
# Beta plot
plot(x = sim_results[(i*9+1):((i+1)*9),]$N, y = sim_results[(i*9+1):((i+1)*9),]$beta_error, type = "b", 
     col = "blue",                    
     pch = 19,                        
     lwd = 2,                         
     xlab = "Number of Individuals (N)", 
     ylab = "Beta Error",
     ylim = c(-0.1, 0.5),
     main = "Beta Error vs. Sample Size (N) for varying rho", 
     cex.axis = 0.9,                  
     cex.lab = 1,                   
     cex.main = 1.2)
grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")
abline(h = 0, col = "black", lty = "dotted", lwd = 2)
colors = c("lightblue", "grey", "green", "yellow", "orange", "red")
for (i in 1:length(colors)){
  points(x = sim_results[(i*9+1):((i+1)*9),]$N, y = sim_results[(i*9+1):((i+1)*9),]$beta_error,
         type = "b", pch = 19, col = colors[i], lwd = 2)
}



# Create the ggplot
ggplot(sim_results_long, aes(x = N, y = error, color = Coefficient, group = Coefficient)) +
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
    title = "Error of Coefficient Estimates Across N",
    x = "Number of Individuals (N)",
    y = "Error",
    color = "Coefficient") #+
  # annotate("text", x = max(sim_results$N) * 0.7, y = max(sim_results_long$error) * 0.8, 
  #          label = paste0("Tfull = ", unique(sim_results$Tfull), 
  #                         "\nrho = ", unique(sim_results$rho), 
  #                         "\nn_reps = ", unique(sim_results$n_reps)),
  #          hjust = 0, color = "black", size = 5, alpha = 0.7)
foo <- sim_results[,1:5] %>% group_by(rho)
ggplot(sim_results_long %>% filter(Coefficient == "beta_error"), aes(x = N, y = error, color = rho)) +
  geom_line(size = 1.2)

# Print the plot
print(error_plot)

# Create the ggplot
#ggplot(sim_results, aes(x = N, y = beta_error)) + geom_line(size = 1)
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
# Step 1 (rho -0.9 to -0.1)
rho = seq(from = -0.9, to = -0.1, by = 0.1)
Tfull = 5
N = 20000 # maybe 20000?
n_reps = 5
timer_error = T
beta_free = 1
omega_var = 1

parameters = create_grid(N, n_reps, Tfull, rho, timer_error)
(parameters)
set.seed(1) # estimation process for rho = 0 & seed = 1 is bugged, thus estimate in steps

# Step 2 (rho 0 to 0.9)
rho = seq(from = 0, to = 0.9, by = 0.1)
parameters = create_grid(N, n_reps, Tfull, rho, timer_error)
(parameters)
set.seed(2)

# Step 3 (rho -0.9 to 0.9 for Tfull 2 to 4)
rho = seq(from = -0.9, to = 0.9, by = 0.1)
Tfull = 2:4
parameters = create_grid(N, n_reps, Tfull, rho, timer_error)
(parameters)
set.seed(3)

# Step 4 (rho -0.9 to 0.9 for Tfull = 10)
rho = seq(from = -0.9, to = 0, by = 0.1)
Tfull = 10
N = 20000 # maybe 20000?
n_reps = 5
timer_error = T
parameters = create_grid(N, n_reps, Tfull, rho, timer_error)
(parameters)
set.seed(4)

# Step 5
rho = seq(from = 0.1, to = 0.7, by = 0.1)
Tfull = 10
parameters = create_grid(N, n_reps, Tfull, rho, timer_error)
(parameters)
set.seed(5)

# Step 6
rho = seq(from = 0.8, to = 0.9, by = 0.1)
parameters = create_grid(N, n_reps, Tfull, rho, timer_error)
(parameters)
set.seed(6)

# Fit
TrueVarInfo <- pmap(parameters, get_error_and_averages)
TrueVarResults = map_dfr(TrueVarInfo, as.data.frame) %>% select(beta_sd_avg, omega_sd_avg, rho_sd_avg)

# Extract the SD entries from the list
beta_sd_avg <- sapply(TrueVarInfo, function(x) x$beta_sd_avg)
omega_sd_avg <- sapply(TrueVarInfo, function(x) x$omega_sd_avg)
rho_sd_avg <- sapply(TrueVarInfo, function(x) x$rho_sd_avg)

# Combine them into data frame
sd_dataframe <- data.frame(
  rho = parameters$rho,
  beta_coef = parameters$beta_free,
  Tfull = parameters$Tfull,
  beta_sd_avg = beta_sd_avg,
  omega_sd_avg = omega_sd_avg,
  rho_sd_avg = rho_sd_avg,
  N = parameters$N,
  omega_var = parameters$omega_var
)

sd_dataframe
#old_path = getwd()
#setwd(dir = paste0(getwd(), "/GitHub/MasterThesis/TrueVarianceInfoCollection"))
#sd_dataframe2 = read.csv(file = paste0(getwd(),"/TrueVarianceDataframe.csv"))
#sd_dataframe = dplyr::bind_rows(sd_dataframe2, sd_dataframe)
#write.csv(sd_dataframe, file = "TrueVarianceDataframe.csv", row.names = F)
## ----- extracting from dfs
# old_path = getwd()
# setwd(dir = paste0(getwd(), "/GitHub/MasterThesis/ErrAvgSD_Folder"))
# rho = seq(from = 0.1, to = 0.7, by = 0.1)
# sd_missing <- data.frame(
#   rho = rho,
#   beta_coef = 1,
#   Tfull = 10,
#   beta_sd_avg = rep.int(NA, length(rho)),
#   omega_sd_avg = rep.int(NA, length(rho)),
#   rho_sd_avg = rep.int(NA, length(rho))
# )
# for (i in 1:length(rho)){
#   foo = read.csv(file = paste0(getwd(),"/ErrAvgSD_N=20000rho=", rho[i], "Tfull=10nreps=5.csv"))
#   AvgSd = c(mean(foo$beta_sd, na.rm = T),
#             mean(foo$omega_sd, na.rm = T),
#             mean(foo$rho_sd, na.rm = T))
#   sd_missing[i,4:6] = AvgSd
# }
# sd_missing
# setwd(dir = old_path)
# setwd(dir = paste0(getwd(), "/GitHub/MasterThesis/TrueVarianceInfoCollection"))
# sd_dataframe = read.csv(file = paste0(getwd(),"/TrueVarianceDataframe.csv"))
# sd_dataframe = dplyr::bind_rows(sd_dataframe, sd_missing)
# write.csv(sd_dataframe, file = "TrueVarianceDataframe.csv", row.names = F)
# setwd(dir = old_path)

# Obtain Confidence Interval stats ----------------------------------------
get_CI_results <- function(n_reps, Tfull, N, quest = 1, beta_free, rho, omega_var, alpha = 0.05, 
                           return_val = "", timer_error = T, timer_data = F, DontSkipFit = T, saveEstimates = T){
  
  if(timer_error){ 
    parameterMappedTo <- paste0("N=", N, ", Tfull=", Tfull, ", beta=", beta_free, ", rho=", rho, ", omega=", omega_var)
    msg = paste0("Fitting of Rprobit Models: ", parameterMappedTo, " , Finished: [", Sys.time(), "]")
    tic(msg = msg)
  }
  sim_replications <- replicate(n_reps, get_simulation_estimates(Tfull, # number of choice occasions (time)
                                                                 N, # individuals
                                                                 quest, # num of questions
                                                                 beta_free, # beta coefs; here: (free, fixed)
                                                                 rho, # rho in error AR(1)-process
                                                                 omega_var = omega_var,
                                                                 return_val = return_val, # return coefs & sds
                                                                 timer = timer_data,
                                                                 DontSkipFit = DontSkipFit)) 
  if(timer_error){
    # toc()
    foo = toc(quiet = TRUE)
    elapsed_time = foo$toc - foo$tic
    elapsed_minutes = round(elapsed_time/60, 2)
    cat(paste0("Elapsed fitting time (n_reps=", n_reps, "): ", elapsed_minutes, " min, Parameters: ", parameterMappedTo, "\n"))
    #cat(sprintf("Elapsed time: %.2f minutes", elapsed_minutes, "; Parameters: ", parameterMappedTo, "\n"))
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
    msg <- paste0("All nlm models (n_reps=", n_reps ,") were successfully estimated (N=", N, ", Tfull=", Tfull, ", rho=", rho, ", beta=", beta_free, ", omega=", omega_var, ")")
  } else {
    # Some (or all) models failed
    if (!DontSkipFit) {
      msg <- paste0("Skipped nlm model iterations (", n_reps, ") (N=", N, ", Tfull=", Tfull, ", rho=", rho, ", beta=", beta_free, ", omega=", omega_var, ")")
    } else {
      msg <- paste0("Some nlm models (", length(conv_failed_index), "/", n_reps, 
                    ") failed to converge (N=", N, ", Tfull=", Tfull, ", rho=", rho, ", beta=", beta_free, ", omega=", omega_var, ")")
    }
    # Adjust sim_replications to exclude failed indices
    sim_replications <- sim_replications[, -conv_failed_index]
  }
  print(msg)
  cat("\n")
  
  # Save true parameter
  true_parameters = c(beta = beta_free, omega = omega_var, rho = rho)
  
  # If all estimations failed/were skipped, return NA
  if(length(conv_success_index) == 0){
    
    results_table <- data.frame(
      Param = c("beta", "omega", "rho"),
      n_reps = c(rep.int(0, 3)),
      Success_Count = c(rep.int(NA, 3)),
      Success_rate = c(rep.int(NA, 3)),
      CI_Lower = c(rep.int(NA, 3)),
      True_param = true_parameters,
      CI_Upper = c(rep.int(NA, 3)),
      alpha = alpha,
      N = N,
      Tfull = Tfull,
      rho = rho,
      beta = beta_free,
      omega_var = omega_var,
      row.names = NULL
    )
    
    return(results_table)
  }
  
  # Get estimates
  beta_estimates <- sapply(sim_replications[1,], function(x) x["b_coef"])
  omega_estimates <- sapply(sim_replications[2,], function(x) x["omega_coef"])
  rho_estimates <- sapply(sim_replications[3,], function(x) x["rho_coef"])
  
  # Calculate the proportion of model estimates that lie within (1-alpha)%-CI
  
  # First, load in correct true variance file that was calculated on large N
  rho_to_filter <- rho
  beta_to_filter <- beta_free
  omega_to_filter <- omega_var
  Tfull_to_filter <- Tfull
  
  if(saveEstimates == T){
    # Save sds as well:
    beta_sds = sapply(sim_replications[1,], function(x) x["b_sd"])
    omega_sds <- sapply(sim_replications[2,], function(x) x["omega_sd"])
    rho_sds <- sapply(sim_replications[3,], function(x) x["rho_sd"])
    
    filename <- paste0("ConInt_", "N=", N, "Tfull=", Tfull, "nreps=", n_reps, "rho=", rho, "beta=", beta_free, "omega=", omega_var, ".csv")
    
    # Create df
    df = data.frame(beta = beta_estimates, 
                    beta_sd = beta_sds,
                    omega = omega_estimates,
                    omega_sd = omega_sds,
                    rho = rho_estimates,
                    rho_sd = rho_sds)
    for (i in 1:6){
      correct_tau = paste0("tau_coef", i)
      correct_tau_sd = paste0("tau_sd", i)
      tau_i_estimates = sapply(sim_replications[5,], function(x) x[correct_tau])
      tau_i_sd = sapply(sim_replications[5,], function(x) x[correct_tau_sd])
      
      old_row_count = dim(df)[2]
      df = cbind(df, tau_i_estimates, tau_i_sd)
      names(df)[(old_row_count+1):(old_row_count+2)] = c(paste0("tau", i), paste0("tau", i, "_sd"))
    }
    
    if(file.exists(paste0(getwd(),"/GitHub/MasterThesis/TrueVarianceInfoCollection"))){
      old_path = getwd()
      setwd(dir = paste0(getwd(),"/GitHub/MasterThesis/TrueVarianceInfoCollection"))
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
  
  # Loading true variance dataframe that contains all necessary info
  if(file.exists(paste0(getwd(),"/GitHub/MasterThesis/TrueVarianceInfoCollection"))){
    old_path = getwd()
    setwd(dir = paste0(getwd(),"/GitHub/MasterThesis/TrueVarianceInfoCollection"))
    true_var = read.csv(file = paste0(getwd(),"/TrueVarianceDataframe.csv"))
    setwd(dir = old_path)
  } else {
    msg = paste0("Error: TrueVariance DF file path inaccurate; File not found") 
    print(msg)
  }
  cat("\n")
  
  # Extract the correct info out of true var df
  true_var = dplyr::filter(true_var, rho == rho_to_filter,
                           beta_coef == beta_to_filter,
                           Tfull == Tfull_to_filter,
                           omega_var == omega_to_filter)
  if(dim(true_var)[1] == 0){
    msg = paste0("Error: Didn't find required true var (rho=", rho_to_filter,
                 ", beta=", beta_to_filter, ", Tfull=", Tfull_to_filter, ", omega=", omega_to_filter, "); Returning NA value instead")
    print(msg)
    cat("\n")
    
    results_table <- data.frame(
      Param = c("beta", "omega", "rho"),
      n_reps = c(rep.int(length(beta_estimates), 3)),
      Success_Count = c(rep.int(NA, 3)),
      Success_rate = c(rep.int(NA, 3)),
      CI_Lower = c(rep.int(NA, 3)),
      True_param = true_parameters,
      CI_Upper = c(rep.int(NA, 3)),
      alpha = alpha,
      N = N,
      Tfull = Tfull_to_filter,
      rho = rho_to_filter,
      beta = beta_to_filter,
      omega_var = omega_to_filter,
      row.names = NULL
    )
    return(results_table)
  }
  
  # Adjust true var accordingly to sample size 
  adjust_var_to_N <- function(unadjusted_true_var, N_new, N_old){
    return((N_old*unadjusted_true_var)/N_new)
  }
  
  N_to_load = true_var$N
  true_var <- true_var[1,4:6] 
  true_var <- true_var^2 # because info is SD so far
  names(true_var) = c("beta_true_var", "omega_true_var", "rho_true_var")
  
  adjusted_true_var <- adjust_var_to_N(true_var, N_new = N, N_old = N_to_load)
  
  # Calculate (1-alpha)%-CI boundaries
  if(alpha > 0 & alpha < 1){
    z_value = qnorm(1-alpha/2) # For different CIs
  } else {
    print("Error: Alpha value invalid")
  }
  
  ci_lower <- true_parameters - z_value * sqrt(adjusted_true_var)
  ci_lower <- unlist(ci_lower) # beta_true_var omega_true_var rho_true_var
  
  ci_upper <- true_parameters + z_value * sqrt(adjusted_true_var)
  ci_upper <- unlist(ci_upper)
  
  # Check cuccess or failure for each parameter
  success_beta <- sum(beta_estimates >= ci_lower[1] & beta_estimates <= ci_upper[1])
  failure_beta <- length(beta_estimates) - success_beta
  
  success_omega <- sum(omega_estimates >= ci_lower[2] & omega_estimates <= ci_upper[2])
  failure_omega <- length(omega_estimates) - success_omega
  
  success_rho <- sum(rho_estimates >= ci_lower[3] & rho_estimates <= ci_upper[3])
  failure_rho <- length(rho_estimates) - success_rho
  
  # Display results in table
  results_table <- data.frame(
    Param = c("beta", "omega", "rho"),
    n_reps = c(rep.int(length(beta_estimates), 3)),
    Success_Count = c(success_beta, success_omega, success_rho),
    Success_rate = c(success_beta, success_omega, success_rho)/length(beta_estimates),
    CI_Lower = ci_lower,
    True_param = true_parameters,
    CI_Upper = ci_upper,
    alpha = alpha,
    N = N,
    Tfull = Tfull_to_filter,
    rho = rho_to_filter,
    beta = beta_to_filter,
    omega_var = omega_to_filter,
    row.names = NULL
  )
  
  return(results_table)
}

get_CI_results(n_reps = 20, Tfull = 2, N = 100, quest = 1, beta_free = 1,
               rho = 0.55, omega_var = 1, alpha = 0.05, timer_error = T, timer_data = T, DontSkipFit = F, saveEstimates = T)

# ----- Calculate first set (hold everything but N)
n_reps = 200
Tfull = 2
N = c(10, 20, 30, 40, 50, 100)
rho = 0.5
omega_var = 1
beta_free = 1


parameters = create_grid(N, n_reps, Tfull, rho, beta_free, omega_var)
(parameters)

# Fit models

CI_overview <- pmap(parameters, get_CI_results)

# Extract info
#(CI_results <- dplyr::select(dplyr::bind_rows(CI_overview, .id = "Source"), -Source))
CI_results = map_dfr(CI_overview, as.data.frame)

# Save info
old_path = getwd()
setwd(dir = paste0(getwd(), "/GitHub/MasterThesis/TrueVarianceInfoCollection"))
write.csv(CI_results, file = "CIresultDataFrame.csv", row.names = F)
setwd(dir = old_path)


# Graphics for confidence intervals ---------------------------------------
library(tidyr)
library(ggplot2)

CI_results_beta = CI_results %>% dplyr::filter(Param == "beta")














# old stuff ---------------------------------------------------------------
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

