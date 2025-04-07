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
    }
    cat("\n")
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

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20000)
n_reps = c(1)
beta_free = c(seq(from = 0, to = 0.9, by = 0.1))
omega_var = 1
Tfull = c(2:4)
set.seed(31) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(1.8, 1.9, 2)
omega_var = 1
Tfull = c(2:4)
set.seed(37) # set.seed(40) starting from row 856, adjust later for prior rows

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(2.4, 2.5)
omega_var = 1
Tfull = c(2:4)
set.seed(41) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(1.1, 1.2, 1.3)
omega_var = 0.25^2
Tfull = c(2:4)
set.seed(45)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(1.4, 1.5, 1.6)
omega_var = 0.25^2
Tfull = c(2:4)
set.seed(47) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(2.1, 2.2)
omega_var = 0.25^2
Tfull = c(2:4)
set.seed(49) # set.seed(50) starting from row 589 (first entry of parameters with N = 1000), adjust entries later 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(2.3, 2.4, 2.5)
omega_var = 0.25^2
Tfull = c(2:4)
set.seed(51)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) # N = 20, 30, 40 missing due to instability 
n_reps = c(500)
beta_free = c(0, 0.1, 0.2)
omega_var = 4^2
Tfull = c(2:4)
set.seed(52)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) # N = 20, 30, 40 missing due to instability 
n_reps = c(500)
beta_free = c(0.3, 0.4)
omega_var = 4^2
Tfull = c(2:4)
set.seed(53)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) # N = 20, 30, 40 missing due to instability 
n_reps = c(500)
beta_free = c(0.5, 0.6, 0.7)
omega_var = 4^2
Tfull = c(2:4)
set.seed(54) # set.seed(56) starting from row 217, adjust later

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) # N = 20, 30, 40 missing due to instability 
n_reps = c(500)
beta_free = c(1.7)
omega_var = 4^2
Tfull = c(2:4)
set.seed(58) # set.seed(59) starting from row 123, adjust later

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) # N = 20, 30, 40 missing due to instability 
n_reps = c(500)
beta_free = c(2, 2.1, 2.2)
omega_var = 4^2
Tfull = c(2:4)
set.seed(60) # set.seed(62) starting from row 182 (=ind+1), adjust later

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(3, 3.5)
omega_var = 1
Tfull = c(2:4)
set.seed(61)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) # 20:40 missing, adjust later
n_reps = c(500)
beta_free = c(2.6, 2.7)
omega_var = 4^2
Tfull = c(2:4)
set.seed(70)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) 
n_reps = c(500)
beta_free = c(2.6, 2.7)
omega_var = c(1, 0.25^2)
Tfull = c(2:4)
set.seed(73)

# rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
# N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) 
# n_reps = c(2500)
# beta_free = c(1)
# omega_var = c(1)
# Tfull = c(2:4)
# set.seed(77) # set.seed(81) starting from row 281 (=ind+1), adjust later

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(2500)
beta_free = c(1)
omega_var = c(1)
Tfull = c(2:4)
set.seed(84) # recalced

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(2500)
beta_free = c(1)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(82)

## On Work PC --------------------------------------------------------------
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

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(1.6, 1.7)
omega_var = 1
Tfull = c(2:4) # c(2:5)
set.seed(38) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(2.7)
omega_var = 1
Tfull = c(2:4) # c(2:5)
set.seed(42) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(0, 0.1, 0.2, 0.3, 0.4)
omega_var = 0.25^2
Tfull = c(2:4) # c(2:5)
set.seed(43) # set.seed(46) starting from row 234

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(0.8, 0.9, 1.1)
omega_var = 4^2
Tfull = c(2:4)
set.seed(55) 


rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(1.8, 1.9)
omega_var = 4^2
Tfull = c(2:4)
set.seed(59) # set.seed(61) starting from row 89 (=ind+1), adjust later
(parameters = create_grid(rho, n_reps, beta_free, omega_var, Tfull, N))
(ind = which(parameters$N == 70 & parameters$Tfull == 2 & parameters$beta == 1.8 & parameters$rho == 0))
parameters = parameters[(ind+1):nrow(parameters),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(2.6, 2.7, 2.8)
omega_var = 4^2
Tfull = c(2:4)
set.seed(69) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(2.8, 2.9)
omega_var = c(1, 0.25^2)
Tfull = c(2:4)
set.seed(74) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(1000) # adjust later 
n_reps = c(2500)
beta_free = c(1)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(79)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(3)
omega_var = c(1)
Tfull = c(2:4)
set.seed(189)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200)
n_reps = c(500)
beta_free = c(3)
omega_var = c(0.25^2)
Tfull = c(2:4)
set.seed(196)

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

## Omega Variance = 16 low sample calcs ------------------------------------
rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40)
n_reps = c(500)
beta_free = c(0)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(85) # set.seed(87) starting from row 7, set.seed(88) starting from row 35
parameters = parameters[-c(1:35),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40)
n_reps = c(500)
beta_free = c(0.1)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(89) # set.seed(90) starting from row 7, set.seed(91) starting from row 14
parameters = parameters[-c(1:13),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40)
n_reps = c(500)
beta_free = c(0.2)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(92) # set.seed(94) starting from row 6, set.seed(95) starting from row 28
parameters = parameters[-c(1:27),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30)
n_reps = c(500)
beta_free = c(0.3)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(100) # set.seed(102) starting from row 14, set.seed(104) starting from row 27
parameters = parameters[-c(1:27),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30)
n_reps = c(500)
beta_free = c(0.4)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(106) # set.seed(108) starting from row 6
parameters = parameters[-c(1:5),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30)
n_reps = c(500)
beta_free = c(0.5)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(112) # set.seed(114) starting from row 7, set.seed(116) starting from row 21
parameters = parameters[-c(1:20),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30)
n_reps = c(500)
beta_free = c(0.6)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(121) # set.seed(122) starting from row 7, set.seed(124) starting from row 11
parameters = parameters[-c(1:10),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30)
n_reps = c(500)
beta_free = c(0.7)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(126) # set.seed(127) starting from row 6, set.seed(128) starting from row 7
parameters = parameters[-c(1:6),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30)
n_reps = c(500)
beta_free = c(0.8)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(132) # set.seed(135) starting from row 7, set.seed(139) starting from row 14
parameters = parameters[-c(1:13),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30)
n_reps = c(500)
beta_free = c(0.9)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(146) # set.seed(148) starting from row 7
parameters = parameters[-c(1:6),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30)
n_reps = c(500)
beta_free = c(1)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(151) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30)
n_reps = c(500)
beta_free = c(1.1)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(152) # set.seed(153) starting from row 7
parameters = parameters[-c(1:6),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30)
n_reps = c(500)
beta_free = c(1.2)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(159)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30)
n_reps = c(500)
beta_free = c(1.3)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(162) # set.seed(164) starting from row 13, set.seed(168) starting from row 28
parameters = parameters[-c(1:27),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(40)
n_reps = c(2500)
beta_free = c(1)
omega_var = c(4^2)
Tfull = c(2:4)
set.seed(168) # set.seed(164) starting from row 13, set.seed(168) starting from row 28
parameters = parameters[-c(1:27),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40)
n_reps = c(500)
beta_free = c(1)
omega_var = c(0.25^2)
Tfull = c(2:4)
set.seed(168) # set.seed(164) starting from row 13, set.seed(168) starting from row 28
parameters = parameters[-c(1:27),]

# The last remaining parameter sets missing from the estimate collection
rho = c(0.9, 0.9, 0.9, 0.6, 0.9, 0.6)
N = c(30, 20, 40, 40, 40, 20)
n_reps = c(500)
beta_free = c(0, 0.5, 0.6, 0.7, 1.8, 2.1)
omega_var = c(4^2)
Tfull = c(3, 2, 2, 2, 2, 4)
set.seed(202) # set.seed(203) starting from row 6
(parameters = data.frame(rho, n_reps, beta_free, omega_var, Tfull, N)) # <- not create_grid here
parameters = parameters[-c(1:5),]

rho = c(0.6, 0.9, 0.6, 0.9)
N = c(20, 30, 20, 20)
n_reps = c(500)
beta_free = c(2.9, 2.9, 3.0, 3.0)
omega_var = c(4^2)
Tfull = c(2, 2, 3, 4)
set.seed(205)
(parameters = data.frame(rho, n_reps, beta_free, omega_var, Tfull, N)) # <- not create_grid here

rho = c(0.9, -0.9, 0.9, -0.9)
N = c(20000)
n_reps = c(1)
beta_free = c(0, 0.4, 0.4, 0.6)
omega_var = c(0.25^2, 4^2, 4^2, 4^2)
Tfull = c(4, 3, 3, 4)
set.seed(207)
(parameters = data.frame(rho, n_reps, beta_free, omega_var, Tfull, N)) # <- not create_grid here

rho = c(0.9, -0.3, 0.9)
N = c(20000)
n_reps = c(1)
beta_free = c(1.1, 1.2, 1.3)
omega_var = c(1)
Tfull = c(4, 3, 4)
set.seed(208)
(parameters = data.frame(rho, n_reps, beta_free, omega_var, Tfull, N)) # <- not create_grid here

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20000)
n_reps = c(1)
beta_free = c(1.6, 1.7, 1.8, 1.9)
omega_var = c(1)
Tfull = c(2:4)
set.seed(209)

rho = c(-0.9, 0.6, 0.3, 0)
N = c(20000)
n_reps = c(1)
beta_free = c(1.7, 1.9, 2.7, 2.8)
omega_var = c(1,4^2, 1, 0.25^2)
Tfull = c(3, 3, 3, 3)
set.seed(214)
(parameters = data.frame(rho, n_reps, beta_free, omega_var, Tfull, N)) # <- not create_grid here

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70) 
n_reps = c(500)
beta_free = c(4)
omega_var = c(0.25^2, 1, 4^2)
Tfull = c(2:4)
set.seed(217) # set.seed(219) starting from row 4, set.seed(220) starting from row 39
# set.seed(221) starting from row 45, set.seed(222) starting from row 49
# set.seed(223) starting from row 89
which(parameters$N == 30 & parameters$Tfull == 3 & parameters$rho == 0 & parameters$omega_var == 0.0625)
parameters = parameters[-c(1:88),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(100, 125, 150, 200, 300, 500) 
n_reps = c(2500)
beta_free = c(0)
omega_var = c(0.25^2, 1, 4^2)
Tfull = c(2:4)
set.seed(224) # set.seed(230) starting from row 288, set.seed(241) starting from row 351
# set.seed(242) starting from row 364
(ind = which(parameters$N == 500 & parameters$Tfull == 4 & parameters$rho == 0.6 & parameters$omega_var == 0.0625))
parameters = parameters[-c(1:363),]


rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(1000, 500, 300, 200, 150, 125, 100, 90, 80, 70, 60, 50) 
n_reps = c(500)
beta_free = c(3.2, 3.4)
omega_var = c(0.25^2, 1, 4^2)
Tfull = c(2:4)
set.seed(251) # continue at first entry of 500

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(150, 125, 100) 
n_reps = c(2500)
beta_free = c(4)
omega_var = c(0.25^2, 1, 4^2)
Tfull = c(2:4)
set.seed(263)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(500) 
n_reps = c(2500)
beta_free = c(4)
omega_var = c(0.25^2, 1, 4^2)
Tfull = c(4)
set.seed(267)

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40) 
n_reps = c(2500)
beta_free = c(0)
omega_var = c(0.25^2, 1, 4^2)
Tfull = c(2:4)
set.seed(271) # set.seed(274) starting from row 15, set.seed(279) starting from row 20, set.seed(280) starting from row 21
# set.seed(288) starting from row 62, set.seed(295) starting from row 75
(ind = which(parameters$N == 30 & parameters$Tfull == 2 & parameters$rho == 0 & parameters$omega_var == 1))
parameters = parameters[-c(1:74),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20) 
n_reps = c(2500)
beta_free = c(3)
omega_var = c(0.25^2, 1, 4^2)
Tfull = c(2:3)
set.seed(309) # set.seed(312) starting from row 14, set.seed(315) starting from row 21, set.seed(317) starting from row 22, N=20, Tfull=3, beta=3, rho=0.3, omega=0.0625
(ind = which(parameters$N == 20 & parameters$Tfull == 2 & parameters$rho == 0.6 & parameters$omega_var == 16))
parameters = parameters[-c(1:21),]

#
# True Variance calcs -----------------------------------------------------
rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20000) 
n_reps = c(1)
beta_free = c(seq(from = 0, to = 0.9, by = 0.1))
omega_var = 1
Tfull = c(2:4) 
set.seed(31) # set.seed(35) starting from row 159

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20000) 
n_reps = c(1)
beta_free = c(seq(from = 0, to = 1, by = 0.1))
omega_var = c(0.25^2, 4^2)
Tfull = c(2:4) 
set.seed(64) # set.seed(67) starting from row 328 (=ind+1), no adjustment needed

(ind = which(parameters$rho == 0.3 & parameters$Tfull == 4 & parameters$beta_free == 0.2 & parameters$omega_var == 0.0625))
parameters = parameters[(ind+1):nrow(parameters),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20000)
n_reps = c(1)
beta_free = c(seq(from = 1, to = 2, by = 0.1))
omega_var = c(0.25^2, 4^2)
Tfull = c(2:4)
set.seed(65) # set.seed(66) starting from row 222 (=ind+1), no adjustment needed
ind = which(parameters$Tfull == 3 & parameters$rho == 0 & parameters$beta == 1.9 & parameters$omega == 0.0625)
parameters = parameters[(ind+1):nrow(parameters),]

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20000)
n_reps = c(1)
beta_free = c(seq(from = 2.1, to = 2.5, by = 0.1))
omega_var = c(0.25^2, 4^2)
Tfull = c(2:4)
set.seed(68) 

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20000)
n_reps = c(1)
beta_free = c(seq(from = 2.6, to = 3, by = 0.1))
omega_var = c(0.25^2, 1, 4^2)
Tfull = c(2:4)
set.seed(80) # there is a note on Home PC with seed 76, but 80 was used here according to history

# Fit grid and obtain results ---------------------------------------------

(parameters = create_grid(rho, n_reps, beta_free, omega_var, Tfull, N))

# Fit models to grid
error <- pmap(parameters, get_error_and_averages)
(error_results = map_dfr(error, as.data.frame))
(error_results = cbind(error_results, parameters))

# Save calculated dataframe for later
file.exists(paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"))
previous_error_results = read.csv(file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"))
error_results = dplyr::bind_rows(previous_error_results, error_results)
write.csv(error_results, file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"), row.names = F)


## Recalc results from estimate collection ---------------------------------

# Parameter sets that need to be included in results:
rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = 1
omega_var = c(16, 25)
Tfull = c(2:4) 
(parameters = create_grid(rho, n_reps, beta_free, omega_var, Tfull, N))
set.seed(12) # set.seed(14) starting from ro 34: parameters = parameters[34:nrow(parameters),]
# set.seed(15) starting from row 63, set.seed(16) starting from row 92

rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(1.8, 1.9, 2)
omega_var = 1
Tfull = c(2:4)
(parameters = create_grid(rho, n_reps, beta_free, omega_var, Tfull, N))
set.seed(37) # set.seed(40) starting from row 856, adjust later for prior rows


# Recalc code
list_of_dataframes <- list()

for (i in 1:nrow(parameters)){
  N = parameters$N[i]
  Tfull = parameters$Tfull[i]
  n_reps = parameters$n_reps[i]
  rho = parameters$rho[i]
  beta_free = parameters$beta_free[i]
  omega_var = parameters$omega_var[i]
  
  filename <- paste0("ErrAvgSD_N=", N, "_Tfull=", Tfull, "_nreps=", n_reps, 
                     "_rho=", rho, "_beta=", beta_free, "_omega=", omega_var, ".csv")
  
  if(!file.exists(paste0(getwd(), "/GitHub/MasterThesis/Estimate_Collection/", filename))){
    msg = paste0("Error occured in iteration:", i)
    print(msg)
    break
  }
  foo = read.csv(paste0(getwd(), "/GitHub/MasterThesis/Estimate_Collection/", filename))
  
  # Get estimates
  beta_estimates <- foo$beta
  omega_estimates <- foo$omega
  rho_estimates <- foo$rho
  
  true_parameters <- c(beta = beta_free, omega = omega_var, rho = rho)
  
  # Get mean error
  beta_error <- mean(beta_estimates - true_parameters["beta"], na.rm = T)
  omega_error <- mean(omega_estimates - true_parameters["omega"], na.rm = T) # we are comparing omega_var entry with estimated omega_var parameter
  rho_error <- mean(rho_estimates - true_parameters["rho"], na.rm = T)
  
  # Now, calculate average standard deviation for each estimate
  beta_sd_avg <- mean(foo$beta_sd, na.rm = T)
  omega_sd_avg <- mean(foo$omega_sd, na.rm = T) # these are sd's of omega_var
  rho_sd_avg <- mean(foo$rho_sd, na.rm = T)
  
  list_of_dataframes[[i]] = data.frame(beta_error = beta_error, 
                                 beta_sd_avg = beta_sd_avg,
                                 omega_error = omega_error,
                                 omega_sd_avg = omega_sd_avg,
                                 rho_error = rho_error,
                                 rho_sd_avg = rho_sd_avg,
                                 rho = rho,
                                 N = N,
                                 n_reps = n_reps,
                                 beta_free = beta_free,
                                 omega_var = omega_var,
                                 Tfull = Tfull,
                                 beta_estimate_sd = NA,
                                 omega_estimate_sd = NA,
                                 rho_estimate_sd = NA)
  
}

combined_df <- do.call(rbind, list_of_dataframes)

file.exists(paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"))
previous_error_results = read.csv(file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv")) # previous_error_results <- foo %>% mutate(across(where(is.character), as.numeric))
error_results = dplyr::bind_rows(previous_error_results, combined_df)
write.csv(error_results, file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"), row.names = F)

## Add in results from prof pc ---------------------------------------------
file.exists(paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"))
previous_error_results = read.csv(file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv")) # previous_error_results <- foo %>% mutate(across(where(is.character), as.numeric))
prof_pc_results = read.csv(file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_large_Nreps_omega_var16.csv"))
error_results = dplyr::bind_rows(previous_error_results, prof_pc_results)
write.csv(error_results, file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"), row.names = F)


# --------------------------------------------------------------------------------
# Obtain Results from Source files (build entire results.csv from scratch --------
library(dplyr)
library(tidyverse)

# All (relevant) parameter combinations that were estimated prior to this: 
rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500)
beta_free = c(seq(from = 0, to = 3, by = 0.1))
omega_var = c(0.25^2, 1, 4^2)
Tfull = c(2:4) 
(parameters = create_grid(rho, n_reps, beta_free, omega_var, Tfull, N))

# Check which models are still missing/need to be recalced due to NA entry (i.e. true var df with 1 nreps)
list_of_dataframes <- list()
vector_of_missings <- c()
counter = 0

for (i in 1:nrow(parameters)){
  N = parameters$N[i]
  Tfull = parameters$Tfull[i]
  n_reps = parameters$n_reps[i]
  rho = parameters$rho[i]
  beta_free = parameters$beta_free[i]
  omega_var = parameters$omega_var[i]
  
  filename <- paste0("ErrAvgSD_N=", N, "_Tfull=", Tfull, "_nreps=", n_reps, 
                     "_rho=", rho, "_beta=", beta_free, "_omega=", omega_var, ".csv")
  
  if(!file.exists(paste0(getwd(), "/GitHub/MasterThesis/Estimate_Collection/", filename))){
    msg = paste0("Error occured in iteration:", i)
    print(msg)
    counter = 1
    list_of_missing[[counter]] = i
    next
  }
  foo = read.csv(paste0(getwd(), "/GitHub/MasterThesis/Estimate_Collection/", filename))
  
  # Get estimates
  beta_estimates <- foo$beta
  omega_estimates <- foo$omega
  rho_estimates <- foo$rho
  
  true_parameters <- c(beta = beta_free, omega = omega_var, rho = rho)
  
  # Get mean error
  beta_error <- mean(beta_estimates - true_parameters["beta"], na.rm = T)
  omega_error <- mean(omega_estimates - true_parameters["omega"], na.rm = T) # we are comparing omega_var entry with estimated omega_var parameter
  rho_error <- mean(rho_estimates - true_parameters["rho"], na.rm = T)
  
  # Now, calculate average standard deviation for each estimate
  beta_sd_avg <- mean(foo$beta_sd, na.rm = T)
  omega_sd_avg <- mean(foo$omega_sd, na.rm = T) # these are sd's of omega_var
  rho_sd_avg <- mean(foo$rho_sd, na.rm = T)
  
  list_of_dataframes[[i]] = data.frame(beta_error = beta_error, 
                                       beta_sd_avg = beta_sd_avg,
                                       omega_error = omega_error,
                                       omega_sd_avg = omega_sd_avg,
                                       rho_error = rho_error,
                                       rho_sd_avg = rho_sd_avg,
                                       rho = rho,
                                       N = N,
                                       n_reps = n_reps,
                                       beta_free = beta_free,
                                       omega_var = omega_var,
                                       Tfull = Tfull,
                                       beta_estimate_sd = NA,
                                       omega_estimate_sd = NA,
                                       rho_estimate_sd = NA)
  
}
combined_df <- do.call(rbind, list_of_dataframes)

# Save calculated dataframe
write.csv(combined_df, file = paste0(getwd(), "/GitHub/MasterThesis/DGP1_500reps_Results.csv"), row.names = F)


# Normal models with 500 reps
for (i in 1:nrow(parameters)){
  N = parameters$N[i]
  Tfull = parameters$Tfull[i]
  n_reps = parameters$n_reps[i]
  rho = parameters$rho[i]
  beta_free = parameters$beta_free[i]
  omega_var = parameters$omega_var[i]
  
  filename <- paste0("ErrAvgSD_N=", N, "_Tfull=", Tfull, "_nreps=", n_reps, 
                     "_rho=", rho, "_beta=", beta_free, "_omega=", omega_var, ".csv")
  
  if(!file.exists(paste0(getwd(), "/GitHub/MasterThesis/Estimate_Collection/", filename))){
    msg = paste0("Error occured in iteration:", i)
    print(msg)
    counter = counter + 1
    vector_of_missings[counter] = i
    next
  }
}
vector_of_missings = as.vector(unlist(vector_of_missings))
parameters = parameters[vector_of_missings,]
View(parameters %>% arrange(beta_free))

# True variance NA filtering (searching which results need to be recalced due
# to them being calc'ed for one model only, but with large sample which may 
# be crashed/not converged at the time/had NA result for whatever reason)
rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20000)
n_reps = c(1)
beta_free = c(seq(from = 0, to = 3, by = 0.1))
omega_var = c(0.25^2, 1, 4^2)
Tfull = c(2:4) 
(parameters = create_grid(rho, n_reps, beta_free, omega_var, Tfull, N))

# True variance models with 1 rep
list_of_dataframes <- list()
vector_of_missings <- c()
vector_of_NA_rows <- c()
counter = counter_NA = 0

for (i in 1:nrow(parameters)){
  N = parameters$N[i]
  Tfull = parameters$Tfull[i]
  n_reps = parameters$n_reps[i]
  rho = parameters$rho[i]
  beta_free = parameters$beta_free[i]
  omega_var = parameters$omega_var[i]
  
  filename <- paste0("ErrAvgSD_N=", N, "_Tfull=", Tfull, "_nreps=", n_reps, 
                     "_rho=", rho, "_beta=", beta_free, "_omega=", omega_var, ".csv")
  
  if(!file.exists(paste0(getwd(), "/GitHub/MasterThesis/True_Variance_Estimate_Collection/", filename))){
    msg = paste0("Error occured in iteration:", i)
    print(msg)
    counter = counter + 1
    vector_of_missings[counter] = i
    next
  }
  
  df <- read.csv(paste0(getwd(), "/GitHub/MasterThesis/True_Variance_Estimate_Collection/", filename))
  
  if(nrow(df) == 1 && all(is.na(df))) {
    msg = paste0("NA occured in iteration:", i)
    print(msg)
    counter_NA = counter_NA + 1
    vector_of_NA_rows[counter_NA] = i
    next
  }
}
vector_of_missings = as.vector(unlist(vector_of_missings))
vector_of_NA_rows = as.vector(unlist(vector_of_NA_rows))
parameters = parameters[vector_of_missings,]
View(parameters %>% arrange(beta_free))

file.exists(paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"))
previous_error_results = read.csv(file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv")) # previous_error_results <- foo %>% mutate(across(where(is.character), as.numeric))
error_results = dplyr::bind_rows(previous_error_results, combined_df)
write.csv(error_results, file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"), row.names = F)

# Testing for the appropriate replication sample size for beta = 1 value:
rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000)
n_reps = c(500, 2500)
beta_free = c(0, 1, 2, 3, 4)
omega_var = c(0.25^2, 1, 4^2)
Tfull = c(2:4) 
nrep_samplesizes = c(100, 200, 300, 400, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000)
(parameters = create_grid(n_reps, rho, nrep_samplesizes, beta_free, omega_var, Tfull, N))




# Recalc code where the 2500 and 500 (independent) model results are first combined,
# then the first nrep_samplesizes number's are taken out of the 3000 (i.e. first 
# 500 from 3000, and so on...), and the results are calculated on that basis
list_of_dataframes <- list()
vector_of_missings <- c()
counter = 0

for (i in (seq(from = 0, to = dim(parameters)[1], by = 2) + 1)){
  flag_1 = flag_2 = 0
  N = parameters$N[i]
  Tfull = parameters$Tfull[i]
  n_reps = parameters$n_reps[i]
  rho = parameters$rho[i]
  beta_free = parameters$beta_free[i]
  omega_var = parameters$omega_var[i]
  nreps_samplesize_number = parameters$nrep_samplesizes[i]
  
  filename1 <- paste0("ErrAvgSD_N=", N, "_Tfull=", Tfull, "_nreps=", n_reps, 
                     "_rho=", rho, "_beta=", beta_free, "_omega=", omega_var, ".csv")
  
  if(!file.exists(paste0(getwd(), "/GitHub/MasterThesis/Estimate_Collection/", filename1))){
    msg = paste0("Error occured in iteration:", i, " (filename1)")
    print(msg)
    counter = counter + 1
    vector_of_missings[counter] = i
    flag_1 = 1
  }
  
  N = parameters$N[i+1]
  Tfull = parameters$Tfull[i+1]
  n_reps = parameters$n_reps[i+1]
  rho = parameters$rho[i+1]
  beta_free = parameters$beta_free[i+1]
  omega_var = parameters$omega_var[i+1]
  
  filename2 <- paste0("ErrAvgSD_N=", N, "_Tfull=", Tfull, "_nreps=", n_reps, 
                      "_rho=", rho, "_beta=", beta_free, "_omega=", omega_var, ".csv")
  
  if(!file.exists(paste0(getwd(), "/GitHub/MasterThesis/Estimate_Collection/", filename2))){
    msg = paste0("Error occured in iteration:", i, " (filename2)")
    print(msg)
    counter = counter + 1
    vector_of_missings[counter] = i
    flag_2 = 1
  }
  if(parameters$nrep_samplesizes[i] != parameters$nrep_samplesizes[i+1]){
    msg = c("2 different nreps samplesize numbers after each other in 'parameters'")
    print(msg)
    break
  }
  
  if(flag_1 != 1 & flag_2 != 1){
    foo1 = read.csv(paste0(getwd(), "/GitHub/MasterThesis/Estimate_Collection/", filename1)) # 500 nreps one
    foo2 = read.csv(paste0(getwd(), "/GitHub/MasterThesis/Estimate_Collection/", filename2)) # 2500 nreps one
    foo = rbind(foo1, foo2) # combine
    foo = foo[1:nreps_samplesize_number,]# take only the first nreps_samplesize number of rows
    foo = na.omit(foo)
  } else {
      if(flag_1 == 1){
        foo = read.csv(paste0(getwd(), "/GitHub/MasterThesis/Estimate_Collection/", filename2))
        if(nreps_samplesize_number <= dim(foo)[1]){
          foo = foo[1:nreps_samplesize_number,]
        }
      } else {
        foo = read.csv(paste0(getwd(), "/GitHub/MasterThesis/Estimate_Collection/", filename1)) # 500 nreps one
        if(nreps_samplesize_number <= dim(foo)[1]){
          foo = foo[1:nreps_samplesize_number,]
        }
      }
    foo = na.omit(foo)
  }
  
  
  # Get estimates
  beta_estimates <- foo$beta
  omega_estimates <- foo$omega
  rho_estimates <- foo$rho
  
  true_parameters <- c(beta = beta_free, omega = omega_var, rho = rho)
  
  # Get mean error
  beta_error <- mean(beta_estimates - true_parameters["beta"], na.rm = T)
  omega_error <- mean(omega_estimates - true_parameters["omega"], na.rm = T) # we are comparing omega_var entry with estimated omega_var parameter
  rho_error <- mean(rho_estimates - true_parameters["rho"], na.rm = T)
  
  # Now, calculate average standard deviation for each estimate
  beta_sd_avg <- mean(foo$beta_sd, na.rm = T)
  omega_sd_avg <- mean(foo$omega_sd, na.rm = T) # these are sd's of omega_var
  rho_sd_avg <- mean(foo$rho_sd, na.rm = T)
  
  # Also, calculate standard deviation for standard deviation estimates of 
  # parameter estimates that will be used as a CI for avg estimate SD plots
  sd_of_beta_sds <- sd(foo$beta_sd, na.rm = T)
  sd_of_omega_sds <- sd(foo$omega_sd, na.rm = T) # these are sd's of omega_var
  sd_of_rho_sds <- sd(foo$rho_sd, na.rm = T)
  
  # Calculate SD of the Sample of Estimates (so on beta, omega and rho estimates)
  # This uses the code in CalcSampleEstimateSD, but directly implemented
  
  beta_estimate_sd = sd(beta_estimates, na.rm = T)
  omega_estimate_sd = sd(omega_estimates, na.rm = T)
  rho_estimate_sd = sd(rho_estimates, na.rm = T)
  
  # Calculate the Confidence Intervals + the success rate of estimates lying
  # within the interval (same as code in CalcConfidenceIntervalOverlaps, but here)
  
  # Calculate the proportion of model estimates that lie within (1-alpha)%-CI
  
  # First, load in correct true variance file that was calculated on large N
  N_to_load = 20000
  truevar_file <- paste0("ErrAvgSD_N=", N_to_load, "_Tfull=", Tfull, "_nreps=", 1, 
                         "_rho=", rho, "_beta=", beta_free, "_omega=", omega_var, ".csv")
  
  # Loading true variance dataframe that contains all necessary info
  if( file.exists(paste0(getwd(),"/GitHub/MasterThesis/True_Variance_Estimate_Collection/", truevar_file)) ){
    true_var_large_N = read.csv(file = paste0(getwd(),"/GitHub/MasterThesis/True_Variance_Estimate_Collection/", truevar_file))
    true_var_large_N = true_var_large_N %>% select(beta_sd, omega_sd, rho_sd) # <----------------------------------------------------------------------------------
    true_var_large_N = true_var_large_N^2 # because info is SD so far
    names(true_var_large_N) = c("beta_true_var", "omega_true_var", "rho_true_var")
  } else {
    msg = paste0("Error: TrueVariance DF file path inaccurate and/or file not found") 
    print(msg)
    cat("\n")
    next
    
    results_table <- data.frame(
      n_reps = n_reps,
      N = N,
      Tfull = Tfull,
      rho = rho,
      beta = beta_free,
      omega_var = omega_var,
      
      ci_lower_beta = NA,
      ci_upper_beta = NA,
      success_beta = NA,
      failure_beta = NA,
      p_beta = NA,
      bino_ci_norm_lower_beta = NA,
      bino_ci_exact_lower_beta = NA,
      bino_ci_wilson_lower_beta = NA,
      bino_ci_norm_upper_beta = NA,
      bino_ci_exact_upper_beta = NA,
      bino_ci_wilson_upper_beta = NA,
      
      ci_lower_omega = NA,
      ci_upper_omega = NA,
      success_omega = NA,
      failure_omega = NA,
      p_omega = NA,
      bino_ci_norm_lower_omega = NA,
      bino_ci_exact_lower_omega = NA,
      bino_ci_wilson_lower_omega = NA,
      bino_ci_norm_upper_omega = NA,
      bino_ci_exact_upper_omega = NA,
      bino_ci_wilson_upper_omega = NA,
      
      
      ci_lower_rho = NA,
      ci_upper_rho = NA,
      success_rho = NA,
      failure_rho = NA,
      p_rho = NA,
      bino_ci_norm_lower_rho = NA,
      bino_ci_exact_lower_rho = NA,
      bino_ci_wilson_lower_rho = NA,
      bino_ci_norm_upper_rho = NA,
      bino_ci_exact_upper_rho = NA,
      bino_ci_wilson_upper_rho = NA,
      row.names = NULL
    )
    
    list_of_dataframes[[i]] <- results_table
    
    next
  }
  
  
  # # Extract the correct info out of true var df
  # true_var = dplyr::filter(true_var, rho == rho_to_filter,
  #                          beta_coef == beta_to_filter,
  #                          Tfull == Tfull_to_filter,
  #                          omega_var == omega_to_filter)
  
  adjust_var_to_N <- function(unadjusted_true_var, N_new, N_old){
    return((N_old*unadjusted_true_var)/N_new)
  }
  
  adjusted_true_var_large_N <- adjust_var_to_N(true_var_large_N, N_new = N, N_old = N_to_load)
  alpha = 0.05
  
  if(alpha > 0 & alpha < 1){
    z_value = qnorm(1-alpha/2) # For different CIs
  } else {
    print("Error: Alpha value invalid")
  }
  
  ci_lower <- true_parameters - z_value * sqrt(adjusted_true_var_large_N)
  ci_lower <- unlist(ci_lower) # beta_true_var omega_true_var rho_true_var
  
  ci_upper <- true_parameters + z_value * sqrt(adjusted_true_var_large_N)
  ci_upper <- unlist(ci_upper)
  
  # Check cuccess or failure for each parameter
  success_beta <- sum(beta_estimates >= ci_lower[1] & beta_estimates <= ci_upper[1])
  failure_beta <- length(beta_estimates) - success_beta
  
  success_omega <- sum(omega_estimates >= ci_lower[2] & omega_estimates <= ci_upper[2])
  failure_omega <- length(omega_estimates) - success_omega
  
  success_rho <- sum(rho_estimates >= ci_lower[3] & rho_estimates <= ci_upper[3])
  failure_rho <- length(rho_estimates) - success_rho
  
  # Obtain uncertainty (upper and lower) of CI-overlap result by normal approximation
  eps = 1e-10
  p_beta = success_beta/length(beta_estimates)
  p_omega = success_omega/length(omega_estimates)
  p_rho = success_rho/length(rho_estimates)
  
  binomial_norm_upper_beta = p_beta + z_value * sqrt( (p_beta*(1-p_beta) + eps)/length(beta_estimates) )
  binomial_norm_upper_omega = p_omega + z_value * sqrt( (p_omega*(1-p_omega) + eps)/length(omega_estimates) )
  binomial_norm_upper_rho = p_rho + z_value * sqrt( (p_rho*(1-p_rho) + eps)/length(rho_estimates) )
  
  binomial_norm_lower_beta = p_beta - z_value * sqrt( (p_beta*(1-p_beta) + eps)/length(beta_estimates) )
  binomial_norm_lower_omega = p_omega - z_value * sqrt( (p_omega*(1-p_omega) + eps)/length(omega_estimates) )
  binomial_norm_lower_rho = p_rho - z_value * sqrt( (p_rho*(1-p_rho) + eps)/length(rho_estimates) )
  
  # Now, do the same, but to obtain exact ones (by clopper and pearson)
  binomial_exact_upper_beta = 1 - qbeta(p = alpha/2, shape1 = length(beta_estimates) - success_beta, shape2 = success_beta + 1)
  binomial_exact_upper_omega = 1 - qbeta(p = alpha/2, shape1 = length(omega_estimates) - success_omega, shape2 = success_omega + 1)
  binomial_exact_upper_rho = 1 - qbeta(p = alpha/2, shape1 = length(rho_estimates) - success_rho, shape2 = success_rho + 1)
  
  binomial_exact_lower_beta = 1 - qbeta(p = 1 - alpha/2, shape1 = length(beta_estimates) - success_beta + 1, shape2 = success_beta)
  binomial_exact_lower_omega = 1 - qbeta(p = 1 - alpha/2, shape1 = length(omega_estimates) - success_omega + 1, shape2 = success_omega)
  binomial_exact_lower_rho = 1 - qbeta(p = 1 - alpha/2, shape1 = length(rho_estimates) - success_rho + 1, shape2 = success_rho)
  
  # Compute Wilson Score Interval
  wilson_ci <- function(p, n, z) {
    denominator <- 1 + (z^2 / n)
    center <- p + (z^2 / (2 * n))
    margin <- z * sqrt((p * (1 - p) / n) + (z^2 / (4 * n^2)))
    
    lower <- (center - margin) / denominator
    upper <- (center + margin) / denominator
    
    return(c(lower, upper))
  }
  
  # Compute Wilson CI for all estimates
  wilson_beta <- wilson_ci(p_beta, length(beta_estimates), z_value)
  wilson_omega <- wilson_ci(p_omega, length(omega_estimates), z_value)
  wilson_rho <- wilson_ci(p_rho, length(rho_estimates), z_value)
  
  
  list_of_dataframes[[i]] = data.frame(beta_error = beta_error, 
                                       beta_sd_avg = beta_sd_avg,
                                       sd_of_beta_sds = sd_of_beta_sds,
                                       omega_error = omega_error,
                                       omega_sd_avg = omega_sd_avg,
                                       sd_of_omega_sds = sd_of_omega_sds,
                                       rho_error = rho_error,
                                       rho_sd_avg = rho_sd_avg,
                                       sd_of_rho_sds = sd_of_rho_sds,
                                       rho = rho,
                                       N = N,
                                       n_reps = nreps_samplesize_number,
                                       beta_free = beta_free,
                                       omega_var = omega_var,
                                       Tfull = Tfull,
                                       
                                       beta_estimate_sd = beta_estimate_sd,
                                       omega_estimate_sd = omega_estimate_sd,
                                       rho_estimate_sd = rho_estimate_sd,
                                       
                                       ci_lower_beta = ci_lower[1],
                                       ci_upper_beta = ci_upper[1],
                                       success_beta = success_beta,
                                       failure_beta = failure_beta,
                                       p_beta = p_beta,
                                       bino_ci_norm_lower_beta = binomial_norm_lower_beta,
                                       bino_ci_exact_lower_beta = binomial_exact_lower_beta,
                                       bino_ci_wilson_lower_beta = wilson_beta[1],
                                       bino_ci_norm_upper_beta = binomial_norm_upper_beta,
                                       bino_ci_exact_upper_beta = binomial_exact_upper_beta,
                                       bino_ci_wilson_upper_beta = wilson_beta[2],
                                       
                                       ci_lower_omega = ci_lower[2],
                                       ci_upper_omega = ci_upper[2],
                                       success_omega = success_omega,
                                       failure_omega = failure_omega,
                                       p_omega = p_omega,
                                       bino_ci_norm_lower_omega = binomial_norm_lower_omega,
                                       bino_ci_exact_lower_omega = binomial_exact_lower_omega,
                                       bino_ci_wilson_lower_omega = wilson_omega[1],
                                       bino_ci_norm_upper_omega = binomial_norm_upper_omega,
                                       bino_ci_exact_upper_omega = binomial_exact_upper_omega,
                                       bino_ci_wilson_upper_omega = wilson_omega[2],
                                       
                                       
                                       ci_lower_rho = ci_lower[3],
                                       ci_upper_rho = ci_upper[3],
                                       success_rho = success_rho,
                                       failure_rho = failure_rho,
                                       p_rho = p_rho,
                                       bino_ci_norm_lower_rho = binomial_norm_lower_rho,
                                       bino_ci_exact_lower_rho = binomial_exact_lower_rho,
                                       bino_ci_wilson_lower_rho = wilson_rho[1],
                                       bino_ci_norm_upper_rho = binomial_norm_upper_rho,
                                       bino_ci_exact_upper_rho = binomial_exact_upper_rho,
                                       bino_ci_wilson_upper_rho = wilson_rho[2],
                                       row.names = NULL)
  
}
combined_df <- do.call(rbind, list_of_dataframes)
View(combined_df) 
View(combined_df[seq(from = 1, to = 57, by = 7),]) # compare various sample sizes for N = 20, rho = -0.9, T = 2
View(combined_df[seq(from = 5230, to = 5286, by = 7),]) # compare various sample sizes for N = 125, rho = -0.9, T = 3, omega = 16

# Save calculated dataframe
write.csv(combined_df, file = paste0(getwd(), "/GitHub/MasterThesis/DGP1_3000reps_Results.csv"), row.names = F)



