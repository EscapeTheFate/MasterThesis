
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

