
# Load package
library(Rprobit)

# -------------------------------------------------------------------------
# Exemplified data generating process and package usage

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
                                     quest = 1, beta_coefs = c(1, 1), rho = 0.5, return_val = ""){
  
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
  
  # fit the system with cml_pair_type = 2 (only adjacent pairs + first and last decision)
  Rprob_mod <- fit_Rprobit(Rprobit_obj, init_method = "theta", cml_pair_type = 0)
  
  if(nzchar(return_val)){
    switch(return_val,
           "summary" = return(summary(Rprob_mod)),
           "var" = return(diag(Rprob_mod$vv)),
           "sd" = return(sqrt(diag(Rprob_mod$vv)))
    )
  } else {
    # Omega is Var-Cov, d.h. für Vergleich entweder sd quadrieren oder var wurzeln 
    b = c(b_coef = Rprob_mod$theta[1], b_sd = sqrt(diag(Rprob_mod$vv))[1])
    omega = c(omega_coef = Rprob_mod$theta[2], omega_sd = sqrt(diag(Rprob_mod$vv))[2]) # coef = omega sd entry, sd = sd of omega sd entry
    rho = c(rho_coef = Rprob_mod$theta[9], rho_sd = sqrt(diag(Rprob_mod$vv))[9])
    conv = c(nlm_code = Rprob_mod$info$nlm_info$code)
    #Rprob_conv = c(H_conv = T, J_conv = T)
    return(list(b, omega, rho, conv))
  }
}

# Obtain Confidence Interval stats ----------------------------------------
# Get true sd by fitting model to quasi-infinite sample (N extremely large; other parameters fixed)
rho = seq(from = -0.9, to = 0.9, by = 0.1)
N = 20000
Tfull = 5

for (r in rho){
  
  set.seed(1)
  true_var = get_simulation_estimates(Tfull = Tfull, # number of choice occasions (time)
                                      N = N, # individuals
                                      quest = 1, # num of questions
                                      beta_coefs = c(1, 1), # beta coefs; here: (free, fixed)
                                      rho = r, # rho in error AR(1)-process
                                      return_val = "var") # if return should be something else
  # Dynamic file name
  filename <- paste0("trueVarN=", N, "rho=", r, "Tfull=", Tfull, ".csv")
  write.csv(true_var, file = filename, row.names = F)
}

