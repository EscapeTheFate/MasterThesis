# -------------------------------------------------------------------------
# Master thesis - Code Space for Data Generating Process (DGP) of
# Ordered MNP-based Discrete Choice Models in Rprobit environment
# -------------------------------------------------------------------------
# Install package from source
#file.exists("~/Rprobit_0.3.2.tar.gz")
#install.packages("~/Rprobit_0.3.2.tar.gz", repos = NULL, type = "source")

# Load package
library(Rprobit)

# -------------------------------------------------------------------------
# Examplified data generating process and package usage

check_if_stationary <- function(phi){
  if(!is.vector(phi)){
    phi = as.vector(phi)
  }
  dim = length(phi)
  if(dim == 1){
    if(abs(phi[1]) < 1){
      return(T)
    }
    else{
      return(F)
    }
  }
  if(dim == 2){
    if(abs(phi[2]) < 1 & (phi[1] + phi[2]) < 1 & (phi[2] - phi[1]) < 1){
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

initialise_model_specification <- function(){
  
}

draw_simulation_data <- function(Tfull, N, quest, AR_coef, system){
  
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
  gam = LO %*% rnorm(lRE, sd = 0.5) # mixed part here
  
  get_AR_epsilon <- function(time, AR_coef, previous_error, Sigma_stationary = c(1), Sigma = c(1)){
    
    get_stationary_epsilon <- function(Sigma = c(1)){
      epsilon = rnorm(n=1, mean = 0, sd = Sigma)
      return(epsilon)
    }
    
    if(time == 1){
      epsilon_1 = get_stationary_epsilon(Sigma_stationary)
      return(epsilon_1)
    }
    else{
      u = get_stationary_epsilon(Sigma)
      epsilon_new = AR_coef * previous_error + u
      return(epsilon_new)
    }
  }
  
  # Obtain utilities for each choice occasion for individual 1
  U = matrix(0, nrow = Tp)
  # for (t in 1:Tp){
  #   if(t == 1){
  #     epsilon_old = get_AR_epsilon(time = t)
  #     U[t] = X[t,] %*% b + X[t,1] %*% beta[t] + epsilon_old
  #   }
  #   else{
  #     epsilon_new = get_AR_epsilon(time = t, AR_coef = system$factor, previous_error = epsilon_old)
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
    gam = LO %*% rnorm(lRE, sd = 0.5)
    
    # Obtain utilities for each choice occasion for individual 1
    U = matrix(0, nrow = Tp)
    # for (t in 1:Tp){
    #   if(t == 1){
    #     epsilon_old = get_AR_epsilon(time = t)
    #     U[t] = X[t,] %*% b + X[t,1] %*% beta[t] + epsilon_old
    #   }
    #   else{
    #     epsilon_new = get_AR_epsilon(time = t, AR_coef = system$factor, previous_error = epsilon_old)
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

#test = draw_simulation_data(Tfull = 5, N = 100, quest = 1, system = build_system_from_model_AR_R(theta_0, mod, c(1:Tfull)))
#table(test$y)
#View(test)

Tfull = 5
N = 100
quest = 1
AR_coef = 0.5

get_simulation_estimates <- function(Tfull = 5, # number of choice occasions (time)
                                     N = 100, # individuals
                                     quest = 1, AR_coef = 0.5){
  
  mod <- mod_AR_cl$new(alt  = 7,
                       Hb   = diag(2)[,-2,drop=FALSE],
                       fb   = matrix(0,2,1),
                       HO   = matrix(1,1,1),
                       fO   = matrix(0,1,1),
                       HL   = matrix(0,1,0),
                       fL   = matrix(1,1,1),
                       ordered = TRUE
  )
  
  mod$fb[2] = 2 # <- Since b = [b_1 (free), 2 (fixed)], we need to change fb entry
  mod$lag_length = 1 # The p of the error AR(p)-process
  mod$stationary <- TRUE # indicates whether the stationary distribution shall be used to start the state process
  
  tau = c(0,rep.int(log(2),5))
  theta_0 = c(1, 1, tau, 0.5) # (beta (-coef1), Omega, Sigma (-all), tau, Phi)
  
  # Advance
  time = c(1:Tfull) # choice occasion ID
  Tp = length(time)*quest # number of questions per person 
  
  # Setup dataframe and system for DGP:
  system <- build_system_from_model_AR_R(theta_0, mod, time)
  df = draw_simulation_data(Tfull = Tfull, N = N, quest = quest, AR_coef = AR_coef, system = system)
  
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
  summary(Rprob_mod)
  return(summary(Rprob_mod))
  
}

get_simulation_estimates(AR_coef = 0.5)



# -------------------------------------------------------------------------
# To-Do:
# - Get things to work
# - Extract parameters of interest of Rprobit mod output
# - Do the same for all grid values
# - Improve get_simultion_estimates to only depend on draw_simulation_data
# - Slot in stability check infront of parameter grid
get_bias <- function(n_reps = 500)

rho_grid <- seq(from = -0.95, to = 0.95, by = 0.5)
Tp_grid <- seq(from = 2, to = 20, by = 5)
parameters = expand.grid(rho = rho_grid, Tp = Tp_grid)

bias <- Map(get_bias, n = parameters$n, Tp = parameters$Tp)








