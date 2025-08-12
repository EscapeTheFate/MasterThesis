
# Load package
library(Rprobit)
library(tictoc)
library(purrr)
library(future.apply)



# This works --------------------------------------------------------------
# Preliminary settings:
Tfull = 3
N = 300
quest = 1
rho = 0.5
beta_free = c(1, 2) # Coefs of beta_1 and beta_2, beta_3 is fixed for normalisation/scale
omega_var = 2

# Modeling 
draw_simulation_data <- function(Tfull, N, quest, rho, system){
  
  time = c(1:Tfull)
  Tp = length(time)*quest # number of questions per person 
  
  # Setup dataframe by first individual
  latent_class = sample(x = c(1,2,3), size = 1, prob = c(rep.int(1/3,3)))
  latent_means = c(2,4,7)
  X = matrix(0,Tp,3)
  X[,1] = rnorm(n = Tp, mean = latent_means[latent_class], sd = 1) # b_1 is free
  X[,2] = rnorm(n = Tp, mean = 0, sd = 1)
  X[,3] = rnorm(n = Tp, mean = 0, sd = 1) # b_3 is fixed 
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
  df = data.frame(ID = 1, X1 = X[,1], X2 = X[,2], X3 = X[,3], y = y, time = time_df, quest = quest_df)
  
  # cycle over deciders 
  for (j in 2:N){
    latent_class = sample(x = c(1,2,3), size = 1, prob = c(rep.int(1/3,3)))
    X = matrix(0,Tp,3)
    X[,1] = rnorm(n = Tp, mean = latent_means[latent_class], sd = 1) # b_1 is free
    X[,2] = rnorm(n = Tp, mean = 0, sd = 1)
    X[,3] = rnorm(n = Tp, mean = 0, sd = 1) # b_3 is fixed 
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
    df = rbind(df, data.frame(ID = j, X1 = X[,1], X2 = X[,2], X3 = X[,3], y = y, time = time_df, quest = quest_df))
  } 
  
  return(df)
}

get_simulation_estimates <- function(Tfull, # number of choice occasions (time)
                                     N, # individuals
                                     quest = 1, beta_free, rho, omega_var, return_val = "", timer = F, DontSkipFit = T){
  mod <- mod_AR_cl$new(alt  = 7,
                       Hb   = diag(3)[,-3,drop=FALSE],
                       fb   = diag(3)[,-c(1:2),drop=FALSE],
                       HO   = matrix(1,1,1),
                       fO   = matrix(0,1,1),
                       HL   = matrix(0,1,0),
                       fL   = matrix(1,1,1),
                       ordered = TRUE
  )
  
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
  data_raw <- data_raw_StSp_cl$new(df = df, alt_names = c("1","2","3","4","5","6","7"), id = "ID", choice = "y", ordered = TRUE, varying = "", dec_char = c("X1","X2","X3"))
  data_raw$set_time_col("time")
  data_raw$set_quest_col("quest")
  data_raw$alt_names <- alt_names <- sprintf('%d',1:7)
  
  # set up Rprobit_obj: only type 2 regressors as usual in ordered probit situations. 
  form = y ~ 0 | X1 + X2 + X3 | 0
   
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
  
  Rprob_mod <- suppressMessages(fit_Rprobit(Rprobit_obj, init_method = "theta", cml_pair_type = 0))
  summary(Rprob_mod)
  
  return(summary(Rprob_mod))
}



# This doesn't work (now it does) -------------------------------------------------------

# Preliminary settings:
Tfull = 3
N = 300
quest = 1
rho = 0.5
beta_free = c(1, 2) # Coefs of beta_1 and beta_2, beta_3 is fixed for normalisation/scale
alpha_1 = 2 # var(beta_1)
alpha_3 = 1 # var(beta_2)
rho_beta = 0.2
# then, cov(beta_1, beta_2) = rho_beta * sqrt(Var(beta_1)) * sqrt(Var(beta_2))
# Thus, we set alpha_2 = Covariance entry of Omega for beta_1, beta_2 to be:

build_omega_var <- function(alpha_1, alpha_3, rho_beta){
  
  # Safety checks
  if(alpha_1 <= 0 || alpha_3 <= 0){
    stop("Variances alpha_1 and alpha_3 must be positive.")
  }
  
  if(abs(rho_beta) > 1){
    stop("Correlation rho must be between -1 and 1.")
  }
  
  alpha_2 <- rho_beta * sqrt(alpha_1 * alpha_3)
  omega_var <- matrix(c(alpha_1, alpha_2, alpha_2, alpha_3), nrow = 2, byrow = TRUE)
  
  return(omega_var)
}

omega_var = build_omega_var(alpha_1, alpha_3, rho_beta) # Omega matrix regarding beta_1 and beta_2
# Requirement: if omega_var matrix is {{alpha_1, alpha_2},{alpha_2, alpha_3}}, then alpha_1*alpha_3 > alpha_2^2 for chol-decomp. to work  

check_cholDecomp_cond <- function(omega_var){
  # Check if input is a 2x2-matrix
  if(!is.matrix(omega_var) || !all(dim(omega_var) == c(2, 2))){
    stop("omega_var must be a 2x2 matrix.")
  }
  
  # Check symmetry
  if(!isTRUE(all.equal(omega_var, t(omega_var)))){
    stop("omega_var must be symmetric.")
  }
  
  # Check positive definiteness (leading minor > 0 and determinant > 0)
  alpha1 <- omega_var[1, 1]
  alpha2 <- omega_var[1, 2]  # == omega_var[2, 1]
  alpha3 <- omega_var[2, 2]
  
  if(alpha1 <= 0){
    stop("Top-left entry must be positive: omega_var[1,1] > 0.")
  }
  
  determinant <- alpha1 * alpha3 - alpha2^2
  if(determinant <= 0){
    stop("Matrix is not positive definite: determinant <= 0.")
  }
  
  # All checks passed
  invisible(TRUE)
}

# Modeling 
draw_simulation_data <- function(Tfull, N, quest, rho, system){
  
  time = c(1:Tfull)
  Tp = length(time)*quest # number of questions per person 
  
  # Setup dataframe by first individual
  latent_class = sample(x = c(1,2,3), size = 1, prob = c(rep.int(1/3,3)))
  latent_means = c(2,4,7)
  X = matrix(0,Tp,3)
  X[,1] = rnorm(n = Tp, mean = latent_means[latent_class], sd = 1) # b_1 is free
  X[,2] = rnorm(n = Tp, mean = 0, sd = 1)
  X[,3] = rnorm(n = Tp, mean = 0, sd = 1) # b_3 is fixed 
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
  df = data.frame(ID = 1, X1 = X[,1], X2 = X[,2], X3 = X[,3], y = y, time = time_df, quest = quest_df)
  
  # cycle over deciders 
  for (j in 2:N){
    latent_class = sample(x = c(1,2,3), size = 1, prob = c(rep.int(1/3,3)))
    X = matrix(0,Tp,3)
    X[,1] = rnorm(n = Tp, mean = latent_means[latent_class], sd = 1) # b_1 is free
    X[,2] = rnorm(n = Tp, mean = 0, sd = 1)
    X[,3] = rnorm(n = Tp, mean = 0, sd = 1) # b_3 is fixed 
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
    df = rbind(df, data.frame(ID = j, X1 = X[,1], X2 = X[,2], X3 = X[,3], y = y, time = time_df, quest = quest_df))
  } 
  
  return(df)
}

get_simulation_estimates <- function(Tfull, # number of choice occasions (time)
                                     N, # individuals
                                     quest = 1, beta_free, rho, omega_var, return_val = "", timer = F, DontSkipFit = T){
  mod <- mod_AR_cl$new(alt  = 7,
                       Hb   = diag(3)[,-3,drop=FALSE],
                       fb   = diag(3)[,-c(1:2),drop=FALSE],
                       HO   = diag(3),
                       fO   = matrix(0,3,1),
                       HL   = matrix(0,1,0),
                       fL   = matrix(1,1,1),
                       ordered = TRUE
  )
  
  mod$lag_length = 1 # The p of the error AR(p)-process
  mod$stationary <- TRUE # indicates whether the stationary distribution shall be used to start the state process
  
  # Check if omega_var is correctly specified
  check_cholDecomp_cond(omega_var)
  
  # Setup for theta_0 true parameter value
  tau = c(0,rep.int(log(2),5))
  
  l11 <- sqrt(omega_var[1, 1])
  l21 <- omega_var[2, 1] / l11
  l22 <- sqrt(omega_var[2, 2] - l21^2)
  omega <- c(l11, l21, l22) # these are the starting values after all transformation/parametrization are done
  
  theta_0 = c(beta_free, omega, tau, rho) # (beta (-coef1), Omega (affin param), Sigma (-all), tau, rho)
  
  # Advance
  time = c(1:Tfull) # choice occasion ID
  Tp = length(time)*quest # number of questions per person 
  
  # Setup dataframe and system for DGP:
  system <- build_system_from_model_AR_R(theta_0, mod, time)
  df = draw_simulation_data(Tfull = Tfull, N = N, quest = quest, rho = rho, system = system)
  
  # Convert dataframe to data_raw_StSp_cl for estimation:
  data_raw <- data_raw_StSp_cl$new(df = df, alt_names = c("1","2","3","4","5","6","7"), id = "ID", choice = "y", ordered = TRUE, varying = "", dec_char = c("X1","X2","X3"))
  data_raw$set_time_col("time")
  data_raw$set_quest_col("quest")
  data_raw$alt_names <- alt_names <- sprintf('%d',1:7)
  
  # set up Rprobit_obj: only type 2 regressors as usual in ordered probit situations. 
  form = y ~ 0 | X1 + X2 + X3 | 0
  
  control_simul <- list(Tp = Tp)
  control <- list(control_simulation = control_simul)
  re <- c("X1", "X2")
  
  # generate Rprobit_obj 
  Rprobit_obj <- setup_Rprobit(form = form, data_raw = data_raw,
                               mod = mod,
                               re = re,
                               seed = 17,
                               theta_0 = theta_0,
                               control = control_simul
  )
  
  Rprob_mod <- suppressMessages(fit_Rprobit(Rprobit_obj, init_method = "theta", cml_pair_type = 0))
  summary(Rprob_mod)
  
  return(summary(Rprob_mod))
}

get_simulation_estimates(Tfull, N, quest = 1, beta_free, rho, omega_var = build_omega_var(alpha_1, alpha_3, rho_beta), return_val = "summary", timer = F, DontSkipFit = T)

# Running some tests
Tfull = 2:4
N = 500
quest = 1
rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
beta_1 = 1
beta_2 = 1
alpha_1 = 1 # var(beta_1)
alpha_3 = 1 # var(beta_2)
rho_beta = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8)
n_reps = 100

(parameters = create_grid(rho, n_reps, beta_1, beta_2, alpha_1, alpha_3, rho_beta, Tfull, N))

plan(multisession)
for (i in 1:5){
  plan(multisession)
  sim_replications <- future_replicate(n_reps, get_simulation_estimates(
    Tfull, N, quest, beta_free, rho, omega_var, return_val, timer_data, DontSkipFit),
    simplify = FALSE) 
}







