#
# Master thesis - Code Space for Data Generating Process (DGP) of
# Ordered MNP-based Discrete Choice Models in Rprobit environment
#

## Install package from source ---------------------------------------------
#file.exists("~/Rprobit_0.3.2update.tar.gz")
#install.packages("~/Rprobit_0.3.2update.tar.gz", repos = NULL, type = "source")

## Load package
library(Rprobit)
library(tictoc)
library(purrr)
library(future.apply)

## Testing set
# Tfull = 3
# N = 300
# quest = 1
# rho = 0.5
# beta_free = c(1, 2) # Coefs of beta_1 and beta_2, beta_3 is fixed for normalisation/scale
# alpha_1 = 2 # var(beta_1)
# alpha_3 = 1 # var(beta_2)
# rho_beta = 0.2
# omega_var = build_omega_var(alpha_1, alpha_3, rho_beta)

## Data-Generating-Process Structure

create_grid <- function(...) {
  vars <- list(...)
  var_names <- sapply(substitute(list(...))[-1], deparse)
  setNames(expand.grid(vars), var_names)
}

build_omega_var <- function(alpha_1, alpha_3, rho_beta){
  
  # Here, we assume that Omega is 2x2-matrix, so that it has the following 
  # structure of: alpha_1 = ω_11, alpha_2 = ω_12 = ω_21, alpha_3 = ω_22
  
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

extract_rho_vector <- function(alpha_1_vec, alpha_2_vec, alpha_3_vec){
  mapply(function(a1, a2, a3) a2 / sqrt(a1 * a3),
         a1 = alpha_1_vec,
         a2 = alpha_2_vec,
         a3 = alpha_3_vec)
}

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

draw_simulation_data_dgp2 <- function(Tfull, N, quest, rho, system){
  
  time = c(1:Tfull)
  Tp = length(time)*quest # number of questions per person 
  
  # Setup dataframe by first individual
  latent_class = sample(x = c(1,2,3), size = 1, prob = c(rep.int(1/3,3)))
  latent_means = c(2,4,7)
  X = matrix(0,Tp,3)
  X[,1] = rnorm(n = Tp, mean = latent_means[latent_class], sd = 1) # b_1 is free
  X[,2] = rnorm(n = Tp, mean = 0, sd = 1) # b_2 also free
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

get_simulation_estimates_dgp2 <- function(Tfull, # number of choice occasions (time)
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
  
  l11 <- sqrt(omega_var[1,1])
  l21 <- omega_var[2,1]/l11
  l22 <- sqrt(omega_var[2,2]-l21^2)
  omega <- c(l11, l21, l22) # these are the starting values after all transformation/parametrization are done
  
  theta_0 = c(beta_free, omega, tau, rho) # (beta (-coef1), Omega (affin param), Sigma (-all), tau, rho)
  
  # Advance
  time = c(1:Tfull) # choice occasion ID
  Tp = length(time)*quest # number of questions per person 
  
  # Setup dataframe and system for DGP:
  system <- build_system_from_model_AR_R(theta_0, mod, time)
  df = draw_simulation_data_dgp2(Tfull = Tfull, N = N, quest = quest, rho = rho, system = system)
  
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
                 "theta" = replace(Rprob_mod$theta, 3:5, # replace omega entries by delta method, yielding omega entries:
                                   c(Rprob_mod$theta[3]^2, # Var(X2) in omega (not SE or Var of param)
                                     Rprob_mod$theta[3]*Rprob_mod$theta[4], # Cov(X1,X2) in omega (")
                                     Rprob_mod$theta[4]^2+Rprob_mod$theta[5]^2) # Var(X1) in omega (")
                                   ), 
                 "var" = replace(diag(Rprob_mod$vv), 3:5, # replace second entry (omega) by delta method
                                 c(4*Rprob_mod$theta[3]^2*diag(Rprob_mod$vv)[3], # <-- these are variance entries, so no need to square (we can take sqrt for SE)
                                   Rprob_mod$theta[4]^2 * Rprob_mod$vv[3,3] + Rprob_mod$theta[3]^2 * Rprob_mod$vv[4,4] + 2 * Rprob_mod$theta[3] * Rprob_mod$theta[4] * Rprob_mod$vv[3,4],
                                   4*Rprob_mod$theta[4]^2 * Rprob_mod$vv[4,4] + 4*Rprob_mod$theta[5]^2 * Rprob_mod$vv[5,5] + 8*Rprob_mod$theta[4] * Rprob_mod$theta[5] * Rprob_mod$vv[4,5])
                                 ), 
                 "sd" = replace(sqrt(diag(Rprob_mod$vv)), 3:5, # replace second entry (omega) by delta method
                                sqrt(c(4*Rprob_mod$theta[3]^2*diag(Rprob_mod$vv)[3],
                                  Rprob_mod$theta[4]^2 * Rprob_mod$vv[3,3] + Rprob_mod$theta[3]^2 * Rprob_mod$vv[4,4] + 2 * Rprob_mod$theta[3] * Rprob_mod$theta[4] * Rprob_mod$vv[3,4],
                                  4*Rprob_mod$theta[4]^2 * Rprob_mod$vv[4,4] + 4*Rprob_mod$theta[5]^2 * Rprob_mod$vv[5,5] + 8*Rprob_mod$theta[4] * Rprob_mod$theta[5] * Rprob_mod$vv[4,5]))
                               )
                 )
        } else {
          # Omega is Var-Cov, d.h. für Vergleich entweder sd quadrieren oder var wurzeln 
          b2 = c(b2_coef = Rprob_mod$theta[1], b2_sd = sqrt(diag(Rprob_mod$vv))[1])
          b1 = c(b1_coef = Rprob_mod$theta[2], b1_sd = sqrt(diag(Rprob_mod$vv))[2])
          omega_b1 = c(omega_b1_coef = Rprob_mod$theta[4]^2+Rprob_mod$theta[5]^2, omega_b1_sd = sqrt(4*Rprob_mod$theta[4]^2 * Rprob_mod$vv[4,4] + 4*Rprob_mod$theta[5]^2 * Rprob_mod$vv[5,5] + 8*Rprob_mod$theta[4] * Rprob_mod$theta[5] * Rprob_mod$vv[4,5]))
          omega_b2 = c(omega_b2_coef = Rprob_mod$theta[3]^2, omega_b2_sd = sqrt(c(4*Rprob_mod$theta[3]^2*diag(Rprob_mod$vv)[3]))) # coef_omega = theta_omega^2 (omega matrix entry), sd = delta method (sd of omega mat entry)
          omega_cov = c(omega_cov_coef = Rprob_mod$theta[3]*Rprob_mod$theta[4], omega_cov_sd = sqrt(Rprob_mod$theta[4]^2 * Rprob_mod$vv[3,3] + Rprob_mod$theta[3]^2 * Rprob_mod$vv[4,4] + 2 * Rprob_mod$theta[3] * Rprob_mod$theta[4] * Rprob_mod$vv[3,4]))
          rho = c(rho_coef = Rprob_mod$theta[12], rho_sd = sqrt(diag(Rprob_mod$vv))[12])                                
          conv = c(nlm_code = Rprob_mod$info$nlm_info$code)
          tau_raw = c(tau_coef = Rprob_mod$theta[6:11], tau_sd = sqrt(diag(Rprob_mod$vv)[6:11]))
          list(b1, b2, omega_b1, omega_b2, omega_cov, rho, conv, tau_raw)
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
          b2 = c(b2_coef = NA, b2_sd = NA)
          b1 = c(b1_coef = NA, b1_sd = NA)
          omega_b1 = c(omega_b1_coef = NA, omega_b1_sd = NA)
          omega_b2 = c(omega_b2_coef = NA, omega_b2_sd = NA) # coef_omega = theta_omega^2 (omega matrix entry), sd = delta method (sd of omega mat entry)
          omega_cov = c(omega_cov_coef = NA, omega_cov_sd = NA)
          rho = c(rho_coef = NA, rho_sd = NA)                                
          conv = c(nlm_code = 99)
          tau_raw = c(tau_coef = rep.int(NA, 6), tau_sd = rep.int(NA, 6))
          list(b1, b2, omega_b1, omega_b2, omega_cov, rho, conv, tau_raw)
          
          message("Model fitting failed: ", conditionMessage(e))
          
          # Return the same structure with NA values
          list(b = b, omega = omega, rho = rho, conv = conv, tau_raw = tau_raw)
        }
      }
    )
  } else { # if skipped, return NA same as if error occured
    b2 = c(b2_coef = NA, b2_sd = NA)
    b1 = c(b1_coef = NA, b1_sd = NA)
    omega_b1 = c(omega_b1_coef = NA, omega_b1_sd = NA)
    omega_b2 = c(omega_b2_coef = NA, omega_b2_sd = NA) # coef_omega = theta_omega^2 (omega matrix entry), sd = delta method (sd of omega mat entry)
    omega_cov = c(omega_cov_coef = NA, omega_cov_sd = NA)
    rho = c(rho_coef = NA, rho_sd = NA)                                
    conv = c(nlm_code = 99)
    tau_raw = c(tau_coef = rep.int(NA, 6), tau_sd = rep.int(NA, 6))
    result = list(b1, b2, omega_b1, omega_b2, omega_cov, rho, conv, tau_raw)
  }
  if(timer == T){
    toc()
  }
  # Return results 
  return(result)
}

#get_simulation_estimates_dgp2(Tfull, N, quest = 1, beta_free = c(1, 1), rho, omega_var, return_val = "summary", timer = F, DontSkipFit = T)





# Obtain error stats ------------------------------------------------------

get_error_and_averages_dgp2 <- function(n_reps, Tfull, N, quest = 1, beta_free, rho, omega_var,
                                   return_val = "", timer_error = T, timer_data = F, DontSkipFit = T, saveEstimates = T){
  
  if(timer_error){ 
    parameterMappedTo <- paste0("N=", N, ", Tfull=", Tfull, ", beta=c(", toString(beta_free), ")", ", rho=", rho, ", omega=c(", toString(omega_var[lower.tri(omega_var, diag = TRUE)]),
                                ", cor=", extract_rho_vector(omega_var[1,1], omega_var[2,1], omega_var[2,2]) ,"))")
    msg = paste0("Fitting of Rprobit Models: ", parameterMappedTo)
    tic(msg = msg)
  }
  
  plan(multisession)
  sim_replications <- future_replicate(n_reps, get_simulation_estimates_dgp2(
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
  conv_status <- sapply(sim_replications, function(x) x[[7]]["nlm_code"])
  conv_success_index <- conv_status %in% c(1, 2, 3)
  conv_failed_index <- conv_status %in% c(4, 5, 99) # code 99 for nlm crash
  
  if (sum(conv_failed_index) == 0 & DontSkipFit == TRUE) {
    msg <- paste0(
      "All nlm models (n_reps=", n_reps, 
      ") were successfully estimated (N=", N, 
      ", Tfull=", Tfull, 
      ", rho=", rho, 
      ", beta=c(", toString(beta_free), ")", 
      ", omega=c(", toString(omega_var[lower.tri(omega_var, diag = TRUE)]),
      ", cor=", extract_rho_vector(omega_var[1,1], omega_var[2,1], omega_var[2,2]) ,"))"
    )
  } else {
    if (!DontSkipFit) {
      msg <- paste0(
        "Skipped nlm model iterations (", n_reps, 
        ") (N=", N, 
        ", Tfull=", Tfull, 
        ", rho=", rho, 
        ", beta=c(", toString(beta_free), ")",
        ", omega=c(", toString(omega_var[lower.tri(omega_var, diag = TRUE)]),
        ", cor=", extract_rho_vector(omega_var[1,1], omega_var[2,1], omega_var[2,2]) ,"))"
      )
    } else {
      msg <- paste0(
        "Some nlm models (", sum(conv_failed_index), "/", n_reps, 
        ") failed to converge (N=", N, 
        ", Tfull=", Tfull, 
        ", rho=", rho, 
        ", beta=c(", toString(beta_free), ")", 
        ", omega=c(", toString(omega_var[lower.tri(omega_var, diag = TRUE)]),
        ", cor=", extract_rho_vector(omega_var[1,1], omega_var[2,1], omega_var[2,2]) ,"))"
      )
    }
  }
  print(msg)
  cat("\n")
  
  if(saveEstimates){
    filename <- paste0("ErrAvgSD_N=", N, "_Tfull=", Tfull, "_nreps=", n_reps, 
                       "_rho=", rho, "_beta1=", beta_free[1], "_beta2=", beta_free[2],
                       "_omegab1=", omega_var[1,1], "_omegab2=", omega_var[2,2],
                       "_omegacor=", extract_rho_vector(omega_var[1,1], omega_var[2,1], omega_var[2,2]), ".csv")
    
    # Get estimates and build dataframe at the same time (not as before), calculate metrics later
    df <- data.frame(
      beta1       = sapply(sim_replications, function(x) x[[1]]["b1_coef"]),
      beta1_sd    = sapply(sim_replications, function(x) x[[1]]["b1_sd"]),
      beta2       = sapply(sim_replications, function(x) x[[2]]["b2_coef"]),
      beta2_sd    = sapply(sim_replications, function(x) x[[2]]["b2_sd"]),
      omega_b1    = sapply(sim_replications, function(x) x[[3]]["omega_b1_coef"]),
      omega_b1_sd = sapply(sim_replications, function(x) x[[3]]["omega_b1_sd"]),
      omega_b2    = sapply(sim_replications, function(x) x[[4]]["omega_b2_coef"]),
      omega_b2_sd = sapply(sim_replications, function(x) x[[4]]["omega_b2_sd"]),
      omega_cov   = sapply(sim_replications, function(x) x[[5]]["omega_cov_coef"]),
      omega_cov_sd= sapply(sim_replications, function(x) x[[5]]["omega_cov_sd"]),
      rho         = sapply(sim_replications, function(x) x[[6]]["rho_coef"]),
      rho_sd      = sapply(sim_replications, function(x) x[[6]]["rho_sd"]),
      conv_code   = conv_status
    )
    
    omega_cor_estimates <- extract_rho_vector(
      df$omega_b1, df$omega_cov, df$omega_b2
    )
    df$omega_cor <- extract_rho_vector(df$omega_b1, df$omega_cov, df$omega_b2)
    
    for (i in 1:6){
      # Extract (raw) tau_i + it's (raw) sd (so, just sqrt of diagonal varcov-mat) in iteration i
      correct_tau = paste0("tau_coef", i)
      correct_tau_sd = paste0("tau_sd", i)
      
      tau_i_estimates <- sapply(sim_replications, function(x) x[[8]][correct_tau])
      tau_i_sd <- sapply(sim_replications, function(x) x[[8]][correct_tau_sd])
      
      old_col_count = dim(df)[2]
      df = cbind(df, tau_i_estimates, tau_i_sd)
      names(df)[(old_col_count+1):(old_col_count+2)] = c(paste0("tau", i), paste0("tau", i, "_sd"))
    }
    
    if(N == 10000 | N == 20000){
      if(file.exists(paste0(getwd(),"/GitHub/MasterThesis/DGP2/True_Variance_Estimate_Collection"))){
        old_path = getwd()
        setwd(dir = paste0(getwd(),"/GitHub/MasterThesis/DGP2/True_Variance_Estimate_Collection"))
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
      if(file.exists(paste0(getwd(),"/GitHub/MasterThesis/DGP2/Estimate_Collection"))){
        old_path = getwd()
        setwd(dir = paste0(getwd(),"/GitHub/MasterThesis/DGP2/Estimate_Collection"))
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
}

# get_error_and_averages_dgp2(n_reps = 100, Tfull = 3, N = 100, quest = 1, beta_free = c(1, 1), rho = 0.3, omega_var = build_omega_var(alpha_1, alpha_3, rho_beta),
#                        return_val = "", timer_error = T, timer_data = T,
#                        DontSkipFit = T, saveEstimates = T)


## Extract error and average sd for different models --------------------------

# Define grid level to map get_error on (before on all pcs):
rho = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
N = c(500, 300, 200, 100) # c(20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000) 
n_reps = c(3000)
beta1 = 1
beta2 = 1
Tfull = c(2:4)
alpha_1 = c(0.25^2, 0.5^2, 0.75^2, 1) # var(beta_1)
alpha_3 = c(0.25^2, 0.5^2, 0.75^2, 1) # var(beta_2)
rho_beta = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
set.seed(500) 

# Fit grid and obtain results ---------------------------------------------
(parameters = create_grid(rho, n_reps, beta1, beta2, alpha_1, alpha_3, rho_beta, Tfull, N))

# Iterate over parameter df dimension (different from DGP1 which was based on 
# pmap due to legacy reasons of calc'ing metrics immediately and not afterwards)
for (i in 1:5){ # dim(parameters)[1]
  get_error_and_averages_dgp2(n_reps = 100,
                              Tfull = parameters$Tfull[i],
                              N = parameters$N[i],
                              quest = 1,
                              beta_free = c(parameters$beta1[i], parameters$beta2[i]),
                              rho = parameters$rho[i],
                              omega_var = build_omega_var(parameters$alpha_1[i],
                                                          parameters$alpha_3[i],
                                                          parameters$rho_beta[i]),
                              return_val = "",
                              timer_error = T,
                              timer_data = F,
                              DontSkipFit = T,
                              saveEstimates = T)
}
error <- pmap(parameters, get_error_and_averages)


# -------------------------------------------------------------------------


