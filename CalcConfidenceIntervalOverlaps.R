
library(tidyverse)
library(dplyr)

error_results = read.csv(file = paste0(getwd(), "/GitHub/MasterThesis/DGP1_500reps_Results.csv"))
vars_of_importance = c("rho", "N", "n_reps", "beta_free", "omega_var", "Tfull")
results = error_results %>% select(all_of(vars_of_importance))

failed_rows = c(NULL)
list_of_dataframes <- list()

for(i in 1:nrow(results)){
  N = results$N[i]
  rho = results$rho[i]
  n_reps = results$n_reps[i]
  Tfull = results$Tfull[i]
  beta_free = results$beta_free[i]
  omega_var = results$omega_var[i]
  
  true_parameters = c(beta = beta_free, omega = omega_var, rho = rho)
  
  filename <- paste0("ErrAvgSD_N=", N, "_Tfull=", Tfull, "_nreps=", n_reps, 
                     "_rho=", rho, "_beta=", beta_free, "_omega=", omega_var, ".csv")
  
  if(file.exists(paste0(getwd(), "/GitHub/MasterThesis/Estimate_Collection/", filename))){
    foo = read.csv(file = paste0(getwd(), "/GitHub/MasterThesis/Estimate_Collection/", filename))
    foo = na.omit(foo)
    
  } else {
    msg = paste0("Failed row: ", i, ", Not found: ", filename)
    print(msg)
    
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
    
    if(length(failed_rows) == 0){
      failed_rows[1] = i
    } else {
      old_largest = length(failed_rows)
      failed_rows[old_largest+1] = i
    }
    next
  }
  
  
  
  beta_estimates <- foo$beta
  omega_estimates <- foo$omega
  rho_estimates <- foo$rho
  
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
  
  results_table <- data.frame(
    n_reps = n_reps,
    N = N,
    Tfull = Tfull,
    rho = rho,
    beta = beta_free,
    omega_var = omega_var,
    
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
    row.names = NULL
  )
  
  list_of_dataframes[[i]] <- results_table
}

combined_df <- do.call(rbind, list_of_dataframes)
View(combined_df)
foo = combined_df %>% select(-c(n_reps,
                                N,
                                Tfull,
                                rho,
                                beta,
                                omega_var))
combined_df = cbind(error_results, foo)
write.csv(combined_df, file = paste0(getwd(), "/GitHub/MasterThesis/DGP1_500reps_Results.csv"), row.names = F)






