
library(tidyverse)
library(dplyr)

error_results = read.csv(file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"))
results = error_results %>% select(names(error_results)[7:ncol(error_results)])
#results = cbind(results, conv_1_2 = NA, conv_3 = NA, conv_4_5 = NA, crash = NA)

failed_rows = c(NULL)
list_of_dataframes <- list()

for(i in 1:nrow(results)){
  N = results$N[i]
  rho = results$rho[i]
  n_reps = results$n_reps[i]
  Tfull = results$Tfull[i]
  beta_free = results$beta_free[i]
  omega_var = results$omega_var[i]
  
  filename <- paste0("ErrAvgSD_N=", N, "_Tfull=", Tfull, "_nreps=", n_reps, 
                     "_rho=", rho, "_beta=", beta_free, "_omega=", omega_var, ".csv")
  
  if(file.exists(paste0(getwd(), "/GitHub/MasterThesis/Estimate_Collection/", filename))){
    foo = read.csv(file = paste0(getwd(), "/GitHub/MasterThesis/Estimate_Collection/", filename))
    foo = na.omit(foo)
    
  } else {
    msg = paste0("Failed row: ", i, ", Not found: ", filename)
    print(msg)
    
    if(length(failed_rows) == 0){
      failed_rows[1] = i
    } else {
      old_largest = length(failed_rows)
      failed_rows[old_largest+1] = i
    }
    next
  }
  
  true_parameters = c(beta = beta_free, omega = omega_var, rho = rho)
  
  beta_estimates <- foo$beta
  omega_estimates <- foo$omega
  rho_estimates <- foo$rho
  
  # Calculate the proportion of model estimates that lie within (1-alpha)%-CI
  
  # First, load in correct true variance file that was calculated on large N
  # rho_to_filter <- rho
  # beta_to_filter <- beta_free
  # omega_to_filter <- omega_var
  # Tfull_to_filter <- Tfull
  
  truevar_file <- paste0("ErrAvgSD_N=", 10000, "_Tfull=", Tfull, "_nreps=", 1, 
                         "_rho=", rho, "_beta=", beta_free, "_omega=", omega_var, ".csv")
  
  # Loading true variance dataframe that contains all necessary info
  if(file.exists(paste0(getwd(),"/GitHub/MasterThesis/True_Variance_Estimate_Collection")) & 
     file.exists(paste0(getwd(),"/GitHub/MasterThesis/True_Variance_Estimate_Collection/", truevar_file))){
    true_var = read.csv(file = paste0(getwd(),"/GitHub/MasterThesis/True_Variance_Estimate_Collection/", truevar_file))
    N_to_load = 10000
    true_var = true_var %>% select(beta_sd, omega_sd, rho_sd)
    true_var = true_var^2 # because info is SD so far
    names(true_var) = c("beta_true_var", "omega_true_var", "rho_true_var")
  } else {
    msg = paste0("Error: TrueVariance DF file path inaccurate and/or file not found") 
    print(msg)
    cat("\n")
  }
  
  
  # # Extract the correct info out of true var df
  # true_var = dplyr::filter(true_var, rho == rho_to_filter,
  #                          beta_coef == beta_to_filter,
  #                          Tfull == Tfull_to_filter,
  #                          omega_var == omega_to_filter)
  
  adjust_var_to_N <- function(unadjusted_true_var, N_new, N_old){
    return((N_old*unadjusted_true_var)/N_new)
  }
  
  adjusted_true_var <- adjust_var_to_N(true_var, N_new = N, N_old = N_to_load)
  alpha = 0.05
  
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
  
  results_table <- data.frame(
    Param = c("beta", "omega", "rho"),
    n_reps = c(rep.int(n_reps, 3)),
    Success_Count = c(success_beta, success_omega, success_rho),
    Success_rate = c(success_beta, success_omega, success_rho)/length(beta_estimates),
    CI_Lower = ci_lower,
    True_param = true_parameters,
    CI_Upper = ci_upper,
    alpha = alpha,
    N = N,
    Tfull = Tfull,
    rho = rho,
    beta = beta_free,
    omega_var = omega_var,
    row.names = NULL
  )
  
  list_of_dataframes[[i]] <- results_table
}

combined_df <- do.call(rbind, list_of_dataframes)
View(combined_df)
write.csv(combined_df, file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/CI_Overlap_Results.csv"), row.names = F)






