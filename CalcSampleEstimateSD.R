
library(tidyverse)
library(dplyr)
error_results = read.csv(file = paste0(getwd(), "/GitHub/MasterThesis/DGP1_500reps_Results.csv"))
vars_of_importance = c("rho", "N", "n_reps", "beta_free", "omega_var", "Tfull")
results = error_results %>% select(all_of(vars_of_importance))
#results = cbind(results, conv_1_2 = NA, conv_3 = NA, conv_4_5 = NA, crash = NA)

failed_rows = c()
beta_estimate_sd = omega_estimate_sd = rho_estimate_sd = c(rep.int(NA, nrow(results)))

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
    beta_estimate_sd[i] = omega_estimate_sd[i] = rho_estimate_sd[i] = NA
    next
  }
  
  beta_estimates <- foo$beta
  omega_estimates <- foo$omega
  rho_estimates <- foo$rho
  
  beta_estimate_sd[i] = sd(beta_estimates)
  omega_estimate_sd[i] = sd(omega_estimates)
  rho_estimate_sd[i] = sd(rho_estimates)
  
}

# Delete previous entries as new one's added are NA in there
error_results = error_results %>% select(-beta_estimate_sd, -rho_estimate_sd, -omega_estimate_sd) 
error_results = error_results %>% mutate(beta_estimate_sd, omega_estimate_sd, rho_estimate_sd)
write.csv(error_results, file = paste0(getwd(), "/GitHub/MasterThesis/DGP1_500reps_Results.csv"), row.names = F)


print(failed_rows)
which(is.na(error_results$beta_estimate_sd))
which(is.na(error_results$omega_estimate_sd))
which(is.na(error_results$rho_estimate_sd))

