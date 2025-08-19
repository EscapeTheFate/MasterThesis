
library(tidyverse)
library(dplyr)

error_results = read.csv(file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ErrorAvgSd_Results.csv"))
results = error_results %>% select(names(error_results)[7:ncol(error_results)])
results = cbind(results, conv_1_2 = NA, conv_3 = NA, conv_4_5 = NA, crash = NA)

failed_rows = c(NULL)

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
  
  results$conv_1_2[i] <- sum(foo$conv_code == 1 | foo$conv_code == 2, na.rm = TRUE)
  results$conv_3[i] <- sum(foo$conv_code == 3, na.rm = TRUE)
  results$conv_4_5[i] <- sum(foo$conv_code == 4 | foo$conv_code == 5, na.rm = TRUE)
  results$crash[i] <- sum(foo$conv_code == 99, na.rm = TRUE)
}

summary(results$conv_1_2)
summary(results$conv_3)
summary(results$conv_4_5)
summary(results$crash)

write.csv(results, file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ConvergenceOverview.csv"), row.names = F)
results = read.csv(file = paste0(getwd(), "/GitHub/MasterThesis/Results_Collection/ConvergenceOverview.csv"))
