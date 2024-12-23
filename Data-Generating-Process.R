# -------------------------------------------------------------------------
# Master thesis - Code Space for Data Generating Process (DGP) of
# MNP-based Discrete Choice Models
# -------------------------------------------------------------------------
# Preliminary settings for simulative data set creation
J = 3 # Number of alternatives (choice set, could also depend on N & T)
N = 100 # Number of individuals
T = 20 # Number of choice occasions

# Other:
library(MASS)
vec = c(1,0,0,0,1,0,0,0,1)*1/3


# VAR(1)-Error component creation -----------------------------------------

create_Phi <- function(vector){
  return(Phi=matrix(vector, nrow = J))
}

get_error_t <- function(t, sd, Phi, err_tminus, Const = c(rep.int(0, J))){
  
  if(t == 1){
    error_t <- as.matrix(rnorm(n = J, mean = 0, sd = sd))
     return(error_t)
  }
  
  error_t <- Const + Phi %*% err_tminus + rnorm(n = J, mean = 0, sd = sd)
  return(error_t)
}

# For wide data format ----------------------------------------------------
# Columns: ID, ID_Choice_ocassion, choice, vars.....

J = 3
info_i = 100000
info = c(rep.int(NA,info_i))
for (i in 1:info_i){
  asc = c(0, rep.int(0, J-1))
  #b = c(5,3,8,-2,0)[1:J]
  #beta = MASS::mvrnorm(1, mu = c(rep.int(0, J)), Sigma = diag(c(5,10,2,8,12)[1:J]))
  b = c(rep.int(0, J))
  sign = sample(c(1,-1), size = 1, prob = c(0.75, 0.25))
  beta = c(rep.int(sign,3))
  X = diag(MASS::mvrnorm(1, mu = c(rep.int(0,J)), Sigma = diag(c(rep.int(1,J))))) + matrix(0, J, J)
  #epsilon = get_error_t(t = 1, sd = c(9,4,7,3,10)[1:J], Phi = create_Phi(vec), test)
  epsilon = get_error_t(t = 1, sd = c(rep.int(1, J)), Phi = create_Phi(9*3*vec), test)
  V = asc + X %*% beta + X %*% b
  U = V + epsilon
  info[i] = which(U == max(U))
  if(i==info_i){ print(table(info)/info_i)}
}
for (i in 1:info_i){
  asc = c(0, rep.int(0, J-1))
  #b = c(5,3,8,-2,0)[1:J]
  #beta = MASS::mvrnorm(1, mu = c(rep.int(0, J)), Sigma = diag(c(5,10,2,8,12)[1:J]))
  b = c(rep.int(1, J))
  beta = MASS::mvrnorm(1, mu = c(rep.int(0, J)), Sigma = diag(c(rep.int(1,J))))
  X = diag(MASS::mvrnorm(1, mu = c(rep.int(0,J)), Sigma = diag(c(rep.int(1,J))))) + matrix(0, J, J)
  #epsilon = get_error_t(t = 1, sd = c(9,4,7,3,10)[1:J], Phi = create_Phi(vec), test)
  epsilon = get_error_t(t = 1, sd = c(rep.int(1, J)), Phi = create_Phi(9*3*vec), test)
  V = asc + X %*% beta + X %*% b
  U = V + epsilon
  info[i] = which(U == max(U))
  if(i==info_i){ print(table(info)/info_i)}
}


# Processes for balanced Panel --------------------------------------------
get_individual_i_df <- function(N, J, T_max, person_i, T_properties = list("type" = "balanced",
                                                                           "dist" = "uniform"),
                                b = T, beta = F){
  
  # Decide whether dataframe is a balanced panel, and if not what distribution/process is used
  panel_type = T_properties$type
  if(panel_type != "balanced"){
    
    T_distribution = T_properties$dist
    T_i = switch(  
      T_distribution,  
      "uniform"= sample(1:T_max, size = 1),
      "truncated_normal" = round(truncnorm::rtruncnorm(n = 1, a = 1, b = T_max, mean = T_max/2, sd = T_max/6)),
      "equal_beta" = round(1 + (T_max - 1) * rbeta(1, shape1 = 2, shape2 = 2)),
      "reliable_beta" = round(1 + (T_max - 1) * rbeta(1, shape1 = 5, shape2 = 1)),
      FALSE
    )
    if(T_i == F){
      print("Distribution is incorrect - Use one of: uniform, truncated_normal, equal_beta, or reliable_beta")
      stop()
    }
  } else {
    T_i = T_max
  }
  
  # i-th iteration data frame for individual i, with  fixed J and T_i
  df = data.frame(id = c(rep.int(NA,T_i)),
                  idc = c(rep.int(NA, T_i)),
                  choice = c(rep.int(NA, T_i)))
  
  # Fill in individual i's dataframe
  df$id = i = person_i
  df$idc = 1:T_i
  
  alphabet <- "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  entries <- strsplit(alphabet, "")[[1]][1:J]
  
  for (v in 1:J) {
    col_name <- paste0("x_", substr(alphabet, v, v))  
    df[[col_name]] <- NA  
  }
  
  # Apply underlying utility choice process
  for (t in 1:T_i){
    
    asc = c(0, rep.int(0, J-1))
    b = c(5,3,8,-2,0)[1:J]
    beta = MASS::mvrnorm(1, mu = c(rep.int(0, J)), Sigma = diag(c(5,10,2,8,12)[1:J]))
    X = diag(MASS::mvrnorm(1, mu = c(rep.int(0,J)), Sigma = diag(c(rep.int(1,J))))) + matrix(0, J, J)
    epsilon = get_error_t(t = 1, sd = c(9,4,7,3,10)[1:J], Phi = create_Phi(vec), test)
    V = asc + X %*% beta + X %*% b
    U = V + epsilon
    
    # Extract & fill in choice by maximum utility
    choice = which(U == max(U))
    df$choice[t] = substr(alphabet, choice, choice)
    start_off_var_cols = which(names(df) == paste0("x_", substr(alphabet, 1, 1)))
    df[t,start_off_var_cols:dim(df)[2]] = diag(X)
  }
  return(df)
}

simulate_individual_df <- function(N, J, T_max, T_properties = list("type" = "balanced",
                                                                    "dist" = "uniform")){
  for (i in 1:N){
    if(i == 1){
      df <- get_individual_i_df(N, J, T_max, person_i = i, T_properties)
    } else {
      df <- rbind(df, get_individual_i_df(N, J, T_max, person_i = i, T_properties))
    }
  }
 
  return(df)                                                                                   
}  


get_individual_i_df(N = 100, J = 3, T_max = 20, person_i = 2, T_properties = list("type" = "unbalanced",
                                                                                  "dist" = "uniform"))
simulate_individual_df(N = 100, J = 3, T_max = 20, T_properties = list("type" = "balanced",
                                                                       "dist" = "uniform"))

# Fitting a probit model on dataset ---------------------------------------
