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

J = 5
info_i = 100000
info = c(rep.int(NA,info_i))
for (i in 1:info_i){
  mu = c(28,-10,24,-2,0)[1:J]
  sigma = diag(c(5,1,2,8,12)[1:J])
  beta = MASS::mvrnorm(1, mu, sigma)
  X = diag(MASS::mvrnorm(1, mu = c(rep.int(0,J)), Sigma = diag(c(rep.int(1,J))))) + matrix(0, J, J)
  epsilon = get_error_t(t=1, sd=c(9,4,7,3,10)[1:J], Phi = create_Phi(vec), test)
  V = X %*% beta
  U = V + epsilon
  info[i] = which(U == max(U))
  if(i==info_i){ print(table(info)/info_i)}
}



# Processes for balanced Panel --------------------------------------------
J = 3
N = 100
T = 20
alphabet <- "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
entries <- strsplit(alphabet, "")[[1]][1:J]
# Initial data frame, on which joins are repeatedly performed on
df1 = data.frame(id = c(rep.int(NA,T)),
                 idc = c(rep.int(NA, T)),
                 choice = c(rep.int(NA, T)))

for (i in 1:N){
  # Iteration data frame that is repeatedly created and filled for each individual
  df = data.frame(id = c(rep.int(NA,T)),
                   idc = c(rep.int(NA, T)),
                   choice = c(rep.int(NA, T)))
  
  # Fill in individual i's dataframe
  df$id = i
  df$idc = 1:T
  for (v in 1:J) {
    col_name <- paste0("x_", substr(alphabet, v, v))  
    df[[col_name]] <- NA  
  }
  
  # Apply underlying utility choice process
  for (t in 1:T){
    
    beta = MASS::mvrnorm(1, mu = c(5,3,8,-2,0)[1:J], sigma = diag(c(5,10,2,8,12)[1:J]))
    X = diag(MASS::mvrnorm(1, mu = c(rep.int(0,J)), Sigma = diag(c(rep.int(1,J))))) + matrix(0, J, J)
    epsilon = get_error_t(t = 1, sd = c(9,4,7,3,10)[1:J], Phi = create_Phi(vec), test)
    V = X %*% beta
    U = V + epsilon
    
    # Extract & fill in choice by maximum utility
    choice = which(U == max(U))
    df$choice[t] = substr(alphabet, choice, choice)
    start_off_var_cols = which(names(df) == paste0("x_", substr(alphabet, 1, 1)))
    df[t,start_off_var_cols:dim(df)[2]] = diag(X)
  }
  
  if (i == 1){
    df1 = df
  } else {
    df1 <- rbind(df1, df)
  }
}


# Fitting a probit model on dataset ---------------------------------------







