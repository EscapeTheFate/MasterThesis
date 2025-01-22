# -------------------------------------------------------------------------
# Master thesis - Code Space for Data Generating Process (DGP) of
# MNP-based Discrete Choice Models
# -------------------------------------------------------------------------
# Preliminary settings for simulative data set creation
J = 5 # Number of alternatives (choice set, could also depend on N & T)
N = 100 # Number of individuals
Time = 20 # Number of choice occasions

# Other:
library(MASS)
library(Rprobit)

# VAR(1)-Error component creation -----------------------------------------

vec2matDimJ <- function(vector){
  J <- get("J", envir = parent.frame()) # a bit dangerous, as it retrieves J from parent-env.
  return(Phi=matrix(vector, nrow = J, byrow = T))
}

get_stationary_error <- function(Sigma, mean = c(rep.int(0,J)), diag = F){
  J <- get("J", envir = parent.frame()) # a bit dangerous, as it retrieves J from parent-env.
  if(diag == T & length(Sigma) == J){
    error_t_zero = MASS::mvrnorm(n = 1, mu = mean, Sigma = diag(Sigma))
    return(error_t_zero)
  }
  if(diag == F & dim(Sigma)[1] == J & dim(Sigma)[2] == J){
    error_t_zero = MASS::mvrnorm(n = 1, mu = mean, Sigma = Sigma)
    return(error_t_zero)
  }
  else{
    print("Error in specification (=get_stationary)")
  }
}

get_VAR_error <- function(time, Phi, Sigma_tilde, previous_error, Sigma_stationary, C = c(rep.int(0,J)), diag = F){
  if(time == 1){
    epsilon_new = get_stationary_error(Sigma = Sigma_stationary, diag = diag)
    return(epsilon_new)
  }
  if(time > 1){
    epsilon_old = previous_error
    epsilon_tilde = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), Sigma = Sigma_tilde)
    epsilon_new = C + Phi %*% epsilon_old + epsilon_tilde
    return(as.vector(epsilon_new))
  }
  else{
    print("Error in specification (=get_VAR_error)")
  }
}

#
# Obtaining simulated probabilities ---------------------------------------

# No random effects for coefficients - No errors
sim_probs_no_rng_no_err <- function(type){
  info_i = 100000
  info = c(rep.int(NA,info_i))
  J = 5
  C = 3
  for (i in 1:info_i){
    # Create ASCs
    asc = c(rep.int(0, J))
    for (j in 1:length(J)){
      asc[j] = C*(j-1)/j
    }
    
    if(type == 1){
      b_1 = c(rep.int(1, J))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
      x_2 = diag(rexp(n = J, rate = 0.5))
      b_3 = c(rep.int(0,J))
      x_3 = diag(b_3)
    }
    if(type == 2){
      b_1 = c(rep.int(1, J))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      b_2 = c(rep.int(0,J))
      x_2 = diag(b_2)
      b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
    
      
    }
    if(type == 3){
      b_1 = c(rep.int(0,J))
      x_1 = diag(b_1)
      b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
      x_2 = diag(rexp(n = J, rate = 0.5))
      b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
    }
    if(type == 4){
      b_1 = c(rep.int(1, J))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
      x_2 = diag(rexp(n = J, rate = 0.5))
      b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
    }
    V = asc + b_1 %*% x_1 + b_2 %*% x_2 + b_3 %*% x_3
    U = V #+ epsilon
    info[i] = which(U == max(U))
    if(i==info_i){ print(table(info)/info_i)}
  }
}
sim_probs_no_rng_no_err(type = 3)

# No random effects for coefficients - IID errors
sim_probs_no_rng_iid_err <- function(type){
  info_i = 100000
  info = c(rep.int(NA,info_i))
  J = 5
  C = 3
  for (i in 1:info_i){
    # Create ASCs
    asc = c(rep.int(0, J))
    for (j in 1:length(J)){
      asc[j] = C*(j-1)/j
    }
    
    if(type == 1){
      b_1 = c(rep.int(1, J))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
      x_2 = diag(rexp(n = J, rate = 0.5))
      b_3 = c(rep.int(0,J))
      x_3 = diag(b_3)
    }
    if(type == 2){
      b_1 = c(rep.int(1, J))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      b_2 = c(rep.int(0,J))
      x_2 = diag(b_2)
      b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
    
      
    }
    if(type == 3){
      b_1 = c(rep.int(0,J))
      x_1 = diag(b_1)
      b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
      x_2 = diag(rexp(n = J, rate = 0.5))
      b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
    }
    if(type == 4){
      b_1 = c(rep.int(1, J))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
      x_2 = diag(rexp(n = J, rate = 0.5))
      b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
    }
    
    # Obtain systematic utility
    V = asc + b_1 %*% x_1 + b_2 %*% x_2 + b_3 %*% x_3
    # Obtain random utility
    epsilon = get_stationary_error(Sigma = c(3,2,9,4,5), diag = T)
    # Put both together
    U = V + epsilon
    # Extract info on final choice by maximum utility entry
    info[i] = which(U == max(U))
    if(i==info_i){ print(table(info)/info_i)}
  }
}
sim_probs_no_rng_iid_err(type = 1)

# No random effects for coefficients - (V)AR errors
sim_probs_no_rng_var_err <- function(type){
  info_i = 10000
  info = c(rep.int(NA_integer_,info_i))
  J = 5
  Time = 10
  C = 3
  rho = 1/3
  for (t in 1:Time){
    for (i in 1:info_i){
      # Create ASCs
      asc = c(rep.int(0, J))
      for (j in 1:length(J)){
        asc[j] = C*(j-1)/j
      }
      
      if(type == 1){
        b_1 = c(rep.int(1, J))
        x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
        b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
        x_2 = diag(rexp(n = J, rate = 0.5))
        b_3 = c(rep.int(0,J))
        x_3 = diag(b_3)
      }
      if(type == 2){
        b_1 = c(rep.int(1, J))
        x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
        b_2 = c(rep.int(0,J))
        x_2 = diag(b_2)
        b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
        x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
        
        
      }
      if(type == 3){
        b_1 = c(rep.int(0,J))
        x_1 = diag(b_1)
        b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
        x_2 = diag(rexp(n = J, rate = 0.5))
        b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
        x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
      }
      if(type == 4){
        b_1 = c(rep.int(1, J))
        x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
        b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
        x_2 = diag(rexp(n = J, rate = 0.5))
        b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
        x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
      }
      
      # Obtain systematic utility
      V = asc + b_1 %*% x_1 + b_2 %*% x_2 + b_3 %*% x_3
      # Obtain random utility
      if(t == 1){
        epsilon = epsilon_tminus = get_VAR_error(time = t, 
                                Phi = rho*diag(1, nrow = J, ncol = J),
                                Sigma_tilde = diag(5, nrow = J, ncol = J),
                                previous_error = epsilon,
                                Sigma_stationary = c(3,2,9,4,5), 
                                diag = T)
      } else {
        epsilon = get_VAR_error(time = t, 
                                Phi = rho*diag(1, nrow = J, ncol = J),
                                Sigma_tilde = diag(5, nrow = J, ncol = J),
                                previous_error = epsilon_tminus,
                                Sigma_stationary = c(3,2,9,4,5), 
                                diag = T)
      }
      # Put both together
      U = V + epsilon
      
      
      # Extract info on final choice by maximum utility entry
      info[i] = which(U == max(U))
      if(i==info_i){
        print(table(info)/info_i)
        print(paste("time: t=", t))
      }
    }
    # Lastly, transfer old epsilon
    epsilon_tminus = epsilon
  }
}
sim_probs_no_rng_var_err(type = 4)

# Random effects for coefficients - No errors
sim_probs_rng_no_err <- function(type){
  info_i = 100000
  info = c(rep.int(NA,info_i))
  J = 5
  C = 3
  for (i in 1:info_i){
    # Create ASCs
    asc = c(rep.int(0, J))
    for (j in 1:length(J)){
      asc[j] = C*(j-1)/j
    }
    
    if(type == 1){
      b_1 = c(rep.int(1, J))
      beta_1 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                             Sigma = vec2matDimJ(c(3,  0,0.5,0.2,0.1,
                                                   0,  2,  0,0.5,0.2,
                                                 0.5,  0,  1,  0,0.5,
                                                 0.2,0.5,  0,  2,  0,
                                                 0.1,0.2,0.5,  0,  3)))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      # ------------------------------------------------------------------------
      b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
      beta_2 = c(rep.int(rnorm(n=1, mean=0, sd = 0.5),J))
      x_2 = diag(rexp(n = J, rate = 0.5))
      # ------------------------------------------------------------------------
      beta_3 = b_3 = c(rep.int(0,J))
      x_3 = diag(b_3)
    }
    if(type == 2){
      b_1 = c(rep.int(1, J))
      beta_1 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                             Sigma = vec2matDimJ(c(3,  0,0.5,0.2,0.1,
                                                   0,  2,  0,0.5,0.2,
                                                   0.5,  0,  1,  0,0.5,
                                                   0.2,0.5,  0,  2,  0,
                                                   0.1,0.2,0.5,  0,  3)))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      # ------------------------------------------------------------------------
      beta_2 = b_2 = c(rep.int(0,J))
      x_2 = diag(b_2)
      # ------------------------------------------------------------------------
      b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
      beta_3 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                             Sigma = vec2matDimJ(c(3,  0,-0.5,-0.2,-0.1,
                                                   0,  2,  0,-0.5,-0.2,
                                                   -0.5,  0,  1,  0,0.5,
                                                   -0.2,-0.5,  0,  2,  0,
                                                   -0.1,-0.2,-0.5,  0,  3)))
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
      
    }
    if(type == 3){
      beta_1 = b_1 = c(rep.int(0,J))
      x_1 = diag(b_1)
      b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
      beta_2 = c(rep.int(rnorm(n=1, mean=0, sd = 0.5),J))
      x_2 = diag(rexp(n = J, rate = 0.5))
      b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
      beta_3 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                             Sigma = vec2matDimJ(c(3,  0,-0.5,-0.2,-0.1,
                                                   0,  2,  0,-0.5,-0.2,
                                                   -0.5,  0,  1,  0,0.5,
                                                   -0.2,-0.5,  0,  2,  0,
                                                   -0.1,-0.2,-0.5,  0,  3)))
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
    }
    if(type == 4){
      b_1 = c(rep.int(1, J))
      beta_1 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                             Sigma = vec2matDimJ(c(3,  0,0.5,0.2,0.1,
                                                   0,  2,  0,0.5,0.2,
                                                   0.5,  0,  1,  0,0.5,
                                                   0.2,0.5,  0,  2,  0,
                                                   0.1,0.2,0.5,  0,  3)))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
      beta_2 = c(rep.int(rnorm(n=1, mean=0, sd = 0.5),J))
      x_2 = diag(rexp(n = J, rate = 0.5))
      b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
      beta_3 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                             Sigma = vec2matDimJ(c(3,  0,-0.5,-0.2,-0.1,
                                                   0,  2,  0,-0.5,-0.2,
                                                   -0.5,  0,  1,  0,0.5,
                                                   -0.2,-0.5,  0,  2,  0,
                                                   -0.1,-0.2,-0.5,  0,  3)))
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
    }
    V = asc + (b_1 + beta_1) %*% x_1 + (b_2 + beta_2) %*% x_2 + (b_3 + beta_3) %*% x_3
    U = V #+ epsilon
    info[i] = which(U == max(U))
    if(i==info_i){ print(table(info)/info_i)}
  }
}
sim_probs_rng_no_err(type = 4)

# Random effects for coefficients - IID errors
sim_probs_rng_iid_err <- function(type){
  info_i = 100000
  info = c(rep.int(NA,info_i))
  J = 5
  C = 3
  for (i in 1:info_i){
    # Create ASCs
    asc = c(rep.int(0, J))
    for (j in 1:length(J)){
      asc[j] = C*(j-1)/j
    }
    
    if(type == 1){
      b_1 = c(rep.int(1, J))
      beta_1 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                             Sigma = vec2matDimJ(c(3,  0,0.5,0.2,0.1,
                                                   0,  2,  0,0.5,0.2,
                                                   0.5,  0,  1,  0,0.5,
                                                   0.2,0.5,  0,  2,  0,
                                                   0.1,0.2,0.5,  0,  3)))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      # ------------------------------------------------------------------------
      b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
      beta_2 = c(rep.int(rnorm(n=1, mean=0, sd = 0.5),J))
      x_2 = diag(rexp(n = J, rate = 0.5))
      # ------------------------------------------------------------------------
      beta_3 = b_3 = c(rep.int(0,J))
      x_3 = diag(b_3)
    }
    if(type == 2){
      b_1 = c(rep.int(1, J))
      beta_1 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                             Sigma = vec2matDimJ(c(3,  0,0.5,0.2,0.1,
                                                   0,  2,  0,0.5,0.2,
                                                   0.5,  0,  1,  0,0.5,
                                                   0.2,0.5,  0,  2,  0,
                                                   0.1,0.2,0.5,  0,  3)))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      # ------------------------------------------------------------------------
      beta_2 = b_2 = c(rep.int(0,J))
      x_2 = diag(b_2)
      # ------------------------------------------------------------------------
      b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
      beta_3 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                             Sigma = vec2matDimJ(c(3,  0,-0.5,-0.2,-0.1,
                                                   0,  2,  0,-0.5,-0.2,
                                                   -0.5,  0,  1,  0,0.5,
                                                   -0.2,-0.5,  0,  2,  0,
                                                   -0.1,-0.2,-0.5,  0,  3)))
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
      
    }
    if(type == 3){
      beta_1 = b_1 = c(rep.int(0,J))
      x_1 = diag(b_1)
      b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
      beta_2 = c(rep.int(rnorm(n=1, mean=0, sd = 0.5),J))
      x_2 = diag(rexp(n = J, rate = 0.5))
      b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
      beta_3 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                             Sigma = vec2matDimJ(c(3,  0,-0.5,-0.2,-0.1,
                                                   0,  2,  0,-0.5,-0.2,
                                                   -0.5,  0,  1,  0,0.5,
                                                   -0.2,-0.5,  0,  2,  0,
                                                   -0.1,-0.2,-0.5,  0,  3)))
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
    }
    if(type == 4){
      b_1 = c(rep.int(1, J))
      beta_1 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                             Sigma = vec2matDimJ(c(3,  0,0.5,0.2,0.1,
                                                   0,  2,  0,0.5,0.2,
                                                   0.5,  0,  1,  0,0.5,
                                                   0.2,0.5,  0,  2,  0,
                                                   0.1,0.2,0.5,  0,  3)))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
      beta_2 = c(rep.int(rnorm(n=1, mean=0, sd = 0.5),J))
      x_2 = diag(rexp(n = J, rate = 0.5))
      b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
      beta_3 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                             Sigma = vec2matDimJ(c(3,  0,-0.5,-0.2,-0.1,
                                                   0,  2,  0,-0.5,-0.2,
                                                   -0.5,  0,  1,  0,0.5,
                                                   -0.2,-0.5,  0,  2,  0,
                                                   -0.1,-0.2,-0.5,  0,  3)))
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
    }
    # Obtain systematic utility
    V = asc + (b_1 + beta_1) %*% x_1 + (b_2 + beta_2) %*% x_2 + (b_3 + beta_3) %*% x_3
    # Obtain random utility
    epsilon = get_stationary_error(Sigma = c(3,2,9,4,5), diag = T)
    # Put both together into perceived utility
    U = V + epsilon
    # Extract info on final choice by maximum utility entry
    info[i] = which(U == max(U))
    if(i==info_i){ print(table(info)/info_i)}
  }
}
sim_probs_rng_iid_err(type = 4)

# No random effects for coefficients - (V)AR errors
sim_probs_rng_var_err <- function(type){
  info_i = 10000
  info = c(rep.int(NA_integer_,info_i))
  J = 5
  Time = 10
  C = 3
  rho = 1/3
  for (t in 1:Time){
    for (i in 1:info_i){
      # Create ASCs
      asc = c(rep.int(0, J))
      for (j in 1:length(J)){
        asc[j] = C*(j-1)/j
      }
      
      if(type == 1){
        b_1 = c(rep.int(1, J))
        beta_1 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                               Sigma = vec2matDimJ(c(3,  0,0.5,0.2,0.1,
                                                     0,  2,  0,0.5,0.2,
                                                     0.5,  0,  1,  0,0.5,
                                                     0.2,0.5,  0,  2,  0,
                                                     0.1,0.2,0.5,  0,  3)))
        x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
        # ------------------------------------------------------------------------
        b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
        beta_2 = c(rep.int(rnorm(n=1, mean=0, sd = 0.5),J))
        x_2 = diag(rexp(n = J, rate = 0.5))
        # ------------------------------------------------------------------------
        beta_3 = b_3 = c(rep.int(0,J))
        x_3 = diag(b_3)
      }
      if(type == 2){
        b_1 = c(rep.int(1, J))
        beta_1 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                               Sigma = vec2matDimJ(c(3,  0,0.5,0.2,0.1,
                                                     0,  2,  0,0.5,0.2,
                                                     0.5,  0,  1,  0,0.5,
                                                     0.2,0.5,  0,  2,  0,
                                                     0.1,0.2,0.5,  0,  3)))
        x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
        # ------------------------------------------------------------------------
        beta_2 = b_2 = c(rep.int(0,J))
        x_2 = diag(b_2)
        # ------------------------------------------------------------------------
        b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
        beta_3 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                               Sigma = vec2matDimJ(c(3,  0,-0.5,-0.2,-0.1,
                                                     0,  2,  0,-0.5,-0.2,
                                                     -0.5,  0,  1,  0,0.5,
                                                     -0.2,-0.5,  0,  2,  0,
                                                     -0.1,-0.2,-0.5,  0,  3)))
        x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
        
      }
      if(type == 3){
        beta_1 = b_1 = c(rep.int(0,J))
        x_1 = diag(b_1)
        b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
        beta_2 = c(rep.int(rnorm(n=1, mean=0, sd = 0.5),J))
        x_2 = diag(rexp(n = J, rate = 0.5))
        b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
        beta_3 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                               Sigma = vec2matDimJ(c(3,  0,-0.5,-0.2,-0.1,
                                                     0,  2,  0,-0.5,-0.2,
                                                     -0.5,  0,  1,  0,0.5,
                                                     -0.2,-0.5,  0,  2,  0,
                                                     -0.1,-0.2,-0.5,  0,  3)))
        x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
      }
      if(type == 4){
        b_1 = c(rep.int(1, J))
        beta_1 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                               Sigma = vec2matDimJ(c(3,  0,0.5,0.2,0.1,
                                                     0,  2,  0,0.5,0.2,
                                                     0.5,  0,  1,  0,0.5,
                                                     0.2,0.5,  0,  2,  0,
                                                     0.1,0.2,0.5,  0,  3)))
        x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
        b_2 = MASS::mvrnorm(1, mu = c(seq(from = 2, to =-2, length.out = J)), Sigma = diag(c(rep.int(0,J))))
        beta_2 = c(rep.int(rnorm(n=1, mean=0, sd = 0.5),J))
        x_2 = diag(rexp(n = J, rate = 0.5))
        b_3 = MASS::mvrnorm(n = 1, mu = c(4,2,1,2,4), Sigma = diag(c(rep.int(0,J))))
        beta_3 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                               Sigma = vec2matDimJ(c(3,  0,-0.5,-0.2,-0.1,
                                                     0,  2,  0,-0.5,-0.2,
                                                     -0.5,  0,  1,  0,0.5,
                                                     -0.2,-0.5,  0,  2,  0,
                                                     -0.1,-0.2,-0.5,  0,  3)))
        x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
      }
      # Obtain systematic utility
      V = asc + (b_1 + beta_1) %*% x_1 + (b_2 + beta_2) %*% x_2 + (b_3 + beta_3) %*% x_3
      # Obtain random utility
      if(t == 1){
        epsilon = epsilon_tminus = get_VAR_error(time = t, 
                                Phi = rho*diag(1, nrow = J, ncol = J),
                                Sigma_tilde = diag(5, nrow = J, ncol = J),
                                previous_error = epsilon,
                                Sigma_stationary = c(3,2,9,4,5), 
                                diag = T)
      } else {
        epsilon = get_VAR_error(time = t, 
                                Phi = rho*diag(1, nrow = J, ncol = J),
                                Sigma_tilde = diag(5, nrow = J, ncol = J),
                                previous_error = epsilon_tminus,
                                Sigma_stationary = c(3,2,9,4,5), 
                                diag = T)
      }
      # Put both together
      U = V + epsilon
      
      
      # Extract info on final choice by maximum utility entry
      info[i] = which(U == max(U))
      if(i==info_i){
        print(table(info)/info_i)
        print(paste("time: t=", t))
      }
    }
    # Lastly, transfer old epsilon
    epsilon_tminus = epsilon
  }
}
sim_probs_rng_var_err(type = 4)

# Processes for balanced Panel --------------------------------------------
# For wide data format - Columns: ID, ID_Choice_ocassion, choice, vars.....

get_individual_i_df <- function(J, T_max, person_i, type, T_properties = list("type" = "balanced",
                                                                           "dist" = "uniform"),
                                beta = T){
  vec2matDimJ <- function(vector){
    J <- get("J", envir = parent.frame()) # a bit dangerous, as it retrieves J from parent-env.
    return(Phi=matrix(vector, nrow = J, byrow = T))
  }
  
  get_stationary_error <- function(Sigma, mean = c(rep.int(0,J)), diag = F){
    J <- get("J", envir = parent.frame()) # a bit dangerous, as it retrieves J from parent-env.
    if(diag == T & length(Sigma) == J){
      error_t_zero = MASS::mvrnorm(n = 1, mu = mean, Sigma = diag(Sigma))
      return(error_t_zero)
    }
    if(diag == F & dim(Sigma)[1] == J & dim(Sigma)[2] == J){
      error_t_zero = MASS::mvrnorm(n = 1, mu = mean, Sigma = Sigma)
      return(error_t_zero)
    }
    else{
      print("Error in specification (=get_stationary)")
    }
  }
  
  get_VAR_error <- function(time, Phi, Sigma_tilde, previous_error, Sigma_stationary, C = c(rep.int(0,J)), diag = F){
    if(time == 1){
      epsilon_new = get_stationary_error(Sigma = Sigma_stationary, diag = diag)
      return(epsilon_new)
    }
    if(time > 1){
      epsilon_old = previous_error
      epsilon_tilde = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), Sigma = Sigma_tilde)
      epsilon_new = C + Phi %*% epsilon_old + epsilon_tilde
      return(as.vector(epsilon_new))
    }
    else{
      print("Error in specification (=get_VAR_error)")
    }
  }
  
  
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
  
  if(type == 1 | type == 2 | type == 4){
    for (v in 1:J) {
      col_name <- paste0("x_", substr(alphabet, v, v))  
      df[[col_name]] <- NA
    }
  }
  if(type == 1 | type == 3 | type == 4){
    for (v in 1:J) {
      col_name <- paste0("y_", substr(alphabet, v, v))  
      df[[col_name]] <- NA  
    }
  }
  if(type == 2 | type == 3 | type == 4){
    for (v in 1:J){
      col_name <- paste0("z_", substr(alphabet, v, v))  
      df[[col_name]] <- NA  
    }
  }
  
  C = 3
  rho = 1/3
  
  # Apply underlying utility choice process
  for (t in 1:T_i){
    
    
    # Create ASCs
    asc = c(rep.int(0, J))
    for (j in 1:length(J)){
      asc[j] = C*(j-1)/j
    }
    
    if(type == 1){
      b_1 = c(rep.int(1, J))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      b_2 = c(seq(from = 2, to =-2, length.out = J))
      x_2 = diag(rep.int(rexp(n = 1, rate = 0.5), J)) # diag(rexp(n = J, rate = 0.5))
      b_3 = c(rep.int(0,J))
      x_3 = diag(b_3)
    }
    if(type == 2){
      b_1 = c(rep.int(1, J))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      b_2 = c(rep.int(0,J))
      x_2 = diag(b_2)
      b_3 = c(4,2,1,2,4)
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
      
    }
    if(type == 3){
      b_1 = c(rep.int(0,J))
      x_1 = diag(b_1)
      b_2 = c(seq(from = 2, to =-2, length.out = J))
      x_2 = diag(rep.int(rexp(n = 1, rate = 0.5), J)) # diag(rexp(n = J, rate = 0.5))
      b_3 = c(4,2,1,2,4)
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
    }
    if(type == 4){
      b_1 = c(rep.int(1, J))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      b_2 = c(seq(from = 2, to =-2, length.out = J))
      x_2 = diag(rep.int(rexp(n = 1, rate = 0.5), J)) # diag(rexp(n = J, rate = 0.5))
      b_3 = c(4,2,1,2,4)
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
    }
    
    if(beta == T){
      if(type == 1 | type == 2 | type == 4){
        # beta_1 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
        #                        Sigma = vec2matDimJ(c(3,  0,0.5,0.2,0.1,
        #                                              0,  2,  0,0.5,0.2,
        #                                              0.5,  0,  1,  0,0.5,
        #                                              0.2,0.5,  0,  2,  0,
        #                                              0.1,0.2,0.5,  0,  3)))
        beta_1 = c(rep.int(0, J))
      }
      if(type == 1 | type == 3 | type == 4){
        beta_2 = c(rep.int(rnorm(n=1, mean=0, sd = 0.5),J))
      }
      if(type == 2 | type == 3 | type == 4){
        beta_3 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                               Sigma = vec2matDimJ(c(3,  0,-0.5,-0.2,-0.1,
                                                     0,  2,  0,-0.5,-0.2,
                                                     -0.5,  0,  1,  0,0.5,
                                                     -0.2,-0.5,  0,  2,  0,
                                                     -0.1,-0.2,-0.5,  0,  3)))
      }
    } else {
      beta_1 = beta_2 = beta_3 = c(rep.int(0,J))
    }
    # Obtain systematic utility
    V = asc + (b_1 + beta_1) %*% x_1 + (b_2 + beta_2) %*% x_2 + (b_3 + beta_3) %*% x_3
    # Obtain random utility
    if(t == 1){
      epsilon = epsilon_tminus = get_VAR_error(time = t, 
                                               Phi = rho*diag(1, nrow = J, ncol = J),
                                               Sigma_tilde = diag(5, nrow = J, ncol = J),
                                               previous_error = epsilon,
                                               Sigma_stationary = c(3,2,9,4,5), 
                                               diag = T)
    } else {
      epsilon = get_VAR_error(time = t, 
                              Phi = rho*diag(1, nrow = J, ncol = J),
                              Sigma_tilde = diag(5, nrow = J, ncol = J),
                              previous_error = epsilon_tminus,
                              Sigma_stationary = c(3,2,9,4,5), 
                              diag = T)
    }
    # Put both together
    U = V + epsilon
    
    # Extract & fill in choice by maximum utility
    choice = which(U == max(U))
    df$choice[t] = substr(alphabet, choice, choice)
    # Fill in variable values by combining x1, x2 and x3 into once vector
    X = c(diag(x_1), diag(x_2), diag(x_3))
    start_off_var_cols = which(names(df) == paste0("x_", substr(alphabet, 1, 1)))
    df[t,start_off_var_cols:dim(df)[2]] = X
  }
  return(df)
}

get_alternative_i_df <- function(J, T_max, person_i, type, T_properties = list("type" = "balanced",
                                                                              "dist" = "uniform"),
                                beta = T){
  vec2matDimJ <- function(vector){
    J <- get("J", envir = parent.frame()) # a bit dangerous, as it retrieves J from parent-env.
    return(Phi=matrix(vector, nrow = J, byrow = T))
  }
  
  get_stationary_error <- function(Sigma, mean = c(rep.int(0,J)), diag = F){
    J <- get("J", envir = parent.frame()) # a bit dangerous, as it retrieves J from parent-env.
    if(diag == T & length(Sigma) == J){
      error_t_zero = MASS::mvrnorm(n = 1, mu = mean, Sigma = diag(Sigma))
      return(error_t_zero)
    }
    if(diag == F & dim(Sigma)[1] == J & dim(Sigma)[2] == J){
      error_t_zero = MASS::mvrnorm(n = 1, mu = mean, Sigma = Sigma)
      return(error_t_zero)
    }
    else{
      print("Error in specification (=get_stationary)")
    }
  }
  
  get_VAR_error <- function(time, Phi, Sigma_tilde, previous_error, Sigma_stationary, C = c(rep.int(0,J)), diag = F){
    if(time == 1){
      epsilon_new = get_stationary_error(Sigma = Sigma_stationary, diag = diag)
      return(epsilon_new)
    }
    if(time > 1){
      epsilon_old = previous_error
      epsilon_tilde = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), Sigma = Sigma_tilde)
      epsilon_new = C + Phi %*% epsilon_old + epsilon_tilde
      return(as.vector(epsilon_new))
    }
    else{
      print("Error in specification (=get_VAR_error)")
    }
  }
  
  
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
  
  if(type == 1 | type == 2 | type == 4){
    for (v in 1:J) {
      col_name <- paste0("x_", substr(alphabet, v, v))  
      df[[col_name]] <- NA
    }
  }
  if(type == 1 | type == 3 | type == 4){
    for (v in 1:J) {
      col_name <- paste0("y_", substr(alphabet, v, v))  
      df[[col_name]] <- NA  
    }
  }
  if(type == 2 | type == 3 | type == 4){
    for (v in 1:J){
      col_name <- paste0("z_", substr(alphabet, v, v))  
      df[[col_name]] <- NA  
    }
  }
  
  C = 3
  rho = 1/3
  
  # Apply underlying utility choice process
  for (t in 1:T_i){
    
    
    # Create ASCs
    asc = c(rep.int(0, J))
    for (j in 1:length(J)){
      asc[j] = C*(j-1)/j
    }
    
    if(type == 1){
      b_1 = c(rep.int(1, J))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      b_2 = c(seq(from = 2, to =-2, length.out = J))
      x_2 = diag(rexp(n = J, rate = 0.5))
      b_3 = c(rep.int(0,J))
      x_3 = diag(b_3)
    }
    if(type == 2){
      b_1 = c(rep.int(1, J))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      b_2 = c(rep.int(0,J))
      x_2 = diag(b_2)
      b_3 = c(4,2,1,2,4)
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
      
    }
    if(type == 3){
      b_1 = c(rep.int(0,J))
      x_1 = diag(b_1)
      b_2 = c(seq(from = 2, to =-2, length.out = J))
      x_2 = diag(rexp(n = J, rate = 0.5))
      b_3 = c(4,2,1,2,4)
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
    }
    if(type == 4){
      b_1 = c(rep.int(1, J))
      x_1 = diag(MASS::mvrnorm(n = 1, mu = c(-2:2), Sigma = diag(c(8,6,3,2.5,2))))
      b_2 = c(seq(from = 2, to =-2, length.out = J))
      x_2 = diag(rexp(n = J, rate = 0.5))
      b_3 = c(4,2,1,2,4)
      x_3 = diag(rlnorm(n = 5, meanlog = 0, sdlog = 2))
    }
    
    if(beta == T){
      if(type == 1 | type == 2 | type == 4){
        beta_1 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)),
                               Sigma = vec2matDimJ(c(3,  0,0.5,0.2,0.1,
                                                     0,  2,  0,0.5,0.2,
                                                     0.5,  0,  1,  0,0.5,
                                                     0.2,0.5,  0,  2,  0,
                                                     0.1,0.2,0.5,  0,  3)))
      }
      if(type == 1 | type == 3 | type == 4){
        beta_2 = c(rep.int(rnorm(n=1, mean=0, sd = 0.5),J))
      }
      if(type == 2 | type == 3 | type == 4){
        beta_3 = MASS::mvrnorm(n = 1, mu = c(rep.int(0,J)), 
                               Sigma = vec2matDimJ(c(3,  0,-0.5,-0.2,-0.1,
                                                     0,  2,  0,-0.5,-0.2,
                                                     -0.5,  0,  1,  0,0.5,
                                                     -0.2,-0.5,  0,  2,  0,
                                                     -0.1,-0.2,-0.5,  0,  3)))
      }
    } else {
      beta_1 = beta_2 = beta_3 = c(rep.int(0,J))
    }
    # Obtain systematic utility
    V = asc + (b_1 + beta_1) %*% x_1 + (b_2 + beta_2) %*% x_2 + (b_3 + beta_3) %*% x_3
    # Obtain random utility
    if(t == 1){
      epsilon = epsilon_tminus = get_VAR_error(time = t, 
                                               Phi = rho*diag(1, nrow = J, ncol = J),
                                               Sigma_tilde = diag(5, nrow = J, ncol = J),
                                               previous_error = epsilon,
                                               Sigma_stationary = c(3,2,9,4,5), 
                                               diag = T)
    } else {
      epsilon = get_VAR_error(time = t, 
                              Phi = rho*diag(1, nrow = J, ncol = J),
                              Sigma_tilde = diag(5, nrow = J, ncol = J),
                              previous_error = epsilon_tminus,
                              Sigma_stationary = c(3,2,9,4,5), 
                              diag = T)
    }
    # Put both together
    U = V + epsilon
    
    # Extract & fill in choice by maximum utility
    choice = which(U == max(U))
    df$choice[t] = substr(alphabet, choice, choice)
    # Fill in variable values by combining x1, x2 and x3 into once vector
    X = c(diag(x_1), diag(x_2), diag(x_3))
    start_off_var_cols = which(names(df) == paste0("x_", substr(alphabet, 1, 1)))
    df[t,start_off_var_cols:dim(df)[2]] = X
  }
  return(df)
}

simulate_panel_df <- function(N, J, T_max, type, T_properties = list("type" = "balanced",
                                                                    "dist" = "uniform"), beta = T){
  for (i in 1:N){
    if(i == 1){
      df <- get_individual_i_df(J, T_max, person_i = i, type, T_properties, beta)
    } else {
      df <- rbind(df, get_individual_i_df(J, T_max, person_i = i, type, T_properties, beta))
    }
  }
  #df$choice <- factor(df$choice, levels = c("A", "B", "C", "D", "E"), labels = c(1, 2, 3, 4, 5))
  df$choice <- as.factor(df$choice)
  return(df)                                                                                   
} 

simulate_alt_panel_df <- function(N, J, T_max, type, T_properties = list("type" = "balanced",
                                                                     "dist" = "uniform"), beta = T){
  for (i in 1:N){
    if(i == 1){
      df <- get_alternative_i_df(J, T_max, person_i = i, type, T_properties, beta)
    } else {
      df <- rbind(df, get_alternative_i_df(J, T_max, person_i = i, type, T_properties, beta))
    }
  }
  #df$choice <- factor(df$choice, levels = c("A", "B", "C", "D", "E"), labels = c(1, 2, 3, 4, 5))
  df$choice <- as.factor(df$choice)
  return(df)                                                                                   
} 


# Testing:
get_individual_i_df(J = 5, T_max = 20, person_i = 99, type = 4, T_properties = list("type" = "unbalanced",
                                                                                  "dist" = "uniform"))
df = simulate_panel_df(N = 500, J = 5, T_max = 1, type = 4, T_properties = list("type" = "balanced",
                                                                                 "dist" = "uniform"), beta = F)
df2 = simulate_alt_panel_df(N = 500, J = 5, T_max = 1, type = 4, T_properties = list("type" = "balanced",
                                                                                "dist" = "uniform"), beta = F)


# Fitting a probit model on dataset ---------------------------------------
library(Rprobit)
data("Train", package = "mlogit")
# choice ~ v1 | v2 | v3, TD: Change x2 to type 2
sim_mod <- setup_Rprobit(form = choice ~ x | y | z, data_raw = df, ids = c("id"))
sim_mod <- setup_Rprobit(form = choice ~ x + y + z | 0 | 0, data_raw = df2, ids = c("id"))
sim_mod <- setup_Rprobit(form = choice ~ 0 | 0 | x + y + z, data_raw = df2, ids = c("id"))
check_identifiability(sim_mod)

mod1 <- mod_cl$new(
  Hb = sim_mod$mod$Hb,
  fb  = sim_mod$mod$fb,
  HO  = sim_mod$mod$HO,
  fO  = sim_mod$mod$fO,
  HL  = sim_mod$mod$HL,
  fL  = sim_mod$mod$fL,
  alt = 5)
check_identifiability(mod1)

sim_fit <- fit_Rprobit(Rprobit_obj = sim_mod)
summary(sim_fit)
# for mixing of coefs of type 2 or type 3: re <- c("V1", "V3") --> sim_mod$mod$lRE changes
print(sim_mod)




# Change scale fixation (for identification purpose)
#simulation_mod$mod$lthL = 1 (is originally 1?)
sim_mod$mod$HL = matrix(0,nrow = J, ncol = 1) 
sim_mod$mod$HL[5,1] = 1 # <- last coef to 1
sim_mod$mod$fL = matrix(0,nrow = J, ncol = 1)

#sim_mod$mod$lthb = 3 (is originally 7?)
sim_mod$mod$Hb = diag(15)[,-1] # diag(4)[,-1]
sim_mod$mod$fb = as.matrix(c(-0.06,0,0,0),ncol=1)

sim_mod$theta_0 <- c()

sim_fit <- fit_Rprobit(Rprobit_obj = sim_mod)
summary(sim_fit)


# -------------------------------------------------------------------------
data("Train", package = "mlogit")
Train$choiceid <- 1:nrow(Train)
train_mod <- setup_Rprobit(form = choice ~ price + time + change + comfort | 0 | 0,  data_raw = Train,  ids = c("id"))
train_fit <- fit_Rprobit(train_mod)
summary(train_fit)

sim_mod <- setup_Rprobit(form = choice ~ x + y + z | 0 | 0, data_raw = df, ids = c("id"))
sim_fit <- fit_Rprobit(sim_mod)
summary(sim_fit)