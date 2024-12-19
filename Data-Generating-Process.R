# -------------------------------------------------------------------------
# Master thesis - Code Space for Data Generating Process (DGP) of
# MNP-based Discrete Choice Models
# -------------------------------------------------------------------------
# Preliminary settings for simulative data set creation
J = 3 # Number of alternatives
N = 100 # Number of individuals
T = 20 # Number of choice occasions

# Other:
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

# # Test for t=1
# get_error_t(t=1, sd=c(1,1,1), Phi = create_Phi(vec), test)
# # Test for t in {1,..., 10}
# for (t in 1:10){
#   err <- get_error_t(t, sd=c(2,3,4), Phi=create_Phi(vec), old_err, Const = c(100,100,100))
#   err <- as.matrix(err)
#   print(err)
#   old_err <- err
# }

