# Conditional interactions functions used in the Gibbs functions in rcpp



conditional_interactions_mu <- function(Omega, interactions, observations){
  cond_mu <- c()
  interactions  <- interactions[upper.tri(interactions)]
  
  for(i in 1:length(interactions)){
    cond_mu[i] <- Omega[i, -i] %*% solve(Omega[-i,-i]) %*% interactions[-i]
  }
  
  cond_mu_mat <- matrix(0, nrow = ncol(observations), ncol = ncol(observations))
  cond_mu_mat[upper.tri(cond_mu_mat, diag = FALSE)] <- cond_mu
  cond_mu_mat[lower.tri(cond_mu_mat, diag = FALSE)] <- cond_mu
  
  
  return(cond_mu_mat)
}


conditional_interactions_sigma <- function(Omega, interactions, observations){
  cond_sigma_sq <- c()
  interactions  <- interactions[upper.tri(interactions)]
  
  for(i in 1:length(interactions)){
    cond_sigma_sq[i] <- Omega[i, i] - Omega[i, -i] %*% solve(Omega[-i,-i]) %*%  Omega[-i, i]
  }
  
  cond_sigma <- sqrt(cond_sigma_sq)
  cond_sigma_mat <- matrix(0, nrow = ncol(observations), ncol = ncol(observations))
  cond_sigma_mat[upper.tri(cond_sigma_mat, diag = FALSE)] <- cond_sigma
  cond_sigma_mat[lower.tri(cond_sigma_mat, diag = FALSE)] <- cond_sigma
  
  return(cond_sigma_mat)
}