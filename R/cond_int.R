conditional_interactions_mu <- function(Omega, interactions){
  cond_mu <- c()
  cond_mu_mat <- matrix(0, nrow = ncol(interactions), ncol = ncol(interactions))  
  interactions  <- interactions[lower.tri(interactions)]
  # interactions <- matrix(interactions)
  
  for(i in 1:length(interactions)){
    cond_mu[i] <- Omega[i, -i] %*% solve(Omega[-i,-i]) %*% interactions[-i]
  }
  
  cond_mu_mat[upper.tri(cond_mu_mat, diag = FALSE)] <- cond_mu
  cond_mu_mat[lower.tri(cond_mu_mat, diag = FALSE)] <- cond_mu
  
  
  return(cond_mu_mat)
}


conditional_interactions_sigma <- function(Omega, interactions){
  cond_sigma_sq <- c()
  cond_sigma_mat <- matrix(0, nrow = ncol(interactions), ncol = ncol(interactions))  
  interactions  <- interactions[lower.tri(interactions)]
  
  for(i in 1:length(interactions)){
    cond_sigma_sq[i] <- Omega[i, i] - Omega[i, -i] %*% solve(Omega[-i,-i]) %*%  Omega[-i, i]
  }
  
  cond_sigma <- sqrt(cond_sigma_sq)  #obtain SDs
  cond_sigma_mat[upper.tri(cond_sigma_mat, diag = FALSE)] <- cond_sigma
  cond_sigma_mat[lower.tri(cond_sigma_mat, diag = FALSE)] <- cond_sigma
  
  return(cond_sigma_mat)
}
