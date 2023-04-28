#include <Rcpp.h>
#include "gibbs_graphical_model.h"
#include "utility_functions.h"
using namespace Rcpp;


// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of an edge + interaction
//  pair (using a cauchy prior)
// ----------------------------------------------------------------------------|
List metropolis_edge_interaction_pair_cauchy(NumericMatrix interactions,
                                             NumericMatrix thresholds,
                                             IntegerMatrix gamma,
                                             IntegerMatrix observations,
                                             IntegerVector no_categories,
                                             NumericMatrix proposal_sd,
                                             double cauchy_scale,
                                             IntegerMatrix index,
                                             int no_interactions,
                                             int no_persons,
                                             NumericMatrix rest_matrix,
                                             NumericMatrix inclusion) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  int node1;
  int node2;

  for(int cntr = 0; cntr < no_interactions; cntr ++) {
    node1 = index(cntr, 1) - 1;
    node2 = index(cntr, 2) - 1;

    current_state = interactions(node1, node2);

    if(gamma(node1, node2) == 0) {
      proposed_state = R::rnorm(current_state, proposal_sd(node1, node2));
    } else {
      proposed_state = 0.0;
    }

    log_prob = log_pseudolikelihood_ratio(interactions,
                                          thresholds,
                                          observations,
                                          no_categories,
                                          no_persons,
                                          node1,
                                          node2,
                                          proposed_state,
                                          current_state,
                                          rest_matrix);

    if(gamma(node1, node2) == 0) {
      log_prob += R::dcauchy(proposed_state, 0.0, cauchy_scale, true);
      log_prob -= R::dnorm(proposed_state,
                           current_state,
                           proposal_sd(node1, node2),
                           true);
      log_prob += log(inclusion(node1, node2) /
        (1 - inclusion(node1, node2)));
    } else {
      log_prob -= R::dcauchy(current_state, 0.0, cauchy_scale,  true);
      log_prob += R::dnorm(current_state,
                           proposed_state,
                           proposal_sd(node1, node2),
                           true);
      log_prob -= log(inclusion(node1, node2) /
        (1 - inclusion(node1, node2)));
    }

    U = R::unif_rand();
    if(std::log(U) < log_prob) {
      gamma(node1, node2) = 1 - gamma(node1, node2);
      gamma(node2, node1) = 1 - gamma(node2, node1);

      interactions(node1, node2) = proposed_state;
      interactions(node2, node1) = proposed_state;

      //Update the matrix of rest scores
      double state_diff = proposed_state - current_state;
      for(int person = 0; person < no_persons; person++) {
        rest_matrix(person, node1) += observations(person, node2) *
          state_diff;
        rest_matrix(person, node2) += observations(person, node1) *
          state_diff;
      }
    }
  }
  return List::create(Named("interactions") = interactions,
                      Named("gamma") = gamma,
                      Named("rest_matrix") = rest_matrix);
}



// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of an edge + interaction 
//  pair (using a laplace prior)
// ----------------------------------------------------------------------------|
List metropolis_edge_interaction_pair_laplace(NumericMatrix interactions,
                                              NumericMatrix thresholds,
                                              IntegerMatrix gamma,
                                              IntegerMatrix observations,
                                              IntegerVector no_categories,
                                              NumericMatrix proposal_sd,
                                              double cauchy_scale,
                                              IntegerMatrix index,
                                              int no_interactions,
                                              int no_persons,
                                              NumericMatrix rest_matrix,
                                              NumericMatrix inclusion) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;
  
  int node1;
  int node2;
  
  for(int cntr = 0; cntr < no_interactions; cntr ++) {
    node1 = index(cntr, 1) - 1;
    node2 = index(cntr, 2) - 1;
    
    current_state = interactions(node1, node2);
    
    if(gamma(node1, node2) == 0) {
      proposed_state = R::rnorm(current_state,
                                sd_approx_lap(proposal_sd(node1, node2)));
    } else {
      proposed_state = 0.0;
    }
    
    log_prob = log_pseudolikelihood_ratio(interactions,
                                          thresholds,
                                          observations,
                                          no_categories,
                                          no_persons,
                                          node1,
                                          node2,
                                          proposed_state,
                                          current_state,
                                          rest_matrix);
    
    if(gamma(node1, node2) == 0) {
      log_prob += dlap_1(proposed_state, 0.0, cauchy_scale, true);
      log_prob -= R::dnorm(proposed_state,
                           current_state,
                           sd_approx_lap(proposal_sd(node1, node2)),
                           true);
      log_prob += log(inclusion(node1, node2) /
        (1 - inclusion(node1, node2)));
    } else {
      log_prob -= dlap_1(current_state, 0.0, cauchy_scale,  true);
      log_prob += R::dnorm(current_state,
                           proposed_state,
                           sd_approx_lap(proposal_sd(node1, node2)),
                           true);
      log_prob -= log(inclusion(node1, node2) /
        (1 - inclusion(node1, node2)));
    }
    
    U = R::unif_rand();
    if(std::log(U) < log_prob) {
      gamma(node1, node2) = 1 - gamma(node1, node2);
      gamma(node2, node1) = 1 - gamma(node2, node1);
      
      interactions(node1, node2) = proposed_state;
      interactions(node2, node1) = proposed_state;
      
      //Update the matrix of rest scores
      double state_diff = proposed_state - current_state;
      for(int person = 0; person < no_persons; person++) {
        rest_matrix(person, node1) += observations(person, node2) *
          state_diff;
        rest_matrix(person, node2) += observations(person, node1) *
          state_diff;
      }
    }
  }
  return List::create(Named("interactions") = interactions,
                      Named("gamma") = gamma,
                      Named("rest_matrix") = rest_matrix);
}

// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of an edge + interaction
//  pair (using a unit information prior)
// ----------------------------------------------------------------------------|
List metropolis_edge_interaction_pair_unitinfo(NumericMatrix interactions,
                                               NumericMatrix thresholds,
                                               IntegerMatrix gamma,
                                               IntegerMatrix observations,
                                               IntegerVector no_categories,
                                               NumericMatrix proposal_sd,
                                               NumericMatrix unit_info,
                                               IntegerMatrix index,
                                               int no_interactions,
                                               int no_persons,
                                               NumericMatrix rest_matrix,
                                               NumericMatrix inclusion) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;

  int node1;
  int node2;

  for(int cntr = 0; cntr < no_interactions; cntr ++) {
    node1 = index(cntr, 1) - 1;
    node2 = index(cntr, 2) - 1;

    current_state = interactions(node1, node2);

    if(gamma(node1, node2) == 0) {
      proposed_state = R::rnorm(current_state, proposal_sd(node1, node2));
    } else {
      proposed_state = 0.0;
    }

    log_prob = log_pseudolikelihood_ratio(interactions,
                                          thresholds,
                                          observations,
                                          no_categories,
                                          no_persons,
                                          node1,
                                          node2,
                                          proposed_state,
                                          current_state,
                                          rest_matrix);

    if(gamma(node1, node2) == 0) {
      log_prob += R::dnorm(proposed_state,
                           0.0,
                           unit_info(node1, node2),
                           true);
      log_prob -= R::dnorm(proposed_state,
                           current_state,
                           proposal_sd(node1, node2),
                           true);
      log_prob += log(inclusion(node1, node2) /
        (1 - inclusion(node1, node2)));
    } else {
      log_prob += R::dnorm(current_state,
                           proposed_state,
                           proposal_sd(node1, node2),
                           true);
      log_prob -= R::dnorm(current_state,
                           0.0,
                           unit_info(node1, node2),
                           true);
      log_prob -= log(inclusion(node1, node2) /
        (1 - inclusion(node1, node2)));
    }

    U = R::unif_rand();
    if(std::log(U) < log_prob) {
      gamma(node1, node2) = 1 - gamma(node1, node2);
      gamma(node2, node1) = 1 - gamma(node2, node1);
      interactions(node1, node2) = proposed_state;
      interactions(node2, node1) = proposed_state;

      //Update the matrix of rest scores
      double state_diff = proposed_state - current_state;
      for(int person = 0; person < no_persons; person++) {
        rest_matrix(person, node1) += observations(person, node2) *
          state_diff;
        rest_matrix(person, node2) += observations(person, node1) *
          state_diff;
      }
    }
  }
  return List::create(Named("interactions") = interactions,
                      Named("gamma") = gamma,
                      Named("rest_matrix") = rest_matrix);
}


// ----------------------------------------------------------------------------|
// MH algorithm to sample from the cull-conditional of an edge + interaction
//  pair (using a extended unit information prior)
// ----------------------------------------------------------------------------|
List metropolis_edge_interaction_pair_unitinfo_plus(NumericMatrix interactions,
                                               NumericMatrix thresholds,
                                               IntegerMatrix gamma,
                                               IntegerMatrix observations,
                                               IntegerVector no_categories,
                                               NumericMatrix proposal_sd,
                                               NumericMatrix unit_info,
                                               IntegerMatrix index,
                                               int no_interactions,
                                               int no_persons,
                                               NumericMatrix rest_matrix,
                                               NumericMatrix inclusion) {
  double proposed_state;
  double current_state;
  double log_prob;
  double U;
  
  Function conditional_interactions_mu("conditional_interactions_mu");
  NumericMatrix cond_mu = conditional_interactions_mu(unit_info, interactions);
  Function conditional_interactions_sigma("conditional_interactions_sigma");
  NumericMatrix cond_sigma = conditional_interactions_sigma(unit_info, interactions);
  
  int node1;
  int node2;
  
  for(int cntr = 0; cntr < no_interactions; cntr ++) {
    node1 = index(cntr, 1) - 1;
    node2 = index(cntr, 2) - 1;
    
    current_state = interactions(node1, node2);
    
    if(gamma(node1, node2) == 0) {
      proposed_state = R::rnorm(current_state, proposal_sd(node1, node2));
    } else {
      proposed_state = 0.0;
    }
    
    log_prob = log_pseudolikelihood_ratio(interactions,
                                          thresholds,
                                          observations,
                                          no_categories,
                                          no_persons,
                                          node1,
                                          node2,
                                          proposed_state,
                                          current_state,
                                          rest_matrix);
    
    if(gamma(node1, node2) == 0) {
      log_prob += R::dnorm(proposed_state,
                           cond_mu(node1, node2), 
                           cond_sigma(node1, node2),
                           true);
      log_prob -= R::dnorm(proposed_state,
                           current_state,
                           proposal_sd(node1, node2),
                           true);
      log_prob += log(inclusion(node1, node2) /
        (1 - inclusion(node1, node2)));
    } else {
      log_prob += R::dnorm(current_state,
                           proposed_state,
                           proposal_sd(node1, node2),
                           true);
      log_prob -= R::dnorm(current_state,
                           cond_mu(node1, node2), 
                           cond_sigma(node1, node2),
                           true);
      log_prob -= log(inclusion(node1, node2) /
        (1 - inclusion(node1, node2)));
    }
    
    U = R::unif_rand();
    if(std::log(U) < log_prob) {
      gamma(node1, node2) = 1 - gamma(node1, node2);
      gamma(node2, node1) = 1 - gamma(node2, node1);
      interactions(node1, node2) = proposed_state;
      interactions(node2, node1) = proposed_state;
      
      //Update the matrix of rest scores
      double state_diff = proposed_state - current_state;
      for(int person = 0; person < no_persons; person++) {
        rest_matrix(person, node1) += observations(person, node2) *
          state_diff;
        rest_matrix(person, node2) += observations(person, node1) *
          state_diff;
      }
    }
  }
  return List::create(Named("interactions") = interactions,
                      Named("gamma") = gamma,
                      Named("rest_matrix") = rest_matrix);
}