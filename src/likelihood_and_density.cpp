#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double log_pseudolikelihood(NumericMatrix interactions,
                            NumericMatrix thresholds,
                            IntegerMatrix observations,
                            IntegerVector no_categories) {
  int no_nodes = observations.ncol();
  int no_persons = observations.nrow();
  double rest_score = 0.0;
  double bound =  0.0;
  double log_pseudolikelihood = 0.0;
  double denominator = 0.0;
  double exponent = 0.0;
  int score = 0;
  
  //Contributions of the full-conditional of nodes (pseudolikelihoods) ---------
  for(int s = 0; s <  no_nodes; s++) {
    //Numerator of full-conditional of node s (pseudolikelihood) ---------------
    for(int person = 0; person < no_persons; person++) {
      rest_score = 0.0;
      for(int node = 0; node < no_nodes; node++) {
        rest_score += observations(person, node) * interactions(node, s);  
      }
      log_pseudolikelihood += observations(person,s) * rest_score;
      bound = no_categories[s] * rest_score;
      log_pseudolikelihood -= bound;
      denominator = std::exp(-bound);
      for(int category = 0; category < no_categories[s]; category++) {
        if(observations(person, s) == category + 1) {
          log_pseudolikelihood += thresholds(s, category);
        }
        score = category + 1;
        exponent = thresholds(s, category) + 
          score * rest_score - 
          bound;
        denominator += std::exp(exponent);
      }
      //Denominator of full-conditional of node s (pseudolikelihood) -----------
      log_pseudolikelihood -= log(denominator);
    }
  }
  return log_pseudolikelihood;
}

// [[Rcpp::export]]
double log_unnormalized_pseudoposterior_normal(NumericMatrix interactions,
                                               NumericMatrix thresholds,
                                               IntegerMatrix observations,
                                               IntegerVector no_categories,
                                               NumericMatrix interaction_var,
                                               double threshold_alpha = 1.0,
                                               double threshold_beta = 1.0) {
  int no_nodes = observations.ncol();

  double unn_pseudo_post = log_pseudolikelihood(interactions,
                                                thresholds,
                                                observations,
                                                no_categories);
  
  //Contribution of the prior densities (interactions) -------------------------
  for(int s = 0; s < no_nodes - 1; s++) {
    for(int t = s + 1; t < no_nodes; t++) {
        unn_pseudo_post += R::dnorm(interactions(s, t), 
                            0.0, 
                            std::sqrt(interaction_var(s, t)), 
                            true);
    }
  }
  
  //Contribution of the prior densities (thresholds) ---------------------------
  for(int s = 0; s < no_nodes; s++) {
    for(int category = 0; category < no_categories[s]; category++) {
      unn_pseudo_post -= R::lbeta(threshold_alpha, threshold_beta);
      unn_pseudo_post += threshold_alpha * thresholds(s, category);
      unn_pseudo_post -= (threshold_alpha + threshold_beta) * 
        std::log(1 + std::exp(thresholds(s, category)));
    }
  }
  return unn_pseudo_post;
}

// [[Rcpp::export]]
double log_unnormalized_pseudoposterior_cauchy(NumericMatrix interactions,
                                               NumericMatrix thresholds,
                                               IntegerMatrix observations,
                                               double scale,
                                               IntegerVector no_categories,
                                               double threshold_alpha = 1.0,
                                               double threshold_beta = 1.0) {
  int no_nodes = observations.ncol();
  double unn_pseudo_post = log_pseudolikelihood(interactions,
                                                thresholds,
                                                observations,
                                                no_categories);
  
  //Contribution of the prior densities (interactions) -------------------------
  for(int s = 0; s < no_nodes - 1; s++) {
    for(int t = s + 1; t < no_nodes; t++) {
      unn_pseudo_post += R::dcauchy(interactions(s, t), 
                            0.0, 
                            scale, 
                            true);
    }
  }
  
  //Contribution of the prior densities (thresholds) ---------------------------
  for(int s = 0; s < no_nodes; s++) {
    for(int category = 0; category < no_categories[s]; category++) {
      unn_pseudo_post -= R::lbeta(threshold_alpha, threshold_beta);
      unn_pseudo_post += threshold_alpha * thresholds(s, category);
      unn_pseudo_post -= (threshold_alpha + threshold_beta) * 
        std::log(1 + std::exp(thresholds(s, category)));
    }
  }
  return unn_pseudo_post;
}


// [[Rcpp::export]]
double dlap(double interaction,
            double mu = 0.0,
            double b = 1.0,
            bool log = true) {
  
  
  double logden = - std::log(2) - std::log(b) - std::abs(interaction-mu)/b;
  
  if(log == true){
    return(logden);
  }
  
  else{
    return(std::exp(logden));
  }
}


// [[Rcpp::export]]
double log_unnormalized_pseudoposterior_laplace(NumericMatrix interactions,
                                            NumericMatrix thresholds,
                                            IntegerMatrix observations,
                                            double scale,
                                            IntegerVector no_categories,
                                            double threshold_alpha = 1.0,
                                            double threshold_beta = 1.0) {
  int no_nodes = observations.ncol();
  double unn_pseudo_post = log_pseudolikelihood(interactions,
                                                thresholds,
                                                observations,
                                                no_categories);
  
  //Contribution of the prior densities (interactions) -------------------------
  for(int s = 0; s < no_nodes - 1; s++) {
    for(int t = s + 1; t < no_nodes; t++) {
      unn_pseudo_post += dlap(interactions(s, t), 
                              0.0, 
                              scale, 
                              true);
    }
  }
  
  //Contribution of the prior densities (thresholds) ---------------------------
  for(int s = 0; s < no_nodes; s++) {
    for(int category = 0; category < no_categories[s]; category++) {
      unn_pseudo_post -= R::lbeta(threshold_alpha, threshold_beta);
      unn_pseudo_post += threshold_alpha * thresholds(s, category);
      unn_pseudo_post -= (threshold_alpha + threshold_beta) * 
        std::log(1 + std::exp(thresholds(s, category)));
    }
  }
  return unn_pseudo_post;
}


// [[Rcpp::export]]
double dh(double interaction,
            double scale = 1, 
            double tau = 1, 
            double prop_rel_edges = 1,         // it doesn't allow me not to give it default arguments 
            bool log = false, 
            IntegerVector no_categories = 1,
            double no_persons = 1, 
            double no_interactions = 1) {
  
  double lambda = abs(R::rcauchy(0,1));   //I made the scale to always equal 1
  
  double rel_edges =  prop_rel_edges*no_interactions;
  
  int max_no_categories = max(no_categories);
  double tau_0 = (rel_edges / (no_interactions - rel_edges)) * (max_no_categories / sqrt(no_persons));
  
  if (tau == 1) {
    tau = abs(R::rnorm(0, tau_0));
  } else if (tau == 2) {
    tau = abs(R::rcauchy(0, tau_0));
  } else if (tau == 3) {
    tau = abs(R::rcauchy(0,1));
  }
  
  double sd =  lambda* tau;
  double logden = R::dnorm(interaction, 0, sd, log = log);
  return logden;
}

// [[Rcpp::export]]
double log_unnormalized_pseudoposterior_horseshoe(NumericMatrix interactions,
                                                  NumericMatrix thresholds,
                                                  IntegerMatrix observations,
                                                  double scale,
                                                  double tau,
                                                  double prop_rel_edges,
                                                  IntegerVector no_categories,
                                                  double threshold_alpha = 1.0,
                                                  double threshold_beta = 1.0,
                                                  double no_persons = 1, 
                                                  double no_interactions = 1) {
  int no_nodes = observations.ncol();
  double unn_pseudo_post = log_pseudolikelihood(interactions,
                                                thresholds,
                                                observations,
                                                no_categories);
  
  //Contribution of the prior densities (interactions) -------------------------
  for(int s = 0; s < no_nodes - 1; s++) {
    for(int t = s + 1; t < no_nodes; t++) {
      unn_pseudo_post += dh(interactions(s, t),
                            scale,
                            tau,
                            prop_rel_edges, 
                            true, 
                            no_categories, 
                            no_persons, 
                            no_interactions);
    }
  }
  
  //Contribution of the prior densities (thresholds) ---------------------------
  for(int s = 0; s < no_nodes; s++) {
    for(int category = 0; category < no_categories[s]; category++) {
      unn_pseudo_post -= R::lbeta(threshold_alpha, threshold_beta);
      unn_pseudo_post += threshold_alpha * thresholds(s, category);
      unn_pseudo_post -= (threshold_alpha + threshold_beta) * 
        std::log(1 + std::exp(thresholds(s, category)));
    }
  }
  return unn_pseudo_post;
}

// [[Rcpp::export]]
double emvs_log_unnormalized_pseudoposterior(NumericMatrix interactions,
                                             NumericMatrix thresholds,
                                             IntegerMatrix observations,
                                             IntegerVector no_categories,
                                             double xi,
                                             NumericMatrix slab_var,
                                             double theta = 0.5,
                                             bool hierarchical = false,
                                             double indicator_alpha = 1.0,
                                             double indicator_beta = 1.0,
                                             double threshold_alpha = 1.0,
                                             double threshold_beta = 1.0) {
  int no_nodes = observations.ncol();
  int no_persons= observations.nrow();
  double unn_pseudo_post = 0.0;

  unn_pseudo_post = log_pseudolikelihood(interactions,
                                         thresholds,
                                         observations,
                                         no_categories);
  
  //Contribution of the prior densities (interactions) -------------------------
  for(int s = 0; s < no_nodes - 1; s++) {
    for(int t = s + 1; t < no_nodes; t++) {
      unn_pseudo_post += std::log(
        theta * R::dnorm(interactions(s, t), 
                         0.0, 
                         std::sqrt(slab_var(s, t)), 
                         false) + 
                           (1 - theta) * R::dnorm(interactions(s, t), 
                            0.0, 
                            std::sqrt(slab_var(s, t) * xi / no_persons), 
                            false));
    }
  }
  
  //Contribution of the prior densities (thresholds) ---------------------------
  for(int s = 0; s < no_nodes; s++) {
    for(int category = 0; category < no_categories[s]; category++) {
      unn_pseudo_post -= R::lbeta(threshold_alpha, threshold_beta);
      unn_pseudo_post += threshold_alpha * thresholds(s, category);
      unn_pseudo_post -= (threshold_alpha + threshold_beta) * 
        std::log(1 + std::exp(thresholds(s, category)));
    }
  }
  
  if(hierarchical == true)
    unn_pseudo_post += R::dbeta(theta, indicator_alpha, indicator_beta, true);
  
  return unn_pseudo_post;
}