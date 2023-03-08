#' Bayesian structure learning in Markov Random Fields of mixed binary and 
#' ordinal variables using MCMC. 
#'
#' The function \code{bgm} explores the joint pseudoposterior distribution of 
#' structures and parameters in a Markov Random Field for mixed binary and 
#' ordinal variables. 
#' 
#' collapsed into other categories after recoding. See \code{reformat_data} for
#' details.
#' 
#' @param iter The number of iterations of the Gibbs sampler. Defaults to 
#' \code{1e4}. This usually gives a good idea. But for good estimates it is 
#' recommended to run the procedure for \code{1e5} iterations. 
#' 
#' @param interaction_prior The prior distribution for the interaction effects. 
#' Currently, two prior densities are implemented: The Unit Information prior
#' \code{caching = TRUE} these terms are cached in an \code{n} by \code{p} 
#' matrix of real numbers. This matrix is computed once and updated only when 
#' its values change. This decreases the run time of the algorithm, but requires 
#' more memory. If \code{caching = FALSE}, the procedure does not cache terms. 
#' Defaults to \code{TRUE}. 
#'
#' @param display_progress Should the function show a progress bar 
#' (\code{display_progress = TRUE})? Or not (\code{display_progress = FALSE})?
#' Defauls to \code{TRUE}.
#' 
#' @param burnin The number of burnin iterations. The output of the Gibbs 
#' sampler is stored after \code{burnin} iterations.   
#' are numeric matrices that contain the (model or structure-averaged) posterior 
#' means (EAP estimates) of the pairwise associations, and category thresholds, 
#' respectively. If \code{save = TRUE}, a list containing the 
#' \code{iter} by \code{p *  (p - 1) / 2} matrices \code{samples.gamma} 
#' and \code{samples.interactions}, and the \code{iter} by 
#' \code{sum(no_categories)} matrix \code{samples.thresholds}. These contain the 
#' parameter states at every iteration of the Gibbs sampler. Column averages 
#' offer the EAP estimates.
#' 
#' @examples 
#' \dontrun{
#'  ##Analyse the Wenchuan dataset
#'    
#'  # Here, we use 1e4 iterations, for an actual analysis please use at least 
#'  # 1e5 iterations.
#'  fit = bgm(x = Wenchuan)
#'
#'   
#'  #------------------------------------------------------------------------------|
#'  # INCLUSION - EDGE WEIGHT PLOT
#'  #------------------------------------------------------------------------------|
#'   
#'  par(mar = c(6, 5, 1, 1))
#'  plot(x = fit$interactions[lower.tri(fit$interactions)], 
#'       y = fit$gamma[lower.tri(fit$gamma)], ylim = c(0, 1), 
#'       xlab = "", ylab = "", axes = FALSE, pch = 21, bg = "gray", cex = 1.3)
#'  abline(h = 0, lty = 2, col = "gray")
#'  abline(h = 1, lty = 2, col = "gray")
#'  abline(h = .5, lty = 2, col = "gray")
#'  mtext("Posterior Inclusion Probability", side = 1, line = 3, cex = 1.7)
#'  mtext("Posterior Mode Edge Weight", side = 2, line = 3, cex = 1.7)
#'  axis(1)
#'  axis(2, las = 1)
#'   
#'   
#'  #------------------------------------------------------------------------------|
#'  # EVIDENCE - EDGE WEIGHT PLOT
#'  #------------------------------------------------------------------------------|
#'  
#'  #The bgms package currently assumes that the prior odds are 1:
#'  prior.odds = 1
#'  posterior.inclusion = fit$gamma[lower.tri(fit$gamma)]
#'  posterior.odds = posterior.inclusion / (1 - posterior.inclusion)
#'  log.bayesfactor = log(posterior.odds / prior.odds)
#'  log.bayesfactor[log.bayesfactor > 5] = 5
#' 
#'  par(mar = c(5, 5, 1, 1) + 0.1)
#'  plot(fit$interactions[lower.tri(fit$interactions)], log.bayesfactor, pch = 21, bg = "#bfbfbf", 
#'       cex = 1.3, axes = FALSE, xlab = "", ylab = "", ylim = c(-5, 5.5),
#'       xlim = c(-0.5, 1.5))
#'  axis(1)
#'  axis(2, las = 1)
#'  abline(h = log(1/10), lwd = 2, col = "#bfbfbf")
#'  abline(h = log(10), lwd = 2, col = "#bfbfbf")
#' 
#'  text(x = 1, y = log(1 / 10), labels = "Evidence for Exclusion", pos = 1,
#'       cex = 1.7)
#'  text(x = 1, y = log(10), labels = "Evidence for Inclusion", pos = 3, cex = 1.7)
#'  text(x = 1, y = 0, labels = "Absence of Evidence", cex = 1.7)
#'  mtext("Log-Inclusion Bayes Factor", side = 2, line = 3, cex = 1.5, las = 0)
#'  mtext("Posterior Mean Interactions ", side = 1, line = 3.7, cex = 1.5, las = 0)
#'  
#'  
#'  #------------------------------------------------------------------------------|
#'  # THE LOCAL MEDIAN PROBABILITY NETWORK
#'  #------------------------------------------------------------------------------|
#'   
#'  tmp = fit$interactions[lower.tri(fit$interactions)]
#'  tmp[posterior.inclusion < 0.5] = 0
#'   
#'  median.prob.model = matrix(0, nrow = ncol(Wenchuan), ncol = ncol(Wenchuan))
#'  median.prob.model[lower.tri(median.prob.model)] = tmp
#'  median.prob.model = median.prob.model + t(median.prob.model)
#'  
#'  rownames(median.prob.model) = colnames(Wenchuan)
#'  colnames(median.prob.model) = colnames(Wenchuan)
#'   
#'  library(qgraph)
#'  qgraph(median.prob.model, 
#'         theme = "TeamFortress", 
#'         maximum = .5,
#'         fade = FALSE,
#'         color = c("#f0ae0e"), vsize = 10, repulsion = .9, 
#'         label.cex = 1.1, label.scale = "FALSE", 
#'         labels = colnames(Wenchuan))
#'  }
bgm = function(x,
                iter = 1e5,
                burnin = 1e3,
                interaction_prior = c("UnitInfo", "Cauchy", "Laplace", "Horseshoe", "UnitInfo+"),
                scale = 2.5,
                threshold_alpha = 1,
                threshold_beta = 1,
                tau = 1,
                prop_rel_edges = 0.5,
                save = FALSE,
                caching = TRUE,
                display_progress = FALSE) {
  
  #Check Gibbs input -----------------------------------------------------------
  if(abs(iter - round(iter)) > sqrt(.Machine$double.eps)) 
    stop("Parameter ``iter'' needs to be a positive integer.")
  if(abs(burnin - round(burnin)) > sqrt(.Machine$double.eps) || burnin < 0) 
    stop("Parameter ``burnin'' needs to be a non-negative integer.")
  
  #Check prior set-up for the interaction parameters ---------------------------
  interaction_prior = match.arg(interaction_prior)
  if(interaction_prior == "Cauchy" | interaction_prior == "Laplace" 
     | interaction_prior == "Horseshoe") {
    if(scale <= 0)
      stop("The scale of the Cauchy prior needs to be positive.")
  }  
  
  #Check prior set-up for the threshold parameters -----------------------------
  if(threshold_alpha <= 0  | !is.finite(threshold_alpha))
    stop("Parameter ``threshold_alpha'' needs to be positive.")
  if(threshold_beta <= 0  | !is.finite(threshold_beta))
    stop("Parameter ``threshold_beta'' needs to be positive.")
  
  #Check data input ------------------------------------------------------------
  if(!inherits(x, what = "matrix"))
    stop("The input x is supposed to be a matrix.")
  
  if(ncol(x) < 2)
    stop("The matrix x should have more than one variable (columns).")
  if(nrow(x) < 2)
    stop("The matrix x should have more than one observation (rows).")
  
  
  #Check HS arguments------------------------------------------------------------
  if (!(tau == 1 | tau == 2 | tau == 3)) 
    stop("You have chosen an invalid option for the computation of the global 
         shrinkage paramater for the Horseshoe slab. Currently, there are three options: 
         1 = Half Normal distribution with a scale defined by the total number of relevant 
         edges, suppiled in the prop_rel_edges argument; 2 = Half Cauchy distribution with 
         a scale defined by the total number of relevant edges, suppiled in the 
         prop_rel_edges argument 3 = Half Cauchy distribution with a scale of 1.")
  if(prop_rel_edges > 1 | prop_rel_edges < 0)
  stop("The supplied value for the argument prop_rel_edges should be between 0 and 1, 
        indicating the proportion of relevant edges that the users a priori believes 
        should be included in the network.")
  
  #Format the data input -------------------------------------------------------
  data = reformat_data(x = x)
  x = data$x
  no_categories = data$no_categories
  
  no_nodes = ncol(x)
  no_interactions = no_nodes * (no_nodes - 1) / 2
  no_thresholds = sum(no_categories)
  
  #Proposal set-up for the interaction parameters ------------------------------
  if(interaction_prior == "Cauchy" | interaction_prior == "Laplace" 
     | interaction_prior == "Horseshoe") {
    pps = try(mppe(x = x, 
                   no_categories = no_categories, 
                   interaction_prior = interaction_prior, 
                   scale = scale,
                   tau = tau,
                   prop_rel_edges = prop_rel_edges), 
              silent = TRUE)
  } else {
    pps = try(mppe(x = x, 
                   no_categories = no_categories, 
                   interaction_prior = interaction_prior), 
              silent = TRUE)
  }
  if(inherits(pps, what = "try-error"))
    stop("We use normal approximations to the posterior as proposal distribution 
  in a Metropolis within Gibbs approach. The normal approximation is based on 
  the curvature around the mode, estimated from optimizing the full
  pseudoposterior. For your data the pseudoposterior could not be optimized. 
  Please check your data for missing categories, or low category counts. Please 
  contact the package author if the data checks out and the data do not explain 
  why optimization failed.")
  
  if(interaction_prior == "UnitInfo") {
    unit_info = sqrt(pps$unit_info)
  } 
  
  else if(interaction_prior == "UnitInfo+") {
    unit_info = pps$unit_info
  } 
  
  else {
    unit_info = matrix(data = NA, nrow = 1, ncol = 1)
  }
  
  #Set up the variance of the (normal) proposal distribution
  hessian = pps$hessian[-c(1:no_thresholds), -c(1:no_thresholds)]
  
  proposal_sd = matrix(0, 
                       nrow = no_nodes,
                       ncol = no_nodes)
  cntr = 0
  for(node1 in 1:(no_nodes - 1)) {
    for(node2 in (node1 + 1):no_nodes) {
      cntr = cntr + 1
      proposal_sd[node1, node2] = sqrt(-1 / hessian[cntr, cntr])
      proposal_sd[node2, node1] = proposal_sd[node1, node2]
    }
  }
  
  # # Starting value of model matrix:
  gamma = matrix(1,
                 nrow = no_nodes,
                 ncol = no_nodes)

  #Starting values of interactions and thresholds (posterior mode)
  interactions = pps$interactions
  thresholds = pps$thresholds
  
  #Precomputing number of observations per category for each node.
  n_cat_obs = matrix(0, 
                     nrow = max(no_categories) + 1, 
                     ncol = no_nodes)
  for(node in 1:no_nodes) {
    for(category in 0:no_categories[node]) {
      n_cat_obs[category + 1, node] = sum(x[, node] == category)
    }
  }
  
  # Index vector used to sample interactions in a random order.
  Index = matrix(0, 
                 nrow = no_nodes * (no_nodes - 1) / 2, 
                 ncol = 3)
  cntr = 0
  for(node1 in 1:(no_nodes - 1)) {
    for(node2 in (node1 + 1):no_nodes) {
      cntr =  cntr + 1
      Index[cntr, 1] = cntr
      Index[cntr, 2] = node1
      Index[cntr, 3] = node2
    }
  }
  
  #The Metropolis within Gibbs sampler -----------------------------------------
  out = gibbs_sampler(observations = x,
                      gamma = gamma,
                      interactions = interactions,
                      thresholds = thresholds,
                      no_categories  = no_categories,
                      interaction_prior = interaction_prior,
                      scale = scale,
                      tau,
                      prop_rel_edges = prop_rel_edges,
                      unit_info = unit_info,
                      proposal_sd = proposal_sd,
                      Index = Index,
                      iter = iter,
                      burnin = burnin,
                      n_cat_obs = n_cat_obs, 
                      threshold_alpha = threshold_alpha,
                      threshold_beta = threshold_beta,
                      save = save,
                      caching = caching,
                      display_progress = display_progress)
  
  #Preparing the output --------------------------------------------------------
  if(save == FALSE) {
    gamma = out$gamma
    interactions = out$interactions
    tresholds = out$thresholds
    
    colnames(interactions) = paste0("node ", 1:no_nodes)
    rownames(interactions) = paste0("node ", 1:no_nodes)
    colnames(gamma) = paste0("node ", 1:no_nodes)
    rownames(gamma) = paste0("node ", 1:no_nodes)
    colnames(tresholds) = paste0("category ", 1:max(no_categories))
    rownames(tresholds) = paste0("node ", 1:no_nodes)
    
    return(list(gamma = gamma, 
                interactions = interactions,
                thresholds = tresholds))
  } else {
    gamma = out$gamma
    interactions = out$interactions
    thresholds = out$thresholds
    
    names1 = names2 = character(length = no_nodes * (no_nodes - 1) / 2)
    cntr = 0
    for(node in 1:(no_nodes - 1)) {
      for(node_2 in (node + 1):no_nodes) {
        cntr = cntr + 1
        names1[cntr] = paste0("sigma(",node, ", ",node_2,")")
        names2[cntr] = paste0("gamma(",node, ", ",node_2,")")
      }
    }
    colnames(gamma) = names2
    colnames(interactions) = names1
    
    names = character(length = sum(no_categories))
    cntr = 0
    for(node in 1:no_nodes) {
      for(category in 1:no_categories[node]) {
        cntr = cntr + 1
        names[cntr] = paste0("threshold(",node, ", ",category,")")
      }
    }
    colnames(thresholds) = names
    
    rownames(gamma) = paste0("Iter. ", 1:iter)
    rownames(interactions) = paste0("Iter. ", 1:iter)
    rownames(thresholds) = paste0("Iter. ", 1:iter)
    
    return(list(gamma = gamma, 
                interactions = interactions,
                thresholds = thresholds))
  }
}