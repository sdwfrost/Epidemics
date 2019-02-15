#'
#' Psuedo-Marginal MCMC For Inference on Epidemic Data
#'


#' Takes multiple discrete observations of an epidemic
#' and the complete observation and evaluates the
#' pdf of observing the partial information given
#' the complete information.

#' function to simulate an epidemic with parameters theta and calculate pi(y|x)
#' @param N, size of complete population
#' @param beta
#' @param gamma
#' @param Y, observed subset of the complete epidemic data
#' @param k, noumber of observations of epidemic
#' @param T_obs, amount of time epidemic is observed for.
#' @param kernel,


pi_Y_given_X = function(N, Y, k, T_obs, beta, par_gamma, kernel){

  no_sampled = sum(Y[[1]])

  Epi_sim = Epidemic_Simulation(N, par_gamma, beta, kernel)

  X = transitions_between_observations(T_obs, k, cbind(Epi_sim$t_inf, Epi_sim$t_rem))

  pi_Y_X = sum(mapply(dmvhyper, Y, X, list(k = no_sampled,  log = TRUE)))

  return(pi_Y_X)
}


Pseudo_Marginal_MCMC = function(Y, k, T_obs, beta0, gamma0, theta_beta, theta_gamma, lambda, no_sims,
                                no_its, burn_in,
                                lag_max = NA, thinning_factor = 1){

  beta = beta0
  gamma = gamma0

  #' Intialise Missing Data by Simulating no_sims epidemics with
  #' parameters beta and gamma

  pi_Y_given_X_hat = replicate(no_sims, )

  for(i in 1:no_its){

    #' Propose new beta and gamma using Multiplicative RW propsal


  }




}




























