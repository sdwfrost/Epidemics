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


pi_Y_given_X = function(N, Y, k, T_obs, a, beta, par_gamma, kernel){

  no_sampled = sum(Y[[1]])

  Epi_sim = Epidemic_Simulation(N, a, par_gamma, beta, kernel)

  X = transitions_between_observations(T_obs, k, cbind(Epi_sim$t_inf, Epi_sim$t_rem))

  pi_Y_X = prod(mapply(extraDistr::dmvhyper, Y, X, MoreArgs = list(k = no_sampled)))

  return(pi_Y_X)
}


Pseudo_Marginal_MCMC = function(Y, k, T_obs, kernel, par_beta0, par_gamma0, theta_beta, theta_gamma, lambda, no_sims,
                                no_its, burn_in,
                                lag_max = NA, thinning_factor = 1){

  theta = c(par_beta0, par_gamma0[1])

  #' Intialise Missing Data by Simulating no_sims epidemics with
  #' parameters beta and gamma

  pi_Y_given_theta_hat = mean(replicate(no_sims, pi_Y_given_X(N, Y, k, T_obs, a, par_beta, par_gamma, kernel)))

  for(i in 1:no_its){

    #' Propose new beta and gamma using Multiplicative RW propsal
    log_gamma = log(gamma)
    log_beta = log(beta)
    log_gamma_prop = log_theta + rmvn(1, mu = 0, sigma = V)
    theta_prop = exp(log_theta_prop)

    pi_Y_given_theta_hat_prop = mean(replicate(no_sims, pi_Y_given_X(N, Y, k, T_obs, a, par_beta, par_gamma, kernel)))

    log_u = runif(1)

    log_a = (log(pi_Y_given_theta_hat_prop) + sum(dexp(c(beta_prop, gamma_prop), rate = prior_rate, log = TRUE)) + sum(beta) + sum(gamma)) -
            (log(pi_Y_given_theta_hat) + sum(dexp(c(beta, gamma), rate = prior_rate, log = T)) + sum(beta_prop) + sum(gamma_prop))



  }




}




























