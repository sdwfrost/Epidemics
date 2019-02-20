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


pi_Y_given_X = function(N, Y, k, T_obs, initial_infective, beta, gamma){

  no_sampled = sum(Y[[1]])

  X = Epidemic_Gillespie(N, initial_infective, gamma, beta, k, T_obs)$panel_data



  print(X)
  pi_Y_X = prod(mapply(extraDistr::dmvhyper, Y, X, MoreArgs = list(k = no_sampled)))

  return(pi_Y_X)
}


Pseudo_Marginal_MCMC = function(N, Y, k, T_obs, beta0, gamma0, prior_rate, initial_infective, lambda, V, no_sims,
                                no_its, burn_in, lag_max = NA, thinning_factor = 1){

  theta = c(beta0, gamma0)

  #' Intialise Missing Data by Simulating no_sims epidemics with
  #' parameters beta and gamma


  pi_Y_given_theta_hat = mean(replicate(no_sims, pi_Y_given_X(N, Y, k, T_obs, initial_infective,
                                                                   theta[1], theta[2])))


  #' Create Storage Matrix
  draws = matrix(NA, nrow = no_its + 1, ncol = length(theta) + 1)
  draws[1,] = c(theta, pi_Y_given_theta_hat)
  accept = 0
  for(i in 1:no_its){

    #' Propose new beta and gamma using Multiplicative RW propsal
    log_theta = log(theta)
    log_theta_prop = log_theta + mvnfast::rmvn(1, mu = c(0,0), sigma = lambda*V)
    theta_prop = exp(log_theta_prop)
    pi_Y_given_theta_hat_prop = mean(replicate(no_sims, pi_Y_given_X(N, Y, k, T_obs, initial_infective,
                                                                     theta_prop[1], theta_prop[2])))

    log_u = runif(1)

    log_a = (log(pi_Y_given_theta_hat_prop) + sum(dexp(theta_prop, rate = prior_rate, log = TRUE)) + sum(theta)) -
            (log(pi_Y_given_theta_hat) + sum(dexp(theta, rate = prior_rate, log = T)) + sum(theta_prop))

    if(log_a < log_u){
      pi_Y_given_theta_hat = pi_Y_given_theta_hat_prop
      theta = theta_prop
      accept = accept + 1
    }
    #' Store State
    draws[i + 1, ] = c(theta, pi_Y_given_theta_hat)
  }
  return(draws)
}


























