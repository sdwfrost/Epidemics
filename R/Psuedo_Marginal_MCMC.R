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

  pi_Y_X = prod(mapply(extraDistr::dmvhyper, Y, X, MoreArgs = list(k = no_sampled)))

  return(pi_Y_X)
}

pi_Y_given_theta_hat = function(no_sims, N, Y, k, T_obs, initial_infective, beta, gamma, parallel = FALSE, no_cores = 2){
  if(parallel){
    pi_Y_X_i = parallel::mcmapply(1:no_sims,
                                 FUN =  function(X) return(pi_Y_given_X(N, Y, k, T_obs, initial_infective, beta, gamma)),
                                 mc.cores = no_cores)

  } else{
    pi_Y_X_i = replicate(no_sims, pi_Y_given_X(N, Y, k, T_obs, initial_infective, beta, gamma))
  }
  return(mean(pi_Y_X_i))
}






Pseudo_Marginal_MCMC = function(N, Y, k, T_obs, beta0, gamma0, prior_rate, initial_infective, lambda, V, no_sims,
                                no_its, burn_in, lag_max = NA, thinning_factor = 1, parallel = FALSE, no_cores){

  Start = as.numeric(Sys.time())

  theta = c(beta0, gamma0)

  #' Variable List to Pass to Cores


  #' Intialise Missing Data by Simulating no_sims epidemics with
  #' parameters beta and gamma
  pi_Y_given_theta_hat_current = pi_Y_given_theta_hat(no_sims, N, Y, k, T_obs, initial_infective, theta[1], theta[2], parallel, no_cores)

  #' Create Storage Matrix
  draws = matrix(NA, nrow = no_its + 1, ncol = length(theta) + 1)
  draws[1,] = c(theta, pi_Y_given_theta_hat_current)

  #' Proposal Acceptance Counter
  accept = 0

  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = no_its)

  for(i in 1:no_its){
    pb$tick()
    #' Propose new beta and gamma using Multiplicative RW propsal
    log_theta = log(theta)
    log_theta_prop = log_theta + mvnfast::rmvn(1, mu = c(0,0), sigma = lambda*V)
    theta_prop = exp(log_theta_prop)
    pi_Y_given_theta_hat_prop = as.numeric(pi_Y_given_theta_hat(no_sims, N, Y, k, T_obs, initial_infective, theta_prop[1], theta_prop[2], parallel, no_cores))

    log_u = log(runif(1))

    log_a = (log(pi_Y_given_theta_hat_prop) + sum(dexp(theta_prop, rate = prior_rate, log = TRUE)) + sum(theta)) -
            (log(pi_Y_given_theta_hat_current) + sum(dexp(theta, rate = prior_rate, log = T)) + sum(theta_prop))

    if(log_u < log_a){
      pi_Y_given_theta_hat_current = pi_Y_given_theta_hat_prop
      theta = theta_prop
      accept = accept + 1
    }
    #' Store State
    draws[i + 1, ] = c(theta, pi_Y_given_theta_hat_current)
  }

  End <- as.numeric(Sys.time())

  time_taken <- End - Start

  # Thin the samples
  draws <- draws[seq(from = burn_in + 1, to = (no_its + 1) - burn_in, by = thinning_factor),]

  # Calculate R_0 sample values
  R0_samples = (N*draws[,1])/draws[,2]

  # Calculate Effective Sample Sizes (and Per Second) and Acceptance Rates
  ESS <- c(coda::effectiveSize(draws[, 1:2]), coda::effectiveSize(R0_samples))
  ESS_sec <- ESS/time_taken
  accept_rate <-  accept/no_its

  # = Plots =
  par(mfrow = c(2,2))

  # Plot Beta Samples and Sample Auto-Corrolation Function
  if(is.na(lag_max)){
    # Beta
    plot(draws[, 1], type = 'l')
    acf(draws[, 1])

    # Gamma
    plot(draws[, 2], type = 'l')
    acf(draws[, 2])
  } else{
    # Beta
    plot(draws[, 1], type = 'l')
    acf(draws[, 1], lag_max)

    # Gamma
    plot(draws[, 2], type = 'l')
    acf(draws[, 2], lag_max)
  }

  #' Calculating Summary Statistics for samples
  beta_summary = c(mean(draws[,1]), sd(draws[,1]))
  gamma_summary = c(mean(draws[,2]), sd(draws[,2]))
  R0_summary = c(mean(R0_samples), sd(R0_samples))

  printed_output(rinf_dist = "Exp", no_proposals = NA, no_its, ESS, time_taken, ESS_sec, accept_rate)
  return(list(draws = draws, accept_rate = accept_rate, ESS = ESS, ESS_sec = ESS_sec, beta_summary = beta_summary,
              gamma_summary = gamma_summary, R0_summary = R0_summary, time_taken = time_taken))
}


























