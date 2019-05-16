#'
#' Endemic (SIS) fsMCMC
#'
#'



#' Epidemic_fsMCMC MCMC algorithm to make inference on panel data of an epidemic through forward simulation
#' @param N Total size of closed population
#' @param initial_infective Number of individuals who are initially infected in the population
#' @param x panel data observed. Follows a random sample of n individuals and observes them at k timepoints
#' @param beta0 beta starting value (Infection Rate Parameter)
#' @param gamma0 gamma starting value (Removal rate parameter)
#' @param no_draws How many rexp(1) and runif(1) draws to make for use in the Gillespie algorithm.
#' @param T_obs the period for which the epidemic is observed.
#' @param k how many equally spaced observations take place
#' @param lambda RWM proposal parameter
#' @param V RWM proposal covariance matrix
#' @param no_its The number of MCMC iterations
#' @param burn_in How many of the MCMC iterations are thrown away as burn in (Convergence to Stationary Distn)
#' @param lag_max When plotting the estimated ACF of samples, what will be the maximum lag estimated/plotted
#' @param thinning_factor
#'
Epidemic_fsMCMC = function(N, a, x, beta0, gamma0, kernel = NULL, no_draws, s, T_obs, k, lambda, V, no_its,
                           burn_in = 0, lag_max = NA, thinning_factor = 1){

  Start = as.numeric(Sys.time())

  no_sampled = sum(x[[1]])

  beta_indices = 1:length(beta0)
  gamma_indices = (length(beta) + 1):(length(beta0) + length(gamma0))
  theta_curr = c(beta0, gamma0)
  logP_curr = -Inf

  while(logP_curr == -Inf){
    U_curr = runif(no_draws)
    E_curr = rexp(no_draws)
    if(is.null(kernel)){
      Y_curr = Deterministic_Gillespie1(N, a, theta_curr[beta_indices], theta_curr[gamma_indices], E_curr, U_curr, T_obs, k, store = FALSE)$panel_data
    } else{
      Y_curr = SIS_Deterministic_Gillespie(N, a, theta_curr[beta_indices], theta_curr[gamma_indices], E_curr, U_curr, T_obs, k, kernel,
                                           obs_end = T_obs[2], store = FALSE)$panel_data
    }
    logP_curr = dHyperGeom(x, Y_curr, no_sampled, log = TRUE)
  }

  #' Create Storage Matrix
  draws = matrix(NA, nrow = no_its + 1, ncol = length(theta_curr) + 1)
  draws[1,] = c(theta_curr, logP_curr)

  #' Proposal Acceptance Counter
  accept_theta = 0
  accept_RVs = 0
  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = no_its)

  for(i in 1:no_its){
    pb$tick()
    #' Propose new beta and gamma using Multiplicative RW propsal
    log_theta_curr = log(theta_curr)
    log_theta_prop = log_theta_curr + mvnfast::rmvn(1, mu = c(0,0), sigma = lambda*V)
    theta_prop = exp(log_theta_prop)

    if(is.null(kernel)){
      Y_prop = Deterministic_Gillespie1(N, a, theta_prop[beta_indices], theta_prop[gamma_indices], E_curr, U_curr, T_obs, k, store = FALSE)$panel_data
    } else{
      Y_prop = SIS_Deterministic_Gillespie(N, a, theta_prop[beta_indices], theta_prop[gamma_indices], E_curr, U_curr, T_obs, k, kernel, T_obs[2], store = FALSE)$panel_data
    }

    logP_prop = dHyperGeom(x, Y_prop, no_sampled, log = TRUE)

    log_a = (logP_prop + sum(theta_curr)) - (logP_curr + sum(theta_prop))

    log_u = log(runif(1))

    if(log_u < log_a){
      logP_curr = logP_prop
      theta_curr = theta_prop
      accept_theta = accept_theta + 1
    }

    #' Draw New Random Variables
    E_prop = E_curr
    U_prop = U_curr

    proposal_set = sample(1:no_draws, size = s)

    E_prop[proposal_set] = rexp(s)
    U_prop[proposal_set] = runif(s)

    if(is.null(kernel)){
      Y_prop = Deterministic_Gillespie1(N, a, theta_curr[beta_indices], theta_curr[gamma_indices], E_prop, U_prop, T_obs, k, store = FALSE)$panel_data
    } else{
      Y_prop = SIS_Deterministic_Gillespie(N, a, theta_curr[beta_indices], theta_curr[gamma_indices], E_prop, U_prop, T_obs, k, kernel, T_obs[2], store = FALSE)$panel_data
    }

    logP_prop = dHyperGeom(x, Y_prop, no_sampled, log = TRUE)


    log_u = log(runif(1))

    log_a = (logP_prop) - (logP_curr)

    if(log_u < log_a){
      logP_curr = logP_prop
      E_curr = E_prop
      U_curr = U_prop
      accept_RVs = accept_RVs + 1
    }


    #' Store State
    draws[i+1, ] = c(theta_curr, logP_curr)
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
  accept_rate <-  c(accept_theta, accept_RVs)/no_its

  # = Plots =
  par(mfrow = c(2,2))

  # Plot Beta Samples and Sample Auto-Corrolation Function
  if(is.na(lag_max)){
    # Beta
    plot(draws[, 1], type = 'l', ylab = expression(beta))
    acf(draws[, 1], main = "")

    # Gamma
    plot(draws[, 2], type = 'l', ylab = expression(gamma))
    acf(draws[, 2], main = "")
  } else{
    # Beta
    plot(draws[, 1], type = 'l', ylab = expression(beta))
    acf(draws[, 1], lag_max, main = "")

    # Gamma
    plot(draws[, 2], type = 'l',  ylab = expression(gamma))
    acf(draws[, 2], lag_max, main = "")
  }

  #' Calculating Summary Statistics for samples
  beta_summary = c(mean(draws[,1]), sd(draws[,1]))
  gamma_summary = c(mean(draws[,2]), sd(draws[,2]))
  R0_summary = c(mean(R0_samples), sd(R0_samples))

  printed_output(rinf_dist = "Exp", no_proposals = NA, no_its, ESS, time_taken, ESS_sec, accept_rate)
  return(list(draws = draws, accept_rate = accept_rate, ESS = ESS, ESS_sec = ESS_sec, beta_summary = beta_summary,
              gamma_summary = gamma_summary, R0_summary = R0_summary, time_taken = time_taken))
}
