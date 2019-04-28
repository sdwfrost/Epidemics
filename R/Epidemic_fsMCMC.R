#'
#' Forward Simulation MCMC (fsMCMC, Neal & Huang 2015)
#'
#'

#'
#' Seperate simulation of epidemic into parameters (infection, removal rate) and random elements
#' (event times, event type)
#'
#' Initialise with a simulation such that P(theta; E, U) > 0
#'
#'
#' 1. Propose new theta
#'
#' 2. Calculate probability of observing x* given simulated Y(theta; E, U), P(theta; E, U)
#'
#' 3. Metropolis-Hastings Accept/Reject Step
#'
#' 4. Simulate New E & U
#'
#' 5. Calculate P(theta; E, U)
#'
#' 6. Metropolis-Hastings Accept/Reject Step
#'
#' 7. Back to 1
#'

dHyperGeom = function(x, Y, k, log = TRUE){
  return(sum(mapply(extraDistr::dmvhyper, x, Y, MoreArgs = list(k, log))))
}

#' lambda, initial_infected, T_obs, k
#'
#' @param current_state A list of the current state the Continuous Markov Chain of MCMC
#'                      is at, i.e parameter values and log density value.
#' @param MCMC_parameters A list of any extra parameters which are needed to make proposals
#'                        and/or calculate the log density
fsMCMC_theta_proposal = function(current_state, lambda){

  #' Current State
  theta_curr = current_state[[1]]
  E_curr = current_state[[2]]
  U_curr = current_state[[3]]
  logP_curr = current_state[[4]]
  accept_theta = current_state[[5]]
  #' Define MCMC Parameters
  #lambda = MCMC_parameters[[1]]
  V = current_state[[6]]
  x = current_state[[7]]
  N = current_state[[8]]
  T_obs = current_state[[9]]
  k = current_state[[10]]
  initial_infectives = current_state[[11]]
  no_sampled = current_state[[12]]

  log_theta = log(theta)
  log_theta_prop = log_theta + mvnfast::rmvn(1, mu = c(0,0), sigma = lambda*V)
  theta_prop = exp(log_theta_prop)

  Y_prop = Deterministic_Gillespie1(N, initial_infective, theta_prop[1], theta_prop[2], E_curr, U_curr, T_obs, k, store = FALSE)$panel_data

  logP_prop = dHyperGeom(x, Y_prop, no_sampled, log = TRUE)

  log_a = (logP_prop + sum(theta)) - (logP_curr + sum(theta_prop))

  if(log(runif(1)) < log_a){
    logP_curr = logP_prop
    theta = theta_prop
    accept_theta = accept_theta + 1
  }

  return(list(theta_curr, E_curr, U_curr, logP_curr, accept_theta, V, x, N, T_obs, k,
              initial_infectives, no_sampled))

}


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
  theta_curr = c(beta0, gamma0)
  logP_curr = -Inf

  while(logP_curr == -Inf){
    U_curr = runif(no_draws)
    E_curr = rexp(no_draws)
    if(is.null(kernel)){
      Y_curr = Deterministic_Gillespie1(N, a, theta_curr[1], theta_curr[2], E_curr, U_curr, T_obs, k, store = FALSE)$panel_data
    } else{
      Y_curr = Kernel_Deterministic_Gillespie(N, a, theta_curr[1], theta_curr[2], E_curr, U_curr, T_obs, k, kernel, store = FALSE)$panel_data
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
      Y_prop = Deterministic_Gillespie1(N, a, theta_prop[1], theta_prop[2], E_curr, U_curr, T_obs, k, store = FALSE)$panel_data
    } else{
      Y_prop = Kernel_Deterministic_Gillespie(N, a, theta_prop[1], theta_prop[2], E_curr, U_curr, T_obs, k, kernel, store = FALSE)$panel_data
    }

    logP_prop = dHyperGeom(x, Y_prop, no_sampled, log = TRUE)

    log_a = (logP_prop + sum(theta_curr)) - (logP_curr + sum(theta_prop))

    if(log(runif(1)) < log_a){
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
      Y_prop = Deterministic_Gillespie1(N, a, theta_curr[1], theta_curr[2], E_prop, U_prop, T_obs, k, store = FALSE)$panel_data
    } else{
      Y_prop = Kernel_Deterministic_Gillespie(N, a, theta_curr[1], theta_curr[2], E_prop, U_prop, T_obs, k, kernel, store = FALSE)$panel_data

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




