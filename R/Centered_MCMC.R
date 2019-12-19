#'
#'
#' PhD
#' A Gibbs Sampler for a partially observed SIR epidemic.
#'
#' A Centered Parameterisation of Bayesian Modelling of the Infection and Removal times
#' of an Epidemic Centered
#'
#' ==== Model ====
#'
#' theta(beta, gamma) --> Infection Times (I, Unobserved) --> Removal Times (R, Observed)
#'
#'
#' Approach is based on that used by Neal & Roberts (2003).
#'
#'

# ==== Preamble ====

#' 1. Take observed data and simulate a set of valid infection times using inital values
#'    for gamma and beta.
#'
#'    Steps for one iteration:
#'
#' 2. Calculate the parameters of the conditional posteriors
#'    to be sampled from.
#'
#' 3. Sample from the conditional distributions for the infection and removal rate
#'    in turn
#'
#' 4. Sample an Infection Time to change. Propose the change and Accept/Reject using
#'    Metropolis-Hastings Acceptance probability
#'

# ==== Documentation ====

# === Epidemic Parameters ===

# N := Number of Individuals in the closed Population of the Epidemic
# I := infection times
# R := removal times
# rinf.dist := Chosen Distribution of the infectious period of individuals
# alpha := Shape Parameter for Infectious Period Distribution (If Gamma or Weibull)

# === Bayesian Construction Parameters ===

# gamma.prior := Prior distribution for gamma parameter in Bayesian Construction
# beta.prior := See 'gamma.prior'

# theta.gamma := hyperparameters for gamma prior
# theta.beta := '                   ' beta prior




# === MCMC Parameters ===

# gamma0 := intial value of gamma parameter
# beta0 := intial value of beta parameter
# no.proposals := How many infection times are to be proposed per iteration
# no.its := Number of iterations of Sampling to be carried out
# burn.in := The number of iterations to be thrown away as a period of convergence to the
#            target distribution
# thinning.factor := The factor by which samples are thinned to reduce lag in the simulated Markov Chain
# lag.max := A parameter to be passed to the acf() function. Controls the maximum lag of acf to be calculated/plotted

# ==== Metropolis Within Gibbs MCMC Algoritm for a Centered Parameterisation ====


inf_time_proposal = function(t_inf, t_rem, no_proposals, par_rem){
  which_infected = which(t_inf < Inf)
  t_inf_prop = t_inf
  proposed_infected  = sample(which_infected, no_proposals)
  t_inf_prop[proposed_infected] = t_rem[proposed_infected] - rgamma(no_proposals, rate = par_rem[1],
                                                                    shape = par_rem[2])
  return(t_inf_prop)
}

Centered_MCMC = function(N, a, t_rem, gamma0, theta_gamma, beta0, theta_beta, kernel, no_proposals, no_its, burn_in,
                         PLOT = TRUE, thinning_factor = 1, lag_max = NA){
  Start <- as.numeric(Sys.time())
  # == Initialise ==

  # Number of people who were removed
  n_R <- sum(t_rem < Inf)
  n_I <- n_R - a # For finished Epidemic n_I = n_R

  # Which individuals got infected?
  which_infected <- which(t_rem < Inf)

  # Rescale
  t_rem <- t_rem - min(t_rem)

  # Initialise Beta and Gamma Values
  beta <- beta0
  gamma <- gamma0

  inf_period <- rep(0, N)
  t_inf <- rep(Inf, N)

  # Draw infectious periods
  inf_period[which_infected] <- rexp(n_R, rate = gamma)

  # Calculate infection times
  t_inf <- t_rem - inf_period

  # == Checking validity of drawn infection times ==

  waifw = sapply(t_inf[which_infected], function(t) cbind(t_inf, t_rem)[which_infected, 1] < t & t < cbind(t_inf, t_rem)[which_infected, 2])

  while(sum(colSums(waifw) > 0) != n_R - 1){

    # Widen infectious periods
    inf_period <- inf_period*1.1
    #inf_period[which_infected] <- rexp(n_R, rate = gamma)
    # Calculate new infection times
    t_inf <- t_rem - inf_period

    waifw = sapply(t_inf[which_infected], function(t) cbind(t_inf, t_rem)[which_infected, 1] < t & t < cbind(t_inf, t_rem)[which_infected, 2])
    print(sum(colSums(waifw) > 0))
  }

  # Calculate components of posterior parameters for beta and the likelihood
  # associated with these infection times and removal times

  # Infectious Pressure Integral
  IP_integral <- integral_part_inf(events = cbind(t_inf, t_rem), t_inf_j = t_inf, with.beta = FALSE)

  # Removal Pressure Integral
  RP_integral <- sum(t_rem[which_infected] - t_inf[which_infected])

  # Empty matrix to store samples
  draws <- matrix(NA, nrow = no_its + 1, ncol = N + 2)

  # First entry is the intial draws
  draws[1,] <- c(beta, gamma, t_inf)


  # Acceptance Counter
  accept <- 0

  # == The Sampler ==

  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = no_its)
  for(i in 1:no_its){
    pb$tick()
    # == Drawing Beta and Gamma (Gibbs Step) ==

    # Update beta and gamma by drawing from their conditional posteriors

    # Exponential Prior ---> Gamma Posterior
    gamma <-  rgamma(1, shape = n_R + theta_gamma[1], rate = theta_gamma[2] + RP_integral)
    beta  <-  rgamma(1, shape = n_I + theta_beta[1], rate = theta_beta[2] + IP_integral)

    loglikelihood_inf = log_likelihood_inf(beta, t_inf, t_rem, kernel)

    # ==== Updating Infection time(s) ====

    #t_inf_prop = inf_time_proposal(t_inf, t_rem, no_proposals, c(gamma, alpha))

    # Choose an infected individual at random
    proposed_infected <- sample(which_infected, no_proposals)

    # Propose a new infectious period
    inf_period_prop <- rexp(no_proposals, rate = gamma)

    # Use this to calculate the new infection time
    t_inf_prop <- t_inf
    t_inf_prop[proposed_infected] <- t_rem[proposed_infected] - inf_period_prop

    # == Calculate the likelihood components with proposed infection time ==
    loglikelihood_inf_prop = log_likelihood_inf(beta, t_inf_prop, t_rem, kernel)

    log_u <- log(runif(1))
    log_a <- (loglikelihood_inf_prop) -
                (loglikelihood_inf)
    # + sum(dgamma(inf_period[proposed_infected], shape = alpha, rate = gamma, log = TRUE))
    # + sum(dgamma(inf_period_prop, shape = alpha, rate = gamma, log = TRUE))

    if(log_u < log_a){
      t_inf <- t_inf_prop
      IP_integral <- integral_part_inf(events = cbind(t_inf, t_rem), t_inf_j = t_inf, with.beta = FALSE)
      RP_integral <- sum(t_rem[which_infected] - t_inf[which_infected])
      accept <- accept + 1
    }
    draws[i+1,] <- c(beta, gamma, t_inf)
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
  accept_rate <- accept/no_its

  # = Plots =
  par(mfrow = c(2,2))

  if(PLOT){
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

  }

  # Calculating Summary Statistics for samples
  beta_summary = c(mean(draws[,1]), var(draws[,1]))
  gamma_summary = c(mean(draws[,2]), var(draws[,2]))
  R0_summary = c(mean(R0_samples), var(R0_samples))

  printed_output(rinf_dist = "Exp", no_proposals, no_its, ESS, time_taken, ESS_sec, accept_rate)
  return(list(draws = draws, accept_rate = accept_rate, ESS = ESS, ESS_sec = ESS_sec, beta_summary = beta_summary,
              gamma_summary = gamma_summary, R0_summary = R0_summary, time_taken = time_taken))
}





