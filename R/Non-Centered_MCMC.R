#'
#'
#' MCMC Sampler for a non-centered Bayesian Parameterisation of 
#' Epidemics
#'
#'

# ==== Preamble ====

# ==== Documentation ====

# === Epidemic Parameters ===

# N := Number of Individuals in the closed Population of the Epidemic
# I := infection times
# R := removal times
# rinf.dist := Chosen Distribution of the infectious period of individuals

# === Bayesian Construction Parameters ===

# gamma.prior := Prior distribution for gamma parameter in Bayesian Construction
# beta.prior := See 'gamma.prior'

# theta.gamma := hyperparameters for gamma prior
# theta.beta := '                   ' beta prior


# === MCMC Parameters ===

# gamma0 := intial value of gamma parameter
# beta0 := intial value of beta parameter
# no.proposals := How many infection times are to be proposed per iteration
# lambda := Tuning parameter of RWM Proposal Jumps for gamma parameter
# no.its := Number of iterations of Sampling to be carried out
# burn.in := The number of iterations to be thrown away as a period of convergence to the
#            target distribution
# thinning.factor := The factor by which samples are thinned to reduce lag in the simulated Markov Chain
# lag.max := A parameter to be passed to the acf() function. Controls the maximum lag of acf to be calculated/plotted

# ==== Non-Centered Algorithm ====

NonCentered_MCMC <- function(N, t_rem, gamma0, alpha = 1, theta_gamma, beta0, theta_beta, kernel, no_its,
                         burn_in = 0, no_proposals, lambda, thinning_factor = 1, lag_max = NA){
  
  # ==== Initialise ====
  
  Start <- as.numeric(Sys.time())
  Start <- as.numeric(Sys.time())
  # == Initialise ==
  
  # Number of people who were removed 
  n_R <- sum(t_rem < Inf) 
  n_I <- n_R # For finished Epidemic n_I = n_R
  
  # Which individuals got infected?
  which_infected <- which(t_rem < Inf)
  
  # Rescale 
  t_rem <- t_rem - min(t_rem)
  
  # Initialise Beta and Gamma Values
  beta <- beta0
  B = kernel(beta0)
  gamma <- gamma0
  
  inf_period <- rep(0, N)
  t_inf <- rep(Inf, N)
  
  # Draw infectious periods
  inf_period[which_infected] <- rgamma(n_I, rate = gamma, shape = alpha)
  
  # Calculate infection times
  t_inf <- t_rem - inf_period
  
  # == Checking validity of drawn infection times ==
  
  waifw = sapply(t_inf[which_infected], function(t) cbind(t_inf, t_rem)[which_infected, 1] < t & t < cbind(t_inf, t_rem)[which_infected, 2])
  
  while(sum(colSums(waifw) > 0) != n_I - 1){
    
    # Widen infectious periods
    inf_period <- inf_period*1.1
    
    # Calculate new infection times
    t_inf <- t_rem - inf_period
    
    waifw = sapply(t_inf[which_infected], function(t) cbind(t_inf, t_rem)[which_infected, 1] < t & t < cbind(t_inf, t_rem)[which_infected, 2])
    print(sum(colSums(waifw) > 0))
  }
  
  
  
  # Calculate U from the intial Infection Times 
  U <- gamma*inf_period
  
  # Calculate components of posterior parameters for beta and the likelihood
  # associated with these infection times and removal times
  
  # Infectious Pressure Integral
  IP_integral <- integral_part_inf(events = cbind(t_inf, t_rem), t_inf_j = t_inf, with.beta = FALSE)
  
  # Empty matrix to store samples
  draws <- matrix(NA, nrow = no_its + 1, ncol = N + 2)
  
  
  # First entry is the intial draws
  draws[1,] <- c(beta, gamma, t_rem)
  
  # Counter for acceptance of proposals
  accept_gamma <- 0 
  accept_U <- 0
  
  # ==== Sampling ====
  print("Sampling Progress")
  pb <- progress_bar$new(total = no_its)
  for(i in 1:no_its){
    pb$tick()
    
    # == Beta ==
    beta  <-  rgamma(1, shape = (n_I - 1) + 1, rate = theta_beta + IP_integral)
    
    # Log-Likelihood
    loglikelihood = epidemic_loglikelihood(c(gamma, alpha), beta, t_inf, t_rem, kernel)
    
    
    # == Gamma ==
    
    log_gamma_prop <- log(gamma) + lambda*rnorm(1, mean = 0, sd = 1)
    gamma_prop <- exp(log_gamma_prop)
    
    # = Likelihood Components with proposed gamma = 
    t_inf_prop <- t_rem - (1/gamma_prop)*U 
    
    loglikelihood_prop = epidemic_loglikelihood(c(gamma_prop, alpha), beta, t_inf_prop, t_rem, kernel)
    
    # = Accept/Reject Step = 

    
    log_a <- (loglikelihood_prop + dexp(gamma_prop, rate = theta_gamma, log = TRUE) + log_gamma_prop) -
             (loglikelihood + dexp(gamma, rate = theta_gamma, log = TRUE) + log(gamma))
    
    log_u <- log(runif(1))

    if(log_u < log_a){
      
      # Update Parameters/Likelihood Components
      gamma <- gamma_prop
      t_inf <- t_inf_prop
      loglikelihood = loglikelihood_prop
      IP_integral <- integral_part_inf(events = cbind(t_inf, t_rem), t_inf_j = t_inf, with.beta = FALSE)
      
      # Sync t_inf and U
      U[which_infected] <- gamma*(t_rem[which_infected] - t_inf[which_infected])
      
      accept_gamma <- accept_gamma + 1
    }
    
    # == Augmented Data U ==
    
    # Randomly select infectives to propose a new U for
    proposed_infected <- sample(which_infected, no_proposals)
    
    # Draw the new Us using Exp(1)
    
    U_prop <- U
    U_prop[proposed_infected] <- rgamma(no_proposals, rate = 1, shape = alpha)
    
    # Sync t_inf_prop and U_prop
    t_inf_prop <- t_rem - (1/gamma)*U_prop
    
    # Proposed Log-likelihood
    loglikelihood_prop = epidemic_loglikelihood(c(gamma, alpha), beta, t_inf_prop, t_rem, kernel)

    # = Accept/Reject =
    
    log_a <- (loglikelihood_prop + sum(dgamma(U[proposed_infected], rate = gamma, shape = alpha, log = TRUE))) -
      (loglikelihood + sum(dgamma(U_prop[proposed_infected], rate = gamma, shape = alpha, log = TRUE)))
    
    log_u <- log(runif(1))
    
    if(log_u < log_a){
      U <- U_prop
      t_inf <- t_inf_prop
      IP_integral <- integral_part_inf(events = cbind(t_inf, t_rem), t_inf_j = t_inf, with.beta = FALSE)
      accept_U <- accept_U + 1
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
  ESS <- c(effectiveSize(draws[, 1:2]), effectiveSize(R0_samples))  
  ESS_sec <- ESS/time_taken
  accept_rate <- c(accept_gamma, accept_U)/no_its
  
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
  
  # Calculating Summary Statistics for samples
  beta_summary = c(mean(draws[,1]), sd(draws[,1]))
  gamma_summary = c(mean(draws[,2]), sd(draws[,2]))
  R0_summary = c(mean(R0_samples), sd(R0_samples))
  
  # = Summary Information of MCMC Performance =
  printed.output("NA",no_proposals, no_its, ESS, time_taken, ESS_sec, accept_rate)
  
  return(list(draws = draws, accept_rate = accept_rate, ESS = ESS, ESS_sec = ESS_sec,
              beta_summary = beta_summary, gamma_summary = gamma_summary, R0_summary = R0_summary,
              time_taken = time_taken))
}

