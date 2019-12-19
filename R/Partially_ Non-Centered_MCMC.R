#'
#' PhD
#' Partially Non-Centered Sampling Algorithm for Epidemics
#'
#'

# ==== Preamble ====

# ==== Partially Non-Centered Sampling Algorithm ====

PNC.Algorithm <- function(N, a, t_rem, gamma0, beta0, theta.gamma, theta.beta,
                          rinf.dist, no_proposals, lambda, centering_prob, no_its, burn_in, thinning_factor = 1,
                          lag_max = NA){

  # Starting Time
  Start <- as.numeric(Sys.time())

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
  inf_period[which_infected] <- rexp(n_I, rate = gamma)

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
  }

  # Use these infectious periods and the known removal times to
  # to calculate infection times.
  t_inf <- t_rem - inf_period

  # == Checking validity of drawn infection times ==

  # Each infection time must lie within the infectious period of another
  # individual (excluding the initial infected individual)
  valid.infected <- 0
  while(valid.infected < n_I - 1){
    valid.infected <- 0
    # Valid infection times check
    for(i in infected){
      if(sum(I[-i] < I[i] & I[i] < R[-i]) > 0){
        valid.infected <- valid.infected + 1
      }
    }
    #print(valid.infected)
    # If the infection times do not create a valid epidemic, we can
    # widen the length of the infectious periods so that there is
    # a higher chance of crossover in infectious periods
    if(valid.infected < n_I - 1){

      # Widen infectious periods
      Q <- Q*1.1

      # Calculate new infection times
      I <- R - Q
    }
  }

  # Calculate Us for Non-Centered Parameterisation
  U <- gamma*Q

  # Initial Choice off Centered and Non-Centered Infection Times

  # We center each Infection Time with probability 'cenetering.prob'
  Z <- rep(NA, N)
  Z[infected] <- rbinom(1, 1, prob = centering.prob)

  # Set of Centered Individuals
  centered = which(Z == 1)

  # Calculate U from the intial Infection Times
  U <- gamma*Q

  # Calculate components of posterior parameters for beta and the likelihood
  # associated with these infection times and removal times

  # Infectious Pressure Integral
  IP.int <- IP.integral(I, R)

  # Infectious Process Product
  LIP <- Likelihood.Infection.Product(I, R, log = TRUE)

  # Removal Pressure Integral
  RP.int.centered <- RP.integral(I[centered], R[centered])
  #RP.int.centered <- RP.integral(I, R)

  # Empty matrix to store samples
  draws <- matrix(NA, nrow = no.its + 1, ncol = 2*N + 2)


  # First entry is the intial draws
  draws[1,] <- c(beta, gamma, I, Z)

  # Gamma Proposal Counter
  no.gamma.props <- 0

  # Acceptance Counters
  accept.inf <- 0
  accept.gamma <- 0

  # Number of Accepted Proposals which were centered
  accept.cenetered <- 0

  # ==== Sampling ====

  print("Sampling Progress")
  pb <- progress_bar$new(total = no.its)

  for(i in 1:no.its){

    pb$tick()

    # == Beta Sample ==

    # Exponential Prior --> Gamma Posterior
    beta  <-  rgamma(1, shape = (n_I - 1) + 1, rate = theta.beta + IP.int)

    # == Gamma Sample ==

    # Some Infection Times Centered, some non-centered.
    no.gamma.props <- no.gamma.props + 1

    # gamma.prop = gamma + lambda*rnorm(1)
    log.gamma.prop <- log(gamma) + lambda*rnorm(1, mean = 0, sd = 1)
    gamma.prop <- exp(log.gamma.prop)

    # = Likelihood Components with proposed gamma =
    I.prop <- I
    I.prop[Z == 0] <- R - (1/gamma.prop)*U[Z == 0]

    LIP.prop <- Likelihood.Infection.Product(I.prop, R, log = TRUE)
    IP.int.prop <- IP.integral(I.prop, R)

    # = Accept/Reject Step =

    # log.alpha = ( LIP.prop - beta*IP.int.prop - gamma.prop*theta.gamma) -
    #  (LIP - beta*IP.int - gamma*theta.gamma)

    log.alpha <- (LIP.prop - beta*IP.int.prop - gamma.prop*(theta.gamma + RP.int.centered) +
                    log.gamma.prop) -
      (LIP - beta*IP.int - gamma*(theta.gamma + RP.int.centered) + log(gamma))

    log.u <- log(runif(1))

    if(log.u < log.alpha){
      gamma <- gamma.prop
      I <- I.prop
      U <- (R - I)*gamma.prop
      LIP <- LIP.prop
      IP.int <- IP.int.prop
      accept.gamma <- accept.gamma + 1
      #print(accept.gamma)
    }


    # == Infection Time Sample ==

    # First Choose which Infection time to update
    proposed.infected <- sample(infected, no.proposals)


    # Infection Time proposals are equivalent in both the centered
    # and non-centered parameterisations of the likelihood.
    # We will work in terms of U here, but make sure to update I
    # as Us change
    # == Non-Centered Proposal ==

    U.prop <- U
    U.prop[proposed.infected] <- rexp(no.proposals, rate = 1)
    I.prop <- R - (1/gamma)*U.prop

    # = Likelihood Components with Proposed Us =

    IP.int.prop <- IP.integral(I.prop, R)
    LIP.prop <- Likelihood.Infection.Product(I.prop, R, log = TRUE)

    # = Accept/Reject =

    log.alpha <- (LIP.prop - beta*IP.int.prop) -
                 (LIP - beta*IP.int)

    log.u <- log(runif(1))

    if(log.u < log.alpha){
      U <- U.prop
      I <- I.prop
      LIP <- LIP.prop
      IP.int <- IP.int.prop
      accept.inf <- accept.inf + 1
    }

    # ==== Resample Zs ====

    Z[infected] <- rbinom(n_I, 1, prob = centering.prob)
    centered = which(Z == 1)
    RP.int.centered = RP.integral(I[centered], R[centered])
    #RP.int.centered = RP.integral(I, R)
    # Store Samples
    draws[i+1,] <- c(beta, gamma, I, Z)
  }

  End <- as.numeric(Sys.time())
  time.taken <- End - Start

  # Calculate R_0 sample values
  R_0.samples = (N*draws[,1])/draws[,2]

  # Calculate Effective Sample Sizes (and Per Second) and Acceptance Rates
  ESS <- c(effectiveSize(draws[, 1:2]), effectiveSize(R_0.samples))
  ESS.sec <- ESS/time.taken
  accept.rate <- c(accept.gamma/no.gamma.props, accept.inf/no.its)

  # = Plots =
  par(mfrow = c(2,2))

  # Plot Beta Samples and Sample Auto-Corrolation Function
  if(is.na(lag.max)){
    # Beta
    plot(draws[, 1], type = 'l')
    acf(draws[, 1])

    # Gamma
    plot(draws[, 2], type = 'l')
    acf(draws[, 2])
  } else{
    # Beta
    plot(draws[, 1], type = 'l')
    acf(draws[, 1], lag.max)

    # Gamma
    plot(draws[, 2], type = 'l')
    acf(draws[, 2], lag.max)
  }

  # Calculating Summary Statistics for samples
  beta.summary = c(mean(draws[,1]), sd(draws[,1]))
  gamma.summary = c(mean(draws[,2]), sd(draws[,2]))
  R_0.summary = c(mean(R_0.samples), sd(R_0.samples))

  # = Summary Information of MCMC Performance =
  printed.output(rinf.dist, no.proposals, no.its, ESS, time.taken, ESS.sec, accept.rate)

  # = Output =
  return(list(draws = draws, accept.rate = accept.rate, ESS = ESS, ESS.sec = ESS.sec, beta.summary = beta.summary,
              gamma.summary = gamma.summary, R_0.summary = R_0.summary, time.taken = time.taken))
}


