#'
#'
#' Ancillarity–Sufficiency Interweaving Strategy (ASIS) MCMC
#'
#'
#'
#'

# The ASIS MCMC algorithm makes back-to-back Centred and Non-Centred proposals
# in a bid to switch between two different types of MCMC updates

# ==== Preamble ====


# ==== ASIS Sampler ====

ASIS.sampler <- function(N, R, gamma0, beta0, gamma.prior, theta.gamma, beta.prior, theta.beta, rinf.dist, alpha, no.its,
                       burn.in = 0, no.prop.CP, no.prop.NCP, lambda, thinning.factor = 1, lag.max = NA){

  # ==== Initialise ====

  Start <- as.numeric(Sys.time())

  if(((gamma.prior == "gamma" & length(theta.gamma) != 2) |
      (gamma.prior == "exp" & length(theta.gamma) != 1))|
     ((beta.prior == "gamma" & length(theta.beta) != 2) |
      (beta.prior == "exp" & length(theta.beta) != 1))){
    stop("Number of prior hyperparameters specified does not match
         the number required for the prior assigned to beta/gamma")
  }


  #beta = rprior(1, dist = beta.prior, theta = theta.beta)
  #gamma = rprior(1, dist = gamma.prior, theta = theta.gamma)

  beta <- beta0
  gamma <- gamma0

  # Number of people who were removed
  n_R <- sum(R < 10000)

  # For finished epidemics, the number of people who were infected
  # is equal to the number who were removed.
  #n_I = sum(I < 10000)
  n_I <- n_R

  # Which individuals got infected?
  infected <- which(R < 10000)

  # Rescale so that first removal time is at t = 0
  # This time is when the epidemic would
  # begin to be observed, in reality
  R[infected] <- R[infected] - min(R[infected])

  Q <- rep(0, N)
  I <- rep(10000, N)
  # Draw infectious periods
  Q[infected] <- rinf.period(n = n_I, dist = rinf.dist,
                            theta = switch(rinf.dist,
                                           exp = gamma, gamma = c(gamma, alpha),
                                           weibull = c(gamma, alpha)))

  # Use these infectious periods and the known removal times to
  # to calculate infection times.
  I <- R - Q

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

  # Calculate U from the intial Infection Times
  U <- gamma*Q

  # Calculate components of posterior parameters for beta and the likelihood
  # associated with these infection times and removal times

  # Infectious Pressure Integral
  IP.int <- IP.integral(I, R)

  # Infectious Process Product
  LIP <- Likelihood.Infection.Product(I, R, log = TRUE)

  # Removal Pressure Integral
  RP.int <- RP.integral(I, R)

  # Empty matrix to store samples
  draws <- matrix(NA, nrow = no.its + 1, ncol = N + 2)


  # First entry is the intial draws
  draws[1,] <- c(beta, gamma, I)

  # Counter for acceptance of proposals
  accept.gamma <- 0
  accept.infection <- 0

  # ==== Sampling ====
  print("Sampling Progress")
  pb <- progress_bar$new(total = no.its/2)
  for(i in 1:(no.its/2)){
    pb$tick()
    # === Centred Sampling Step ===

    # == Drawing Beta and Gamma (Gibbs Step) ==

    # Update beta and gamma by drawing from their conditional posteriors

    # Gamma Prior ---> Gamma Posterior
    #gamma <-  rgamma(1, shape = theta.gamma[1] + n_R, rate = theta.gamma[2] + sum(R - I, na.rm = TRUE) )
    #beta  <-  rgamma(1, shape = theta.beta[1] + (n_I - 1) + 1, rate = theta.beta[2] + inf.pressure )

    # Exponential Prior ---> Gamma Posterior
    gamma <-  rgamma(1, shape = n_R + 1, rate = theta.gamma + RP.int)
    beta  <-  rgamma(1, shape = (n_I - 1) + 1, rate = theta.beta + IP.int)

    # ==== Updating Infection time(s) ====

    # Choose an infected individual at random
    proposed.infected <- sample(infected, no.prop.CP)

    # Propose a new infectious period
    Q.proposal <- rinf.period(n = no.prop.CP, dist = rinf.dist,
                             theta = switch(rinf.dist,
                                            exp = gamma, gamma = c(gamma, alpha),
                                            weibull = c(gamma, alpha)))

    # Use this to calculate the new infection time
    I.prop <- I
    I.prop[proposed.infected] <- R[proposed.infected] - Q.proposal

    # == Calculate the likelihood components with proposed infection time ==

    IP.int.prop <- IP.integral(I.prop, R)
    LIP.prop <- Likelihood.Infection.Product(I.prop, R)
    RP.int.prop <- RP.integral(I.prop, R)

    log.u <- log(runif(1))
    log.a <- (LIP.prop - beta*IP.int.prop) - (LIP - beta*IP.int)

    if(log.u < log.a){
      I <- I.prop
      LIP <- LIP.prop
      IP.int <- IP.int.prop
      RP.int <- RP.int.prop
      U <- gamma*(R - I.prop)
      accept.infection <- accept.infection + 1
    }
    #print(accept.infection)
    draws[(2*i + 1) - 1,] <- c(beta, gamma, I)


    # === Non-Centred Sampling Step ===

    # == Beta ==

    # Conditional Posterior will be of a simple form
    # Therefore we can just Gibbs Sample Beta

    # Exponential Prior --> Gamma Posterior
    beta  <-  rgamma(1, shape = (n_I - 1) + 1, rate = theta.beta + IP.int)

    # == Gamma ==

    # Conditional Posterior of Gamma will be more complicated
    # and Gibbs is not possible

    # Therefore, an alternative MH proposal must be used

    # RWM (Optimal Acceptance Prob.y will be around 23.4%)

    #gamma.prop <- gamma + lambda*rnorm(1)
    log.gamma.prop <- log(gamma) + lambda*rnorm(1, mean = 0, sd = 1)
    gamma.prop <- exp(log.gamma.prop)

    # = Likelihood Components with proposed gamma =
    I.prop <- R - (1/gamma.prop)*U
    LIP.prop <- Likelihood.Infection.Product(I.prop, R, log = TRUE)
    IP.int.prop <- IP.integral(I.prop, R)

    # = Accept/Reject Step =

    #log.alpha <- ( LIP.prop - beta*IP.int.prop - gamma.prop*theta.gamma) -
    #  (LIP - beta*IP.int - gamma*theta.gamma)

    log.alpha <- (LIP.prop - beta*IP.int.prop - gamma.prop*theta.gamma +
                   log.gamma.prop) -
      (LIP - beta*IP.int - gamma*theta.gamma + log(gamma))

    log.u <- log(runif(1))

    if(log.u < log.alpha){
      gamma <- gamma.prop
      I <- I.prop
      LIP <- LIP.prop
      IP.int <- IP.int.prop
      gamma.accepted <- TRUE
      accept.gamma <- accept.gamma + 1
    } else{
      gamma.accepted <- FALSE
    }

    # == Augmented Data U ==

    # Randomly select infectives to propose a new U for
    proposed.infected <- sample(infected, no.prop.NCP)

    # Draw the new Us using Exp(1)

    U.prop <- U
    U.prop[proposed.infected] <- rexp(no.prop.NCP, rate = 1)
    I.prop <- R - (1/gamma)*U.prop

    # = Likelihood Components with Proposed Us =

    IP.int.prop <- IP.integral(I.prop, R)
    LIP.prop <- Likelihood.Infection.Product(I.prop, R, log = TRUE)

    # = Accept/Reject =

    log.alpha <- (LIP.prop - beta*IP.int.prop) -
      (LIP - beta*IP.int)

    log.u <- log(runif(1))

    if(log.u < log.alpha){
      infection.accepted <- TRUE
      U <- U.prop
      I <- I.prop
      LIP <- LIP.prop
      IP.int <- IP.int.prop
      RP.int <- RP.integral(I, R)
      accept.infection <- accept.infection + 1
    } else{
      infection.accepted <- FALSE
    }

    # RP.integral does not need to be computed for the U proposal step.
    # However, infection times will change if this proposal is accepted
    # If both the Gamma and U propsals are accepted with RP.int being
    # calculated after each acceptance, then RP.int will
    # be recalculated uneccesarily at the Gamma acceptance. Therefore we
    # wait to see whether the U proposal is accepted to recalculate.
    if(gamma.accepted == TRUE & infection.accepted == FALSE){
      RP.int <- RP.integral(I, R)
    }
    draws[2*i + 1,] <- c(beta, gamma, I)
    #print(accept.infection)
    #print(accept.gamma)
  }

  End <- as.numeric(Sys.time())
  time.taken <- End - Start

  # == Diagnostics ==

  # Remove Burn In
  draws <- draws[-(1:burn.in),]

  # Thin the samples
  draws <- draws[seq(from = burn.in + 1, to = (no.its + 1) - burn.in, by = thinning.factor),]


  # Calculate R_0 sample values
  R_0.samples = (N*draws[,1])/draws[,2]

  # Calculate Effective Sample Sizes (and Per Second) and Acceptance Rates
  ESS <- c(effectiveSize(draws[, 1:2]), effectiveSize(R_0.samples))
  ESS.sec <- ESS/time.taken
  accept.rate <- c(2*accept.gamma, accept.infection)/no.its

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
  printed.output(rinf.dist, c(no.prop.CP, no.prop.NCP), no.its, ESS, time.taken, ESS.sec, accept.rate)

  # = Output =
  return(list(draws = draws, accept.rate = accept.rate, ESS = ESS, ESS.sec = ESS.sec,
              gamma.summary = gamma.summary, beta.summary = beta.summary, R_0.summary = R_0.summary,
              time.taken = time.taken))
}



