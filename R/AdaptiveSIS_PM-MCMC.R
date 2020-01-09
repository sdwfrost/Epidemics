
#' Adaptive Pseudo-Marginal MCMC for SIS Epidemic Panel Data
#'
#' Adapts proposal parameters for with a view of optimal of a target using Pseudo-Marginal MCMC scheme.

#' @family Panel Data MCMC
#' @param obsTransData Interpanel transition data.
#' @param I_0 Initial number of infectives in the population.
#' @param obsTimes Times at which epidemic cohort were followed up.
#' @param N Population size.
#' @param beta0 Starting value for infectious process parameter.
#' @param gamma0 Starting value for removal/recovery process parameter.
#' @param lambda0 Starting value for RWM proposal parameter which is to be adapted.
#' @param V0 Starting state for RWM proposal Covariance matrix which is to be adapted.
#' @param delta Probability that a proposal is made based on the starting proposal conditions.
#' @param noSims Number of epidemic simulations to be used per iteration for estimation of likelihood.
#' @param noIts Number of MCMC iterations.
#' @param lagMax Plotting parameter for acf() function.
#' @param thinningFactor Controls the factor by which MCMC samples are thinned, to reduce dependency.
#' @param parallel Are epidemic simulations run in parallel?
#' @param noCores If epidemic simulations are run in parallel, this is the number of cores utilised.
#'
#' @return Proposal parameters which can be used for more optimal exploration of target distribution (plus MCMC summary).


adaptiveSIS_PseudoMarginalMCMC = function(obsTransData, I_0, obsTimes, N, beta0, gamma0, lambda0 = 2.38/sqrt(2), V0 = diag(c(1/N, 1)), noSims,
                                      noIts, burnIn, lagMax = NA, thinningFactor = 1, parallel = FALSE, noCores,
                                      delta = 0.05){

  Start = as.numeric(Sys.time())

  lambda = lambda0
  noSampled = sum(obsTransData[[1]])
  thetaCurr = c(beta0, gamma0)
  logPCurr = -Inf

  while(logPCurr == -Inf){
    panelDataSim = replicate(noSims, homogeneousPanelDataSIS_Gillespie(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr[1],
                                                                       thetaCurr[2], obsTimes)$panelData, simplify = FALSE)
    transDataSim = sapply(X = panelDataSim, function(X) transitionData(X, states = 1:2), simplify = FALSE)
    logPCurr = mean(sapply(X = transDataSim, function(X) dHyperGeom(obsTransData, X, noSampled, log = T)))
  }

  #' Create Storage Matrix
  draws = matrix(NA, nrow = noIts + 1, ncol = length(thetaCurr) + 1)
  draws[1,] = c(thetaCurr, logPCurr)

  #' Proposal Acceptance Counter
  accept = 0

  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = noIts)

  for(i in 1:noIts){
    pb$tick()
    #' Propose new beta and gamma using Multiplicative RW propsal
    u1 = runif(1, 0, 1)
    if(u1 > delta & accept > 10){
      Vi = var(draws[,1:2], na.rm = T)
      thetaProp = abs(thetaCurr + lambda*mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = Vi))
    } else{
      thetaProp = abs(thetaCurr + lambda0*mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = V0))
    }

    #print(sim$noDraws)
    #print(thetaCurr)
    #print(accept)

    panelDataSim = replicate(noSims, homogeneousPanelDataSIS_Gillespie(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaProp[1],
                                                                       thetaProp[2], obsTimes)$panelData, simplify = FALSE)
    transDataSim = sapply(X = panelDataSim, function(X) transitionData(X, states = 1:2), simplify = FALSE)
    logPProp = mean(sapply(X = transDataSim, function(X) dHyperGeom(obsTransData, X, noSampled, log = T)))

    logU = log(runif(1))

    logA = (logPProp + sum(dexp(thetaProp, rate = 0.001, log = TRUE))) -
      (logPCurr + sum(dexp(thetaCurr, rate = 0.001, log = T)))

    if(logU < logA){
      logPCurr = logPProp
      thetaCurr = thetaProp
      accept = accept + 1
      if(u1 > delta & accept > 10){
        lambda = lambda + 0.93*(lambda/sqrt(i))
      }
    } else{
      if(u1 > delta & accept > 10){
        lambda = lambda - 0.07*(lambda/sqrt((i)))
      }
    }
    #' Store State
    draws[i + 1, ] = c(thetaCurr, logPCurr)
    # print(c(thetaCurr))
    # print(lambda)
  }

  End <- as.numeric(Sys.time())

  timeTaken <- End - Start

  # Thin the samples
  draws <- draws[seq(from = burnIn + 1, to = (noIts + 1) - burnIn, by = thinningFactor),]

  # Calculate Effective Sample Sizes (and Per Second) and Acceptance Rates
  ESS <- coda::effectiveSize(draws[, 1:2])
  ESS.sec <- ESS/timeTaken
  acceptRate <-  accept/noIts

  # = Plots =
  par(mfrow = c(2,2))

  # Plot Beta Samples and Sample Auto-Corrolation Function
  if(is.na(lagMax)){
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
  betaSummary = c(mean(draws[,1]), sd(draws[,1]))
  gammaSummary = c(mean(draws[,2]), sd(draws[,2]))
  if(accept < 10){
    Vi = NULL
  }

  printed_output(rinf_dist = "Exp", no_proposals = NA, noIts, ESS, timeTaken, ESS.sec, acceptRate)
  return(list(draws = draws, acceptRate = acceptRate, lambda = lambda, V = Vi, ESS = ESS, ESS.sec = ESS.sec, betaSummary = betaSummary,
              gammaSummary = gammaSummary, timeTaken = timeTaken))
}
