#'
#' (Multiplicative) Random-Walk Metropolis for Heterogeneous SIR Epidemic
#'
#'

#' Full information on epidemic
#'
#' Heterogeneous Kernel set up
#'
#' Infectious model is parameterised by a predefined kernel.
#'
#' Removal process is independent of infectious process a posteriori
#'
#' Hence no need to sample removal parameter,
#' given we have chosen a suitable prior (i.e Exponential/Gamma)
#'
#' Paired Infection-Removal times.
#'


RWMSampler =  function(beta0, betaPrior = c(1, 0.001), gammaPrior = c(1,0.001), N,
                       eventTimes, kernel, lambda, V = 1, noIts, burnIn, thinFactor = 1, lagMax = NULL){
  Start <- as.numeric(Sys.time())
  no_beta = length(beta0)

  #' Initialise
  betaCurr = beta0
  #betaCurr = rgamma(no_beta, betaPrior[,1], betaPrior[,2])

  #' Calculate Infectious Process part of the likelihood (Removal Process will cancel out).
  logLikelihoodCurr = log_likelihood_inf(betaCurr, eventTimes[,1], eventTimes[,2], kernel)

  if(no_beta > 1){
    betaCurrPrior = sum(mapply(dgamma, betaCurr, shape = betaPrior[,1], rate = betaPrior[,2],
                               MoreArgs = list(log = TRUE)))
  } else{
    betaCurrPrior = dgamma(betaCurr, betaPrior[1], betaPrior[2], log = TRUE)
  }

  i_removed = eventTimes[,1] < Inf
  i_infected = eventTimes[,2] < Inf
  #' Gamma Posterior Parameters
  gammaPosterior = c(gammaPrior[1] + sum(eventTimes[i_removed,2] - eventTimes[i_infected,1]),
                     gammaPrior[2] + sum(eventTimes[,2] < Inf))

  gammaCurr = rgamma(1, gammaPosterior[1], gammaPosterior[2])
  # betaPosterior = c(betaPrior[1] + integral_part_inf(eventTimes, eventTimes[,1], with.beta = FALSE),
  #                   betaPrior[2] + prod_part_inf(eventTimes[,1], eventTimes, B = matrix(1, nrow = N, ncol = N)))


  #' Data
  results = matrix(NA, nrow = noIts + 1, ncol = no_beta + 1)
  results[1, ] = c(betaCurr, gammaCurr)
  acceptCounter = 0

  for(i in 1:noIts){

    logBetaCurr = log(betaCurr)

    logBetaProposal = logBetaCurr + mvtnorm::rmvnorm(1, mean = rep(0, no_beta),
                                                    sigma = lambda*V)

    betaProposal = as.numeric(exp(logBetaProposal))

    #' Accept/Reject
    log_u = log(runif(1))

    logLikelihoodProp = log_likelihood_inf(betaProposal, eventTimes[,1], eventTimes[,2], kernel)

    #betaProposalPrior = dgamma(betaProposal, betaPrior[,1], betaPrior[,2], log = TRUE)


    if(no_beta > 1){
      betaProposalPrior = sum(mapply(dgamma, betaProposal, shape = betaPrior[,1], rate = betaPrior[,2],
                                     MoreArgs = list(log = TRUE)))
    } else{
      betaProposalPrior = dgamma(betaProposal, betaPrior[1], betaPrior[2], log = TRUE)
    }


    logAlpha = (logLikelihoodProp + betaProposalPrior + sum(logBetaCurr)) -
            (logLikelihoodCurr + betaCurrPrior + sum(logBetaProposal))
    print(i)
    print(betaCurr)
    if(i == 3334){
       #print("ERROR HERE")
    }
    if(log_u < logAlpha){
      betaCurr = betaProposal
      logLikelihoodCurr = logLikelihoodProp
      betaCurrprior = betaProposalPrior
      acceptCounter = acceptCounter + 1
    }

  gammaCurr = rgamma(1, gammaPosterior[1], gammaPosterior[2])
  results[i + 1, ] = c(betaCurr, gammaCurr)
  }




  End <- as.numeric(Sys.time())

  timeTaken <- End - Start

  # Thin the samples
  results <- results[seq(from = burnIn + 1, to = (noIts + 1) - burnIn, by = thinFactor),]

  # Calculate R_0 sample values
  #R0_samples = (N*draws[,1])/draws[,2]

  # Calculate Effective Sample Sizes (and Per Second) and Acceptance Rates
  ESS <- coda::effectiveSize(results[, 1:no_beta])
  ESS_sec <- ESS/timeTaken
  acceptRate <-  acceptCounter/noIts

  # = Plots =
  par(mfrow = c(2, no_beta))

  #' MCMC Plots
  for(i in 1:no_beta){
    MCMCplots(results[, i], lagMax, ylab = expression(beta[as.expression(i)]))
  }

  return(list(timeTaken = timeTaken, results = results, gammaPosterior = gammaPosterior,
              ESS_sec = ESS_sec))
}


