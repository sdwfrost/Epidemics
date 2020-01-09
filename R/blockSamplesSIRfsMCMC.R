
#' fsMCMC MCMC for SIR Epidemic Panel Data (with blocked proposals of underlying Random Variables)
#'
#' Adapts proposal parameters for with a view of optimal of a target using forward simulation MCMC scheme.

#' @family Panel Data MCMC
#' @param obsTransData Interpanel transition data.
#' @param I_0 Initial number of infectives in the population.
#' @param obsTimes Times at which epidemic cohort were followed up.
#' @param N Population size.
#' @param beta0 Starting value for infectious process parameter.
#' @param gamma0 Starting value for removal/recovery process parameter.
#' @param lambda Starting value for RWM proposal parameter which is to be adapted.
#' @param V Starting state for RWM proposal Covariance matrix which is to be adapted.
#' @param noDraws Number of underlying variables to be drawn (i.e max. number of events in observation period)
#' @param blockSize Number of underlying random variables of the process to refresh.
#' @param noIts Number of MCMC iterations.
#' @param lagMax Plotting parameter for acf() function.
#' @param thinningFactor Controls the factor by which MCMC samples are thinned, to reduce dependency.
#'
#' @return MCMC summary.

SIR_fsMCMC_blockedIS = function(obsTransData, I_0, obsTimes, N, beta0, gamma0, lambda, V, noDraws = 2*N - I_0, blockSize, noIts,
                                burnIn = 0, lagMax = NA, thinningFactor = 1){

  # Calculate the amount of blockSamples which need to be made
  # To accomodate for a potentially different final block size (blockSamples is not integer), calculate the residual noDraws
  blockSamples = ceiling(noDraws/blockSize)
  indicies = lapply(X = 0:(blockSamples - 1), function(X){c(X*blockSize + 1, min((X+1)*blockSize, noDraws))})

  Start = as.numeric(Sys.time())


  noSampled = sum(obsTransData[[1]])
  thetaCurr = c(beta0, gamma0)
  logPCurr = -Inf

  while(logPCurr == -Inf){
    ECurr = rexp(noDraws)
    UCurr = runif(noDraws)
    sim = homogeneousPanelDataSIR_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr[1],
                                              thetaCurr[2], obsTimes,
                                              ECurr, UCurr)
    panelDataSim = sim$panelData
    transDataSim = transitionData(panelDataSim, states = 1:3)
    logPCurr = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)

  }

  # Create Storage Matrix
  draws = matrix(NA, nrow = noIts + 1, ncol = length(thetaCurr) + 1)
  draws[1,] = c(thetaCurr, logPCurr)

  # Proposal Acceptance Counter
  acceptTheta = 0
  acceptEU = 0
  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = noIts)

  for(i in 1:noIts){
    pb$tick()
    # ==== Beta and Gamma Proposal ====

    # Folded Normal
    thetaProp = abs(thetaCurr + mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = lambda*V))


    newSim = homogeneousPanelDataSIR_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaProp[1],
                                                 thetaProp[2], obsTimes, ECurr, UCurr)
    transDataSim = transitionData(newSim$panelData, states = 1:3)
    logPProp = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)

    logA = (logPProp + sum(dexp(thetaProp, rate = 0.001, log = T)) ) -
      (logPCurr + sum(dexp(thetaCurr, rate = 0.001, log = T)) )

    logU = log(runif(1))

    if(logU < logA){
      logPCurr = logPProp
      thetaCurr = thetaProp
      acceptTheta = acceptTheta + 1
    }

    # ==== E and U proposal ====
    proposalSet = sapply(X = indicies, function(X) sample(X[1]:X[2], size = 1))

    EProp = ECurr
    UProp = UCurr

    EProp[proposalSet] = rexp(blockSamples)
    UProp[proposalSet] = runif(blockSamples)

    newSim = homogeneousPanelDataSIR_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr[1],
                                                 thetaCurr[2], obsTimes, EProp, UProp)
    transDataSim = transitionData(newSim$panelData, states = 1:3)
    logPProp = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)

    logA = (logPProp) - (logPCurr)
    logU = log(runif(1))

    if(logU < logA){
      logPCurr = logPProp
      ECurr = EProp
      UCurr = UProp
      acceptEU = acceptEU + 1
    }

    # Store State
    draws[i+1, ] = c(thetaCurr, logPCurr)
  }

  End <- as.numeric(Sys.time())

  timeTaken <- End - Start

  # Thin the samples
  draws <- draws[seq(from = burnIn + 1, to = (noIts + 1) - burnIn, by = thinningFactor),]

  # Calculate Effective Sample Sizes (and Per Second) and Acceptance Rates
  ESS <- coda::effectiveSize(draws[,1:2])
  ESS.sec <- ESS/timeTaken
  acceptRate <-  c(acceptTheta, acceptEU)/noIts

  # = Plots =
  par(mfrow = c(2,2))

  # Plot Beta Samples and Sample Auto-Corrolation Function
  if(is.na(lagMax)){
    # Beta
    plot(draws[, 1], type = 'l', ylab = expression(beta))
    acf(draws[, 1], main = "")
    # Gamma
    plot(draws[, 2], type = 'l', ylab = expression(gamma))
    acf(draws[, 2], main = "")
  } else{
    # Beta
    plot(draws[, 1], type = 'l', ylab = expression(beta))
    acf(draws[, 1], lagMax, main = "")
    # Gamma
    plot(draws[, 2], type = 'l', ylab = expression(gamma))
    acf(draws[, 2], lagMax, main = "")
  }

  # Calculating Summary Statistics for samples
  betaSummary = c(mean(draws[,1]), sd(draws[,1]))
  gammaSummary = c(mean(draws[,2]), sd(draws[,2]))


  printed_output(rinf_dist = "Exp", no_proposals = NA, noIts, ESS, timeTaken, ESS.sec, acceptRate)
  return(list(draws = draws, acceptRate = acceptRate, ESS = ESS, ESS.sec = ESS.sec,
              betaSummary = betaSummary, gammaSummary = gammaSummary, timeTaken = timeTaken))


}


# ==== Adaptive MCMC ====

# Adapt RWM Jump Proposal Parameter \lambda and Covariance Matrix V

# Find a sufficient amount of draws

# Adapt blockSize parameter which changes the number of RVs changed


adaptiveSIRfsMCMC7_blockedIS = function(obsTransData, I_0, obsTimes, N, beta0, gamma0, thetaLim, lambda0, V,
                                        firstRun = TRUE, noDraws = 2*N - I_0,
                                        blockSize, noIts, burnIn = 0, lagMax = NA, thinningFactor = 1, delta){

  Start = as.numeric(Sys.time())

  #' Calculate the amount of blockSamples which need to be made
  #' To accomodate for a potentially different final block size (blockSamples is not integer), calculate the residual noDraws
  blockSamples = ceiling(noDraws/blockSize)
  indicies = lapply(X = 0:(blockSamples - 1), function(X){ c(X*blockSize + 1, min((X+1)*blockSize, noDraws))})

  if(missing(V)){
    V = diag(c(1/N, 1))
  }

  if(firstRun){
    lambda = 2.38/sqrt(2)
    V = diag(c(1/N, 1))
  } else{
    lambda = lambda0
  }
  Vi = diag(c(1/N, 1))
  noSampled = sum(obsTransData[[1]])
  thetaCurr = c(beta0, gamma0)
  logPCurr = -Inf
  while(logPCurr == -Inf){
    ECurr = rexp(noDraws)
    UCurr = runif(noDraws)
    sim = homogeneousPanelDataSIR_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr[1],
                                              thetaCurr[2], obsTimes,
                                              ECurr, UCurr)

    panelDataSim = sim$panelData
    transDataSim = transitionData(panelDataSim, states = 1:3)
    logPCurr = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)
    }


  #' Create Storage Matrix
  draws = matrix(NA, nrow = noIts + 1, ncol = length(thetaCurr) + 1)
  draws[1,] = c(thetaCurr, logPCurr)

  #' Proposal Acceptance Counter
  acceptTheta = 0
  accProbs = c()
  acceptEU = 0
  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = noIts)

  for(i in 1:noIts){
    pb$tick()
    # ==== Beta and Gamma Proposal ====
    u1 = runif(1, 0, 1)
    if(u1 > delta & acceptTheta > 10){
      # save  = Vi
      Vi = var(draws[,1:2], na.rm = T)
      # print(sum(save - Vi))
      thetaProp = abs(thetaCurr + lambda*mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = Vi))
    } else{
      thetaProp = abs(thetaCurr + lambda0*mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = V))
    }

    newSim = homogeneousPanelDataSIR_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaProp[1],
                                                 thetaProp[2], obsTimes, ECurr, UCurr)
    transDataSim = transitionData(newSim$panelData, states = 1:3)
    logPProp = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)

    if(sum(thetaProp < thetaLim) == 2){
      logA = (logPProp + sum(dexp(thetaCurr, rate = 0.001, log = T)) ) -
        (logPCurr + sum(dexp(thetaProp, rate = 0.001, log = T)) )
    } else{
      logA = -Inf
    }

    logU = log(runif(1))

    if(logU < logA){
      logPCurr = logPProp
      thetaCurr = thetaProp
      acceptTheta = acceptTheta + 1
      if(u1 > delta){
        lambda = lambda + 0.93*(lambda/(sqrt(i)))
      }
    } else{
      if(u1 > delta){
        lambda = lambda - 0.07*(lambda/(sqrt(i)))
      }
    }
    # ==== E and U proposal ====
    proposalSet = sapply(X = indicies, function(X) sample(X[1]:X[2], size = 1))

    EProp = ECurr
    UProp = UCurr

    EProp[proposalSet] = rexp(blockSamples)
    UProp[proposalSet] = runif(blockSamples)


    newSim = homogeneousPanelDataSIR_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr[1],
                                                 thetaCurr[2], obsTimes, EProp, UProp)

    transDataSim = transitionData(newSim$panelData, states = 1:3)

    logPProp = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)

    logA = (logPProp) - (logPCurr)
    logU = log(runif(1))

    if(logU < logA){
      logPCurr = logPProp
      ECurr = EProp
      UCurr = UProp
      acceptEU = acceptEU + 1
      if(blockSize > 3){
        blockSize = blockSize - 3*rbinom(1, size = 1, prob = 1/i^(0.25))
      }
    } else{
      blockSize = blockSize + 1*rbinom(1, size = 1, prob = 1/i^(0.25))
    }

    #' Calculate the amount of blockSamples which need to be made
    #' To accomodate for a potentially different final block size (blockSamples is not integer), calculate the residual noDraws
    blockSamples = ceiling(noDraws/blockSize)
    indicies = lapply(X = 0:(blockSamples - 1), function(X){c(X*blockSize + 1, min((X+1)*blockSize, noDraws))})



    #' Store State
    draws[i+1, ] = c(thetaCurr, logPCurr)
  }

  End <- as.numeric(Sys.time())

  timeTaken <- End - Start

  # Thin the samples
  draws <- draws[seq(from = burnIn + 1, to = (noIts + 1) - burnIn, by = thinningFactor),]

  # Calculate Effective Sample Sizes (and Per Second) and Acceptance Rates
  ESS <- coda::effectiveSize(draws[,1:2])
  ESS.sec <- ESS/timeTaken
  acceptRate <-  c(acceptTheta, acceptEU)/noIts

  printed_output(rinf_dist = "Exp", no_proposals = NA, noIts, ESS, timeTaken, ESS.sec, acceptRate)

  # = Plots =
  par(mfrow = c(2,2))

  # Plot Beta Samples and Sample Auto-Corrolation Function
  if(is.na(lagMax)){
    # Beta
    plot(draws[, 1], type = 'l', ylab = expression(beta))
    acf(draws[, 1], main = "")
    # Gamma
    plot(draws[, 2], type = 'l', ylab = expression(gamma))
    acf(draws[, 2], main = "")
  } else{
    # Beta
    plot(draws[, 1], type = 'l', ylab = expression(beta))
    acf(draws[, 1], lagMax, main = "")
    # Gamma
    plot(draws[, 2], type = 'l', ylab = expression(gamma))
    acf(draws[, 2], lagMax, main = "")
  }

  #' Calculating Summary Statistics for samples
  betaSummary = c(mean(draws[,1]), sd(draws[,1]))
  gammaSummary = c(mean(draws[,2]), sd(draws[,2]))

  return(list(draws = draws, acceptRate = acceptRate, lambda = lambda, V = Vi, blockSize = blockSize, ESS = ESS, ESS.sec = ESS.sec,
              betaSummary = betaSummary, gammaSummary = gammaSummary, timeTaken = timeTaken))
}


