
# SIS fsMCMC with a more pragmatic independence sampler.

# U and E and organised in sequential blocks of blockSize. Sample one member of each bloack to be updated
# Tune k for approx. 23% accept rate.

# ==== SIS fsMCMC Blocked Independence Sampler ====

SIS_fsMCMC_blockedIS = function(obsTransData, I_0, obsTimes, N, beta0, gamma0, thetaLim, lambda, V, noDraws, blockSize, noIts,
                                burnIn = 0, lagMax = NA, thinningFactor = 1){

  #' Calculate the amount of blockSamples which need to be made
  #' To accomodate for a potentially different final block size (blockSamples is not integer), calculate the residual noDraws
  blockSamples = ceiling(noDraws/blockSize)
  indicies = lapply(X = 0:(blockSamples - 1), function(X){ c(X*blockSize + 1, min((X+1)*blockSize, noDraws))})

  Start = as.numeric(Sys.time())


  noSampled = sum(obsTransData[[1]])
  thetaCurr = c(beta0, gamma0)
  logPCurr = -Inf

  saved = noDraws
  while(logPCurr == -Inf){
    ECurr = rexp(noDraws)
    UCurr = runif(noDraws)
    sim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr[1],
                                              thetaCurr[2], obsTimes,
                                              ECurr, UCurr)
    if(!is.list(sim)){
      noDraws = noDraws + 1000
    } else{
      panelDataSim = sim$panelData
      transDataSim = transitionData(panelDataSim, states = 1:2)
      logPCurr = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)
    }
  }

  if(noDraws > saved){
    blockSamples = ceiling(noDraws/blockSize)
    indicies = lapply(X = 0:(blockSamples - 1), function(X){ c(X*blockSize + 1, min((X+1)*blockSize, noDraws))})
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

    #' Folded Normal
    thetaProp = abs(thetaCurr + mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = lambda*V))


    newSim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaProp[1],
                                                 thetaProp[2], obsTimes, ECurr, UCurr)

    save = noDraws
    while(!is.list(newSim)){
      noDraws = noDraws + 1000
      UCurr = c(UCurr, runif(1000))
      ECurr = c(ECurr, rexp(1000))
      newSim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaProp[1],
                                                   thetaProp[2], obsTimes, ECurr, UCurr)
      print("Number of Draws Increased")
    }

    if(noDraws > save){
      blockSamples = ceiling(noDraws/blockSize)
      indicies = lapply(X = 0:(blockSamples - 1), function(X){ c(X*blockSize + 1, min((X+1)*blockSize, noDraws))})
    }
    transDataSim = transitionData(newSim$panelData, states = 1:2)

    logPProp = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)

    # + sum(thetaCurr)
    if(sum(thetaProp < thetaLim) == 2){
      logA = (logPProp + sum(dexp(thetaProp, rate = 0.001, log = T)) ) -
        (logPCurr + sum(dexp(thetaCurr, rate = 0.001, log = T)) )
    } else{
      logA = -Inf
    }
    accProbs[i] = exp(logA)
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

    newSim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr[1],
                                                       thetaCurr[2], obsTimes, EProp, UProp)

    save = noDraws
    while(!is.list(newSim)){
      noDraws = noDraws + 1000
      UProp = c(UProp, runif(1000))
      EProp = c(EProp, rexp(1000))
      newSim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr[1],
                                                   thetaCurr[2], obsTimes, EProp, UProp)
      print("Number of Draws Increased")
    }

    if(noDraws > save){
      blockSamples = ceiling(noDraws/blockSize)
      indicies = lapply(X = 0:(blockSamples - 1), function(X){ c(X*blockSize + 1, min((X+1)*blockSize, noDraws))})
    }
    transDataSim = transitionData(newSim$panelData, states = 1:2)

    logPProp = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)

    logA = (logPProp) - (logPCurr)
    logU = log(runif(1))

    if(logU < logA){
      logPCurr = logPProp
      ECurr = EProp
      UCurr = UProp
      acceptEU = acceptEU + 1
    }

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


  printed_output(rinf_dist = "Exp", no_proposals = NA, noIts, ESS, timeTaken, ESS.sec, acceptRate)
  print(mean(accProbs))
  print(sum(accProbs > 1))
  return(list(draws = draws, acceptRate = acceptRate, ESS = ESS, ESS.sec = ESS.sec,
              betaSummary = betaSummary, gammaSummary = gammaSummary, timeTaken = timeTaken))


}

# ==== Adaptive MCMC ====

# Adapt RWM Jump Proposal Parameter \lambda and Covariance Matrix V

# Find a sufficient amount of draws

# Adapt blockSize parameter which changes the number of RVs changed


adaptiveSISfsMCMC7_blockedIS = function(obsTransData, I_0, obsTimes, N, beta0, gamma0, thetaLim, lambda0, V, firstRun = TRUE, noDraws,
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
  noSampled = sum(obsTransData[[1]])
  thetaCurr = c(beta0, gamma0)
  logPCurr = -Inf

  saved = noDraws
  while(logPCurr == -Inf){
    ECurr = rexp(noDraws)
    UCurr = runif(noDraws)
    sim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr[1],
                                              thetaCurr[2], obsTimes,
                                              ECurr, UCurr)
    if(!is.list(sim)){
      noDraws = noDraws + 1000
    } else{
      panelDataSim = sim$panelData
      transDataSim = transitionData(panelDataSim, states = 1:2)
      logPCurr = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)
    }
  }

  if(noDraws > saved){
    blockSamples = ceiling(noDraws/blockSize)
    indicies = lapply(X = 0:(blockSamples - 1), function(X){ c(X*blockSize + 1, min((X+1)*blockSize, noDraws))})
  }

  maxNoEvents = sim$noEvents

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
      Vi = var(draws[,1:2], na.rm = T)
      thetaProp = abs(thetaCurr + lambda*mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = Vi))
    } else{
      thetaProp = abs(thetaCurr + lambda0*mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = V))
    }

    newSim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaProp[1],
                                                 thetaProp[2], obsTimes, ECurr, UCurr)

    save = noDraws
    while(!is.list(newSim)){
      noDraws = noDraws + 1000
      UCurr = c(UCurr, runif(1000))
      ECurr = c(ECurr, rexp(1000))
      newSim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaProp[1],
                                                 thetaProp[2], obsTimes, ECurr, UCurr)
    }

    if(noDraws > save){
      blockSamples = ceiling(noDraws/blockSize)
      indicies = lapply(X = 0:(blockSamples - 1), function(X){ c(X*blockSize + 1, min((X+1)*blockSize, noDraws))})
    }

    maxNoEvents = max(newSim$noEvents, maxNoEvents)
    transDataSim = transitionData(newSim$panelData, states = 1:2)

    logPProp = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)

    if(sum(thetaProp < thetaLim) == 2){
      logA = (logPProp + sum(dexp(thetaCurr, rate = 0.001, log = T)) ) -
        (logPCurr + sum(dexp(thetaProp, rate = 0.001, log = T)) )
    } else{
      logA = -Inf
    }
    accProbs[i] = exp(logA)
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


    newSim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr[1],
                                                       thetaCurr[2], obsTimes, EProp, UProp)

    save = noDraws
    while(!is.list(newSim)){
      noDraws = noDraws + 1000
      UProp = c(UProp, runif(1000))
      EProp = c(EProp, rexp(1000))
      newSim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr[1],
                                                   thetaCurr[2], obsTimes, EProp, UProp)
      print(noDraws)
    }

    if(noDraws > save){
      blockSamples = ceiling(noDraws/blockSize)
      indicies = lapply(X = 0:(blockSamples - 1), function(X){ c(X*blockSize + 1, min((X+1)*blockSize, noDraws))})
    }
    maxNoEvents = max(newSim$noEvents, maxNoEvents)
    transDataSim = transitionData(newSim$panelData, states = 1:2)

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
  print(c("No. Draws Unused: ", noDraws - maxNoEvents))
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

  return(list(draws = draws, acceptRate = acceptRate, lambda = lambda, V = Vi, noDraws = noDraws, blockSize = blockSize, ESS = ESS, ESS.sec = ESS.sec,
              betaSummary = betaSummary, gammaSummary = gammaSummary, timeTaken = timeTaken))
}
