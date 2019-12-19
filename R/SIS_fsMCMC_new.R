

# ==== SIS_fsMCMC_beta ====
SIS_fsMCMC_beta = function(obsTransData, I_0, obsTimes, N, beta0, gamma, thetaLim, lambda, noDraws, s, noIts,
                      burnIn = 0, lagMax = NA, thinningFactor = 1){

  Start = as.numeric(Sys.time())

  noSampled = sum(obsTransData[[1]])
  thetaCurr = beta0
  logPCurr = -Inf

  while(logPCurr == -Inf){
    ECurr = rexp(noDraws)
    UCurr = runif(noDraws)
    sim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr[1],
                                              gamma, obsTimes,
                                              ECurr, UCurr)
    panelDataSim = sim$panelData
    transDataSim = transitionData(panelDataSim, states = 1:2)
    logPCurr = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)
  }

  #' Create Storage Matrix
  draws = matrix(NA, nrow = noIts + 1, ncol = length(thetaCurr) + 1)
  draws[1,] = c(thetaCurr, logPCurr)

  #' Proposal Acceptance Counter
  acceptTheta = 0
  acceptEU = 0
  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = noIts)

  for(i in 1:noIts){
    pb$tick()
    # ==== Beta and Gamma Proposal ====

    #' Folded Normal
    thetaProp = abs(thetaCurr + rnorm(1, 0, lambda))
    # logThetaCurr = log(thetaCurr)
    # logThetaProp = logThetaCurr + mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = lambda*V)
    # thetaProp = exp(logThetaProp)

    newSim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaProp,
                                                 gamma, obsTimes, ECurr, UCurr)

    transDataSim = transitionData(newSim$panelData, states = 1:2)

    logPProp = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)

    # + sum(thetaCurr)
    if(sum(thetaProp < thetaLim) == length(thetaCurr)){
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
    }

    # ==== E and U proposal ====
    proposalSet = sample(1:noDraws, size = s, replace = F)

    EProp = ECurr
    UProp = UCurr

    EProp[proposalSet] = rexp(s)
    UProp[proposalSet] = runif(s)

    panelDataSim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr,
                                                       gamma, obsTimes, EProp, UProp)$panelData
    transDataSim = transitionData(panelDataSim, states = 1:2)

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
  par(mfrow = c(1,2))

  # Plot Beta Samples and Sample Auto-Corrolation Function
  if(is.na(lagMax)){
    # Beta
    plot(draws[, 1], type = 'l', ylab = expression(beta))
    acf(draws[, 1], main = "")
  } else{
    # Beta
    plot(draws[, 1], type = 'l', ylab = expression(beta))
    acf(draws[, 1], lagMax, main = "")
  }

  #' Calculating Summary Statistics for samples
  betaSummary = c(mean(draws[,1]), sd(draws[,1]))

  printed_output(rinf_dist = "Exp", no_proposals = NA, noIts, ESS, timeTaken, ESS.sec, acceptRate)
  return(list(draws = draws, acceptRate = acceptRate, ESS = ESS, ESS.sec = ESS.sec,
              betaSummary = betaSummary, timeTaken = timeTaken))


}


# ==== SIS_fsMCMC_oneProposal ====
SIS_fsMCMC_oneProposal = function(obsTransData, I_0, obsTimes, N, beta0, gamma0, thetaLim, lambda, V, noDraws, s, noIts,
                                  burnIn = 0, lagMax = NA, thinningFactor = 1){

  Start = as.numeric(Sys.time())

  noSampled = sum(obsTransData[[1]])
  thetaCurr = c(beta0, gamma0)
  logPCurr = -Inf

  while(logPCurr == -Inf){
    ECurr = rexp(noDraws)
    UCurr = runif(noDraws)
    sim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr[1],
                                                     thetaCurr[2], obsTimes,
                                                     ECurr, UCurr)
    panelDataSim = sim$panelData
    transDataSim = transitionData(panelDataSim, states = 1:2)
    logPCurr = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)
  }

  #' Create Storage Matrix
  draws = matrix(NA, nrow = noIts + 1, ncol = length(thetaCurr) + 1)
  draws[1,] = c(thetaCurr, logPCurr)

  #' Proposal Acceptance Counter
  acceptTheta = 0
  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = noIts)

  for(i in 1:noIts){
    pb$tick()
    # ==== Beta and Gamma Proposal ====

    #' Folded Normal
    thetaProp = abs(thetaCurr + mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = lambda*V))
    # logThetaCurr = log(thetaCurr)
    # logThetaProp = logThetaCurr + mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = lambda*V)
    # thetaProp = exp(logThetaProp)

    # ==== E and U proposal ====
    proposalSet = sample(1:noDraws, size = s, replace = F)

    EProp = ECurr
    UProp = UCurr

    EProp[proposalSet] = rexp(s)
    UProp[proposalSet] = runif(s)


    newSim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaProp[1],
                                                     thetaProp[2], obsTimes, ECurr, UCurr)

    transDataSim = transitionData(newSim$panelData, states = 1:2)

    logPProp = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)

    # + sum(thetaCurr)
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
      UCurr = UProp
      ECurr = EProp
      acceptTheta = acceptTheta + 1
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
  acceptRate <-  acceptTheta/noIts

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

  return(list(draws = draws, acceptRate = acceptRate, ESS = ESS, ESS.sec = ESS.sec,
              betaSummary = betaSummary, gammaSummary = gammaSummary, timeTaken = timeTaken))

}



# ==== SIS_fsMCMC ====

SIS_fsMCMC = function(obsTransData, I_0, obsTimes, N, beta0, gamma0, thetaLim, lambda, V, noDraws, s, noIts,
                      burnIn = 0, lagMax = NA, thinningFactor = 1){

  Start = as.numeric(Sys.time())

  noSampled = sum(obsTransData[[1]])
  thetaCurr = c(beta0, gamma0)
  logPCurr = -Inf

  while(logPCurr == -Inf){
    ECurr = rexp(noDraws)
    UCurr = runif(noDraws)
    sim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr[1],
                                              thetaCurr[2], obsTimes,
                                              ECurr, UCurr)
    panelDataSim = sim$panelData
    transDataSim = transitionData(panelDataSim, states = 1:2)
    logPCurr = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)
  }

  #' Create Storage Matrix
  draws = matrix(NA, nrow = noIts + 1, ncol = length(thetaCurr) + 1)
  draws[1,] = c(thetaCurr, logPCurr)

  #' Proposal Acceptance Counter
  acceptTheta = 0
  accProbs = c()
  acceptEU = 0
  #print("Sampling Progress")
  #pb <- progress::progress_bar$new(total = noIts)

  for(i in 1:noIts){
    #pb$tick()
    # ==== Beta and Gamma Proposal ====

    #' Folded Normal
    thetaProp = abs(thetaCurr + mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = lambda*V))
    # logThetaCurr = log(thetaCurr)
    # logThetaProp = logThetaCurr + mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = lambda*V)
    # thetaProp = exp(logThetaProp)

    newSim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaProp[1],
                                                 thetaProp[2], obsTimes, ECurr, UCurr)

    transDataSim = transitionData(newSim$panelData, states = 1:2)

    logPProp = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)

    # + sum(thetaCurr)
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
    }

    # ==== E and U proposal ====
    proposalSet = sample(1:noDraws, size = s, replace = F)

    EProp = ECurr
    UProp = UCurr

    EProp[proposalSet] = rexp(s)
    UProp[proposalSet] = runif(s)

    panelDataSim = homogeneousPanelDataSIS_GillespieEU(initialState = c(rep(1, N - I_0), rep(2, I_0)), thetaCurr[1],
                                                       thetaCurr[2], obsTimes, EProp, UProp)$panelData
    transDataSim = transitionData(panelDataSim, states = 1:2)

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

