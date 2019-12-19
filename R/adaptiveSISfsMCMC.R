# ==== Adaptive fsMCMC SIS, aim for 25% Accept Rate ====

adaptiveSISfsMCMC25 = function(obsTransData, I_0, obsTimes, N, beta0, gamma0, thetaLim, lambda0, V, firstRun = TRUE, noDraws, s, noIts,
                      burnIn = 0, lagMax = NA, thinningFactor = 1, delta){

  Start = as.numeric(Sys.time())

  if(missing(V)){
    V = diag(c(1/N, 1))
  }

  if(firstRun){
    lambda = 2.38/sqrt(2)
    #lambda = lambda1
    V = diag(c(1/N, 1))
  } else{
    lambda = lambda0
  }
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
    u1 = runif(1, 0, 1)
    if(u1 > delta & acceptTheta > 10){
      thetaProp = abs(thetaCurr + lambda*mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = var(draws[,1:2], na.rm = T)))
    } else{
      thetaProp = abs(thetaCurr + lambda0*mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = V))
    }
    #' Folded Normal
    #thetaProp = abs(thetaCurr + lambda*mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = V))
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
      if(u1 > delta){
        lambda = lambda + 3*(lambda/(sqrt(i)))
      }
    } else{
      if(u1 > delta){
        lambda = lambda - 1*(lambda/(sqrt(i)))
      }
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

  printed_output(rinf_dist = "Exp", no_proposals = NA, noIts, ESS, timeTaken, ESS.sec, acceptRate)
  if(acceptRate[1] > 0.25 | acceptRate[1] < 0.20){
    print("Further Adaptation Advised")
  } else{
    print("Jump Size Adequate")
  }
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



  return(list(draws = draws, acceptRate = acceptRate, lambda = lambda, V = V, ESS = ESS, ESS.sec = ESS.sec,
              betaSummary = betaSummary, gammaSummary = gammaSummary, timeTaken = timeTaken))


}


# ==== Adaptive fsMCMC SIS, aim for 7% Accept Rate ====


adaptiveSISfsMCMC7 = function(obsTransData, I_0, obsTimes, N, beta0, gamma0, thetaLim, lambda0, V, firstRun = TRUE, noDraws, s, noIts,
                               burnIn = 0, lagMax = NA, thinningFactor = 1, delta){

  Start = as.numeric(Sys.time())

  if(missing(V)){
    V = diag(c(1/N, 1))
  }

  if(firstRun){
    lambda = 2.38/sqrt(2)
    #lambda = lambda1
    V = diag(c(1/N, 1))
  } else{
    lambda = lambda0
  }
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
    #' Folded Normal
    #thetaProp = abs(thetaCurr + lambda*mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = V))
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
      if(u1 > delta){
        lambda = lambda + 0.93*(lambda/(sqrt(i)))
      }
    } else{
      if(u1 > delta){
        lambda = lambda - 0.07*(lambda/(sqrt(i)))
      }
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

  printed_output(rinf_dist = "Exp", no_proposals = NA, noIts, ESS, timeTaken, ESS.sec, acceptRate)
  if(acceptRate[1] > 0.25 | acceptRate[1] < 0.20){
    print("Further Adaptation Advised")
  } else{
    print("Jump Size Adequate")
  }
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



  return(list(draws = draws, acceptRate = acceptRate, optLambda = lambda, optV = Vi, ESS = ESS, ESS.sec = ESS.sec,
              betaSummary = betaSummary, gammaSummary = gammaSummary, timeTaken = timeTaken))


}
