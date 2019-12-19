


oneParticlePM_SIS_beta = function(obsTransData, obsTimes, N, beta0, gamma, lambda, noIts,
                            burnIn = 0, lagMax = NA, thinningFactor = 1){

  Start = as.numeric(Sys.time())

  noSampled = sum(obsTransData[[1]])
  betaCurr = beta0
  logPCurr = -Inf

  while(logPCurr == -Inf){
    panelDataSim = homogeneousPanelDataSIS_Gillespie(initialState = c(rep(1, N - 1), 2), betaCurr, gamma,
                                                     obsTimes)$panelData
    transDataSim = transitionData(panelDataSim, states = 1:2)
    logPCurr = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)
  }

  #' Create Storage Matrix
  draws = matrix(NA, nrow = noIts + 1, ncol = length(betaCurr) + 1)
  draws[1,] = c(betaCurr, logPCurr)

  #' Proposal Acceptance Counter
  acceptBeta = 0
  print("Sampling Progress")
  pb <- progress::progress_bar$new(total = noIts)

  for(i in 1:noIts){
    pb$tick()
    #' Propose new beta and gamma using Multiplicative RW propsal
    logBetaCurr = log(betaCurr)
    logBetaProp = logBetaCurr + lambda*rnorm(1, 0, 1)
    betaProp = exp(logBetaProp)

    panelDataSim = homogeneousPanelDataSIS_Gillespie(initialState = c(rep(1, N - 1), 2), betaProp, gamma,
                                                     obsTimes)$panelData
    transDataSim = transitionData(panelDataSim, states = 1:2)

    logPProp = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)

    loga = (logPProp + sum(betaCurr)) - (logPCurr + sum(betaProp))

    logu = log(runif(1))

    if(logu < loga){
      logPCurr = logPProp
      betaCurr = betaProp
      acceptBeta = acceptBeta + 1
    }

    #' Store State
    draws[i+1, ] = c(betaCurr, logPCurr)
  }

  End <- as.numeric(Sys.time())

  timeTaken <- End - Start

  # Thin the samples
  draws <- draws[seq(from = burnIn + 1, to = (noIts + 1) - burnIn, by = thinningFactor),]

  # Calculate Effective Sample Sizes (and Per Second) and Acceptance Rates
  ESS <- coda::effectiveSize(draws[,1])
  ESS.sec <- ESS/timeTaken
  acceptRate <-  acceptBeta/noIts

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



