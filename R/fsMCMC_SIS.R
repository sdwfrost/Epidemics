#'


fsMCMC_SIS = function(obsTransData, obsTimes, N, beta0, gamma0, lambda, s, noIts,
                             burnIn = 0, lagMax = NA, thinningFactor = 1){

  Start = as.numeric(Sys.time())

  noSampled = sum(obsTransData[[1]])
  thetaCurr = c(beta0, gamma0)
  uCurr = runif()
  logPCurr = -Inf


  while(logPCurr == -Inf){
    panelDataSim = homogeneousPanelDataSIS_Gillespie(initialState = c(rep(1, N - 1), 2), thetaCurr[1],
                                                     thetaCurr[2],
                                                     obsTimes)$panelData
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
    #' Propose new beta and gamma using Multiplicative RW propsal
    logThetaCurr = log(thetaCurr)
    logThetaProp = logThetaCurr + mvtnorm::rmvnorm(1, mean = rep(0, 2), sigma = lambda*diag(1, 2))
    thetaProp = exp(logThetaProp)

    panelDataSim = homogeneousPanelDataSIS_Gillespie(initialState = c(rep(1, N - 1), 2), thetaProp[1],
                                                     thetaProp[2], obsTimes)$panelData
    transDataSim = transitionData(panelDataSim, states = 1:2)

    logPProp = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)

    loga = (logPProp + sum(thetaCurr)) - (logPCurr + sum(thetaProp))

    logu = log(runif(1))

    if(logu < loga){
      logPCurr = logPProp
      thetaCurr = thetaProp
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
    plot(draws[, 1], type = 'l', ylab = expression(gamma))
    acf(draws[, 1], main = "")
  } else{
    # Beta
    plot(draws[, 1], type = 'l', ylab = expression(beta))
    acf(draws[, 1], lagMax, main = "")
    # Gamma
    plot(draws[, 1], type = 'l', ylab = expression(gamma))
    acf(draws[, 1], lagMax, main = "")
  }

  #' Calculating Summary Statistics for samples
  betaSummary = c(mean(draws[,1]), sd(draws[,1]))
  gammaSummary = c(mean(draws[,2]), sd(draws[,2]))

  printed_output(rinf_dist = "Exp", no_proposals = NA, noIts, ESS, timeTaken, ESS.sec, acceptRate)
  return(list(draws = draws, acceptRate = acceptRate, ESS = ESS, ESS.sec = ESS.sec,
              betaSummary = betaSummary, gammaSummary = gammaSummary, timeTaken = timeTaken))

}
