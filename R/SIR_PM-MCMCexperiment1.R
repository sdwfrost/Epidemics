#' SIR Test Epidemic 1 PM-MCMC Experiment


SIR_PM_MCMCexperiment1 = function(noPanels, noSims, lambda0, V0, adapt = T, noIts = 1e6){
  #' Test data
  N = 200
  I_0 = 0.05*N
  gamma = 1
  R0 = 1.4
  beta = gamma*R0/N

  obsTimes = seq(0, 5, length = noPanels)

  set.seed(1)
  testData = homogeneousPanelDataSIR_GillespieEU(c(rep(1, N - I_0), rep(2, I_0)), beta, gamma, obsTimes, E = rexp(2*N - I_0),
                                                 U = runif(2*N - I_0))
  testDataSample = panelDataSample(testData$panelData, m = 0.1*N)
  Y = transitionData(testDataSample, states = 1:3)

  if(adapt){
    print(paste(c("==== Adaptive Step (", noPanels, "Panels,", noSims, "Sims )",   "===="), collapse = " "))

    adaptRun = adaptiveSIR_PseudoMarginalMCMC(Y, I_0, obsTimes, N, beta, gamma, lambda0, V0,
                                              noSims, noIts = 10000, burnIn = 0)
    lambda = adaptRun$lambda
    V = adaptRun$V
  } else{
    adaptRun = NULL
    lambda = lambda0
    V = V0
  }

  print(paste(c("==== MCMC Step (", noPanels, "Panels,", noSims, "Sims )",   "===="), collapse = " "))
  MCMCrun = SIR_PseudoMarginalMCMC(Y, I_0, obsTimes, N, beta, gamma, lambda = lambda, V = V,
                                    noSims, noIts = noIts, burnIn = 0)

  return(list(MCMCrun = MCMCrun, adaptRun = adaptRun))
}

