#' Testing New blocked Independent Sampler fsMCMC for SIR (+ Adaptive Step)


#' SIR fsMCMC experiment

#' Multiple Datasets, Different speed outbreaks.

SIRfsMCMCexperiment1 = function(noPanels, lastObs, lambda0 = 1e-3, V0, blockSize, adapt = T, noIts = 1e6){

  #' Test data
  N = 200
  I_0 = 0.05*N
  gamma = 1
  R0 = 1.4
  beta = gamma*R0/N

  k = 5
  obsTimes = seq(0, 5, length = k)

  set.seed(1)
  testData = homogeneousPanelDataSIR_GillespieEU(c(rep(1, N - I_0), rep(2, I_0)), beta, gamma, obsTimes, E = rexp(2*N - I_0),
                                                 U = runif(2*N - I_0))

  testDataSample = panelDataSample(testData$panelData, m = 0.1*N)
  Y = transitionData(testDataSample, states = 1:3)

  if(adapt){
    print(paste(c("==== Adaptive Step (", noPanels, "Panels )",   "===="), collapse = " "))
    adaptRun = adaptiveSIRfsMCMC7_blockedIS(Y, I_0, obsTimes, N, beta, gamma, thetaLim = rep(Inf, 2), lambda0 = 1e-3, V = V0,
                                            blockSize = 2*N - I_0, noIts = 10000 , delta  = 0.05)

    lambda = adaptRun$lambda
    V = adaptRun$V
    BS = adaptRun$blockSize
  } else{
    adaptRun = NULL
    lambda = lambda0
    V = V0
    BS = blockSize
  }

  print(paste(c("==== MCMC Step (", noPanels, "Panels )",   "===="), collapse = " "))

  MCMCrun = SIR_fsMCMC_blockedIS(Y, I_0, obsTimes, N, beta, gamma, thetaLim = rep(Inf, 2), lambda = lambda, V = V,
                                  blockSize = BS, noIts = noIts)

  return(list(MCMCrun = MCMCrun, adaptRun = adaptRun))
}

