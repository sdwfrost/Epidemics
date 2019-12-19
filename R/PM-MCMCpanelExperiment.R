
#' PM-MCMC Panel Experiment


PM_MCMCexperiment = function(noSims, noPanels, lastObs, lambda0, V0 = diag(c(1/N, 1)), adapt = TRUE, noIts){

  # ==== Sample Panel Data ====
  set.seed(1)
  N = 200
  gamma = 1
  R_0 = 1.25
  beta = gamma*R_0/(N)
  I_0 = 40
  initialState = c(rep(1, N - I_0), rep(2, I_0))

  m = 0.1*N
  obsTimes = seq(0, lastObs, length = noPanels)
  SIS_sim = homogeneousPanelDataSIS_Gillespie(initialState, beta, gamma, obsTimes)

  #' Sample Panel Data
  panelData = panelDataSample(SIS_sim$panelData, m = m)
  Y = transitionData(panelData, states = 1:2)

  if(adapt){
    print(paste(c("==== Adaptive Step (", noPanels, "Panels,", noSims, "Sims )",   "===="), collapse = " "))
    adaptRun = adaptiveSIS_PseudoMarginalMCMC(Y, I_0, obsTimes, N, beta0 = beta, gamma0 = gamma, lambda0 = lambda0, V0 = V0, noSims = noSims, noIts = 10000,
                                          burnIn = 0)
    lambda = adaptRun$lambda
    V = adaptRun$V
  } else{
    adaptRun = NULL
    lambda = lambda0
    V = V0
  }

  print("==== MCMC Step ====")
  print(paste(c("==== MCMC Step (", noPanels, "Panels,", noSims, "Sims )",   "===="), collapse = " "))
  MCMCrun = SIS_PseudoMarginalMCMC(Y, I_0, obsTimes, N, beta0 = beta, gamma0 = gamma, lambda, V, noSims, noIts, burnIn = 0)

  return(list(MCMCrun = MCMCrun, adaptRun = adaptRun))
}
