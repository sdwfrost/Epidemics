#' Variance of beta, gamma for varying number of panels.
#'
#' Same SIS epidemic which occurs in population size N, (set.seed(1))
#'
#' We observe the same m = 0.1*N individuals, but vary how often they are observed.
#'
#' In this simulated situation, the average infectious period is a unit length of time (\gamma = 1).
#'
#' R_0 is chosen to be above 1 and relitively small, but large enough that the disease will reach an endemic state.
#'
#' \beta is then calculated from \gamma, R_0 and N
#'


# ==== Hypothesis ====

# There is some indication that having more panels leads to a decrease in posterior variance (for gamma and beta)

# This would make sense, as we have more information about the epidemic process.

# One problem that is anticipated, however, is that mixing will degrade as no. panels increases. This is due to it
# becoming harder to simulate an epidemic which will line up with the sample panels. Proposals which fail to do this
# will be thrown away as the likelihood estimate will evaluate to zero. Doing multiple simulations may help with
# this.

# ==== Set up ====

SISfsMCMCexperiment = function(noPanels, lastObs, noDraws, lambda0, V0, blockSize, adapt = TRUE, noIts = 10000, burnIn = 1000){

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

  # ==== Adapt ====
  if(adapt){
    print(paste(c("==== Adaptive Step (", noPanels, "Panels )",   "===="), collapse = " "))
    adaptiveBlockedIS_run = adaptiveSISfsMCMC7_blockedIS(Y, I_0, SIS_sim$obsTimes, N, beta0 = beta, gamma0 = gamma, thetaLim = c(Inf, Inf),
                                                         lambda0 = lambda0, noDraws = noDraws, blockSize = blockSize, noIts = 10000, delta = 0.05)

    lambda = adaptiveBlockedIS_run$lambda
    V = adaptiveBlockedIS_run$V
    noDraws = adaptiveBlockedIS_run$noDraws
    blockSize = adaptiveBlockedIS_run$blockSize
  } else{
    adaptiveBlockedIS_run = NULL
  }

  # ==== Final Run ====
  print(paste(c("==== MCMC Step (", noPanels, "Panels )",   "===="), collapse = " "))
  blockedIS_run = SIS_fsMCMC_blockedIS(Y, I_0, SIS_sim$obsTimes, N, beta0 = beta, gamma0 = gamma, thetaLim = c(Inf, Inf),
                                       lambda = lambda, V = V, noDraws = noDraws,
                                       blockSize = blockSize, noIts = noIts)

  return(list(adaptStep = adaptiveBlockedIS_run, MCMCstep = blockedIS_run))
}

