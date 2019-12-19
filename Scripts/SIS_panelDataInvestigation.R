#' Simulate SIS
#'
#' Sample Panel Data
#'

#
# R0 = 4
# I_0 = 1
#' 1. Simulate SIS


#' One Particle Psuedo-Marginal MCMC
lambda = 0.01
noIts = 1000
burnIn = 0
run = oneParticlePM_SIS(Y, SIS_sim$obsTimes, N, beta0 = beta, gamma0 = gamma, lambda, noIts, burnIn)


#' Forward Simulation MCMC
noIts = 1000
noDraws = 2000
lambda = 0.001
V = diag(c(1/N, 1))
s = 500
thetaLim = 2*c(beta, gamma)
run2 = SIS_fsMCMC(Y, I_0, SIS_sim$obsTimes, N, beta0 = beta, gamma0 = gamma, thetaLim, lambda, V, noDraws, s, noIts)


#' Epidemic Gillespie simulation template

#' Individual level gillepsie
#' A transition is assumed to move an individual from one state to another

SIR_gillespieSetup = function(N){

  #' Event Table Template
  eventTableTemplate <<- matrix(NA, nrow = 2*N, ncol = 5)

  #' Stochiometry Matrix, the effect of an event occuring
  stochMatrix <<- matrix(c(-1, 1, 0,
                         0, -1, 1), ncol = 3, byrow = T)

  #' Rate equations
  rateEquations <<- function(beta, gamma, popState){
    infectionRate = popState[1]*popState[2]*beta
    removalRate = popState[2]*gamma
    return(c(infectionRate, removalRate))
  }

  whatEvent <<- function(U, infectionRate, totalRate){
    if(missing(U)){
      U = runif(1, 0, totalRate)
    }
    if(U < infectionRate){
      return(1)
    } else{
      return(2)
    }
  }
  #'return(list(stochMatrix, rateEquations, whatEvent))
}



SIR_gillespie = function(beta, gamma, initialState){
  currentTime = 0
  popState = initialState
  SIR_gillespieSetup(sum(initialState))

  eventTable = eventTableTemplate
  eventTable[1, ] = c(currentTime, NA, popState)
  noEvent = 1
  while(popState[2] > 0){
    rates = rateEquations(beta, gamma, popState)
    Totalrate = sum(rates)
    #' Time til next event
    oldTime = currentTime
    currentTime = oldTime + rexp(1, rate = Totalrate)

    #' Which event
    u1 = runif(1, 0, Totalrate)
    if(u1 < rates[1]){
      event = 1
      popState = popState + stochMatrix[event, ]
    } else{
      event = 2
      popState = popState + stochMatrix[event, ]
      print(popState)
    }
    eventTable[noEvent + 1, ] = c(currentTime, event, popState)
  }
  return(eventTable)
}

# rateEquations = function(reducedB, gamma, individualStates, popState){
#   dimReducedB = dim(reducedB)
#   infectionRate = .Internal(colSums(reducedB, dimReducedB[1], dimReducedB[2], na.rm = FALSE))
#   removalRate = popState[2]*theta[2]
# }

#' Gillespie naming convention

#' The state space

#' NoStates (Redundent be)
noStates = 3
#' Directed Graph to show possible transitions

#' S --> I --> R
SIRgraph = matrix(c(0, 1, 0,
                    0, 0, 1,
                    0, 0, 0), byrow = T, nrow = 3, ncol = 3)

#' Connection List
SIRgraph = list(2, 3, NULL)

#' Define the nature of transitions between states.


#' Population Struture (Homogeneous/Heterogeneous)


# ==== Distribution of number of events ====

noEventsSIS = sapply(1:10000,
       function(X) nrow(individualHomogeneousSIS_Gillespie(initialState, beta, gamma, obsTimes)$eventTable) - 1)

par(mfrow = c(1,1))
hist(noEventsSIS)





