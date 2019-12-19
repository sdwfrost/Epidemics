#' SIR Gillespie *NEW*

#' Split into seperate cases


# ==== Homogeneous SIR ====

#' State is (X, Y, Z)

#' Input
#' Initial State (X_0, Y_0, Z_0)

homogeneousSIR_Gillespie = function(state0, beta, gamma){

  #' Initialise state
  state = state0
  currentTime = 0
  #' Possible events are infection and removal
  infection = c(-1, 1, 0)
  removal   = c(0, -1, 1)

  #' Construct event matrix
  eventMatrix = matrix(c(infection, removal), nrow = 2, ncol = 3, byrow = T)

  #' Create event table
  eventTable = matrix(NA, nrow = 2*sum(state) - state[2] - 2*state[3] + 1, ncol = 5)

  eventTable[1, ] = c(currentTime, NA, state)

  eventNo = 1
  while(state[2] > 0){
    #' Calculate rate of each possible event
    infectionRate = state[1]*state[2]*beta
    removalRate   = state[2]*gamma

    #' Total rate will dictate the waiting time
    #' distribution.
    totalRate     = infectionRate + removalRate

    waitingTime = rexp(1, rate = totalRate)
    currentTime = currentTime + waitingTime

    #' Which event occured?
    u1 = runif(1, 0, totalEventRate)
    if(u1 < infectionRate){
      event = 1
    } else{
      event = 2
    }

    #'event = sample(c(1, 2), size = 1, prob = c(infectionRate, removalRate))

    #' Alter state of epidemic
    state = state + eventMatrix[event, ]

    #' Update Event table
    eventTable[eventNo + 1, ] = c(currentTime, event, state)
    eventNo = eventNo + 1
  }
  return(eventTable)
}

# ==== Individual level tracking in a homogenous SIR Gillespie ====

#' Incorporate the ideas of the event matrix
#' and a matrix of TRUE and FALSE indicating
#' which individual is in which state
#'
#' e.g for SIR
#'
#' 3 states, N individuals
#'
#' Nx3 matrix
#'
#' ith row pertains to the ith individual
#' One element of the ith row will be TRUE
#' which will indicate which infectious state
#' the ith individual is in. The other enteries
#' will have value FALSE.
#'
#' Due to TRUE/FALSE behaving in the exact same
#' way as 1/0 repsectively in R, if needed the
#' summary of number of individuals in each
#' state can be obtained by summing over the
#' respective columns.
#'
#' Hence, this matrix tells us everything about the
#' current state of the epidemic.

#' Need to know the set of individuals how are in each state
#' Also need summary of how many people (this will save operations)
#'
#' Not knowing how many people are in a state will mean everytime we want to
#' calculate the total rate we have to sum over the columns of the first two
#' columns of the stateMatrix
#'
#' This will mean there are 2N more operations per event.
#'
#' So it is more efficient to keep a similar state vector to
#' the homogeneous population level Gillespie. This takes
#' only one extra operation per event (updating the state vector)


individualHomogeneousSIR_Gillespie = function(initialState, beta, gamma){

  currentTime = 0
  #' Possible states and their respective
  #' logical representations.
  susceptible = c(TRUE, FALSE, FALSE)
  infected    = c(FALSE, TRUE, FALSE)
  removed     = c(FALSE, FALSE, TRUE)

  #' Summarise in matrix
  stateMatrix = matrix(c(susceptible, infected, removed), ncol = 3, byrow = T)

  #' Initial state of each individual
  N = length(initialState)
  individualStates = matrix(stateMatrix[initialState, ], ncol = 3, nrow = N)

  #' Initial state of epidemic at population level
  popState = .Internal(colSums(individualStates, dim(individualStates)[1], dim(individualStates)[2], FALSE))

  #' Possible Events
  infection = c(-1, 1, 0)
  removal = c(0, -1, 1)

  #' Summarise in a matrix
  eventMatrix = matrix(c(infection, removal), ncol = 3, byrow = TRUE)

  #' event table
  eventTable = matrix(NA, nrow = 2*N - popState[2] - 2*popState[3] + 1, ncol = 6)
  eventTable[1, ] = c(currentTime, NA, NA, popState)

  eventNo = 1
  while(popState[2] > 0){

    infectionRate = popState[1]*popState[2]*beta
    removalRate   = popState[2]*gamma
    totalEventRate = infectionRate + removalRate

    #' Waiting Time
    waitingTime = rexp(1, rate = totalEventRate)
    currentTime = currentTime + waitingTime

    #' #' Which event happens?
    u1 = runif(1, 0, totalEventRate)
    if(u1 < infectionRate){
      event = 1
    } else{
      event = 2
    }

    cs = cumsum(individualStates[,event])
    U  = runif(1, min = 0, max = cs[length(cs)])
    individual = sum(U > cs) + 1

    #' Update Population State
    popState = popState + eventMatrix[event, ]

    #' Update Individual State
    individualStates[individual, ] = stateMatrix[event + 1, ]

    #' Record event

    eventTable[eventNo + 1, ] = c(currentTime, event, individual, popState)
    eventNo = eventNo + 1
  }
  return(eventTable)
}



# ==== Heterogeneous (using Kernel) ====

heterogeneousSIR_Gillespie = function(initialState, beta, gamma, kernel){

  currentTime = 0
  B = kernel(beta)
  #' Possible states and their respective
  #' logical representations.
  susceptible = c(TRUE, FALSE, FALSE)
  infected    = c(FALSE, TRUE, FALSE)
  removed     = c(FALSE, FALSE, TRUE)

  #' Summarise in matrix
  stateMatrix = matrix(c(susceptible, infected, removed), ncol = 3, byrow = T)

  #' Initial state of each individual
  N = length(initialState)
  individualStates = matrix(stateMatrix[initialState, ], ncol = 3, nrow = N)

  #' Initial state of epidemic at population level
  popState = .Internal(colSums(individualStates, dim(individualStates)[1], dim(individualStates)[2], FALSE))
  #' Possible Events
  infection = c(-1, 1, 0)
  removal = c(0, -1, 1)


  #' Summarise in a matrix
  eventMatrix = matrix(c(infection, removal), ncol = 3, byrow = TRUE)

  #' event table
  eventTable = matrix(NA, nrow = 2*N - popState[2] - 2*popState[3] + 1, ncol = 6)
  eventTable[1, ] = c(currentTime, NA, NA, popState)

  eventNo = 1
  while(popState[2] > 0){
    reducedB = B[individualStates[,2], individualStates[,1], drop = F]
    dimReducedB = dim(reducedB)
    infectionRate = sum(reducedB)
    removalRate   = popState[2]*gamma
    totalEventRate = infectionRate + removalRate

    #' Waiting Time
    waitingTime = rexp(1, rate = totalEventRate)
    currentTime = currentTime + waitingTime

    #' Which event happens?
    #'event = sample(c(1,2), size = 1, prob = c(infectionRate, removalRate))

    #' To who?
    #' Construct rate vector
    who = c(which(individualStates[,1]), which(individualStates[,2]))
    rates = c(.Internal(colSums(reducedB, dimReducedB[1], dimReducedB[2], FALSE)), rep(gamma, popState[2]))
    cs = cumsum(rates)
    U  = runif(1, min = 0, max = cs[length(cs)])
    index = sum(U > cs) + 1

    event = if(index <= popState[1]){
      1
      } else {2}
    individual = who[sum(U > cs) + 1]

    # if(event == 2){
    #   individual = sample(1:N, size = 1, prob = individualStates[, 2], replace = T)
    # } else{
    #   reducedB = B[individualStates[,2], individualStates[,1], drop = F]
    #   dimReducedB = dim(reducedB)
    #   individual = sample(c(0, which(individualStates[,1])), size = 1,
    #                       prob = c(0, .Internal(colSums(reducedB, dimReducedB[1], dimReducedB[2], FALSE))),
    #                       replace = T)
    # }

    #' Update Population State
    popState = popState + eventMatrix[event, ]

    #' Update Individual State
    individualStates[individual, ] = stateMatrix[event + 1, ]

    #' Record event
    eventTable[eventNo + 1, ] = c(currentTime, event, individual, popState)
    eventNo = eventNo + 1
  }
  return(eventTable)
}


#' #heterogeneousSIR_Gillespie2 = function(initialState, beta, gamma, kernel){
#'
#' currentTime = 0
#' B = kernel(beta)
#' #' Possible states and their respective
#' #' logical representations.
#' susceptible = c(TRUE, FALSE, FALSE)
#' infected    = c(FALSE, TRUE, FALSE)
#' removed     = c(FALSE, FALSE, TRUE)
#'
#' #' Summarise in matrix
#' stateMatrix = matrix(c(susceptible, infected, removed), ncol = 3, byrow = T)
#'
#' #' Initial state of each individual
#' N = length(initialState)
#' individualStates = matrix(stateMatrix[initialState, ], ncol = 3, nrow = N)
#'
#' #' Initial state of epidemic at population level
#' popState = .Internal(colSums(individualStates, dim(individualStates)[1], dim(individualStates)[2], FALSE))
#' #' Possible Events
#' infection = c(-1, 1, 0)
#' removal = c(0, -1, 1)
#'
#'
#' #' Summarise in a matrix
#' eventMatrix = matrix(c(infection, removal), ncol = 3, byrow = TRUE)
#'
#' #' event table
#' eventTable = matrix(NA, nrow = 2*N - popState[2] - 2*popState[3] + 1, ncol = 6)
#' eventTable[1, ] = c(currentTime, NA, NA, popState)
#'
#' eventNo = 1
#' while(popState[2] > 0){
#'
#'   infectionRate = sum(B[individualStates[,2], individualStates[,1]])
#'   removalRate   = popState[2]*gamma
#'   totalEventRate = infectionRate + removalRate
#'
#'   #' Waiting Time
#'   waitingTime = rexp(1, rate = totalEventRate)
#'   currentTime = currentTime + waitingTime
#'
#'   #' Which event happens?
#'   event = sample(c(1,2), size = 1, prob = c(infectionRate, removalRate))
#'
#'   #' To who?
#'   # print(event)
#'   # print(popState)
#'   # print(sum(individualStates[,event]))
#'   if(event == 2){
#'     individual = sample(1:N, size = 1, prob = individualStates[, 2], replace = T)
#'   } else{
#'     reducedB = B[individualStates[,2], individualStates[,1], drop = F]
#'     dimReducedB = dim(reducedB)
#'     individual = sample(c(0, which(individualStates[,1])), size = 1,
#'                         prob = c(0, .Internal(colSums(reducedB, dimReducedB[1], dimReducedB[2], FALSE))),
#'                         replace = T)
#'   }
#'
#'   #' Update Population State
#'   popState = popState + eventMatrix[event, ]
#'
#'   #' Update Individual State
#'   individualStates[individual, ] = stateMatrix[event + 1, ]
#'
#'   #' Record event
#'   eventTable[eventNo + 1, ] = c(currentTime, event, individual, popState)
#'   eventNo = eventNo + 1
#' }
#' return(eventTable)
#' }

# ==== Output Panel Data ====

#' individual homogeneous

homogeneousPanelDataSIR_Gillespie = function(initialState, beta, gamma, obsTimes){

  currentTime = 0
  #' Possible states and their respective
  #' logical representations.
  susceptible = c(TRUE, FALSE, FALSE)
  infected    = c(FALSE, TRUE, FALSE)
  removed     = c(FALSE, FALSE, TRUE)
  rates = rep(NA, N)

  #' Summarise in matrix
  stateMatrix = matrix(c(susceptible, infected, removed), ncol = 3, byrow = T)
  stateLabels = c(1, 2, 3) #' 1 Susceptible, 2 Infected, 3 Removed


  #' Initial state of each individual
  N = length(initialState)
  individualStates = matrix(stateMatrix[initialState, ], ncol = 3, nrow = N)

  #' Initial state of epidemic at population level
  popState = .Internal(colSums(individualStates, dim(individualStates)[1], dim(individualStates)[2], FALSE))

  #' Possible Events
  infection = c(-1, 1, 0)
  removal = c(0, -1, 1)

  #' Summarise in a matrix
  eventMatrix = matrix(c(infection, removal), ncol = 3, byrow = TRUE)

  #' event table
  panelData = lapply(X = 1:length(obsTimes), function(X) NA)
  panelData[[1]] = matrix(c(1:N, initialState), nrow = N, ncol = 2, byrow = F)

  eventNo = 1
  while(popState[2] > 0){

    infectionRate = popState[1]*popState[2]*beta
    removalRate   = popState[2]*gamma
    totalEventRate = infectionRate + removalRate

    #' Waiting Time
    waitingTime = rexp(1, rate = totalEventRate)
    oldTime = currentTime
    currentTime = currentTime + waitingTime

    #' Panel Data
    #' Panel Data
    obsTimesPassed = currentTime > obsTimes & obsTimes > oldTime
    if(sum(obsTimesPassed) > 0){
      #' Convert individualStates into panelData
      currentState = sapply(X = 1:N, function(X) stateLabels[individualStates[X,]]) # THIS LINE (AVOID WHICH)
      currentPanelData = matrix(c(1:N, currentState), nrow = N, ncol = 2, byrow = F)
      panelData[obsTimesPassed] = list(currentPanelData)
    }

    #' #' Which event happens?
    u1 = runif(1, 0, totalEventRate)
    if(u1 < infectionRate){
      event = 1
    } else{
      event = 2
    }
    #event = sample(c(1,2), size = 1, prob = c(infectionRate, removalRate), replace = T)

    #' #' To who?
    #individual = sample(1:N, size = 1, prob = individualStates[,event], replace = T)
    cs = cumsum(individualStates[,event])
    u2 = runif(1, min = 0, max = cs[length(cs)])
    individual = sum(u2 > cs) + 1


    #' Update Population State
    popState = popState + eventMatrix[event, ]

    #' Update Individual State
    individualStates[individual, ] = stateMatrix[event + 1, ]

    eventNo = eventNo + 1
  }

  #' Filling in Observation times not passed
  obsTimesNotPassed = obsTimes > currentTime
  if(sum(obsTimesNotPassed) > 0){
    #' Convert individualStates into panelData
    currentState = sapply(X = 1:N, function(X) stateLabels[individualStates[X,]])
    currentPanelData = matrix(c(1:N, currentState), nrow = N, ncol = 2, byrow = F)
    panelData[obsTimesNotPassed] = list(currentPanelData)
  }
  return(list(panelData = panelData, obsTimes = obsTimes))
}



homogeneousPanelDataSIR_GillespieEU = function(initialState, beta, gamma, obsTimes, E, U){

  N = length(initialState)
  currentTime = 0
  lastObs = tail(obsTimes, n = 1)
  #' Possible states and their respective
  #' logical representations.
  #'
  susceptible = c(TRUE, FALSE, FALSE)
  infected    = c(FALSE, TRUE, FALSE)
  removed     = c(FALSE, FALSE, TRUE)
  rates = rep(NA, N)
  #' Summarise in matrix
  stateMatrix = matrix(c(susceptible, infected, removed), ncol = 3, byrow = T)
  stateLabels = c(1, 2, 3) #' 1 Susceptible, 2 Infected, 3 Removed

  #' Initial state of each individual
  N = length(initialState)
  individualStates = matrix(stateMatrix[initialState, ], ncol = 3, nrow = N)

  #' Initial state of epidemic at population level
  popState = .Internal(colSums(individualStates, dim(individualStates)[1], dim(individualStates)[2], FALSE))

  #' Possible Events
  infection = c(-1, 1, 0)
  removal = c(0, -1, 1)

  #' Summarise in a matrix
  eventMatrix = matrix(c(infection, removal), ncol = 3, byrow = TRUE)

  #' event table
  panelData = lapply(X = 1:length(obsTimes), function(X) NA)
  panelData[[1]] = matrix(c(1:N, initialState), nrow = N, ncol = 2, byrow = F)

  eventNo = 1
  while(popState[2] > 0 & currentTime < lastObs){

    if(is.na(U[eventNo])){
      return("Not Enough Draws")
    }

    infectionRate = popState[1]*popState[2]*beta
    removalRate   = popState[2]*gamma
    totalEventRate = infectionRate + removalRate

    #' Waiting Time
    waitingTime = E[eventNo]/totalEventRate
    oldTime = currentTime
    currentTime = currentTime + waitingTime

    #' Panel Data
    obsTimesPassed = currentTime > obsTimes & obsTimes > oldTime
    if(sum(obsTimesPassed) > 0){
      #' Convert individualStates into panelData
      currentState = sapply(X = 1:N, function(X) stateLabels[individualStates[X,]]) # THIS LINE (AVOID WHICH)
      currentPanelData = matrix(c(1:N, currentState), nrow = N, ncol = 2, byrow = F)
      panelData[obsTimesPassed] = list(currentPanelData)
    }

    rates[individualStates[,1]] = popState[2]*beta
    rates[individualStates[,2]] = gamma
    rates[individualStates[,3]] = 0

    cs = cumsum(rates)/totalEventRate
    individual = sum(U[eventNo] > cs) + 1
    if(individualStates[individual, 1]){
      event = 1
    } else{
      event = 2
    }

    #' Update Population State
    popState = popState + eventMatrix[event, ]

    #' Update Individual State
    individualStates[individual, ] = stateMatrix[event + 1, ]

    eventNo = eventNo + 1
  }

  #' Filling in Observation times not passed
  obsTimesNotPassed = obsTimes > currentTime
  if(sum(obsTimesNotPassed) > 0){
    #' Convert individualStates into panelData
    currentState = sapply(X = 1:N, function(X) stateLabels[individualStates[X,]])
    currentPanelData = matrix(c(1:N, currentState), nrow = N, ncol = 2, byrow = F)
    panelData[obsTimesNotPassed] = list(currentPanelData)
  }
  return(list(panelData = panelData, obsTimes = obsTimes, noEvents = eventNo))
}

panelDataSIR_Gillespie = function(initialState, beta, gamma, obsTimes, kernel, replace = F){
  currentTime = 0
  B = kernel(beta)
  #' Possible states and their respective
  #' logical representations.
  susceptible = c(TRUE, FALSE, FALSE)
  infected    = c(FALSE, TRUE, FALSE)
  removed     = c(FALSE, FALSE, TRUE)

  stateEnum = matrix(c(1,2,3), byrow = T)
  #' Summarise in matrix
  stateMatrix = matrix(c(susceptible, infected, removed), ncol = 3, byrow = T)

  #' Initial state of each individual
  N = length(initialState)
  individualStates = matrix(stateMatrix[initialState, ], ncol = 3, nrow = N)

  #' Initial state of epidemic at population level
  popState = .Internal(colSums(individualStates, N, 3, FALSE))
  N = sum(popState)
  #' Possible Events
  infection = c(-1, 1, 0)
  removal = c(0, -1, 1)

  #' Summarise in a matrix
  eventMatrix = matrix(c(infection, removal), ncol = 3, byrow = TRUE)

  #' Panel Data
  panelTemplate = matrix(c(1:N, rep(NA, N)), nrow = N, ncol = 2, byrow = F)
  panelData = lapply(X = 1:length(obsTimes), function(X) panelTemplate)
  panelData[[1]][, 2] = initialState
  eventNo = 1

  while(popState[2] > 0){

    reducedB = B[individualStates[,2], individualStates[,1], drop = F]
    dimReducedB = dim(reducedB)
    infectionRate = sum(reducedB)
    removalRate   = popState[2]*gamma
    totalEventRate = infectionRate + removalRate

    #' Waiting Time
    waitingTime = rexp(1, rate = totalEventRate)
    oldTime = currentTime
    currentTime = currentTime + waitingTime

    #' Panel Data
    obsTimesPassed = which(currentTime > obsTimes & obsTimes > oldTime)
    if(length(obsTimesPassed) > 0){
      #' Convert individualStates into panelData
      #' currentState = sapply(X = 1:N, function(X) which(individualStates[X,]))
      currentState = apply(individualStates, 1, function(X) stateEnum[X])

      currentState = matrix(c(1:N, currentState), nrow = N, ncol = 2, byrow = F)
      for(i in obsTimesPassed){
        panelData[[i]] = currentState
      }
    }

    #' Which event happens?
    event = sample(c(1,2), size = 1, prob = c(infectionRate, removalRate))

    #' To who?
    who = c(which(individualStates[,1]), which(individualStates[,2]))
    rates = c(.Internal(colSums(reducedB, dimReducedB[1], dimReducedB[2], FALSE)), rep(gamma, popState[2]))
    cs = cumsum(rates)
    U  = runif(1, min = 0, max = cs[length(cs)])
    index = sum(U > cs) + 1

    event = if(index <= popState[1]){
      1
    } else {2}
    individual = who[sum(U > cs) + 1]


    #' Update Population State
    popState = popState + eventMatrix[event, ]

    #' Update Individual State
    individualStates[individual, ] = stateMatrix[event + 1, ]


    eventNo = eventNo + 1
  }
  obsTimesNotPassed = which(obsTimes > currentTime)
  if(sum(obsTimesNotPassed) > 0){
    #' Convert individualStates into panelData
    currentState = sapply(X = 1:N, function(X) which(individualStates[X,]))
    currentState = matrix(c(1:N, currentState), nrow = N, ncol = 2, byrow = F)
    for(i in obsTimesNotPassed){
      panelData[[i]] = currentState
    }
  }
  return(list(panelData = panelData, obsTimes = obsTimes))
}

panelDataSample = function(panelData, m){

  samples = sample(1:length(panelData[[1]][,2]), size = m, replace = F)

  panelDataSample = lapply(X = panelData, function(X) X[samples, ])

  return(panelDataSample)
}


