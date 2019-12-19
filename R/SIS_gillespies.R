#' ==== SIS Gillespie with Panel Data Output ====

homogeneousPanelDataSIS_Gillespie = function(initialState, beta, gamma, obsTimes){

  N = length(initialState)
  currentTime = 0
  #' Possible states and their respective
  #' logical representations.
  susceptible = c(TRUE, FALSE)
  infected    = c(FALSE, TRUE)

  #' Summarise in matrix
  stateMatrix = matrix(c(susceptible, infected), ncol = 2, byrow = T)
  stateLabels = c(1, 2) #' 1 Susceptible, 2 Infected

  #' Initial state of each individual
  N = length(initialState)
  individualStates = matrix(stateMatrix[initialState, ], ncol = 2, nrow = N)

  #' Initial state of epidemic at population level
  popState = .Internal(colSums(individualStates, dim(individualStates)[1], dim(individualStates)[2], FALSE))

  #' Possible Events
  infection = c(-1, 1)
  recovery = c(1, -1)

  #' Summarise in a matrix
  eventMatrix = matrix(c(infection, recovery), ncol = 2, byrow = TRUE)

  #' event table
  #panelTemplate = matrix(nrow = N, ncol = 2)
  panelData = lapply(X = 1:length(obsTimes), function(X) NA)
  panelData[[1]] = matrix(c(1:N, initialState), nrow = N, ncol = 2, byrow = F)
  eventNo = 1
  while(popState[2] > 0 & currentTime < tail(obsTimes, n = 1)){
    #print(currentTime)
    infectionRate = popState[1]*popState[2]*beta
    recoveryRate   = popState[2]*gamma
    totalEventRate = infectionRate + recoveryRate

    #' Waiting Time
    waitingTime = rexp(1, rate = totalEventRate)
    oldTime = currentTime
    currentTime = currentTime + waitingTime

    #' Panel Data
    obsTimesPassed = currentTime > obsTimes & obsTimes > oldTime
    if(sum(obsTimesPassed) > 0){
      #' Convert individualStates into panelData
      #' This which does take up significant time
      #currentState = sapply(X = 1:N, function(X) which(individualStates[X,])) # THIS LINE (AVOID WHICH)
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
    individualStates[individual, ] = stateMatrix[abs(event - 3), ]

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
  return(list(panelData = panelData, obsTimes = obsTimes, noDraws = eventNo))
}

#' ==== Non-Centered SIS Gillespie with Panel Data Output ====

homogeneousPanelDataSIS_GillespieEU = function(initialState, beta, gamma, obsTimes, E, U){

  N = length(initialState)
  noRVs = length(E)
  currentTime = 0
  #' Possible states and their respective
  #' logical representations.
  susceptible = c(TRUE, FALSE)
  infected    = c(FALSE, TRUE)
  rates = rep(NA, N)
  #' Summarise in matrix
  stateMatrix = matrix(c(susceptible, infected), ncol = 2, byrow = T)
  stateLabels = c(1, 2) #' 1 Susceptible, 2 Infected

  #' Initial state of each individual
  N = length(initialState)
  individualStates = matrix(stateMatrix[initialState, ], ncol = 2, nrow = N)

  #' Initial state of epidemic at population level
  popState = .Internal(colSums(individualStates, dim(individualStates)[1], dim(individualStates)[2], FALSE))

  #' Possible Events
  infection = c(-1, 1)
  recovery = c(1, -1)

  #' Summarise in a matrix
  eventMatrix = matrix(c(infection, recovery), ncol = 2, byrow = TRUE)

  #' event table
  #panelTemplate = matrix(nrow = N, ncol = 2)
  panelData = lapply(X = 1:length(obsTimes), function(X) NA)
  panelData[[1]] = matrix(c(1:N, initialState), nrow = N, ncol = 2, byrow = F)
  eventNo = 1
  while(popState[2] > 0 & currentTime < tail(obsTimes, n = 1)){

    infectionRate = popState[1]*popState[2]*beta
    recoveryRate   = popState[2]*gamma
    totalEventRate = infectionRate + recoveryRate

    if(is.na(U[eventNo])){
      return("Not Enough Draws")
    }

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


    #'To who?
    #individual = sample(1:N, size = 1, prob = individualStates[,event], replace = T)
    rates[individualStates[,1]] = popState[2]*beta
    rates[individualStates[,2]] = gamma

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
    individualStates[individual, ] = stateMatrix[abs(event - 3), ]

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


individualHomogeneousSIS_Gillespie = function(initialState, beta, gamma, obsTimes){

  currentTime = 0
  #' Possible states and their respective
  #' logical representations.
  susceptible = c(TRUE, FALSE)
  infected    = c(FALSE, TRUE)

  #' Summarise in matrix
  stateMatrix = matrix(c(susceptible, infected), ncol = 2, byrow = T)

  #' Initial state of each individual
  N = length(initialState)
  individualStates = matrix(stateMatrix[initialState, ], ncol = 2, nrow = N)

  #' Initial state of epidemic at population level
  popState = .Internal(colSums(individualStates, dim(individualStates)[1], dim(individualStates)[2], FALSE))

  #' Possible Events
  infection = c(-1, 1)
  recovery = c(1, -1)

  #' Summarise in a matrix
  eventMatrix = matrix(c(infection, recovery), ncol = 2, byrow = TRUE)

  #' event table
  eventTable = matrix(NA, nrow = 100*N, ncol = 5)
  eventTable[1, ] = c(currentTime, NA, NA, popState)

  eventNo = 1
  while(popState[2] > 0 & currentTime < tail(obsTimes, n = 1)){

    infectionRate = popState[1]*popState[2]*beta
    recoveryRate   = popState[2]*gamma
    totalEventRate = infectionRate + recoveryRate

    #' Waiting Time
    waitingTime = rexp(1, rate = totalEventRate)
    oldTime = currentTime
    currentTime = currentTime + waitingTime

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
    individualStates[individual, ] = stateMatrix[abs(event - 3), ]

    eventTable[eventNo + 1, ] = c(currentTime, event, individual, popState)
    eventNo = eventNo + 1
  }

  eventTable = rbind(eventTable[1, ], na.omit(eventTable))
  return(list(eventTable = eventTable, obsTimes = obsTimes))
}

#' homogeneousPanelDataSIS_Gillespie = function(initialState, beta, gamma, obsTimes, U, E){
#'
#'   currentTime = 0
#'   #' Possible states and their respective
#'   #' logical representations.
#'   susceptible = c(TRUE, FALSE)
#'   infected    = c(FALSE, TRUE)
#'
#'   #' Summarise in matrix
#'   stateMatrix = matrix(c(susceptible, infected), ncol = 2, byrow = T)
#'
#'   #' Initial state of each individual
#'   N = length(initialState)
#'   individualStates = matrix(stateMatrix[initialState, ], ncol = 2, nrow = N)
#'
#'   #' Initial state of epidemic at population level
#'   popState = .Internal(colSums(individualStates, dim(individualStates)[1], dim(individualStates)[2], FALSE))
#'
#'   #' Possible Events
#'   infection = c(-1, 1)
#'   recovery = c(1, -1)
#'
#'   #' Summarise in a matrix
#'   eventMatrix = matrix(c(infection, recovery), ncol = 2, byrow = TRUE)
#'
#'   #' event table
#'   panelTemplate = matrix(c(1:N, rep(NA, N)), nrow = N, ncol = 2, byrow = F)
#'   panelData = lapply(X = 1:length(obsTimes), function(X) panelTemplate)
#'   panelData[[1]] = matrix(c(1:N, initialState), nrow = N, ncol = 2, byrow = F)
#'   eventNo = 1
#'   while(popState[2] > 0 & currentTime < tail(obsTimes, n = 1)){
#'
#'     infectionRate = popState[1]*popState[2]*beta
#'     recoveryRate   = popState[2]*gamma
#'     totalEventRate = infectionRate + recoveryRate
#'
#'     #' Waiting Time
#'     waitingTime = E[no_event]/totalEventRate
#'     oldTime = currentTime
#'     currentTime = currentTime + waitingTime
#'
#'     #' Panel Data
#'     obsTimesPassed = which(currentTime > obsTimes & obsTimes > oldTime)
#'     if(length(obsTimesPassed) > 0){
#'       #' Convert individualStates into panelData
#'       currentState = sapply(X = 1:N, function(X) which(individualStates[X,]))
#'       currentState = matrix(c(1:N, currentState), nrow = N, ncol = 2, byrow = F)
#'       for(i in obsTimesPassed){
#'         panelData[[i]] = currentState
#'       }
#'     }
#'
#'     #' #' Which event happens?
#'     #'u1 = runif(1, 0, totalEventRate)
#'     if(U[event_no] < infectionRate){
#'       event = 1
#'     } else{
#'       event = 2
#'     }
#'     #event = sample(c(1,2), size = 1, prob = c(infectionRate, removalRate), replace = T)
#'
#'     #' #' To who?
#'     #individual = sample(1:N, size = 1, prob = individualStates[,event], replace = T)
#'     cs = cumsum(individualStates[,event])
#'     u2 = runif(1, min = 0, max = cs[length(cs)])
#'     individual = sum(u2 > cs) + 1
#'
#'     #' Update Population State
#'     popState = popState + eventMatrix[event, ]
#'
#'     #' Update Individual State
#'     individualStates[individual, ] = stateMatrix[abs(event - 3), ]
#'
#'     eventNo = eventNo + 1
#'   }
#'
#'   #' Filling in Observation times not passed
#'   obsTimesNotPassed = which(obsTimes > currentTime)
#'   if(length(obsTimesNotPassed) > 0){
#'     #' Convert individualStates into panelData
#'     currentState = sapply(X = 1:N, function(X) which(individualStates[X,]))
#'     currentState = matrix(c(1:N, currentState), nrow = N, ncol = 2, byrow = F)
#'     for(i in obsTimesNotPassed){
#'       panelData[[i]] = currentState
#'     }
#'   }
#'   return(list(panelData = panelData, obsTimes = obsTimes))
#' }
#'
#'
#' individualHomogeneousSIS_Gillespie = function(initialState, beta, gamma, obsTimes){
#'
#'   currentTime = 0
#'   #' Possible states and their respective
#'   #' logical representations.
#'   susceptible = c(TRUE, FALSE)
#'   infected    = c(FALSE, TRUE)
#'
#'   #' Summarise in matrix
#'   stateMatrix = matrix(c(susceptible, infected), ncol = 2, byrow = T)
#'
#'   #' Initial state of each individual
#'   N = length(initialState)
#'   individualStates = matrix(stateMatrix[initialState, ], ncol = 2, nrow = N)
#'
#'   #' Initial state of epidemic at population level
#'   popState = .Internal(colSums(individualStates, dim(individualStates)[1], dim(individualStates)[2], FALSE))
#'
#'   #' Possible Events
#'   infection = c(-1, 1)
#'   recovery = c(1, -1)
#'
#'   #' Summarise in a matrix
#'   eventMatrix = matrix(c(infection, recovery), ncol = 2, byrow = TRUE)
#'
#'   #' event table
#'   eventTable = matrix(NA, nrow = 100*N, ncol = 5)
#'   eventTable[1, ] = c(currentTime, NA, NA, popState)
#'
#'   eventNo = 1
#'   while(popState[2] > 0 & currentTime < tail(obsTimes, n = 1)){
#'
#'     infectionRate = popState[1]*popState[2]*beta
#'     recoveryRate   = popState[2]*gamma
#'     totalEventRate = infectionRate + recoveryRate
#'
#'     #' Waiting Time
#'     waitingTime = rexp(1, rate = totalEventRate)
#'     oldTime = currentTime
#'     currentTime = currentTime + waitingTime
#'
#'     #' #' Which event happens?
#'     u1 = runif(1, 0, totalEventRate)
#'     if(u1 < infectionRate){
#'       event = 1
#'     } else{
#'       event = 2
#'     }
#'     #event = sample(c(1,2), size = 1, prob = c(infectionRate, removalRate), replace = T)
#'
#'     #' #' To who?
#'     #individual = sample(1:N, size = 1, prob = individualStates[,event], replace = T)
#'     cs = cumsum(individualStates[,event])
#'     u2 = runif(1, min = 0, max = cs[length(cs)])
#'     individual = sum(u2 > cs) + 1
#'
#'     #' Update Population State
#'     popState = popState + eventMatrix[event, ]
#'
#'     #' Update Individual State
#'     individualStates[individual, ] = stateMatrix[abs(event - 3), ]
#'
#'     eventTable[eventNo + 1, ] = c(currentTime, event, individual, popState)
#'     eventNo = eventNo + 1
#'   }
#'
#'   eventTable = rbind(eventTable[1, ], na.omit(eventTable))
#'   return(list(eventTable = eventTable, obsTimes = obsTimes))
#' }
#'
#'
#'
#'
#'
#'
