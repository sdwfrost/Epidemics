#' SIR_Gillespie function
#'
#' Takes parameters for an SIR Infectious Process and returns either an event table or
#' panel data.

SIR_Gillespie = function(N, a, beta, gamma, trans = possible_transitions(1:3),
                         obs_end = Inf, kernel = NULL,
                         U, E, output = "event", obs_times = NULL){

  current_time = 0
  X = N - a
  Y = a
  Z = N - X - Y

  #' Track specific people
  S = 1:(N - a) #' Set of Susceptibles
  I = (N - a + 1):N #' Set of Infectives
  R = NULL #' Set of removed

  current_state = state(S, I, R)
  #' Infection Matrix
  if(!is.null(kernel)){
    B = kernel(beta)
  }

  if(output == "panel"){
    event_table = NULL
    panel_data = lapply(rep(NA, length(obs_times)), function(X) return(X))
  } else{
    panel_data = NULL
    event_table = data.frame(matrix(NA, nrow = 2*N + 1, ncol = 4))

    event_table[1:(N-a), ] = matrix(c(1:(N-a), rep(0, N - a), rep(1, N - a), rep(1, N - a)),
                                    nrow = N - a, ncol = 4, byrow = FALSE)
    event_table[(N - a + 1):N, ] = matrix(c((N - a + 1):N, rep(0, a), rep(2, a), rep(2, a)),
                                          nrow = a, ncol = 4, byrow = FALSE)
    colnames(event_table) = c("ID", "times", "state", "prev_state")
  }

  no_event = 1
  while(Y[no_event] > 0 & current_time < obs_end){

    old_time = current_time

    if(!is.null(kernel)){
      reduced_B = B[I,S, drop = F]
    }

    if(X[no_event] == 0){
      individual_inf_rate = numeric(0)
    } else{
      if(!is.null(kernel)){
        individual_inf_rate = Matrix::colSums(reduced_B)
      } else{
        individual_inf_rate = rep(Y[no_event]*beta, X[no_event])
      }
    }

    total_inf_rate = sum(individual_inf_rate)
    total_rem_rate = Y[no_event]*gamma
    rate_next_event = total_rem_rate + total_inf_rate
    if(missing(E)){
      time_to_next_event = rexp(1, rate_next_event)
    } else{
      time_to_next_event = E[no_event]/(rate_next_event)
    }
    current_time = current_time + time_to_next_event

    if(output == "panel"){
      obs_times_passed = which(old_time <= obs_times & obs_times <= current_time)

      if(length(obs_times_passed) > 0){

        prev_state = current_state
        current_state = state(S, I, R)

        transitions = state_table(prev_state, current_state, trans)

        panel_data[[obs_times_passed[1]]] = transitions
        if(length(obs_times_passed) > 1){
          no_events_transitions = state_table(current_state, current_state, trans)
        }
        for(i in obs_times_passed[-1]){
          panel_data[[i]] = no_events_transitions
        }
      }
    }

    event = event.epidemics(individual_inf_rate, gamma, Y[no_event], if(!missing(U)){
                                                             U[no_event]} else{NULL})
    if(event$event == 1){
      individual = I[event$ID_index]
      I = I[-event$ID_index]
      R = c(R, individual)

      X[no_event + 1] = X[no_event]
      Y[no_event + 1] = Y[no_event] - 1
      Z[no_event + 1] = Z[no_event] + 1

      prev_state = 2
      state = 3
    } else{
      individual = S[event$ID_index]
      S = S[-event$ID_index]
      I = c(I, individual)

      X[no_event + 1] = X[no_event] - 1
      Y[no_event + 1] = Y[no_event] + 1
      Z[no_event + 1] = Z[no_event]

      prev_state = 1
      state = 2
    }
    if(is.na(individual)){
      print(no_event)
    }
    if(output != "panel"){
      event_table[N + no_event, ] = c(individual, current_time, state, prev_state)
    }
    no_event = no_event + 1
  }
  #' Last and unused panels
  if(output == "panel" & sum(is.na(panel_data)) > 0){
    NA_panels = which(is.na(panel_data))
    if(current_time < obs_times[NA_panels[1]]){
      old_state = current_state
      current_state = state(S, I, R)
      panel_data[[NA_panels[1]]] = state_table(old_state, current_state, trans)
      no_events_transitions = state_table(current_state, current_state, trans)
      for(i in NA_panels[-1]){
        panel_data[[i]] = no_events_transitions
      }
    } else{
      no_events_transitions = state_table(current_state, current_state, trans)
      for(i in NA_panels){
        panel_data[[i]] = no_events_transitions
      }
    }

  }

  event_table = na.omit(event_table)
  return(list(event_table = event_table, X = X, Y = Y, Z = Z, kernel = kernel,
              panel_data = panel_data))
}

