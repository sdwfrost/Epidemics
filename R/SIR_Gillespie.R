#'
#' Epidemic Gillespie w/ kernel set-up
#'

state = function(S, I, R){
  state = c(rep(1, length(S)), rep(2, length(I)), rep(3, length(R)))
  state_data = data.frame(ID = c(S, I, R), state = state)
  state_data = dplyr::arrange(state_data, ID)
  return(state_data$state)
}

possible_transitions = function(states){
  possible_trans = expand.grid(states, states)
  possible_trans = lapply(X = 1:nrow(possible_trans), function(X) as.numeric(possible_trans[X,]))
  return(possible_trans)
}

state_table = function(x, y, trans){
  #if(is.null(states)){
  #  states = levels
  #}
  #if(missing(levels)){
  #  min(x,y):max(x,y)
  #}
  #if(length(states) != length(levels) & !is.null(states)){
  #  stop("Number of states and levels must be equal")
  #}
  #return(table(factor(x, levels, states), factor(y, levels, states)))
  freq = sapply(X = trans, function(X) sum(X[1] == x & X[2] == y))
}

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
    colnames(event_table) = c("ID", "time", "state", "prev_state")
  }

  no_event = 1
  while(Y[no_event] > 0 & current_time < obs_end){

    old_time = current_time

    if(!is.null(kernel)){
      reduced_B = B[I,S, drop = F]
    }

    if(X[no_event] == 0){
      individual_inf_rate = 0
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

Kernel_Deterministic_Gillespie = function(N, a, beta, gamma, E, U, T_obs, k, kernel, store = TRUE){

  #' Initialise
  current_time = 0
  X = N - a
  Y = a
  Z = N - X - Y

  S = 1:(N - a)
  I_i = (N - a + 1):N
  I_s = NULL

  B = kernel(beta)

  #' Panel Data
  if(length(T_obs) == 1){
    T_obs = c(0, T_obs)
  }

  obs_times = seq(T_obs[1], T_obs[2], length = k)

  #' Tracking Variables for Panel Data
  X_t0 = X
  Y_t0 = Y
  Y_s = 0
  Y_i = Y

  #' Data Storage
  if(store){
    sim_data = matrix(NA, nrow = 2*N + 1, ncol = 5)
    sim_data[1, ] = c(current_time, X, Y, Z, NA)
  } else{
    sim_data = NULL
  }


  event_no = 1

  while(Y > 0){
    old_time = current_time

    # Rate of Infection

    reduced_B = B[c(I_s, I_i),S, drop = F]

    if(X == 0){
      individual_inf_rate = 0
    } else{
      individual_inf_rate = Matrix::colSums(reduced_B)
    }

    total_inf_pressure = sum(individual_inf_rate)

    # Rate of Removal
    rem_rate_from_S = Y_s*gamma
    rem_rate_from_I = Y_i*gamma
    total_rem_rate =  rem_rate_from_S + rem_rate_from_I

    # Time of Next Event
    rate_next_event = total_rem_rate + total_inf_pressure
    time_to_next_event = (1/rate_next_event)*E[event_no]
    current_time = current_time + time_to_next_event

    #' Recording Panel Data
    #' Record which observation times have been passed
    obs_times_passed = which(old_time <= obs_times & obs_times <= current_time)

    if(length(obs_times_passed) > 0){
      # Record Panel
      panel_data[[obs_times_passed[1]]] = c(n_ss = X, n_si = Y_s, n_sr = X_t0 - (X + Y_s), n_ii = Y_i, n_ir = Y_t0 - Y_i, n_rr = N - X_t0 - Y_t0)

      # Reset
      X_t0 = X
      Y_t0 = Y
      Y_i = Y
      Y_s = 0

      I_i = c(I_i, I_s)
      I_s = NULL

      for(i in obs_times_passed[-1]){
        panel_data[[i]] = c(n_ss = X, n_si = Y_s, n_sr = X_t0 - (X + Y_s), n_ii = Y_i, n_ir = Y_t0 - Y_i, n_rr = N - X_t0 - Y_t0)
      }
    }

    which_event = Heterogeneous_Event_Deterministic(individual_inf_rate, rep(gamma, Y), U[event_no])

    if(which_event <= X){ #' Infection
      individual = S[which_event]
      #' Susceptibles
      X = X - 1
      S = S[-which_event]

      #' Infected
      Y_s = Y_s + 1
      Y = Y + 1
      I_s = c(I_s, individual)

    } else if(which_event <= X + Y_s){

      individual = I_s[which_event - X]

      Y_s = Y_s - 1
      Y = Y - 1
      I_s = I_s[!(I_s == individual)]

    } else{
      individual = I_i[which_event - X - Y_s]
      Y_i = Y_i - 1
      Y = Y - 1
      I_i = I_i[!(I_i == individual)]
    }

    #' Store old time
    #' Update extra parameters
    Z = N - X - Y

    #sim_data = rbind(sim_data, c(current_time, X, Y, Z))
    if(store){
      sim_data[event_no + 1,] = c(current_time, X, Y, Z, which_event)
    }
    event_no = event_no + 1
  }

  NA_panels = which(is.na(panel_data))

  for(i in NA_panels){
    panel_data[[i]] = c(n_ss = X, n_si = Y_s, n_sr = X_t0 - (X + Y_s), n_ii = Y_i, n_ir = Y_t0 - Y_i, n_rr = N - X_t0 - Y_t0)
  }

  return(list(sim_data = sim_data, panel_data = panel_data, no_events = event_no))
}


Heterogeneous_Event_Deterministic = function(individual_inf_rate, removal_rates, U){
  if(length(individual_inf_rate) == 1){
    if(individual_inf_rate == 0){
      probs = c(removal_rates)/sum(c(removal_rates))
    } else{
      probs = c(individual_inf_rate, removal_rates)/sum(c(individual_inf_rate, removal_rates))
    }
  } else{
    probs = c(individual_inf_rate, removal_rates)/sum(c(individual_inf_rate, removal_rates))
  }

  probs_cumsum = cumsum(probs)

  index =  sum(probs_cumsum < U) + 1

  return(index)
}
