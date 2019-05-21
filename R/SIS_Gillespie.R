#'
#' SIS Epidemic Simulation
#'
#'

SIS_Gillespie = function(N, a, gamma, beta, obs_end = Inf, kernel, U, E){

  current_time = 0
  X = N - a
  Y = a

  #' Track specific people
  S = 1:(N - a) #' Set of Susceptibles
  I = (N - a + 1):N #' Set of Infectives

  #' Infection Matrix
  if(!missing(kernel)){
    B = kernel(beta)
  }

  event_table = data.frame(matrix(NA, nrow = 2*N + 1, ncol = 4))
  colnames(event_table) = c("ID", "time", "state", "prev_state")

  event_table[1,] = c(NA, current_time, NA, NA)

  event_table[1:(N-a), ] = matrix(c(1:(N-a), rep(0, N - a), rep(1, N - a), rep(1, N - a)),
                                  nrow = N - a, ncol = 4, byrow = FALSE)
  event_table[(N - a + 1):N, ] = matrix(c((N - a + 1):N, rep(0, a), rep(2, a), rep(2, a)),
                                         nrow = a, ncol = 4, byrow = FALSE)

  no_event = 1
  while(Y[no_event] > 0 & current_time < obs_end){

    old_time = current_time

    if(!missing(kernel)){
      reduced_B = B[I,S, drop = F]
    }

    if(X[no_event] == 0){
      individual_inf_rate = 0
    } else{
      if(!missing(kernel)){
        individual_inf_rate = Matrix::colSums(reduced_B)
      } else{
        individual_inf_rate = rep(Y[no_event]*beta, X[no_event])
      }
    }

    total_inf_pressure = sum(individual_inf_rate)
    total_rem_rate = Y[no_event]*gamma
    rate_next_event = total_rem_rate + total_inf_pressure

    if(missing(E)){
      time_to_next_event = rexp(1, rate_next_event)
    } else{
      time_to_next_event = E[no_event]/rate_next_event
    }
    current_time = current_time + time_to_next_event

    event = event.epidemics(individual_inf_rate, gamma, Y[no_event], U = if(!missing(U)){
      U[no_event]} else{NULL})

    if(event$event == 1){

      individual = I[event$ID_index]

      I = I[-event$ID_index]
      S = c(S, individual)

      Y[no_event + 1] = Y[no_event] - 1
      X[no_event + 1] = X[no_event] + 1

      prev_state = 2
      state = 1
    } else{
      individual = S[event$ID_index]
      S = S[-event$ID_index]
      I = c(I, individual)

      X[no_event + 1] = X[no_event] - 1
      Y[no_event + 1] = Y[no_event] + 1

      prev_state = 1
      state = 2
    }
    event_table[N + no_event, ] = c(individual, current_time, state, prev_state)
    no_event = no_event + 1
    #print(current_time)
  }
  event_table = na.omit(event_table)
  return(list(event_table = event_table, X = X, Y = Y, kernel = if(!missing(kernel)){
                                                                  kernel} else{NULL}))
}

# ==== Deterministic ====


SIS_Deterministic_Gillespie = function(N, a, beta, gamma, E, U, T_obs, k, kernel, obs_end, store = TRUE){

  #' Initialise
  current_time = 0
  X = N - a
  Y = a

  S_s = 1:(N - a)
  S_i = NULL

  I_i = (N - a + 1):N
  I_s = NULL

  B = kernel(beta)

  #' Panel Data
  if(length(T_obs) == 1){
    T_obs = c(0, T_obs)
  }

  obs_times = seq(T_obs[1], T_obs[2], length = k)

  #' Tracking Variables for Panel Data

  Y_s = 0
  Y_i = Y
  X_s = X
  X_i = 0


  #' Data Storage
  if(store){
    sim_data = matrix(NA, nrow = length(U), ncol = 4)
    sim_data[1, ] = c(current_time, X, Y, NA)
  } else{
    sim_data = NULL
  }

  panel_data = lapply(rep(NA, k), function(X) return(X))

  event_no = 1

  while(Y > 0 & current_time < obs_end){
    old_time = current_time

    # Rate of Infection

    reduced_B = B[c(I_s, I_i),c(S_s, S_i), drop = F]


    if(X == 0){
      individual_inf_rate = 0
    } else{
      # First X_s rates will be those corresponding to susceptibles
      # which started in the susceptible state.
      # The next X_i rates will correspond to those who started in the
      # infected state
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
      panel_data[[obs_times_passed[1]]] = c(n_ss = X_s, n_si = Y_s, n_ii = Y_i, n_is = X_i)

      # Reset
      X_t0 = X
      Y_t0 = Y

      Y_i = Y
      Y_s = 0
      X_s = X
      X_i = 0

      S_s = c(S_s, S_i)
      S_i = NULL

      I_i = c(I_i, I_s)
      I_s = NULL

      for(i in obs_times_passed[-1]){
        panel_data[[i]] = c(n_ss = X_s, n_si = Y_s, n_ii = Y_i, n_is = X_i)
      }
    }

    which_event = Heterogeneous_Event_Deterministic(individual_inf_rate, rep(gamma, Y), U[event_no])

    if(which_event <= X_s){ #' Infection
      individual = S_s[which_event]
      #' Susceptibles
      X = X - 1
      X_s = X_s - 1
      S_s = S_s[!(individual == S_s)]

      #' Infected
      Y = Y + 1
      Y_s = Y_s + 1
      I_s = c(I_s, individual)

    } else if(which_event <= X_s + X_i){
      individual = S_i[which_event - X_s]
      #' Susceptibles
      X = X - 1
      X_i = X_i - 1
      S_i = S_i[!(individual == S_i)]

      #' Infected
      Y = Y + 1
      Y_i = Y_i + 1
      I_i = c(I_i, individual)


    } else if(which_event <= X + Y_s){

      individual = I_s[which_event - X]

      #' Infected
      Y = Y - 1
      Y_s = Y_s - 1
      I_s = I_s[!(I_s == individual)]

      #' Susceptibles
      X = X + 1
      X_s = X_s + 1
      S_s = c(S_s, individual)

    } else{
      individual = I_i[which_event - X - Y_s]

      #' Infected
      Y_i = Y_i - 1
      Y = Y - 1
      I_i = I_i[!(I_i == individual)]

      #' Susceptibles
      X = X + 1
      X_i = X_i + 1
      S_i = c(S_i, individual)

    }

    if(store){
      sim_data[event_no + 1,] = c(current_time, X, Y, which_event)
    }
    event_no = event_no + 1
  }

  NA_panels = which(is.na(panel_data))

  for(i in NA_panels){
    panel_data[[i]] = c(n_ss = X_s, n_si = Y_s, n_ii = Y_i, n_is = X_i)
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
  #print(c(probs_cumsum[100], U))
  index = sum(probs_cumsum < U) + 1
  #print(c("Index", index))
  return(index)
}

