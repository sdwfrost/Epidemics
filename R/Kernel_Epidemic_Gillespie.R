#'
#' Epidemic Gillespie w/ kernel set-up
#'

Kernel_Epidemic_Gillespie = function(N, a, gamma, beta, kernel){

  current_time = 0
  X = N - a
  Y = a
  Z = N - X - Y

  #' Track specific people
  S = 1:(N - a) #' Set of Susceptibles
  I = (N - a + 1):N #' Set of Infectives
  R = NULL #' Set of removed

  #' INfection Matrix
  B = kernel(beta)

  sim_data = matrix(NA, nrow = 2*N + 1, ncol = 6)
  sim_data[1, ] = c(current_time, X, Y, Z, NA, NA)
  no_event = 1
  while(Y > 0){

    old_time = current_time

    reduced_B = B[I,S]
    total_inf_pressure = sum(reduced_B)
    total_rem_rate = Y*gamma
    rate_next_event = total_rem_rate + total_inf_pressure
    time_to_next_event = rexp(1, rate_next_event)
    current_time = current_time + time_to_next_event

    #which_event = sample(c(0,1), size = 1, prob = c(total_inf_pressure, total_rem_rate)/rate_next_event)

    which_event = Heterogeneous_Event(reduced_B, total_rem_rate)

    if(which_event == 0){
      if(length(I) == 1){
        individual = I
      } else{
        individual = sample(I, size = 1)
      }
      Y = Y - 1
      Z = Z + 1
      I = I[I != individual]
      R = c(R, individual)
    } else{
      individual = S[which_event]
      S = S[-which_event]
      I = c(I, individual)
      X = X - 1
      Y = Y + 1
    }
    which_event = 1*(which_event == 0)
    sim_data[no_event + 1, ] = c(current_time, X, Y, Z, which_event, individual)
    no_event = no_event + 1
  }
  return(sim_data)
  #return(list(X = sim_data[,2], Y = sim_data[,3], Z = sim_data[,4], event_times = sim_data[,1], event_type = sim_data[,5], individual = sim_data[,6]))
}



Kerenel_Deterministic_Gillespie = function(N, initial_infective, beta, gamma, E, U, T_obs, k, store = TRUE, kernel){

  #' Initialise
  current_time = 0
  X = N - initial_infective
  Y = initial_infective
  Z = N - X - Y

  #' Panel Data
  if(length(T_obs) == 1){
    T_obs = c(0, T_obs)
  }

  obs_times = seq(T_obs[1], T_obs[2], length = k)

  #' Tracking Variables for Panel Data
  X_t0 = X
  Y_t0 = Y
  I_s = 0
  I_i = Y

  #' Data Storage
  if(store){
    sim_data = matrix(NA, nrow = 2*N + 1, ncol = 5)
    sim_data[1, ] = c(current_time, X, Y, Z, NA)
  } else{
    sim_data = NULL
  }

  panel_data = lapply(rep(NA, k), function(X) return(X))

  event_no = 1

  while(Y > 0){
    old_time = current_time
    reduced_B = B[I,S]
    total_inf_pressure = sum(reduced_B)
    total_rem_rate = Y*gamma
    rate_next_event = total_rem_rate + total_inf_pressure
    time_to_next_event = (1/rate_next_event)*E[event_no]
    current_time = current_time + time_to_next_event

    #' Recording Panel Data
    #' Record which observation times have been passed

    obs_times_passed = which(old_time <= obs_times & obs_times <= current_time)

    if(length(obs_times_passed) > 0){
      # Record Panel
      panel_data[[obs_times_passed[1]]] = c(n_ss = X, n_si = I_s, n_sr = X_t0 - (X + I_s), n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)

      # Reset
      X_t0 = X
      Y_t0 = Y
      I_i = Y
      I_s = 0

      for(i in obs_times_passed[-1]){
        panel_data[[i]] = c(n_ss = X, n_si = I_s, n_sr = X_t0 - (X + I_s), n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)
      }
    }

    which_event = Heterogeneous_Event(reduced_B, total_rem_rate, U[event_no])

    if(which_event == 0){
      X = X - 1
      I_s = I_s + 1
      Y = Y + 1
    } else if(which_event == 1){
      I_s = I_s - 1
      Y = Y - 1
    } else{
      I_i = I_i - 1
      Y = Y - 1
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
    panel_data[[i]] = c(n_ss = X, n_si = I_s, n_sr = X_t0 - (X + I_s), n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)
  }

  return(list(sim_data = sim_data, panel_data = panel_data, no_events = event_no))
}


#' Samples the event which occurs in an heterogeneous epidemic.
#' If a removal occurs, the funciton returns 0, which represents
#' the occurance of a removal
#' If an infection occurs, the index of the individual who becomes
#' infected is returned.

Heterogeneous_Event = function(infection_rate_matrix, removal_rate){
  if(length(infection_rate_matrix) == 0){
    return(0)
  } else if(is.null(dim(infection_rate_matrix))){
    individual_inf_rate = infection_rate_matrix
  } else{
    individual_inf_rate = colSums(infection_rate_matrix)
  }
  event = sample(c(1:length(individual_inf_rate), 0), size = 1,
                 prob = c(individual_inf_rate, removal_rate))
  return(event)
}

Heterogeneous_Event_Deterministic = function(infection_rate_matrix, removal_rate, U){
  if(length(infection_rate_matrix) == 0){
    return(0)
  } else if(is.null(dim(infection_rate_matrix))){
    individual_inf_rate = infection_rate_matrix
  } else{
    individual_inf_rate = colSums(infection_rate_matrix)
  }

  probs = c(individual_inf_rate, removal_rate)/sum(c(individual_inf_rate, removal_rate))

  probs_cumsum = cumsum(probs)
  more = TRUE
  index = 0
  while(more){
    index = index + 1
    more =  probs_cumsum[index] < U
  }
  event = index*(index != length(probs_cumsum))
  return(event)
}

