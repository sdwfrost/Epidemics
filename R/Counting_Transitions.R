#'
#'
#'


#' @param interval, a time interval (e.g (0, 2)) to
#'                  look at changes of state.
#' @param event_times, a N by (no.states) matrix denoting
#'                     the time at which each individual
#'                     transitions to another state.
#' @param states, labels for the states (optional)
#'
interval_transitions = function(interval, event_times){
  # Transitions from susceptible state

  susc_at_t1 = event_times[,1] > interval[1]
  S_t1 = sum(susc_at_t1)

  susc_at_t2 = event_times[,1] > interval[2]
  S_t2 = sum(susc_at_t2)

  infected_before_t2 =  event_times[,1] < interval[2]
  removed_after_t2 = event_times[,2] > interval[2]
  n_si = sum(susc_at_t1 & infected_before_t2 &
             removed_after_t2)

  n_sr = (S_t1 - S_t2) - n_si # The number of individual moving out of S, not including those who just went to I
  n_ss = S_t1 - n_si - n_sr # Those moving from S to I and R between t_1 and t_2 are removed from S.

  # Transitions from infectious state
  infectious_at_t1 = event_times[,1] <= interval[1] & event_times[,2] > interval[1]
  I_t1 = sum(infectious_at_t1)

  removed_before_t2 = event_times[,2] < interval[2]
  n_ir = sum(infectious_at_t1 & removed_before_t2)

  n_ii = I_t1 - n_ir

  # Transitions from removal state
  n_rr = sum(event_times[,2] < interval[1])

  return(c(n_ss = n_ss, n_is = 0, n_rs = 0, n_si = n_si, n_ii = n_ii, n_ri = 0, n_sr = n_sr, n_ir = n_ir, n_rr = n_rr))
}



#' @param k, number of evenly distributed observations
#' @param T_obs, total time which an epidemic is observed for.
#' @param event_times, the times at which transitions between
#'                     epidemic states occur for each individual.

transitions_between_observations = function(T_obs, k, event_times){
  if(length(T_obs) == 1){
    obs_times = seq(0, T_obs, length = k)
  } else{
    obs_times = seq(T_obs[1], T_obs[2], length = k)
  }
  intervals = lapply(1:(length(obs_times)-1), function(i) cbind(obs_times[i], obs_times[i+1]))
  X = lapply(intervals, interval_transitions, event_times)
  return(X)
}

#' Function for getting transition data from the output of a Gillespie style algorithm
#'
#'
#'

gillespie_transitions = function(X_0, Y_0, gillespie_events, T_obs, k, subset = NULL){

  X = X_0
  Y = Y_0
  X_t0 = X_0
  Y_t0 = Y_0
  I_i = Y_0
  I_s = 0


  if(is.null(subset)){
    events = gillespie_events$event_type
    times = gillespie_events$event_times
  } else{
    relavent_events = sort(unlist(sapply(subset, function(X) return(which(X == gillespie_events$individual)))))
    events = gillespie_events$event_type[relavent_events]
    times = gillespie_events$event_times[relavent_events]
  }


  N = X_0 + Y_0
  current_time = 0

  obs_times = seq(T_obs[1], T_obs[2], length = k)

  panel_data = lapply(rep(NA, k), function(X) return(X))
  panel_data[[1]] = c(n_ss = X_0, n_si = I_s, n_sr = X_t0 - (X_0 + I_s), n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)

  event_no = 1
  while(Y > 0){
    old_time = current_time
    current_time = times[event_no]

    obs_times_passed = which(old_time <= obs_times & obs_times <= current_time)

    if((length(obs_times_passed) > 0)){
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
    if(events[event_no] == 0){
      I_s = I_s + 1
      Y_t0 = Y_t0 + 1
      X = X - 1
    } else if(events[event_no] == 1){
      I_s = I_s - 1
      Y_t0 = Y - 1
    } else{
      I_i = I_i - 1
      Y = Y - 1
    }
    event_no = event_no + 1
  }

  NA_panels = which(is.na(panel_data))

  for(i in NA_panels){
    panel_data[[i]] = c(n_ss = X, n_si = I_s, n_sr = X_t0 - (X + I_s), n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)
  }

  return(panel_data)
}


# Observation 1
gillespie_transitions_simplify = function(X_0, Y_0, gillespie_events, T_obs, k, subset = NULL){

  X = X_0
  Y = Y_0
  X_t0 = X_0
  Y_t0 = Y_0
  I_i = Y_0
  I_s = 0


  if(is.null(subset)){
    events = gillespie_events$event_type
    times = gillespie_events$event_times
  } else{
    relavent_events = sort(unlist(sapply(subset, function(X) return(which(X == gillespie_events$individual)))))
    events = gillespie_events$event_type[relavent_events]
    times = gillespie_events$event_times[relavent_events]
  }

  N = X_0 + Y_0
  obs_times = seq(T_obs[1], T_obs[2], length = k)

  panel_data = lapply(rep(NA, k), function(X) return(X))
  panel_data[[1]] = c(n_ss = X_0, n_si = I_s, n_sr = X_t0 - (X_0 + I_s), n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)

  for(i in 2:length(obs_times)){
    time_indices = which(obs_times[i-1] < times & times < obs_times[i])
    I_s = sum(events[time_indices] == 0) - sum(events[time_indices] == 1)
    I_i = I_i - sum(events[time_indices] == 2)
    Y = Y_t0 - sum(events[time_indices] == 1) - sum(events[time_indices] == 2)
    X = X_t0 - sum(events[time_indices] == 0)
    panel_data[[i]] = c(n_ss = X, n_si = I_s, n_sr = X_t0 - (X + I_s), n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)

    #' Reset
    Y_t0 = Y
    X_t0 = X
    I_s = 0
    I_i = Y
  }
  return(panel_data)
}


