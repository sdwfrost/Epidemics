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

  return(c(n_ss = n_ss, n_si = n_si, n_sr = n_sr, n_ii = n_ii, n_ir = n_ir, n_rr = n_rr))
}



#' @param k, number of evenly distributed observations
#' @param T_obs, total time which an epidemic is observed for.
#' @param event_times, the times at which transitions between
#'                     epidemic states occur for each individual.

transitions_between_observations = function(T_obs, k, event_times){
  obs_times = seq(0, T_obs, length = k)
  intervals = lapply(1:(length(obs_times)-1), function(i) cbind(obs_times[i], obs_times[i+1]))
  X = lapply(intervals, interval_transitions, event_times)
  return(X)
}

