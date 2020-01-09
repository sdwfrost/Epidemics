#'
#' Gillespie Algorithm for SIR Epidemic
#'

#' Homogenous SIR Simulator
#' @param N size of closed population
#' @param initial infective, Number of individuals who are initially infected
#' @param par_gamma, vector of rate of removal of infectives and shape parameter for removal distribution
#' @param beta, rate of infection of individuals
#' @param obs, if TRUE,
#'
#'


Epidemic_Gillespie = function(N, a, gamma, beta){

  current_time = 0
  X = N - a
  Y = a
  Z = N - X - Y

  #' Track specific people
  S = 1:(N - a) #' Set of Susceptibles
  I = (N - a + 1):N #' Set of Infectives
  R = NULL #' Set of removed

  sim_data = matrix(NA, nrow = 2*N + 1, ncol = 6)
  sim_data[1, ] = c(current_time, X, Y, Z, NA, NA)
  no_event = 1
  while(Y > 0){

    old_time = current_time
    rate_next_event = Y*gamma + X*Y*beta
    time_to_next_event = rexp(1, rate_next_event)
    current_time = current_time + time_to_next_event

    which_event = sample(c(0,1), size = 1, prob = c(X*Y*beta, Y*gamma)/rate_next_event)

    if(which_event == 0){
      if(length(S) == 1){
        individual = S
      } else{
        individual = sample(S, size = 1)
      }
      X = X - 1
      Y = Y + 1
      S = S[S != individual]
      I = c(I, individual)
    }
    if(which_event == 1){
      if(length(I) == 1){
        individual = I
      } else{
        individual = sample(I, size = 1)
      }
      Y = Y - 1
      Z = Z + 1
      I = I[I != individual]
      R = c(R, individual)
    }

    sim_data[no_event + 1, ] = c(current_time, X, Y, Z, which_event, individual)
    no_event = no_event + 1
  }
  return(sim_data)
  #return(list(X = sim_data[,2], Y = sim_data[,3], Z = sim_data[,4], event_times = sim_data[,1], event_type = sim_data[,5], individual = sim_data[,6]))
}

gillespie_panel_transitions = function(X_0, Y_0, events, times, individuals, T_obs, k, subset = NULL){

  X_t0 = X_0
  Y_t0 = Y_0
  N = X_0 + Y_0
  I_s = 0
  I_i = Y_t0

  times = times[!is.na(times) & times != 0]
  events = events[!is.na(events)]
  individuals = individuals[!is.na(individuals)]
  if(1 - is.null(subset)){
    relevent_events = sort(unlist(sapply(subset, function(X) return(which(X == individuals)))))
    events = events[relevent_events]
    times = times[relevent_events]
    individuals = individuals[relevent_events]
  }


  #' Panel Data
  if(length(T_obs) == 1){
    T_obs = c(0, T_obs)
  }

  obs_times = seq(T_obs[1], T_obs[2], length = k)

  panel_data = lapply(rep(NA, k), function(X) return(X))

  #if(obs_times[1] == 0){
  #  panel_data[[1]] = c(n_ss = X_t0, n_si = I_s, n_sr = X_t0 - (X_0 + I_s), n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)
  #}

  for(i in 1:length(obs_times)){
    if(i == 1){
      time_indicies = which(times < obs_times[i])
    } else{
      time_indicies = which(times < obs_times[i] & times > obs_times[i-1])
    }


    #' At the start of the period, Individual is susceptible but becomes infected (but
    #' not removed) before the next observation
    infections = sum(events[time_indicies] == 0)

    #' At the start of the period, Individual is infected and is removed before
    #' next observation.
    no_events_individuals = sapply(subset, function(X) sum(individuals[time_indicies] == X))
    removals_from_susc = sum(no_events_individuals == 2)

    #' Everyone who is removed in this period
    removals = sum(events[time_indicies] == 1)

    #' At the start of the period, Individual is susceptible but is removed before
    #' next observation
    #'

    removals_from_inf = removals - removals_from_susc

    X_t1 = X_t0 - infections
    Y_t1 = Y_t0 + infections - removals
    I_s = infections - removals_from_susc
    I_i = Y_t0 - removals_from_inf

    panel_data[[i]] = c(n_ss = X_t1, n_si = I_s, n_sr = removals_from_susc, n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)

    X_t0 = X_t1
    Y_t0 = Y_t1
  }
  return(panel_data)
}



Epidemic_Gillespie_transisitions = function(N, initial_infective, gamma, beta, k, T_obs, store = TRUE){

  #' Steps
  #' 1. Draw Exponential waiting time to next event
  #' 2. Sample which event it is to be (with appropriate weighting) (S -> I, I -> R (from S), I -> R (from I))
  #' 3. Update numbers

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

  #' Track specific people
  S = 1:(N - initial_infective) #' Set of Susceptibles
  I = (N - initial_infective + 1):N #' Set of Infectives
  R = NULL #' Set of removed
  I_from_S = NULL
  I_from_I = I
  #' Data Storage
  if(store){
    sim_data = matrix(NA, nrow = 2*N + 1, ncol = 6)
    sim_data[1, ] = c(current_time, X, Y, Z, NA, NA)
  } else{
    sim_data = NULL
  }

  panel_data = lapply(rep(NA, k), function(X) return(X))

  panel_data[[1]] = c(n_ss = X, n_si = I_s, n_sr = X_t0 - (X + I_s), n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)
  #' sim_data = c(current_time, X, Y, Z)
  #'if(panel){
  #'
  #'} else{
  #'  panel_data = NULL
  #'}
  #' Observation Counter (Corresponds to the current observation time we're looking out for)
  event_no = 1

  while(Y > 0){
    old_time = current_time
    rate_next_event = Y*gamma + X*Y*beta
    time_to_next_event = rexp(1, rate_next_event)
    current_time = current_time + time_to_next_event

    #' Recording Panel Data
    #' Record which observation times have been passed

    obs_times_passed = which(old_time <= obs_times & obs_times <= current_time)

    if((length(obs_times_passed) > 0)){
      # Record Panel
      panel_data[[obs_times_passed[1]]] = c(n_ss = X, n_si = I_s, n_sr = X_t0 - (X + I_s), n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)

      # Reset
      X_t0 = X
      Y_t0 = Y
      I_i = Y
      I_s = 0
      I_from_S = NULL
      I_from_I = I

      for(i in obs_times_passed[-1]){
        panel_data[[i]] = c(n_ss = X, n_si = I_s, n_sr = X_t0 - (X + I_s), n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)
      }
    }

    which_event = sample(c(0,1,2), size = 1, prob = c(X*Y*beta, I_s*gamma, I_i*gamma)/rate_next_event)


    if(which_event == 0){

      if(X == 1){
        individual = S
      } else{
        individual = sample(c(S), size = 1) # Which Susceptible Indv is infected
      }
      S = S[!(S == individual)] # Remove from Susceptible set
      X = X - 1  # Decrease Number of Susceptibles

      I_from_S = c(I_from_S, individual)
      I_s = I_s + 1 # Increase number of Infected which have come from susceptible

      I = c(I, individual)
      Y = Y + 1
    } else if(which_event == 1){

      if(I_s == 1){
        individual = I_from_S
      }  else{
        individual = sample(I_from_S, size = 1)
      }

      I_from_S = I_from_S[!(I_from_S == individual)]
      I_s = I_s - 1

      I = I[!(I == individual)]
      Y = Y - 1

      R = c(R, individual)
      Z = Z + 1
    } else{

      if(I_i == 1){
        individual = I_from_I
      } else{
        individual = sample(c(I_from_I), size = 1)
      }
      I_from_I = I_from_I[!(I_from_I == individual)]
      I_i = I_i - 1

      I = I[!(I == individual)]
      Y = Y - 1

      R = c(R, individual)
      Z = Z + 1
    }

    #' Store old time
    #' Update extra parameters
    Z = N - X - Y

    if(store){
      sim_data[event_no + 1,] = c(current_time, X, Y, Z, which_event, individual)
    }
    event_no = event_no + 1
  }

  #if(panel){

  NA_panels = which(is.na(panel_data))

  for(i in NA_panels){
    panel_data[[i]] = c(n_ss = X, n_si = I_s, n_sr = X_t0 - (X + I_s), n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)
  }

  return(list(X = sim_data[,2], Y = sim_data[,3], Z = sim_data[,4], event_times = sim_data[,1], event_type = sim_data[,5], individual = sim_data[,6], panel_data = panel_data))
}


