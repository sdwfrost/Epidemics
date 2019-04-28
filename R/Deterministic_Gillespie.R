#'
#'  Deterministic Gillespie
#'


#' Two Options
#'
#' 1. Draw a predetermined amount of exp(1) and unif(0,1) variables
#'
#' Can predict how many you are likely to need in simulation based settings
#' but still might not have enough. Maybe draw excessive amounts but then wasting
#' time and memory
#'
#'
#' 2. Select random seed(s) to follow a particular stream of draws.
#'
#' Don't need to worry about how many draws you need. Less applicable to MCMC updates
#'
#' Could do MCMC updates by doing chuncks of RVs generated from different random seeds
#' then redrawing a seed for one chunck. Would this work?
#'


#' Deterministic_Gillespie1 A Gillespie algorithm which is deterministic once parameters and
#'                          and RVs have been supplied
#' @param N Size of closed population
#' @param initial_infective How many of the population are initially infected
#' @param beta Infection Paramete, the rate at which an infectious individual infects a susceptible individual
#' @param gamma Removal Parameter, the rate at which an individual is removed once being infected
#' @param E, Exponential Random Variables ith mean 1, which determine the time to the next event when combined
#'           parameters and current state of the epidemic.
#' @param U, Uniform(0,1) random variables, which determine what type of event occurs when combined with
#'           parameters and current state of the epidemic.
#' @param T_obs, Interval of time for which the epidemic is observed for in practice. If a length of time is
#'               provided, it is assumed observation starts at time 0.
#' @param k, How many equally spaced panel observations are made between the start and end of the observation
#'           interval.

Deterministic_Gillespie1 = function(N, a, beta, gamma, E, U, T_obs, k, store = TRUE){

  #' Initialise
  current_time = 0
  X = N - a
  Y = a
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
    rate_next_event = Y*gamma + X*Y*beta
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

    P_1 = X*Y*beta/rate_next_event
    P_2 = I_s*gamma/rate_next_event

    which_event = sum(U[event_no] > c(P_1, P_1 + P_2))

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


#' Deterministice_Gillespie2 A function which constructs an epidemic based on parameters and Random Variables
#'                           predetermined by a seed.
#'

Deterministic_Gillespie2 = function(N, initial_infective, beta, gamma, seeds, no_draws_per_seed, k, T_obs){

  #' Initialise
  current_time = 0
  X = N - initial_infective
  Y = initial_infective
  Z = N - X - Y

  #' Panel Data
  if(length(T_obs) == 1){
    T_obs = c(0, T_obs)
  }

  obs_times = seq(T_obs[1], T_obs[2], length = k)[-1]

  #' Tracking Variables for Panel Data
  X_t0 = X
  Y_t0 = Y
  I_s = 0
  I_i = Y

  #' Data Storage
  sim_data = c(current_time, X, Y, Z)
  panel_data = lapply(rep(NA, k - 1), function(X) return(X))
  #' Observation Counter (Corresponds to the current observation time we're looking out for)
  i = 1
  event_no = 0
  seed_no = 1
  while(Y > 0){
    event_no = event_no + 1
    if((event_no - 1)%%no_draws_per_seed == 0){
      set.seed(seeds[seed_no])
      seed_no = seed_no + 1
    }

    old_time = current_time
    rate_next_event = Y*gamma + X*Y*beta
    E = rexp(1)
    time_to_next_event = (1/rate_next_event)*E
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

    U = runif(1)
    P_1 = X*Y*beta/rate_next_event
    P_2 = I_s*gamma/rate_next_event

    which_event = sum(U > c(P_1, P_1 + P_2))

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

    sim_data = rbind(sim_data, c(current_time, X, Y, Z))

  }

  NA_panels = which(is.na(panel_data))

  for(i in NA_panels){
    panel_data[[i]] = c(n_ss = X, n_si = I_s, n_sr = X_t0 - (X + I_s), n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)
  }
  rm(.Random.seed, envir=.GlobalEnv)
  return(list(sim_data = sim_data, panel_data = panel_data, no_events = event_no))
}

