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

Epidemic_Gillespie = function(N, initial_infective, gamma, beta, k, T_obs){

  #' Steps
  #' 1. Draw Exponential waiting time to next event
  #' 2. Sample which event it is to be (with appropriate weighting) (S -> I, I -> R (from S), I -> R (from I))
  #' 3. Update numbers

  #' Initialise
  time = 0
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
  sim_data = c(time, X, Y, Z)
  panel_data = lapply(rep(NA, k), function(X) return(X))
  #' Observation Counter (Corresponds to the current observation time we're looking out for)
  i = 1

  while(Y > 0){
    old_time = time
    rate_next_event = Y*gamma + X*Y*beta
    time_to_next_event = rexp(1, rate_next_event)
    time = time + time_to_next_event

    print(time)
    print(old_time)
    # Recording Panel Data
    if(old_time <= obs_times[i] & time >= obs_times[i] & i <= k){
      # Record Panel
      panel_data[[i]] = c(time = obs_times[i], n_ss = X, n_si = I_s, n_sr = X_t0 - (X + I_s), n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)

      # Reset
      X_t0 = X
      Y_t0 = Y
      I_i = Y
      I_s = 0

      # Look for next observation.
      i = i + 1
    }

    which_event = sample(c(0,1,2), size = 1, prob = c(X*Y*beta, I_s*gamma, I_i*gamma)/rate_next_event)

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



    sim_data = rbind(sim_data, c(time, X, Y, Z))

  }

  for(i in which(is.na(panel_data))){
    # Record Panel
    panel_data[[i]] = c(time = obs_times[i], n_ss = X, n_si = I_s, n_sr = X_t0 - (X + I_s), n_ii = I_i, n_ir = Y_t0 - I_i, n_rr = N - X_t0 - Y_t0)

    # Reset
    X_t0 = X
    Y_t0 = Y
    I_i = Y
    I_s = 0
  }
  return(list(sim_data = sim_data, panel_data = panel_data))
}











