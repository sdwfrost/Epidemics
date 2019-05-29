#'
#' fsMCMC_sim
#' Simulating a realisation of infectious process with current state
#' of the MCMC and calculating the probability of the observing the
#' partial panel data x given the full panel data Y
#'

#' @param x observed panel data
#' @param N size of population
#' @param a initial size of infectives
#' @param beta infectious process parameters
#' @param gamma removal/recovery process parameters
#' @param obs_end When observation of the infectious process ends
#' @param kernel function specifying the dynamics of the infectious process
#' @param U set of Uniform(0,1) random varibales which determine the nature of events
#' @param E set of Exponential(1) random variables which determine when events happen
#' @param T_obs Period for which the infectious process is observed.
#' @param k How many point observations are made during T_obs

fsMCMC_sim = function(type, x, N, a, beta, gamma, kernel, U, E, obs_times){
  prop_obs = sum(x[[1]])
  #Y = SIR_Gillespie(N, a, beta, gamma, obs_end = tail(obs_times, n = 1), kernel = kernel, U = U, E = E,
  #                  obs_times = obs_times, output = "panel")$panel_data
  Y = Infectious_Disease_Gillespie(type, N, a, beta, gamma, obs_end = tail(obs_times, n = 1), kernel = kernel, U = U, E = E,
                                   output = "panel", obs_times = obs_times)$panel_data
  logP = dHyperGeom(x, Y, prop_obs, log = T)
  return(logP)
}


