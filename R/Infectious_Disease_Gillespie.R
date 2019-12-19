#'
#' Infectious_Disease_Gillespie function
#'
#' Wrapper function for various Infectious Disease Process Gillespie algorithms.
#'
#' Takes the usual arguments plus an extra argument to choose which process you
#' want to run i.e SIR/SIS/SEIR
#'

Infectious_Disease_Gillespie = function(type = "SIR", N, a, beta, gamma, obs_end = Inf, kernel = NULL, U, E,
                                        output = "event", obs_times = NULL){
  if(type == "SIR"){
    return(SIR_Gillespie(N, a, beta, gamma, trans = possible_transitions(1:3),
                         obs_end = obs_end, kernel, U, E, output, obs_times))
  } else if(type == "SIS"){
    return(SIS_Gillespie(N, a, beta, gamma, trans = possible_transitions(1:2),
                         obs_end, kernel, U, E, output, obs_times))
  }
}
