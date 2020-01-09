#'
#' PhD
#' General Stochastic Epidemic Simulation
#'
#'

# ==== Preamble ====

#' @func next_infection returns the individual in a current epidemic to be infected next and
#'                      the individual who infected them and the time at which it happened.
#' @param inf_period,  The time which infection lasts for individual i.
#' @param contact_times, time at which and infectious individual i makes infectious contact
#'                       with a susceptiable individual j.
#' @param t_inf, time at which infection takes place for each individual.
#' @param state, State which each individual is in (S = 0, I = 1, R = 2)

next_infection = function(inf_period, contact_times, t_inf, which_infected, which_susceptible){

  # Infection times matrix
  t_inf_mat <- matrix(t_inf[which_infected], nrow = length(which_infected), ncol = length(which_susceptible), byrow = FALSE)

  # Potential times of next infections
  future_t_inf <- t_inf_mat + contact_times[which_infected, which_susceptible]


  next_t_inf = min(future_t_inf)
  indices = which(future_t_inf == next_t_inf, arr.ind = TRUE)

  return(list(next_t_inf = next_t_inf, infector = which_infected[indices[1,1]], infected = which_susceptible[indices[1,2]]))

}


new_states = function(state, next_t_inf, infector, infected, contact_times){
  if(next_t_inf  == Inf){
    state[which(state == 1)] = 2
  } else{
    state[infector]  =  2 - 1*(sum(contact_times[infector, which(state == 0)] < Inf) > 1)
    state[infected] = 1
  }
  return(state)
}

#' @func Epidemic_Simulation, simulates an epidemic
#' @param N, size of the closed population
#' @param a, Number of initial infectious individuals
#' @param par_beta, parameters to be passed to the kernel in order to calculate infection process parameters
#' @param par_gamma, parameters to be passed to functions concerning the infectious period distribution
#' @param kernel, model for the infectious process of the epidemic

Epidemic_Simulation <- function(N, a, par_gamma, par_beta, kernel){

  #' Choose at random who is initially infected
  #initialInfected <- sample(1:N, a, replace = F)

  state <- c(rep(1, a), rep(0, N - a))
  t_inf <- c(rep(0 , a), rep(Inf, N - a))
  #state = rep(0 , N)
  #t_inf <- rep(0, N) # Time at which an individual is infected.

  #state[initialInfected] <- 1
  #t_inf [initialInfected] <- NA

  inf_period <- rgamma(N, rate = par_gamma, shape = 1) # Period of time for which an individual is infected

  B = kernel(par_beta)

  contact_times <- matrix(rexp(N^2, rate = B), nrow = N, ncol = N, byrow = FALSE) # Time, relative to time of infection, an individual i has infectious contact with individual j

  diag(contact_times) <- Inf # Cannot have infectious contact with self

  # Handle contact times which will not be meaningful.
  inf_period_mat = matrix(rep(inf_period, N), nrow = N, ncol = N, byrow = FALSE)
  contact_times = ifelse(contact_times > inf_period_mat, Inf, contact_times)

  # Final Size Counter
  final_size <- 0

  # Check if
  if(min(contact_times[1,]) == Inf){
    return(list(N = N, final_size = final_size, t_inf = c(rep(0, a), rep(Inf, N - a)),
                t_rem = c(inf_period[1:a], rep(Inf, N - a))))
  }

  while(sum(state == 0) > 0 & sum(state == 1) > 0 & sum(state == 1) < N){

    # Set of susceptibles
    which_susceptible <- which(state == 0)
    # Set of Infected
    which_infected <- which(state == 1)

    next_inf = next_infection(inf_period, contact_times, t_inf, which_infected, which_susceptible)



    t_inf[next_inf$infected] = next_inf$next_t_inf
    state = new_states(state, next_inf$next_t_inf, next_inf$infector, next_inf$infected, contact_times)
    final_size = final_size + 1*(next_inf$next_t_inf < Inf)

  }

  t_rem <- t_inf + inf_period
  t_inf[!(t_inf >= 0) | is.na(t_inf)] <- Inf
  t_rem[!(t_rem >= 0) | is.na(t_inf)] <- Inf

  return(list(N = N, final_size = final_size, t_inf = t_inf, t_rem = t_rem,
              par_beta = par_beta, par_gamma = par_gamma, kernel = kernel))
}



