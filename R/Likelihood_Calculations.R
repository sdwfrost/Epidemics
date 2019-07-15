#'
#' Infection component of Heterogeneous Likelihood Functions
#' Original Author: Chris Jewell (2019)
#' Link: http://fhm-chicas-code.lancs.ac.uk/jewell/epilikelihoods/tree/master/R
#'


# ==== Infection Process ====

# == Contact Network Example ==


# == Network with Distances ==

# Assume and x and y coordinate for each individual

# coordinates will be a N x 2 matrix giving the x and y coordinates
# for each individual.

# The spdistsN1 function will calculate the distances between each



#  ==== Infection Product ====

# The product component of infection deals with who is
# exerting infectious pressure onto an individual at
# the time of infection and how much infectious pressure there is.

# WAIFW := Who aquired infection from whom
# A matrix which tells us who was exerting pressure on an individual
# when they became infected.

prod_part_inf = function(t_inf_j, events, B, log = TRUE){
  is_infected = t_inf_j < Inf
  waifw = sapply(t_inf_j, function(t) events[is_infected, 1] < t & t < events[is_infected, 2])
  lambdaj = colSums(B[is_infected, is_infected] * waifw[,is_infected])
  t_inf_0 = which(t_inf_j[is_infected] == min(t_inf_j[is_infected]))
  if(log){
    sum(log(lambdaj[-t_inf_0]))
  } else{
    prod(lambdaj[-t_inf_0])
  }
}

# ==== Infection Integral ====

#' Relies on the period of time for which individual j exerts pressure on individual j
#' Denote this as the exposure time.

#' @func interval_intersect, calculates the llength of intersection of two intervals
#' @param interval_i
#' @param interval_j
#' @return The length of intersection between the respective intervals in interval_i and interval_j

interval_intersect = function(interval_i, interval_j){
  interval_start = sapply(interval_j[,1], function(x) pmax(x, interval_i[,1]))
  interval_end = sapply(interval_j[,2], function(x) pmin(x, interval_i[,2]))
  pmax(interval_end - interval_start, 0)
}


#' @func integral_part, calculates the integral of infectious pressure portion of an epidemic
#'                      likeihood function
#'
#'
#'
#'
#'

integral_part_inf = function(B, events, t_inf_j, with.beta = TRUE){
  i_infected = events[,1] < Inf
  E = interval_intersect(events[i_infected,], cbind(min(t_inf_j), t_inf_j))
  if(with.beta){
    integral  = E*B[i_infected,]
    sum(integral)
  } else{
    integral = E
    sum(integral)
  }
}

#' @func log_likelihood_inf, calculates the infectious process part of an epidemic likelihood

log_likelihood_inf = function(par, t_inf, t_rem, kernel){
  B = kernel(par)
  prod = prod_part_inf(t_inf, cbind(t_inf, t_rem), B)
  integral = integral_part_inf(B, cbind(t_inf, t_rem), t_inf)
  return(prod - integral)
}

# ==== Removal Process (Gamma Infectious Period) ====

#' Computes the Removal portion of an epidemic likelihood
#'
#' @param par parameter values to be passed to gamma density function (par[1] rate parameter, par[2] shape parameter)
#' @param t_inf infection times of individuals involved in the epidemic
#' @param t_rem removal times of individuals involved in the epidemic

log_likelihood_rem = function(par, t_inf, t_rem){
  i_removed = which(t_rem < Inf)
  llh = sum(dgamma(t_rem[i_removed] - t_inf[i_removed], rate = par, shape = 1, log = TRUE))
  return(llh)
}


# ==== Full Likelihood ====

#' Computes the full epidemic likelihood

epidemic_loglikelihood = function(par_inf, par_rem, t_inf, t_rem, kernel){
  return(log_likelihood_rem(par_rem, t_inf, t_rem) + log_likelihood_inf(par_inf, t_inf, t_rem, kernel))
}



