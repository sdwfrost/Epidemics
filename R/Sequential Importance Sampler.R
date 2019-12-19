#' Sequential Importance Sampler


#' Vanilla Importance Sampler samples N_p random particles
#' and extracts panel data which lines up with the observed
#' panel data.
#'
#' Then (normalised) Importance weights are calculated.
#' Ideal value of each weight is 1/N_p, giving an
#' effective sample size of N_p.
#'
#' The magnitude of the importance weights deteriotes
#' as more panels are added to the observation process
#' as the simulated particles must align with more observed
#' data. However, we would like more panels as we gain more
#' information about the epidemic process. Hence, we need a
#' method which reduces the number of zero-valued
#' importance weights for a higher number a panels. Or
#' better yet, an importance sampling scheme which is
#' independent of panel number.
#'
#' The goal of sequential Importance Sampling is to reduce this
#' burden by simulating particles up to each observation point,
#' calculating importance wieghts, then simulating each particle forward
#' and adjusting the importance weight accordingly until all observation
#' times have been simulated too.
#'
#' Although this seems to reduce the problem to k, 1 panel simulations, which are the most "efficient",
#' the problem has not been fixed. This is because weak particles will still
#' fail to be valid.
#'
#'
#'
#' epidemicSequentialImportanceSampler <- function(panelData, obsTimes, theta, Nparticles,
#'                                                 resample = F){
#'
#'   m = sum(panelData[1])
#'
#'   #' 1. Simulate N particles up to the first observation time. $x^{i}$ \sim q() (This is a particle)
#'   #'
#'   #' First Simulation
#'
#'   i = 1
#'   particles = lapply(X = 1:N_p, function(X){
#'     homogeneousPanelDataSIR_Gillespie(initialState, theta[1], theta[2])
#'     SIR_Gillespie(N, initial_state = , beta = theta[1], gamma = theta[2], output = "panel",
#'                   obs_times = obs_times[(i-1):i], obs_end = obs_times[i])$panel_data
#'   })
#'
#'   #' Calculate (Normalised) Importance weights for these simulations
#'
#'   #' extract transition data from particles
#'   particle_trans = lapply(particles, transition_data)
#'   y_given_x = sapply(X = particle_trans, function(X) extraDistr::dmvhyper(y[[1]], X[[1]], k = m))
#'
#'   ISweights = y_given_x/sum(y_given_x)
#'
#'   return(list(ESS = ESS, ISweights = ISweights, particles = particles, y = y))
#'
#'
#'   #' 4. Use this weighted sample to estimate integrals involving
#'   #'    posterior distribution.
#'
#'
#'
#'
#'
#' }
#'

