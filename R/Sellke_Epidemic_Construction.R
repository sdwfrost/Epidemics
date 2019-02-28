#'
#'
#' fsMCMC for Epidemic Panel Data
#'

#' We have a partial panel data observations x* of an epidemic
#'
#' The epidemic was observed in a subpopulation of m people, at k equaly spaced timepoints
#' between time 0 and time T
#'
#' The number of Susceptible, Infected and Removed are observed at each timepoint
#' From this the following can be deduced: The number of individuals who ...
#'    ... remained susceptible
#'    ... became infected being susceptible last time they were observed
#'    ... became removed  '                                             '
#'
#'    ... remained infected
#'    ... became removed after being infected last time they were observed
#'
#'    ... remained removed

#' fsMCMC will use data augmentation to make inference on the parameters of the Epidemic
#' process
#'
#' Non-Centered Selke Construction of an epidemic.
#' A combination of parameter values theta and augmented data Z give the observation of an epidemic
#' Then the probability of observing x* given the constructed empidemic can be calculated and used in
#' an Accept/Reject step
#'

#' SELLKE CONSTRUCTION
#'
#' Infectious period Q_i ~ Exp(gamma) = (in D) U_i/gamma ~ Exp(1) i = 1,...,N
#'
#' Each individual imposes total infectious pressure beta*Q_i over their whole
#' infectious period.
#'
#' Infectious Thresholds T_i ~ Exp(1/N), N = closed popn size i = 1,...,N
#'
#' Once this Infectious Threshold is exceeded individual i becomes infected
#'
#' However, by ordering these threshold times, they can be split up as follows
#'
#'  ordered_T_i = sum_1^(i-1) L_j, L_j ~ Exp((N-j)/N)
#'
#' So Epidemic is decided by simulating L_1,...,L_N and U_1,...,U_N and beta and gamma.
#'
#' Neal & Huang used this to construct an epidemic and calculate its final size m
#'
#' I want the Infection Times, which can then be combined with the infectious periods Q = U/gamma
#' to give the removal times. Then can use transitions_between_observations to obtain N_ss, N_si,..., N_rr.
#'
#' This is definitely possible.
#'
