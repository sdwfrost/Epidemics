#'
#' Epidemics of varying popn size but constant R_0
#'
#'

lambda = 1 # Normalised infection rate
p = 0.05 # Proportion Infected
N = c(10^2, 10^3, 10^4, 10^5, 10^6) # Population
a = p*N
beta = lambda/N # INfection parameter
gamma = 0.15 # Removal Parameter
k = 12
T_obs = c(0,10)

Epidemic1 = Epidemic_Gillespie(N[1], a[1], gamma, beta[1], k, T_obs)

Epidemic2 = Epidemic_Gillespie(N[2], a[2], gamma, beta[2], k, T_obs)

Epidemic3 = Epidemic_Gillespie(N[3], a[3], gamma, beta[3], k, T_obs)

Epidemic4 = Epidemic_Gillespie(N[4], a[4], gamma, beta[4], k, T_obs)

Epidemic5 = Epidemic_Gillespie(N[5], a[5], gamma, beta[5], k, T_obs)












