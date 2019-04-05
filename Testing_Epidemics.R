#'
#' Epidemics of varying popn size but constant R_0
#'
#'

psi = 1 # Normalised infection rate
p = 0.05 # Proportion Infected
N = c(2*10^2, 10^3, 10^4, 10^5, 10^6) # Population
a = p*N # Number of Initial Infectives
beta = psi/N # Infection parameter
gamma = 0.15 # Removal Parameter

Epidemic1 = Epidemic_Gillespie(N[1], a[1], gamma, beta[1])
usethis::use_data(Epidemic1, overwrite = T)

Epidemic2 = Epidemic_Gillespie(N[2], a[2], gamma, beta[2])
usethis::use_data(Epidemic2, overwrite = T)

Epidemic3 = Epidemic_Gillespie(N[3], a[3], gamma, beta[3])
usethis::use_data(Epidemic3, overwrite = T)

Epidemic4 = Epidemic_Gillespie(N[4], a[4], gamma, beta[4])
usethis::use_data(Epidemic4, overwrite = T)

#Epidemic5 = Epidemic_Gillespie(N[5], a[5], gamma, beta[5])
#usethis::use_data(Epidemic5, overwrite = T)










