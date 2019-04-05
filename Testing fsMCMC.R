#'
#' Testing Epidemic1 fsMCMC
#'


psi = 1 # Normalised infection rate
p = 0.05 # Proportion Infected
N = c(2*10^2, 10^3, 10^4, 10^5) # Population
a = p*N
beta = psi/N # Infection parameter
gamma = 0.15 # Removal Parameter
prop_obs = 0.1

#max(Epidemic4[,1], na.rm = T)

i = 1 # Which Epidemic
rep_sample = c(sample(1:(N[i] - a[i]), size = prop_obs*(N[i] - a[i])),
               sample((N[i] - a[i] + 1 ):N[i], size = prop_obs*a[i]))


T_obs = c(0,60)
k = 12
x_rep = gillespie_panel_transitions(X_0 = prop_obs*(N[i] - a[i]), Y_0 = prop_obs*a[i],
                                    Epidemic1[,5], Epidemic1[,1], Epidemic1[,6],
                                    T_obs, k, subset = rep_sample)

V = diag(1, 2)
lambda = 0.05
s = 400
no_its = 10000
burn_in = 1000
Test_rep = Epidemic_fsMCMC(N[i], a[i], x_rep, beta[i], gamma, no_draws = 2*N[i], s, T_obs, k, lambda, V,
                              no_its, burn_in)

random_sample = sample(1:N[i], size = prop_obs*N[i])


x_random = gillespie_panel_transitions(X_0 = sum(random_sample <= N[i] - a[i]),
                                       Y_0 = sum(random_sample > N[i] - a[i]),
                                       Epidemic1[,5], Epidemic1[,1], Epidemic1[,6],
                                       T_obs, k, subset = random_sample)

V = diag(1, 2)
lambda = 0.005
s = 200
no_its = 10000
burn_in = 1000
Test_random = Epidemic_fsMCMC(N[i], a[i], x_random, beta[i], gamma, no_draws = 2*N[i], s, T_obs, k, lambda, V,
                       no_its, burn_in)


T_obs = c(0,1.6)
k = 6
seq(T_obs[1], T_obs[2], length = k)
bad_sample_200 = random_sample
x_bad = gillespie_panel_transitions(X_0 = sum(random_sample <= N[i] - a[i]), Y_0 = sum(random_sample > N[i] - a[i]),
                                     Epidemic1[,5], Epidemic1[,1], Epidemic1[,6],
                                     T_obs, k, subset = bad_sample_200)


V = diag(1, 2)
lambda = 0.05
s = 250
no_its = 10000
burn_in = 1000
Test_bad = Epidemic_fsMCMC(N[i], a[i], x_bad, beta[i], gamma, no_draws = 2*N[i], s, T_obs, k, lambda, V,
                              no_its, burn_in)

