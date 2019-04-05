#'
#'
#' Varying number of (uniform) observations, k
#'
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

set.seed(2)
rep_sample = c(sample(1:(N[i] - a[i]), size = prop_obs*(N[i] - a[i])),
               sample((N[i] - a[i] + 1 ):N[i], size = prop_obs*a[i]))




T_obs = c(0,60)

# ==== Working Beta Samples ====
jpeg("Enough_Infectious_Process_Obs.jpeg", width = 960, height = 960)
k = 12
x_rep = gillespie_panel_transitions(X_0 = prop_obs*(N[i] - a[i]), Y_0 = prop_obs*a[i],
                                    Epidemic1[,5], Epidemic1[,1], Epidemic1[,6],
                                    T_obs, k, subset = rep_sample)

V = diag(1, 2)
lambda = 0.03
s = 200
no_its = 10000
burn_in = 1000
Test_rep = Epidemic_fsMCMC(N[i], a[i], x_rep, beta[i], gamma, no_draws = 2*N[i], s, T_obs, k, lambda, V,
                           no_its, burn_in)
dev.off()

# ==== Not Working Beta Samples ====

jpeg("Not_Enough_Epidemic_Process_Obs.jpeg", width = 960, height = 960)
k = 2
x_rep = gillespie_panel_transitions(X_0 = prop_obs*(N[i] - a[i]), Y_0 = prop_obs*a[i],
                                    Epidemic1[,5], Epidemic1[,1], Epidemic1[,6],
                                    T_obs, k, subset = rep_sample)

V = diag(1, 2)
lambda = 0.2
s = 400
no_its = 10000
burn_in = 1000
Test_rep = Epidemic_fsMCMC(N[i], a[i], x_rep, beta[i], gamma, no_draws = 2*N[i], s, T_obs, k, lambda, V,
                           no_its, burn_in)
dev.off()


# ==== How is mixing affected by increasing/decreasing k ====

