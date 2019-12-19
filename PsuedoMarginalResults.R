#' Psuedo Marginal Results

#' Simulate Epidemic

N = c(200, 500, 1000)
beta0 = 1
gamma = 0.5
beta = beta0/N
a = 0.05*N
p = c(0.1, 0.2, 0.3, 0.4, 0.5)
k = c(2, 4, 6, 8, 10, 15, 20)


#' ==== Epidemic N = 200 ====

set.seed(1)
U = runif(2*N[1] - a[1])
E = rexp(2*N[1] - a[1])
epidemic200 = SIR_Gillespie(N[1], a[1], beta[1], gamma = 0.5, obs_end = Inf, output = "not_panel",
                          U = U, E = E)
last_event_time = max(epidemic200$event_table$times)

#' Sample
i = 1
j = 4
rep_sample = c(sample(1:a[1], p[i]*a[1]), sample((a[1] + 1):N[1], p[i]*(N[1]- a[1])))
obs_times = seq(0, last_event_time, length = k[j])

Pseudo_Marginal_MCMC(N[1], panelData200, obs_times, beta[1], gamma, 0.001,
                     a[1], lambda = 0.01, V = diag(1, 2), no_sims = 50, no_its = 1000, burn_in = 0, parallel = T,
                     no_cores = 4)




#' ==== Epidemic N = 500 ====


kernel500 = Contact_Kernel(matrix(1, nrow = N[2], ncol = N[2]))

epidemic500 = Epidemic_Simulation(N[2], a[2], gamma, beta[2], kernel500)

#' ==== Epidemic N = 1000 ====

kernel1000 = Contact_Kernel(matrix(1, nrow = N[3], ncol = N[3]))

epidemic1000 = Epidemic_Simulation(N[3], a[3], gamma, beta[3], kernel1000)

#' ==== Panel Data & P-M MCMC N = 200 ====

#' Panel Data
panel_data.epidemics()











