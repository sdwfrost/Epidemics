#' fsMCMC Results

#' Simulate Epidemic

N = c(200, 500, 1000)
beta0 = 0.1
gamma = 0.05
beta = beta0/N
a = 0.05*N
p = c(0.1, 0.2, 0.3, 0.4, 0.5)
k = c(2, 4, 6, 8, 10, 15, 20)


#' ==== Epidemic N = 200 ====

set.seed(1)
U = runif(2*N[1] - a[1])
E = rexp(2*N[1] - a[1])
epidemic200 = SIR_Gillespie(N[1], a[1], beta[1], gamma = gamma, obs_end = Inf, output = "not_panel",
                            U = U, E = E)
last_event_time = max(epidemic200$event_table$times)

#' Sample
i = 1
rep_sample = c(sample(1:a[1], p[i]*a[1]), sample((a[1] + 1):N[1], p[i]*(N[1]- a[1])))

j = 2
obs_times = seq(0, last_event_time, length = k[j])
panelData200 = with(epidemic200$event_table, panel_data_new.epidemics(ID, times, state, prev_state,
                                                                      subset = rep_sample, obs_times = obs_times))

s = 35
lambda = 0.01
no_its = 10000
burn_in = 0.1*no_its
run = fsMCMC.epidemics(type = "SIR", panelData200, obs_times, N[1], a[1], beta[1], gamma, kernel = NULL, no_draws = 2*N[1] - a[1],
                 s, lambda, V = diag(1, 2), no_its, burn_in)

