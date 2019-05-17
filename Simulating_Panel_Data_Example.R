#'
#' SIS: Simulation, transition and panel data example
#'
#'

#' Run line below to install the package in its current iteration
devtools::install_github("JMacDonaldPhD/Epidemics", ref = "general_epidemics")

library(Epidemics)

# ==== Simulating Epidemic (SIR process) ====

N = 100
a = 5
gamma = 0.1
beta = 0.001
kernel = Contact_Kernel(matrix(1, N, N))
obs_end = Inf
epidemic = SIR_Gillespie(N, a, gamma, beta, kernel, obs_end)

attach(epidemic$event_table)

#' Plot evolution of the epidemic
plot(time, epidemic$Y, type = 'l', col = 'red', ylim = c(0,100))
lines(time, epidemic$X, type = 'l', col = 'blue')

transition_table.epidemics(ID, time, state, prev_state, period = c(0.0001, 51), 1:3, c("S", "I", "R"))

panel_data.epidemics(ID, time, state, prev_state, state_names = c("S", "I", "R"),
                     T_obs = c(0,10), k = 10)

detach(epidemic$event_table)

# ==== Simulating Endemic (SIS process) ====
set.seed(1)
N = 10
a = 1
gamma = 0.1
beta = 0.001
kernel = Contact_Kernel(matrix(1, N, N))
obs_end = 50

epidemic = SIS_Gillespie(N, a , gamma, beta, kernel, obs_end)
attach(epidemic$event_table)

#' Plot evolution of the epidemic
plot(time, epidemic$Y, type = 'l', col = 'red', ylim = c(0,100))
lines(time, epidemic$X, type = 'l', col = 'blue')

transition_table.epidemics(ID, time, state, prev_state, period = c(0, 51), 1:2, c("S", "I"))

panel_data.epidemics(ID, time, state, prev_state, state_names = c("S", "I"),
                     T_obs = c(0,51), k = 3)

detach(epidemic$event_table)

