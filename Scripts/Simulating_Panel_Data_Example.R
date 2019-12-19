#'
#' SIS: Simulation, transition and panel data example
#'
#'

#' Run line below to install the package in its current iteration
devtools::install_github("JMacDonaldPhD/Epidemics", ref = "general_epidemics")

library(Epidemics)

# ==== Simulating Epidemic (SIR process) ====

N = 200
a = 5
gamma = 0.5
beta = 0.1
kernel = Contact_Kernel(matrix(1, N, N))
obs_end = Inf

set.seed(1)
epidemic = SIR_Gillespie(N, a, gamma, beta)

set.seed(1)
U = NULL
E = NULL
for(i in 1:(2*N - a)){
  E[i] = rexp(1)
  U[i] = runif(1)
}
epidemic.pre = SIR_Gillespie(N, a, gamma, beta, U = U, E = E)


set.seed(1)
epidemic_kernel = SIR_Gillespie(N, a, gamma, beta, kernel, obs_end)

#' Plot evolution of the epidemic
par(mfrow = c(1,2))
plot.epidemic(epidemic$event_table$time, N, epidemic$X, epidemic$Y, epidemic$Z)

plot.epidemic(epidemic.pre$event_table$time, N, epidemic.pre$X, epidemic.pre$Y,
              epidemic.pre$Z)

lines(c(0,epidemic$event_table$time[epidemic$event_table$time !=0]), epidemic$Y)

attach(epidemic$event_table)
transition_table.epidemics(ID, time, state, prev_state, period = c(0.0001, 51), 1:3, c("S", "I", "R"), output = "freq")$Freq

as.data.frame(panel_data.epidemics(ID, time, state, prev_state, state_names = c("S", "I", "R"),
                     T_obs = c(0,10), k = 10))

detach(epidemic$event_table)

# ==== Simulating Endemic (SIS process) ====

N = 1000
a = 5
gamma = 0.5
beta = 0.01
kernel = Contact_Kernel(matrix(1, N, N))
obs_end = 10

set.seed(1)
epidemic = SIS_Gillespie(N, a , gamma, beta, obs_end)

set.seed(1)
E = NULL
U = NULL

for(i in 1:10000){
  E[i] = rexp(1)
  U[i] = runif(1)
}
epidemic.pre = SIS_Gillespie(N, a, gamma, beta, obs_end , U = U, E = E)

set.seed(1)
epidemic_kernel = SIS_Gillespie(N, a , gamma, beta, obs_end, kernel)
#attach(epidemic$event_table)

par(mfrow = c(1,2))
plot.epidemic(epidemic$event_table$time, N, epidemic$X, epidemic$Y)

plot.epidemic(epidemic.pre$event_table$time, N, epidemic.pre$X, epidemic.pre$Y)

transition_table.epidemics(ID, time, state, prev_state, period = c(0, 2), 1:3, c("S", "I", "R"), output = "freq")

panel_data.epidemics(ID, time, state, prev_state, state_names = c("S", "I", "R"),
                     T_obs = c(0,10), k = 10)

detach(epidemic$event_table)

