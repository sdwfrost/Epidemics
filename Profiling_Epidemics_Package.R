#'
#'
#'
tmp = tempfile()
Rprof(tmp, interval = 0.01)
set.seed(3)
N = 200
a = 0.05*N
beta = 0.02
gamma = 0.5

obs_times = seq(0, 10, length = 10)

U = runif(2*N - a)
E = rexp(2*N - a)
sim = Infectious_Disease_Gillespie(type = "SIR", N, a, beta, gamma, obs_end = tail(obs_times, n = 1), output = "event",
                                   U = U, E = E)
attach(sim)
par(mfrow  = c(1,1))
plot.epidemic(event_table$times, N, X, Y, Z)
detach(sim)
x = with(sim$event_table, panel_data_new.epidemics(ID, times, state, prev_state,
                                                   subset = c(1:19, 191), obs_times = obs_times))

s = 2*N - a
lambda = 3
V = diag(1, 2)

run = fsMCMC.epidemics(type = "SIR", x, obs_times, N, a, beta, gamma, kernel = NULL, no_draws = 2*N - a, s,
                       lambda, V, no_its = 10000)
Rprof(NULL)
summaryRprof(tmp)$by.self
summaryRprof(tmp)$sampling.time


tmp = tempfile()
Rprof(tmp, interval = 0.01)
X = 95
Y = 5
individual_inf_rate = rep(0.01*Y, X)

gamma = 0.5
for(i in 1:1e6){
  #event_epidemicsC(individual_inf_rate, gamma, Y, U = NaN)
  event.epidemics(individual_inf_rate, gamma, Y)
}
#bench::mark(event_epidemicsC(individual_inf_rate, gamma, Y, U),
#            event.epidemics(individual_inf_rate, gamma, Y, U))

Rprof(NULL)
summaryRprof(tmp)$by.self
summaryRprof(tmp)$sampling.time


individual_inf_rate = rep(beta, 190)
for(j in 1:no_reps){
  if(is.null(U[j])){
    U[j] = runif(1)
  }
  X = length(individual_inf_rate)
  rates = c(individual_inf_rate, rep(gamma,Y))
  total_rate = sum(rates)
  ID_index = sum(cumsum(rates)/total_rate < U[j]) + 1

  if(ID_index <= X){
    event = 0
  } else{
    event = 1
    ID_index = ID_index - X
  }
}

microbenchmark::microbenchmark(scan(text = "1 2 3", quiet = T))

beta  = 0.01
individual_inf_rate = rep(beta, 190)
Y = 10
gamma = 0.5

microbenchmark::microbenchmark(event.epidemics(individual_inf_rate, gamma, Y), event.epidemics2(),
                               times = 10000)

tmp = tempfile()
Rprof(tmp, interval = 0.001)

no_reps = 1:10000000
U = NULL
Y = 10
beta = 0.02
individual_inf_rate = rep(beta, 190)
gamma = 0.5

for(i in no_reps){
  #event.epidemics(individual_inf_rate, gamma, Y)
  #event.epidemics2()
  Y = 10
  beta = 0.02
  individual_inf_rate = rep(beta, 190)
  gamma = 0.5
}

Rprof(NULL)
summaryRprof(tmp)$by.self
summaryRprof(tmp)$sampling.time
S = 1:95
I = 96:100
R = NULL
for(i in 1:50000){
  prev_state = c(rep(1, 3), rep(2, 3), rep(3, 3), rep(c(1,2,3,1), 50))
  state = c(rep(c(1, 2, 3), 3), rep(c(2,3,3,1), 50))
  trans = data.frame(prev_state, state)

  #state_table(x = c(rep(1, 3), rep(2, 3), rep(3, 3), rep(c(1,2,3,1), 50)),
  #            y = c(rep(c(1, 2, 3), 3), rep(c(2,3,3,1), 50)),
  #            levels = 1:3)
  print(i)
}



library(dplyr)
freq = dplyr::count(mtcars, c(`cyl`, `gear`))

mtcars
state_table(x, y, levels = 1:3)


before = table(x[-(1:9)], y[-(1:9)])
after = table(x, y) - 1
before
after
obs_times = seq(T_obs[1], T_obs[2], length = k)

rep_sample = c(sample(N[1] - a[1], size = 0.1*(N[1] - a[1])),
               sample((N[1] - a[1] + 1 ):N[1], size = 0.1*a[1]))

x_sim = with(obs_sim$event_table, panel_data.epidemics(ID, time, state, prev_state, state_names = c("S", "I", "R"),
                                           levels = 1:3, rep_sample, obs_times = obs_times))
kernel = Contact_Kernel(matrix(1, nrow = N, ncol = N))
run = fsMCMC.epidemics(x_sim, obs_times, N[1], a[1], beta[1], gamma, kernel,
                       no_draws = 2*N[1] - a[1],
                       s = 180, lambda = 0.03, V = diag(1, 2), no_its = 100)

U = runif(2*N[1] - a[1])
E = rexp(2*N[1] - a[1])


max_events = 2*N[1] - a[1]
x = 1:max_events

end_state(subject, state)

state[sapply(X = ID_events, function(X) head(which(subject == X), n = 1))]
start_state(subject, state)

head(which(subject == 14), n = 1)
sum(subject == 14)
ID_events = 1:30
subject = rep(1:30, 4)
relevant_events = 1:length(subject)
state = rep(1, length(subject))
state[sapply(X = ID_events, function(X) head(which(subject == X), n = 1))]

tmp = tempfile()
Rprof(tmp, interval = 0.01)
for(i in 1:10000){
  Y = with(obs_sim$event_table, panel_data.epidemics(ID, time, state, prev_state,
                                                     state_names = c("S", "I", "R"),
                                                     levels = 1:3, T_obs = T_obs, k = k))
  print(i)
}
Rprof(NULL)
summaryRprof(tmp)

sum(Rprof.summary$by.total$self.time)
detach(obs_sim$event_table)

# Finding the last element of a subset
tail(which(subject[previous_events] == X), n = 1))

which(1 == x)
x[x %in% 1]

x[which(x == 1)]


ID_event = 1:3
subject = c(1, 2, 3, 2, 1, 2, 3)
state = c(1, 1, 1, 2, 2, 3, 2)
state[sapply(X = ID_event, function(X) tail(which(subject == X), n = 1))]

end_state(subject, state)



SIR_Gillespie(N[1], a[1], beta[1], gamma, U = U, E = E)
