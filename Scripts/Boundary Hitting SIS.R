#'
#' Boundary Hitting in a Stochastic SIS model
#'
#'


#' Find SIS endimic which has a steady at x, y != 0 |  x, y != N
#'
#'

par(mfrow = c(1,2))
N = 100
psi = 0.02
gamma = 0.018
beta = psi/N
y_0 = N - gamma/beta
x_0 = N - y_0


Det_SIS = SIS_y_soln(N, x_0, y_0, beta, gamma, t_lim = 10, t_res = 0.0001)
plot.epidemic(times = Det_SIS$t, N, Det_SIS$x, Det_SIS$y)

final_time = NULL
for(i in 1:100){
  N = 1000
  y_0 = 0.05*N
  x_0 = N - y_0
  beta = psi/N
  Stoch_SIS = SIS_Gillespie(N, y_0, gamma, beta, obs_end = 2000)
  final_time[i] = tail(Stoch_SIS$event_table$time, n = 1)
  print(i)
}

extinction = sum(final_time < 2000)/100

final_time_200 = NULL

for(i in 1:100){
  N = 500
  y_0 = 0.05*N
  x_0 = N - y_0
  beta = psi/N
  Stoch_SIS = SIS_Gillespie(N, y_0, gamma, beta, obs_end = 2000)
  final_time_200[i] = tail(Stoch_SIS$event_table$time, n = 1)
  print(i)
}
extinction = sum(final_time_200 < 2000)/100


par(mfrow = c(2,2))

N = 200
psi = 0.02
gamma = 0.018
beta = psi/N
y_0 = 0.05*N
x_0 = N - y_0

Det_SIS = SIS_y_soln(N, x_0, y_0, beta, gamma, t_lim = 10, t_res = 0.0001)
plot.epidemic(times = Det_SIS$t, N, Det_SIS$x, Det_SIS$y)

Stoch_SIS = SIS_Gillespie(N, y_0, gamma, beta, obs_end = 2000)
plot.epidemic(Stoch_SIS$event_table$time, N, Stoch_SIS$X, Stoch_SIS$Y)

N = 1000
psi = 0.02
gamma = 0.018
beta = psi/N
y_0 = 0.05*N
x_0 = N - y_0

Det_SIS = SIS_y_soln(N, x_0, y_0, beta, gamma, t_lim = 10, t_res = 0.0001)
plot.epidemic(times = Det_SIS$t, N, Det_SIS$x, Det_SIS$y)

Stoch_SIS = SIS_Gillespie(N, y_0, gamma, beta, obs_end = 2000)
plot.epidemic(Stoch_SIS$event_table$time, N, Stoch_SIS$X, Stoch_SIS$Y)


