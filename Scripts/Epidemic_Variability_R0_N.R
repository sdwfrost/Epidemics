#' Variance of an Epidemic

#' Variance of Binomial Random Variable
#'
p = seq(0, 1, length = 1000)
plot(p, p*(1 - p), type = 'l', col = 'blue')

N = c(200, 500, 1000)
a = 0.05*N
R0 = c(1.2, 2, 5)
beta0 = 1
beta = beta0/N
gamma = beta0/R0

#' Plot epidmeics side by side
jpeg("R0highvsR0low.jpeg", width = 1200, height = 800)
par(mfrow = c(1,2))
epidemic1 = SIR_Gillespie(N[1], a[1], beta[1], gamma[1], U = runif(2*N[1] - a[1]), E = rexp(2*N[1] - a[1]))
with(epidemic1, plot.epidemic(event_table$times, N[1], X, Y, Z, newPlot = T))

epidemic2 = SIR_Gillespie(N[1], a[1], beta[1], gamma[3], U = runif(2*N[1] - a[1]), E = rexp(2*N[1] - a[1]))
with(epidemic2, plot.epidemic(event_table$times, N[1], X, Y, Z, newPlot = T))
dev.off()

jpeg("VarOfEpidemicIncreasingR0.jpeg", width = 1200, height = 800)
par(mfrow = c(1, 3))
#' R0 = 1.2

epidemic = SIR_Gillespie(N[1], a[1], beta[1], gamma[1], U = runif(2*N[1] - a[1]), E = rexp(2*N[1] - a[1]))
with(epidemic, plot.epidemic(event_table$times,N[1], X, Y, Z, newPlot = T))

# obs_times = seq(0, max(epidemic$event_table$times), length = 10)
#
# for(i in 1:length(obs_times)){
#   abline(v = obs_times[i], lty = 2, col = 'black')
# }

for(i in 1:50){
  epidemic = SIR_Gillespie(N[1], a[1], beta[1], gamma[2], U = runif(2*N[1] - a[1]), E = rexp(2*N[1] - a[1]))
  with(epidemic, plot.epidemic(event_table$times, N[1], X, Y, Z, newPlot = F))
}

#' R0 = 2

epidemic = SIR_Gillespie(N[1], a[1], beta[1], gamma[3], U = runif(2*N[1] - a[1]), E = rexp(2*N[1] - a[1]))
with(epidemic, plot.epidemic(event_table$times, N[1], X, Y, Z, newPlot = T))

# obs_times = seq(0, max(epidemic$event_table$times), length = 10)
#
# for(i in 1:length(obs_times)){
#   abline(v = obs_times[i], lty = 2, col = 'black')
# }

for(i in 1:50){
  epidemic = SIR_Gillespie(N[1], a[1], beta[1], gamma[2], U = runif(2*N[1] - a[1]), E = rexp(2*N[1] - a[1]))
  with(epidemic, plot.epidemic(event_table$times, N[1], X, Y, Z, newPlot = F))
}

#' R0 = 5

epidemic = SIR_Gillespie(N[1], a[1], beta[1], gamma[3], U = runif(2*N[1] - a[1]), E = rexp(2*N[1] - a[1]))
with(epidemic, plot.epidemic(event_table$times, N[1], X, Y, Z, newPlot = T))

# obs_times = seq(0, max(epidemic$event_table$times), length = 10)
#
# for(i in 1:length(obs_times)){
#   abline(v = obs_times[i], lty = 2, col = 'black')
# }

for(i in 1:50){
  epidemic = SIR_Gillespie(N[1], a[1], beta[1], gamma[3], U = runif(2*N[1] - a[1]), E = rexp(2*N[1] - a[1]))
  with(epidemic, plot.epidemic(event_table$times,N[1], X, Y, Z, newPlot = F))
}

dev.off()

jpeg("VarOfEpidemicIncreasingN.jpeg", width = 1200, height = 800)

par(mfrow = c(1, 3))
#' N = 200

epidemic = SIR_Gillespie(N[1], a[1], beta[1], gamma[2], U = runif(2*N[1] - a[1]), E = rexp(2*N[1] - a[1]))
with(epidemic, plot.epidemic(event_table$times, N[1], X, Y, Z, newPlot = T))

# obs_times = seq(0, max(epidemic$event_table$times), length = 10)
#
# for(i in 1:length(obs_times)){
#   abline(v = obs_times[i], lty = 2, col = 'black')
# }

for(i in 1:50){
  epidemic = SIR_Gillespie(N[1], a[1], beta[1], gamma[2], U = runif(2*N[1] - a[1]), E = rexp(2*N[1] - a[1]))
  with(epidemic, plot.epidemic(event_table$times,N[1], X, Y, Z, newPlot = F))
}

#' N = 500
epidemic = SIR_Gillespie(N[2], a[2], beta[2], gamma[2], U = runif(2*N[2] - a[2]), E = rexp(2*N[2] - a[2]))
with(epidemic, plot.epidemic(event_table$times, N[2], X, Y, Z, newPlot = T))

for(i in 1:50){
  epidemic = SIR_Gillespie(N[2], a[2], beta[2], gamma[2], U = runif(2*N[2] - a[2]), E = rexp(2*N[2] - a[2]))
  with(epidemic, plot.epidemic(event_table$times, N[2], X, Y, Z, newPlot = F))
}

#' N = 1000
epidemic = SIR_Gillespie(N[3], a[3], beta[3], gamma[2], U = runif(2*N[3] - a[3]), E = rexp(2*N[3] - a[3]))
with(epidemic, plot.epidemic(event_table$times, N[3], X, Y, Z, newPlot = T))

for(i in 1:50){
  epidemic = SIR_Gillespie(N[3], a[3], beta[3], gamma[2], U = runif(2*N[3] - a[3]), E = rexp(2*N[3] - a[3]))
  with(epidemic, plot.epidemic(event_table$times, N[3], X, Y, Z, newPlot = F))
}

dev.off()

#' Changing the epidemic


#' \pi(y | X) over X

epidemic = SIR_Gillespie(N, a, beta0/N, gamma[2], U = runif(2*N - a), E = rexp(2*N - a))

#' y
obs_times = seq(0, max(epidemic$event_table$times), length = 10)
panel_sample = c(sample(a, size = 0.1*a), sample((a+1):N, size = 0.1*(N - a)))
panel_data = with(epidemic$event_table, panel_data_new.epidemics(ID, times, state, prev_state,
                                                                      subset = panel_sample, obs_times = obs_times))


reject = 0
for(i in 1:1000){
  if(pi_Y_given_X(N, panel_data, obs_times, a, beta0/N, gamma[2]) == 0){
    reject = reject + 1
  }
  print(i)
}

plot(density(`piY|X`))

par(mfrow = c(1,1))
with(epidemic, plot.epidemic(event_table$times,N, X, Y, Z, newPlot = T))





