#' Effect of epidemic variability on Pseudo-Marginal MCMC

#N = seq(200, 1000, by = 80)
N = c(200, 500, 1000)
a = 0.05*N
R_0 = seq(1.2, 5, length = 10)
s = seq(1, 50, by = 5)
beta0 = 1
beta = beta0/N
gamma = beta0/R_0

#' For every N
#' Simulate an Epidemic with theta
#' sample panel data y
#' Observe 10% of the population
#' obs_times = seq(0, max(epidemic$event_table$times), k = 10)
#' simulate s epidemics calculate pi_y_X. sum(pi_y_X > 0)

i = 1
j = 1
epidemic = SIR_Gillespie(N[i], a[i], beta[i], gamma[j], U = runif(2*N[i] - a[i]), E = rexp(2*N[i] - a[i]))

panelSample = sample(1:N[i], size = 0.1*N[i])
obsTimes = seq(0, max(epidemic$event_table$times), length = 20)
panelData = with(epidemic$event_table, panel_data_new.epidemics(ID, times, state, prev_state,
                                                                subset = panelSample, obs_times = obsTimes))
noIts = 100

goodItsperSec = c()
for(k in 1:length(s)){
  start = as.numeric(Sys.time())
  goodIts = sum(sapply(1:noIts, function(X) pi_Y_given_theta_hat(no_sims = s[k], N[i], panelData, obsTimes, a[i], beta[i], gamma[j]) > 0))
  timeTaken = as.numeric(Sys.time()) - start
  print(timeTaken)
  goodItsperSec[k] = goodIts/timeTaken
  print(k)
}

par(mfrow = c(1,1))
plot(s, goodItsperSec)
#sum(sapply(1, function(s) pi_Y_given_X(N[i], panelData, obsTimes, a[i], beta[i], gamma[j])) > 0)

noIts = 1000
vartime = c()
for(k in 1:5){
  start = as.numeric(Sys.time())
  varRun = var(sapply(1:noIts, function(X) pi_Y_given_theta_hat(no_sims = s[k], N[i], panelData, obsTimes, a[i], beta[i], gamma[j])))
  timeTaken = as.numeric(Sys.time()) - start
  print(timeTaken)
  vartime[k] = varRun*timeTaken
  print(k)
}

plot(s[1:5], vartime)



pi_Y_given_X(N[i], panelData, obsTimes, a[i], beta[i], gamma[j])


