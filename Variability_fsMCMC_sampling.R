#'
#' Variablility of Inference
#' Random Samples
#'
#'
#'
#'

N = test_epidemic$N
beta0 = test_epidemic$par_beta
gamma0 = test_epidemic$par_gamma[1]
initial_infective = 5
t_inf = test_epidemic$t_inf
t_rem = test_epidemic$t_rem

T_obs = c(0,10)
k = 6



V = diag(1, 2)
lambda = 0.04
s = 50
no_its = 10000
burn_in = 1000

#' Storage

Variability = matrix(NA, nrow = 1000, ncol = 9)
x_store = list()

for(i in 11:20){
  random_sample = sample(1:N, 0.1*N)
  x = transitions_between_observations(T_obs, k, event_times = cbind(t_inf[random_sample], t_rem[random_sample]))
  Test = Epidemic_fsMCMC(N, initial_infective, x, beta0, gamma0, no_draws = 2*N, s, T_obs, k, lambda, V,
                         no_its, burn_in)

  Variability[i,] = c(Test$beta_summary, Test$gamma_summary, Test$R0_summary, min(Test$ESS), Test$accept_rate)
  x_store[[i]] = x
}


par(mfrow = c(1,2))

boxplot(Variability[,1])
boxplot(Variability[,3])


# ==== Big Epidemic ====

N = big_epidemic$N
beta0 = big_epidemic$par_beta
gamma0 = big_epidemic$par_gamma[1]
initial_infective = 1
t_inf = big_epidemic$t_inf
t_rem = big_epidemic$t_rem

T_obs = c(0,60)
k = 60




V = diag(1, 2)
lambda = 0.0005
s = 75
no_its = 10000
burn_in = 1000

random_sample = sample(1:N, 0.1*N)
x = transitions_between_observations(T_obs, k, event_times = cbind(t_inf[random_sample], t_rem[random_sample]))
Test = Epidemic_fsMCMC(N, initial_infective, x, beta0, gamma0, no_draws = 2*N, s, T_obs, k, lambda, V,
                       no_its, burn_in)


#' Storage

Variability = matrix(NA, nrow = 1000, ncol = 9)
x_store = list()
i = 21
for(i in 11:20){
  random_sample = sample(1:N, 0.1*N)
  x = transitions_between_observations(T_obs, k, event_times = cbind(t_inf[random_sample], t_rem[random_sample]))
  Test = Epidemic_fsMCMC(N, initial_infective, x, beta0, gamma0, no_draws = 2*N, s, T_obs, k, lambda, V,
                         no_its, burn_in)

  Variability[i,] = c(Test$beta_summary, Test$gamma_summary, Test$R0_summary, min(Test$ESS), Test$accept_rate)
  x_store[[i]] = x
}


par(mfrow = c(1,2))
boxplot(Variability[,1])
boxplot(Variability[,3])




