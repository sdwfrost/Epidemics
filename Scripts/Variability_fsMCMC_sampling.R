#'
#' Variablility of Inference
#' Random Samples
#'
#'
#'
#'

psi = 1 # Normalised infection rate
p = 0.05 # Proportion Infected
N = c(2*10^2, 10^3, 10^4, 10^5) # Population
a = p*N
beta = psi/N # Infection parameter
gamma = 0.15 # Removal Parameter
prop_obs = 0.1
T_obs = c(0,60)
k = 12
#' Storage
Variability = matrix(NA, nrow = 400, ncol = 11)
sample_type = c("random", "representative")

V = diag(1,2)
no_its = 10000
burn_in = 1000

i = 1
lambda = 0.03
s = 180
for(j in 1:100){
  rep_sample = c(sample(N[i] - a[i], size = prop_obs*(N[i] - a[i])),
                 sample((N[i] - a[i] + 1 ):N[i], size = prop_obs*a[i]))
  x_rep = gillespie_panel_transitions(X_0 = prop_obs*(N[i] - a[i]), Y_0 = prop_obs*a[i],
                                      Epidemic1[,5], Epidemic1[,1], Epidemic1[,6],
                                      T_obs, k, subset = rep_sample)
  Test = Epidemic_fsMCMC(N[i], a[i], x_rep, beta[i], gamma, no_draws = 2*N[i], s, T_obs, k, lambda, V,
                         no_its, burn_in)
  Variability[j,] = c(sample_type[2], N[i], Test$beta_summary, Test$gamma_summary, Test$R0_summary, min(Test$ESS), Test$accept_rate)
  print(j)
}

for(j in 101:200){
  random_sample = sample(1:N[i], size = prop_obs*N[i])
  x_random = gillespie_panel_transitions(X_0 = sum(random_sample <= N[i] - a[i]),
                                         Y_0 = sum(random_sample > N[i] - a[i]),
                                         Epidemic1[,5], Epidemic1[,1], Epidemic1[,6],
                                         T_obs, k, subset = random_sample)
  Test = Epidemic_fsMCMC(N[i], a[i], x_random, beta[i], gamma, no_draws = 2*N[i], s, T_obs, k, lambda, V,
                         no_its, burn_in)
  Variability[j,] = c(sample_type[1], N[i], Test$beta_summary, Test$gamma_summary, Test$R0_summary, min(Test$ESS), Test$accept_rate)
  print(j)
}

i = 2
lambda = 0.005
s = 550

for(j in 201:300){
  rep_sample = c(sample(N[i] - a[i], size = prop_obs*(N[i] - a[i])),
                 sample((N[i] - a[i] + 1 ):N[i], size = prop_obs*a[i]))
  x_rep = gillespie_panel_transitions(X_0 = prop_obs*(N[i] - a[i]), Y_0 = prop_obs*a[i],
                                      Epidemic2[,5], Epidemic2[,1], Epidemic2[,6],
                                      T_obs, k, subset = rep_sample)
  Test = Epidemic_fsMCMC(N[i], a[i], x_rep, beta[i], gamma, no_draws = 2*N[i], s, T_obs, k, lambda, V,
                         no_its, burn_in)
  Variability[j,] = c(sample_type[2], N[i], Test$beta_summary, Test$gamma_summary, Test$R0_summary, min(Test$ESS), Test$accept_rate)
  print(j)
}

for(j in 301:400){
  random_sample = sample(1:N[i], size = prop_obs*N[i])
  x_random = gillespie_panel_transitions(X_0 = sum(random_sample <= N[i] - a[i]),
                                         Y_0 = sum(random_sample > N[i] - a[i]),
                                         Epidemic2[,5], Epidemic2[,1], Epidemic2[,6],
                                         T_obs, k, subset = random_sample)
  Test = Epidemic_fsMCMC(N[i], a[i], x_random, beta[i], gamma, no_draws = 2*N[i], s, T_obs, k, lambda, V,
                         no_its, burn_in)
  Variability[j,] = c(sample_type[1], N[i], Test$beta_summary, Test$gamma_summary, Test$R0_summary, min(Test$ESS_sec), Test$accept_rate)
  print(j)
}

usethis::use_data(Variability, overwrite = T)




