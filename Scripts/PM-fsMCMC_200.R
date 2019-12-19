#' Pseudo-Marginal/fsMCMC simdata200

#' Set m = 0.1*N

m = with(simdata200, 0.1*N)
panelSample = readRDS(file = "test_panelSample.rds")
eventTimesSample = with(simdata200, cbind(t_inf, t_rem)[panelSample,])

#' k = 6 works well with current sample, slight underestimation of gamma
#' Seems to converge at about 10000 iterations
#' k = 12
#' Works less well
#'
#' k = 4
#' Seems to work less well as well
#' lambda = 0.02, s = 120
#'

#' Observe second half of epidemic only?
#' k = 6
#' Works really well
#' s = 90, lambda = 0.012
#'
#' Observe first half, over estimates drastically
#' k = 6
#'
#'

observe_param = list()
k_list = c(6, 12, 4, 6, 6)
final_event_time = with(simdata200, max(t_rem[t_rem < Inf]))
T_obs_list = rbind(c(0, final_event_time), c(0, final_event_time), c(0, final_event_time),
                   c(0.5*final_event_time, final_event_time), c(0, 0.5*final_event_time))
lambda_list = c(0.011, 0.009, 0.02, 0.012 , 75)
s_list = c(75, 60, 120, 90, 75 )
for(i in 1:length(k_list)){
  observe_param[[i]] = list(k = k_list[i], T_obs = T_obs_list[i,],
                            obs_times = seq(T_obs_list[i, 1], T_obs_list[i, 2], length = k_list[i]),
                            paneldata200 = with(simdata200, transitions_between_observations(T_obs_list[i,], k_list[i], event_times = eventTimesSample)),
                            s = s_list[i], lambda = lambda_list[i])
}

no_its = 50000
burn_in = 0
RUNS = lapply(X = 1:5, function(X) with(observe_param[[X]],
                                 with(simdata200, fsMCMC.epidemics(type = "SIR", paneldata200, obs_times = obs_times,
                                                                   N, 0.05*N, no_draws = 2*N - 0.05*N, par_beta, par_gamma,
                                                                   s = 75, lambda = 0.012, V = diag(1, 2),
                                                                   no_its = no_its, burn_in = burn_in))))


par(mfrow = c(1, 2))
for(j in 1:2){
  if(j == 1){
    ylim = c(0,1000)
    v = simdata200$par_beta
  } else{
    ylim = c(0, 10)
    v = simdata200$par_gamma
  }
  with(RUNS[[1]], plot(density(draws[,j]), ylim = ylim, type = 'l'))
  for(i in 1:length(RUNS)){
    with(RUNS[[i]], lines(density(draws[,j]), col = i))
  }
  abline(v = v, col = 'blue', lty = 2)
}



save.image(file = "fsMCMC200Results.rds")
#
# load("fsMCMC200Results.rds")


#'
#' k = 6
#' #T_obs = c(0, with(simdata200, max(t_rem[t_rem < Inf])))
#'
#' T_obs = c(0, 0.5*with(simdata200, max(t_rem[t_rem < Inf])))
#' #'T_obs = c(0.5*with(simdata200, max(t_rem[t_rem < Inf])), with(simdata200, max(t_rem[t_rem < Inf])))
#'
#' obs_times = seq(T_obs[1], T_obs[2], length = k)
#'
#'
#' paneldata200 = with(simdata200, transitions_between_observations(T_obs, k, event_times = eventTimesSample))
#'
#' # = fsMCMC
#'
#' no_its = 10000
#' burn_in = 0
#'
#' opt_run = with(simdata200, fsMCMC.epidemics(type = "SIR", paneldata200, obs_times = obs_times,
#'                                             N, 0.05*N, no_draws = 2*N - 0.05*N, par_beta, par_gamma,
#'                                             s = 75, lambda = 0.012, V = diag(1, 2),
#'                                             no_its = no_its, burn_in = burn_in))
#'
#' par(mfrow = c(1, 2))
#' with(opt_run, plot(density(draws[,1])))
#' with(simdata200, abline(v = par_beta, col = 'blue', lty = 2))
#' with(opt_run, plot(density(draws[,2])))
#' with(simdata200, abline(v = par_gamma, col = 'blue', lty = 2))
#'
#' RUNS = lapply(X = 1:2, function(X) with(simdata200, fsMCMC.epidemics(type = "SIR", paneldata200, obs_times = obs_times,
#'                                                                       N, 0.05*N, no_draws = 2*N - 0.05*N, par_beta, par_gamma,
#'                                                                       s = 90, lambda = 0.012, V = diag(1, 2),
#'                                                                       no_its = no_its, burn_in = burn_in)))
#'



#'
#'
#'
#' #' Pseudo-Marginal MCMC
#' #'
#' k = 6
#' T_obs = c(0, with(simdata200, max(t_rem[t_rem < Inf])))
#'
#' #'T_obs = c(0, 0.5*with(simdata200, max(t_rem[t_rem < Inf])))
#' #'T_obs = c(0.5*with(simdata200, max(t_rem[t_rem < Inf])), with(simdata200, max(t_rem[t_rem < Inf])))
#' obs_times = seq(T_obs[1], T_obs[2], length = k)
#' paneldata200 = with(simdata200, transitions_between_observations(T_obs, k, event_times = eventTimesSample))
#'
#'
#' no_its = 1000
#' burn_in = 0
#'
#' opt_run = with(simdata200, Pseudo_Marginal_MCMC(N, paneldata200, obs_times = seq(0, T_obs[2], length = k),
#'                                       par_beta, par_gamma, prior_rate = 0.001, initial_infective = 0.05*N,
#'                                       lambda = 0.1, V = diag(1, 2), no_sims = 50,
#'                                       no_its = no_its, burn_in = burn_in, parallel = T, no_cores = 4))
#'
#' par(mfrow = c(1, 2))
#' with(opt_run, plot(density(draws[,1])))
#' with(simdata200, abline(v = par_beta, col = 'blue', lty = 2))
#' with(opt_run, plot(density(draws[,2])))
#' with(simdata200, abline(v = par_gamma, col = 'blue', lty = 2))
#'
#'
#'
#' #' Check Convergence
#'
#' no_its = 50000
#' burn_in = 0
#' convCheckRun = with(simdata200, fsMCMC.epidemics(type = "SIR", paneldata200, obs_times = seq(0, T_obs[2], length = k),
#'                                             N, 0.05*N, no_draws = 2*N - 0.05*N, 0.1, 0.5,
#'                                             s = 70, lambda = 0.01, V = diag(1, 2),
#'                                             no_its = no_its, burn_in = burn_in))
#'
#' with(opt_run, plot(density(draws[,1])))
#'
#'
#'
#'

