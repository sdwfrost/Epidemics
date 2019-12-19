#' Importance Sampler Performance (Homogeneously Mixing SIR)

m = with(simdata200, 0.1*N)
panelSample = readRDS(file = "test_panelSample.rds")
eventTimesSample = with(simdata200, cbind(t_inf, t_rem)[panelSample,])

final_event_time = with(simdata200, max(t_rem[t_rem < Inf]))
k = 6
T_obs = c(0, final_event_time)
#' Increasing number panels
# for(k in seq(2, 12, by = 1)){
#   obs_times = seq(T_obs[1], T_obs[2], length = k)
#   paneldata200 = with(simdata200, transitions_between_observations(T_obs, k, event_times = eventTimesSample))
#   ISrun = with(simdata200,
#                epidemicImportanceSampler(N, a = N*0.05, paneldata200[1], obs_times[1:2], theta = c(par_beta, par_gamma),
#                                          N_p = 100))
# }


ESS = lapply(X = 2:12, function(X){
  print(X)
  obs_times = seq(T_obs[1], T_obs[2], length = X)
  paneldata200 = with(simdata200, transitions_between_observations(T_obs, X, event_times = eventTimesSample))
  ISrun = with(simdata200,
               epidemicImportanceSampler(N, a = N*0.05, paneldata200, obs_times, theta = c(par_beta, par_gamma),
                                         N_p = 1000))$ESS
})



save.image(file = "IS_ESS_Increasing_panel_no.rds")

plot(2:12, ESS)
