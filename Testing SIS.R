#'
#' Testing SIS Endemic
#'
#'
R_0 = 1.2 #' psi/gamma
gamma = 0.5
psi = R_0*gamma
N = 100
beta = psi/N
kernel = Contact_Kernel(Matrix::Matrix(1, nrow = N, ncol = N))
a  = 5
T_obs = c(0,30)
obs_end = T_obs[2]

no_events = c()
time_taken = c()
for(i in 1:10000){
  print(i)
  start_time = as.numeric(Sys.time())
  no_events[i] = SIS_Deterministic_Gillespie(N, a , beta, gamma,
                              E = rexp(20000), U = runif(20000), T_obs = c(0, 30),
                              k = 6, kernel, obs_end = 30)$no_events
  time_taken[i] = as.numeric(Sys.time()) - start_time

}

max(no_events)
mean(no_events)

max(time_taken)
mean(time_taken)
min(time_taken)

boxplot(time_taken, ylim = c(0, 0.25))

2*sum(time_taken)/60

# ==== SIS fsMCMC ====









