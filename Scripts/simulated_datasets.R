#' Experiment Simulated Datasets

#' Simulate 3 density dependent SIR epidemics
#' Same R_0, gamma and beta_0
#' Vary the population size N
#' These datasets will be used to invesigate the scalability of the methods
#' and how methods deal with increasing missing data.

N = c(200, 500, 1000)
a = 0.05*N
gamma = 0.5
beta = 1/N

#' N  = 200

kernel200 = Contact_Kernel(matrix(1, nrow = N[1], ncol = N[1]))
FINALSIZE = 0
while(FINALSIZE == 0){
  simdata200 = Epidemic_Simulation(N[1], a[1], gamma, beta[1], kernel200)
  FINALSIZE = simdata200$final_size
}

#' N = 500

kernel500 = Contact_Kernel(matrix(1, nrow = N[2], ncol = N[2]))
FINALSIZE = 0
while(FINALSIZE == 0){
  simdata500 = Epidemic_Simulation(N[2], a[2], gamma, beta[2], kernel500)
  FINALSIZE = simdata500$final_size
}


#' N = 1000

kernel1000 = Contact_Kernel(matrix(1, nrow = N[3], ncol = N[3]))
FINALSIZE = 0
while(FINALSIZE == 0){
  simdata1000 = Epidemic_Simulation(N[3], a[3], gamma, beta[3], kernel1000)
  FINALSIZE = simdata1000$final_size
}

usethis::use_data(simdata200, simdata500, simdata1000, overwrite  = T)







