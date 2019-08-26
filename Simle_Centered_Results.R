#'
#' Simple Centered Results
#'
#'



#' ===== Homogeneous Epidemic, N = 200 =====
#'
#'


#' Simulate

N = 100
a = 0.05*N
R_0 = 2
gamma = 0.5
beta = gamma*R_0/N

kernel = Contact_Kernel(matrix(1, nrow = N, ncol = N))

#' Gives Final Size 144
i = 0
final_size = 0
while(final_size == 0){
  set.seed(i)
  sim = Epidemic_Simulation(N, a, gamma, beta, kernel)
  final_size = sim_high$final_size
  i = i + 1
}

set.seed(i - 1)
sim = Epidemic_Simulation(N, a, gamma, beta, kernel)

#' ===== Inference: 1 infection time changing =====

run = Centered_MCMC(N, t_rem = sim$t_rem, gamma, theta_gamma = c(1, 0.001), alpha = 1, beta, theta_beta = c(1, 0.001), kernel, no_proposals = 15,
                    no_its = 10000, burn_in = 1000)

#' Density Plots

par(mfrow = c(1,3))
plot(density(run$draws[,1]), type = 'l', col = "blue", xlab = expression(beta))
abline(v = beta, col = "red", lty = 2)

plot(density(run$draws[,2]), type = 'l', col = "blue", xlab = expression(gamma))
abline(v = gamma, col = "red", lty = 2)

plot(density(run$draws[,1]*N/run$draws[,2]), type= 'l', col = "blue", xlab = expression(R[0]))
abline(v = R_0, col = "red", lty = 2)



run2 = NonCentered_MCMC(N, t_rem = sim$t_rem, gamma, alpha = 1, theta_gamma = c(1, 0.001), beta, theta_beta = c(1,0.001),
                        kernel, no_its = 10000, burn_in = 1000, no_proposals = 1, lambda = 0.025)






