
# =========================================== Testing Centered MCMC

N = 200
beta0 = 1
beta = beta0/N
gamma = 0.5
a = 0.05*N
kernel = Contact_Kernel(matrix(1, nrow = N, ncol = N))

final_size = 0
while(final_size == 0){
  epidemic = Epidemic_Simulation(N, a, gamma, beta, kernel)
  final_size = epidemic$final_size
}

run = Centered_MCMC(N, a, epidemic$t_rem, gamma, c(1, 0.001), beta, c(1, 0.001), kernel, no_proposals = 5, no_its = 10000,
              burn_in = 1000)

# =========================================== Testing Non-centered MCMC



