#'
#' Profile Likelihoods of Toy Epidemic Example
#'



N = toy_epidemic$N
t_inf = toy_epidemic$t_inf
t_rem = toy_epidemic$t_rem
par_gamma = toy_epidemic$par_gamma
par_beta = toy_epidemic$par_beta
kernel = Contact_Kernel(matrix(1, nrow = N, ncol = N) - diag(1, N))

no_evals = 1000

individual_inf_time = seq(0, 50, length = no_evals)
epidemic_loglikelihood(par_rem = par_gamma,
                       par_inf = par_beta, t_inf_new,
                       t_rem, kernel)
infected_individuals  =  which(t_inf < Inf)
profile_log_likelihood_inf_matrix = matrix(NA, nrow = length(infected_individuals),
                                           ncol = no_evals)
for(j in 1:length(infected_individuals)){
  for(i in 1:no_evals){
    t_inf_new = t_inf
    t_inf_new[infected_individuals[j]] = individual_inf_time[i]
    profile_log_likelihood_inf_matrix[j,i] = epidemic_loglikelihood(par_rem = par_gamma,
                                                           par_inf = par_beta, t_inf_new,
                                                           t_rem, kernel)
  }
}

plot(individual_inf_time, exp(profile_log_likelihood_inf_matrix[3,]), type = 'l',
     xlim = c(0, 10), ylim = c(0, 1*10^-25))

for(i in 2:6){
  lines(individual_inf_time, exp(profile_log_likelihood_inf_matrix[i,]), type = 'l',
        col = i)
}

# ================ Removal Times ========================


no_evals = 1000

individual_rem_time = seq(0, 50, length = no_evals)

removed_individuals  =  which(t_rem < Inf)
profile_log_likelihood_rem_matrix = matrix(NA, nrow = length(removed_individuals),
                                           ncol = no_evals)

for(j in 1:length(removed_individuals)){
  for(i in 1:no_evals){
    t_rem_new = t_rem
    t_rem_new[removed_individuals[j]] = individual_rem_time[i]
    profile_log_likelihood_rem_matrix[j,i] = epidemic_loglikelihood(par_rem = par_gamma,
                                                                    par_inf = par_beta, t_inf,
                                                                    t_rem_new, kernel)
  }
}

plot(individual_rem_time, exp(profile_log_likelihood_rem_matrix[6,]), type = 'l')

for(i in 2:6){
  lines(individual_inf_time, exp(profile_log_likelihood_inf_matrix[i,]), type = 'l',
        col = i)
}

#================= Profile Beta ================


beta_seq = seq(0, 0.1, length = no_evals)
for(i in 1:no_evals){
  profile_log_likelihood_beta[i] = epidemic_loglikelihood(par_rem = par_gamma,
                                                          par_inf = beta_seq[i], t_inf,
                                                          t_rem, kernel)
}



plot(beta_seq, exp(profile_log_likelihood_beta), type = 'l')

exp(epidemic_loglikelihood(par_rem = par_gamma, par_inf = par_beta, t_inf, t_rem, kernel))















