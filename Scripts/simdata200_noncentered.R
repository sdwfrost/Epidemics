#' simdata200_noncentered

#' Centered and Non-Centered MCMC on Simulated datasets

flatPrior = c(1, 0.001)
no_its = 10000
burn_in = 0
# ==== Centered ====

#' N = 200

#' Tune
opt_run = with(simdata200, NonCentered_MCMC(N, 0.05*N, t_rem, par_gamma, theta_gamma = flatPrior, par_beta,
                                         theta_beta = flatPrior, no_proposals = 5, lambda = 0.02,
                                         no_its = no_its, burn_in = burn_in))

#' Optimal:
#' no_proposals = 12
#' lambda = 0.02

#' Density
par(mfrow = c(1, 2))
with(opt_run, plot(density(draws[-(1:1000), 1]*200)))
with(opt_run, plot(density(draws[-(1:1000), 2])))

#' Check Convergence
# Choose a range of values for beta and gamma
# Do runs with these values, and so if they converge to the same distribution
no_its = 1000
burn_in = 0
beta_values = runif(10, 0, 0.1)
gamma_values = runif(10, 0, 5)

convCheckRuns = lapply(X = 1:10, FUN =  function(X) with(simdata200, NonCentered_MCMC(N, 0.05*N, t_rem, gamma_values[X], theta_gamma = flatPrior, beta_values[X],
                                                                                   theta_beta = flatPrior, no_proposals = 12, lambda = 0.02,
                                                                                   no_its = no_its, burn_in = burn_in)))

with(convCheckRuns[[1]], plot(draws[,1], type = 'l'))

with(convCheckRuns[[2]], lines(draws[,1], type = 'l'))
#' Trace Plots showing they all converge

par(mfrow = c(1, 2))
for(j in 1:2){
  if(j == 1){
    ylim = c(0, 0.1)
  } else{
    ylim = c(0, 10)
  }
  with(convCheckRuns[[1]], plot(draws[1:1000,j], ylim = ylim, type = 'l'))

  for(i in 1:length(convCheckRuns)){
    with(convCheckRuns[[i]], lines(draws[1:1000,j], col = i))
  }
}

#' Convergence Issues
#' Do these convergence issues subside as N grows?


#' Final Long Run
no_its = 50000
burn_in = 5000

finalRun = with(simdata200,  NonCentered_MCMC(N, 0.05*N, t_rem, par_gamma, theta_gamma = flatPrior, par_beta,
                                           theta_beta = flatPrior, no_proposals = 5, lambda = 0.02,
                                           no_its = no_its, burn_in = burn_in))

#' Density Plots

par(mfrow = c(1, 2))

with(finalRun, plot(density(draws[-(1:20000), 1]), col = 'blue', xlab = expression(beta[0])))
with(simdata200, abline(v = par_beta, col = 'red', lty = 2))

with(finalRun, plot(density(draws[-(1:20000), 2]), col = 'blue', xlab = expression(gamma)))
with(simdata200, abline(v = par_gamma, col = 'red', lty = 2))

with(finalRun, plot(density(200*draws[-(1:20000), 1]/draws[-(1:20000), 2]), col = 'blue', xlab = expression(gamma)))
with(simdata200, abline(v = N*par_beta/par_gamma, col = 'red', lty = 2))



#' Comments
#' Beta estimated well
#' Gamma Underestimated Slightly

save.image(file = "simdata200_noncentered.RData")

