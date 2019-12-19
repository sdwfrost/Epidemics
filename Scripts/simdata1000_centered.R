#' simdata1000_centered

#' simdata1000_Centered

#' Centered and Non-Centered MCMC on Simulated datasets

flatPrior = c(1, 0.001)
no_its = 10000
burn_in = 0
# ==== Centered ====

#' N = 200


#' Tune
opt_run = with(simdata1000, Centered_MCMC(N, 0.05*N, t_rem, par_gamma, theta_gamma = flatPrior, par_beta,
                                         theta_beta = flatPrior, no_proposals = 10,
                                         no_its = no_its, burn_in = burn_in))



#' Optimal:
#' no_proposals = 7

#' Density
par(mfrow = c(1, 2))
with(opt_run, plot(contour(draws[-(1:1000), 1]*500, draws[-(1:1000), 2])))
with(opt_run, plot(density(draws[-(1:1000), 2])))


#' Check Convergence
# Choose a range of values for beta and gamma
# Do runs with these values, and so if they converge to the same distribution
no_its = 1000
burn_in = 0
beta_values = runif(10, 0, 0.1)
gamma_values = runif(10, 0, 5)

convCheckRuns = lapply(X = 1:10, FUN =  function(X) with(simdata1000, Centered_MCMC(N, 0.05*N, t_rem, gamma_values[i], theta_gamma = flatPrior, beta_values[i],
                                                                                   theta_beta = flatPrior, kernel, no_proposals = 7,
                                                                                   no_its, burn_in)))

#' Trace Plots showing they all converge

par(mfrow = c(1, 2))
for(j in 1:2){
  if(j == 1){
    ylim = c(0, 0.1)
  } else{
    ylim = c(0, 10)
  }
  with(convCheckRuns[[1]], plot(draws[,j], ylim = ylim, type = 'l'))

  for(i in 1:length(convCheckRuns)){
    with(convCheckRuns[[1]], lines(draws[,j], col = i))
  }
}




#' Final Long Run
no_its = 50000
burn_in = 5000

finalRun = with(simdata1000,  Centered_MCMC(N, 0.05*N, t_rem, par_gamma, theta_gamma = flatPrior, par_beta,
                                           theta_beta = flatPrior, kernel, no_proposals = 7,
                                           no_its, burn_in))

#' Density Plots

par(mfrow = c(1, 2))

with(finalRun, plot(density(draws[, 1]), col = 'blue', xlab = expression(beta[0])))
with(simdata1000, abline(v = par_beta, col = 'red', lty = 2))

with(finalRun, plot(density(draws[, 1]), col = 'blue', xlab = expression(gamma)))
with(simdata1000, abline(v = par_gamma, col = 'red', lty = 2))

#' Comments
#' Beta estimated well
#' Gamma Underestimated Slightly

save.image(file = "simdata1000_centered.RData")
