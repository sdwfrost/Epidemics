#'
#'
#' Test MCMC using Reduce function
#'

Normal_MCMC_proposal = function(current_state, lambda){

  x_curr = current_state[1]
  log_pi_curr = current_state[2]
  accept_x = current_state[3]

  #' Propose new beta and gamma using Multiplicative RW propsal
  x_prop = x_curr + lambda*rnorm(1)

  log_pi_prop = dnorm(x_prop, log = TRUE)

  log_a = (log_pi_prop) - (log_pi_curr)

  if(log(runif(1)) < log_a){
    log_pi_curr = log_pi_prop
    x_curr = x_prop
    accept_x = accept_x + 1
  }
  return(c(x_curr, log_pi_curr, accept_x))
}

Normal_MCMC = function(x_0, lambda, no_its, burn_in, lag_max = NA, thinning_factor = 1, reduce = F){
  Start = as.numeric(Sys.time())


  log_pi_0 = dnorm(x_0, log = TRUE)

  #' Proposal Acceptance Counter


  #' MCMC
  if(reduce){
    draws = matrix(unlist(Reduce(Normal_MCMC_proposal, init = c(x_0, log_pi_0, 0), x = rep(lambda, no_its),
                                 accumulate = TRUE)), ncol = 3, byrow = T)
    accept_x = tail(draws[,3], n = 1)
  } else{
    accept_x = 0
    draws = matrix(NA, nrow = no_its + 1, ncol = 2)
    draws[1, ] = c(x_0, log_pi_0)
    x_curr = x_0
    log_pi_curr = log_pi_0
    #print("Sampling Progress")
    #pb <- progress::progress_bar$new(total = no_its)
    for(i in 1:no_its){
      #pb$tick()
      #' Propose new beta and gamma using Multiplicative RW propsal
      x_prop = x_curr + lambda*rnorm(1)

      log_pi_prop = dnorm(x_prop, log = TRUE)

      log_a = (log_pi_prop) - (log_pi_curr)

      if(log(runif(1)) < log_a){
        log_pi_curr = log_pi_prop
        x_curr = x_prop
        accept_x = accept_x + 1
      }

      draws[i + 1, ] = c(x_curr, log_pi_curr)
    }
  }

  End <- as.numeric(Sys.time())

  time_taken <- End - Start

  draws <- draws[seq(from = burn_in + 1, to = (no_its + 1) - burn_in, by = thinning_factor), 1:2]

  #' Thin the samples

  # Calculate Effective Sample Sizes (and Per Second) and Acceptance Rates
  ESS <- coda::effectiveSize(draws[, 1])
  ESS_sec <- ESS/time_taken
  accept_rate <-  accept_x/no_its

  # = Plots =
  par(mfrow = c(1,2))

  # Plot Beta Samples and Sample Auto-Corrolation Function
  if(is.na(lag_max)){
    # x plots
    plot(draws[, 1], type = 'l', ylab = expression(x))
    acf(draws[, 1], main = "")
  } else{
    # x plots
    plot(draws[, 1], type = 'l', ylab = expression(x))
    acf(draws[, 1], lag_max, main = "")
  }

  #' Calculating Summary Statistics for samples
  x_summary = c(mean(draws[,1]), sd(draws[,1]))
  printed_output(rinf_dist = "-", no_proposals = "-", no_its, ESS, time_taken, ESS_sec, accept_rate)
  return(list(draws = draws, accept_rate = accept_rate, ESS = ESS, ESS_sec = ESS_sec, x_summary = x_summary,
              time_taken = time_taken))
}

start = as.numeric(Sys.time())
set.seed(1)
Test = Normal_MCMC(x_0 = rnorm(1), lambda = 4.1, no_its = 500000, burn_in = 1000, reduce = T)
time = as.numeric(Sys.time()) - start

plot(Test$draws[,1], exp(Test$draws[,2]))


