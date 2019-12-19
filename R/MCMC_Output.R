#'
#'
#' Functions for MCMC Output
#'
#'

# ==== Preamble ====



#' Trace and ACF plotting function

MCMCplots = function(samples, lagMax = NA, ylab, trueValue = NULL, set.mfrow = FALSE){

  if(set.mfrow){
    par(mfrow = c(1,2))
  }

  #' Trace Plot
  plot(samples, type = 'l', ylab = ylab)

  # ACF
  acf(samples, lagMax, ylab = ylab)

  #' Density Plot

  plot(density(samples), type = 'l', col = 'blue')

  if(!is.null(trueValue)){
    abline(v = trueValue)
  }

}




# ==== Printed Output for an MCMC run ====

printed_output <- function(rinf_dist, no_proposals, no_its, ESS, time_taken, ESS_sec, accept_rate){
  print(paste(c("Infectious Period Distribution:", rinf_dist), sep = "", collapse = ""))
  print(paste(c("No. Infection time proposals:", no_proposals)), sep = "", collapse = "")
  print(paste(c("Number of Iterations:", no_its)), sep = "", collapse = "")
  print(paste(c("Effective Sample Size:", min(ESS))), sep = "", collapse = "")
  print(paste(c("Time Taken:", time_taken)), sep = "", collapse = "")
  print(paste(c("Effective Sample Size per Second:", min(ESS_sec))), sep = "", collapse = "")
  print(paste(c("Acceptance Rate:", accept_rate)), sep = "", collapse = "")
}



