#' fsMCMC Results

load("fsMCMC200Results.rds")

#' Take one example to show that algorithm works
#' and converges

par(mfrow = c(2, 2))
with(RUNS[[4]], plot(draws[,1], type = 'l', ylab = expression(beta), main = 'k = 6, 2nd half of epidemic'))
with(RUNS[[4]], acf(draws[,1], main = ''))

with(RUNS[[4]], plot(draws[,2], type = 'l', ylab = expression(gamma), main = 'k = 6, 2nd half of epidemic'))
with(RUNS[[4]], acf(draws[,2], main = ''))

#' Presenting the problem
#' Different observations of the same sample lead to different inferences.
#' Without knowing the whole process, we wouldn't know what observation times
#' would give the 'best' inference. What constitutes as 'best' in this case?

#' Interested in learning about the process but also want good parameter estimates
#' and good prediction. The optimal experimental design may be different for these
#' two proposals.
par(mfrow = c(1, 2))
for(j in 1:2){
  if(j == 1){
    xlab = expression(beta)
    ylim = c(0,1000)
    v = simdata200$par_beta
  } else{
    xlab = expression(gamma)
    ylim = c(0, 10)
    v = simdata200$par_gamma
  }
  with(RUNS[[1]], plot(density(draws[,j]), ylim = ylim, type = 'l', main = "", xlab = xlab))
  for(i in 2:length(RUNS)){
    with(RUNS[[i]], lines(density(draws[,j]), col = i, lty = i))
  }

  #' TRUE PARAMETER VALUE
  abline(v = v, col = 'black', lty = 2)

  #' LEGEND
  legend("topright", c("k = 6, Whole", "k = 12, Whole", "k = 4, Whole",
                       "k = 6, 2nd Half", "k = 6, 1st Half"), col = 1:length(RUNS), lty = 1:length(RUNS),
         cex = 0.6, xjust = 1)
}














