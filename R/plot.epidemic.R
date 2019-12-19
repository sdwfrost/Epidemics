#'
#' Function for plotting the evolution of epidemic
#'


plot.epidemic = function(times, N, X, Y, Z, newPlot = TRUE){

  times = c(0, times[times != 0])

  if(newPlot){
    plot(times, Y, type = 'l', col = 'red', ylim = c(0, N),
         ylab = "No. Individuals", xlab = "time")

    lines(times, X, type = 'l', col = 'blue')

    if(!missing(Z)){
      lines(times, Z, type = 'l', col = 'purple')
      legend("topright", legend = c("X", "Y", "Z"), lty = c(1,1,1), col = c("blue", "red", "purple"), cex = 1)
    } else{
      legend("topright", legend = c("X", "Y"), lty = c(1,1,1), col = c("blue", "red"), cex = 1)
    }

  }else{
    lines(times, Y, type = 'l', col = 'red', ylim = c(0, N))

    lines(times, X, type = 'l', col = 'blue')

    if(!missing(Z)){
      lines(times, Z, type = 'l', col = 'purple')
    }
  }


}



