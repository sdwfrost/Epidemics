#'
#' Function for plotting the evolution of epidemic
#'


plot.epidemic = function(times, N, X, Y, Z){

  times = c(0, times[times != 0])

  plot(times, Y, type = 'l', col = 'red', ylim = c(0, N))

  lines(times, X, type = 'l', col = 'blue')

  if(!missing(Z)){
    lines(times, Z, type = 'l', col = 'purple')
  }
}



