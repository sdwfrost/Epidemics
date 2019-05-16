#'
#' Deterministic Epidemic Solutions
#'


# ==== Logisitic Growth Epidemic ====

Logistic_y_t = function(N, y_0, beta, t){
  if(t < 0){
    stop("Time t must be non-negative.")
  }
  return(N*y_0/(y_0 + (N - y_0)*exp(-beta*N*t)))
}

Logistic_end = function(N, y_0, beta){
  return((1/(beta*N))*log((N-1)*(N - y_0)/y_0))
}

Logistic_Epidemic = function(N, y_0, beta, time_res, t_lim, PLOT = TRUE){

  times = seq(t_lim[1], t_lim[2], by = time_res)

  y_t = sapply(X = times, function(X) Logistic_y_t(N, y_0, beta, X))

  x_t = N - y_t

  full_infection = min(times[y_t > N - 1])

  if(PLOT){
    plot(times, y_t, ylim = c(0,N), col = "red", type = 'l')
    lines(times, x_t, col = "blue", type = 'l', lty = 2)
    abline(v = full_infection)
  }

  return(list(x_t = x_t, y_t = y_t, t = times))
}

Logistic_EC = function(N, y_0, beta, time_res, t_lim, PLOT = TRUE){

  LE = Logistic_Epidemic(N, y_0, beta, time_res, t_lim, PLOT = FALSE)

  dy_dt = beta*LE$y_t*LE$x_t

  if(PLOT){
    max_rate = max(dy_dt)
    plot(LE$t, dy_dt, ylim = c(-max_rate,max_rate), col = "red", type = 'l')
    lines(LE$t, -dy_dt, col = "blue", type = 'l', lty = 2)
  }
  return(list(dy_dt = dy_dt, t = LE$t))
}

Logistic_max_EC = function(N, y_0, beta){
  return(1/(beta*N)*log((N - y_0)/y_0))
}

# ==== Interacting Groups ====

#' A general solution for a logistic growth model with interacting
#' groups is only possible with numerical integration methods (Eulers, Runge-Kutta)



Logistic_group_epidemic = function(N, y_0, beta, kappa, time_res, time_lim){




}



