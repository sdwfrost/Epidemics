#'
#' SIS Functions
#'
#'


SIS_y_soln = function(N, x_0, y_0, beta, gamma, t_lim, t_res){
  t_seq = seq(0, t_lim, by = t_res)

  y = (N*beta - gamma)*y_0/((beta*x_0 - gamma)*exp(-t_seq) + beta*y_0)

  y[y < 0] = 0
  y[y > N] = N

  x = N - y

  return(list(x = x, y = y, t = t_seq))
}




