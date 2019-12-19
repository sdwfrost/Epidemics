#' dHyperGeom function
#'
#' Multivariate Hypergeometric Density Function
#'
#'

dHyperGeom = function(x, Y, k, log = TRUE){
  return(sum(mapply(extraDistr::dmvhyper, x, Y, MoreArgs = list(k, log))))
}

