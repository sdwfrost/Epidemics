#' state_table function
#'
#' Takes two vectors of states and the possible transitions,
#' then returns counts of unique transitions


state_table = function(x, y, trans){
  #if(is.null(states)){
  #  states = levels
  #}
  #if(missing(levels)){
  #  min(x,y):max(x,y)
  #}
  #if(length(states) != length(levels) & !is.null(states)){
  #  stop("Number of states and levels must be equal")
  #}
  #return(table(factor(x, levels, states), factor(y, levels, states)))
  freq = sapply(X = trans, function(X) sum(X[1] == x & X[2] == y))
}

