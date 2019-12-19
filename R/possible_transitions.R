#' possible_transitions function
#'
#' Given a vector of states
#'

possible_transitions = function(states){
  possible_trans = expand.grid(states, states)
  possible_trans = lapply(X = 1:nrow(possible_trans), function(X) as.numeric(possible_trans[X,]))
  return(possible_trans)
}

