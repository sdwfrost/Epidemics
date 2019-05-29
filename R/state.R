#' State Function
#'
#' Take multiple sets of individuals and convert into a vector of states



state = function(...){
  args = list(...)
  i = 1:length(args)
  #state = c(rep(1, length(S)), rep(2, length(I)), rep(3, length(R)))
  state = mapply(function(args, i) rep(i, length(args)), args, i, SIMPLIFY = F)
  state_data = data.frame(ID = unlist(args), state = unlist(state))
  state_data = dplyr::arrange(state_data, ID)
  return(state_data$state)
}
