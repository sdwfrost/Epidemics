#'
#' Which Event Occurs next in an Epidemic (and to who)?
#'
#'



#' Samples the event which occurs in an epidemic.
#'
#' Two out outputs.
#'
#' event, equals 0 if infection and 1 if removal/recovery
#'
#' ID_index, the index of the individual the event occurs to.


event.epidemics = function(individual_inf_rate, gamma, Y, U = NULL){
  if(is.null(U) | is.logical(is.na(U))){
    U = runif(1, 0, 1)
  }
  X = length(individual_inf_rate)
  rates = c(individual_inf_rate, rep(gamma,Y))
  total_rate = sum(rates)
  ID_index = sum(cumsum(rates)/total_rate < U) + 1

  if(ID_index <= X){
    event = 0
  } else{
    event = 1
    ID_index = ID_index - X
  }
  return(list(event = event, ID_index = ID_index))
}
