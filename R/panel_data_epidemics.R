#'
#' Transition Panel Data for a General Epidemic Event Table
#'

#' Function aim:
#'
#' Take a event table, which is output from an Epidemic Gillespie algorithm,
#' and return panel data on the number of each possible transitions that happened
#' between the specified observation points.
#'
#' An event table will have a few things. For each event, there will be ...
#'
#'
#'  Time of event
#'  ID of individual
#'  State of Individual after event
#'  State of Individual before event
#'  Covariates of individual
#'

#' Then we can give panel_data.epidemics the event table
#' perhaps tell it a subset of individuals we want to observe
#' also, and got transition panel data out.
#'

# One time period

# Starting State of Individual
# Ending State of Individual

#' function to keep transitions with zero transition count in table


state_table = function(x, y, states = NULL, levels){
  if(is.null(states)){
    states = levels
  }
  if(missing(levels)){
    min(x,y):max(x,y)
  }
  if(length(states) != length(levels) & !is.null(states)){
    stop("Number of states and levels must be equal")
  }
  return(table(factor(x, levels, states), factor(y, levels, states)))
}


transition_table.epidemics = function(subject, state, prev_state, levels, state_names = NULL, output = "table"){
  ID = unique(subject)
  start_state = prev_state[sapply(ID, function(ID) head(which(subject == ID), n = 1))]
  end_state = state[sapply(ID, function(ID) tail(which(subject == ID), n = 1))]
  trans = state_table(start_state, end_state, state_names, levels)
  #trans = table(start_state, end_state)

  names(dimnames(trans)) = c("t_1", "t_0")

  if(output == "table"){
    return(trans)
  } else if(output == "freq"){
    trans = as.data.frame(trans)
    #trans_names =  sapply(X = 1:nrow(trans), function(X) paste(trans[X,1:2], sep = "", collapse = ""))
    return(trans)
  }
}

panel_data.epidemics = function(times, subject, state, prev_state, state_names, T_obs, k, obs_times = NULL, data = NULL){

  if(!is.null(data)){
    data = as.data.frame(data)
    times = eval(substitute(times), data, parent.frame())
    state = eval(substitute(state), data, parent.frame())
    prev_state = eval(substitute(prev_state), data, parent.frame())
  }

  n = length(times)

  if(!is.null(data)){
    subject = if(missing(subject)){
              rep(1, n)
            } else{
              eval(substitute(subject), data, parent.frame())
            }
  }

  if(is.null(obs_times)){
    obs_times = seq(T_obs[1], T_obs[2], length = k)
  }

  i = 1:length(obs_times)

  time_indices = lapply(i, function(i) which(obs_times[i] < times & times < obs_times[i + 1]))

  panel_data = lapply(X = time_indices,
                      function(X) transition_table.epidemics(subject[X], state[X], prev_state[X]))

}


