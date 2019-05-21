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


transition_table.epidemics = function(subject, times, state, prev_state, period, levels, state_names = NULL, output = "table"){

  relevant_events = period[1] <= times & times <= period[2]
  previous_events = times < period[1]

  ID = unique(subject)
  ID_r_e = unique(subject[relevant_events])

  ID_p_e = ID[!(ID %in% ID_r_e)]

  if(length(ID_p_e) > 0){
    state_ID_p_e = state[sapply(X = ID_p_e, function(X) tail(which(subject[previous_events] == X), n = 1))]
  } else{
    state_ID_p_e = NULL
  }

  if(length(ID_r_e) > 0){
    start_state = c(prev_state[relevant_events][sapply(X = ID_r_e, function(X) head(which(subject[relevant_events] == X), n = 1))],
                    state_ID_p_e)

    end_state = c(state[relevant_events][sapply(X = ID_r_e, function(X) tail(which(subject[relevant_events] == X), n = 1))],
                  state_ID_p_e)
  } else{
    start_state = state_ID_p_e

    end_state = state_ID_p_e
  }

  trans = state_table(start_state, end_state, state_names, levels)
  #trans = table(start_state, end_state)

  names(dimnames(trans)) = c("t_0", "t_1")

  if(output == "table"){
    return(trans)
  } else if(output == "freq"){
    trans = as.data.frame(trans)
    #trans_names =  sapply(X = 1:nrow(trans), function(X) paste(trans[X,1:2], sep = "", collapse = ""))
    return(trans$Freq)
  }
}

panel_data.epidemics = function(subject, times, state, prev_state, state_names, levels, T_obs, k, obs_times = NULL, data = NULL){
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
  if(missing(levels)){
    levels = 1:length(state_names)
  }
  if(is.null(obs_times)){
    obs_times = seq(T_obs[1], T_obs[2], length = k)
  }

  panel_data = lapply(X = 1:(length(obs_times) - 1),
                      function(X) transition_table.epidemics(subject, times, state,
                                                             prev_state, obs_times[c(X, X + 1)],
                                                             levels, state_names, output = "freq"))

  return(panel_data)
}


