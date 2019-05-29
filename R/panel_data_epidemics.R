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

start_state = function(subject, prev_state){
  no_events = length(subject)
  first_event = 1

  if(no_events == 0){
    return(NULL)
  } else if(no_events > 1){
    for(i in 2:no_events){
      if(!(subject[i] %in% subject[1:(i-1)])){
        first_event = c(first_event, i)
      }
    }
  }
  return(prev_state[first_event])
}

end_state = function(subject, state){
  no_events = length(subject)
  last_event = NULL

  if(no_events == 0){
    return(NULL)
  } else if(no_events > 1){
    for(i in 1:(no_events - 1)){
      if(!(subject[i] %in% subject[(i+1):no_events])){
        last_event = c(last_event, i)
      }
    }
  }
  last_event = c(last_event, no_events)
  return(state[last_event])
}

transition_table.epidemics = function(subject, times, state, prev_state, period, levels, state_names = NULL, output = "table"){

  relevant_events = period[1] <= times & times <= period[2]


  ID = unique(subject)
  ID_relevant_events = unique(subject[relevant_events])
  previous_events = (times < period[1]) & !(subject %in% ID_relevant_events)
  ID_prev_events = ID[!(ID %in% ID_relevant_events)]

  if(sum(previous_events) > 0){
    #state_prev_events = end_state(subject[previous_events], state[previous_events])
    state_prev_events = state[sort(sapply(X = ID_prev_events,
                                       function(X) tail(which(subject[previous_events] == X),
                                                         n = 1)))]
  } else{
    state_prev_events = NULL
  }

  if(sum(relevant_events) > 0){
    #start_state = c(start_state(subject[relevant_events], prev_state[relevant_events]),
    #                state_prev_events)
    start_state = c(prev_state[relevant_events][sort(sapply(X = ID_relevant_events,
                                                       function(X) head(which(subject[relevant_events] == X),
                                                                        n = 1)))],
                    state_prev_events)

    #end_state = c(end_state(subject[relevant_events], state[relevant_events]),
    #              state_prev_events)
    end_state = c(state[relevant_events][sort(sapply(X = ID_relevant_events,
                                                     function(X) tail(which(subject[relevant_events] == X), n = 1)))],
                  state_prev_events)
  } else{
    start_state = state_prev_events

    end_state = state_prev_events
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

panel_data.epidemics = function(subject, times, state, prev_state, state_names, levels, subset = NULL, T_obs, k, obs_times = NULL, data = NULL){
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
  if(!is.null(subset)){
    subset_index = subject %in% subset
    subject = subject[subset_index]
    times = times[subset_index]
    state = state[subset_index]
    prev_state = prev_state[subset_index]
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


