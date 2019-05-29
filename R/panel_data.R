#'
#' Temporary Panel Data Generator
#'
#'


panel_data_new.epidemics = function(subject, times, state, prev_state,
                                    trans = possible_transitions(1:3), subset = NULL, obs_times){

  if(!is.null(subset)){
    subset_index = subject %in% subset
    subject = subject[subset_index]
    times = times[subset_index]
    state = state[subset_index]
    prev_state = prev_state[subset_index]
  }

  # Initial State
  initial_events = times == 0
  curr_state = state[initial_events]
  S = subject[initial_events][curr_state == 1]
  I = subject[initial_events][curr_state == 2]
  R = subject[initial_events][curr_state == 3]

  #' Panel Data
  panel_data = lapply(rep(NA, length(obs_times)), function(X) X)

  panel_data[[1]] = state_table(curr_state, curr_state, trans)
  for(i in 1:(length(obs_times) - 1)){

    time_index = which(obs_times[i] < times & times < obs_times[i + 1])

    if(length(time_index) > 0){
      next_subject = subject[time_index]
      next_state = state[time_index]
      next_prev_state = prev_state[time_index]

      for(j in 1:length(time_index)){
        if(next_state[j] == 3){
          individual = next_subject[j]
          I = I[!(individual == I)]
          R = c(R, individual)

        } else{
          individual = next_subject[j]
          S = S[!(individual == S)]
          I = c(I, individual)
        }
      }
    }

    old_state = curr_state
    curr_state = state(S, I, R)
    panel_data[[i + 1]] = state_table(old_state, curr_state, trans)
  }
  NA_panels = which(is.na(panel_data))
  no_events_transitions = state_table(curr_state, curr_state, trans)
  for(i in NA_panels){
    panel_data[[i]] = no_events_transitions
  }
  return(panel_data)
}
