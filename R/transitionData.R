#' Function which transforms panel data into Transition Data


transitionData = function(panelData, states){
  lapply(X = 2:length(panelData), FUN = function(X){state_table(panelData[[X-1]][,2], panelData[[X]][,2], trans = possible_transitions(states))})
}
