#'
#' General Epidemic Process
#'
#'

# ==== Purpose ====

#' To have one algorithm which deals with all constructions of an Epidemic.
#' For example, to be able to deal with;
#'  Multiple States
#'  The mechanics of transitioning between these states
#'  Track when events happen and what the nature of the event was.

# ==== Initial Thoughts ====

#' How many states?
#' How do individuals transition between these states?
#' Is it affected by the make up of other states?
#' How do you choose which event is next?
#'
#'
#' Choosing when next event is easy, E ~ exp(Total_rate_of_events)
#'
#' Choosing which event occurs is also okay, use sample with the individual transition rates and
#' use the index thrown out to chose which event occurs and which individual it happens to.
#'
#'

# ==== What we would need ====

#' (Homogeneous SIR context)
#' There are k states. (3 states)
#' Want to be able to name these states to give them meaning. (Susceptible, Infectious, Removed)
#'
#' Specify the mechanics of transitioning between these states (S ---> (beta) I ---> (gamma) R)
#'
#' Be able to draw the time of the next event (E ~ exp(total_inf_rate + total_rem_rate))
#'
#' Specify the probabilities of an event occuring and possibly which individual this event happens
#' to. (sample(c(S, I), prob = c(rep(length(I)*beta, length(S)), rep(gamma, length(I)))))
#'


# ==== Time to next event and which event functions ====


# == "population" object ==

#' Simply an object that will tell the algorithm what state an individual is in
#' at this particular time point.

# == "parameters" object ==

#' A list object which will hold the parameters for each transition process
#'
#'


# == "rates" object ==

#' A list of objects which will hold all the information about the rates at which an
#' individual transitions from one state to another.
#'
#' If the transition from one state to another is homogenous (i.e the same for all individuals,
#' the the object will be the parameters of the process)
#'
#' If the process depends on the relationship between individuals and aspects of themselves,
#' the object will be a (Sparse) Matrix.
#'
#' Use kernels to assist computing the rates
#'

rates = function(parameters, kernels){
  rates = mapply(function(kernels, parameters) kernel(parameters), kernels, parameters)
  return(rates)
}

#' Seems like drawing the time to next event and which event will occur will use
#' the same information, but in different ways. Either way, both functions will
#' have to weedle out which rates are relevant given what state people are in.
#' We do not want to have to decide this more than once, therefore there should
#' be a function which does this, and then we can pass the result to the functions.
#'

#' relevant_rates A function which will decided which rates are important, in terms of
#'                deciding when the next event occurs and which event occurs. The aim of
#'                this function is to avoid extracting this information twice.
#' @param rates Holds information about the rates ate which each individual transition
#'              between each state
#' @param population Holds information about what state each individual in the population
#'                   is in
#' @return relevant_rates
relevant_rates = function(rates, population){

}


#' time_to_next_event Draws the time to next event in the Infectious Disease process
#' @param rates Holds information about the rates ate which each individual transition
#'              between each state
#' @param population Holds information about what state each individual in the population is in
#' @param dist What distribution does the time_to_next_event follow (Exponential for a
#'             First Order Markov Process)
time_to_next_event = function(rates, population, dist = "Exp"){

}

#' which_event Decides which event will occur next in the Infectious Disease process
#' @param rates Holds information about the rates ate which each individual transition
#'              between each state.
#' @param population Holds information about what state each individual in the
#'                   population is in.
#'
#'
#'
which_event = function(rates, population){

}




# ==== The Algorithm ====

# == Rate Structure and Parameters ==

#' The parameters and the rate structure together will determine the rates of the
#' process for each individual
#'
#' States should be ordered in such a way which follows the flow of the epidemic.
#'
#' Therefore the rate structure should be ordered in the same way
#'
#' SIR Example:
#'
#' State 1 ---> S, State 2 ---> I, State 3 ---> R
#'
#' 1 ---> 2,  Beta Parameters
#'
#' 2 ---> 3,  Gamma Parameters
#'

General_Epidemic_Process = function(parameters, no_states, state_names, rate_structure,
                                    initial_population){

}





