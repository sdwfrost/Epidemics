#' Importance Sampling with known epidemic parameters (SIR)

#' Epidemic process {x_t ; t \in (0, T)}
#' T is the end of the epidemic
#' x_t is the infectious state of every individual
#' in the population at time t
#'
#' The process is assumed to evolve according to
#' parameters \theta
#'
#' Observed Panel Data {y_t ; t \in (t_1, t_2, ..., t_k)}
#' y_t is the infectious state of a sample of size m from the
#' population

#' y_t are independent conditional on the epidemic process at
#' time t.

# = Distributions =

#' We would like to receive samples from the posterior
#' distribution
#'
#' \pi(x_{0:t}| y_{1:t}) \propto \pi(y_{1:t}| x_{0:t}) \pi(x_{0:t})
#'

#' Do this via importance sampling
#'
#' 1. Propose $x^{i}$ \sim q() (This is a particle)
#' 2. Calculate importance weight \omega(x^{i})
#' (Although this is not possible is cannot calculate \pi(y_{1:t}))
#' 3. Repeat 1-2 N times
#' 4. Calculate Normalised weights \tilde{\omega}(x^{i})
#' (This is possible without \pi(y_{1:t}))
#' 5. Use this weighted sample to estimate integrals involving
#'    posterior distribution.

#' weights \omega are the ratio of the posterior and proposal,
#' measuring the difference telling us how representative the
#' sample is to the posterior. If proposals were made directly
#' from the posterior, all weights would be 1.

#' Diagnostic of the Importance Sampler
#' Effective Sample Size
#'
#' $$ ESS := 1/sum_{i = 1}^{N} \tilde{\omega}(x^{i})^{2} $$
#'
#' How many direct posterior samples the sample
#' obtained through importance sampling is worth
#'
#' If sampled directly from the posterior, all
#' weights would be 1
#'
#'

#' The epidemic process is easily simulated using the Gillespie
#' algorithm. This means using a proposal distribution of the
#' form q() = \pi(x_(0:t)|\theta) is trival. Although this is
#' not really informed by the data
#' = (Can we bias simulation of particle x^{i} according to y?) =

epidemicImportanceSampler <- function(transData, N, a, obsTimes, theta, Nparticles){

  start = as.numeric(Sys.time())
  # 1. Simulate N particles $x^{i}$ \sim q() (This is a particle)
  m = sum(transData[[1]])

  particles = lapply(X = 1:Nparticles, function(X){
    homogeneousPanelDataSIR_Gillespie(initialState = c(rep(1, N - a), rep(2, a)), beta = theta[1], gamma = theta[2],
                  obsTimes = obsTimes)$panelData
  })

  # 2. Calculate Likelihood (which is the main component of the weights)
  # extract transition data from particles
  particleTrans = lapply(particles, transitionData, states = 1:3)
  if(length(obsTimes) > 2){
    YgivenX = sapply(X = particleTrans, function(X) prod(mapply(extraDistr::dmvhyper, transData, X, MoreArgs = list(k = m))))
    } else{
      YgivenX = sapply(X = particleTrans, function(X) extraDistr::dmvhyper(transData[[1]], X[[1]], k = m))
    }

  # 3. Calculate Normalised weights
  ISweights = YgivenX/sum(YgivenX)

  # 4. Resample
  resample = sample(1:Nparticles, prob = ISweights, replace = T)
  particles = particles[resample]
  #ISweights = ISweights[resample]

  ESS = 1/sum(ISweights^2)
  timeTaken = as.numeric(Sys.time()) - start
  ESS.sec = ESS/timeTaken
  return(list(ESS = ESS, ESS.sec = ESS.sec, ISweights = ISweights, particles = particles, transData = transData,
              timeTaken = timeTaken))


  #' 4. Use this weighted sample to estimate integrals involving
  #'    posterior distribution.


}

sequentialImportanceSampler = function(transData, N, a, obsTimes, theta, Nparticles){

  # empty vector for components of marginal estimate
  marg = c()

  # Timestep 1, simulation of N particles
  start = as.numeric(Sys.time())
  m = sum(transData[[1]])
  # Simulate
  particles = lapply(X = 1:Nparticles, function(X){
    homogeneousPanelDataSIR_Gillespie(initialState = c(rep(1, N - a), rep(2, a)), beta = theta[1], gamma = theta[2],
                                      obsTimes = obsTimes[1:2])$panelData
  })

  # Calculate Transition data and conditional observation densities
  particleTrans = lapply(particles, transitionData, states = 1:3)

  ISweights = sapply(X = particleTrans, function(X) extraDistr::dmvhyper(transData[[1]], X[[1]], k = m))

  # 3. Calculate Normalised weights
  ISweightsNorm = ISweights/sum(ISweights)

  # 4. Resample
  resample = sample(1:Nparticles, prob = ISweightsNorm, replace = T)
  particles = particles[resample]

  for(i in 2:(length(obsTimes) - 1)){
    # 1. Simulate/Propogate Particles to the next panel observation time



    # a) Use Gillespie simulation which outputs panel data, run between ith and (i+1)th
    #    observation time.
    nextObs = lapply(X = particles , function(X){
      homogeneousPanelDataSIR_Gillespie(initialState = tail(X, n = 1)[[1]][,2], beta = theta[1], gamma = theta[2],
                                        obsTimes = obsTimes[i:(i+1)])$panelData[[2]]
    })
    particles = lapply(X = 1:Nparticles, function(X){
      particles[[X]][[i + 1]] = nextObs[[X]]
      return(particles[[X]])
    })

    # 2. Calculate Normalised Weights
    # a) Particle transition data
    particleTrans = lapply(particles, transitionData, states = 1:3)

    # b) Hypergeometric Calculation
    ISweights = sapply(X = particleTrans, function(X) extraDistr::dmvhyper(transData[[i]], X[[i]], k = m))

    # c) Normalise weights
    ISweightsNorm = ISweights/sum(ISweights)

    # 3. Estimate Marginal Density
    marg[i] = mean(ISweights)

    # 3. Resample
    resample = sample(1:Nparticles, prob = ISweights, replace = T)
    particles = particles[resample]
    #ISweights = ISweights[resample]
  }

  margEst = prod(marg)
  ESS = 1/sum(ISweights^2)

  timeTaken = as.numeric(Sys.time()) - start
  ESS.sec = ESS/timeTaken
  return(list(ESS = ESS, ESS.sec = ESS.sec, particles = particles, ISweights = ISweights,
              timeTaken = timeTaken))
}


