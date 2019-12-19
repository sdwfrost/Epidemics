set.seed(1)
#' Epidemic Initial Conditions
N = 200
gamma = 1
R0 = 1.25
beta = gamma*R0/N
I_0 = 40
initialState = c(rep(1, N - I_0), rep(2, I_0))


#' Observation Parameters
m = 0.1*N
k = 10
lastObs = 5


set.seed(1)
obsTimes = seq(0, lastObs, length = k)
SIS_sim = homogeneousPanelDataSIS_Gillespie(initialState, beta, gamma, obsTimes)


#' Sample Panel Data
panelData = panelDataSample(SIS_sim$panelData, m = m)
Y = transitionData(panelData, states = 1:2)


noIts = 10000
lambda = 5*10^(-5)
V = diag(c(1/N, 1))
thetaLim= rep(Inf, 2)
noDraws = 1000


adaptiveBlockedIS_run = adaptiveSISfsMCMC7_blockedIS(Y, I_0, SIS_sim$obsTimes, N, beta0 = beta, gamma0 = gamma, thetaLim, lambda0 = lambda,
                                                     noDraws = noDraws, blockSize = 30, noIts = 10000, delta = 0.05)

blockedIS_run = SIS_fsMCMC_blockedIS(Y, I_0, SIS_sim$obsTimes, N, beta0 = beta, gamma0 = gamma, thetaLim, lambda, V = adaptiveBlockedIS_run$V, noDraws,
                                     blockSize = adaptiveBlockedIS_run$blockSize, noIts)


#' Not Blocked


adaptiveRun = adaptiveSISfsMCMC7(Y, I_0, SIS_sim$obsTimes, N, beta0 = beta, gamma0 = gamma, thetaLim, lambda0 = lambda,
                                                     noDraws = noDraws, s = 150, noIts = 10000, delta = 0.05)

run = SIS_fsMCMC(Y, I_0, SIS_sim$obsTimes, N, beta0 = beta, gamma0 = gamma, thetaLim, lambda = adaptiveRun$optLambda,
                 V = adaptiveRun$optV, noDraws, s = 150, noIts)




tmp = tempfile()
Rprof(tmp, interval = 0.01)
adaptiveRun = adaptiveSISfsMCMC25(Y, I_0, SIS_sim$obsTimes, N, beta0 = beta, gamma0 = gamma, thetaLim, lambda0 = lambda,
                                  noDraws = noDraws,s = s, noIts = 10000, delta = 0.05)

Rprof(NULL)
summaryRprof(tmp)$by.self
summaryRprof(tmp)$sampling.time
E = rexp(10000)
U = runif(10000)
tmp2 = tempfile()
Rprof(tmp2, interval = 0.01)
adaptiveRun2 = adaptiveSISfsMCMC10(Y, I_0, SIS_sim$obsTimes, N, beta0 = beta, gamma0 = gamma, thetaLim, lambda0 = lambda,
                                   noDraws = noDraws,s = s, noIts = 1000, delta = 0.05)
Rprof(NULL)
head(summaryRprof(tmp2)$by.self)
summaryRprof(tmp2)$sampling.time

head(summaryRprof(tmp)$by.self)
head(summaryRprof(tmp2)$by.self)












