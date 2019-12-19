#' PROFILING SIS GILLESPIES
devtools::load_all(".")

N = 200
gamma = 1
R0 = 4
beta = gamma*R0/N
initialState = c(rep(1, N - 1), rep(2, 1))
obsTimes = seq(0, 5, length = 6)

no_draws = 100
s = 10
UCurr = runif(no_draws)
ECurr = rexp(no_draws)
proposalSet = sample(1:no_draws, size = s)

EProp = ECurr
UProp = UCurr
EProp[proposalSet] = rexp(s)
UProp[proposalSet] = runif(s)

sum(dexp(UProp, rate = 1, log = TRUE)) + sum(dexp(UCurr[proposalSet], rate = 1, log = TRUE)) -
  (sum(dexp(UCurr, rate = 1, log = TRUE)) + sum(dexp(UProp[proposalSet], rate = 1, log = TRUE)))




set.seed(1)
x = homogeneousPanelDataSIS_Gillespie(initialState, beta, gamma, obsTimes)

set.seed(1)
x = homogeneousPanelDataSIS_GillespieEU(initialState, beta, gamma, obsTimes, E = rexp(2000), U = runif(2000))

x$panelData[[6]] == x1$panelData[[6]]

panelTemplate = matrix(nrow = N, ncol = 2)

tmp = tempfile()
Rprof(tmp, interval = 0.01)
replicate(1000, homogeneousPanelDataSIS_Gillespie(initialState, beta, gamma, obsTimes))
Rprof(NULL)
summaryRprof(tmp)$by.self
summaryRprof(tmp)$sampling.time
E = rexp(10000)
U = runif(10000)
tmp2 = tempfile()
Rprof(tmp2, interval = 0.01)
replicate(1000, homogeneousPanelDataSIS_GillespieEU(initialState, beta, gamma, obsTimes, E , U))
Rprof(NULL)
summaryRprof(tmp2)$by.self
summaryRprof(tmp2)$sampling.time
