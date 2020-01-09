#' Testing PM-MCMC + Adaptive
library(ggplot2)

# ==== Sample Panel Data ====

noPanels = 10
lastObs = 5

set.seed(1)
N = 200
gamma = 1
R_0 = 1.25
beta = gamma*R_0/(N)
I_0 = 40
initialState = c(rep(1, N - I_0), rep(2, I_0))

m = 0.1*N
obsTimes = seq(0, lastObs, length = noPanels)
SIS_sim = homogeneousPanelDataSIS_Gillespie(initialState, beta, gamma, obsTimes)

#' Sample Panel Data
panelData = panelDataSample(SIS_sim$panelData, m = m)
Y = transitionData(panelData, states = 1:2)


# ==== PM-MCMC ====

#' Do PM-MCMC runs to investigate its behaviour
lambda = 2.85/sqrt(2)
V = diag(c(1/N, 1))
noIts = 10000

newPseudoMarginalMCMC(Y, I_0, obsTimes, N, beta0 = beta, gamma0 = gamma, lambda = lambda, V = V, noSims = 2, noIts = noIts,
                      burnIn = 0)

adaptivePseudoMarginalMCMC(Y, I_0, obsTimes, N, beta0 = beta, gamma0 = gamma, noSims = 1, noIts = noIts,
                           burnIn = 0)

PM_MCMCexperiment(noSims = 1, noPanels= 5, lastObs = 5, lambda0 = 1e-5, noIts = 10000)
panelDataSim = sim$panelData
transDataSim = transitionData(panelDataSim, states = 1:2)
logPProp = dHyperGeom(obsTransData, transDataSim, noSampled, log = T)



# ==== Variance of No Draws when simulating Epidemic ====

betaSeq = seq(0.001, 0.015, length = 100)
gammaSeq = seq(0.5, 2, length = 100)

betaSeq[1]*N/gammaSeq[1]
betaSeq[100]*N/gammaSeq[100]


noDrawsMat = matrix(NA, nrow = 100, ncol = 100)
rownames(noDrawsMat) = betaSeq
colnames(noDrawsMat) = gammaSeq
for(i in 1:length(betaSeq)){
  noDrawsMat[i, ] = sapply(X = gammaSeq, function(X) homogeneousPanelDataSIS_Gillespie(initialState, betaSeq[i], X, obsTimes)$noDraws)
  print(i)
  }

library(ggplot2)
library(reshape2)

#' Puts Matrix A into long data format (i.e V1, V2, value)
longData <- melt(noDrawsMat)
View(longData)
#' Removes 0 value entries
longData <- longData[longData$value!=0, ]

ggplot(longData, aes(x = Var1, y = Var2)) +
  geom_raster(aes(fill = value)) +
  scale_fill_gradient(low = "lightblue", high = "red") +
  labs(x = expression(beta), y = expression(gamma), title = "No Events") +
  theme_bw() + theme(axis.text.x = element_text(size = 9, angle = 0, vjust = 0.3),
                     axis.text.y = element_text(size = 9),
                     plot.title = element_text(size = 11, hjust = 0.5))



