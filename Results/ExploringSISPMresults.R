#' Exploring SIS/SIR


# ==== Preamble ====
library(ggplot2)
library(gridExtra)

load("SIS_PMnoSimsExperiment.Rdata")
load("bigRuns_SISIncreasingPanels.Rdata")



#' Plot of Test Epidemic
set.seed(1)
N = 200
gamma = 1
R_0 = 1.25
beta = gamma*R_0/(N)
I_0 = 40
initialState = c(rep(1, N - I_0), rep(2, I_0))

m = 0.1*N
noPanels = 10
obsTimes = seq(0, lastObs, length = noPanels)
set.seed(1)
SIS_sim = homogeneousPanelDataSIS_Gillespie(initialState, beta, gamma, obsTimes)

set.seed(1)
SIS_sim_events = individualHomogeneousSIS_Gillespie(initialState, beta, gamma, obsTimes = 5)
jpeg("SIStestEpidemic.jpeg", width = 960, height = 960)
par(mfrow = c(1,1))
with(SIS_sim_events, plot(eventTable[,c(1,5)], type = 'l', col = "red", ylim = c(0, 100),
                          xlab = "Time", ylab = "No. Infectives"))
dev.off()

#' Mixing Improving with No. Simulations
for(i in 1:5){
  jpeg(paste(c('PM_MCMC_SIS', i, '.jpeg'), sep = "", collapse = ""), width = 1.25*960, height = 960)
  par(mfrow = c(1, 2))
  plot(SIS_PMexperiment10[[i]]$MCMCrun$draws[, 1], type = 'l', ylab = expression(beta))
  acf(SIS_PMexperiment10[[i]]$MCMCrun$draws[, 1], main = "")
  dev.off()
}


jpeg("PM_MCMC_PosteriorSamples10.jpeg", width = 1.25*960, height = 960)
par(mfrow = c(1,1))
with(SIS_PMexperiment10[[5]]$MCMCrun, plot(draws[,1:2], pch = ".", cex = 2,
                                           xlab = expression(beta), ylab = expression(gamma),
                                           main = "PM-MCMC Posterior Samples (10 Panels)"))
with(SIS_PMexperiment10[[5]]$MCMCrun, points(betaSummary[1], gammaSummary[1], col = "blue", cex = 5, pch = "*"))
points(0.007, 1, col = "red", cex = 5, pch = "*")

legend("topleft", legend =c("Posterior Mean", "Truth"), col = c("blue", "red"), pch = "*",
       cex = 1)
dev.off()



#' ESS + ESS/sec
jpeg("PM_MCMC_ESS10.jpeg", width = 1.25*960, height = 960)
ESS = c()
ESS.sec = c()
noSims = 1:5
for(i in 1:5){
  ESS[i] = min(SIS_PMexperiment10[[i]]$MCMCrun$ESS)
  ESS.sec[i] = min(SIS_PMexperiment10[[i]]$MCMCrun$ESS.sec)
}
par(mfrow = c(1,2))
plot(noSims, ESS, xlab = "No. Simulations", ylab = "Effective Sample Size (N = 10000)",
     pch = "x", cex  = 3)
plot(noSims, ESS.sec, xlab = "No. Simulations", ylab = "ESS/sec", pch = "x", cex = 3)
dev.off()



#' fsMCMC mixing

jpeg("fsMCMCmixing.jpeg", width = 1.25*960, height = 960)
par(mfrow = c(1,2))
plot(finalRun[[1]]$MCMCstep$draws[1:10000, 1], type = 'l', ylab = expression(beta))
acf(finalRun[[1]]$MCMCstep$draws[1:10000, 1], main = "")
dev.off()

#' fsMCMC Posterior Samples

jpeg("fsMCMC_PosteriorSamples10.jpeg", width = 1.25*960, height = 960)
par(mfrow = c(1,1))
with(finalRun[[1]]$MCMCstep, plot(draws[1:10000,1:2], pch = ".", cex = 2,
                                           xlab = expression(beta), ylab = expression(gamma),
                                           main = "fsMCMC Posterior Samples (10 Panels)"))
with(finalRun[[1]]$MCMCstep, points(betaSummary[1], gammaSummary[1], col = "blue", cex = 5, pch = "*"))
points(0.007, 1, col = "red", cex = 5, pch = "*")
with(SIS_PMexperiment10[[5]]$MCMCrun, points(betaSummary[1], gammaSummary[1], col = "purple", cex = 5, pch = "*"))
legend("topleft", legend =c("fsMCMC Mean", "Truth", "PM-MCMC Mean"), col = c("blue", "red", "purple"), pch = "*",
       cex = 1.2)
dev.off()

#' ESS + ESS/sec + fsMCMC ESS and ESS/sec
jpeg("fsMCMC_ESS10.jpeg", width = 1.25*960, height = 960)
ESS = c()
ESS.sec = c()
noSims = 1:5
for(i in 1:5){
  ESS[i] = min(SIS_PMexperiment10[[i]]$MCMCrun$ESS)
  ESS.sec[i] = min(SIS_PMexperiment10[[i]]$MCMCrun$ESS.sec)
}
par(mfrow = c(1,1))
plot(noSims, ESS.sec, xlab = "No. Simulations", ylab = "Effective Sample Size (N = 10000)",
     pch = "x", cex  = 2, xlim = c(1, 6), col = "red")
points(6, min(finalRun[[1]]$MCMCstep$ESS.sec), col = "blue", pch = "x", cex = 2)
legend("topright", legend = c("PM-MCMC", "fsMCMC"), col = c("red", "blue"), pch = "x",
       cex = 1.2)
dev.off()

# ==== Scatterplots ====
par(mfrow = c(2, 3))

SIS_PMexperiment10final = list()
SIS_PMexperiment11final = list()
SIS_PMexperiment12final = list()
SIS_PMexperiment13final = list()
SIS_PMexperiment14final = list()
SIS_PMexperiment15final = list()

for(i in 1:length(SIS_PMexperiment10)){
  SIS_PMexperiment10final[[i]] = SIS_PMexperiment10[[i]]$MCMCrun
  SIS_PMexperiment11final[[i]] = SIS_PMexperiment11[[i]]$MCMCrun
  SIS_PMexperiment12final[[i]] = SIS_PMexperiment12[[i]]$MCMCrun
  SIS_PMexperiment13final[[i]] = SIS_PMexperiment13[[i]]$MCMCrun
  SIS_PMexperiment14final[[i]] = SIS_PMexperiment14[[i]]$MCMCrun
  SIS_PMexperiment15final[[i]] = SIS_PMexperiment15[[i]]$MCMCrun
}

for(i in 1:length(SIS_PMexperiment10final)){
  SIS_PMexperiment10final[[i]]$noPanels = 10
  SIS_PMexperiment10final[[i]]$noSims = i
}
for(i in 1:length(SIS_PMexperiment11final)){
  SIS_PMexperiment11final[[i]]$noPanels = 11
  SIS_PMexperiment11final[[i]]$noSims = i
}
for(i in 1:length(SIS_PMexperiment12final)){
  SIS_PMexperiment12final[[i]]$noPanels = 12
  SIS_PMexperiment12final[[i]]$noSims = i
}
for(i in 1:length(SIS_PMexperiment13final)){
  SIS_PMexperiment13final[[i]]$noPanels = 13
  SIS_PMexperiment13final[[i]]$noSims = i
}
for(i in 1:length(SIS_PMexperiment14final)){
  SIS_PMexperiment14final[[i]]$noPanels = 14
  SIS_PMexperiment14final[[i]]$noSims = i
}
for(i in 1:length(SIS_PMexperiment15final)){
  SIS_PMexperiment15final[[i]]$noPanels = 15
  SIS_PMexperiment15final[[i]]$noSims = i
}

# ==== Results Dataset ====

PMresults = data.frame(matrix(nrow = 6*length(SIS_PMexperiment10final), ncol = 11))

colnames(PMresults) = c("betaESS", "gammaESS", "betaESS.sec", "gammaESS.sec", "betaMean", "betaVar", "gammaMean", "gammaVar",
                      "timeTaken", "noPanels", "noSims")

for(i in 1:length(SIS_PMexperiment10final)){
  #' 10 Panels
  PMresults[i, ] = unlist(SIS_PMexperiment10final[[i]][-(1:2)])
  #' 11 Panels
  PMresults[i + 5, ] = unlist(SIS_PMexperiment11final[[i]][-(1:2)])
  #' 12 Panels
  PMresults[i + 10, ] = unlist(SIS_PMexperiment12final[[i]][-(1:2)])
  #' 13 Panels
  PMresults[i + 15, ] = unlist(SIS_PMexperiment13final[[i]][-(1:2)])
  #' 14 Panels
  PMresults[i + 20, ] = unlist(SIS_PMexperiment14final[[i]][-(1:2)])
  #' 15 Panels
  PMresults[i + 25, ] = unlist(SIS_PMexperiment15final[[i]][-(1:2)])
}

rm(list = ls()[!(ls() == "PMresults")])

# ==== ESS & ESS/sec ====

par(mfrow = c(2, 3))
for(i in 1:6){
  with(PMresults, plot(noSims[1:5 + 5*(i - 1)], pmin(betaESS.sec[1:5 + 5*(i - 1)], gammaESS.sec[1:5 + 5*(i - 1)]) , xlab = "No Sims", ylab = "Min ESS",
                     main = paste(c(9 + i, "Panels"), sep = "", collapse = " ")))
}

par(mfrow = c(2, 3))
for(i in 1:6){
  with(PMresults, plot(noSims[1:5 + 5*(i - 1)], pmin(betaESS.sec[1:5 + 5*(i - 1)], gammaESS.sec[1:5 + 5*(i - 1)]) , xlab = "No Sims", ylab = "Min ESS/sec",
                     main = paste(c(9 + i, "Panels"), sep = "", collapse = " ")))
}


ggplot(data = PMresults, aes(x = noPanels, y = noSims)) + geom_raster(aes(fill = gammaESS))
ggplot(data = PMresults, aes(x = noPanels, y = noSims)) + geom_raster(aes(fill = pmin(betaESS.sec, gammaESS.sec)))

# Compare to SIS fsMCMC
load("bigRuns_SISIncreasingPanels.Rdata")

fsMCMCresults

fsMCMCresults = data.frame(matrix(nrow = 6, ncol = 10))
colnames(fsMCMCresults) = c("betaESS", "gammaESS", "betaESS.sec", "gammaESS.sec", "betaMean", "betaVar", "gammaMean", "gammaVar",
                            "timeTaken", "noPanels")

for(i in 1:length(finalRun)){
  #' 10 Panels
  fsMCMCresults[i, ] = c(unlist(finalRun[[i]]$MCMCstep[-(1:2)]), i)
}




#' 10 Panels
par(mfrow = c(2, 3))
for(i in 1:6){
  with(PMresults, plot(noSims[1:5 + 5*(i - 1)], pmin(betaESS.sec[1:5 + 5*(i - 1)], gammaESS.sec[1:5 + 5*(i - 1)]) , xlab = "No Sims", ylab = "Min ESS",
                       main = paste(c(9 + i, "Panels"), sep = "", collapse = " ")))
  with(fsMCMCresults, points(, pmin(betaESS.sec[i], gammaESS.sec[i]), pch = 4, col = "red"))
}

par(mfrow = c(1,1))
with(PMresults,plot(10, max(pmin(betaESS.sec[1:5], gammaESS.sec[1:5])), xlab = "No Panels", ylab = "Min ESS/sec",
                    xlim = c(10,15),ylim = c(0, max(betaESS.sec, gammaESS.sec))))

for(i in 1:6){
  with(PMresults,points(9 + i, max(pmin(betaESS.sec[1:5 + 5*(i-1)], gammaESS.sec[1:5 + 5*(i-1)]))))
  with(fsMCMCresults, points(9 + i, max(pmin(betaESS.sec[i], gammaESS.sec[i])), pch = 4, col ="red"))
}

with(PMresults,  max(pmin(betaESS.sec[1:5 + 5*(i-1)], gammaESS.sec[1:5 + 5*(i-1)])))

with(PMresults,  max(betaESS.sec, gammaESS.sec))



par(mfrow = c(2,2))
#' Beta
plot(noPanels, betaESS)
plot(noPanels, betaESS.sec)

#' Gamma
plot(noPanels, gammaESS)
plot(noPanels, gammaESS.sec)

par(mfrow = c(1, 2))
#' Min ESS & ESS.sec
plot(noPanels, minESS)
plot(noPanels, minESS.sec)


# ==== 2D Densities ====
grid.arrange(plot1, plot2, ncol=2)

xlabel = xlab(expression(beta))
ylabel = ylab(expression(gamma))

data1 = as.data.frame(finalRun[[1]]$MCMCstep$draws)
plot1 <- ggplot(data1, aes(x= V1, y= V2)) +
  geom_density_2d() + title("10 panels")

data2 = as.data.frame(finalRun[[2]]$MCMCstep$draws)
plot2 <- ggplot(data2, aes(x= V1, y= V2)) +
  geom_density_2d() + xlabel + ylabel

data3 = as.data.frame(finalRun[[3]]$MCMCstep$draws)
plot3 <- ggplot(data3, aes(x= V1, y= V2)) +
  geom_density_2d() + xlabel + ylabel

data4 = as.data.frame(finalRun[[4]]$MCMCstep$draws)
plot4 <- ggplot(data4, aes(x= V1, y= V2)) +
  geom_density_2d() + xlabel + ylabel

data5 = as.data.frame(finalRun[[5]]$MCMCstep$draws)
plot5 <- ggplot(data5, aes(x= V1, y= V2)) +
  geom_density_2d() + xlabel + ylabel

data6 = as.data.frame(finalRun[[6]]$MCMCstep$draws)
plot6 <- ggplot(data6, aes(x= V1, y= V2)) +
  geom_density_2d()

grid.arrange(plot1, plot2,plot3, plot4, plot5, plot6, ncol=3)

SIRepidemic = homogeneousPanelDataSIR_GillespieEU(c(rep(1, 190), rep(2, 10)), beta = 0.06, gamma = 3,
                                                  obsTimes = seq(0, 10, length = 6),
                                                  E = rexp(2*200 - 10), U = runif(2*200 - 10))

# ==== Marginal Density Plots ====

par(mfrow = c(2, 3))

plot(density(finalRun[[1]]$MCMCstep$draws[,1]))
abline(v = 1.4/200, lty = 2, col = 2)
for(i in 2:length(finalRun)){
  plot(density(finalRun[[i]]$MCMCstep$draws[,1]))
  abline(v = 1.4/200, lty = 2, col = 2)
}

plot(density(finalRun[[1]]$MCMCstep$draws[,2]))
abline(v = 1, lty = 2, col = 2)
for(i in 2:length(finalRun)){
  plot(density(finalRun[[i]]$MCMCstep$draws[,2]))
  abline(v = 1, lty = 2, col = 2)
}

# ==== Boxplots ====

noPanels = 10:15
par(mfrow = c(1,1))

boxplot(finalRun[[1]]$MCMCstep$draws[-(1:1000),1], finalRun[[2]]$MCMCstep$draws[,1], finalRun[[3]]$MCMCstep$draws[-(1:10000),1],
        finalRun[[4]]$MCMCstep$draws[-(1:1000),1], finalRun[[5]]$MCMCstep$draws[-(1:10000),1], finalRun[[6]]$MCMCstep$draws[-(1:10000),1],
        xlab = "No Panels", ylab = expression(beta), names = noPanels)
abline(h = 1.4/200, lty  = 2, col =2)

par(mfrow = c(1,1))
boxplot(finalRun[[1]]$MCMCstep$draws[,2], finalRun[[2]]$MCMCstep$draws[,2], finalRun[[3]]$MCMCstep$draws[,2],
        finalRun[[4]]$MCMCstep$draws[,2], finalRun[[5]]$MCMCstep$draws[,2], finalRun[[6]]$MCMCstep$draws[,2],
        xlab = "No Panels", ylab = expression(gamma))
abline(h = 1, lty  = 2, col =2)

# ==== Mean & Variance ====

betaMean = c()
betaVar = c()

gammaMean = c()
gammaVar = c()


for(i in 1:length(finalRun)){
  betaMean[i] = finalRun[[i]]$MCMCstep$betaSummary[1]
  betaVar[i] = finalRun[[i]]$MCMCstep$betaSummary[2]

  gammaMean[i] = finalRun[[i]]$MCMCstep$gammaSummary[1]
  gammaVar[i] = finalRun[[i]]$MCMCstep$gammaSummary[2]
}

par(mfrow = c(2,2))
noPanels = 10:15

plot(noPanels, betaMean, ylab = expression(hat(beta)))
abline(h = 1.4/200, lty = 2, col = "red")

plot(noPanels, betaVar, ylab = expression(Var(beta)))

plot(noPanels, gammaMean, xlab = expression(hat(gamma)))
abline(h = 1, lty = 2, col = "red")

plot(noPanels, gammaVar, ylab = expression(Var(gamma)))

