#' Exploring SIS fsMCMC experiment results
#'

# ==== Preamble ====
library(ggplot2)
library(gridExtra)

load("bigRuns_SISIncreasingPanels.Rdata")

# ==== Scatterplots ====
par(mfrow = c(2, 3))

for(i in 1:length(finalRun)){
  plot(finalRun[[i]]$MCMCstep$draws[,-3], xlab = expression(beta), ylab = expression(gamma), pch = 2,
       main = paste(c(i + 9, "Panels", sep = " ")))
}

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

# ==== ESS & ESS/sec ====
betaESS = c()
betaESS.sec = c()

gammaESS = c()
gammaESS.sec = c()

minESS = c()
minESS.sec = c()
for(i in 1:length(finalRun)){
  betaESS[i] = finalRun[[i]]$MCMCstep$ESS[1]
  betaESS.sec[i] = finalRun[[i]]$MCMCstep$ESS.sec[1]

  gammaESS[i] = finalRun[[i]]$MCMCstep$ESS[2]
  gammaESS.sec[i] = finalRun[[i]]$MCMCstep$ESS.sec[2]

  minESS[i] = min(finalRun[[i]]$MCMCstep$ESS)
  minESS.sec[i] = min(finalRun[[i]]$MCMCstep$ESS.sec)
}

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

# ==== Boxplots ====
noPanels = 10:15
par(mfrow = c(1,1))
boxplot(finalRun[[1]]$MCMCstep$draws[-(1:10000),1], finalRun[[2]]$MCMCstep$draws[,1], finalRun[[3]]$MCMCstep$draws[-(1:10000),1],
        finalRun[[4]]$MCMCstep$draws[-(1:10000),1], finalRun[[5]]$MCMCstep$draws[-(1:10000),1], finalRun[[6]]$MCMCstep$draws[-(1:10000),1],
        xlab = "No Panels", ylab = expression(beta), names = noPanels)
abline(h = 1.4/200, lty  = 2, col =2)

par(mfrow = c(1,1))
boxplot(finalRun[[1]]$MCMCstep$draws[,2], finalRun[[2]]$MCMCstep$draws[,2], finalRun[[3]]$MCMCstep$draws[,2],
        finalRun[[4]]$MCMCstep$draws[,2], finalRun[[5]]$MCMCstep$draws[,2], finalRun[[6]]$MCMCstep$draws[,2],
        xlab = "No Panels", ylab = expression(gamma))
abline(h = 1, lty  = 2, col =2)

