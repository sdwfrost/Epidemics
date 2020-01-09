#' Final Runs for SIS Experiment

noIts = 100000
adapt = FALSE

# ==== load data ====

load("SISIncreasingPanels.RData")
rm(SISpanelExperiment)
adaptStep = results[seq(1, 11, by = 2)]

for(i in 1:6){
  adaptStep[[i]]$noPanels = i + 9
}

finalRun = lapply(X = adaptStep,
                  function(X) {SISpanelExperiment(noPanels = X$noPanels, lastObs = 5, noDraws = X$noDraws, lambda = X$lambda, V = X$V,
                                                  blockSize = X$blockSize, adapt = F, noIts = 1e5)})

rm(adaptStep, results)
save.image("bigRuns_SISIncreasingPanels.Rdata")

