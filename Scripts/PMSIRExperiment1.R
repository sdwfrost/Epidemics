#' PM SIR Experiment

# ==== PM SIR Experiment 1 (NOT PRIORITY) ====
noIts = 100000
adapt = TRUE
N = 200
noSims = 1:5

PMSIRexperiment5panels = lapply(X = noSims,
                        function(X) {SIR_PM_MCMCexperiment1(noSims = X, noPanels = 5, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})

PMSIRexperiment6panels = lapply(X = noSims,
                        function(X) {SIR_PM_MCMCexperiment1(noSims = X, noPanels = 6,lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})

PMSIRexperiment7panels = lapply(X = noSims,
                        function(X) {SIR_PM_MCMCexperiment1(noSims = X, noPanels = 7, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})
PMSIRexperiment8panels = lapply(X = noSims,
                        function(X) {SIR_PM_MCMCexperiment1(noSims = X, noPanels = 8, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})
PMSIRexperiment9panels = lapply(X = noSims,
                        function(X) {SIR_PM_MCMCexperiment1(noSims = X, noPanels = 9, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})
PMSIRexperiment10panels = lapply(X = noSims,
                              function(X) {SIR_PM_MCMCexperiment1(noSims = X, noPanels = 10, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})

save.image("SIRPMnoSimsExperiment5Panels.Rdata")

rm(list = ls())

# ==== fsMCMC SIR Experiment 1 ====

noIts = 100000
adapt = TRUE
N = 200
noPanels = 5:15

fsMCMCSIRexperiment5panels = lapply(X = noPanels,
                                function(X) {SIRfsMCMCexperiment1(noPanels = X, lastObs = 5, lambda0 = 1e-3,
                                                                    V0 = diag(c(1/N, 1)), blockSize = 100, adapt = T, noIts = noIts)})

save.image("SIRfsMCMCnoPanelsExperiment.Rdata")

rm(list = ls())

# ==== PM SIS Experiment (Run one for 100000 to exhibit behaviour) ====


noIts = 10000
adapt = TRUE
N = 200
noSims = 1:5

SIS_PMexperiment5 = lapply(X = noSims,
                        function(X) {PM_MCMCexperiment(noSims = X, noPanels = 5, lastObs = 5, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})

SIS_PMexperiment6 = lapply(X = noSims,
                        function(X) {PM_MCMCexperiment(noSims = X, noPanels = 6, lastObs = 5, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})

SIS_PMexperiment7 = lapply(X = noSims,
                        function(X) {PM_MCMCexperiment(noSims = X, noPanels = 7, lastObs = 5, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})
SIS_PMexperiment8 = lapply(X = noSims,
                        function(X) {PM_MCMCexperiment(noSims = X, noPanels = 8, lastObs = 5, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})
SIS_PMexperiment9 = lapply(X = noSims,
                        function(X) {PM_MCMCexperiment(noSims = X, noPanels = 9, lastObs = 5, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})

SIS_PMexperiment10 = lapply(X = noSims,
                        function(X) {PM_MCMCexperiment(noSims = X, noPanels = 10, lastObs = 5, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})

SIS_PMexperiment11 = lapply(X = noSims,
                        function(X) {PM_MCMCexperiment(noSims = X, noPanels = 11, lastObs = 5, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})

SIS_PMexperiment12 = lapply(X = noSims,
                        function(X) {PM_MCMCexperiment(noSims = X, noPanels = 12, lastObs = 5, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})
SIS_PMexperiment13 = lapply(X = noSims,
                        function(X) {PM_MCMCexperiment(noSims = X, noPanels = 13, lastObs = 5, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})
SIS_PMexperiment14 = lapply(X = noSims,
                        function(X) {PM_MCMCexperiment(noSims = X, noPanels = 14, lastObs = 5, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})
SIS_PMexperiment15 = lapply(X = noSims,
                        function(X) {PM_MCMCexperiment(noSims = X, noPanels = 15, lastObs = 5, lambda0 = 1e-3,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})

save.image("SIS_PMnoSimsExperiment.Rdata")

rm(list = ls())

# ==== fsMCMC SIS Experiment ====

noIts = 100000
adapt = TRUE
N = 200
noPanels = 5:15

fsMCMCSISexperiment = lapply(X = noPanels,
                                    function(X) {SISfsMCMCexperiment(noPanels = X, lastObs = 5, noDraws = 1000, lambda0 = 1e-3,
                                                                      V0 = diag(c(1/N, 1)), blockSize = 100,  adapt = T, noIts = noIts)})

save.image("SISfsMCMCnoPanelsExperiment.Rdata")

rm(list = ls())

ls()[!(ls() == "results")]
