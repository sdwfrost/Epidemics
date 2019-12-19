#' PM SIS experiment

noIts = 10000
adapt = TRUE
N = 200
noSims = 1:5

PMexperiment10 = lapply(X = noSims,
                      function(X) {PM_MCMCexperiment(noSims = X, noPanels = 10, lastObs = 5, lambda0 = 1e-4,
                                                      V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})

PMexperiment11 = lapply(X = noSims,
                      function(X) {PM_MCMCexperiment(noSims = X, noPanels = 11, lastObs = 5, lambda0 = 1e-4,
                                                     V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})

PMexperiment12 = lapply(X = noSims,
                        function(X) {PM_MCMCexperiment(noSims = X, noPanels = 12, lastObs = 5, lambda0 = 1e-4,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})
PMexperiment13 = lapply(X = noSims,
                        function(X) {PM_MCMCexperiment(noSims = X, noPanels = 13, lastObs = 5, lambda0 = 1e-4,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})
PMexperiment14 = lapply(X = noSims,
                        function(X) {PM_MCMCexperiment(noSims = X, noPanels = 14, lastObs = 5, lambda0 = 1e-4,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})
PMexperiment15 = lapply(X = noSims,
                        function(X) {PM_MCMCexperiment(noSims = X, noPanels = 15, lastObs = 5, lambda0 = 1e-4,
                                                       V0 = diag(c(1/N, 1)), adapt = T, noIts = noIts)})

save.image("PMnoSimsExperiment.Rdata")




