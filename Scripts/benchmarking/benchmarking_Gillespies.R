#' Benchmarking purposed Gillespies against "God function" Gillespie

N = 200
beta = 0.02
gamma = 0.5
kernel = Contact_Kernel(matrix(1, nrow = N, ncol = N))
initialState = c(rep(1, N - 1), 2)
obsTimes = seq(0, 10, length = 6)

results = matrix(NA, nrow = 5, ncol = 2)

# ==== Homogeneous Population Level Gillespie ====

results[1, 1] = system.time(replicate(1000, homogeneousSIR_Gillespie(c(199, 1, 0), beta, gamma)))[3]

results[1, 2] = system.time(replicate(1000, SIR_Gillespie(N, beta = beta, gamma = gamma)))[3]

# ==== Homogeneous Individual Level Gillespie ====

results[2, 1] = system.time(replicate(1000, individualHomogeneousSIR_Gillespie(initialState, beta, gamma)))[3]
results[2, 2] = system.time(replicate(1000, SIR_Gillespie(N, beta = beta, gamma = gamma)))[3]

# ==== Heterogeneous kernel-based Gillespie ====

results[3, 1] = system.time(replicate(1000, heterogeneousSIR_Gillespie(initialState, beta, gamma, kernel)))[3]
results[3, 2] = system.time(replicate(1000, SIR_Gillespie(N, beta = beta, gamma = gamma, kernel = kernel)))[3]

# ==== Extracting Panel Data ====

#' No Kernel
results[4, 1] = system.time(replicate(1000,
                                      individualHomogeneousPanelDataSIR_Gillespie(initialState, beta, gamma, obsTimes)))[3]
results[4, 2] = system.time(replicate(1000, SIR_Gillespie(N, beta = beta, gamma = gamma,
                                                          obs_times = obsTimes, output = "panel")))[3]



#' Kernel
results[5, 1] = system.time(replicate(1000,
                                      panelDataSIR_Gillespie(initialState, beta, gamma, obsTimes, kernel)))[3]
results[5, 2] = system.time(replicate(1000, SIR_Gillespie(N, beta = beta, gamma = gamma, kernel = kernel,
                                                          obs_times = obsTimes, output = "panel")))[3]

# ==== Plot Results ====

barplot(t(results), beside = T, names.arg = c("Pop. Homogeneous", "Ind. Homogeneous",
                                              "Heterogeneous (K)", "Panel Data", "Panel Data (K)"),
        col = c("lightblue", "coral3"), ylab = "Time (s) (1000 reps)", main = "Run time of purposed Gillespie Algorithms vs. an All purpose Gillespie")
legend("topleft", legend = c("Purposed", "God"), fill = c("lightblue", "coral3"))

#' Profile SIR_Gillespie for Extracting Panel Data and Heterogenoeus Kernel-based

tmp = tempfile()
Rprof(tmp, interval = 0.01)
x = replicate(1000, SIR_Gillespie(N, beta = beta, gamma = gamma))
Rprof(NULL)
summaryRprof(tmp)$by.self
summaryRprof(tmp)$sampling.time

# amending data.frames is expensive compared to a vanilla matrix.
tmp = tempfile()
Rprof(tmp, interval = 0.01)
x = replicate(1000, panelDataSIR_Gillespie(initialState, beta, gamma, obsTimes,  kernel))
Rprof(NULL)
summaryRprof(tmp)$by.self
summaryRprof(tmp)$sampling.time


# ==== FINDINGS ====

# Within SIR_Gillespie:
#'1. amending data.frames is expensive compared to a vanilla matrix.

#' SIR_Gillespie vs. Purposed

#' 1. Similar performance for individual based methods.
#'
#' Purposed methods waste time calling sample.int twice.
#'
#' Where as SIR_Gillespie wastes time calling event.epidemics
#' which uses a quick method (I think), but creates overhead
#' through calling this user defined function a lot. (FURTHER TESTING)

#' Within Purposed:
#'  Faster to use .Internal(colSums()) to avoid unnecessary checks
#'  Anymore Internal functions that can be used?
#'
#' Cumsum and runif is faster than two samples
#' Find a cleverer way to generate list of individuals and rates
#' to sample from.

X = matrix(0, nrow = 200, ncol = 200)
X_nrow = dim(X)[1]
X_ncol = dim(X)[2]
x = .Internal(colSums(X, dim(X)[1], dim(X)[2], F))

system.time(replicate(400, .Internal(colSums(X, X_nrow, X_ncol, F))))
system.time(replicate(400, .Internal(colSums(X, dim(X)[1], dim(X)[2], F))))
system.time(replicate(400, colSums(X)))
#' This is not faster

system.time(replicate(100, colSums(X)))
system.time(replicate(100, colSums.real.matrix(X)))


p  = runif(10e3)
cs = cumsum(p)
U  = runif(1, min = 0, max = cs[length(cs)])











