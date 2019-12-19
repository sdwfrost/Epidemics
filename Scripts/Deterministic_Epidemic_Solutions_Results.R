#'
#' Deterministic Epidemic Solution Results
#'

Logistic_y_t(N = 100, y_0 = 5, beta = 0.001, t = 0.01)



LE = Logistic_Epidemic(N = 100, y_0 = 5, beta = 0.01, time_res = 0.001, t_lim = c(0,100))

Logistic_end(N = 100, y_0 = 5, beta = 0.01)

N = 1:1000
beta = 1/N


T_1_order = log(N)/(beta*N)

T_1 = sapply(X = N, function(X) Logistic_end(X, y_0 = 1, beta = 1/X))

plot(N, T_1, type = 'l', col = "darkgreen")
abline(v = 50, lty = 2)


Logistic_EC(N = 100, y_0 = 5, beta = 0.01, time_res = 0.001, t_lim = c(0,100))
x_t = 1:99
x_0 = 99
ln_x = log(x_t/x_0)

plot(x_t, ln_x)

object.size(library(coda))



