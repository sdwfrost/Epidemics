#'
#' Testing Household Kernel
#'

households = sample(1:3, size = 10, replace = TRUE)
household_matrix = 1*sapply(1:length(households), function(X) households[X] == households)

isSymmetric(household_matrix)


lower.tri(household_matrix) ==  t(upper.tri(household_matrix))

household_matrix[lower.tri(household_matrix, diag = T)] ==
  household_matrix[upper.tri(household_matrix, diag = T)]

kernel = Household_Kernel(households)


s = Matrix::Matrix(rbinom(100^2,size = 1, prob = 0.5), nrow = 100, ncol = 100)
s[lower.tri(s)] = Matrix::t(s)[lower.tri(s)]
Matrix::isSymmetric(s)
