#'
#' Epidemic Gillespie w/ Kernel Test
#'
#'


# Set seed then carry out Epidemic Gillespie

# Set same seed then do Kernel Epidemic Gillespie with full contact matrix
# (Everyone in contact with each other)

R_0 = 1.5
gamma = 0.5
psi = R_0*gamma
N = 10^(2:5)
beta = psi/N

#' Which Population
i = 1

# ==== Without Kernel ====
susceptibles_left = c()
for(j in 1:2000){
  print(j)
  Test_Epidemic_norm = Epidemic_Gillespie(N[i], a = 5, gamma, beta[i])
  susceptibles_left[j] = min(Test_Epidemic_norm[,2], na.rm = T)
}

# ==== With Kernel ====

kernel = Contact_Kernel(matrix(rbinom(N[i], size = 1, prob = 0.99),
                               nrow = N[i], ncol = N[i]))

susceptibles_left_split = c()
for(j in 1:2000){
  print(j)
  Test_Epidemic_kernel = Kernel_Epidemic_Gillespie(N[i], a = 5, gamma, beta[i],
                                                   kernel)
  susceptibles_left_split[j] = min(Test_Epidemic_kernel[,2], na.rm = T)
}

susceptibles_left_det = c()
for(j in 1:2000){
  print(j)
  U = runif(2*N[i])
  E = rexp(2*N[i])
  Det_Epidemic = Deterministic_Gillespie1(N[i], initial_infective = 5, beta[1], gamma, E, U,
                                          T_obs = c(0,10), k = 1, store = TRUE)
  susceptibles_left_det[j] = min(Det_Epidemic$sim_data[,2], na.rm = T)
}

boxplot(susceptibles_left, susceptibles_left_det)

plot(final_size)
plot(susceptibles_left, ylim = c(0, 100))


sum(samples == 0)
sum(split_samples == 0)
X = 95
Y = 5
B = kernel(beta[1])
reduced_B = B[1:5, 1:95]
inf_rate = beta[1]*X*Y
rem_rate = gamma*Y
rate_next_event = rem_rate + inf_rate
#' 1 Removal
#' 0 Infection




for(j in 1:200000){
  samples[j] = sample(c(0,1), size = 1, prob = c(X*Y*beta[1], Y*gamma)/rate_next_event)
}

for(j in 1:200000){
  split_samples[j] = 1*(Heterogeneous_Event(integer(0), rem_rate) == 0)
}

samples = c()
split_samples = c()

X= 95
for(j in 1:10000){
  #set.seed(j)
  samples[j] = sample(c(0,1), size = 1, prob = c(X*Y*beta[1], Y*gamma)/rate_next_event)

  #set.seed(j)
  split_samples[j] = 1*(Heterogeneous_Event(reduced_B, rem_rate) == 0)

}


sum(samples == 1)/length(samples)
sum(split_samples == 1)/length(samples)

final_size = c()
final_size_split = c()
first_event = c()
for(j in 10001:20000){
  Split_Epidemic = Kernel_Epidemic_Gillespie(N[1], a = 5, gamma, beta[1], kernel)
  final_size_split[j] = max(Split_Epidemic[,4], na.rm = T) - 5
  Epidemic = Epidemic_Gillespie(N[1], 5, gamma, beta[1])
  final_size[j] = max(Epidemic[,4], na.rm = T) - 5
  #first_event[j] = Split_Epidemic[2,5]
}
sum(final_size )
sum(first_event)
plot(final_size_split)

results = data.frame(c(rep("A", 1000), rep("B", 1000)), c(final_size, final_size_split))

plot.coords(xlabel = c(rep("A", 1000), rep("B", 1000)), c(final_size, final_size_split))

boxplot(final_size, final_size_split)
sum(final_size == 0)
