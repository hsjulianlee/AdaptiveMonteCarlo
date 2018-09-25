source("HMC.R")

# HMC run
L = 150 # the number of steps in each proposal (T-segment) 
N = 1000 # the number of T-segments 
k = 100 # the dimension of vectors (d in the paper) 

Sigma = diag(seq(0.01, k * 0.01, by = 0.01)^2)
mu = rep(0, k)
q = mvrnorm(n = 1, mu, Sigma)

HMC_sample = matrix(nrow = N, ncol = k)
epsilon = runif(1, 0.0104, 0.0156)
HMC_sample[1, ] = HMC(U, grad_U, epsilon, L, q)
for (i in 2:N) {
    epsilon = runif(1, 0.0104, 0.0156)
    HMC_sample[i, ] = HMC(U, grad_U, epsilon, L, HMC_sample[i-1, ])
}

# Figure 5.6
plot(HMC_sample[, k], pch = 16, ylim = c(-3, 3),
     las = 1,
     main = "Hamiltonian Monte Carlo",
     xlab = "Iteration",
     ylab = "Last position coordinate")

# Figure 5.7 means
means <- apply(HMC_sample, 2, mean)
plot(seq(0.01, k * 0.01, by = 0.01), means, pch = 16, ylim = c(-0.7, 0.7),
     las = 1,
     main = "Hamiltonian Monte Carlo",
     xlab = "Standard deviation of coordinate",
     ylab = "Sample mean of coordinate")
abline(0, 0)

# Figure 5.7 standard deviations
sdevs <- apply(HMC_sample, 2, sd)
plot(seq(0.01, k * 0.01, by = 0.01), sdevs, pch = 16, ylim = c(0, 1.2), xlim = c(0, 1),
     xaxs = "i", yaxs = "i", las = 1,
     main = "Hamiltonian Monte Carlo",
     xlab = "Standard deviation of coordinate",
     ylab = "Sample standard deviation of coordinate")
abline(0, 1)
