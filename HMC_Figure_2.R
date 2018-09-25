source("HMC.R")

library(R.utils) 

#' Run HMC for t seconds. 
#' The vector vars defining the standard deviations of coordinates should be defined outside the function. 
#'
#' @param q the initial point 
#'
#' @return the sample where the number of rows is variable depending on how many samples were obtained in the given time. 
#'
runtsec = function(t, q, epsilon) {
  HMC_sample = matrix(q, nrow = 1) 
  withTimeout({ 
    repeat {
      q <- HMC(U, grad_U, epsilon, L, q) 
      HMC_sample <- rbind(HMC_sample, q) 
    }
  }, substitute = FALSE, timeout = t, onTimeout = "silent") 
  return ( HMC_sample )  
}

# HMC run
L = 100 # the number of steps in each proposal (T-segment) 
k = 100 # the dimension of vectors (d in the paper) 

t = 1 

mu = rep(0, k) 
vars = seq(0.01, k * 0.01, by = 0.01)^2 
q = mvrnorm(n = 1, mu, diag(vars)) 

epsilon = 0.013 

HMC_sample <- runtsec(t, q, epsilon) 

# Figure 5.7 means
means <- apply(HMC_sample, 2, mean)
plot(sqrt(vars), means, pch = 16, ylim = c(-0.7, 0.7),
     las = 1,
     main = "Hamiltonian Monte Carlo",
     xlab = "Standard deviation of coordinate",
     ylab = "Sample mean of coordinate")
abline(0, 0)

# Figure 5.7 standard deviations
sdevs <- apply(HMC_sample, 2, sd)
plot(sqrt(vars), sdevs, pch = 16, ylim = c(0, 1.2), xlim = c(0, 1),
     xaxs = "i", yaxs = "i", las = 1,
     main = "Hamiltonian Monte Carlo",
     xlab = "Standard deviation of coordinate",
     ylab = "Sample standard deviation of coordinate")
abline(0, 1)
