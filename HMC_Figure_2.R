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
    accepted <- 0 
    repeat {
      qnew <- HMC(U, grad_U, epsilon, L, q) 
      if (any(qnew != q)) accepted <- accepted + 1 
      q <- qnew 
      HMC_sample <- rbind(HMC_sample, q) 
    }
  }, substitute = FALSE, timeout = t, onTimeout = "silent") 
  nT <- length(HMC_sample[,1]) - 1 
  return ( list( nT = nT, sample = HMC_sample , ar = accepted / nT ) ) 
}

# HMC run
L = 100 # the number of steps in each proposal (T-segment) 
k = 100 # the dimension of vectors (d in the paper) 
n = 10 # the number of simulated trajectories 

t = 1 

mu = rep(0, k) 
vars = rep(1, k)^2 

epsilon = 0.7 

nTs <- rep(0, n) 
ses <- rep(0, n) 
ars <- rep(0, n) 
for (i in 1:n) {
  q = mvrnorm(n = 1, mu, diag(vars)) 
  results <- runtsec(t, q, epsilon) 
  nTs[i] <- results$nT 
  ses[i] <- mean(results$sample[,1])^2 
  ars[i] <- results$ar 
}

print(median(nTs)) 
print(median(ars)) 
boxplot(ses, xlab = as.character(median(ars))) 
