source("HMC.R")

library(dplyr)
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
n = 120 # the number of simulated trajectories 

t = 3 

mu = rep(0, k) 
vars = rep(1, k)^2 

allses <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("mar", "se")) 

for (epsilon in seq(0.85, 1.25, by = 0.1)) {

  print(epsilon) 

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
  mar = median(ars) 

  for (se in ses) {
    allses <- add_row(allses, mar = round(mar, 3), se = se) 
  }
} 

# plot but first remove what boxplot would think is an outlier! 
boxplot(se~mar, allses, outline = FALSE) 
