source("AdaptiveMetropolis.R")

d = 100
N <- 10000
b <- 0.05

O = rep(0, d)
M = matrix(rnorm(d * d, 0, 1), ncol = d)
MM = M %*% t(M)

targetNMM = function(x, M = O, S = MM) {
  return(dmvnorm(x, M, S))
}

start_time <- Sys.time()
trajectory = AMetro(target = targetNMM, trajLength = N, x = O, beta = b, burnin = 0)
print(Sys.time() - start_time)

plot(1:length(trajectory[,1]), trajectory[,1], "l")
