require(mvtnorm)

AMetro = function(target, trajLength, x, beta, burnin = 0) {
  
  d = length(x) # the dimension of the target distribution
  
  trajectory = matrix(0, trajLength, d)
  trajectory[1,] = x # initial point
  
  nAccepted = 0 # to count the number of accepted proposals
  nRejected = 0 # to count the number of rejected proposals
  
  c1 <- (1 - beta)^2 * 2.38^2 / d
  C0 <- (0.1)^2 * diag(d) / d # covariance matrix for t <= 2d
  c2 <- beta^2 * C0
  d2 = d * 2
  
  currentPosition <- x
  density <- target(x)
  sumx <- rep(0, d)
  
  for (i in 1:(trajLength - 1)) {
    sumx_1 <- sumx
    sumx <- sumx + currentPosition
    if (i < d2) {
      # the proposal distribution given at iteration i <= 2d
      proposedJump = rmvnorm(1, currentPosition, C0)
    } else {
      # the proposal distribution for i > 2d
      if (i == d2) {
        Ci_1 <- (i - 1) * cov(trajectory[1:i,])
      }
      else {
        # C[i] = 1/(i-1) * (sum[j=1..i] xj xj^t - i * x_[i] x_[i]^t)
        # (i-1)C[i] = sum[j=1..i] xj xj^t - sumx[i] sumx[i]^t / i
        # (i-2)C[i-1] = sum[j=1..i-1] xj xj^t - sumx[i-1] sumx[i-1]^t / (i-1)
        # difference = xi xi^t - sumx[i] sumx[i]^t / i + sumx[i-1] sumx[i-1]^t / (i-1) 
        Ci_1 <- Ci_1 + (t(currentPosition) %*% currentPosition) - (t(sumx) %*% sumx) / i + (t(sumx_1) %*% sumx_1) / (i - 1)
      }
      proposedJump = rmvnorm(1, currentPosition, c1 * Ci_1 / (i - 1) + c2)
    }
    
    densityNew <- target(proposedJump)
    
    if (density <= densityNew || runif(1) * density < densityNew) {
      trajectory[i + 1,] = proposedJump
      currentPosition = proposedJump
      density <- densityNew
      if (i >= burnin) {
        nAccepted = nAccepted + 1
      }
    } else {
      trajectory[i + 1,] = currentPosition
      if (i >= burnin) {
        nRejected = nRejected + 1
      }
    }
  }
  
  AcceptanceRate = nAccepted / (trajLength - burnin)
  trajectory[(burnin + 1):trajLength,]
  
}
