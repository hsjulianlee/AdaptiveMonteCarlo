AMetro = function(target, trajLength, x, beta, burnin = 0) {
  require(MASS) 
  
  d = length(x) # the dimension of the target distribution
  
  trajectory = matrix(0, trajLength, d)
  trajectory[1,] = x # initial point
  
  nAccepted = 0 # to count the number of accepted proposals
  nRejected = 0 # to count the number of rejected proposals
  
  C0 = (0.1)^2 * diag(d) / d # covariance matrix for t <= 2d
  
  currentPosition <- x
  density <- target(x)
  
  for (i in 1:(trajLength - 1)) {
    if (i <= 2*d) {
      # the proposal distribution given at iteration i <= 2d
      proposedJump = mvrnorm(1, currentPosition, C0)
    } else {
      # the proposal distribution for i > 2d
      C = cov(trajectory[1:i,])
      proposedJump = mvrnorm(1, currentPosition, (1 - beta)^2 * 2.38^2 * C / d + beta^2 * C0)
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
