require(mvtnorm)

RWMetro = function(target, trajLength, x, COV, burnin=0)  {
  
  d = length(x) # dimension of the target distribution
  trajectory = matrix(0,trajLength, d)
  trajectory[1,] = x # initial point
  
  nAccepted = 0 #to count the number of accepted proposals
  nRejected = 0 #to count the number of rejected proposals
  
  currentPosition = x
  density = target(x)
  
  for(i in 1:(trajLength - 1)) {
    
    proposedJump = rmvnorm(1, currentPosition, COV)
    densityNew = target(proposedJump)
    
    if(density <= densityNew || runif(1) * density < densityNew) {
      trajectory[i + 1,] = proposedJump
      currentPosition = proposedJump
      density = densityNew
      if(i >= burnin) {
        nAccepted = nAccepted + 1
      }
      
    } else {
      trajectory[i + 1,] = currentPosition
      if(i >= burnin) {
        nRejected = nRejected + 1
      }
    }
    
  }
  
  AcceptanceRate = nAccepted / (trajLength - burnin)
  trajectory[(burnin + 1):trajLength,]
}
