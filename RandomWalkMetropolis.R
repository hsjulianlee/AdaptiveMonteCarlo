RWMetro = function(target, trajLength, x, COV, burnin=0)
{
  require(MASS)
  trajectory = rep(0,trajLength)
  trajectory[1] = x
  
  nAccepted = 0
  nRejected = 0
  
  for(t in 1:(trajLength-1)){
    currentPosition = trajectory[t]
    
    proposedJump = mvrnorm(n=1, currentPosition, COV)
    probAccept = min(1, target(proposedJump)/target(currentPosition))
    
    if(runif(1)<probAccept){
      trajectory[t+1] = proposedJump
      
      if(t > burnin){
        nAccepted = nAccepted + 1
      }
    } else{
      trajectory[t] = currentPosition
      
      if(t > burnin){
        nRejected = nRejected + 1
      }
    }
  }
  
  acceptedTraj = trajectory[(burnin+1):length(trajectory)]
  AcceptanceRate = signif(nAccepted/length(acceptedTraj),3)
  
}