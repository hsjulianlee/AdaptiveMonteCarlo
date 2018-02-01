RWMetro = function(target, trajLength, x, beta, burnin=0)
{
  require(MASS) 
  trajectory = rep(0,trajLength) 
  trajectory[1] = x # initial point
  
  d = length(x) # the dimension of the target distribution
  
  nAccepted = 0 # to count the number of accepted proposals
  nRejected = 0 # to count the number of rejected proposals
  
  for(t in 1:(trajLength-1)){
    
    if(t <= 2*d){
      currentPosition = trajectory[t]
      #the proposal distribution given at iteration t<=2d 
      proposedJump = mvrnorm(currentPosition, (0.1)^2*diag(x=1, d, d)/d)
      
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
      
    } else{
      currentPosition = trajectory[t]
      
      #proposal distribution for t>2d
      proposedJump = (1-beta)*mvrnorm(currentPosition, (2.38)*2)
      
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