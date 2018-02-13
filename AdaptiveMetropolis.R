RWMetro = function(target, trajLength, x, beta, burnin=0)
{
  require(MASS) 
  trajectory = rep(0,trajLength) 
  trajectory[1] = x # initial point
  
  d = length(x) # the dimension of the target distribution
  
  nAccepted = 0 # to count the number of accepted proposals
  nRejected = 0 # to count the number of rejected proposals
  
  C0 = (0.1)^2 * diag(d)/d # covariance matrix for t <=2d
  mx = rep(0,d) # initial mean value of the trajectory
  
  
  for(t in 1:(trajLength-1)){
    
    if(t <= 2*d){
      currentPosition = trajectory[t]
      #the proposal distribution given at iteration t<=2d 
      proposedJump = mvrnorm(currentPosition, (0.1)^2*diag(x=1, d, d)/d)
      
      probAccept = min(1, target(proposedJump)/target(currentPosition))
      
      if(runif(1)<probAccept){
        trajectory[t+1] = proposedJump
        mxold = mx #mean up to time t
        mx = (mx*t + trajectory[t+1]) / (t+1) #updated mean
        C = (1/t)*(trajectory[1:(t+1)]%*%t(trajectory[1:(t+1)])
                   - (t+1)*mx%*%t(mx)) #empirical covariance
        
        if(t > burnin){
          nAccepted = nAccepted + 1
        }
      } else{
        trajectory[t+1] = currentPosition
        mxold = mx #mean up to time t
        mx = (mx*t + trajectory[t+1]) / (t+1) #updated mean
        C = (1/t)*(trajectory[1:(t+1)]%*%t(trajectory[1:(t+1)])
                   - (t+1)*mx%*%t(mx)) #empirical covariance
        if(t > burnin){
          nRejected = nRejected + 1
        }
      }
      
    } else{
      currentPosition = trajectory[t]
      
      #proposal distribution for t>2d
      proposedJump = (1-beta)*mvrnorm(currentPosition, (2.38)*2*C/d) +
        beta * mvrnorm(currentPosition, (0.1)^2*diag(d)/d)
      
      probAccept = min(1, target(proposedJump)/target(currentPosition))
      
      sd = (1-beta)*2.38^2 / d
      sde = beta * (0.1)^2 / d
      
      if(runif(1)<probAccept){
        trajectory[t+1] = proposedJump
        mxold = mx #mean up to time t
        mx = (mx*t + trajectory[t+1]) / (t+1) #updated mean
        C = t/(t+1)C + sd/(t+1) * ((t+1)*mxold %*% t(mxold) - 
                                 (t+2)*mx %*% t(mx) + 
                                   trajectory[1:(t+1)] %*% 
                                   t(trajectory[1:(t+1)] + 
                                       e * diag(d)))
        #empirical covariance based on (3) Haario.
        
        if(t > burnin){
          nAccepted = nAccepted + 1
        }
      } else{
        trajectory[t+1] = currentPosition
        mxold = mx #mean up to time t
        mx = (mx*t + trajectory[t+1]) / (t+1) #updated mean
        C = t/(t+1)C + sd/(t+1) * ((t+1)*mxold %*% t(mxold) - 
                                     (t+2)*mx %*% t(mx) + 
                                     trajectory[1:(t+1)] %*% 
                                     t(trajectory[1:(t+1)] + 
                                         e * diag(d)))
        #empirical covariance based on (3) Haario.
        
        if(t > burnin){
          nRejected = nRejected + 1
        }
      }
    }
    
    
  
  acceptedTraj = trajectory[(burnin+1):length(trajectory)]
  AcceptanceRate = signif(nAccepted/length(acceptedTraj),3)
  
}