library(MASS)

#' Single iteration of the HMC algorithm
#'
#' @param U a function which returns the potential energy given a value for q
#' @param grad_U a function which returns the vector of partial derivatives of U given q
#' @param epsilon the stepsize for leapfrog steps
#' @param L the number of leapfrog steps in the trajectory
#' @param current_q the current position that the trajectory starts from
#'
#' @return the position of the trajectory after L steps
#'
HMC = function (U, grad_U, epsilon, L, current_q)
{
  q = current_q
  p = rnorm(length(q),0,1) # independent standard normal variates
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  for (i in 1:L) {
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * grad_U(q)
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)) {
    return (q) # accept
  } else {
    return (current_q) # reject
  }
}

########################################################################################
### 
###       Sampling from a 100-dimensional multivariate Gaussian distribution in which
###       the variables are independent, with means of zero, and standard deviations
###       of 0.01, 0.02, ..., 0.99, 1.00.
### 
########################################################################################

#' A function to evaluate minus the log of the normal density, i.e. the "potential energy". 
#' The vector vars defining the standard deviations of coordinates should be defined outside the function. 
#'
#' @param q the point at which to estimate the function 
#'
#' @return minus log of density at q 
#'
U = function (q) 
{
  return ( 0.5 * sum(q^2 / vars) ) 
}

#' A function which returns the vector of partial derivatives of U at q 
#' The vector vars defining the standard deviations of coordinates should be defined outside the function. 
#'
#' @param q the point at which to estimate the function 
#'
#' @return the vector of partial derivatives of U at q 
#'
grad_U = function (q)
{
  return ( q / vars ) 
}
