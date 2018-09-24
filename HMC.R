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

#' A function to evaluate minus the log of the density of the distribution to be sampled, i.e. the "potential energy".
#'
#' @param q the point at which to estimate the function 
#'
#' @return minus log of density at q 
#'
U = function (q)
{
  k = 100
  # Covariance matrix
  Sigma = diag(seq(0.01, 1, by = 0.01)^2)
  # Determinant of covariance matrix
  detSigma = prod(diag(Sigma))
  # Diagonal of inverse of covariance matrix
  invSigmad = diag(solve(Sigma))
  # Density of the distribution
  f = exp(-0.5 * sum(invSigmad * q^2)) / sqrt((2 * pi)^k * detSigma)
  # Potential energy
  return(-log(f))
}

#' A function which returns the vector of partial derivatives of U given q
#'
#' @param q the point at which to estimate the function 
#'
#' @return the vector of partial derivatives of U at q 
#'
grad_U = function (q)
{
  k = 100
  # Covariance matrix
  Sigma = diag(seq(0.01, 1, by = 0.01)^2)
  # Determinant of covariance matrix
  detSigma = prod(diag(Sigma))
  # Diagonal of inverse of covariance matrix
  invSigmad = diag(solve(Sigma))
  # Density of the distribution
  f = exp(-0.5 * sum(invSigmad * q^2)) / sqrt((2 * pi)^k * detSigma)
  # Partial derivatives of density of the distribution
  df = -invSigmad * q * exp(-0.5 * sum(invSigmad * q^2)) / sqrt((2 * pi)^k * detSigma)
  # Partial derivatives of U
  return(-df/f)
}

# HMC run
L = 150
N = 1000
k = 100

Sigma = diag(seq(0.01, 1, by = 0.01)^2)
mu = rep(0, k)
q = mvrnorm(n = 1, mu, Sigma)

HMC_sample = matrix(nrow = N, ncol = 100)
epsilon = runif(1, 0.0104, 0.0156)
HMC_sample[1, ] = HMC(U, grad_U, epsilon, L, q)
for (i in 2:N) {
  epsilon = runif(1, 0.0104, 0.0156)
  HMC_sample[i, ] = HMC(U, grad_U, epsilon, L, HMC_sample[i-1, ])
}

# Figure 5.6
plot(HMC_sample[, 100], pch = 16, ylim = c(-3, 3),
     las = 1,
     main = "Hamiltonian Monte Carlo",
     xlab = "Iteration",
     ylab = "Last position coordinate")

# Figure 5.7 means
means <- apply(HMC_sample, 2, mean)
plot(seq(0.01, 1, by = 0.01), means, pch = 16, ylim = c(-0.7, 0.7),
     las = 1,
     main = "Hamiltonian Monte Carlo",
     xlab = "Standard deviation of coordinate",
     ylab = "Sample mean of coordinate")
abline(0, 0)

# Figure 5.7 standard deviations
sdevs <- apply(HMC_sample, 2, sd)
plot(seq(0.01, 1, by = 0.01), sdevs, pch = 16, ylim = c(0, 1.2), xlim = c(0, 1),
     xaxs = "i", yaxs = "i", las = 1,
     main = "Hamiltonian Monte Carlo",
     xlab = "Standard deviation of coordinate",
     ylab = "Sample standard deviation of coordinate")
abline(0, 1)
