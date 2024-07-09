# Simulate using essentially Kock & Callot [2015] Experimental Design B and
# retrieve coefficient matrices (as one wide matrix). 
# Note: Since our dimensions (p) are all even, we use blocks of four instead of
# the five in Kock & Callot.

## Dependencies
library("Matrix") # for Diagonal

## MAIN FUNCTION

# INPUTS:
#   n:          Effective sample size, integer
#   p:          System dimension, divisible by four.
#   sigma_eps:  Std.dev of eps_0,i (which are indep. gaussian)
#   seed:       Seed for RNG
#   r:          MC iter
#   nburn:      No. of burn-in periods
# OUTPUT:
#   yminus3ton: (4 + n) x p outcome matrix (after burn-in)
sim_data <- function(n = 100, p = 4, sigma_eps = 0.1,
                     seed = NULL, r = NULL, nburn = 10000) {
  set.seed(seed + r)    # set seed for reproducibility
  q <- 4                # number of lags (hardcoded)
  n_tot <- nburn + q + n          # total no. of periods
  yinit <- matrix(0, p, q)        # p x q initial values (zeros)
  ylong <- matrix(NA, p, n_tot)   # p x n_tot outcome matrix
  ylong[, 1:q] <- yinit           # initiate from all zeros
  normals <- matrix(rnorm(p * n_tot), p, n_tot) # p x n_tot indep. normals
  eps <- sigma_eps * normals      # p x n_tot indep. gaussian innovations
  # Note: first 4 epsilons not actually in use: kept for simple indexing.
  theta <- coef_mat(p)                      # p x 4p coef matrix
  theta1 <- theta[, 1:p]                    # p x p block diagonal
  theta4 <- theta[, (3 * p + 1):(4 * p)]    # p x p block diagonal
  # Simulate VAR(4) process
  for (t in 5:n_tot) {
    ylong[, t] <- as.matrix(theta1 %*% ylong[, t - 1] +
                            theta4 %*% ylong[, t - 4] +
                            eps[, t])
  }
  # Note: as.matrix used to allow storage
  yminus3ton <- t(ylong[, (nburn + 1):n_tot]) # p x (q + n) after burn-in
  return(yminus3ton)
}

## HELPER FUNCTION
# Coefficient matrix for the VAR(4) process
# INPUT:
#   p:          System dimension, divisible by four
# OUTPUT:
#   theta:      p x 4p coefficient matrix
coef_mat <- function(p) {
  if (p %% 4 != 0) {
    warning("p must be divisible by four.")
  }
  pover4 <- p / 4                           # p/4 (num blocks)
  block1 <- matrix(0.15, 4, 4)              # block in theta1
  ipover4 <- Matrix::Diagonal(pover4)       # sparse identity I_{p/4}
  theta1 <- kronecker(ipover4, block1)      # p x p block diagonal
  block4 <- matrix(-0.1, 4, 4)              # block in theta4
  theta4 <- kronecker(ipover4, block4)      # p x p block diagonal
  theta23 <- matrix(0, p, 2 * p)            # zeros in theta2 and theta3
  theta <- cbind(theta1, theta23, theta4)   # p x 4p coefficient matrix
  return(theta)
}