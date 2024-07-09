# Simulate using essentially Kock & Callot [2015] Experimental Design D (which
# we here renamed) and retrieve coefficient matrix.

## Dependencies
library("Matrix") # for triu and Diagonal

## MAIN FUNCTION

# INPUTS:
#   n:          Effective sample size, positive integer
#   p:          System dimension, positive integer
#   sigma_eps:  Std.dev of eps_0,i (which are indep. gaussian)
#   seed:       Seed for RNG
#   r:          MC iter, integer
#   nburn:      No. of burn-in periods, integer
# OUTPUT:
#   y0ton: (1+n) x p outcome matrix (after burn-in)
sim_data <- function(n = 100, p = 4, sigma_eps = 0.1,
                     seed = NULL, r = NULL, nburn = 10000) {
  set.seed(seed + r)
  n_tot <- nburn + 1 + n
  yinit <- matrix(0, p, 1)      # p x 1 initial values (zeros)
  ylong <- matrix(NA, p, n_tot) # p x n_tot outcome matrix
  ylong[, 1] <- yinit           # initiate from all zeros
  normals <- matrix(rnorm(p * n_tot), p, n_tot) # p x n_tot std. normals
  eps <- sigma_eps * normals    # p x n_tot indep. innovations
  # Note: first epsilon not actually in use: kept for simple indexing.
  theta <- coef_mat(p)
  for (t in 2:n_tot) {
    ylong[, t] <- as.matrix(theta %*% ylong[, t - 1] + eps[, t])
  }
  # Note: as.matrix used to allow storage
  y0ton <- t(ylong[, (nburn + 1):n_tot])
  return(y0ton)
}

## HELPER FUNCTION

# Coefficient matrix for the VAR(1) process
# INPUT:
#   p:          System dimension, positive integer
# OUTPUT:
#   theta:      p x p coefficient matrix
coef_mat <- function(p) {
  a <- 0.4                          # parameter (hardcoded)
  theta_temp <- matrix(NA, p, p)    # p x p temporary matrix
  for (i in 1:(p - 1)) {            # fill in strict upper triangle
    for (j in 2:p) {
      theta_temp[i, j] <- ((-1)^(abs(i - j))) * (a^(abs(i - j) + 1))
    }
  }
  tri_up <- Matrix::triu(theta_temp, 1) # strict upper triangle
  di <- Matrix::Diagonal(n = p, x = a)  # diagonal (sparse)
  theta <- di + tri_up + t(tri_up)      # symmetrize
  return(theta)
}