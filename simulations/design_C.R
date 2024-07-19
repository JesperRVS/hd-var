# Simulate using essentially Kock & Callot [2015] Experimental Design D (which
# we here renamed) and retrieve coefficient matrix.

## IMPORT COEFFICIENT MATRIX FUNCTION
source("simulations/coef_mats.R") # for coef_mat_c
# coef_mat_c <- function(p) {0}
# insertSource("coef_mats.R", functions = "coef_mat_c")

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
sim_data_c <- function(n = 100, p = 4, sigma_eps = 0.1, nburn = 10000) {
  n_tot <- nburn + 1 + n
  yinit <- matrix(0, p, 1)      # p x 1 initial values (zeros)
  ylong <- matrix(NA, p, n_tot) # p x n_tot outcome matrix
  ylong[, 1] <- yinit           # initiate from all zeros
  normals <- matrix(rnorm(p * n_tot), p, n_tot) # p x n_tot std. normals
  eps <- sigma_eps * normals    # p x n_tot indep. innovations
  # Note: first epsilon not actually in use: kept for simple indexing.
  theta <- coef_mat_c(p)
  for (t in 2:n_tot) {
    ylong[, t] <- as.matrix(theta %*% ylong[, t - 1] + eps[, t])
  }
  # Note: as.matrix used to allow storage
  y0ton <- t(ylong[, (nburn + 1):n_tot])
  return(y0ton)
}