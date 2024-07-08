# Simulate data using Kock & Callot's [2015 J Economet] Experimental Design A
# and retrieve coefficient matrix.

# TODO:
# - Draw directly from stationary distribution of VAR(1) process instead

## Dependencies
library("Matrix") # for triu (upper triangle of matrix)

## MAIN FUNCTION
sim_data <- function(n = 100, p = 4, theta = 0.5,
                     family = "gaussian", sigma_eps = 0.1, rho = 0, df = 3,
                     seed = NULL, r = NULL, nburn = 10000) {
  # INPUTS:
  #   n:          Effective sample size, integer
  #   p:          System dimension, divisible by four
  #   sigma:      Std.dev of eps_0,i (if gaussian)
  #   rho:        Correlation between eps_0,i and eps_0,j, i != j
  #               (if "gaussian")
  #   family:     Distribution of eps_0i, "gaussian" or "student"
  #   df:         Degrees of freedom for t dist (if "student")
  #   seed:       Seed for RNG
  #   r:          MC iter
  #   nburn:      No. of burn-in periods
  # OUTPUT:
  #   y0ton:      (1+n) x p outcome matrix
  set.seed(seed + r)                            # set seed for reproducibility
  n_tot <- nburn + 1 + n                        # total no. of periods
  # Simulate n_tot x p innovations
  eps <- sim_eps(n = n_tot, p = p, family = family,
                 sigma_eps = sigma_eps, rho = rho, df = df)
  # Note: first epsilon not actually in use: kept for simple indexing.
  eps <- t(eps)                                 # p x n_tot innovations
  yinit <- matrix(0, p, 1)                      # p x 1 initial values (zeros)
  ylong <- matrix(NA, p, n_tot)                 # p x n_tot outcome matrix
  ylong[, 1] <- yinit                           # initiate from all zeros
  theta <- coef_matrix(theta, p)    # specify coef matrix (diagonal)
  for (t in 2:n_tot) {
    ylong[, t] <- theta %*% ylong[, t - 1] + eps[, t] # VAR(1) process
  }
  y0ton <- ylong[, (nburn + 1):n_tot] # p x (1 + n) after burn-in
  y0ton <- t(y0ton)                   # return as (1 + n) x p matrix
  return(y0ton)
}

## HELPER FUNCTIONS

# Simulate either possibly correlated Gaussian innovations or independent
# Student-t innovations
sim_eps <- function(n = 100, p = 4, family = "gaussian",
                    sigma_eps = 0.1, rho = 0, df = NULL) {
  # INPUTS:
  #   n:          Sample size, integer
  #   p:          System dimension, integer
  #   sigma_eps:  Std.dev of eps_0,i (if gaussian)
  #   rho:        Correlation between eps_0,i and eps_0,j, i != j
  #               (if "gaussian")
  #   family:     Distribution of eps_0i, "gaussian" or "student"
  #   df:         Degrees of freedom for t dist (if "student")
  # OUTPUT:
  #   eps:        n x p innovation matrix
  switch(family,
    "gaussian" = {
      cor_eps <- cor_mat(p, rho)                    # innovation covariance
      normals <- matrix(rnorm(n * p), n, p)         # n x p indep. normals
      eps <- sigma_eps * normals %*% chol(cor_eps)  # n x p correlated normals
    },
    "student" = {
      nrmlzr <- sqrt(df / (df - 2))                 # normalizing constant
      students <- matrix(rt(n * p, df), n, p)       # n x p indep. t(df) dist
      eps <- (sigma_eps / nrmlzr) * students        # => var = sigma_eps^2
    }
  )
  return(eps)
}

# Innovation correlation matrix, Toeplitz structure
# INPUTS:
#   p:          System dimension, integer
#   rho:        Level of correlation between eps_0,i and eps_0,j, i != j
# OUTPUT:
#   sigma:      p x p innovation correlation matrix
cor_mat <- function(p, rho) {
  sigma_temp <- matrix(0, p, p)
  for (j in 1:(p - 1)) {
    for (k in (j + 1):p) {
      sigma_temp[j, k] <- rho^(abs(j - k))
    }
  }
  tri_up <- triu(sigma_temp, 1) # strict upper triangle (from Matrix)
  sigma <- diag(p) + tri_up + t(tri_up)
  return(sigma)
}

# Coefficient matrix
# INPUTS:
#   p:          System dimension, integer
# OUTPUT:
#   theta:      p x p VAR(1) coefficient matrix (diagonal)
coef_matrix <- function(theta, p) {
  m <- theta * diag(p)  # convert theta to diagonal matrix
  return(m)
}
