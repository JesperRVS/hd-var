# Simulate data using Kock & Callot's [2015 J Economet] Experimental Design A
# and retrieve coefficient matrix.

sim_data <- function(n, p, seed = NULL, r = NULL) {
  # INPUTS:
  #   n:          Effective sample size, integer
  #   p:          System dimension, divisible by four
  #   seed:       Seed for RNG
  #   r:          MC iter
  # OUTPUT:
  #   y0ton:      (1+n) x p outcome matrix
  if (missing(n) || missing(p)) {
    stop("Not enough input arguments.")
  }
  if (!is.null(seed) && !is.null(r)) {
    set.seed(seed + r)
  } else {
    set.seed(NULL)  # Reset RNG to default state
  }
  sigma_eps <- 0.1                          # std.dev of eps_0i [common]
  theta <- 0.5                              # AR(1) coef [common]
  sigma_y <- sigma_eps / sqrt(1 - theta^2)  # => std(y_0i) [common]
  yinit <- sigma_y * rnorm(p)               # sample from stationary dist.
  y0ton <- matrix(NA, nrow = p, ncol = 1 + n)
  y0ton[, 1] <- yinit                       # initiate from stationary dist.
  # independent eps_ti~N(0,sigma_eps^2) --v
  eps <- sigma_eps * matrix(rnorm(p * (1 + n)), nrow = p)
  # Note: 1st epsilon (i.e. for t=0) not actually in use
  theta <- coef_matrix(p)      # specify coef matrix (diagonal)
  for (t in 2:(1 + n)) {
    y0ton[, t] <- theta %*% y0ton[, t - 1] + eps[, t]
  }
  y0ton <- t(y0ton)  # return as (1+n) x p matrix
  return(y0ton)
}

# Coefficient matrix
coef_matrix <- function(p) {
  theta <- 0.5              # AR(1) coef [common]
  theta <- theta * diag(p)  # convert to diagonal matrix
  return(theta)
}
