library("Matrix")
library("mvtnorm")

#' Simulate data according to a given design
#'
#' @param n Effective sample size
#' @param p System dimension
#' @param design Design type
#' @param sigma_eps Standard deviation of eps_0,i (which are indep. gaussian)
#' @param nburn No. of burn-in periods
#' @return data (q + n) x p outcome matrix (after burn-in)
#' @export
#' @examples
#' sim_data_by_design(n = 100, p = 4, design = "Diagonal", sigma_eps = 0.1,
#'                   seed = 1234, r = 1, nburn = 10000)
#' sim_data_by_design(n = 100, p = 4, design = "Correlated", sigma_eps = 0.1,
#'                  seed = 1234, r = 1, nburn = 10000)
#' sim_data_by_design(n = 100, p = 4, design = "HeavyTailed", sigma_eps = 0.1,
#'                 seed = 1234, r = 1, nburn = 10000)
#' sim_data_by_design(n = 100, p = 4, design = "BlockDiag", sigma_eps = 0.1,
#'                seed = 1234, r = 1, nburn = 10000)
#' sim_data_by_design(n = 100, p = 4, design = "NearBand", sigma_eps = 0.1,
#'               seed = 1234, r = 1, nburn = 10000)
sim_data_by_design <- function(n = 100, p = 4, design, nburn = 10000) {
  switch(design,
    # Diagonal design w/ independent Gaussian innovations
    "Diagonal" = {
      data <- sim_data_diag(n = n, p = p, coef = 0.5, family = "gaussian",
                            sigma_eps = 0.1, rho = 0, nburn = nburn)
    },
    # Near-band design w/ independent Gaussian innovations
    "NearBand" = {
      data <- sim_data_toep(n = n, p = p, sigma_eps = 0.1, nburn = nburn)
    },
    # Block diagonal design w/ independent Gaussian innovations
    "BlockDiag" = {
      data <- sim_data_block(n = n, p = p, sigma_eps = 0.1, nburn = nburn)
    },
    # Diagonal design w/ correlated Gaussian innovations
    "Correlated" = {
      data <- sim_data_diag(n = n, p = p, coef = 0.5, family = "gaussian",
                            sigma_eps = 0.1, rho = 0.9, nburn = nburn)
    },
    # Diagonal design w/ correlated Student-t innovations
    "HeavyTailed" = {
      data <- sim_data_diag(n = n, p = p, coef = 0.5, family = "student",
                            sigma_eps = 0.1, rho = 0.9, df = 5, nburn = nburn)
    },
    # Diagonal design w/ heteroskedastic innovations
    "Heteroskedastic" = {
      data <- sim_data_hetero(n = n, p = p, sigma_eps = 0.1, nburn = nburn)
    },
    # Near-unity design w/ independent Gaussian innovations
    "NearUnity" = {
      data <- sim_data_diag(n = n, p = p, coef = 1 - (5 / n),
                            family = "gaussian", sigma_eps = 0.1, rho = 0,
                            nburn = nburn)
    },
    stop("Design not recognized.")
  ) # end switch
  return(data)
}

### == DESIGNS == ###

## == DIAGONAL DESIGN(S) == ##
# Simulate data using variations of Kock & Callot's [2015 J Economet]
# Experimental Design A, which has a diagonal coefficient matrix.
# INPUTS:
#   n:          Effective sample size, integer
#   p:          System dimension, divisible by four
#   sigma:      Std.dev of eps_0,i (if gaussian)
#   rho:        Correlation between eps_0,i and eps_0,j, i != j
#               (if "gaussian")
#   family:     Distribution of eps_0i, "gaussian" or "student"
#   df:         Degrees of freedom for t dist (only relevant if "student")
#   nburn:      No. of burn-in periods
# OUTPUT:
#   y0ton:      (1 + n) x p outcome matrix (after burn-in)
sim_data_diag <- function(n = 100, p = 4, coef = 0.5, family = "gaussian",
                          sigma_eps = 0.1, rho = 0, df = 3, nburn = 10000) {
  n_tot <- nburn + 1 + n          # total no. of periods
  # Simulate n_tot x p innovations
  eps <- sim_eps(n = n_tot, p = p, family = family,
                 sigma_eps = sigma_eps, rho = rho, df = df)
  # Note: first epsilon not actually in use: kept for simple indexing.
  eps <- t(eps)                   # p x n_tot innovations
  yinit <- matrix(0, p, 1)        # p x 1 initial values (zeros)
  ylong <- matrix(NA, p, n_tot)   # p x n_tot outcome matrix
  ylong[, 1] <- yinit             # initiate from all zeros
  theta <- coef_mat_diag(p = p, coef = coef) # coef matrix (diagonal, sparse)
  for (t in 2:n_tot) {
    ylong[, t] <- as.matrix(theta %*% ylong[, t - 1] + eps[, t])
  }
  y0ton <- t(ylong[, (nburn + 1):n_tot])  # (1 + n) x p after burn-in
  return(y0ton)
}

## == DIAGONAL DESIGN(S) HELPER FUNCTIONS == ##

# Coefficient matrix for the VAR(1) process
# INPUTS:
#   p:          System dimension, integer
#   coef:       VAR(1) coefficient, scalar
# OUTPUT:
#   theta:      p x p VAR(1) coefficient matrix (diagonal, sparse)
coef_mat_diag <- function(p, coef) {
  m <- Matrix::Diagonal(n = p, x = coef)
  return(m)
}

# Simulate either possibly correlated Gaussian/Student innovations
# INPUTS:
#   n:          Sample size, integer
#   p:          System dimension, integer
#   sigma_eps:  Std.dev of eps_0,i (if gaussian)
#   rho:        Correlation between eps_0,i and eps_0,j, i != j
#   family:     Distribution of eps_0i, "gaussian" or "student"
#   df:         Degrees of freedom for t dist (only relevant if "student")
# OUTPUT:
#   eps:        n x p innovation matrix
sim_eps <- function(n = 100, p = 4, family = "gaussian",
                    sigma_eps = 0.1, rho = 0, df = NULL) {
  cor_eps <- cor_mat(p, rho)                        # innovation covariance
  switch(family,
    "gaussian" = {
      normals <- matrix(rnorm(n * p), n, p)         # n x p indep. normals
      eps <- sigma_eps * normals %*% chol(cor_eps)  # n x p correlated normals
    },
    "student" = {
      students <- mvtnorm::rmvt(n, sigma = cor_eps, df = df) # n x p
      # ^-- draws from multivariate (central) t w/ df=nu, scale=sigma
      nrmlzr <- sqrt(df / (df - 2))                 # normalizing constant
      eps <- (sigma_eps / nrmlzr) * students        # n x p scaled mult. t
      # => cov(eps) = sigma_eps^2 * cor_eps
    },
    stop("Family not recognized.")
  )
  return(eps)
}

# Innovation correlation matrix, equi-correlation structure
# INPUTS:
#   p:          System dimension, integer
#   rho:        Level of correlation between eps_0,i and eps_0,j, i != j
# OUTPUT:
#   sigma:      p x p innovation correlation matrix
cor_mat <- function(p, rho) {
  sigma <- rho * matrix(1, p, p) + (1 - rho) * diag(p)
  return(sigma)
}

## == NEAR-BAND DESIGN == ##
# Simulate using essentially Kock & Callot [2015] Experimental Design D.
# INPUTS:
#   n:          Effective sample size, positive integer
#   p:          System dimension, positive integer
#   sigma_eps:  Std.dev of eps_0,i (which are indep. gaussian)
#   nburn:      No. of burn-in periods, integer
# OUTPUT:
#   y0ton: (1+n) x p outcome matrix (after burn-in)
sim_data_toep <- function(n = 100, p = 4, sigma_eps = 0.1, nburn = 10000) {
  n_tot <- nburn + 1 + n
  yinit <- matrix(0, p, 1)      # p x 1 initial values (zeros)
  ylong <- matrix(NA, p, n_tot) # p x n_tot outcome matrix
  ylong[, 1] <- yinit           # initiate from all zeros
  normals <- matrix(rnorm(p * n_tot), p, n_tot) # p x n_tot std. normals
  eps <- sigma_eps * normals    # p x n_tot indep. innovations
  # Note: first epsilon not actually in use: kept for simple indexing.
  theta <- coef_mat_toep(p)
  for (t in 2:n_tot) {
    ylong[, t] <- as.matrix(theta %*% ylong[, t - 1] + eps[, t])
  }
  # Note: as.matrix used to allow storage
  y0ton <- t(ylong[, (nburn + 1):n_tot])
  return(y0ton)
}

## == NEAR-BAND DESIGN HELPER FUNCTIONS == ##

# Coefficient matrix for the VAR(1) process
# INPUT:
#   p:          System dimension, positive integer
# OUTPUT:
#   theta:      p x p coefficient matrix
coef_mat_toep <- function(p) {
  a <- 0.4                          # parameter (hardcoded)
  theta_temp <- matrix(NA, p, p)    # p x p temporary matrix
  for (i in 1:(p - 1)) {            # fill in strict upper triangle
    for (j in 2:p) {
      theta_temp[i, j] <- ((-1)^(abs(i - j))) * (a^(abs(i - j) + 1))
    }
  }
  tri_up <- Matrix::triu(theta_temp, 1) # strict upper triangle
  di <- diag(x = a, nrow = p, ncol = p) # diagonal matrix
  theta <- di + tri_up + t(tri_up)      # symmetrize
  return(theta)
}

## == BLOCK-DIAGONAL DESIGN == ##
# Simulate using essentially Kock & Callot [2015] Experimental Design B.
# Note: We use blocks of four instead of the five in Kock & Callot.
# INPUTS:
#   n:          Effective sample size, integer
#   p:          System dimension, divisible by four.
#   sigma_eps:  Std.dev of eps_0,i (which are indep. gaussian)
#   seed:       Seed for RNG
#   r:          MC iter
#   nburn:      No. of burn-in periods
# OUTPUT:
#   yminus3ton: (4 + n) x p outcome matrix (after burn-in)
sim_data_block <- function(n = 100, p = 4, sigma_eps = 0.1, nburn = 10000) {
  q <- 4                          # number of lags (hardcoded)
  n_tot <- nburn + q + n          # total no. of periods
  yinit <- matrix(0, p, q)        # p x q initial values (zeros)
  ylong <- matrix(NA, p, n_tot)   # p x n_tot outcome matrix
  ylong[, 1:q] <- yinit           # initiate from all zeros
  normals <- matrix(rnorm(p * n_tot), p, n_tot) # p x n_tot indep. normals
  eps <- sigma_eps * normals      # p x n_tot indep. Gaussian innovations
  # Note: first 4 epsilons not actually in use: kept for simple indexing.
  theta <- coef_mat_block(p)                # p x 4p coef matrix
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

## == BLOCK-DIAGONAL HELPER FUNCTIONS == ##
# Coefficient matrix for the VAR(4) process
# INPUT:
#   p:          System dimension, divisible by four
# OUTPUT:
#   theta:      p x 4p coefficient matrix
coef_mat_block <- function(p) {
  if (p %% 4 != 0) {
    warning("p should be divisible by four.")
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

## == HETEROSKEDASTIC DESIGN == ##
# Simulate from diagonal designs with heteroskedasticity in the (Gaussian)
# innovations.
sim_data_hetero <- function(n = 100, p = 4, coef = 0.5,
                            sigma_eps = 0.1, h = 1.5, nburn = 10000) {
  n_tot <- nburn + 1 + n
  yinit <- matrix(0, p, 1)      # p x 1 initial values (zeros)
  ylong <- matrix(NA, p, n_tot) # p x n_tot outcome matrix
  ylong[, 1] <- yinit           # initiate from all zeros
  etas <- matrix(rnorm(p * n_tot), p, n_tot) # p x n_tot std. normals
  eps <- matrix(NA, p, n_tot)    # p x n_tot indep. innovations
  # Note: first epsilon not actually in use, but kept for simple indexing.
  theta <- coef_mat_diag(p = p, coef = coef)
  for (t in 2:n_tot) {
    stds_etalag <- stds_eps_etalag(etalag = etas[, t - 1], h = h)
    eps[, t] <- sigma_eps * stds_etalag * etas[, t]
    ylong[, t] <- as.matrix(theta %*% ylong[, t - 1] + eps[, t])
  }
  # Note: as.matrix used to allow storage
  y0ton <- t(ylong[, (nburn + 1):n_tot])
  return(y0ton)
}

# == HETEROSKEDASTIC DESIGN HELPER FUNCTIONS == #
stds_eps_etalag <- function(etalag, h) {
  p <- length(etalag)
  stds_etalag <- numeric(p)
  for (i in 1:(p - 1)) {
    stds_etalag[i] <-
      exp(- h * abs(etalag[i]) + h * abs(etalag[i + 1]))
      # min(exp(- h * abs(etalag[i]) + h * abs(etalag[i + 1])), 100)
  }
  stds_etalag[p] <- exp(- h * abs(etalag[p]) + h * abs(etalag[1]))
    # min(exp(- h * abs(etalag[p]) + h * abs(etalag[1])), 100)
  return(stds_etalag)
}