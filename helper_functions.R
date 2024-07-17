## == Dependencies == ##
library("glmnet", "MASS", "Matrix")
# Notes:
#   glmnet is used for LASSO estimation
#   MASS is used for Moore-Penrose pseudoinverse (ginv)
#   Matrix is used for rankMatrix

## == DATA PROCESSING FUNCTIONS == ##
# Unpack data and construct predictors and response
# INPUTS
#   data:   (q + n) x p matrix of data (n = effective sample size)
#   q:      autoregressive order; default is 1 (i.e. VAR(1) model)
# OUTPUT: list with the following components:
#   x:      n x pq matrix of predictors
#   y:      n x p matrix of responses
unpack <- function(data, q = 1) {
  nplusq <- nrow(data)  # n plus q
  p <- ncol(data)       # p = dim of output
  n <- nplusq - q       # n = effective sample size
  pq <- p * q           # pq = total number of variables
  y <- data[(q + 1):(q + n), ]  # response
  x <- matrix(NA, n, pq)        # predictors
  for (ell in 1:q) {
    # where to store ellth lag --v
    block_ell <- ((ell - 1) * p + 1):(ell * p)
    # periods corresponding to ellth lag --v
    lag_ell <- (q + 1 - ell):(q + n - ell)
    # store 1st lags first, 2nd second,... --v
    x[, block_ell] <- data[lag_ell, ]
  }
  return(list(x = x, y = y))
}

## == LEAST SQUARES FUNCTIONS == ##

# Find solution to minimize ||y - Xb||_{\ell_2} over b
# INPUTS
#   x:          n x p matrix of regressors
#   y:          n x 1 vector of responses
# OUTPUTs fit: list with the following components
#   sol:        p x 1 vector of OLS estimates
#   full_rank:  logical; if TRUE, regressors of full rank
ls_sol <- function(x, y) {
  full_rank <- ncol(x) == Matrix::rankMatrix(x)[1] # check for full column rank
  if (full_rank == TRUE) { # if full rank, produce unique (OLS) solution
    sol <- solve(crossprod(x), crossprod(x, y))
  } else { # o/w produce solution based on pseudoinverse
    sol <- MASS::ginv(x) %*% y # ginv is Moore-Penrose pseudoinverse
    # TODO: Find more clever way to handle rank-deficient case
  }
  fit <- list(sol = sol, full_rank = full_rank)
  return(fit)
}

# Least squares refitting with multiple responses and same regressors
# INPUTS
#   x: n x pq matrix of predictors
#   y: n x p matrix of responses
#   that: p x pq matrix of estimates
# OUTPUTs refit: list with the following components
#   that: p x pq matrix of refitted estimates
#   full_rank: p x 1 logical vector; if TRUE, selected regressors of full rank
mult_refit <- function(x, y, that) {
  p <- ncol(y)                           # num responses
  sel <- that != 0                       # flag selections
  full_rank <- matrix(NA, p, 1)          # full rank flag
  for (i in 1:p) {
    shat_i <- sum(sel[i, ])              # num selected
    if (shat_i > 0) {                    # if something selected...
      xi <- as.matrix(x[, sel[i, ]])     # get active regressors
      # ^-- as.matrix used to keep the column dimension (also when = 1)
      fit_i <- ls_sol(xi, y[, i])        # refit
      that[i, sel[i, ]] <- fit_i$sol     # overwrite (keeping zeros)
      full_rank[i, 1] <- fit_i$full_rank # check rank
    } else {                             # if nothing selected...
      full_rank[i, 1] <- TRUE            # rank considered "full"
    }
  }
  refit <- list(that = that, full_rank = full_rank)
  return(refit)
}

## == LASSO FUNCTIONS == ##

# Multiple responses weighted LASSO using same regressors
# INPUTS
#   x: n x pq matrix of predictors
#   y: n x p matrix of responses
#   lambda_glmnet: penalty level in eyes of glmnet
#   upsilon: p x pq matrix of penalty loadings
#   full_path: logical; if TRUE, use full path of lambda values
#   tol_glmnet: convergence tolerance for glmnet
# OUTPUT
#   that: p x pq matrix of estimates
mult_lasso <- function(x, y, lambda_glmnet, upsilon = NULL,
                       full_path = FALSE, tol_glmnet = 1e-4) {
  if (missing(x) || missing(y) || missing(lambda_glmnet)) {
    stop("Not enough input arguments.")
  }
  pq <- ncol(x)  # columns in x = pq = total number of regressors
  p <- ncol(y)   # q = autoregressive order, p = dim(output)
  if (nrow(x) != nrow(y)) {
    stop("Number of observations (rows in x and y) do not match.")
  }
  if (is.null(upsilon)) {
    upsilon <- matrix(1, p, pq) # default to unit loadings
  }
  that <- matrix(NA, p, pq) # placeholder for estimates
  if (full_path == FALSE) { # if full path *not* requested (= default)...
    for (i in 1:p) {        # use glmnet w/ single lambda
      xtilde <- sweep(x, 2, upsilon[i, ], FUN = "/")  # rescale using loadings
      fit_i <- glmnet(xtilde, y[, i], family = "gaussian",
                      standardize = FALSE, intercept = FALSE,
                      lambda = lambda_glmnet, thresh = tol_glmnet)
      that[i, ] <- coef(fit_i)[-1] / upsilon[i, ] # original scaling
    }
  } else {            # if full path requested (not default)...
    for (i in 1:p) {  # use glmnet's lambda sequence and interpolate
      xtilde <- sweep(x, 2, upsilon[i, ], FUN = "/")  # rescale using loadings
      fit_i <- glmnet(xtilde, y[, i], family = "gaussian",
                      standardize = FALSE, intercept = FALSE,
                      thresh = tol_glmnet)
      # ^-- glmnet constructs its own lambda sequence
      that_i_temp <- coef(fit_i, s = lambda_glmnet)[-1] # interpolate
      that[i, ] <- that_i_temp / upsilon[i, ]           # original scale
    }
  }
  return(that) # return as sparse matrix
}

## == INFORMATION CRITERIA FUNCTIONS == ##

# Function to calculate the information criterion using glmnet
# INPUTS:
#   x:          matrix of predictors (n x p)
#   y:          vector of response (n x 1)
#   criteria:   character vector specifying the information criteria to use
#               (default is c("aic", "bic", "hqic"))
#   ...:        additional arguments to pass to glmnet
# OUTPUT: list with the following components:
#   betas: matrix of betas, with columns named after the criteria
#   lambdas: vector of lambdas that minimize each criterion
ic_lasso <- function(x, y, criteria = c("aic", "bic", "hqic"), ...) {
  criteria <- match.arg(criteria, several.ok = TRUE)
  # Fit the model using glmnet
  fit <- glmnet(x, y, family = "gaussian",
                intercept = FALSE, standardize = FALSE, ...)
  # Note: Intercepts are handled via prior demeaning
  df <- fit$df                          # degrees of freedom w/o intercept
  n_crit <- length(criteria)            # num criteria
  intrs <- numeric(n_crit)              # placeholder intercepts
  names(intrs) <- criteria              # name intercepts
  betas <- matrix(NA, ncol(x), n_crit)  # placeholder slopes
  colnames(betas) <- criteria           # name slopes
  lambdas <- numeric(n_crit)            # placeholder penalties
  names(lambdas) <- criteria            # name penalties
  n <- nrow(x)                          # num observations
  # Calculate mean squared errors and fit minimum information criterion
  res <- y - predict(fit, x, s = fit$lambda)  # residuals n x num(lambda)
  mse <- colMeans(res^2)                      # mse for each lambda candidate
  for (crit in criteria) {
    # Calculate the chosen information criterion
    if (crit == "aic") {
      ic <- n * log(mse) + 2 * df
    } else if (crit == "bic") {
      ic <- n * log(mse) + log(n) * df
    } else if (crit == "hqic") {
      ic <- n * log(mse) + 2 * log(log(n)) * df
    } else {
      stop("Unknown information criterion")
    }
    # Find the minimum information criterion, corresponding lambda and estimates
    id_min <- which.min(ic)                         # identify minimizer
    lambda_ic <- fit$lambda[id_min]                 # minimizer
    beta <- as.matrix(coef(fit, s = lambda_ic)[-1]) # resulting slopes
    betas[, crit] <- beta           # store slopes
    lambdas[crit] <- lambda_ic      # store penalty (in eyes of glmnet)
  }
  return(list(betas = betas, lambdas = lambdas))
}

## == SQUARE-ROOT LASSO FUNCTIONS == ##
#' Square-Root Lasso Estimator
#'
#' This function computes the Square-Root Lasso estimator, which solves
#' \deqn{min_{\beta}  \sqrt{qhat(\beta)} + \frac{\lambda}{n} \sum_{j} |\beta_{j}|} # nolint: line_length_linter.
#'
#' @param x Design matrix, n x p (without a column of 1s)
#' @param y Response vector, n x 1
#' @param lambda Penalty parameter, scalar
#' @param max_iter Maximum number of iterations (default 100)
#' @param opt_tol_norm Stopping tolerance for relative change in the Euclidean
#' norm of \eqn{\beta} (default 1e-4)
#' @return A vector of Square-Root Lasso estimates, p x 1
#' @details The function iteratively updates the \eqn{\beta} estimates using a
#' coordinate descent approach until the relative change in the Euclidean norm
#' of \eqn{\beta} falls below the specified tolerance or the maximum number of
#' iterations is reached.
#' @note The intercept is handled separately via prior demeaning.
#' @examples
#' # Example usage:
#' n <- 100
#' p <- 10
#' set.seed(42)
#' x <- matrix(rnorm(n * p), n, p)
#' y <- rnorm(n)
#' lambda <- 0.1
#' beta_estimates <- sqrt_lasso(x, y, lambda)
#' @export
sqrt_lasso <- function(x, y, lambda, max_iter = 100, opt_tol_norm = 1e-4) {
  # Compute necessary dimensions
  n <- nrow(x)
  p <- ncol(x)
  # Ridge matrix for the initial beta estimate
  ridge_matrix <- Matrix::Diagonal(p, x = lambda)
  beta <- solve(crossprod(x) + ridge_matrix, crossprod(x, y))
  beta <- as.matrix(beta)
  # Precompute quantities for the algorithm
  iter <- 0
  xx <- crossprod(x) / n
  xy <- crossprod(x, y) / n
  xx_diag <- diag(xx)
  x_beta <- x %*% beta
  error <- y - x_beta
  qhat <- mean(error^2)
  # Small constant to avoid division by zero
  epsilon <- .Machine$double.eps
  while (iter < max_iter) {
    iter <- iter + 1
    beta_old <- beta
    # Coordinate descent update for each beta_j
    for (j in 1:p) {
      s0 <- xx[j, ] %*% beta - xx_diag[j] * beta[j] - xy[j]
      beta_j_old <- beta[j]
      if (n^2 < (lambda^2) / xx_diag[j]) {
        beta[j] <- 0
      } else if (s0 > (lambda / n) * sqrt(qhat)) {
        beta[j] <- ((lambda / sqrt(n^2 - lambda^2 / xx_diag[j])) *
                      sqrt(max(qhat - (s0^2 / xx_diag[j]), 0)) - s0) /
          xx_diag[j]
      } else if (s0 < -(lambda / n) * sqrt(qhat)) {
        beta[j] <- (-(lambda / sqrt(n^2 - lambda^2 / xx_diag[j])) *
                      sqrt(max(qhat - (s0^2 / xx_diag[j]), 0)) - s0) /
          xx_diag[j]
      } else {
        beta[j] <- 0
      }
      if (beta[j] != beta_j_old) { # Update mean-square error if change
        x_beta <- x_beta + x[, j] * (beta[j] - beta_j_old)
        error <- y - x_beta
        qhat <- mean(error^2)
      }
    }
    # Check for convergence using the relative change in the Euclidean norm
    norm_diff <- norm(beta - beta_old, type = "2")          # numerator
    norm_beta_old <- norm(beta_old, type = "2") + epsilon   # denominator
    if (norm_diff / norm_beta_old < opt_tol_norm) {         # check tolerance
      break                                                 # break if met
    }
  }
  return(beta)
}

# Function which calculates the equation-by-equation square-root LASSO estimates
# for a VAR model, akin to mult_lasso
# INPUTS:
#   x:          n x pq matrix of predictors
#   y:          n x p matrix of responses
#   lambda:     penalty level (in eyes of sqrt_lasso)
#   upsilon:    p x pq matrix of penalty loadings
#   max_iter:   maximum number of iterations (default 10000)
#   print_out:  boolean; if TRUE, print iteration details (default FALSE)
#   opt_tol_norm: stopping tolerance for ||beta - beta_old|| (default 1e-6)
# OUTPUT: p x pq matrix of estimates
mult_sqrt_lasso <- function(x, y, lambda, upsilon = NULL,
                            max_iter = 10000, print_out = FALSE,
                            opt_tol_norm = 1e-6) {
  if (missing(x) || missing(y) || missing(lambda)) {
    stop("Not enough input arguments.")
  }
  pq <- ncol(x)  # columns in x = pq = total number of regressors
  p <- ncol(y)   # q = autoregressive order, p = dim(output)
  if (nrow(x) != nrow(y)) {
    stop("Number of observations (rows in x and y) do not match.")
  }
  if (is.null(upsilon)) {
    upsilon <- matrix(1, p, pq) # default to unit loadings
  }
  that <- matrix(NA, p, pq) # placeholder for estimates
  for (i in 1:p) {  # use sqrt_lasso
    xtilde <- sweep(x, 2, upsilon[i, ], FUN = "/")  # rescale using loadings
    that_i_temp <- sqrt_lasso(xtilde, y[, i], lambda, 
                              max_iter = max_iter,
                              opt_tol_norm = opt_tol_norm)
    that[i, ] <- that_i_temp / upsilon[i, ]         # original scale
  }
  return(that) # return as matrix
}