## == Dependencies == ##
library("glmnet", "MASS", "Matrix")
# Notes:
#   glmnet is used for LASSO estimation
#   MASS is used for Moore-Penrose pseudoinverse (ginv)
#   Matrix is used for rankMatrix

## == DATA PROCESSING FUNCTIONS == ##
#' Unpacks time series data into predictors and response for VAR(q) model.
#'
#' @param data A (q + n) x p matrix of time series data, where n is the effective
#'             sample size.
#' @param q Autoregressive order; default is 1 (VAR(1) model).
#'
#' @return A list with the following components:
#'   \item{x}{An n x pq matrix of predictors, where pq = p * q.}
#'   \item{y}{An n x p matrix of responses.}
#'
#' @details
#' The function extracts predictors and responses from time series data for
#' fitting a VAR(q) model. It sets up lagged values of the data as predictors
#' and the subsequent period's data as responses.
#'
#' @examples
#' # Create example data
#' data <- matrix(1:12, ncol = 3)
#' # Unpack data assuming q = 1
#' result <- unpack(data, q = 1)
#' # Expected output for x and y
#' expected_x <- matrix(c(1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6), nrow = 4, byrow = TRUE)
#' expected_y <- matrix(c(5, 6, 7, 8), nrow = 4, byrow = TRUE)
#' # Check if the function output matches expected values
#' identical(result$x, expected_x)  # Should return TRUE
#' identical(result$y, expected_y)  # Should return TRUE
#'
#' @export
unpack <- function(data, q = 1) {
  if (missing(data) || is.null(data)) {
    stop("No data provided.")
  }
  nplusq <- nrow(data)  # n plus q
  p <- ncol(data)       # p = number of variables
  n <- nplusq - q       # n = effective sample size
  pq <- p * q           # pq = total number of variables
  y <- data[(q + 1):(q + n), ]  # responses
  # Initialize x with NA for subsequent checks
  x <- matrix(NA, n, pq)
  # Fill x with lagged responses
  for (ell in 1:q) {
    block_ell <- ((ell - 1) * p + 1):(ell * p)
    lag_ell <- (q + 1 - ell):(q + n - ell)
    x[, block_ell] <- data[lag_ell, ]
  }
  return(list(x = x, y = y))
}

## == LEAST SQUARES FUNCTIONS == ##
#' Least Squares Solution
#'
#' Find the solution to minimize ||y - Xb||_{\ell_2} over b.
#'
#' @param x n x p matrix of regressors.
#' @param y n x 1 vector of responses.
#' @return A list with the following components:
#'   \item{sol}{p x 1 vector of OLS estimates.}
#'   \item{full_rank}{logical; TRUE if regressors are of full rank.}
#' @details
#' This function computes the least squares solution for the linear system
#' (X'X)b = X'y. If X'X is of full rank, it uses the solve function for
#' efficiency. If X'X is rank-deficient, it uses singular value decomposition
#' (SVD) to compute the (Moore-Penrose) solution.
#' @note
#' The Matrix package uses the S4 class system (Chambers, 1998) to retain
#' information on the structure of matrices from the intermediate calculations.
#' A general matrix in dense storage, created by the Matrix function, has class
#' "dgeMatrix" but its cross-product has class "dpoMatrix". The solve methods
#' for the "dpoMatrix" class use the Cholesky decomposition.
#'
#' The rank computation itself involves the SVD of x. However, it is still
#' faster to call Matrix::rankMatrix than to compute the SVD directly and deduce
#' the rank.
#' @examples
#' x <- matrix(rnorm(100), ncol = 5)
#' y <- rnorm(20)
#' fit <- ls_sol(x, y)
#' @export
ls_sol <- function(x, y) {
  # Check for full column rank
  full_rank <- ncol(x) == Matrix::rankMatrix(x)[1]
  if (full_rank) {
    # When full rank, produce unique (OLS) solution
    sol <- solve(crossprod(x), crossprod(x, y))
    # This is (X' X)^{-1} X' y (which equals X^+ y)
  } else {
    # Use SVD for rank-deficient case
    svd_x <- svd(x)                             # compute SVD
    d <- svd_x$d                                # singular values
    u <- svd_x$u                                # left singular vectors
    v <- svd_x$v                                # right singular vectors
    tol <- max(dim(x)) * .Machine$double.eps    # scale free tolerance
    d_inv <- ifelse(d > tol * max(d), 1 / d, 0) # inverse singular values
    sol <- v %*% (d_inv * t(u) %*% y)
    # This is X^+ y
  }
  fit <- list(sol = sol, full_rank = full_rank)
  return(fit)
}

#' Least Squares Refitting with Multiple Responses and Same Regressors
#'
#' This function refits least squares models for multiple responses using the
#' same regressors. It returns a list containing the refitted estimates and
#' a logical vector indicating if the selected regressors are of full rank.
#'
#' @param x A matrix (n x pq) of predictors.
#' @param y A matrix (n x p) of responses.
#' @param that A matrix (p x pq) of estimates.
#'
#' @return A list with the following components:
#' \item{that}{A matrix (p x pq) of refitted estimates.}
#' \item{full_rank}{A logical vector (p x 1); TRUE if selected regressors are
#'                  of full rank.}
#'
#' @examples
#' # Example usage:
#' x <- matrix(rnorm(100), 10, 10)
#' y <- matrix(rnorm(20), 10, 2)
#' that <- matrix(0, 2, 10)
#' result <- mult_refit(x, y, that)
mult_refit <- function(x, y, that) {
  p <- ncol(y)                          # num responses
  sel <- that != 0                      # flag selections
  full_rank <- logical(p)               # full rank flag (initialize FALSE)
  num_selections <- rowSums(sel)        # precompute num selections
  for (i in which(num_selections > 0)) {
    xi <- as.matrix(x[, sel[i, ]])      # get active regressors
    fit_i <- ls_sol(xi, y[, i])         # refit
    that[i, sel[i, ]] <- fit_i$sol      # overwrite (keeping zeros)
    full_rank[i] <- fit_i$full_rank     # check rank
  }
  full_rank[num_selections == 0] <- TRUE  # default full rank for no selections
  refit <- list(that = that, full_rank = full_rank)
  return(refit)
}

## == LASSO FUNCTIONS == ##
#' Multiple LASSO using glmnet
#'
#' This function performs multiple LASSO regression using the glmnet package,
#' optionally allowing for variable rescaling and choice between single lambda
#' or full regularization path solutions.
#'
#' @param x A matrix of predictors (observations x variables).
#' @param y A matrix of responses (observations x outcomes).
#' @param lambda_glmnet Numeric. Regularization parameter for glmnet.
#' @param upsilon Optional. Matrix of rescaling factors for predictors.
#'   Defaults to unit loadings.
#' @param full_path Logical. If TRUE, uses full regularization path of glmnet.
#'   If FALSE (default), uses single lambda specified by lambda_glmnet.
#' @param tol_glmnet Numeric. Convergence threshold for glmnet.
#'
#' @return A matrix of coefficient estimates for each response variable.
#'
#' @details
#' - For each response variable, predictors are rescaled using upsilon if provided.
#' - If full_path = FALSE, glmnet is used with a single lambda specified by
#'   lambda_glmnet.
#' - If full_path = TRUE, glmnet's own lambda sequence is used, and coefficients
#'   are interpolated to obtain estimates corresponding to lambda_glmnet.
#'
#' @seealso \code{\link[glmnet]{glmnet}}
#' @keywords regression
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' n <- 100
#' p <- 5
#' q <- 3
#' x <- matrix(rnorm(n * p), n, p)
#' y <- matrix(rnorm(n * q), n, q)
#' # Perform multiple LASSO
#' mult_lasso(x, y, lambda_glmnet = 0.1)
#'
#' @export
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
## OLD VERSION BELOW ##
# # Multiple responses weighted LASSO using same regressors
# # INPUTS
# #   x: n x pq matrix of predictors
# #   y: n x p matrix of responses
# #   lambda_glmnet: penalty level in eyes of glmnet
# #   upsilon: p x pq matrix of penalty loadings
# #   full_path: logical; if TRUE, use full path of lambda values
# #   tol_glmnet: convergence tolerance for glmnet
# # OUTPUT
# #   that: p x pq matrix of estimates
# mult_lasso <- function(x, y, lambda_glmnet, upsilon = NULL,
#                        full_path = FALSE, tol_glmnet = 1e-4) {
#   if (missing(x) || missing(y) || missing(lambda_glmnet)) {
#     stop("Not enough input arguments.")
#   }
#   pq <- ncol(x)  # columns in x = pq = total number of regressors
#   p <- ncol(y)   # q = autoregressive order, p = dim(output)
#   if (nrow(x) != nrow(y)) {
#     stop("Number of observations (rows in x and y) do not match.")
#   }
#   if (is.null(upsilon)) {
#     upsilon <- matrix(1, p, pq) # default to unit loadings
#   }
#   that <- matrix(NA, p, pq) # placeholder for estimates
#   if (full_path == FALSE) { # if full path *not* requested (= default)...
#     for (i in 1:p) {        # use glmnet w/ single lambda
#       xtilde <- sweep(x, 2, upsilon[i, ], FUN = "/")  # rescale using loadings
#       fit_i <- glmnet(xtilde, y[, i], family = "gaussian",
#                       standardize = FALSE, intercept = FALSE,
#                       lambda = lambda_glmnet, thresh = tol_glmnet)
#       that[i, ] <- coef(fit_i)[-1] / upsilon[i, ] # original scaling
#     }
#   } else {            # if full path requested (not default)...
#     for (i in 1:p) {  # use glmnet's lambda sequence and interpolate
#       xtilde <- sweep(x, 2, upsilon[i, ], FUN = "/")  # rescale using loadings
#       fit_i <- glmnet(xtilde, y[, i], family = "gaussian",
#                       standardize = FALSE, intercept = FALSE,
#                       thresh = tol_glmnet)
#       # ^-- glmnet constructs its own lambda sequence
#       that_i_temp <- coef(fit_i, s = lambda_glmnet)[-1] # interpolate
#       that[i, ] <- that_i_temp / upsilon[i, ]           # original scale
#     }
#   }
#   return(that) # return as sparse matrix
# }

## == INFORMATION CRITERIA FUNCTIONS == ##
#' Compute Information Criteria for Lasso Regression using glmnet
#'
#' This function fits a Lasso regression model using glmnet and calculates
#' information criteria (AIC, BIC, HQIC) to select the best model.
#'
#' @param x Matrix of predictors (n x p).
#' @param y Vector of response (n x 1).
#' @param criteria Character vector specifying the information criteria to use.
#'                 Default is c("aic", "bic", "hqic").
#' @return A list with components:
#'   \item{betas}{Matrix of coefficients (betas), with columns named after criteria.}
#'   \item{lambdas}{Vector of penalty parameters (lambdas) that minimize each criterion.}
#' @examples
#' x <- matrix(rnorm(100), ncol = 5)
#' y <- rnorm(20)
#' ic_lasso(x, y)
#' @import glmnet
#' @export
ic_lasso <- function(x, y, criteria = c("aic", "bic", "hqic")) {
  # Match and validate criteria
  criteria <- match.arg(criteria, several.ok = TRUE)
  # Fit the model using glmnet
  fit <- glmnet(x, y, family = "gaussian",
                intercept = FALSE, standardize = FALSE)
  # Extract degrees of freedom
  df <- fit$df
  n_crit <- length(criteria)
  betas <- matrix(NA, ncol(x), n_crit)
  colnames(betas) <- criteria
  lambdas <- numeric(n_crit)
  # Predict outcome for each candidate lambda
  yhats <- x %*% fit$beta # n x num(lambda) matrix of predictions
  minus_res <- sweep(yhats, 1, y, "-") # implied negative residuals
  # Calculate sum of squared residuals for each candidate lambda
  ssr <- colSums(minus_res^2)
  # Calculate number of observations
  n <- length(y)
  n_log_ssr <- log(ssr) * n # n times log(SSR)
  # Compute information criteria for each specified criterion
  for (i in seq_along(criteria)) {
    crit <- criteria[i]
    if (crit == "aic") {
      ic <- n_log_ssr + 2 * df
    } else if (crit == "bic") {
      ic <- n_log_ssr + log(n) * df
    } else if (crit == "hqic") {
      ic <- n_log_ssr + 2 * log(log(n)) * df
    } else {
      stop("Unknown information criterion")
    }
    id_min <- which.min(ic)         # identify minimizer
    lambda_ic <- fit$lambda[id_min] # minimizing lambda
    beta_ic <- fit$beta[, id_min]   # corresponding estimates
    # Store results
    betas[, i] <- beta_ic
    lambdas[i] <- lambda_ic
  }
  # Return results as a list
  list(betas = betas, lambdas = lambdas)
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
sqrt_lasso <- function(x, y, lambda, max_iter = 100, rel_tol_norm = 1e-4) {
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
    if (norm_diff / norm_beta_old < rel_tol_norm) {         # check tolerance
      break                                                 # break if met
    }
  }
  return(beta)
}

#' Estimate equation-by-equation square-root LASSO coefficients for VAR models.
#'
#' This function estimates coefficients for a Vector Autoregressive (VAR) model
#' using the square-root LASSO method, where predictors are rescaled based on
#' penalty loadings.
#'
#' @param x n x pq matrix of predictors.
#' @param y n x p matrix of responses.
#' @param lambda Penalty level for square-root LASSO regularization.
#' @param upsilon p x pq matrix of penalty loadings. Each row corresponds to the
#'        loadings for one equation in the VAR model. If NULL, defaults to unit
#'        loadings.
#' @param max_iter Maximum number of iterations for the square-root LASSO
#'        algorithm. Default is 100.
#' @param rel_tol_norm Stopping tolerance for ||beta - beta_old|| in the
#'        square-root LASSO algorithm. Default is 1e-4.
#'
#' @return p x pq matrix of estimated coefficients for each equation in the VAR
#'         model.
#'
#' @details
#' - The function rescales predictors \code{x} using provided penalty loadings
#'   \code{upsilon} before applying square-root LASSO estimation.
#' - Each equation in the VAR model is estimated separately, iterating through
#'   \code{p} equations.
#' - Estimated coefficients are rescaled back to original scale using inverse of
#'   \code{upsilon}.
#' - Ensure \code{x} and \code{y} have same number of observations (rows).
#'
#' @examples
#' # Example usage:
#' x_data <- matrix(rnorm(100 * 6), 100, 6)
#' y_data <- matrix(rnorm(100 * 2), 100, 2)
#' lambda <- 0.1
#' penalty_loadings <- matrix(runif(6), 2, 3)
#' estimates <- mult_sqrt_lasso(x_data, y_data, lambda, upsilon = penalty_loadings)
#'
#' @export
mult_sqrt_lasso <- function(x, y, lambda, upsilon = NULL,
                            max_iter = 100, rel_tol_norm = 1e-4) {
  # Check input validity
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
  if (nrow(upsilon) != p || ncol(upsilon) != pq) {
    stop("Invalid dimensions for upsilon.")
  }
  that <- matrix(NA, p, pq) # placeholder for estimates
  for (i in 1:p) {
    # Rescale x for the i-th equation
    xtilde_i <- sweep(x, 2, upsilon[i, ], FUN = "/")
    # Perform sqrt_lasso estimation with max_iter and rel_tol_norm
    that_i_temp <- sqrt_lasso(xtilde_i, y[, i], lambda,
                              max_iter = max_iter,
                              rel_tol_norm = rel_tol_norm)
    # Store estimates back to original scale
    that[i, ] <- that_i_temp / upsilon[i, ]
  }
  return(that) # return as matrix
}