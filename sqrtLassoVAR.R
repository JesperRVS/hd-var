#' Square-Root LASSO Vector Autoregression (VAR) Estimation
#'
#' Estimate vector autoregression (VAR) coefficients using the square-root LASSO method.
#'
#' @param data A numeric matrix or data frame where rows correspond to observations and columns to variables.
#' @param q Integer indicating the lag order of the VAR model.
#' @param post Logical indicating whether to refit estimates after selection (default is TRUE).
#' @param intercept Logical indicating whether to include intercepts in the model (default is TRUE).
#' @param c Penalty scaling parameter for square-root LASSO.
#' @param gamma Tuning parameter for square-root LASSO.
#' @param upsilon Optional p x (p*q) matrix of loadings for the square-root LASSO.
#' @param max_iter Maximum number of iterations for the square-root LASSO algorithm.
#' @param rel_tol_norm Relative tolerance for convergence of the square-root LASSO algorithm based on norm.
#'
#' @return A list with components:
#'   \item{intr}{Vector of intercept estimates if intercept = TRUE, otherwise NULL.}
#'   \item{that}{Estimated VAR coefficients matrix.}
#'   \item{full_rank_post}{Logical indicating whether the refitted model is full rank, if post = TRUE, otherwise NULL.}
#'
#' @details
#' This function estimates the coefficients of a VAR model using the square-root LASSO method,
#' which combines L1 penalization with a square-root transformation. The estimation procedure
#' includes options to include intercepts, specify penalty levels, and tune convergence criteria.
#'
#' @references
#' Yuan, M., & Lin, Y. (2007). Model selection and estimation in the Gaussian graphical model.
#' Biometrika, 94(1), 19-35.
#'
#' @seealso \code{\link{mult_sqrt_lasso}}, \code{\link{mult_refit}}, \code{\link{sweep}}
#'
#' @importFrom MASS ginv
#' @importFrom Matrix rankMatrix Diagonal
#' @import helper_functions.R
#'
#' @examples
#' data <- matrix(rnorm(1000), ncol = 5)  # Example data matrix with 5 variables
#' result <- sqrt_lasso_var(data, q = 2)
#'
#' @export
library("MASS", "Matrix")
# Notes:
#   MASS is used for Moore-Penrose pseudoinverse (ginv)
#   Matrix is used for rankMatrix and Diagonal
source("helper_functions.R")  # for mult_sqrt_lasso and refitting
sqrt_lasso_var <- function(data, q = 1, post = TRUE, intercept = TRUE,
                           c = 1.1, gamma = 0.1 / log(max(dim(data))),
                           upsilon = NULL,
                           max_iter = 100, rel_tol_norm = 1e-4) {
  if (missing(data) || is.null(data)) {
    stop("No data provided.")
  }
  # INITIALIZE
  # Unpack data and construct predictors and response
  nplusq <- nrow(data)  # n plus q
  p <- ncol(data)       # p = dim of output
  n <- nplusq - q       # n = effective sample size
  pq <- p * q           # pq = total number of variables
  xy_data <- unpack(data, q = q) # unpack data
  x <- xy_data$x        # predictors
  y <- xy_data$y        # responses
  # Penalty level in eyes of square-root LASSO
  lambda_star_sqrt <- c * sqrt(n) * qnorm(1 - gamma / (2 * p * pq))
  # Note: Sqrt-LASSO objective cancels out the "2" in the KKT conditions
  if (is.null(upsilon)) {
    upsilon <- matrix(1, p, pq) # default to unit loadings
  }
  if (nrow(upsilon) != p || ncol(upsilon) != pq) {
    stop("Invalid dimensions for upsilon.")
  }
  # If intercepts requested, demean before proceeding
  if (intercept) {
    ybar <- colMeans(y)     # response means as p x 1 vector
    xbar <- colMeans(x)     # predictor means as pq x 1 vector
    y <- sweep(y, 2, ybar)  # demean responses
    x <- sweep(x, 2, xbar)  # demean predictors
  }
  # Note: Means are stored to back out intercepts later
  # ESTIMATE
  # Estimate the square-root LASSO VAR
  that <- mult_sqrt_lasso(x, y, lambda = lambda_star_sqrt,
                                    upsilon = upsilon,
                                    max_iter = max_iter, 
                                    rel_tol_norm = rel_tol_norm)
  # POST-PROCESS
  # If requested, refit estimates after selection
  if (post) {
    refit <- mult_refit(x, y, that)
    that <- refit$that
    full_rank_post <- refit$full_rank
  } else {
    full_rank_post <- NULL
  }
  # Back out intercepts if requested
  if (intercept) {
    intr <- ybar - that %*% xbar
  } else {
    intr <- NULL
  }
  return(list(intr = intr, that = that, full_rank_post = full_rank_post))
}


# ## Dependencies
# library("MASS", "Matrix")
# # Notes:
# #   MASS is used for Moore-Penrose pseudoinverse (ginv)
# #   Matrix is used for rankMatrix and Diagonal
# source("helper_functions.R")  # for mult_sqrt_lasso and refitting
# sqrt_lasso_var <- function(data, q = 1, post = TRUE, intercept = TRUE,
#                            c = 1.1, gamma = 0.1 / log(max(dim(data))),
#                            upsilon = NULL,
#                            max_iter = 100, rel_tol_norm = 1e-4) {
#   if (missing(data) || is.null(data)) {
#     stop("No data provided.")
#   }
#   # INITIALIZE
#   # Unpack data and construct predictors and response
#   nplusq <- nrow(data)  # n plus q
#   p <- ncol(data)       # p = dim of output
#   n <- nplusq - q       # n = effective sample size
#   pq <- p * q           # pq = total number of variables
#   xy_data <- unpack(data, q = q) # unpack data
#   x <- xy_data$x        # predictors
#   y <- xy_data$y        # responses
#   # Penalty level in eyes of square-root LASSO
#   lambda_star_sqrt <- c * sqrt(n) * qnorm(1 - gamma / (2 * p * pq))
#   # Note: Sqrt-LASSO objective cancels out the "2" in the KKT conditions
#   if (is.null(upsilon)) {
#     upsilon <- matrix(1, p, pq) # default to unit loadings
#   }
#   if (nrow(upsilon) != p || ncol(upsilon) != pq) {
#     stop("Invalid dimensions for upsilon.")
#   }
#   # If intercepts requested, demean before proceeding
#   if (intercept) {
#     ybar <- colMeans(y)     # response means as p x 1 vector
#     xbar <- colMeans(x)     # predictor means as pq x 1 vector
#     y <- sweep(y, 2, ybar)  # demean responses
#     x <- sweep(x, 2, xbar)  # demean predictors
#   }
#   # Note: Means are stored to back out intercepts later
#   # ESTIMATE
#   # Estimate the square-root LASSO VAR
#   that <- mult_sqrt_lasso(x, y, lambda = lambda_star_sqrt,
#                                     upsilon = upsilon,
#                                     max_iter = max_iter, 
#                                     rel_tol_norm = rel_tol_norm)
#   # POST-PROCESS
#   # If requested, refit estimates after selection
#   if (post) {
#     refit <- mult_refit(x, y, that)
#     that <- refit$that
#     full_rank_post <- refit$full_rank
#   } else {
#     full_rank_post <- NULL
#   }
#   # Back out intercepts if requested
#   if (intercept) {
#     intr <- ybar - that %*% xbar
#   } else {
#     intr <- NULL
#   }
#   return(list(intr = intr, that = that, full_rank_post = full_rank_post))
# }