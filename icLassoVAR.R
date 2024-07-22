#' Estimate a VAR Model with LASSO and Information Criteria Selection
#'
#' This function estimates a VAR model using LASSO and information criteria 
#' selection. It supports post-selection refitting and back out intercepts.
#'
#' @param data A matrix of data with dimensions (q + n) x p where n is the 
#'   effective sample size.
#' @param q Integer. Autoregressive order. Default is 1 (i.e., VAR(1) model).
#' @param criteria A character vector of information criteria to use. Default 
#'   is c("aic", "bic", "hqic").
#' @param post Logical. If TRUE, refit estimates after selection. Default is 
#'   TRUE.
#' @param intercept Logical. If TRUE, back out intercepts. Default is TRUE.
#' @param upsilon A matrix of loadings with dimensions p x pq. Default is 
#'   unit loadings.
#' @return A list with the following components:
#'   \item{intrs}{A matrix of intercepts with dimensions p x num_crit.}
#'   \item{thats}{An array of estimates with dimensions p x pq x num_crit.}
#'   \item{full_rank_post}{A matrix of full rank flags with dimensions 
#'   p x num_crit.}
#' @import glmnet
#' @import MASS
#' @import Matrix
#' @export
#'
#' @examples
#' # Example usage:
#' data <- matrix(rnorm(1000), ncol = 10)
#' result <- ic_lasso_var(data)
source("helper_functions.R")  # for ic_lasso and refitting
ic_lasso_var <- function(data, q = 1, criteria = c("aic", "bic", "hqic"),
                         post = TRUE, intercept = TRUE, upsilon = NULL) {
  # Check if data is provided
  if (missing(data) || is.null(data)) {
    stop("No data provided.")
  }
  # Match the criteria against allowed values
  criteria <- match.arg(criteria, several.ok = TRUE)
  # Get dimensions of data
  p <- ncol(data)
  pq <- p * q
  xy <- unpack(data, q = q) # Unpack data into predictors and responses
  x <- xy$x
  y <- xy$y
  num_crit <- length(criteria)
  # Set default loadings if upsilon is not provided
  if (is.null(upsilon)) {
    upsilon <- matrix(1, p, pq)
  }
  # Check dimensions of upsilon matrix
  if (nrow(upsilon) != p || ncol(upsilon) != pq) {
    stop("Invalid dimensions for upsilon.")
  }
  # Demean data if intercepts are requested
  if (intercept) {
    ybar <- colMeans(y) # Mean of responses
    xbar <- colMeans(x) # Mean of predictors
    y <- sweep(y, 2, ybar, "-") # Demean responses
    x <- sweep(x, 2, xbar, "-") # Demean predictors
    # Convert means to matrices for later use
    ybar <- matrix(ybar, nrow = p, ncol = 1, byrow = TRUE)
    xbar <- matrix(xbar, nrow = pq, ncol = 1, byrow = TRUE)
    # Placeholder for intercepts
    intrs <- matrix(NA, nrow = p, ncol = num_crit)
    colnames(intrs) <- criteria
  }
  # Placeholder for estimates
  thats <- array(NA, dim = c(p, pq, num_crit))
  dimnames(thats) <- list(1:p, 1:pq, criteria)
  # Loop through each response variable to fit LASSO
  for (i in 1:p) { 
    upsilon_i <- upsilon[i, ] # Current loadings for the response
    xtemp <- sweep(x, 2, upsilon_i, "/") # Scale predictors
    fit_ic <- ic_lasso(xtemp, y[, i], criteria = criteria) # Fit LASSO
    thats_i_temp <- fit_ic$betas # Extract coefficients
    thats[i, , ] <- sweep(thats_i_temp, 1, upsilon_i, "/") # Scale back
  }
  # Refitting estimates if requested
  if (post) {
    full_rank_post <- matrix(NA, nrow = p, ncol = num_crit)
    colnames(full_rank_post) <- criteria
    for (crit in criteria) {
      refit <- mult_refit(x, y, thats[, , crit])
      thats[, , crit] <- refit$that # Update estimates
      full_rank_post[, crit] <- refit$full_rank # Check rank
    }
  } else {
    full_rank_post <- NULL
  }
  # Back out intercepts if requested
  if (intercept) {
    for (crit in criteria) {
      intrs[, crit] <- ybar - thats[, , crit] %*% xbar
    }
  } else {
    intrs <- NULL
  }
  # Return results as a list
  return(list(intrs = intrs, thats = thats,
              full_rank_post = full_rank_post))
}
## OLD VERSION BELOW ##
# ## Dependencies
# # library("glmnet", "MASS", "Matrix")
# # Notes:
# #   glmnet is used for LASSO estimation
# #   MASS is used for Moore-Penrose pseudoinverse (ginv)
# #   Matrix is used for rankMatrix
# source("helper_functions.R")  # for ic_lasso and refitting

# # Function to estimate a VAR w/ LASSO and information criteria selection
# # INPUTS
# #   data:           (q + n) x p matrix of data (n = effective sample size)
# #   q:              autoregressive order; default is 1 (i.e. VAR(1) model)
# #   criteria:       information criteria to use; default is AIC, BIC and HQIC
# #   post:           logical; if TRUE, refit estimates after selection
# #   intercept:      logical; if TRUE, back out intercepts
# #   upsilon:        p x pq matrix of loadings; default is unit loadings
# # OUTPUTs a list with the following components:
# #   intrs:          p x num_crit matrix of intercepts
# #   thats:          p x pq x num_crit array of estimates
# #   full_rank_post: p x num_crit matrix of full rank flags
# ic_lasso_var <- function(data, q = 1, criteria = c("aic", "bic", "hqic"),
#                          post = TRUE, intercept = TRUE, upsilon = NULL) {
#   if (missing(data) || is.null(data)) {
#     stop("No data provided.")
#   }
#   criteria <- match.arg(criteria, several.ok = TRUE)
#   # Unpack data and construct predictors and response
#   nplusq <- nrow(data)      # n plus q
#   p <- ncol(data)           # p = dim of output
#   pq <- p * q               # pq = total number of variables
#   xy <- unpack(data, q = q) # unpack data
#   x <- xy$x                 # predictors
#   y <- xy$y                 # responses
#   num_crit <- length(criteria)  # num information criteria
#   if (is.null(upsilon)) {
#     upsilon <- matrix(1, p, pq) # default to unit loadings
#   }
#   if (nrow(upsilon) != p || ncol(upsilon) != pq) {
#     stop("Invalid dimensions for upsilon.")
#   }
#   # if intercept requested, demean before proceeding
#   if (intercept) {
#     ybar <- colMeans(y)           # response means as p-dim array
#     xbar <- colMeans(x)           # predictor means as pq-dim array
#     y <- sweep(y, 2, ybar, "-")   # demean responses
#     x <- sweep(x, 2, xbar, "-")   # demean predictors
#     ybar <- matrix(ybar, nrow = p, ncol = 1, byrow = TRUE)
#     xbar <- matrix(xbar, nrow = pq, ncol = 1, byrow = TRUE)
#     intrs <- matrix(NA, nrow = p, ncol = num_crit) # placeholder intercepts
#     colnames(intrs) <- criteria
#   }
#   # Note: Means are stored to back out intercepts later
#   thats <- array(NA, dim = c(p, pq, num_crit))  # placeholder slopes
#   dimnames(thats) <- list(1:p, 1:pq, criteria)
#   for (i in 1:p) { # equation-by-equation information criteria LASSO
#     upsilon_i <- upsilon[i, ] # ith row of upsilon
#     xtemp <- sweep(x, 2, upsilon_i, "/")              # scale predictors
#     fit_ic <- ic_lasso(xtemp, y[, i], criteria = criteria)
#     thats_i_temp <- fit_ic$betas                      # slopes on equal scale
#     thats_i <- sweep(thats_i_temp, 1, upsilon_i, "/") # and on original scale
#     thats[i, , ] <- thats_i
#   }
#   # Refit estimates for each information criterion
#   if (post) {
#     full_rank_post <- array(NA, dim = c(p, num_crit)) # rank check placeholder
#     dimnames(full_rank_post) <- list(1:p, criteria)
#     for (crit in criteria) {
#       refit <- mult_refit(x, y, thats[, , crit])  # refit selection
#       thats[, , crit] <- refit$that               # overwrite (keeping zeros)
#       full_rank_post[, crit] <- refit$full_rank   # flag full rank
#     }
#   } else {
#     full_rank_post <- NULL
#   }
#   if (intercept) {                            # if intercepts requested...
#     for (crit in criteria) {
#       intrs[, crit] <- ybar - thats[, , crit] %*% xbar # ... back them out
#     }
#   } else {
#     intrs <- NULL
#   }
#   return(list(intrs = intrs, thats = thats, full_rank_post = full_rank_post))
# }