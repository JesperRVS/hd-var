#' Lasso VAR Model Fitting
#'
#' This function fits a vector autoregressive (VAR) model with Lasso
#' regularization and provides options for refining the estimates.
#'
#' @param data A (q + n) x p matrix of data (n = effective sample size).
#' @param q Autoregressive order; default is 1 (i.e., VAR(1) model).
#' @param post Logical; if TRUE, refit estimates after each loading update.
#' @param intercept Logical; if TRUE, estimate intercepts.
#' @param c Tuning parameter for penalty level (score markup). Default is 1.1,
#'   as in Belloni et al. (2012, ECTA).
#' @param gamma Tuning parameter for penalty level (score quantile). Default is
#'   0.1 / log(max(dim(data))), as in Belloni et al. (2012, ECTA).
#' @param k Maximum number of penalty loading updates. Default is 15 , as in
#'   Belloni et al. (2012, ECTA).
#' @param tol_ups Convergence tolerance for penalty loadings (relative), so
#'   that penalty loadings are considered converged if
#'   ||vec(ups_new - ups_old)||_2 / ||vec(ups_old)||_2 <= tol_ups.
#' @param warn Logical; if TRUE, issue warning if max updates reached.
#' @param full_path Logical; if TRUE, use full path of lambda values. Default is
#'   FALSE.
#' @param tol_glmnet Convergence tolerance for glmnet. (Mainly for debugging.)
#'
#' @return A list with the following components:
#' \item{lambda}{Penalty level used}
#' \item{ups_init}{Initial penalty loadings}
#' \item{ups_refi}{Refined penalty loadings path}
#' \item{ups}{Final penalty loadings}
#' \item{intr_init}{Initial intercepts (if requested)}
#' \item{intr_refi}{Refined intercepts path (if requested)}
#' \item{intr}{Final intercepts (if requested)}
#' \item{that_init}{Initial estimates}
#' \item{that_refi}{Refined estimates path}
#' \item{that}{Final estimates}
#' \item{full_rank_post_init}{Logicals; if TRUE, selected regressors full rank}
#' \item{full_rank_post_refi}{Logicals; if TRUE, selected regressors full rank}
#' \item{k_term}{Number of penalty loading updates performed (up to k)}
#' \item{rel_diff_ups_term}{Relative change in penalty loadings at termination}
#' \item{rel_diffs_ups}{Relative changes in penalty loadings along path}
#'
#' @importFrom stats qnorm
#' @importFrom glmnet glmnet
#' @importFrom Matrix crossprod
#'
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(1000), 100, 10)
#' fit <- lasso_var(data)
#' }
#'
#' @export
source("helper_functions.R")
lasso_var <- function(data, q = 1, post = TRUE, intercept = TRUE,
                      c = 1.1, gamma = NULL, k = 15,
                      tol_ups = 1e-3, warn = TRUE, full_path = FALSE,
                      tol_glmnet = 1e-4) {
  # Unpack data and construct predictors and response
  xy_data <- unpack(data = data, q = q)
  x <- xy_data$x
  y <- xy_data$y
  # Call upon mult_lasso_bcch
  fit <- mult_lasso_bcch(x = x, y = y, post = post, intercept = intercept,
                         c = c, gamma = gamma, k = k,
                         tol_ups = tol_ups, warn = warn, full_path = full_path,
                         tol_glmnet = tol_glmnet)
  return(fit)
}


# # OLD VERSION BELOW
# lasso_var <- function(data, q = 1, post = TRUE, intercept = TRUE,
#                       c = 1.1, gamma = NULL, k = 15,
#                       tol_ups = 1e-3,  warn = TRUE, full_path = FALSE,
#                       tol_glmnet = 1e-4) {
#   # INITIALIZE
#   if (missing(data) || is.null(data)) {
#     stop("No data provided.")
#   }
#   # Unpack data and construct predictors and response
#   nplusq <- nrow(data)  # n plus q
#   p <- ncol(data)       # p = dim of output
#   n <- nplusq - q       # n = effective sample size
#   pq <- p * q           # pq = total number of variables
#   xy_data <- unpack(data = data, q = q) # unpack data
#   x <- xy_data$x        # predictors
#   y <- xy_data$y        # responses

#   # Probability tolerance
#   if (is.null(gamma)) {
#     gamma <-  0.1 / log(max(c(n, pq)))
#   }

#   # Penalty level
#   lambda_star <- 2 * c * sqrt(n) * qnorm(1 - gamma / (2 * pq * p))
#   lambda_glmnet <- lambda_star / (2 * n) # penalty in eyes of glmnet
#   # Note: glmnet defines "lambda" based on *half* of *average* square loss,
#   # while our "lambda" stems from *sum* square loss (w/o the 1/2)

#   # If intercepts requested, demean before proceeding
#   if (intercept) {
#     ybar <- colMeans(y)     # response means as p-dim array
#     xbar <- colMeans(x)     # predictor means as pq-dim array
#     y <- sweep(y, 2, ybar)  # demean responses
#     x <- sweep(x, 2, xbar)  # demean predictors
#     ybar <- as.matrix(ybar) # response means as p x 1 matrix
#     xbar <- as.matrix(xbar) # predictor means as pq x 1 matrix
#   }
#   # Note: Means are stored to back out intercepts later

#   # ESTIMATION
#   # INITIAL STEP (k = 0)
#   # Initial penalty loadings and estimates
#   ups_init <- sqrt((1 / n) * crossprod(y^2, x^2))
#   that_init <- mult_lasso(x, y, lambda_glmnet, upsilon = ups_init,
#                           full_path = full_path, tol_glmnet = tol_glmnet)
#   if (post) {                               # if refitting requested...
#     refit <- mult_refit(x, y, that_init)    # refit initial estimates
#     that_init <- refit$that                 # overwrite (keeping zeros)
#     full_rank_post_init <- refit$full_rank  # flag full rank
#   } else { # if no refitting requested...
#     full_rank_post_init <- NULL # then rank deficiency is irrelevant
#   }
#   if (intercept) { # if intercepts requested...
#     intr_init <- ybar - that_init %*% xbar # back them out
#   } else {
#     intr_init <- NULL # o/w don't
#   }

#   # UPDATING (up to K times)
#   # Stop if relative change in penalty loadings is small enough
#   # Create placeholders
#   ups_refi <- array(NA, dim = c(p, pq, k))  # refined penalty loadings
#   that_refi <- array(NA, dim = c(p, pq, k)) # and estimates
#   rel_diffs_ups <- numeric(k)               # relative changes in loadings
#   if (post) { # if refitting requested...
#     full_rank_post_refi <- matrix(NA, p, k) # flag full rank
#   } else { # o/w don't
#     full_rank_post_refi <- NULL
#   }
#   if (intercept) { # if intercepts requested...
#     intr_refi <- matrix(NA, p, k) # create placeholder
#   } else {
#     intr_refi <- NULL # o/w don't
#   }
#   for (l in 1:k) {                       # start updating
#     if (l == 1) {                        # if first iteration...
#       that_old <- that_init              # use initial estimates
#     } else {
#       that_old <- that_refi[, , l - 1]   # o/w use previous (refined) estimates
#     }
#     res_old <- y - x %*% t(that_old) # implied residuals
#     ups_new <- sqrt((1 / n) * crossprod(res_old^2, x^2)) # update penalty load.
#     ups_refi[, , l] <- ups_new # store
#     that_new <- mult_lasso(x, y, lambda_glmnet, ups_new, tol_glmnet) # estimates
#     if (post) {                            # if refitting requested...
#       refit <- mult_refit(x, y, that_new)  # refit estimates
#       that_new <- refit$that               # overwrite (keeping zeros)
#       full_rank_post_refi[, l] <- refit$full_rank # flag full rank
#     }
#     that_refi[, , l] <- that_new           # store
#     if (intercept) {                       # if intercepts requested...
#       intr_refi[, l] <- ybar - that_new %*% xbar  # back them out
#     }
#     if (l == 1) {          # if at first update...
#       ups_old <- ups_init  # then compare w/ initial penalty loadings
#     } else {               # o/w compare w/ previous (refined) penalty loadings
#       ups_old <- ups_refi[, , l - 1]
#     }
#     diff_ups <- ups_new - ups_old # change in penalty loadings
#     # relative change in vectorized ell_2 norm --v
#     rel_diff_ups <- sqrt(sum(diff_ups^2)) /
#       (sqrt(sum(ups_old^2)) + .Machine$double.eps)
#     rel_diffs_ups[l] <- rel_diff_ups # store
#     if (rel_diff_ups <= tol_ups) { # if change is small enough...
#       break # stop updating
#     }
#   }
#   ups <- ups_refi[, , l] # set refined penalty loadings
#   if (intercept) {       # if intercepts requested...
#     intr <- intr_refi[, l] # set intercepts
#   } else { # o/w don't
#     intr <- NULL
#   }
#   that <- that_refi[, , l] # set estimates

#   # OPTIONAL WARNING
#   if (warn && l == k && rel_diff_ups > tol_ups) { # if max updates reached...
#     warning("Maximum number of updates reached. Consider increasing K.")
#     cat(sprintf("Relative change in penalty loadings is %3.1g percent\n", 
#                 100 * rel_diff_ups))
#   }

#   # OUTPUT
#   fit <- list( # store results in list
#     lambda = lambda_star,               # penalty level used
#     ups_init = ups_init,                # initial penalty loadings
#     ups_refi = ups_refi[, , 1:l],       # refined penalty loadings path
#     ups = ups,                          # final penalty loadings
#     intr_init = intr_init,              # initial intercepts
#     intr_refi = intr_refi[, 1:l],       # refined intercepts path
#     intr = intr,                        # final intercepts
#     that_init = that_init,              # initial estimates
#     that_refi = that_refi[, , 1:l],     # refined estimates path
#     that = that,                        # final estimates
#     full_rank_post_init = full_rank_post_init, # full rank flags
#     full_rank_post_refi = full_rank_post_refi[, 1:l], # (if requested)
#     k_term = l,                         # number of updates performed
#     rel_diff_ups_term = rel_diff_ups,   # rel change in loadings at termination
#     rel_diffs_ups = rel_diffs_ups[1:l]  # rel changes in loadings along path
#   )
#   return(fit) # return list
# }