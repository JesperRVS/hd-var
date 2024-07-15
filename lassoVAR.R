### ESTIMATION TOOLS FOR VECTOR AUTOREGRESSIONS WITH LASSO PENALIZATION

## Dependencies
library("glmnet", "MASS", "Matrix")
# Notes:
#   glmnet is used for LASSO estimation
#   MASS is used for Moore-Penrose pseudoinverse (ginv)
#   Matrix is used for rankMatrix
source("helper_functions.R")  # for multiple Lassos and refitting

## MAIN FUNCTION
# INPUTS
#   data: (q + n) x p matrix of data (n = effective sample size)
#   q: autoregressive order; default is 1 (i.e. VAR(1) model)
#   post: logical; if TRUE, refit estimates after *each* loading update
#   intercept: logical; if TRUE, estimate intercepts
#   c: tuning parameter for penalty level (score markup). Default is 1.1, as
#     in Belloni et al. (2012, ECTA).
#   gamma: tuning parameter for penalty level (score quantile). Default is
#     0.1 / log(max(dim(data))), as in Belloni et al. (2012, ECTA).
#   k: maximum number of penalty loading updates. Default is 15.
#   tol_ups: convergence tolerance for penalty loadings (relative), so that
#     penalty loadings are considered converged if ||vec(ups_new -
#     ups_old)||_{\ell_2} / ||vec(ups_old)||_{\ell_2} <= tol_ups
#   warn: logical; if TRUE, issue warning if max updates reached
#   full_path: logical; if TRUE, use full path of lambda values. Default is
#     FALSE.
#   tol_glmnet: convergence tolerance for glmnet. (Mainly for debugging.)
# OUTPUTs fit: list with the following components
#   lambda: penalty level used
#   ups_init: initial penalty loadings
#   ups_refi: refined penalty loadings path
#   ups: final penalty loadings (i.e. after convergence)
#   intr_init: initial intercepts (if requested)
#   intr_refi: refined intercepts path (if requested)
#   intr: final intercepts (i.e. after loading convergence, if requested)
#   that_init: initial estimates
#   that_refi: refined estimates path
#   that: final estimates (i.e. after loading convergence)
#   full_rank_post_init: logicals; if TRUE, selected regressors full rank
#   full_rank_post_refi: logicals; if TRUE, selected regressors full rank
#   k_term: number of penalty loading updates performed (up to k)
#   rel_diff_ups_term: relative change in penalty loadings at termination
#   rel_diffs_ups: relative changes in penalty loadigns along path
lasso_var <- function(data, q = 1, post = TRUE, intercept = TRUE,
                      c = 1.1, gamma = 0.1 / log(max(dim(data))), k = 15,
                      tol_ups = 1e-3,  warn = TRUE,
                      full_path = FALSE, tol_glmnet = 1e-4) {
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
  # Penalty level
  lambda_star <- 2 * c * sqrt(n) * qnorm(1 - gamma / (2 * q * p^2))
  lambda_glmnet <- lambda_star / (2 * n) # penalty in eyes of glmnet
  # Note: glmnet defines "lambda" based on *half* of *average* square loss,
  # while our "lambda" stems from *sum* square loss (w/o the 1/2)
  # If intercepts requested, demean before proceeding
  if (intercept == TRUE) {
    ybar <- colMeans(y)     # response means as p-dim array
    xbar <- colMeans(x)     # predictor means as pq-dim array
    y <- sweep(y, 2, ybar)  # demean responses
    x <- sweep(x, 2, xbar)  # demean predictors
    ybar <- as.matrix(ybar) # response means as p x 1 matrix
    xbar <- as.matrix(xbar) # predictor means as pq x 1 matrix
  }
  # Note: Means are stored to back out intercepts later
  # ESTIMATION
  # INITIAL STEP (k = 0)
  ups_init <- sqrt((1 / n) * crossprod(y^2, x^2)) # initial penalty loadings
  that_init <- mult_lasso(x, y, lambda_glmnet,
                          ups_init, tol_glmnet)   # and estimates
  if (post == TRUE) {                       # if refitting requested...
    refit <- mult_refit(x, y, that_init)    # refit initial estimates
    that_init <- refit$that                 # overwrite (keeping zeros)
    full_rank_post_init <- refit$full_rank  # flag full rank
  } else { # if no refitting requested...
    full_rank_post_init <- NULL # then rank deficiency is irrelevant
  }
  if (intercept == TRUE) { # if intercepts requested...
    intr_init <- ybar - that_init %*% xbar # back them out
  } else{
    intr_init <- NULL # o/w don't
  }
  # UPDATING (UP TO K TIMES)
  # Stop if relative change in penalty loadings is small enough
  # Create placeholders
  ups_refi <- array(NA, dim = c(p, pq, k))  # refined penalty loadings
  that_refi <- array(NA, dim = c(p, pq, k)) # and estimates
  rel_diffs_ups <- matrix(NA, k, 1)         # relative changes in loadings
  if (post == TRUE) { # if refitting requested...
    full_rank_post_refi <- matrix(NA, p, k) # flag full rank
  } else { # o/w don't
    full_rank_post_refi <- NULL
  }
  if (intercept == TRUE) { # if intercepts requested...
    intr_refi <- matrix(NA, p, k) # create placeholder
  } else {
    intr_refi <- NULL # o/w don't
  }
  for (l in 1:k) {                      # start updating
    if (l == 1) {                       # if first iteration...
      that_old <- that_init             # use initial estimates
    } else {
      that_old <- that_refi[, , l - 1]  # o/w use previous (refined) estimates
    }
    res_old <- y - x %*% t(that_old) # implied residuals
    ups_new <- sqrt((1 / n) * t(res_old^2) %*% x^2) # update penalty loadings
    ups_refi[, , l] <- ups_new # store
    that_new <- mult_lasso(x, y, lambda_glmnet, ups_new, tol_glmnet) # estimates
    if (post == TRUE) {                           # if refitting requested...
      refit <- mult_refit(x, y, that_new)         # refit estimates
      that_new <- refit$that                      # overwrite (keeping zeros)
      full_rank_post_refi[, l] <- refit$full_rank # flag full rank
    }
    that_refi[, , l] <- that_new                  # store
    if (intercept) {                              # if intercepts requested...
      intr_refi[, l] <- ybar - that_new %*% xbar  # back them out
    }
    if (l == 1) {           # if at first update...
      ups_old <- ups_init   # then compare w/ initial penalty loadings
    } else {                # o/w compare w/ previous (refined) penalty loadings
      ups_old <- ups_refi[, , l - 1]
    }
    diff_ups <- ups_new - ups_old # change in penalty loadings
    # relative change in vectorized ell_2 norm --v
    rel_diff_ups <- sqrt(sum(diff_ups^2)) /
      (.Machine$double.eps + sqrt(sum(ups_old^2)))
    rel_diffs_ups[l, 1] <- rel_diff_ups # store
    if (rel_diff_ups <= tol_ups) {  # if change is small enough...
      break                         # stop updating
    }
  }
  ups <- ups_refi[, , l]    # report final penalty loadings
  if (intercept == TRUE) {  # if intercepts requested...
    intr <- intr_refi[, l]  # report final intercepts
  } else {
    intr <- NULL            # o/w don't
  }
  that <- that_refi[, , l]  # report final estimates
  if (warn == TRUE && l == k && rel_diff_ups > tol_ups) {
    warning("Maximum number of updates reached. Consider increasing K.")
    cat(sprintf("Relative change in penalty loadings is %3.1g percent\n",
                100 * rel_diff_ups))
  }
  fit <- list(
    lambda = lambda_star,                             # penalty level
    ups_init = ups_init,                              # initial loadings
    ups_refi = ups_refi[, , 1:l],                     # refined loadings path
    ups = ups,                                        # final loadings
    intr_init = intr_init,                            # initial intercepts
    intr_refi = intr_refi[, 1:l],                     # refined intercepts path
    intr = intr,                                      # final intercepts
    that_init = that_init,                            # initial estimates
    that_refi = that_refi[, , 1:l],                   # refined estimates path
    that = that,                                      # final estimates
    full_rank_post_init = full_rank_post_init,        # full rank initial?
    full_rank_post_refi = full_rank_post_refi[, 1:l], # full rank refined path?
    k_term = l,                                       # no. updates performed
    rel_diff_ups_term = rel_diff_ups,                 # rel change at term
    rel_diffs_ups = rel_diffs_ups[1:l, 1]             # rel changes path
  )
  return(fit)
}