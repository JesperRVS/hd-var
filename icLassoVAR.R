## Dependencies
library("glmnet", "MASS", "Matrix")
# Notes:
#   glmnet is used for LASSO estimation
#   MASS is used for Moore-Penrose pseudoinverse (ginv)
#   Matrix is used for rankMatrix
source("helper_functions.R")  # for ic_lasso and refitting

# Function to estimate a VAR w/ LASSO and information criteria selection
# INPUTS
#   data:       (q + n) x p matrix of data (n = effective sample size)
#   q:          autoregressive order; default is 1 (i.e. VAR(1) model)
#   criteria:   information criteria to use; default is AIC, BIC and HQIC
#   post:       logical; if TRUE, refit estimates after selection
#   intercept:  logical; if TRUE, back out intercepts
#   ...:        additional arguments to ic_lasso (passed on to glmnet)
# OUTPUTs fit: list with the following components
#   intrs:      p x num_crit matrix of intercepts
#   thats:      p x pq x num_crit array of estimates
ic_lasso_var <- function(data, q = 1, criteria = c("aic", "bic", "hqic"),
                         post = TRUE, intercept = TRUE, ...) {
  if (missing(data) || is.null(data)) {
    stop("No data provided.")
  }
  criteria <- match.arg(criteria, several.ok = TRUE)
  # Unpack data and construct predictors and response
  nplusq <- nrow(data)      # n plus q
  p <- ncol(data)           # p = dim of output
  pq <- p * q               # pq = total number of variables
  xy <- unpack(data, q = q) # unpack data
  x <- xy$x                 # predictors
  y <- xy$y                 # responses
  num_crit <- length(criteria) # num information criteria
  # if intercept requested, demean before proceeding
  if (intercept == TRUE) {
    ybar <- colMeans(y)     # response means as p-dim array
    xbar <- colMeans(x)     # predictor means as pq-dim array
    y <- sweep(y, 2, ybar)  # demean responses
    x <- sweep(x, 2, xbar)  # demean predictors
    ybar <- as.matrix(ybar) # response means as p x 1 matrix
    xbar <- as.matrix(xbar) # predictor means as pq x 1 matrix
    intrs <- array(NA, dim = c(p, num_crit)) # placeholder intercepts
    dimnames(intrs) <- list(1:p, criteria)
  } else {
    intrs <- NULL
  }
  # Note: Means are stored to back out intercepts later
  thats <- array(NA, dim = c(p, pq, num_crit))  # placeholder slopes
  dimnames(thats) <- list(1:p, 1:pq, criteria)
  for (i in 1:p) { # equation-by-equation information criteria LASSO
    fit_ic <- ic_lasso(x, y[, i], criteria = criteria, ...)
    thats[i, , ] <- fit_ic$betas # store slopes
  }
  # Refit estimates for each information criterion
  if (post == TRUE) {
    full_rank_post <- array(NA, dim = c(p, num_crit)) # rank check placeholder
    dimnames(full_rank_post) <- list(1:p, criteria)
    for (crit in criteria) {
      refit <- mult_refit(x, y, thats[, , crit])  # refit selection
      thats[, , crit] <- refit$that               # overwrite (keeping zeros)
      refit
      full_rank_post[, crit] <- refit$full_rank   # flag full rank
    }
  } else {
    full_rank_post <- NULL
  }
  if (intercept == TRUE) {                    # if intercepts requested...
    for (crit in criteria) {
      intrs[, crit] <- ybar - thats[, , crit] %*% xbar # ... back them out
    }
  }
  return(list(intrs = intrs, thats = thats, full_rank_post = full_rank_post))
}