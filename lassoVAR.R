### ESTIMATION TOOLS FOR VECTOR AUTOREGRESSIONS WITH LASSO PENALIZATION

## Dependencies
library("glmnet", "MASS", "Matrix")
# Notes:
#   glmnet is used for LASSO estimation
#   MASS is used for Moore-Penrose pseudoinverse (ginv)
#   Matrix is used for rankMatrix

## MAIN FUNCTION
lasso <- function(data, q = 1, post = TRUE, intercept = TRUE,
                  c = 1.1, gamma = 0.1 / log(max(dim(data))), k = 15,
                  tol_ups = 1e-3,  warn = TRUE,
                  full_path = FALSE, tol_glmnet = 1e-4) {
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
  if (missing(data) || is.null(data)) {
    stop("No data provided.")
  }
  # INITIALIZE
  # Unpack data and construct predictors and response
  nplusq <- nrow(data)  # n + q
  p <- ncol(data)       # p = dim of output
  n <- nplusq - q       # n = effective sample size
  pq <- p*q             # pq = total number of variables
  y <- data[(q + 1):(q + n), ] # response
  x <- matrix(NA, nrow = n, ncol = pq) # predictors
  for (ell in 1:q) {
    # where to store ellth lag --v
    block_ell <- ((ell - 1) * p + 1):(ell * p)
    # periods corresponding to ellth lag --v
    lag_ell <- (q + 1 - ell):(q + n - ell)
    # store 1st lags first, 2nd second,... --v
    x[, block_ell] <- data[lag_ell, ]
  }
  # Penalty level
  lambda_star <- 2 * c * sqrt(n) * qnorm(1 - gamma / (2 * q * p^2))
  lambda_glmnet <- lambda_star / (2 * n) # penalty in eyes of glmnet
  # Note: glmnet defines "lambda" based on *half* of *average* square loss,
  # while our "lambda" stems from *sum* square loss (w/o the 1/2)
  #  If intercepts requested, demean before proceeding
  if (intercept == TRUE) {
    ybar <- colMeans(y)     # response means as p-dim array
    xbar <- colMeans(x)     # predictor means as pq-dim array
    y <- sweep(y, 2, ybar)  # demean response
    x <- sweep(x, 2, xbar)  # demean predictors
    ybar <- as.matrix(ybar) # response means as p x 1 matrix
    xbar <- as.matrix(xbar) # predictor means as pq x 1 matrix
  }
  # Note: Means are stored to back out intercepts later
  # ESTIMATION
  # INITIAL STEP (k=0)
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

## HELPER FUNCTIONS

# Multiple responses weighted LASSO using same regressors
mult_lasso <- function(x, y, lambda_glmnet, upsilon = NULL,
                       full_path = FALSE, tol_glmnet = 1e-4) {
  # INPUTS
  #   x: n x pq matrix of predictors
  #   y: n x p matrix of responses
  #   lambda_glmnet: penalty level in eyes of glmnet
  #   upsilon: p x pq matrix of penalty loadings
  #   full_path: logical; if TRUE, use full path of lambda values
  #   tol_glmnet: convergence tolerance for glmnet
  # OUTPUT
  #   that: p x pq matrix of estimates
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
  return(that)
}

# Find solution to minimize ||y - Xb||_{\ell_2} over b
ls_sol <- function(x, y) {
  # INPUTS
  #   x:          n x p matrix of regressors
  #   y:          n x 1 vector of responses
  # OUTPUTs fit: list with the following components
  #   sol:        p x 1 vector of OLS estimates
  #   full_rank:  logical; if TRUE, regressors of full rank
  full_rank <- ncol(x) == Matrix::rankMatrix(x)[1] # check for full column rank
  if (full_rank == TRUE) { # if full rank, produce unique (OLS) solution
    sol <- solve(crossprod(x), crossprod(x, y))
  } else { # o/w produce solution based on pseudoinverse
    sol <- MASS::ginv(x) %*% y # ginv is Moore-Penrose pseudoinverse
  }
  fit <- list(sol = sol, full_rank = full_rank)
  return(fit)
}

# Least squares refitting with multiple responses and same regressors
mult_refit <- function(x, y, that) {
  # INPUTS
  #   x: n x pq matrix of predictors
  #   y: n x p matrix of responses
  #   that: p x pq matrix of estimates
  # OUTPUTs refit: list with the following components
  #   that: p x pq matrix of refitted estimates
  #   full_rank: p x 1 logical vector; if TRUE, selected regressors of full rank
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


