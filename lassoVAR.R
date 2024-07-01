# ESTIMATION TOOLS FOR VECTOR AUTOREGRESSIONS WITH LASSO PENALIZATION

# Dependencies
library("glmnet")

# Multiple-equation weighted LASSO
mult_lasso <- function(x, y, lambda_glmnet, upsilon = NULL, tol_glmnet = 1e-4) {
  if (missing(x) || missing(y) || missing(lambda_glmnet)) {
    stop("Not enough input arguments.")
  }
  qp <- ncol(x)  # extract dimensions
  p <- ncol(y)   # q = autoregressive order, p = dim(output)
  if (nrow(x) != nrow(y)) {
    stop("Number of observations in x and y do not match.")
  }
  if (is.null(upsilon)) {
    upsilon <- matrix(1, p, qp)  # default to unit loadings
  }
  if (is.null(tol_glmnet)) {
    tol_glmnet <- 1e-4
  }
  # Options for glmnet
  opts <- list(standardize = FALSE,    # don't standardize (rescaling below)
               intercept = FALSE,      # don't fit an intercept
               lambda = lambda_glmnet, # penalty in eyes of glmnet
               thresh = tol_glmnet)    # tolerance for coordinate descent
  # Estimate
  that <- matrix(NA, p, qp)
  for (i in 1:p) {
    xtilde <- sweep(x, 2, upsilon[i, ], FUN = "/")  # rescale using loadings
    fit_i <- glmnet(xtilde, y[, i], family = "gaussian",
                    standardize = opts$standardize, intercept = opts$intercept,
                    lambda = opts$lambda, thresh = opts$thresh)
    that[i, ] <- as.numeric(coef(fit_i))[-1] / upsilon[i, ]  # original scaling
  }
  return(that)
}

lasso <- function(data, q = 1, post = TRUE, intr = TRUE,
                      c = 1.1, gamma = 0.1 / log(max(dim(data))), k = 15,
                      tol_ups = 1e-3, tol_glmnet = 1e-4, warn = TRUE) {
  if (missing(data) || is.null(data)) {
    stop("No data provided.")
  }
  # Unpack data
  nplusq <- nrow(data)  # n plus q
  p <- ncol(data)       # p = dim of output
  n <- nplusq - q       # n = effective sample size
  qp <- q * p           # qp = total number of variables
  y <- data[(q + 1):(q + n), ] # response
  x <- matrix(NA, nrow = n, ncol = qp) # predictors
  for (ell in 1:q) {
    # where to store ellth lag --v
    block_ell <- ((ell - 1) * p + 1):(ell * p)
    # periods corresponding to ellth lag --v
    lag_ell <- (q + 1 - ell):(q + n - ell)
    # store 1st lags first, 2nd second,...
    x[, block_ell] <- data[lag_ell, ]
  }
  # Penalty level
  lambda_star <- 2 * c * sqrt(n) * qnorm(1 - gamma / (2 * q * p^2))
  lambda_glmnet <- lambda_star / (2 * n) # penalty in eyes of glmnet
  # Note: glmnet defines "lambda" based on *half* of *average* square loss,
  # while our "lambda" stems from *sum* square loss (w/o the 1/2)
  #  If intercepts requested, demean before proceeding
  if (intr == TRUE) {
    ybar <- colMeans(y)
    y <- sweep(y, 2, ybar)
    xbar <- colMeans(x)
    x <- sweep(x, 2, xbar)
  }
  # Note: Means stored to back out intercept estimates later
  # Initial step
  ups_init <- sqrt((1 / n) * t(y^2) %*% (x^2)) # loadings
  that_init <- mult_lasso(x, y, lambda_glmnet, ups_init, tol_glmnet)
  if (post == TRUE) {
    sel <- that_init != 0 # flag active sets of predictors
    # TODO: DO SOMETHING MORE CLEVER HERE
    for (i in 1:p) {
      if (sum(sel[i, ]) > 0) { # if something selected, refit
        xx_sel_i <- t(x[, sel[i, ]]) %*% x[, sel[i, ]]
        xy_sel_i <- t(x[, sel[i, ]]) %*% y[, i]
        that_init[i, sel[i, ]] <- solve(xx_sel_i, xy_sel_i)
      }
    }
  }
  fit <- list(
                lambda = lambda_star,
                ups_init = ups_init,
                that_init = that_init)
  return(fit)
}