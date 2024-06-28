# mlasso: Multiple equation LASSO

# Dependencies
library("glmnet")

mlasso <- function(x, y, lambda_glmnet, upsilon = NULL, tol_glmnet = 1e-4) {
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
  opts <- list(standardize = FALSE,    # don't standardize (normalization below)
               intercept = FALSE,      # don't fit an intercept
               lambda = lambda_glmnet, # penalty in eyes of glmnet
               thresh = tol_glmnet)    # tolerance for coordinate descent
  # Estimate
  that <- matrix(NA, p, qp)
  for (i in 1:p) {
    xtilde <- sweep(x, 2, upsilon[i, ], FUN = "/")  # rescale
    fit_i <- glmnet(xtilde, y[, i], family = "gaussian",
                    standardize = opts$standardize, intercept = opts$intercept,
                    lambda = opts$lambda, thresh = opts$thresh)
    that[i, ] <- as.numeric(coef(fit_i))[-1] / upsilon[i, ]  # original scaling
  }
  return(that)
}
