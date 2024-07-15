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
#  ...:         additional arguments to pass to glmnet
# OUTPUT: list with the following components:
#   betas: matrix of betas, with columns named after the criteria
#   lambdas: vector of lambdas that minimize each criterion
ic_lasso <- function(x, y, criteria = c("aic", "bic", "hqic"), ...) {
  criteria <- match.arg(criteria, several.ok = TRUE)
  # Fit the model using glmnet
  fit <- glmnet(x, y, family = "gaussian",
                intercept = FALSE, standardize = FALSE, ...)
  # Note: Intercept is handled via prior demeaning
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
    lambdas[crit] <- lambda_ic      # stor penalty (in eyes of glmnet)
  }
  return(list(betas = betas, lambdas = lambdas))
}