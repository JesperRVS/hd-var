## Dependencies
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