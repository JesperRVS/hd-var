sqrt_lasso_var <- function(data, q = 1, post = TRUE, intercept = TRUE,
                           c = 1.1, gamma = 0.1 / log(max(dim(data))),
                           upsilon = NULL) {
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
  lambda_star_sqrt <- c * sqrt(n) * qnorm(1 - gamma / (2 * q * p^2))
  # Note: sqrt-LASSO objective cancels out the "2" in the KKT conditions
  if (is.null(upsilon)) {
    upsilon <- matrix(1, p, pq) # default to unit loadings
  }
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
  # ****HERE HERE HERE***
}