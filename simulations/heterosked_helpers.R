# Define the double factorial function
double_factorial <- function(n) {
  # Check if the input is a nonnegative integer
  if (!is.numeric(n) || n < 0 || floor(n) != n) {
    stop("Input must be a nonnegative integer")
  }
  # Handle the base cases
  if (n == 0 || n == 1) {
    return(1)
  }
  # Initialize the result
  result <- 1
  # Calculate the double factorial
  while (n > 1) {
    result <- result * n
    n <- n - 2
  }
  return(result)
}

var_m <- function(m = 0) {
  if (m == 0) {
    return(1)
  } else {
    return(double_factorial(2 * m - 1) / sqrt(4 * m + 1))
  }
}

var_eps_of_x <- function(x, m = 0) {
  if (length(x) < 2) {
    stop("x must have at least 2 columns")
  }
  return(exp(- 2 * m * x[, 1]^2) * abs(x[, 2])^(2 * m))
}

sim_data_hetero <- function(n = 100, p = 4, s = 1, m = 0) {
  x <- matrix(rnorm(n * p), nrow = n, ncol = p) # standard normals
  z <- matrix(rnorm(n), nrow = n, ncol = 1)     # standard normals
  var_eps_x <- var_eps_of_x(x, m = m)
  eps <- sqrt(var_eps_x) * z    # eps | X = x ~ N(0, sigma_eps^2(x))
  var_m <- var_m(m = m)         # variance as function of m
  eps <- eps / sqrt(var_m)      # rescale to eps | X = x ~ N(0, 1)
  theta <- numeric(p)
  theta[1:s] <- 1
  if (s < p) {
    theta[(1 + s):p] <- 0
  }
  theta <- as.matrix(theta)
  y <- x %*% theta + eps
  return(list(x = x, y = y, eps = eps))
}