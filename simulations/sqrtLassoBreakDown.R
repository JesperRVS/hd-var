# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

RNGkind(normal.kind = "Kinderman-Ramage") # faster normal draws (will do many)
seed <- 2345
set.seed(seed)

# SIMULATE DATA
n <- 1000
p <- 2
s <- 1

sim_data_hetero <- function(n = 100, p = 4, s = 1, m = 0){
  x <- matrix(rnorm(n * p), nrow = n, ncol = p) # standard normals
  z <- matrix(rnorm(n), nrow = n, ncol = 1)     # standard normals
  eps <- sqrt(abs(x[, 2])^m) * z                # eps | X = x ~ N(0, |x_2|^m)
  theta <- numeric(p)
  theta[1:s] <- 1
  theta[(s + 1):p] <- 0
  theta <- as.matrix(theta)
  y <- x %*% theta + eps
  return(list(x = x, y = y, eps = eps))
}

# ESTIMATE COEFFICIENTS
source("helper_functions.R")
c0 <- 1.1
gamma <- 0.05
lambda_sqrtl <- c0 * sqrt(n) * qnorm(1 - gamma / (2 * p))
lambda_lasso <- 2 * c0 * sqrt(n) * qnorm(1 - gamma / (2 * p))

mvec <- seq(0, 10, by = 1)
numm <- length(mvec)
thats <- array(NA, dim = c(p, numm, 3))
dimnames(thats) <- list(j = 1:p, m = mvec,
                        method = c("SqrtLasso", "LassoInit", "LassoIdeal"))
for (thism in 1:numm) {
  m <- mvec[thism]
  data <- sim_data_hetero(n = n, p = p, s = s, m = m)
  that_sqrtl <- sqrt_lasso(x = data$x, y =  data$y, lambda = lambda_sqrtl)
  thats[, thism, 1] <- that_sqrtl
  ups_init <- sqrt((1 / n) * crossprod(data$y^2, data$x^2))
  xtil_init <- sweep(data$x, 2, ups_init, "/")
  fit_init_temp <- glmnet(x = xtil_init, y = data$y,
                          lambda = lambda_lasso / (2 * n),
                          standardize = FALSE, intercept = FALSE)
  that_init <- coef(fit_init_temp)[-1] / ups_init
  thats[, thism, 2] <- that_init
#   that_lasso <- mult_lasso(x = data$x, y = data$y,
#                            lambda_glmnet = lambda_lasso / (2 * n),
#                            upsilon = ups_init) / ups_init
#   thats[, thism, 2] <- as.matrix(that_lasso)
  ups_ideal <- sqrt((1 / n) * crossprod(data$eps^2, data$x^2))
  that_ideal <- mult_lasso(x = data$x, y = data$y,
                           lambda_glmnet = lambda_lasso / (2 * n),
                           upsilon = ups_ideal) / ups_ideal
  thats[, thism, 3] <- as.matrix(that_ideal)
}
thats[1,,]
colMeans(thats[1, , ])

## SANDBOX ##
data <- sim_data_hetero(n = n, p = p, s = s, m = 4)
x <- data$x
y <- data$y
ups_init <- sqrt((1 / n) * crossprod(y^2, x^2))
xtil_init <- sweep(data$x, 2, ups_init, FUN = "/")
fit_init_temp <- glmnet(x = xtil_init, y = y, family = "gaussian",
                        lambda = lambda_lasso / (2 * n),
                        standardize = FALSE, intercept = FALSE,
                        thresh = 1e-4)
that_init <- as.matrix(coef(fit_init_temp)[-1] / ups_init)

that_lasso <- mult_lasso(x, y,
                         lambda_glmnet = lambda_lasso / (2 * n),
                         upsilon = ups_init) / ups_init
# that_init <- mult_lasso(x = x, y = y, lambda_glmnet = lambda_lasso / (2 * n),
#                         upsilon = ups_init) / ups_init

## TODO ##
# Q: that_init and that_lasso are *not* the same - why?


# eps <- data$eps
# ups_ideal <- sqrt((1 / n) * crossprod(eps^2, x^2))
# that_ideal <- mult_lasso(x = x, y = y, lambda_glmnet = lambda_lasso  / (2 * n),
#                          upsilon = ups_ideal) / ups_ideal
# glmnet(x, y, family = "gaussian",
#                      lambda = lambda_glmnet,
#                      standardize = FALSE, intercept = FALSE)
