# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

RNGkind(normal.kind = "Kinderman-Ramage") # faster normal draws (will do many)
seed <- 2345
set.seed(seed)

# SIMULATE DATA
n <- 100
p <- 4
source("simulations/simData.R")
nburn <- 1000
data <- sim_data_by_design(n = n, p = p, design = "Diagonal", nburn = nburn)

# ESTIMATE COEFFICIENTS
q <- 1                  # number of lags
post <- FALSE           # refit?
intercept <- FALSE      # intercept?
# LASSO (weighted)
source("lassoVAR.R")
fit_lasso <- lasso_var(data, q = q,
                       post = post, intercept = intercept)
# LASSO (equiweighted)
xy <- unpack(data, q = q)
x <- xy$x
y <- xy$y
c0 <- 1.1
pq <- p * q
gamma <- 0.1 / log(max(c(n, pq)))
lambda_star <- 2 * c0 * sqrt(n) * qnorm(1 - gamma / (2 * p * pq))
that_lasso_naive <- mult_lasso(x = x, y = y, lambda = lambda_star,
                               upsilon = NULL)

# SQRT LASSO
source("sqrtLassoVAR.R")
fit_sqrtl <- sqrt_lasso_var(data, q = q,
                            post = post, intercept = intercept)
fit_sqrtl_naive <- sqrt_lasso_var(data, q = q,
                                  post = post, intercept = intercept,
                                  upsilon = matrix(1, p, pq))