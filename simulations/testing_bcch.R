rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# SIMULATION SETTINGS
n <- 1000
p <- 2
seed <- 2345
set.seed(seed)
# Monte Carlo (MC) settings
RNGkind(normal.kind = "Kinderman-Ramage") # faster normal draws

# SIMULATE DATA
source("simulations/simData.R")
data <- sim_data_by_design(n = n, p = p, design = "Diagonal")

# ESTIMATE COEFFICIENTS

# Lasso
q <- 1
post <- FALSE
intercept <- FALSE
source("lassoVAR.R")
fit_lasso <- lasso_var(data = data, q = q,
                       post = post, intercept = intercept)
fit_lasso$that

# Manual Lasso
source("helper_functions.R")
xy <- unpack(data = data, q = q)
x <- xy$x
y <- xy$y
fit_bcch <- mult_lasso_bcch(x = x, y = y,
                            post = post, intercept = intercept)
fit_bcch$that

# With wrapper
fit_lasso_var_bcch <- lasso_var_bcch(data = data, q = q,
                                     post = post, intercept = intercept)
fit_lasso_var_bcch$that