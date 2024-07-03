# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages
libest <- c("glmnet")
# Auto-installer (checks if installed - if not, installs and loads)
using <- function(...) {
    libs <- unlist(list(...))
    req <- unlist(lapply(libs, require, character.only = TRUE))
    need <- libs[req == FALSE]
    if (length(need) > 0) {
        install.packages(pkgs = need, repos = "https://cran.us.r-project.org")
        lapply(need, require, character.only = TRUE)
    }
}
using(libest)

# Simulate data
seed <- 2345
r <- 0
source("simulations/design_A.R")
n <- 10
p <- 12
y0ton <- sim_data(n, p, seed = seed, r = r)

# Estimate
source("lassoVAR.R")
q <- 1 # autoregressive order
fit_lasso <- lasso(y0ton, q = q, post = FALSE, intr = FALSE)
fit_lasso_intr <- lasso(y0ton, q = q, post = FALSE, intr = TRUE)
fit_post <- lasso(y0ton, q = q, post = TRUE, intr = FALSE) # HERE HERE HERE
fit_post_intr <- lasso(y0ton, q = q, post = TRUE, intr = TRUE)
# fit_post_2 <- lasso(y0ton, q = q, gamma = p, post = TRUE, intr = TRUE)

fit_ls <- lasso(y0ton, q = q, c = 0, post = FALSE, intr = FALSE)
fit_ls_intr <- lasso(y0ton, q = q, c = 0, post = FALSE, intr = TRUE)
fit_ls_post <- lasso(y0ton, q = q, c = 0, post = TRUE, intr = FALSE)
fit_ls_post_intr <- lasso(y0ton, q = q, c = 0, post = TRUE, intr = TRUE)

x <- y0ton[1:n, ]
y <- y0ton[-1, ]

ybar <- as.matrix(colMeans(y))