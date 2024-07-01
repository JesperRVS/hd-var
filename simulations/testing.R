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
source("simulations/design_A.R")
n <- 1000
p <- 4
y0ton <- sim_data(n, p)

# Estimate
source("lassoVAR.R")
q <- 1 # autoregressive order
fit_lasso <- lasso(y0ton, q = q, post = FALSE, intr = FALSE)

fit_post <- lasso(y0ton, q = q, post = TRUE, intr = FALSE)
