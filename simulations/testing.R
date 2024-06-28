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
n <- 10000
p <- 4
y0ton <- sim_data(n, p)

# Estimate
source("mlasso.R")
lambda_glmnet <- 0.0001
that <- mlasso(y0ton[-1, ], y0ton[-n, ], lambda_glmnet)
print(that)