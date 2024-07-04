# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages
libest <- c("glmnet")
libplt <- c("ggplot2") # c("formattable", "ggplot2", "latex2exp", "reshape", "stringr")
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
using(append(libest, libplt))

# Simulate data
seed <- 2345
r <- 0
source("simulations/design_A.R")
n <- 100
p <- 4
y0ton <- sim_data(n, p, seed = seed, r = r)

# Estimate
source("lassoVAR.R")
q <- 1 # autoregressive order
fit_lasso <- lasso(y0ton, q = q, post = FALSE, intercept = FALSE)
refit <- mult_refit(y0ton[1:n, ], y0ton[-1, ], that = fit_lasso$that_init)

fit_lasso_intr <- lasso(y0ton, q = q, post = FALSE, intercept = TRUE)
fit_post <- lasso(y0ton, q = q, post = TRUE, intercept = FALSE)
fit_post_intr <- lasso(y0ton, q = q, post = TRUE, intercept = TRUE)
# fit_post_2 <- lasso(y0ton, q = q, gamma = p, post = TRUE, intercept = TRUE)

fit_ls <- lasso(y0ton, q = q, c = 0, post = FALSE, intercept = FALSE)
fit_ls_intr <- lasso(y0ton, q = q, c = 0, post = FALSE, intercept = TRUE)
fit_ls_post <- lasso(y0ton, q = q, c = 0, post = TRUE, intercept = FALSE)
fit_ls_post_intr <- lasso(y0ton, q = q, c = 0, post = TRUE, intercept = TRUE)

x <- y0ton[1:n, ]
y <- y0ton[-1, ]

ybar <- as.matrix(colMeans(y))

a <- c(1:9)
b <- rev(a)

#array of data
ab <- array(c(a,b), dim = c(3,3,2))

# Plotting

# Illustrating penalty loading convergence
plot(1:fit_lasso$k_term, fit_lasso$rel_diffs_ups, type = "l",
    xlab = "Number of updates, k", ylab = "Relative change in penalty loadings (vectorized ell_2 norm)",
    main = "LASSO: Relative change in penalty loadings"
)


library(ggplot2)
ggplot(data = data.frame(x = 1:fit_lasso$k_term, y = fit_lasso$rel_diffs_ups), aes(x = x, y = y)) +
    geom_line() +
    labs(x = "Number of updates, k", y = "Relative change in penalty loadings (vectorized ell_2 norm)",
         title = "LASSO: Relative change in penalty loadings"
    )
