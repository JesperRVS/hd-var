## TODO:
# 0. Create following functions:
    # [x] AIC-Lasso,
    # [x] Post-AIC-Lasso,
    # [x] BIC-Lasso,
    # [x] Post-BIC-Lasso,
    # [x] Sqrt-Lasso,
    # [x] Post-Sqrt-Lasso
# 1. Create main_sim file
# 2. Use 5 designs 
#   (1) "Diagonal":     Design A as in KC2015
#   (2) "Correlated":   Design A' w/ strongly correlated innovations (hence outcomes)
#   (3) "HeavyTailed":  Design A'' w/ heavy-tailed (here: student-t(5) innovations)
#   (4) "BlockDiag":    Design B as in KC2015
#   (5) "NearBand":     Design C as in KC2015
# 3. Use the following estimation methods (1) Lasso w/ Lasso updating (2)
#   Post-Lasso w/ post-Lasso updating (3) AIC-Lasso (4) Post-AIC-Lasso (5)
#   BIC-Lasso (6) Post-BIC-Lasso (7) Sqrt-Lasso (8) Post-Sqrt-Lasso (9) Least
#   squares w/ Moore-Penrose inverse
# 4. See if you can store all matrices

# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages
libest <- c("glmnet")
libplt <- c("ggpubr", "reshape2", "latex2exp")
# c("formattable", "latex2exp", "stringr")
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

# SIMULATION
RNGkind(normal.kind = "Kinderman-Ramage") # faster normal draws (will do many)
seed <- 2345
r <- 0
source("simulations/design_A.R")
n <- 1000
p <- 4

y0ton <- sim_data(n, p, seed = seed, r = r)

ydata <- y0ton

## ESTIMATION VIA IC-LASSO
source("icLassoVAR.R")
fit_ic_naive <- ic_lasso_var(ydata, criteria = c("aic", "bic", "hqic"),
    post = FALSE, intercept = FALSE)
x <- ydata[1:n, ]
stdx <- sqrt(colMeans(x^2))
upsilon <- matrix(stdx, nrow = p, ncol = p, byrow = TRUE)
fit_ic <- ic_lasso_var(ydata, criteria = c("aic", "bic", "hqic"),
    post = FALSE, intercept = FALSE, upsilon = upsilon)
fit_ic_naive
fit_ic
fit_ic_naive_intrs <- ic_lasso_var(ydata, criteria = c("aic", "bic", "hqic"),
    post = FALSE, intercept = TRUE)
fit_ic_naive_post <- ic_lasso_var(ydata, criteria = c("aic", "bic", "hqic"),
    post = TRUE, intercept = FALSE)
fit_ic_naive_post_intrs <- ic_lasso_var(ydata, criteria = c("aic", "bic", "hqic"),
    post = TRUE, intercept = TRUE)
##
xtemp <- matrix(rnorm(n * p), nrow = n, ncol = p)
sigma_p <- 10
theta_p <- 1 / sigma_p
x <- cbind(xtemp[, 1:(p - 1)], sigma_p * xtemp[, p])
y <- 1 * x[, 1] + theta_p * x[, p] + sqrt(10) * rnorm(n)

fit <- glmnet(x, y, family = "gaussian",
              intercept = FALSE, standardize = TRUE)

source("helper_functions.R")
fit_ic_naive <- ic_lasso(x, y, criteria = c("aic", "bic", "hqic"))
stdx <- sqrt(colMeans(x^2))
xtilde <- sweep(x, 2, stdx, FUN = "/")
fit_ic_proper <- ic_lasso(xtilde, y, criteria = c("aic", "bic", "hqic"))
sweep(fit_ic_proper$betas, 1, stdx, "/")

fit_ic

# ESTIMATION VIA INFORMATION CRITERIA
source("icLassoVAR.R")
fit_ic_naive <- ic_lasso_var(ydata, criteria = c("aic", "bic", "hqic"),
    post = FALSE, intercept = FALSE, standardize = FALSE)
fit_ic <- ic_lasso_var(ydata, criteria = c("aic", "bic", "hqic"),
    post = FALSE, intercept = FALSE, standardize = TRUE)
fit_ic_intr <- ic_lasso_var(ydata, criteria = c("aic", "bic", "hqic"),
    post = FALSE, intercept = TRUE, standardize = TRUE)
fit_post_ic <- ic_lasso_var(ydata, criteria = c("aic", "bic", "hqic"),
    post = TRUE, intercept = FALSE, standardize = TRUE)
fit_post_ic_intr <- ic_lasso_var(ydata, criteria = c("aic", "bic", "hqic"),
    post = TRUE, intercept = TRUE, standardize = TRUE)



## OLD BELOW THIS LINE ##

# ESTIMATION VIA SQRT-LASSO
q <- 1
source("sqrtLassoVAR.R")
fit_sqrt_lasso_naive <- sqrt_lasso_var(ydata, q = q, post = FALSE, intercept = FALSE)
x <- ydata[1:n, ]
fit_sqrt_lasso <- sqrt_lasso_var(ydata, q = q, post = FALSE, intercept = FALSE,
    upsilon = matrix(sqrt(colMeans(x^2)), nrow = p, ncol = p * q))
fit_post_sqrt_lasso <- sqrt_lasso_var(ydata, q = q, post = TRUE, intercept = FALSE,
    upsilon = matrix(sqrt(colMeans(x^2)), nrow = p, ncol = p * q))
fit_sqrt_lasso_intr <- sqrt_lasso_var(ydata, q = q, post = FALSE, intercept = TRUE,
    upsilon = matrix(sqrt(colMeans(x^2)), nrow = p, ncol = p * q))
fit_post_sqrt_lasso_intr <- sqrt_lasso_var(ydata, q = q, post = TRUE, intercept = TRUE,
    upsilon = matrix(sqrt(colMeans(x^2)), nrow = p, ncol = p * q))

## OLD BELOW
source("helper_functions.R")
pq <- p * q
xy <- unpack(ydata, q = q)
x <- xy$x
y <- xy$y




c <- 1.1
gamma <- 0.1 / log(max(c(n, p)))
lambda_star <- c * sqrt(n) * qnorm(1 - gamma / (2 * p * pq))
upsilon <- matrix(sqrt(colMeans(x^2)), nrow = p, ncol = pq)
# upsilon <- matrix(1, nrow = p, ncol = ncol(x))
thats_sqrt_lasso <- mult_sqrt_lasso(x, y, lambda = lambda_star,
    upsilon = upsilon)
thats_sqrt_lasso
# thats_sqrt_lasso <- mult_sqrt_lasso(x, y, lambda = lambda_star,
#     upsilon = upsilon, print_out = FALSE, max_iter = 100)

x <- ydata[1:n, ]
y <- ydata[2:(1 + n), 1]
sqrt_lasso(x, y, lambda = lambda_star)

## OLD BELOW THIS LINE ##

# ESTIMATION VIA SQUARE-ROOT LASSO
source("helper_functions.R")
x <- ydata[1:n, ]
y <- ydata[2:(1 + n), 1]
c <- 1.1
gamma <- 0.1 / log(max(c(n, p)))
lambda_star <- c * sqrt(n) * qnorm(1 - gamma / (2 * p^2))
upsilon <- sqrt(colMeans(x^2))
# upsilon <- rep(1, p)
xtilde <- sweep(x, 2, upsilon, FUN = "/")
that_sqrt <- sqrt_lasso(xtilde, y, lambda = lambda_star, print_out = TRUE, max_iter = 100) / upsilon

that_sqrt_naive <- sqrt_lasso(x, y, lambda = lambda_star, print_out = TRUE, max_iter = 100)

# ESTIMATION VIA INFORMATION CRITERIA
source("icLassoVAR.R")
fit_ic <- ic_lasso_var(ydata, criteria = c("aic", "bic", "hqic"),
    post = FALSE, intercept = FALSE)
fit_ic_intr <- ic_lasso_var(ydata, criteria = c("aic", "bic", "hqic"),
    post = FALSE, intercept = TRUE)
fit_post_ic <- ic_lasso_var(ydata, criteria = c("aic", "bic", "hqic"),
    post = TRUE, intercept = FALSE)
fit_post_ic_intr <- ic_lasso_var(ydata, criteria = c("aic", "bic", "hqic"),
    post = TRUE, intercept = TRUE)

# ESTIMATION
source("lassoVAR.R")
fit_lasso <- lasso_var(ydata, q = 1, post = FALSE, intercept = FALSE)
fit_postl <- lasso_var(ydata, q = 1, post = TRUE, intercept = FALSE)




x <- ydata[1:n, ]
y <- ydata[2:(1 + n), 1]

source("icGlmnet.R")
x <- ydata[1:n, ]
y <- ydata[2:(1 + n), 2]
fit_ic <- ic_glmnet(x, y, criteria = c("aic", "bic", "hqic"),
    intercept = FALSE, standardize = FALSE)


fit_aic <- ic_glmnet(x, y,
    glmnet_options = list(standardize = FALSE, thresh = 1e-8),
    criterion = "bic")

fit <- glmnet(x, y, family = "gaussian", intercept = FALSE, standardize = FALSE)


# ESTIMATION
source("lassoVAR.R")
q <- 1 # autoregressive order
fit_lasso <- lasso(ydata, q = q, post = FALSE, intercept = FALSE,
                   tol_ups = .Machine$double.eps) # minimal tolerance
fit_postl <- lasso(ydata, q = q, post = TRUE, intercept = FALSE,
                   tol_ups = .Machine$double.eps) # minimal tolerance
# ^-- here actually met b/c no change in selection


y0ton <- sim_data(n, p, seed = seed, r = r,
    family = "gaussian", rho = .9)

# y0ton <- sim_data(n, p, seed = seed, r = r, family = "student", df = 5)

# use ggplot to plot all p times in one figure with different colors
time <- 1:(1 + n)
df <- melt(data.frame(cbind(y0ton, time)), id.vars = "time")
ggplot(data = df, aes(x = time, y = value, color = variable)) +
    geom_line() +
    labs(x = "Time", y = "Value", title = "Time series plot of all variables")






rho <- .5
sigma_eps <- 0.1
cor_eps <- cor_mat(p, rho)       # innovation covariance matrix
n_burn <- 10000
n_tot <- 1 + n + n_burn                     # total no. of periods
normals <- matrix(rnorm(n_tot * p), n_tot, p)
eps <- sigma_eps * normals %*% chol(cor_eps)






# y <- y0ton[2:(1 + n), ]
# x <- y0ton[1:n, ]
# source("lassoVAR.R")
# lambda <- 0.001
# library("rbenchmark")
# benchmark("elementwise" = {
#             lasso_elem <- mult_lasso(x, y, lambda_glmnet = lambda, full_path = FALSE)
#           },
#           "full_path" = {
#             lasso_path <- mult_lasso(x, y, lambda_glmnet = lambda, full_path = TRUE)
#           },
#           replications = 1000, order = "relative", 
#           columns = c("test", "replications", "elapsed", "relative")
#           )

# lasso_path <- mult_lasso(x, y, lambda_glmnet = lambda, full_path = TRUE)
# lasso_elem <- mult_lasso(x, y, lambda_glmnet = lambda, full_path = FALSE)
# lambda <- 0.001
# lasso_elem <- glmnet(x, y, standardize = FALSE, intercept = FALSE)
# lasso_coef <- coef(lasso_elem, s = lambda)[-1]

# lasso_coef

# ESTIMATION
source("lassoVAR.R")
q <- 1 # autoregressive order
fit_lasso <- lasso(y0ton, q = q, post = FALSE, intercept = FALSE,
                   tol_ups = .Machine$double.eps)

# PLOTTING
# Illustrating convergence as a function of the number of updates (k) using LASSO
rel_diffs_that <- matrix(NA, 15, 1)
theta <- coef_matrix(p)
errs_lasso_init <- fit_lasso$that_init - theta
errs_lasso_refi <- sweep(fit_lasso$that_refi, c(1, 2), theta, "-")
ell2_vec_errs_init <- sqrt(sum(errs_lasso_init^2))
ell2_vec_errs_refi <- apply(errs_lasso_refi, 3, function(x) sqrt(sum(x^2)))
ell2_vec_errs <- c(ell2_vec_errs_init, ell2_vec_errs_refi)
for (l in 1:15) {
    if (l == 1) {
        that_old <- fit_lasso$that_init
    } else {
        that_old <- fit_lasso$that_refi[, , l - 1]
    }
    that_new <- fit_lasso$that_refi[, , l]
    diff_that <- that_new - that_old
    rel_diff_that <- sqrt(sum(diff_that^2)) /
        (.Machine$double.eps + sqrt(sum(that_old^2)))
    rel_diffs_that[l] <- rel_diff_that
}
# Create a single plot of the relative change in penalty loadings as a function
# of the number of updates with the y-axis in percent
df_ups <- data.frame(k = 1:15, rel_diffs_ups = fit_lasso$rel_diffs_ups)
p_ups <- ggplot(data = df_ups, aes(x = k, y = rel_diffs_ups)) +
    geom_line() +
    geom_point() +
    labs(x = TeX("Update $k$"), y = "",
         title = TeX("(a) $\\frac{\\|vec(\\{\\widehat{\\upsilon}^{(k)}_{i,j}-\\widehat{\\upsilon}^{(k-1)}_{i,j}\\}_{(i,j)})\\|_{l_{2}}}{\\|vec(\\{\\widehat{\\upsilon}^{(k-1)}_{i,j}\\}_{(i,j)})\\|_{l_{2}}}$")
    ) +
    scale_x_continuous(breaks = c(1, 5, 10, 15)) +
    scale_y_continuous(labels = scales::percent_format(scale = 100)) +
    theme(plot.title = element_text(size = 10))
# Create a single plot of the relative change in estimates as a function
# of the number of updates with the y-axis in percent
df_that <- data.frame(k = 1:15, rel_diffs_that = rel_diffs_that)
p_that <- ggplot(data = df_that, aes(x = k, y = rel_diffs_that)) +
    geom_line() +
    geom_point() +
    labs(x = TeX("Update $k$"), y = "",
        title = TeX("(b) $\\frac{\\|vec(\\{\\widehat{\\beta}^{(k)}_{i,j}-\\widehat{\\beta}^{(k-1)}_{i,j}\\}_{(i,j)})\\|_{l_{2}}}{\\|vec(\\{\\widehat{\\beta}^{(k-1)}_{i,j}\\}_{(i,j)})\\|_{l_{2}}}$")
         ) +
    scale_x_continuous(breaks = c(1, 5, 10, 15)) +
    scale_y_continuous(labels = scales::percent_format(scale = 100)) +
    theme(plot.title = element_text(size = 10))
# title = TeX("(b) $\\|vec(\\widehat{\\Theta}^{(k)}-\\widehat{\\Theta}^{(k-1)})\\|_{l_{2}}/\\|vec(\\widehat{\\Theta}^{(k-1)})\\|_{l_{2}}$")

# Create a single plot of the ell_2 norm of vectorized errors as a function
# of the number of updates
df_errs <- data.frame(k = 0:15, ell2_vec_errs = ell2_vec_errs)
p_errs <- ggplot(data = df_errs, aes(x = k, y = ell2_vec_errs)) +
    geom_line() +
    geom_point() +
    labs(x = TeX("Update $k$"), y = "",
         title = TeX("(c) $\\|vec(\\{\\widehat{\\beta}^{(k)}_{i,j}-\\beta_{0i,j}\\}_{(i,j)})\\|_{l_{2}}$")
         ) +
    scale_x_continuous(breaks = c(0, 1, 5, 10, 15)) +
    scale_y_continuous(breaks = seq(.75, 1, by = .05)) +
    theme(plot.title = element_text(size = 10))
ggarrange(p_ups, p_that, p_errs, nrow = 1)
ggsave("Figure_lasso_convergence_updating.pdf", width = 8, height = 3, units = "in", dpi = 300)



# time <- 1:(1 + n)
# df <- melt(data.frame(cbind(y0ton, time)), id.vars = "time")
# ggplot(data = df, aes(x = time, y = value, color = variable)) +
#     geom_line() +
#     labs(x = "Time", y = "Value", title = "Time series plot of all variables")

# ggplot(data = data.frame(x = 1:(1 + n), y = y0ton[, 1]), aes(x = x, y = y)) +
#     geom_line() +
#     labs(x = "Time", y = "Value", title = "Time series plot of first variable")
# ggplot(data = data.frame(x = 1:(1 + n), y = y0ton), aes(x = x, y = y)) +
#     geom_line() +
#     labs(x = "Time", y = "Value", title = "Time series plot of second variable")




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

# PLOTTING

# Illustrating penalty loading convergence as a function of the number of
# updates Lasso w/ Lasso loading updates w/ impossible tolerance (iteration
# reaches limit K=15)
fit_lasso <- lasso(y0ton, q = q, post = FALSE, intercept = FALSE,
                   tol_ups = .Machine$double.eps)
rel_diffs_that <- matrix(NA, 15, 1)
theta <- coef_matrix(p)
errs_lasso <- sweep(fit_lasso$that_refi, c(1, 2), theta, "-")
ell2_vec_errs <- apply(errs_lasso, 3, function(x) sqrt(sum(x^2)))
for (l in 1:15) {
    if (l == 1) {
        that_old <- fit_lasso$that_init
    } else {
        that_old <- fit_lasso$that_refi[, , l - 1]
    }
    that_new <- fit_lasso$that_refi[, , l]
    diff_that <- that_new - that_old
    rel_diff_that <- sqrt(sum(diff_that^2)) /
        (.Machine$double.eps + sqrt(sum(that_old^2)))
    rel_diffs_that[l] <- rel_diff_that
}
# df <- data.frame(k = 1:15, rel_diffs_ups = fit_lasso$rel_diffs_ups,
#     rel_diffs_that = rel_diffs_that, ell2_vec_errs = ell2_vec_errs)

# Gather these three plots side by side using ggplot2
library(ggpubr)


grid.arrange(p_ups, p_that, p_errs, nrow = 1)



ggplot(data = df, aes(x = k, y = value, color = variable)) +
    geom_line() +
    geom_point() +
    labs(x = "Number of updates, k", y = "",
         title = "LASSO: Relative change in penalty loadings and estimates"
    ) +
    scale_x_continuous(breaks = c(1, 5, 10, 15))




# Plotting ell_2 norm of vectorized errors as function of updates
df <- data.frame(k = 1:15, ell2_vec_errs = ell2_vec_errs)
df <- melt(df, id.vars = "k")
ggplot(data = df, aes(x = k, y = value)) +
    geom_line() +
    geom_point() +
    labs(x = "Number of updates, k", y = "ell_2 norm of vectorized errors",
         title = "LASSO: ell_2 norm of vectorized errors"
    ) +
    scale_x_continuous(breaks = c(1, 5, 10, 15))







# Create a plot with three panels side-by-side (1 row, 3 columns) with the
# following:
# 1. The relative change in penalty loadings as a function of the number of
#    updates
# 2. The relative change in estimates as a function of the number of updates
# 3. The ell_2 norm of vectorized errors as a function of the number of updates



ggplot(data = df, aes(x = k, y = value, color = variable)) +
    geom_line() +
    geom_point() +
    facet_wrap(~variable, scales = "free_y") +
    labs(x = "Number of updates, k", y = "",
         title = "LASSO: Relative change in penalty loadings and estimates"
    ) +
    scale_x_continuous(breaks = c(1, 5, 10, 15))
# provide code to format the first two panels as percentages with the third
# panel remaining a decimal number







# Plot the ell_2 norm of vectorized errors as a function of the number of updates
# and the relative change in penalty loadings and estimates, respectively, as function of updates
# with the relative change on a primary y-axis and the ell_2 norm of vectorized errors on a secondary y-axis







df <- data.frame(k = 1:15, rel_diffs_ups = fit_lasso$rel_diffs_ups,
    rel_diffs_that = rel_diffs_that)
df <- melt(df, id.vars = "k")
ggplot(data = df, aes(x = k, y = value, color = variable)) +
    geom_line() +
    geom_point() +
    labs(x = "Number of updates, k", y = "Relative change",
         title = "LASSO: Relative change in penalty loadings and estimates"
    ) +
    scale_x_continuous(breaks = c(1, 5, 10, 15))

df <- data.frame(k = 1:fit_lasso$k_term,
                 rel_diffs_ups = fit_lasso$rel_diffs_ups)
ggplot(data = df, aes(x = k, y = rel_diffs_ups)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 1e-3, linetype = "dashed", color = "red") +
  labs(x = "Number of updates, k",
    y = "Relative change in penalty loadings (vectorized ell_2 norm)",
    title = "LASSO: Relative change in penalty loadings",
    xlabel = 1:fit_lasso$k_term
  ) +
  scale_x_continuous(breaks = c(1, 5, 10, 15)) +
  scale_y_continuous(labels = scales::percent_format(scale = 100))

mat <- matrix(rep(1, 6), nrow = 2, ncol = 3)


# Illustrating penalty loadings convergence as a function of the number of updates
df <- data.frame(k = 1:fit_lasso$k_term,
    thats = t(apply(fit_lasso$that_refi, 3, diag)))
df <- melt(df, id.vars = "k")
ggplot(data = df, aes(x = k, y = value, color = variable)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of updates, k", y = "Penalty Loadings",
    title = "LASSO: Estimate convergence (diagonal entries)"
  )

# Illustrating estimate convergence as a function of the number of updates
df <- data.frame(k = 1:fit_lasso$k_term, thats = t(apply(fit_lasso$that_refi, 3, diag)))
df <- melt(df, id.vars = "k")
ggplot(data = df, aes(x = k, y = value, color = variable)) +
    geom_line() +
    geom_point() +
    labs(x = "Number of updates, k", y = "Estimate",
         title = "LASSO: Estimate convergence (diagonal entries)"
    )
    
# ggplot(data = df, aes(x = k, y = value)) +
#     geom_line() +
#     labs(x = "Number of updates, k", y = "Estimate of the first penalty loading",
#          title = "LASSO: Estimate convergence (entry [1,1])"
#     )
    




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
