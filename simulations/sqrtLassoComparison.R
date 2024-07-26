# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages for parallel computing
library("doRNG")
library("doParallel")
library("foreach")

# SIMULATION SETTINGS
n <- 1000
p <- 2
s <- 1
# mvec <- seq(0, 5, by = 1)
mvec <- seq(0, 3, by = 0.2)
numm <- length(mvec) # number of m values (= DGPs)

# ESTIMATE COEFFICIENTS
c0 <- 1.1
gamma <- 0.05
lambda_sqrtl <- c0 * sqrt(n) * qnorm(1 - gamma / (2 * p))
lambda_lasso <- 2 * c0 * sqrt(n) * qnorm(1 - gamma / (2 * p))
lambda_lasso_glmnet <- lambda_lasso / (2 * n)
methods <- c("SqrtLasso", "LassoInit", 
             "LassoOnceRefined", "LassoTwiceRefined", "LassoIdeal")
nummet <- length(methods)

# Monte Carlo (MC) settings
RNGkind(normal.kind = "Kinderman-Ramage") # faster normal draws
iseed <- 2345
cl <- makeCluster(detectCores())
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = iseed)
nummc <- 1000

# Placeholders for the results
ell2_errors <- array(NA, dim = c(nummet, nummc, numm))
dimnames(ell2_errors) <- list(method = methods, mc = 1:nummc, m = mvec)

# START PARALLEL COMPUTING
theta <- c(rep(1, s), rep(0, p - s))
for (thism in 1:numm) {
  m <- mvec[thism]
  results <- foreach(icount(nummc)) %dopar% {
    # Simulate data
    source("simulations/heterosked_helpers.R", local = TRUE)
    data <- sim_data_hetero(n = n, p = p, s = s, m = m)
    # Estimate coefficients
    source("helper_functions.R", local = TRUE)
    thats_mc <- matrix(NA, nrow = p, ncol = nummet)
    # SqrtLasso
    that_sqrtl <- sqrt_lasso(x = data$x, y =  data$y, lambda = lambda_sqrtl)
    thats_mc[, 1] <- that_sqrtl
    # LassoInit
    ups_init <- sqrt((1 / n) * crossprod(data$y^2, data$x^2))
    that_init <- mult_lasso(x = data$x, y = data$y,
                            lambda_glmnet = lambda_lasso_glmnet,
                            upsilon = ups_init)
    thats_mc[, 2] <- that_init
    # LassoOnceRefined
    res_init <- data$y - data$x %*% t(that_init)
    ups_refi <- sqrt((1 / n) * crossprod(res_init^2, data$x^2))
    that_refi <- mult_lasso(x = data$x, y = data$y,
                            lambda_glmnet = lambda_lasso_glmnet,
                            upsilon = ups_refi)
    thats_mc[, 3] <- that_refi
    # LassoTwiceRefined
    res_refi <- data$y - data$x %*% t(that_refi)
    ups_refi2 <- sqrt((1 / n) * crossprod(res_refi^2, data$x^2))
    that_refi2 <- mult_lasso(x = data$x, y = data$y,
                             lambda_glmnet = lambda_lasso_glmnet,
                             upsilon = ups_refi2)
    thats_mc[, 4] <- that_refi2
    # LassoIdeal
    ups_ideal <- sqrt((1 / n) * crossprod(data$eps^2, data$x^2))
    that_ideal <- mult_lasso(x = data$x, y = data$y,
                             lambda_glmnet = lambda_lasso_glmnet,
                             upsilon = ups_ideal)
    thats_mc[, 5] <- that_ideal
    errors <- sweep(thats_mc, 1, theta, "-")
    ell2 <- apply(errors, 2, function(x) sqrt(sum(x^2)))
    return(ell2)
  } # end mc loop
  ell2_errors[, , thism] <- array(unlist(results), dim = c(nummet, nummc))
} # end m loop
# Stop parallel computing
stopCluster(cl)

# Calculate the average ell2 errors
ell2_avg <- apply(ell2_errors, c(1, 3), mean)
# Put the average ell2 errors relative to that of SqrtLasso for each m
ell2_avg_rel <- sweep(ell2_avg[-1, ], 2, ell2_avg[1, ], "/")
# Plot the average ell2 errors of LassoInit, LassoRefined and LassoIdeal relative to that of SqrtLasso
library("ggplot2")
library("gridExtra")
df <- reshape::melt(ell2_avg_rel[-1, ])
p1 <- ggplot(df, aes(x = m, y = value, colour = method)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "m", y = "Average ell2 error relative to SqrtLasso") +
  ggtitle("Average ell2 errors relative to SqrtLasso for feasible and ideal Lassos") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black") +
  annotate("text", x = 0, y = 1.1, label = "SqrtLasso", hjust = 0, vjust = 0) +
  scale_x_continuous(breaks = mvec, minor_breaks = NULL)
p1

# save the plot as pdf
file_name <- paste0("ell2_errors_2_regressor_heterosked_sqrtLasso_comparison_n", n,
                    "_p", p, "_mean.pdf")

ggsave(file_name, plot = p1, width = 10, height = 5)

# Same plot but for relative median ell2 errors
ell2_med <- apply(ell2_errors, c(1, 3), median)
ell2_med_rel <- sweep(ell2_med[-1, ], 2, ell2_med[1, ], "/")
df <- reshape::melt(ell2_med_rel[-1, ])
p2 <- ggplot(df, aes(x = m, y = value, colour = method)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "m", y = "Median ell2 error relative to SqrtLasso") +
  ggtitle("Median ell2 errors relative to SqrtLasso for feasible and ideal Lassos") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black") +
  annotate("text", x = 0, y = 1.1, label = "SqrtLasso", hjust = 0, vjust = 0) +
  scale_x_continuous(breaks = mvec, minor_breaks = NULL)
p2

# save the plot as pdf
file_name <- paste0("ell2_errors_2_regressor_heterosked_sqrtLasso_comparison_n", n,
                    "_p", p, "_median.pdf")

ggsave(file_name, plot = p2, width = 10, height = 5)

# Same blot but for the relative maximum ell2 errors
ell2_max <- apply(ell2_errors, c(1, 3), max)

ell2_max_rel <- sweep(ell2_max[-1, ], 2, ell2_max[1, ], "/")
df <- reshape::melt(ell2_max_rel[-1, ])
p3 <- ggplot(df, aes(x = m, y = value, colour = method)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "m", y = "Maximum ell2 error relative to SqrtLasso") +
  ggtitle("Maximum ell2 errors relative to SqrtLasso for feasible and ideal Lassos") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black") +
  annotate("text", x = 0, y = 1.1, label = "SqrtLasso", hjust = 0, vjust = 0) +
  scale_x_continuous(breaks = mvec, minor_breaks = NULL)
p3

# save the plot as pdf
file_name <- paste0("ell2_errors_2_regressor_heterosked_sqrtLasso_comparison_n", n,
                    "_p", p, "_maximum.pdf")

ggsave(file_name, plot = p3, width = 10, height = 5)

# # Plot the average ell2 errors of LassoIdeal alone relative to that of SqrtLasso
# library("ggplot2")
# library("gridExtra")
# df <- reshape::melt(ell2_avg_rel)
# p1 <- ggplot(df, aes(x = m, y = value)) +
#   geom_line(color = "blue") +
#   geom_point(color = "blue") +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   labs(x = "m", y = "Average ell2 error relative to SqrtLasso") +
#   ggtitle("Average ell2 errors relative to SqrtLasso for LassoIdeal") +
#   geom_hline(yintercept = 1, linetype = "dashed", colour = "black") +
#   annotate("text", x = 0, y = 1.01, label = "SqrtLasso", hjust = 0, vjust = 0) +
#   # grid lines and ticks only at integer m values
#   scale_x_continuous(breaks = mvec, minor_breaks = NULL) +
#   # ylimits from 0 to 2
#     coord_cartesian(ylim = c(0, 2))

# #   theme(panel.grid.minor = element_blank(),
# #     panel.grid.major.x = element_line(colour = "grey", linetype = "dotted")) 
# p1



# df <- reshape::melt(ell2_avg[c(1, 3), ])
# # Plot the average ell2 errors as a function of m, one line per method
# # in one plot
# p1 <- ggplot(df, aes(x = m, y = value, colour = method)) +
#   geom_line() +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   labs(x = "m", y = "Average ell2 error") +
#   ggtitle("Average ell2 errors for different methods")
# p1


# # Plot the average estimates as a function of m
# library("ggplot2")
# library("gridExtra")
# thats_avg <- apply(thats, c(1, 2, 4), mean)
# df <- reshape::melt(thats_avg)

# # Plot the average estimates as a function of m for j=1
# p1 <- ggplot(df[df$j == 1, ], aes(x = m, y = value, colour = method)) +
#   geom_line() +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   labs(x = "m", y = "Average estimate") +
#   coord_cartesian(ylim = c(-0.1, 1)) +
#   # add horizontal lines at true value
#   geom_hline(yintercept = 1, linetype = "dashed", colour = "black") +
#   ggtitle("j = 1") +
#   # add explanation that the dashed lines are the true values
#   annotate("text", x = 0, y = 1.01, label = "True value", hjust = 0, vjust = 0)

# # Plot the average estimates as a function of m for j=2
# p2 <- ggplot(df[df$j == 2, ], aes(x = m, y = value, colour = method)) +
#   geom_line() +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   labs(x = "m", y = "Average estimate") +
#   coord_cartesian(ylim = c(-0.1, 1)) +
#   # add horizontal lines at true value
#   geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
#   ggtitle("j = 2") +
#   # add explanation that the dashed lines are the true values
#   annotate("text", x = 0, y = 0.01, label = "True value", hjust = 0, vjust = 0)

# # Combine the two plots
# gridp1p2_mean <- grid.arrange(p1, p2, ncol = 2)

# # Save the combined plot as pdf
# file_name <- paste0("mean_2_regresor_heterosked_sqrtLasso_comparison_n", n,
#                     "_p", p, ".pdf")

# ggsave(file_name, plot = gridp1p2_mean, width = 10, height = 5)

# 

## SANDBOX ##
# n <- 1e6
# p <- 2
# set.seed(123)
# RNGkind(normal.kind = "Kinderman-Ramage") # faster normal draws
# x <- matrix(rnorm(n * p), nrow = n, ncol = p) # standard normals
# mvec <- seq(0, 5, by = 1)
# numm <- length(mvec)
# vars_simu <- numeric(numm)
# vars_theo <- numeric(numm)
# for (thism in 1:numm) {
#   m <- mvec[thism]
#   vars_simu[thism] <- mean(var_eps_of_x(x, m = m))
#   vars_theo[thism] <- var_m(m = m)
# }
# cbind(vars_simu, vars_theo)

# cbind(mean(var_eps_x), var_m(m = m))