# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages for parallel computing
library("doRNG")
library("doParallel")
library("foreach")

# SIMULATE DATA
n <- 100
p <- 2
s <- 1
sigma_x <- 1 # standard deviation of irrelevant regressors
sigma_x_irrel <- 10

sim_data_hetero <- function(n = 100, p = 4, s = 1, m = 0) {
  z1 <- matrix(rnorm(n * p), nrow = n, ncol = p) # standard normals
  x <- cbind(z1[, 1], sigma_x_irrel * z1[, 2:p])
  # x <- sigma_x * matrix(rnorm(n * p), nrow = n, ncol = p) # X ~ N(0, sigma_x^2)
  z <- matrix(rnorm(n), nrow = n, ncol = 1)     # standard normals
  eps <- exp(-0.5 * m * x[, 2]^2) * z   # eps | X = x ~ N(0, exp(-m * x_j^2))
  # eps <- sqrt(abs(x[, 1])^m) * z  # eps | X = x ~ N(0, |x_j|^m)
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
lambda_lasso_glmnet <- lambda_lasso / (2 * n)

mvec <- seq(0, 10, by = 1)
numm <- length(mvec)


# Monte Carlo (MC) settings
RNGkind(normal.kind = "Kinderman-Ramage") # faster normal draws
iseed <- 2345
# registerDoRNG(seed = iseed)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = iseed)
nummc <- 1000

thats <- array(NA, dim = c(p, 3, nummc, numm))
dimnames(thats) <- list(j = 1:p,
                        method = c("SqrtLasso", "LassoInit", "LassoIdeal"),
                        mc = 1:nummc,
                        m = mvec)
# scales <- array(NA, dim = c(p, numm))

for (thism in 1:numm) {
  m <- mvec[thism]
  results <- foreach(icount(nummc)) %dopar% {
    # Simulate data
    data <- sim_data_hetero(n = n, p = p, s = s, m = m)
    # Estimate coefficients
    source("helper_functions.R", local = TRUE)
    thats_mc <- matrix(NA, nrow = p, ncol = 3)
    # SqrtLasso
    that_sqrtl <- sqrt_lasso(x = data$x, y =  data$y, lambda = lambda_sqrtl)
    thats_mc[, 1] <- that_sqrtl
    # vars_x <- colMeans(data$x^2)
    # var_eps <- mean(data$eps^2)
    # exps_of_squares <- (1 / n) * crossprod(data$eps^2, data$x^2)
    # scales_mc <- exps_of_squares / (var_eps * vars_x)
    # print(scales_mc)
    # LassoInit
    ups_init <- sqrt((1 / n) * crossprod(data$y^2, data$x^2))
    that_init <- mult_lasso(x = data$x, y = data$y,
                            lambda_glmnet = lambda_lasso_glmnet,
                            upsilon = ups_init)
    thats_mc[, 2] <- that_init
    # LassoIdeal
    ups_ideal <- sqrt((1 / n) * crossprod(data$eps^2, data$x^2))
    that_ideal <- mult_lasso(x = data$x, y = data$y,
                             lambda_glmnet = lambda_lasso_glmnet,
                             upsilon = ups_ideal)
    thats_mc[, 3] <- that_ideal
    return(thats_mc)
  } # end mc loop
  thats[, , , thism] <- array(unlist(results), dim = c(p, 3, nummc))
} # end m loop
# Stop parallel computing
stopCluster(cl)

# Plot the average estimates as a function of m
library("ggplot2")
library("gridExtra")
thats_avg <- apply(thats, c(1, 2, 4), mean)
df <- reshape::melt(thats_avg)

# Plot the average estimates as a function of m for j=1
p1 <- ggplot(df[df$j == 1, ], aes(x = m, y = value, colour = method)) +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "m", y = "Average estimate") +
  coord_cartesian(ylim = c(-0.1, 1)) +
  # add horizontal lines at true value
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black") +
  ggtitle("j = 1") +
  # add explanation that the dashed lines are the true values
  annotate("text", x = 0, y = 1.01, label = "True value", hjust = 0, vjust = 0)

# Plot the average estimates as a function of m for j=2
p2 <- ggplot(df[df$j == 2, ], aes(x = m, y = value, colour = method)) +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "m", y = "Average estimate") +
  coord_cartesian(ylim = c(-0.1, 1)) +
  # add horizontal lines at true value
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  ggtitle("j = 2") +
  # add explanation that the dashed lines are the true values
  annotate("text", x = 0, y = 0.01, label = "True value", hjust = 0, vjust = 0)

# Combine the two plots
gridp1p2_mean <- grid.arrange(p1, p2, ncol = 2)

# Save the combined plot as pdf
file_name <- paste0("mean_exponential_heterosked_sqrtLassoBreakDown_n", n,
                    "_p", p, "_sigma_x", sigma_x, ".pdf")

ggsave(file_name, plot = gridp1p2_mean, width = 10, height = 5)

# Create similar plots but based on the median instead of the mean
thats_med <- apply(thats, c(1, 2, 4), median)
df_med <- reshape::melt(thats_med)

# Plot the median estimates as a function of m
p1_med <- ggplot(df_med[df_med$j == 1, ], aes(x = m, y = value, colour = method)) +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "m", y = "Median estimate") +
  coord_cartesian(ylim = c(-0.1, 1)) +
  # add horizontal lines at true value
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black") +
  ggtitle("j = 1") +
  # add explanation that the dashed lines are the true values
  annotate("text", x = 0, y = 1.01, label = "True value", hjust = 0, vjust = 0)

# Plot the median estimates as a function of m for j=2
p2_med <- ggplot(df_med[df_med$j == 2, ], aes(x = m, y = value, colour = method)) +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "m", y = "Median estimate") +
  coord_cartesian(ylim = c(-0.1, 1)) +
  # add horizontal lines at true value
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  ggtitle("j = 2") +
  # add explanation that the dashed lines are the true values
  annotate("text", x = 0, y = 0.01, label = "True value", hjust = 0, vjust = 0)

# Combine the two plots
gridp1p2_med <- grid.arrange(p1_med, p2_med, ncol = 2)

# Save the combined plot as high quality pdf
file_name_med <- paste0("median_exponential_heterosked_sqrtLassoBreakDown_n", n,
                        "_p", p, "_sigma_x", sigma_x, ".pdf")
ggsave(file_name_med, plot = gridp1p2_med, width = 10, height = 5)






# # Plot the median estimates as a function of m for j=1
# p1_med <- ggplot(df_med[df_med$j == 1, ], aes(x = m, y = value, colour = method)) +
#   geom_line() +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   labs(x = "m", y = "Median estimate") +
#   coord_cartesian(ylim = c(-0.1, 1)) +
#   # add horizontal lines at true value
#   geom_hline(yintercept = 1, linetype = "dashed", colour = "black") +
#   ggtitle("j = 1") +
#   # add explanation that the dashed lines are the true values
#   annotate("text", x = 0, y = 1.01, label = "True value", hjust = 0, vjust = 0)

# # Plot the median estimates as a function of m for j=2
# p2_med <- ggplot(df_med[df_med$j == 2, ], aes(x = m, y = value, colour = method)) +
#   geom_line() +
#   theme_minimal() +
#   theme(legend.position = "bottom") +
#   labs(x = "m", y = "Median estimate") +
#   coord_cartesian(ylim = c(-0.1, 1)) +
#   # add horizontal lines at true value
#   geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
#   ggtitle("j = 2") +
#   # add explanation that the dashed lines are the true values
#   annotate("text", x = 0, y = 0.01, label = "True value", hjust = 0, vjust = 0)

# # Combine the two plots
# gridp1p2_med <- grid.arrange(p1_med, p2_med, ncol = 2)

# # Save the combined plot as high quality pdf
# ggsave("sqrtLassoBreakDown_med.pdf", plot = gridp1p2_med, width = 10, height = 5)


# thats[1,,]
# colMeans(thats[1, , ])

## SANDBOX ##
# data <- sim_data_hetero(n = n, p = p, s = s, m = 4)
# x <- data$x
# y <- data$y
# ups_init <- sqrt((1 / n) * crossprod(y^2, x^2))
# xtil_init <- sweep(x, 2, ups_init, FUN = "/")

# fit_init_temp <- glmnet(x = xtil_init, y = y, family = "gaussian",
#                         lambda = lambda_lasso / (2 * n),
#                         standardize = FALSE, intercept = FALSE,
#                         thresh = 1e-4)
# that_init <- as.matrix(coef(fit_init_temp)[-1] / ups_init)

# that_lasso <- mult_lasso(x, y,
#                          lambda_glmnet = lambda_lasso / (2 * n),
#                          upsilon = ups_init)


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

