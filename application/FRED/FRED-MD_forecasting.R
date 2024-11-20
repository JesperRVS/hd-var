# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages for parallel computing
libpar <- c("doRNG", "doParallel", "foreach")
lapply(libpar, require, character.only = TRUE)

if (Sys.info()[["sysname"]] == "Linux") {
  setwd("../..") # if on Linux server, back up to parent folder
}

# Read data
data <- read.csv("application/FRED/data/FRED-MD_2022-05_preprocessed.csv")
colnames(data) <- NULL
data <- as.matrix(data) # convert from dataframe to matrix

# Forecast settings: Lag length, forecast horizon
n_all <- nrow(data)
p <- ncol(data)

qvec <- 1:4 #12                  # lags
numlag <- length(qvec)        # number of lags
qmax <- max(qvec)             # maximum lag

numfore <- 10 # 120              # number of forecast horizons

n <- n_all - numfore - qmax   # number of observations used for estimation

methods <- c("Lasso", "PostLasso", "SqrtLasso",
             "BICLasso", "PostSqrtLasso", "PostBICLasso")
nummet <- length(methods)

# Placeholders
ivwsfes <- array(NA, dim = c(numfore, numlag, nummet))
dimnames(ivwsfes) <- list(horizon = 1:numfore, q = qvec, method = methods)

# Sample variances and their inverses (using pre-processed sample)
varvec <- as.matrix(apply(data, 2, var))
invvarvec <- 1 / varvec

cl <- makeCluster(detectCores())
registerDoParallel(cl)

# Loop over lags
for (thislag in 1:numlag) {
  q <- qvec[thislag] # current lag
  print(paste("Lag", q, "started."))
  # Loop over periods
  results <- foreach(thisnplus1 = 1:numfore) %dopar% {
    nplus1 <- qmax + n + thisnplus1           # current one-period-ahead
    sample <- (nplus1 - (n + q)):(nplus1 - 1) # current observations indices
    y_est <- data[sample, ]                   # current data for estimation

    # ESTIMATE
    # 1. VAR(q) LASSO
    source("lassoVAR.R")
    fit_lasso <- lasso_var(data = y_est, q = q, post = FALSE, intercept = TRUE)
    # 2. VAR(q) Post-LASSO
    fit_post <- lasso_var(data = y_est, q = q, post = TRUE, intercept = TRUE)
    # 3./5. VAR(q) Sqrt-LASSO and Post-Sqrt-LASSO
    source("sqrtLassoVAR.R")
    fit_sqrt <- sqrt_lasso_var(data = y_est, q = q,
                               post = TRUE, intercept = TRUE)
    # 4./6 VAR(q) BIC LASSO and Post-BIC-LASSO
    source("icLassoVAR.R")
    fit_bic <- ic_lasso_var(data = y_est, q = q, criteria = "bic",
                            post = TRUE, intercept = TRUE)

    # FORECAST
    # Fetch outcome to be forecasted
    y_nplus1 <- as.matrix(data[nplus1, ])
    # Fetch predictors for the forecast
    zn_mat <- matrix(NA, q, p)
    for (lag in 1:q) {
      zn_mat[lag, ] <- data[(nplus1 - lag), ]
    }
    # Stack the predictors into a vector
    z_n <- as.vector(t(zn_mat)) # transpose, then stack over columns (lags)
    z_n <- as.matrix(z_n) # convert to matrix to allow matrix multiplication

    # 1. VAR(q) LASSO prediction, forecast errors, IVWSFE
    y_nplus1_lasso <- fit_lasso$intr + fit_lasso$that %*% z_n   # p x 1
    fe_lasso <- y_nplus1 - y_nplus1_lasso                       # p x 1
    ivwsfe_lasso <- sum(invvarvec * fe_lasso^2)                 # scalar
    # 2. VAR(q) Post-LASSO
    y_nplus1_post <- fit_post$intr + fit_post$that %*% z_n
    fe_post <- y_nplus1 - y_nplus1_post
    ivwsfe_post <- sum(invvarvec * fe_post^2)
    # 3./5. VAR(q) Sqrt-LASSO and Post-Sqrt-LASSO
    # Sqrt-LASSO
    y_nplus1_sqrt <- fit_sqrt$intr + fit_sqrt$that %*% z_n
    fe_sqrt <- y_nplus1 - y_nplus1_sqrt
    ivwsfe_sqrt <- sum(invvarvec * fe_sqrt^2)
    # Post-Sqrt-LASSO
    y_nplus1_sqrt_post <- fit_sqrt$intr_post + fit_sqrt$that_post %*% z_n
    fe_sqrt_post <- y_nplus1 - y_nplus1_sqrt_post
    ivwsfe_sqrt_post <- sum(invvarvec * fe_sqrt_post^2)
    # 4. VAR(q) BIC-LASSO
    y_nplus1_bic <- fit_bic$intrs[, 1] + fit_bic$thats[, , 1] %*% z_n
    fe_bic <- y_nplus1 - y_nplus1_bic
    ivwsfe_bic <- sum(invvarvec * fe_bic^2)
    # 6. VAR(q) Post-BIC-LASSO
    y_nplus1_bic_post <-
      fit_bic$intrs_post[, 1] + fit_bic$thats_post[, , 1] %*% z_n
    fe_bic_post <- y_nplus1 - y_nplus1_bic_post
    ivwsfe_bic_post <- sum(invvarvec * fe_bic_post^2)

    # Return IVWSFEs for current q and period (one for each method)
    c(ivwsfe_lasso, ivwsfe_post, ivwsfe_sqrt,
      ivwsfe_bic, ivwsfe_sqrt_post, ivwsfe_bic_post)
  }
  print(paste("Lag", q, "done."))
  # convert results into a matrix with numfore rows and nummet columns
  ivwsfes_q <- matrix(unlist(results), nrow = numfore, byrow = TRUE)
  ivwsfes[, thislag, ] <- ivwsfes_q
}
stopCluster(cl)

# Save the workspace
if (Sys.info()[["sysname"]] == "Linux") {
  file_name <- paste0("application_workspace", "_N_", numfore, "_qmax_", qmax)
  save.image(file = paste0("application/FRED/", file_name, ".RData"))
  q("no")
}