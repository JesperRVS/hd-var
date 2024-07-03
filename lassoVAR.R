# ESTIMATION TOOLS FOR VECTOR AUTOREGRESSIONS WITH LASSO PENALIZATION

# Dependencies
library("glmnet", "MASS")

# Multiple-equation weighted LASSO
mult_lasso <- function(x, y, lambda_glmnet, upsilon = NULL, tol_glmnet = 1e-4) {
  if (missing(x) || missing(y) || missing(lambda_glmnet)) {
    stop("Not enough input arguments.")
  }
  pq <- ncol(x)  # extract dimensions
  p <- ncol(y)   # q = autoregressive order, p = dim(output)
  if (nrow(x) != nrow(y)) {
    stop("Number of observations in x and y do not match.")
  }
  if (is.null(upsilon)) {
    upsilon <- matrix(1, p, pq)  # default to unit loadings
  }
  # Options for glmnet
  opts <- list(standardize = FALSE,    # don't standardize (rescaling below)
               intercept = FALSE,      # don't fit an intercept
               lambda = lambda_glmnet, # penalty in eyes of glmnet
               thresh = tol_glmnet)    # tolerance for coordinate descent
  # Estimate
  that <- matrix(NA, p, pq)
  for (i in 1:p) {
    xtilde <- sweep(x, 2, upsilon[i, ], FUN = "/")  # rescale using loadings
    fit_i <- glmnet(xtilde, y[, i], family = "gaussian",
                    standardize = opts$standardize, intercept = opts$intercept,
                    lambda = opts$lambda, thresh = opts$thresh)
    that[i, ] <- as.numeric(coef(fit_i))[-1] / upsilon[i, ]  # original scaling
  }
  return(that)
}

lasso <- function(data, q = 1, post = TRUE, intr = TRUE,
                      c = 1.1, gamma = 0.1 / log(max(dim(data))), k = 15,
                      tol_ups = 1e-3, tol_glmnet = 1e-4, warn = TRUE) {
  if (missing(data) || is.null(data)) {
    stop("No data provided.")
  }
  # Unpack data
  nplusq <- nrow(data)  # n + q
  p <- ncol(data)       # p = dim of output
  n <- nplusq - q       # n = effective sample size
  pq <- p*q             # pq = total number of variables
  y <- data[(q + 1):(q + n), ] # response
  x <- matrix(NA, nrow = n, ncol = pq) # predictors
  for (ell in 1:q) {
    # where to store ellth lag --v
    block_ell <- ((ell - 1) * p + 1):(ell * p)
    # periods corresponding to ellth lag --v
    lag_ell <- (q + 1 - ell):(q + n - ell)
    # store 1st lags first, 2nd second,... --v
    x[, block_ell] <- data[lag_ell, ]
  }
  # Penalty level
  lambda_star <- 2 * c * sqrt(n) * qnorm(1 - gamma / (2 * q * p^2))
  lambda_glmnet <- lambda_star / (2 * n) # penalty in eyes of glmnet
  # Note: glmnet defines "lambda" based on *half* of *average* square loss,
  # while our "lambda" stems from *sum* square loss (w/o the 1/2)
  #  If intercepts requested, demean before proceeding
  if (intr == TRUE) {
    ybar <- colMeans(y)     # means as p-dim array
    xbar <- colMeans(x)     # means as pq-dim array
    y <- sweep(y, 2, ybar)  # demean response
    x <- sweep(x, 2, xbar)  # demean predictors
    ybar <- as.matrix(ybar) # means as p x 1 matrix
    xbar <- as.matrix(xbar) # means as pq x 1 matrix
  }
  # Note: Means are stored to back out intercepts later
  # ESTIMATION
  # Initial step
  ups_init <- sqrt((1 / n) * crossprod(y^2, x^2)) # initial penalty loadings
  that_init <- mult_lasso(x, y, lambda_glmnet, ups_init, tol_glmnet) # and estimates
  if (post == TRUE) {                       # if refitting requested...
    refit <- mult_refit(x, y, that_init)    # refit initial estimates
    that_init <- refit$that                 # overwrite (keeping zeros)
    full_rank_post_init <- refit$full_rank  # check full rank
  } else { # if no refitting requested...
    full_rank_post_init <- NULL # then rank deficiency is irrelevant
  }
  if (intr == TRUE) { # if intercepts requested...
    intr_init <- ybar - that_init %*% xbar # back them out
  } else{
    intr_init <- NULL # o/w don't
  }
  # # Updating
  # ups_refi <- array(NA, dim = c(p, pq, k))  # refined penalty loadings
  # that_refi <- array(NA, dim = c(p, pq, k)) # and estimates
  # if (intr == TRUE) { # if intercepts requested...
  #   intr_refi <- matrix(NA, p, k) # create placeholder
  # } else {
  #   intr_refi <- NULL # o/w don't
  # }
  # for (l in 1:k) {
  #   if (l == 1) {
  #     that_old <- that_init
  #   } else {
  #     that_old <- that_refi[, , l - 1]
  #   }
  #   res_old <- y - x %*% t(that_old)
  #   ups_new <- sqrt((1 / n) * crossprod(res_old^2, x^2))
  #   ups_refi[, , l] <- ups_new
  #   that_new <- mult_lasso(x, y, lambda_glmnet, ups_new, tol_glmnet)
  #   if (post == TRUE) {
  #     sel <- that_new != 0
  #     intr_post <- matrix(NA, p, 1)     # p x 1 intercepts
  #     that_post <- matrix(NA, p, pq)    # p x pq slopes
  #     full_rank_post_init <- matrix(NA, p, 1) # p x 1 full rank flag
  #     for (i in 1:p) {
  #       shat_
        
  #     }
  #   }
  #   that_refi[, , k] <- that_new
    
  #   if (intr) {
  #     intr_refi[, k] <- Ybar - that_new %*% Xbar
  #   }
    
  #   if (k == 1) {
  #     UpsOld <- UpsInit
  #   } else {
  #     UpsOld <- ups_refi[, , k - 1]
  #   }
    
  #   UpsNew <- ups_refi[, , k]
  #   dUps <- UpsNew - UpsOld
  #   reldiffUps <- sqrt(sum(dUps^2)) / sqrt(sum(UpsOld^2))
    
  #   if (reldiffUps <= tol_Ups) {
  #     break
  #   }
  # }
  
  # Ups <- ups_refi[, , k]
  # That <- that_refi[, , k]
  # if (intr) {
  #   const <- Ybar - That %*% Xbar
  # } else {
  #   const <- NULL
  # }
  
  # if (nowarn == 0 && k == K && reldiffUps > tol_Ups) {
  #   warning("Maximum number of updates reached. Consider increasing K.")
  #   cat(sprintf("Relative change in penalty loadings is %3.1g percent\n", 100 * reldiffUps))
  # }
  fit <- list(
                lambda = lambda_star,
                ups_init = ups_init,
                intr_init = intr_init,
                that_init = that_init,
                full_rank_post_init = full_rank_post_init
                )
  return(fit)
}

# Helper functions

# Refit for multiple equations with same regressors
mult_refit <- function(x, y, that) {
  n <- nrow(y)                          # num of observations
  p <- ncol(y)                          # num of equations
  sel <- that != 0                      # flag selections
  full_rank <- matrix(NA, p, 1)         # full rank flag
  for (i in 1:p) {
    shat_i <- sum(sel[i, ])             # num selected
    if (shat_i > 0) {                   # if something selected...
      xi <- as.matrix(x[, sel[i, ]])    # get active regressors
      refit_i <- ls_sol(xi, y[, i])     # refit
      that[i, sel[i, ]] <- refit_i$sol  # overwrite (keeping zeros)
      full_rank[i, 1] <- refit_i$full_rank # check rank
    } else {                            # if nothing selected...
      full_rank[i, 1] <- TRUE           # rank consider "full"
    }
  }
  refit <- list(that = that, full_rank = full_rank)
  return(refit)
}

# Find solution to minimize ||y - Xb||_2 over b
ls_sol <- function(x, y) {
  full_rank <- dim(x)[2] == rankMatrix(x)[1] # check for full rank
  if (full_rank == TRUE) { # if full rank, produce unique solution  
    sol <- solve(crossprod(x), crossprod(x, y))
  } else { # o/w produce solution based on pseudoinverse
    sol <- ginv(x) %*% y # ginv is Moore-Penrose pseudoinverse
  }
  fit <- list(sol = sol, full_rank = full_rank)
  return(fit)
}