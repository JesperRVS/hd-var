# ESTIMATION TOOLS FOR VECTOR AUTOREGRESSIONS WITH LASSO PENALIZATION

# Dependencies
library("glmnet")

# Multiple-equation weighted LASSO
mult_lasso <- function(x, y, lambda_glmnet, upsilon = NULL, tol_glmnet = 1e-4) {
  if (missing(x) || missing(y) || missing(lambda_glmnet)) {
    stop("Not enough input arguments.")
  }
  qp <- ncol(x)  # extract dimensions
  p <- ncol(y)   # q = autoregressive order, p = dim(output)
  if (nrow(x) != nrow(y)) {
    stop("Number of observations in x and y do not match.")
  }
  if (is.null(upsilon)) {
    upsilon <- matrix(1, p, qp)  # default to unit loadings
  }
  # Options for glmnet
  opts <- list(standardize = FALSE,    # don't standardize (rescaling below)
               intercept = FALSE,      # don't fit an intercept
               lambda = lambda_glmnet, # penalty in eyes of glmnet
               thresh = tol_glmnet)    # tolerance for coordinate descent
  # Estimate
  that <- matrix(NA, p, qp)
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
  nplusq <- nrow(data)  # n plus q
  p <- ncol(data)       # p = dim of output
  n <- nplusq - q       # n = effective sample size
  qp <- q * p           # qp = total number of variables
  y <- data[(q + 1):(q + n), ] # response
  x <- matrix(NA, nrow = n, ncol = qp) # predictors
  for (ell in 1:q) {
    # where to store ellth lag --v
    block_ell <- ((ell - 1) * p + 1):(ell * p)
    # periods corresponding to ellth lag --v
    lag_ell <- (q + 1 - ell):(q + n - ell)
    # store 1st lags first, 2nd second,...
    x[, block_ell] <- data[lag_ell, ]
  }
  # Penalty level
  lambda_star <- 2 * c * sqrt(n) * qnorm(1 - gamma / (2 * q * p^2))
  lambda_glmnet <- lambda_star / (2 * n) # penalty in eyes of glmnet
  # Note: glmnet defines "lambda" based on *half* of *average* square loss,
  # while our "lambda" stems from *sum* square loss (w/o the 1/2)
  #  If intercepts requested, demean before proceeding
  if (intr == TRUE) {
    ybar <- colMeans(y)
    y <- sweep(y, 2, ybar)
    xbar <- colMeans(x)
    x <- sweep(x, 2, xbar)
  }
  # Note: Means stored to back out intercept estimates later
  # ESTIMATION
  # Initial step
  ups_init <- sqrt((1 / n) * t(y^2) %*% (x^2)) # loadings
  that_init <- mult_lasso(x, y, lambda_glmnet, ups_init, tol_glmnet) # estimates
  # if (post == TRUE) { # Post-selection refitting
  #   sel <- that_init != 0 # flag active sets of predictors
  #   intr_post <- NA # initialize
  #   that_post <- matrix(NA, qp, 1) # initialize
    # for (i in 1:p) {
    #   shat_i <-  sum(sel[i, ])
    #   if (shat_i > 0) { # if something selected, refit
    #     if (shat == rankMatrix(x[, sel[i, ]])[1]) { # if full rank...
    #       xi <- x[, sel[i, ]] # selected regressors
    #       sol <- solve(crossprod(xi), crossprod(xi, y[, i]))
    #     } else { # if rank deficient...

    #     }
    #     if (intr == TRUE) {
    #       xiaug <- cbind(1, x)
    #       if (1 + shat_i == rankMatrix(xiaug)[1]) {
    #         sol <- solve(crossprod(xiaug), crossprod(xiaug, y[, i]))
    #       }
    #     }
        # } else {

        # }
        # # if intercept requested, refit with intercept

        # if (intr == TRUE) {
        #   refit <- lm(y[, i] ~ 1 + x[, sel[i, ]])
        # } else { # otherwise, refit without intercept
        #   refit <- lm(y[, i] ~ 0 + x[, sel[i, ]])
        # }
        # if (length(fit$coefficients) == fit$rank) { # if full rank
          
        # } else { # otherwise proceed
        #   if (warn == TRUE) { # warn if refit fails
        #     warning(paste("Refitting for equation " i " failed."))
        #   }
        # }
      # }
        #   # Overwrite upon convergence; o/w keep as NA
        #   if (converged == TRUE) {
        #     intr_post[i] <- coef(refit)[1] # store intercept
        #     that_post[sel[i, ], i] <- coef(refit)[-1] # store refit
        #     that_post[!sel[i, ], i] <- 0 # keep zeros at zero
        #   } else {
        #     if (warn == TRUE) {
        #       warning(paste("Refit for equation", i, "did not converge."))
        #     }
        #   }
        # } else { # otherwise, refit without intercept
        #   refit <- glm(y[, i] ~ 0 + x[, sel[i, ]])
        # }
        # that_post[sel[i, ], i] <- coef(refit) # store refit
        # that_post[!sel[i, ], i] <- 0 # keep zeros at zero
        # # Calculate LS estimate using Cholesky factorization
        # xy_sel_i <- crossprod(x[, sel[i, ]], y[, i]) # RHS
        # ch_sel_i <- chol(crossprod(x[, sel[i, ]])) # Cholesky factorization
        # that_init[i, sel[i, ]] <- backsolve(ch_sel_i,
        #         forwardsolve(ch_sel_i, xy_sel_i, upper = TRUE, trans = TRUE))
        # # Solve R'Rb = X'y in two steps
        # xi <- x[, sel[i, ]] # selected predictors
        # that_init[i, sel[i, ]] <- solve(crossprod(xi), crossprod(xi, y[, i]))
      # }
  # }

  # # Intercept estimates
  # if (intr == TRUE) {
  #   constInit=Ybar'-ThatInit*Xbar'; % ...back them out
  # } else {
  #   constInit=[];
  # }
  fit <- list(
                lambda = lambda_star,
                ups_init = ups_init,
                that_init = that_init)
  return(fit)
}