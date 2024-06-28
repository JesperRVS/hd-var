lassoVAR <- function(Ydata, q = 1, post = 1, intr = 1, c = 1.1, gamma = NULL, K = 15, 
                     tol_Ups = 1e-3, tol_glmnet = 1e-4, nowarn = 0) {
  
  if (missing(Ydata) || is.null(Ydata)) {
    stop("Not enough input arguments.")
  }
  
  # Unpack data
  nplusq <- nrow(Ydata)
  p <- ncol(Ydata)
  n <- nplusq - q
  qp <- q * p
  
  Y <- Ydata[(q + 1):(q + n), ]
  X <- matrix(NA, nrow = n, ncol = qp)
  
  for (ell in 1:q) {
    block_ell <- (ell - 1) * p + 1:ell * p
    lag_ell <- (q + 1 - ell):(q + n - ell)
    X[, block_ell] <- Ydata[lag_ell, ]
  }
  
  # Penalty level
  if (is.null(gamma)) {
    gamma <- 0.1 / log(max(n, qp))
  }
  mstar <- floor(n^(1/5))
  lstar <- floor(n / mstar)
  lambdastar <- ((2 * c * n) / sqrt(mstar * lstar)) * qnorm(1 - gamma / (2 * q * p^2))
  lambda_glmnet <- lambdastar / (2 * n)
  
  # Set tolerances
  if (missing(K) || is.null(K)) {
    K <- 15
  }
  if (missing(tol_Ups) || is.null(tol_Ups)) {
    tol_Ups <- 1e-3
  }
  if (missing(tol_glmnet) || is.null(tol_glmnet)) {
    tol_glmnet <- 1e-4
  }
  if (missing(nowarn) || is.null(nowarn)) {
    nowarn <- 0
  }
  
  if (intr) {
    Ybar <- colMeans(Y)
    Y <- sweep(Y, 2, Ybar)
    Xbar <- colMeans(X)
    X <- sweep(X, 2, Xbar)
  }
  
  UpsInit <- sqrt(colMeans(Y^2) %*% t(colMeans(X^2)))
  ThatInit <- Design_A$mlasso(X, Y, lambda_glmnet, UpsInit, tol_glmnet)
  
  if (post) {
    sel <- ThatInit != 0
    for (i in 1:p) {
      ThatInit[i, sel[i, ]] <- solve(t(X[, sel[i, ]]) %*% X[, sel[i, ]]) %*% t(X[, sel[i, ]]) %*% Y[, i]
    }
  }
  
  if (intr) {
    constInit <- Ybar - ThatInit %*% Xbar
  } else {
    constInit <- NULL
  }
  
  UpsRefi <- array(NA, dim = c(p, qp, K))
  ThatRefi <- array(NA, dim = c(p, qp, K))
  if (intr) {
    constRefi <- matrix(NA, nrow = p, ncol = K)
  } else {
    constRefi <- NULL
  }
  
  for (k in 1:K) {
    if (k == 1) {
      Res_kminus1 <- Y - X %*% t(ThatInit)
    } else {
      Res_kminus1 <- Y - X %*% t(ThatRefi[, , k - 1])
    }
    
    UpsRefi[, , k] <- sqrt(colMeans(Res_kminus1^2) %*% t(colMeans(X^2)))
    ThatRefi_k <- Design_A$mlasso(X, Y, lambda_glmnet, UpsRefi[, , k], tol_glmnet)
    
    if (post) {
      sel <- ThatRefi_k != 0
      for (i in 1:p) {
        ThatRefi_k[i, sel[i, ]] <- solve(t(X[, sel[i, ]]) %*% X[, sel[i, ]]) %*% t(X[, sel[i, ]]) %*% Y[, i]
      }
    }
    
    ThatRefi[, , k] <- ThatRefi_k
    
    if (intr) {
      constRefi[, k] <- Ybar - ThatRefi_k %*% Xbar
    }
    
    if (k == 1) {
      UpsOld <- UpsInit
    } else {
      UpsOld <- UpsRefi[, , k - 1]
    }
    
    UpsNew <- UpsRefi[, , k]
    dUps <- UpsNew - UpsOld
    reldiffUps <- sqrt(sum(dUps^2)) / sqrt(sum(UpsOld^2))
    
    if (reldiffUps <= tol_Ups) {
      break
    }
  }
  
  Ups <- UpsRefi[, , k]
  That <- ThatRefi[, , k]
  if (intr) {
    const <- Ybar - That %*% Xbar
  } else {
    const <- NULL
  }
  
  if (nowarn == 0 && k == K && reldiffUps > tol_Ups) {
    warning("Maximum number of updates reached. Consider increasing K.")
    cat(sprintf("Relative change in penalty loadings is %3.1g percent\n", 100 * reldiffUps))
  }
  
  fit <- list(
    lambda = lambdastar,
    UpsInit = UpsInit,
    UpsRefi = UpsRefi,
    Ups = Ups,
    constInit = constInit,
    constRefi = constRefi,
    const = const,
    ThatInit = ThatInit,
    ThatRefi = ThatRefi,
    That = That,
    kterm = k,
    reldiffUps = reldiffUps
  )
  
  return(fit)
}
