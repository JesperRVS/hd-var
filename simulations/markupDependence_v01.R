# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# Packages for parallel computing
libpar <- c("doRNG", "doParallel", "foreach")
lapply(libpar, require, character.only = TRUE)

if (Sys.info()[["sysname"]] == "Linux") {
  setwd("..") # if on Linux server, back up to parent folder
}

testrun <- FALSE # whether to run a test simulation

if (testrun) {
  nvec <- seq(from = 100, to = 500, by = 100)
  pvec <- c(4, 8, 16)
  designs <- c("Diagonal")#, "NearBand")
  methods <- c("Lasso", "PostLasso", "SqrtLasso", "PostSqrtLasso")
  nburn <- 100
  nummc <- 80 # no. MC repetitions
} else {
  nvec <- seq(from = 200, to = 1000, by = 200)
  pvec <- c(16, 32, 64, 128)
  designs <- c("Diagonal")#, "NearBand", "BlockDiag",
               #"Correlated", "HeavyTailed", "Heteroskedastic", "NearUnity")
  methods <- c("Lasso", "PostLasso", "SqrtLasso", "PostSqrtLasso")
  nburn <- 1000
  nummc <- 1000 # no. MC repetitions TODO
}
cvec <- c(1, 1.05, 1.1, 2, 10)
numn <- length(nvec)
nump <- length(pvec)
numdes <- length(designs)
nummet <- length(methods)
numc <- length(cvec)

## == SOME HELPER FUNCTIONS == ##

q_switch <- function(design) {
  q <- list(Diagonal = 1, NearBand = 1, BlockDiag = 4,
            Correlated = 1, HeavyTailed = 1, Heteroskedastic = 1,
            NearUnity = 1)
  return(q[[design]])
}

# Function to determine the (correct) coefficient matrix Theta
coef_mat_diag <- function(p) {0}
coef_mat_toep <- function(p) {0}
coef_mat_block <- function(p) {0}
insertSource("simulations/simData.R",
             functions = c("coef_mat_diag", "coef_mat_toep", "coef_mat_block"))
theta_switch <- function(design, p, n) {
  if (design %in% c("Diagonal", "Correlated",
                    "HeavyTailed", "Heteroskedastic")) {
    theta <- coef_mat_diag(p = p, coef = 0.5)
  } else if (design == "NearUnity") {
    theta <- coef_mat_diag(p = p, coef = 1 - (5 / n))
  } else if (design == "NearBand") {
    theta <- coef_mat_toep(p = p)
  } else if (design == "BlockDiag") {
    theta <- coef_mat_block(p = p)
  } else {
    stop("Design not recognized.")
  }
  return(theta)
}

# Function which calculates the maximum ell_2 norm of the rows of a matrix
max_ell2_row <- function(mat) {
  max_ell2 <- max(apply(mat, 1, function(x) sqrt(sum(x^2))))
  return(max_ell2)
}

# Function which calculates the maximum row sparsity of a matrix
max_row_sparsity <- function(mat) {
  max_sparsity <- max(apply(mat, 1, function(x) sum(x != 0)))
  return(max_sparsity)
}

## == SIMULATIONS == ##

# Monte Carlo (MC) settings
RNGkind(normal.kind = "Kinderman-Ramage") # faster normal draws
iseed <- 2345 # set seed for reproducibility
cl <- makeCluster(detectCores())
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = iseed)

# Placeholders
max_ell2_errors <- array(NA, dim = c(nummc, numn, nump, numdes, nummet, numc))
dimnames(max_ell2_errors) <- list(mc = 1:nummc, n = nvec, p = pvec,
                                  design = designs, method = methods, c = cvec)
max_row_sparsities <- array(NA, dim = c(nummc, numn, nump,
                                        numdes, nummet, numc))
dimnames(max_row_sparsities) <- list(mc = 1:nummc, n = nvec, p = pvec,
                                     design = designs, method = methods,
                                     c = cvec)
num_upd <- array(NA, dim = c(nummc, numn, nump, numdes, 2, numc))
dimnames(num_upd) <- list(mc = 1:nummc, n = nvec, p = pvec,
                          design = designs, method = c("Lasso", "PostLasso"),
                          c = cvec)

intercept <- FALSE # whether to include intercept in simulations

for (thisc in seq_along(cvec)) {
  c <- cvec[thisc] # set score markup
  for (this_design in seq_along(designs)) {
    design <- designs[this_design]
    q <- q_switch(design) # set lag length q (correctly)
    for (thisn in seq_along(nvec)) {
      n <- nvec[thisn]
      for (thisp in seq_along(pvec)) {
        p <- pvec[thisp]
        theta <- theta_switch(design, p = p, n = n)
        print(paste("Design:", design, ", n:", n, ", p:", p, ", c:", c))
        results <- foreach(icount(nummc)) %dopar% {
          # GENERATE DATA
          source("simulations/simData.R", local = TRUE)
          data <- sim_data_by_design(n = n, p = p, design = design,
                                     nburn = nburn)
          # ESTIMATE VAR MODEL
          errors <- numeric(nummet)
          shats <- numeric(nummet)
          kterms <- numeric(2)

          # 1-2. LASSO AND POST-LASSO
          # 1. LASSO
          source("lassoVAR.R", local = TRUE) # for lasso_var
          fit_lasso <- lasso_var(data = data, q = q, post = FALSE,
                                 intercept = intercept, c = c) # <- note c
          that_lasso <- fit_lasso$that # extract estimates
          errors[1] <- max_ell2_row(that_lasso - theta) # calculate errors
          shats[1] <- max_row_sparsity(that_lasso) # calculate sparsity
          kterms[1] <- fit_lasso$k_term # number of loading updates
          # 2. POST-LASSO
          fit_postl <- lasso_var(data = data, q = q, post = TRUE,
                                 intercept = intercept, c = c) # <- note c
          that_postl <- fit_postl$that # extract estimates
          errors[2] <- max_ell2_row(that_postl - theta) # calculate errors
          shats[2] <- max_row_sparsity(that_postl) # calculate sparsity
          kterms[2] <- fit_postl$k_term # number of loading updates
          # Note: Have to do separate runs, as Post-Lasso uses refitting in
          # every step of the algorithm.

          # 3-4. SQRT-LASSO AND POST-SQRT-LASSO
          source("sqrtLassoVAR.R", local = TRUE) # for sqrt_lasso_var
          fit_sqrtl <- sqrt_lasso_var(data = data, q = q, post = TRUE,
                                      intercept = intercept, c = c) # <- note c
          # Note: post = TRUE to get also Post-Sqrt-Lasso
          # 3. SQRT-LASSO
          that_sqrtl <- fit_sqrtl$that
          errors[3] <- max_ell2_row(that_sqrtl - theta)
          shats[3] <- max_row_sparsity(that_sqrtl)
          # 4. POST-SQRT-LASSO
          that_postsqrtl <- fit_sqrtl$that_post
          errors[4] <- max_ell2_row(that_postsqrtl - theta)
          shats[4] <- max_row_sparsity(that_postsqrtl)

          list(errors, shats, kterms) # return results as list
        } # mc loop
        # Extract max rowwise ell_2 errors, max row sparsities
        # and number of loading updates (Lasso and PostLasso), respectively
        max_ell2_errors[, thisn, thisp, this_design, , thisc] <-
          t(sapply(lapply(results, "[[", 1), unlist))
        max_row_sparsities[, thisn, thisp, this_design, , thisc] <-
          t(sapply(lapply(results, "[[", 2), unlist))
        num_upd[, thisn, thisp, this_design, , thisc] <-
          t(sapply(lapply(results, "[[", 3), unlist))
      } # p loop
    } # n loop
  } # design loop
} # c loop
stopCluster(cl)

# Save the workspace
if (Sys.info()[["sysname"]] == "Linux") {
  file_name <- paste("markup_dependence_workspace_", nummc, "_MC_",
                     min(nvec), "_to_", max(nvec),
                     "_n_", min(pvec), "_to_", max(pvec), "_p",
                     "_diagonal_only", sep = "")
  save.image(file = paste0("simulations/", file_name, ".RData"))
  q("no")
}