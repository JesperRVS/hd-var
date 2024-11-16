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

# Simulation settings
if (testrun) {
  nvec <- seq(from = 100, to = 500, by = 100)
  pvec <- c(16, 32, 64, 128)
  designs <- c("Diagonal", "Correlated", "HeavyTailed",
               "BlockDiag", "NearBand",
               "Heteroskedastic_y", "Heteroskedastic_eta")
  methods <- c("Lasso", "PostLasso",
               "BICLasso", "PostBICLasso",
               "SqrtLasso", "PostSqrtLasso")
  nburn <- 100
  nummc <- 240 # no. MC repetitions
} else {
  nvec <- seq(from = 200, to = 1000, by = 200)
  # nvec <- seq(from = 100, to = 1000, by = 100) # TODO
  pvec <- c(16, 32, 64, 128)
  designs <- c("Diagonal", "Correlated", "HeavyTailed",
               "BlockDiag", "NearBand",
               "Heteroskedastic_y", "Heteroskedastic_eta")
  methods <- c("Lasso", "PostLasso",
               "BICLasso", "PostBICLasso",
               "SqrtLasso", "PostSqrtLasso")
  nburn <- 1000
  nummc <- 1000 # no. MC repetitions TODO
}
numn <- length(nvec)
nump <- length(pvec)
numdes <- length(designs)
nummet <- length(methods)

## == SOME HELPER FUNCTIONS == ##

q_switch <- function(design) {
  q <- list(Diagonal = 1, Correlated = 1, HeavyTailed = 1,
            BlockDiag = 4, NearBand = 1,
            Heteroskedastic_y = 1, Heteroskedastic_eta = 1)
  return(q[[design]])
}

# Function to determine the (correct) coefficient matrix Theta
coef_mat_a <- function(p) {0}
coef_mat_b <- function(p) {0}
coef_mat_c <- function(p) {0}
insertSource("simulations/simData.R",
             functions = c("coef_mat_a", "coef_mat_b", "coef_mat_c"))
theta_switch <- function(design, p = 4) {
  if (design %in% c("Diagonal", "Correlated", "HeavyTailed",
                    "Heteroskedastic_y", "Heteroskedastic_eta")) {
    theta <- coef_mat_a(p = p)
  } else if (design == "BlockDiag") {
    theta <- coef_mat_b(p = p)
  } else if (design == "NearBand") {
    theta <- coef_mat_c(p = p)
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
# registerDoRNG(seed = iseed)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = iseed)

# Placeholders
max_ell2_errors <- array(NA, dim = c(nummc, numn, nump, numdes, nummet))
dimnames(max_ell2_errors) <- list(mc = 1:nummc, n = nvec, p = pvec,
                                  design = designs, method = methods)
max_row_sparsities <- array(NA, dim = c(nummc, numn, nump, numdes, nummet))
dimnames(max_row_sparsities) <- list(mc = 1:nummc, n = nvec, p = pvec,
                                     design = designs, method = methods)
num_upd <- array(NA, dim = c(nummc, numn, nump, numdes, 2))
dimnames(num_upd) <- list(mc = 1:nummc, n = nvec, p = pvec,
                          design = designs, method = c("Lasso", "PostLasso"))

intercept <- FALSE # whether to include intercept in simulations

for (this_design in seq_along(designs)) {
  design <- designs[this_design]
  q <- q_switch(design) # set lag length q (correctly)
  for (thisn in seq_along(nvec)) {
    n <- nvec[thisn]
    for (thisp in seq_along(pvec)) {
      p <- pvec[thisp]
      theta <- theta_switch(design, p = p)
      print(paste("Design:", design, ", n:", n, ", p:", p))
      results <- foreach(icount(nummc)) %dopar% {
        # GENERATE DATA
        source("simulations/simData.R", local = TRUE)
        data <- sim_data_by_design(n = n, p = p, design = design, nburn = nburn)
        # ESTIMATE VAR MODEL
        errors <- numeric(nummet)
        shats <- numeric(nummet)
        kterms <- numeric(2)

        # 1-2. LASSO AND POST-LASSO
        # 1. LASSO
        source("lassoVAR.R", local = TRUE) # for lasso_var
        fit_lasso <- lasso_var(data = data, q = q,
                               post = FALSE, intercept = intercept)
        that_lasso <- fit_lasso$that # extract estimates
        errors[1] <- max_ell2_row(that_lasso - theta) # calculate errors
        shats[1] <- max_row_sparsity(that_lasso) # calculate sparsity
        kterms[1] <- fit_lasso$k_term # number of loading updates
        # 2. POST-LASSO
        fit_postl <- lasso_var(data = data, q = q,
                               post = TRUE, intercept = intercept)
        that_postl <- fit_postl$that # extract estimates
        errors[2] <- max_ell2_row(that_postl - theta) # calculate errors
        shats[2] <- max_row_sparsity(that_postl) # calculate sparsity
        kterms[2] <- fit_postl$k_term # number of loading updates

        # 3-4. BIC-LASSO AND POST-BIC-LASSO
        # 3. BIC-LASSO
        source("icLassoVAR.R", local = TRUE) # for bic_lasso_var
        fit_bic <- ic_lasso_var(data = data, q = q, criteria = "bic",
                                post = TRUE, intercept = intercept)
        # Note: post = TRUE to get also Post-BIC-Lasso
        that_bic <- fit_bic$thats[, , "bic"] # extract estimates
        errors[3] <- max_ell2_row(that_bic - theta) # calculate errors
        shats[3] <- max_row_sparsity(that_bic) # calculate sparsity
        # 4. POST-BIC-LASSO
        that_postbic <- fit_bic$thats_post[, , 1] # extract estimates
        errors[4] <- max_ell2_row(that_postbic - theta) # calculate errors
        shats[4] <- max_row_sparsity(that_postbic) # calculate sparsity

        # 5-6. SQRT-LASSO AND POST-SQRT-LASSO
        source("sqrtLassoVAR.R", local = TRUE) # for sqrt_lasso_var
        fit_sqrtl <- sqrt_lasso_var(data = data, q = q,
                                    post = TRUE, intercept = intercept)
        # Note: post = TRUE to get also Post-Sqrt-Lasso
        # 5. SQRT-LASSO
        that_sqrtl <- fit_sqrtl$that
        errors[5] <- max_ell2_row(that_sqrtl - theta)
        shats[5] <- max_row_sparsity(that_sqrtl)
        # 6. POST-SQRT-LASSO
        that_postsqrtl <- fit_sqrtl$that_post
        errors[6] <- max_ell2_row(that_postsqrtl - theta)
        shats[6] <- max_row_sparsity(that_postsqrtl)

        list(errors, shats, kterms) # return results as list
      } # mc loop
      # Extract max rowwise ell_2 errors, max row sparsities
      # and number of loading updates (Lasso and PostLasso), respectively
      max_ell2_errors[, thisn, thisp, this_design, ] <-
        t(sapply(lapply(results, "[[", 1), unlist))
      max_row_sparsities[, thisn, thisp, this_design, ] <-
        t(sapply(lapply(results, "[[", 2), unlist))
      num_upd[, thisn, thisp, this_design, ] <-
        t(sapply(lapply(results, "[[", 3), unlist))
    } # p loop
  } # n loop
} # design loop
stopCluster(cl)

# Save the workspace
if (Sys.info()[["sysname"]] == "Linux") {
  file_name <- paste("simulations_workspace_", nummc, "_MC_",
                     min(nvec), "_to_", max(nvec),
                     "_n_", min(pvec), "_to_", max(pvec), "_p",
                     "_with_num_upd", sep = "")
  save.image(file = paste0("simulations/", file_name, ".RData"))
  q("no")
}