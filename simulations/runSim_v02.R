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

h <- 3 # degree of heteroskedasticity (> 0) for heteroskedastic design only

# Simulation settings
if (testrun) {
  nvec <- seq(from = 100, to = 500, by = 100)
  pvec <- c(4, 8, 16)
  # designs <- c("Diagonal", "BlockDiag")
  # designs <- c("Diagonal", "Heteroskedastic")
  designs <- c("Diagonal", "Correlated", "HeavyTailed",
               "BlockDiag", "NearBand", "Heteroskedastic")
  methods <- c("Lasso", "PostLasso",
               "AICLasso", "PostAICLasso",
               "BICLasso", "PostBICLasso",
               "HQICLasso", "PostHQICLasso",
               "SqrtLasso", "PostSqrtLasso")
  nburn <- 100
  nummc <- 10   # no. MC repetitions
} else {
  # nvec <- seq(from = 100, to = 500, by = 100)
  nvec <- seq(from = 200, to = 1000, by = 200)
  # nvec <- seq(from = 100, to = 1000, by = 100) # TODO
  pvec <- c(16, 32, 64, 128)
  # pvec <- c(16, 32, 64, 128, 256) # TODO
  designs <- c("Diagonal", "Correlated", "HeavyTailed",
               "BlockDiag", "NearBand", "Heteroskedastic")
  methods <- c("Lasso", "PostLasso",
               "AICLasso", "PostAICLasso",
               "BICLasso", "PostBICLasso",
               "HQICLasso", "PostHQICLasso",
               "SqrtLasso", "PostSqrtLasso")
  nburn <- 1000
  # nummc <- 80 # no. MC repetitions
  nummc <- 250 # no. MC repetitions TODO
}
numn <- length(nvec)
nump <- length(pvec)
numdes <- length(designs)
nummet <- length(methods)

## == SOME HELPER FUNCTIONS == ##

q_switch <- function(design) {
  q <- list(Diagonal = 1, Correlated = 1, HeavyTailed = 1,
            BlockDiag = 4, NearBand = 1, Heteroskedastic = 1)
  return(q[[design]])
}

# Function to determine the (correct) coefficient matrix Theta
coef_mat_a <- function(p) {0}
coef_mat_b <- function(p) {0}
coef_mat_c <- function(p) {0}
insertSource("simulations/simData.R",
             functions = c("coef_mat_a", "coef_mat_b", "coef_mat_c"))
theta_switch <- function(design, p = 4) {
  if (design %in% c("Diagonal", "Correlated",
                    "HeavyTailed", "Heteroskedastic")) {
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
        data <- sim_data_by_design(n = n, p = p, design = design,
                                   h = h, nburn = nburn)
        # ESTIMATE VAR MODEL
        errors <- numeric(nummet)
        shats <- numeric(nummet)
        kterms <- numeric(2)

        # 1-2. LASSO AND POST-LASSO
        # 1. LASSO
        source("lassoVAR.R", local = TRUE) # for lasso_var
        fit_lasso <- lasso_var(data = data, q = q,
                               post = FALSE, intercept = intercept)
        that_lasso <- fit_lasso$that
        errors[1] <- max_ell2_row(that_lasso - theta)
        shats[1] <- max_row_sparsity(that_lasso)
        kterms[1] <- fit_lasso$k_term
        # 2. POST-LASSO
        fit_postl <- lasso_var(data = data, q = q,
                               post = TRUE, intercept = intercept)
        that_postl <- fit_postl$that
        errors[2] <- max_ell2_row(that_postl - theta)
        shats[2] <- max_row_sparsity(that_postl)
        kterms[2] <- fit_postl$k_term

        # 3-8. IC-LASSOS AND POST-IC-LASSOS
        # 3/5/7. IC-LASSOS
        source("icLassoVAR.R", local = TRUE) # for ic_lasso_var
        fit_ics <- ic_lasso_var(data = data, q = q,
                                crit = c("aic", "bic", "hqic"),
                                post = FALSE, intercept = intercept)
        # 3. AIC
        that_aic <- fit_ics$thats[, , 1]
        errors[3] <- max_ell2_row(that_aic - theta)
        shats[3] <- max_row_sparsity(that_aic)
        # 5. BIC
        that_bic <- fit_ics$thats[, , 2]
        errors[5] <- max_ell2_row(that_bic - theta)
        shats[5] <- max_row_sparsity(that_bic)
        # 7. HQIC
        that_hqic <- fit_ics$thats[, , 3]
        errors[7] <- max_ell2_row(that_hqic - theta)
        shats[7] <- max_row_sparsity(that_hqic)
        # 4/6/8. POST-IC-LASSOS
        fit_postics <- ic_lasso_var(data = data, q = q,
                                    crit = c("aic", "bic", "hqic"),
                                    post = TRUE, intercept = intercept)
        # 4. POST-AIC
        that_postaic <- fit_postics$thats[, , 1]
        errors[4] <- max_ell2_row(that_postaic - theta)
        shats[4] <- max_row_sparsity(that_postaic)
        # 6. POST-BIC
        that_postbic <- fit_postics$thats[, , 2]
        errors[6] <- max_ell2_row(that_postbic - theta)
        shats[6] <- max_row_sparsity(that_postbic)
        # 8. POST-HQIC
        that_posthqic <- fit_postics$thats[, , 3]
        errors[8] <- max_ell2_row(that_posthqic - theta)
        shats[8] <- max_row_sparsity(that_posthqic)

        # 9-10. SQRT-LASSO AND POST-SQRT-LASSO
        # 9. SQRT-LASSO
        source("sqrtLassoVAR.R", local = TRUE) # for sqrt_lasso_var
        fit_sqrtl <- sqrt_lasso_var(data = data, q = q,
                                    post = FALSE, intercept = intercept)
        that_sqrtl <- fit_sqrtl$that
        errors[9] <- max_ell2_row(that_sqrtl - theta)
        shats[9] <- max_row_sparsity(that_sqrtl)
        # 10. POST-SQRT-LASSO
        fit_postsqrtl <- sqrt_lasso_var(data = data, q = q,
                                        post = TRUE, intercept = intercept)
        that_postsqrtl <- fit_postsqrtl$that
        errors[10] <- max_ell2_row(that_postsqrtl - theta)
        shats[10] <- max_row_sparsity(that_postsqrtl)
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
                     "_n_", min(pvec), "_to_", max(pvec), "_p_",
                     h, "_h", sep = "")
  save.image(file = paste0("simulations/", file_name, ".RData"))
  q("no")
}