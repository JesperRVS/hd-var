## TODO:
# Use 5 designs 
#   (1) "Diagonal":     Design A as in KC2015
#   (2) "Correlated":   Design A' w/ strongly correlated innovations (hence outcomes)
#   (3) "HeavyTailed":  Design A'' w/ heavy-tailed (here: student-t(5) innovations)
#   (4) "BlockDiag":    Design B as in KC2015
#   (5) "NearBand":     Design C as in KC2015

# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

# setwd("C:/Users/kzb125/Downloads/hd-var") # set working directory
# Note: All sourcing relative to this working directory

# Packages for parallel computing
libpar <- c("doRNG", "doParallel", "foreach")
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
using(libpar)

testrun <- TRUE

# Simulation settings
if (testrun) {
  # nvec <- 100
  nvec <- c(100, 200)
  # nvec <- seq(from = 100, to = 400, by = 100)
  # pvec <- 4
  pvec <- c(4, 8)
  # designs <- c("Diagonal")
  # designs <- c("Diagonal", "BlockDiag")
  designs <- "HeavyTailed"
  methods <- c("Lasso", "PostLasso",
               "AICLasso", "PostAICLasso",
               "BICLasso", "PostBICLasso",
               "HQICLasso", "PostHQICLasso",
               "SqrtLasso", "PostSqrtLasso")
  nburn <- 100
} else {
  nvec <- seq(from = 100, to = 2000, by = 100)
  pvec <- c(16, 32, 64, 128, 256)
  designs <- c("Diagonal", "Correlated", "HeavyTailed", "BlockDiag", "NearBand")
  methods <- c("Lasso", "PostLasso", "AICLasso", "PostAICLasso",
               "BICLasso", "PostBICLasso", "HQICLasso", "PostHQICLasso",
               "SqrtLasso", "PostSqrtLasso")
  nburn <- 10000
}
numn <- length(nvec)
nump <- length(pvec)
numdes <- length(designs)
nummet <- length(methods)

# TODO:
# [x] 0. Write function which takes (design, n, p) as input and generates data
# [x] 1. Generate data depending using (design, n, p)
# 2. Estimate VAR model using method
# 3. Store results in a suitable list
# 4. Parallel computing

## SOME HELPER FUNCTIONS

q_switch <- function(design) {
  q <- list(Diagonal = 1, Correlated = 1, HeavyTailed = 1,
            BlockDiag = 4, NearBand = 1)
  return(q[[design]])
}

# Function to determine the (correct) coefficient matrix Theta
coef_mat_a <- function(p) {0}
coef_mat_b <- function(p) {0}
coef_mat_c <- function(p) {0}
insertSource("simulations/simData.R",
             functions = c("coef_mat_a", "coef_mat_b", "coef_mat_c"))
theta_switch <- function(design, p = 4) {
  if (design %in% c("Diagonal", "Correlated", "HeavyTailed")) {
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

# Monte Carlo (MC) settings
RNGkind(normal.kind = "Kinderman-Ramage") # faster normal draws
iseed <- 2345
# registerDoRNG(seed = iseed)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
clusterSetRNGStream(cl, iseed = iseed)
if (testrun) {
  nummc <- 10   # no. MC repetitions
} else {
  nummc <- 2000 # no. MC repetitions
}

# Placeholders
max_ell2_errors <- array(NA, dim = c(nummc, numn, nump, numdes, nummet))
dimnames(max_ell2_errors) <- list(mc = 1:nummc, n = nvec, p = pvec,
                                  design = designs, method = methods)
max_row_sparsities <- array(NA, dim = c(nummc, numn, nump, numdes, nummet))
dimnames(max_row_sparsities) <- list(mc = 1:nummc, n = nvec, p = pvec,
                                     design = designs, method = methods)

intercept <- FALSE # whether to include intercept in simulations

## MAIN SIMULATION LOOP
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
                                   nburn = nburn)
        # ESTIMATE VAR MODEL
        errors <- numeric(nummet)
        shats <- numeric(nummet)
        # LASSO
        source("lassoVAR.R", local = TRUE) # for lasso_var
        fit_lasso <- lasso_var(data = data, q = q,
                               post = FALSE, intercept = intercept)
        that_lasso <- fit_lasso$that
        errors[1] <- max_ell2_row(that_lasso - theta)
        shats[1] <- max_row_sparsity(that_lasso)
        # print("lasso done")
        # POST-LASSO
        fit_postl <- lasso_var(data = data, q = q,
                               post = TRUE, intercept = intercept)
        that_postl <- fit_postl$that
        errors[2] <- max_ell2_row(that_postl - theta)
        shats[2] <- max_row_sparsity(that_postl)
        # print("postlasso done")
        # IC-LASSOS AND POST-IC-LASSOS
        # IC-LASSOS
        source("icLassoVAR.R", local = TRUE) # for ic_lasso_var
        fit_ics <- ic_lasso_var(data = data, q = q,
                                crit = c("aic", "bic", "hqic"),
                                post = FALSE, intercept = intercept)
        # AIC
        that_aic <- fit_ics$thats[, , 1]
        errors[3] <- max_ell2_row(that_aic - theta)
        shats[3] <- max_row_sparsity(that_aic)
        # print("aiclasso done")
        # BIC
        that_bic <- fit_ics$thats[, , 2]
        errors[5] <- max_ell2_row(that_bic - theta)
        shats[5] <- max_row_sparsity(that_bic)
        # print("biclasso done")
        # HQIC
        that_hqic <- fit_ics$thats[, , 3]
        errors[7] <- max_ell2_row(that_hqic - theta)
        shats[7] <- max_row_sparsity(that_hqic)
        # print("hqiclasso done")
        # POST-IC-LASSOS
        fit_postics <- ic_lasso_var(data = data, q = q,
                                    crit = c("aic", "bic", "hqic"),
                                    post = TRUE, intercept = intercept)
        # POST-AIC
        that_postaic <- fit_postics$thats[, , 1]
        errors[4] <- max_ell2_row(that_postaic - theta)
        shats[4] <- max_row_sparsity(that_postaic)
        # print("postaiclasso done")
        # POST-BIC
        that_postbic <- fit_postics$thats[, , 2]
        errors[6] <- max_ell2_row(that_postbic - theta)
        shats[6] <- max_row_sparsity(that_postbic)
        # print("postbiclasso done")
        # POST-HQIC
        that_posthqic <- fit_postics$thats[, , 3]
        errors[8] <- max_ell2_row(that_posthqic - theta)
        shats[8] <- max_row_sparsity(that_posthqic)
        # source("sqrtLassoVAR.R", local = TRUE)
        # print("posthqiclasso done")
        # SQRT-LASSO
        source("sqrtLassoVAR.R", local = TRUE) # for sqrt_lasso_var
        fit_sqrtl <- sqrt_lasso_var(data = data, q = q,
                                    post = FALSE, intercept = intercept)
                                    # upsilon = matrix(1, p, p * q))
        that_sqrtl <- fit_sqrtl$that
        errors[9] <- max_ell2_row(that_sqrtl - theta)
        shats[9] <- max_row_sparsity(that_sqrtl)
        # print("sqrtlasso done")
        # POST-SQRT-LASSO
        fit_postsqrtl <- sqrt_lasso_var(data = data, q = q,
                                        post = TRUE, intercept = intercept)
                                        # upsilon = matrix(1, p, p * q))
        that_postsqrtl <- fit_postsqrtl$that
        errors[10] <- max_ell2_row(that_postsqrtl - theta)
        shats[10] <- max_row_sparsity(that_postsqrtl)
        # print("postsqrtlasso done")

        list(errors, shats) # return results as list
      } # mc loop
      max_ell2_errors[, thisn, thisp, this_design, ] <-
        t(sapply(lapply(results, "[[", 1), unlist))
      max_row_sparsities[, thisn, thisp, this_design, ] <-
        t(sapply(lapply(results, "[[", 2), unlist))
    } # p loop
  } # n loop
} # design loop
stopCluster(cl)

max_ell2_errors[, 2, 1, design = "HeavyTailed",
                method = c("Lasso", "BICLasso", "SqrtLasso")]
max_ell2_errors[, 1, 1, design = "Diagonal",
                method = methods]



