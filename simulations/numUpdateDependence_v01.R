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

n <- 500
p <- 16
k <- 15
design <- "Diagonal"
nburn <- 100
nummc <- 1000 # no. MC repetitions

do_save <- TRUE # whether to save the workspace

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

theta <- theta_switch(design, p, n)

# Function which calculates the maximum ell_2 norm of the rows of a matrix
max_ell2_row <- function(mat) {
  max_ell2 <- max(apply(mat, 1, function(x) sqrt(sum(x^2))))
  return(max_ell2)
}

## == SIMULATIONS == ##

# Monte Carlo (MC) settings
RNGkind(normal.kind = "Kinderman-Ramage") # faster normal draws
iseed <- 2345 # set seed for reproducibility
set.seed(iseed)


# Placeholders for results
lasso_k_term <- numeric(nummc)
lasso_rel_change_loadings <- array(NA, dim = c(nummc, k))
lasso_rel_change_estimates <- array(NA, dim = c(nummc, k))
lasso_max_ell2_errors_init <- numeric(nummc)
lasso_max_ell2_errors_refi <- array(NA, dim = c(nummc, k))
postl_k_term <- numeric(nummc)
postl_rel_change_loadings <- array(NA, dim = c(nummc, k))
postl_rel_change_estimates <- array(NA, dim = c(nummc, k))
postl_max_ell2_errors_init <- numeric(nummc)
postl_max_ell2_errors_refi <- array(NA, dim = c(nummc, k))

intercept <- FALSE # whether to include intercept in simulations

for (mc in 1:nummc) {
  # GENERATE DATA
  source("simulations/simData.R", local = TRUE)
  data <- sim_data_by_design(n = n, p = p, design = design, nburn = nburn)

  # 1. LASSO
  source("lassoVAR.R", local = TRUE) # for lasso_var
  fit_lasso <- lasso_var(data = data, q = q_switch(design), post = FALSE,
                         intercept = intercept, k = k)
  lasso_max_ell2_errors_init[mc] <- max_ell2_row(fit_lasso$that_init - theta)
  lasso_k_term[mc] <- fit_lasso$k_term # number of updates
  # (1.a) Relative change in loadings as a function of the number of updates
  lasso_ups_init <- fit_lasso$ups_init
  lasso_ups_refi <- fit_lasso$ups_refi
  diff_init <- lasso_ups_refi[, , 1] - lasso_ups_init
  lasso_rel_change_loadings[mc, 1] <-
    sqrt(sum((diff_init)^2)) / (sqrt(sum(lasso_ups_init^2))
                                + .Machine$double.eps)
  for (l in 2:fit_lasso$k_term) {
    ups_old <- lasso_ups_refi[, , l - 1]
    ups_new <- lasso_ups_refi[, , l]
    diff_ups <- ups_new - ups_old
    lasso_rel_change_loadings[mc, l] <-
      sqrt(sum((diff_ups)^2)) / (sqrt(sum(ups_old^2)) + .Machine$double.eps)
  }
  if (fit_lasso$k_term < k) {
    lasso_rel_change_loadings[mc, (fit_lasso$k_term + 1):k] <- 0
  }
  # (1.b) Relative change in estimates as a function of the number of updates
  lasso_that_init <- fit_lasso$that_init
  lasso_that_refi <- fit_lasso$that_refi
  diff_init <- lasso_that_refi[, , 1] - lasso_that_init
  lasso_rel_change_estimates[mc, 1] <-
    sqrt(sum((diff_init)^2)) / (sqrt(sum(lasso_that_init^2))
                                + .Machine$double.eps)
  for (l in 2:fit_lasso$k_term) {
    that_old <- lasso_that_refi[, , l - 1]
    that_new <- lasso_that_refi[, , l]
    diff_that <- that_new - that_old
    lasso_rel_change_estimates[mc, l] <-
      sqrt(sum((diff_that)^2)) / (sqrt(sum(that_old^2)) + .Machine$double.eps)
  }
  if (fit_lasso$k_term < k) {
    lasso_rel_change_estimates[mc, (fit_lasso$k_term + 1):k] <- 0
  }
  # (1.c) Estimation error as a function of the number of updates
  for (l in 1:fit_lasso$k_term) {
    lasso_max_ell2_errors_refi[mc, l] <-
      max_ell2_row(fit_lasso$that_refi[, , l] - theta)
  }
  if (fit_lasso$k_term < k) {
    lasso_max_ell2_errors_refi[mc, (fit_lasso$k_term + 1):k] <-
      lasso_max_ell2_errors_refi[mc, fit_lasso$k_term]
  }

  # 2. POST-LASSO
  fit_postl <- lasso_var(data = data, q = q_switch(design), post = TRUE,
                         intercept = intercept, k = k)
  postl_max_ell2_errors_init[mc] <- max_ell2_row(fit_postl$that_init - theta)
  postl_k_term[mc] <- fit_postl$k_term # number of updates
  # (2.a) Relative change in loadings as a function of the number of updates
  postl_ups_init <- fit_postl$ups_init
  postl_ups_refi <- fit_postl$ups_refi
  diff_init <- postl_ups_refi[, , 1] - postl_ups_init
  postl_rel_change_loadings[mc, 1] <-
    sqrt(sum((diff_init)^2)) / (sqrt(sum(postl_ups_init^2))
                                + .Machine$double.eps)
  for (l in 2:fit_postl$k_term) {
    ups_old <- postl_ups_refi[, , l - 1]
    ups_new <- postl_ups_refi[, , l]
    diff_ups <- ups_new - ups_old
    postl_rel_change_loadings[mc, l] <-
      sqrt(sum((diff_ups)^2)) / (sqrt(sum(ups_old^2)) + .Machine$double.eps)
  }
  if (fit_postl$k_term < k) {
    postl_rel_change_loadings[mc, (fit_postl$k_term + 1):k] <- 0
  }
  # (2.b) Relative change in estimates as a function of the number of updates
  postl_that_init <- fit_postl$that_init
  postl_that_refi <- fit_postl$that_refi
  diff_init <- postl_that_refi[, , 1] - postl_that_init
  postl_rel_change_estimates[mc, 1] <-
    sqrt(sum((diff_init)^2)) / (sqrt(sum(postl_that_init^2))
                                + .Machine$double.eps)
  for (l in 2:fit_postl$k_term) {
    that_old <- postl_that_refi[, , l - 1]
    that_new <- postl_that_refi[, , l]
    diff_that <- that_new - that_old
    postl_rel_change_estimates[mc, l] <-
      sqrt(sum((diff_that)^2)) / (sqrt(sum(that_old^2)) + .Machine$double.eps)
  }
  if (fit_postl$k_term < k) {
    postl_rel_change_estimates[mc, (fit_postl$k_term + 1):k] <- 0
  }
  # (2.c) Estimation error as a function of the number of updates
  for (l in 1:fit_postl$k_term) {
    postl_max_ell2_errors_refi[mc, l] <-
      max_ell2_row(fit_postl$that_refi[, , l] - theta)
  }
  if (fit_postl$k_term < k) {
    postl_max_ell2_errors_refi[mc, (fit_postl$k_term + 1):k] <-
      postl_max_ell2_errors_refi[mc, fit_postl$k_term]
  }
}

if (do_save) {
  file_name <- paste("numupd_dependence_workspace_", nummc, "_MC_",
                     n, "_n_", p, "_p_", k, "_k_", design, "_design", sep = "")
  save.image(file = paste0("simulations/", file_name, ".RData"))
}