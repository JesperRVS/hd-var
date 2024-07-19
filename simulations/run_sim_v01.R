## TODO:
# 1. Create main_sim file
# 2. Use 5 designs 
#   (1) "Diagonal":     Design A as in KC2015
#   (2) "Correlated":   Design A' w/ strongly correlated innovations (hence outcomes)
#   (3) "HeavyTailed":  Design A'' w/ heavy-tailed (here: student-t(5) innovations)
#   (4) "BlockDiag":    Design B as in KC2015
#   (5) "NearBand":     Design C as in KC2015
# 3. Use the following estimation methods (1) Lasso w/ Lasso updating (2)
#   Post-Lasso w/ post-Lasso updating (3) AIC-Lasso (4) Post-AIC-Lasso (5)
#   BIC-Lasso (6) Post-BIC-Lasso (7) Sqrt-Lasso (8) Post-Sqrt-Lasso (9) Least
#   squares w/ Moore-Penrose inverse
# 4. See if you can store all matrices
# 5. FIGURE OUT HOW TO EXECUTE FROM SIMS DIRECTORY. sourcing is a pain.

# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

setwd("C:/Users/kzb125/Downloads/hd-var") # set working directory
# Note: All sourcing relative to this working directory

# Packages for parallel computing
libpar <- c("foreach", "iterators", "parallel", "doParallel") 
# c("formattable", "latex2exp", "stringr")
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
  nvec <- 100
  # nvec <- seq(from = 100, to = 200, by = 100)
  pvec <- 4
  # pvec <- c(4, 8)
  designs <- c("Diagonal")
  # designs <- c("Diagonal", "BlockDiag")
  methods <- "Lasso"
  # methods <- c("Lasso", "PostLasso")
  nburn <- 100
} else {
  nvec <- seq(from = 100, to = 2000, by = 100)
  pvec <- c(16, 32, 64, 128, 256)
  designs <- c("Diagonal", "Correlated", "HeavyTailed", "BlockDiag", "NearBand")
  methods <- c("Lasso", "PostLasso", "AICLasso", "PostAICLasso",
               "BICLasso", "PostBICLasso", "SqrtLasso", "PostSqrtLasso")
  nburn <- 10000
}
numn <- length(nvec)
nump <- length(pvec)
numdes <- length(designs)
nummet <- length(methods)

# Monte Carlo (MC) settings
RNGkind(normal.kind = "Kinderman-Ramage") # faster normal draws
cl <- makeCluster(detectCores())
registerDoParallel(cl)
opts <- list(preschedule = TRUE)
seed <- 2345
clusterSetRNGStream(cl, seed) # seed (for reproducibility)
if (testrun) {
  nummc <- 10   # no. MC repetitions
} else {
  nummc <- 2000 # no. MC repetitions
}

# TODO:
# [x] 0. Write function which takes (design, n, p) as input and generates data
# [x] 1. Generate data depending using (design, n, p)
# 2. Estimate VAR model using method
# 3. Store results in a suitable list
# Start without parallel computing

## SOME HELPER FUNCTIONS

q_switch <- function(design) {
  q <- list(Diagonal = 1, Correlated = 1, HeavyTailed = 1,
            BlockDiag = 4, NearBand = 1)
  return(q[[design]])
}

# Function to determine the (correct) coefficient matrix Theta
source("simulations/coef_mats.R") # for coef_mat_a, coef_mat_b, coef_mat_c
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


## MAIN SIMULATION LOOP
# source("helper_functions.R") # for all methods
for (design in designs) {
  q <- q_switch(design) # set lag length q (correctly)
  for (n in nvec) {
    for (p in pvec) {
      theta <- theta_switch(design, p = p)
      for (method in methods) {
        print(paste("Design:", design, ", n:", n, ", p:", p,
                    ", Method:", method))
        results <- foreach(icount(nummc), .options.snow = opts) %do% {
          # GENERATE DATA
          source("simulations/sim_data_by_design.R", local = TRUE)
          data <- sim_data_by_design(n = n, p = p, design = design,
                                     nburn = nburn)
          # # ESTIMATE VAR MODEL
          source("lassoVAR.R", local = TRUE) # for lasso_var
          fit_lasso <- lasso_var(data = data, q = q)
          that_lasso <- fit_lasso$that
          max_ell2_lasso <- max_ell2_row(that_lasso - theta)
          # source("icLassoVAR.R", local = TRUE)
          # source("sqrtLassoVAR.R", local = TRUE)
          list(max_ell2_lasso)
        } # mc loop
      } # method loop
    } # p loop
  } # n loop
} # design loop
stopCluster(cl)
sapply(lapply(results, "[[", 1), unlist)


