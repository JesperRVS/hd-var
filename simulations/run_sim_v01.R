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

# Clear
rm(list = ls(all.names = TRUE)) # will clear all (including hidden) objects.
invisible(gc()) #free up memory

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
  nvec <- seq(from = 100, to = 200, by = 100)
  pvec <- c(4, 8)
  designs <- c("Diagonal", "BlockDiag")
  methods <- c("Lasso", "PostLasso")
} else {
  nvec <- seq(from = 100, to = 2000, by = 100)
  pvec <- c(16, 32, 64, 128, 256)
  designs <- c("Diagonal", "Correlated", "HeavyTailed", "BlockDiag", "NearBand")
  methods <- c("Lasso", "PostLasso", "AICLasso", "PostAICLasso",
               "BICLasso", "PostBICLasso", "SqrtLasso", "PostSqrtLasso")
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
# 0. Write function which takes (design, n, p) as input and generates data
# 1. Generate data depending using (design, n, p)
# 2. Estimate VAR model using method
# 3. Store results in a suitable list

# Start without parallel computing

for (design in designs) {
  for (n in nvec) {
    for (p in pvec) {
      for (method in methods) {
        print(paste("Design:", design, "n:", n, "p:", p, "Method:", method))
        results <- foreach(icount(nummc), .options.snow = opts) %dopar% {
          source("helper_functions.R", local = TRUE)
          source("lassoVAR.R", local = TRUE)
          source("icLassoVAR.R", local = TRUE)
          source("sqrtLassoVAR.R", local = TRUE)
        } # mc loop
      } # method loop
    } # p loop
  } # n loop
} # design loop
