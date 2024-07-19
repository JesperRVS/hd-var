#' Simulate data according to a given design
#'
#' @param n Effective sample size
#' @param p System dimension
#' @param design Design type
#' @param sigma_eps Standard deviation of eps_0,i (which are indep. gaussian)
#' @param nburn No. of burn-in periods
#' @return data (q + n) x p outcome matrix (after burn-in)
#' @export
#' @examples
#' sim_data_by_design(n = 100, p = 4, design = "Diagonal", sigma_eps = 0.1,
#'                   seed = 1234, r = 1, nburn = 10000)
#' sim_data_by_design(n = 100, p = 4, design = "Correlated", sigma_eps = 0.1,
#'                  seed = 1234, r = 1, nburn = 10000)
#' sim_data_by_design(n = 100, p = 4, design = "HeavyTailed", sigma_eps = 0.1,
#'                 seed = 1234, r = 1, nburn = 10000)
#' sim_data_by_design(n = 100, p = 4, design = "BlockDiag", sigma_eps = 0.1,
#'                seed = 1234, r = 1, nburn = 10000)
#' sim_data_by_design(n = 100, p = 4, design = "NearBand", sigma_eps = 0.1,
#'               seed = 1234, r = 1, nburn = 10000)
# sim_data_by_design <- function(n = 100, p = 4, design, sigma_eps = 0.1,
#                                nburn = 10000) {
#   switch(design,
#     # Diagonal design w/ independent Gaussian innovations
#     "Diagonal" = {
#       source("design_A.R")
#       data <- sim_data_a(n = n, p = p, family = "gaussian",
#                          sigma_eps = sigma_eps, rho = 0, nburn = nburn)
#     },
#     # Diagonal design w/ correlated Gaussian innovations
#     "Correlated" = {
#       source("design_A.R")
#       data <- sim_data_a(n = n, p = p, family = "gaussian",
#                          sigma_eps = sigma_eps, rho = 0.5, nburn = nburn)
#     },
#     # Diagonal design w/ independent Student-t innovations
#     "HeavyTailed" = {
#       source("design_A.R")
#       data <- sim_data_a(n = n, p = p, family = "student",
#                          sigma_eps = sigma_eps, df = 5, nburn = nburn)
#     },
#     # Block diagonal design w/ independent Gaussian innovations
#     "BlockDiag" = {
#       source("design_B.R")
#       data <- sim_data_b(n = n, p = p, sigma_eps = sigma_eps, nburn = nburn)
#     },
#     # Near-band design w/ independent Gaussian innovations
#     "NearBand" = {
#       source("design_C.R")
#       data <- sim_data_c(n = n, p = p, sigma_eps = sigma_eps, nburn = nburn)
#     }
#   ) # end switch
#   return(data)
# }
# # rewrite the above function using double spaces for indentation and
# # keeping lines to a maximum of 80 characters
sim_data_by_design <- function(n = 100, p = 4, design, sigma_eps = 0.1,
                               nburn = 10000) {
  switch(design,
    # Diagonal design w/ independent Gaussian innovations
    "Diagonal" = {
      source("simulations/design_A.R")
      data <- sim_data_a(n = n, p = p, family = "gaussian",
                         sigma_eps = sigma_eps, rho = 0, nburn = nburn)
    },
    # Diagonal design w/ correlated Gaussian innovations
    "Correlated" = {
      source("simulations/design_A.R")
      data <- sim_data_a(n = n, p = p, family = "gaussian",
                         sigma_eps = sigma_eps, rho = 0.5, nburn = nburn)
    },
    # Diagonal design w/ independent Student-t innovations
    "HeavyTailed" = {
      source("simulations/design_A.R")
      data <- sim_data_a(n = n, p = p, family = "student",
                         sigma_eps = sigma_eps, df = 5, nburn = nburn)
    },
    # Block diagonal design w/ independent Gaussian innovations
    "BlockDiag" = {
      source("simulations/design_B.R")
      data <- sim_data_b(n = n, p = p, sigma_eps = sigma_eps, nburn = nburn)
    },
    # Near-band design w/ independent Gaussian innovations
    "NearBand" = {
      source("simulations/design_C.R")
      data <- sim_data_c(n = n, p = p, sigma_eps = sigma_eps, nburn = nburn)
    }
  ) # end switch
  return(data)
}