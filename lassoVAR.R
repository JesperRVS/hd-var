#' Lasso VAR Model Fitting
#'
#' This function fits a vector autoregressive (VAR) model with Lasso
#' regularization and provides options for refining the estimates.
#'
#' @param data A (q + n) x p matrix of data (n = effective sample size).
#' @param q Autoregressive order; default is 1 (i.e., VAR(1) model).
#' @param post Logical; if TRUE, refit estimates after each loading update.
#' @param intercept Logical; if TRUE, estimate intercepts.
#' @param c Tuning parameter for penalty level (score markup). Default is 1.1,
#'   as in Belloni et al. (2012, ECTA).
#' @param gamma Tuning parameter for penalty level (score quantile). Default is
#'   0.1 / log(max(dim(data))), as in Belloni et al. (2012, ECTA).
#' @param k Maximum number of penalty loading updates. Default is 15, as in
#'   Belloni et al. (2012, ECTA).
#' @param tol_ups Convergence tolerance for penalty loadings (relative), so
#'   that penalty loadings are considered converged if
#'   ||vec(ups_new - ups_old)||_2 / ||vec(ups_old)||_2 <= tol_ups.
#' @param warn Logical; if TRUE, issue warning if max updates reached.
#' @param full_path Logical; if TRUE, use full path of lambda values. Default
#'   is FALSE.
#' @param tol_glmnet Convergence tolerance for glmnet. (Mainly for debugging.)
#'
#' @return A list with the following components:
#' \item{lambda}{Penalty level used}
#' \item{ups_init}{Initial penalty loadings}
#' \item{ups_refi}{Refined penalty loadings path}
#' \item{ups}{Final penalty loadings}
#' \item{intr_init}{Initial intercepts (if requested)}
#' \item{intr_refi}{Refined intercepts path (if requested)}
#' \item{intr}{Final intercepts (if requested)}
#' \item{that_init}{Initial estimates}
#' \item{that_refi}{Refined estimates path}
#' \item{that}{Final estimates}
#' \item{full_rank_post_init}{Logicals; if TRUE, selected regressors full rank}
#' \item{full_rank_post_refi}{Logicals; if TRUE, selected regressors full rank}
#' \item{k_term}{Number of penalty loading updates performed (up to k)}
#' \item{rel_diff_ups_term}{Relative change in penalty loadings at termination}
#' \item{rel_diffs_ups}{Relative changes in penalty loadings along path}
#'
#' @references
#' Belloni, A., Chernozhukov, V., Chen, D., & Hansen, C. (2012).
#' Sparse models and methods for optimal instruments with an application to 
#' eminent domain. Econometrica, 80(6), 2369-2429.
#' \url{https://doi.org/10.3982/ECTA9626}
#'
#' Kock, A. B., Pedersen, R. S., & Sørensen, J. R.-V. (2024). Data-Driven Tuning
#' Parameter Selection for High-Dimensional Vector Autoregressions.
#' \url{https://arxiv.org/abs/2403.06657}
#'
#' @importFrom stats qnorm
#' @importFrom glmnet glmnet
#' @importFrom Matrix crossprod
#'
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(1000), 100, 10)
#' fit <- lasso_var(data)
#' }
#'
#' @export
source("helper_functions.R")
lasso_var <- function(data, q = 1, post = TRUE, intercept = TRUE,
                      c = 1.1, gamma = NULL, k = 15,
                      tol_ups = 1e-3, warn = TRUE, full_path = FALSE,
                      tol_glmnet = 1e-4) {
  # Unpack data and construct predictors and response
  xy_data <- unpack(data = data, q = q)
  x <- xy_data$x
  y <- xy_data$y
  # Call upon mult_lasso_bcch
  fit <- mult_lasso_bcch(x = x, y = y, post = post, intercept = intercept,
                         c = c, gamma = gamma, k = k,
                         tol_ups = tol_ups, warn = warn, full_path = full_path,
                         tol_glmnet = tol_glmnet)
  return(fit)
}