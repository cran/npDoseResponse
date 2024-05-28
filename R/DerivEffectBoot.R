#' Nonparametric bootstrap inference on the derivative effect via our localized derivative estimator.
#'
#' This function implements the nonparametric bootstrap inference on the derivative
#' of a dose-response curve via our localized derivative estimator.
#'
#' @param Y The input n-dimensional outcome variable vector.
#' @param X The input n*(d+1) matrix. The first column of X stores the
#' treatment/exposure variables, while the other d columns are confounding variables.
#' @param t_eval The m-dimensional vector for evaluating the derivative. (Default:
#' t_eval = NULL. Then, t_eval = \code{X[,1]}, which consists of the observed treatment
#' variables.)
#' @param boot_num The number of bootstrapping times. (Default: boot_num = 500.)
#' @param alpha The confidence level of both the uniform confidence band and
#' pointwise confidence interval. (Default: alpha = 0.95.)
#' @param h_bar The bandwidth parameter for the Nadaraya-Watson conditional
#' CDF estimator. (Default: h_bar = NULL. Then, the Silverman's rule of thumb
#' is applied. See Chen et al. (2016) for details.)
#' @param kernT_bar The name of the kernel function for the Nadaraya-Watson
#' conditional CDF estimator. (Default: kernT_bar = "gaussian".)
#' @param h,b The bandwidth parameters for the treatment/exposure variable
#' and confounding variables in the local polynomial regression.
#' (Default: h = NULL, b = NULL. Then, the rule-of-thumb bandwidth selector
#' in Eq. (A1) of Yang and Tschernig (1999) is used with additional scaling
#' factors C_h and C_b, respectively.)
#' @param C_h,C_b The scaling factors for the rule-of-thumb bandwidth parameters.
#' @param print_bw The indicator of whether the current bandwidth parameters
#' should be printed to the console. (Default: print_bw = TRUE.)
#' @param degree Degree of local polynomials. (Default: degree = 2.)
#' @param deriv_ord The order of the estimated derivative of the conditional mean
#' outcome function. (Default: deriv_ord = 1. It shouldn't be changed in most cases.)
#' @param kernT,kernS The names of kernel functions for the treatment/exposure
#' variable and confounding variables. (Default: kernT = "epanechnikov",
#' kernS = "epanechnikov".)
#' @param parallel The indicator of whether the function should be parallel
#' executed. (Default: parallel = TRUE.)
#' @param cores The number of cores for parallel execution. (Default: cores = 6.)
#'
#' @return A list that contains four elements.
#' \item{theta_est}{The estimated derivative of the dose-response curve evaluated at points \code{t_eval}.}
#' \item{theta_est_boot}{The estimated derivative of the dose-response curve evaluated at points \code{t_eval} for all the bootstrap samples.}
#' \item{theta_alpha}{The width of the uniform confidence band.}
#' \item{theta_alpha_var}{The widths of the pointwise confidence bands at evaluation points \code{t_eval}.}
#'
#' @author Yikun Zhang, \email{yikunzhang@@foxmail.com}
#' @references Zhang, Y., Chen, Y.-C., and Giessing, A. (2024)
#' \emph{Nonparametric Inference on Dose-Response Curves Without the Positivity Condition.}
#' \url{https://arxiv.org/abs/2405.09003}.
#'
#' @keywords estimator derivative localized our through curve dose-response a of derivative the on inference bootstrap
#'
#' @examples
#' \donttest{
#'   set.seed(123)
#'   n <- 300
#'
#'   S2 <- cbind(2 * runif(n) - 1, 2 * runif(n) - 1)
#'   Z2 <- 4 * S2[, 1] + S2[, 2]
#'   E2 <- 0.2 * runif(n) - 0.1
#'   T2 <- cos(pi * Z2^3) + Z2 / 4 + E2
#'   Y2 <- T2^2 + T2 + 10 * Z2 + rnorm(n, mean = 0, sd = 1)
#'   X2 <- cbind(T2, S2)
#'
#'   t_qry2 = seq(min(T2) + 0.01, max(T2) - 0.01, length.out = 100)
#'   chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
#'   if (nzchar(chk) && chk == "TRUE") {
#'     # use 2 cores in CRAN/Travis/AppVeyor
#'     num_workers <- 2L
#'   } else {
#'     # use all cores in devtools::test()
#'     num_workers <- parallel::detectCores()
#'   }
#'
#'   # Increase bootstrap times "boot_num" to a larger integer in practice
#'   theta_boot2 = DerivEffectBoot(Y2, X2, t_eval = t_qry2, boot_num = 3, alpha = 0.95,
#'                                h_bar = NULL, kernT_bar = "gaussian", h = NULL,
#'                                b = NULL, C_h = 7, C_b = 3, print_bw = FALSE,
#'                                degree = 2, deriv_ord = 1, kernT = "epanechnikov",
#'                                kernS = "epanechnikov", parallel = TRUE,
#'                                cores = num_workers)
#' }
#'
#' @export
#' @importFrom stats quantile
#'

DerivEffectBoot <- function(Y, X, t_eval = NULL, boot_num = 500, alpha = 0.95,
                            h_bar = NULL, kernT_bar = "gaussian", h = NULL,
                            b = NULL, C_h = 7, C_b = 3, print_bw = TRUE,
                            degree = 2, deriv_ord = 1, kernT = "epanechnikov",
                            kernS = "epanechnikov", parallel = TRUE, cores = 6) {
  if (is.null(t_eval)) {
    t_eval <- X[, 1]
  }

  n <- nrow(X)

  theta_est <- DerivEffect(Y, X, t_eval = t_eval, h_bar = h_bar,
                           kernT_bar = kernT_bar, h = h, b = b, C_h = C_h,
                           C_b = C_b, print_bw = print_bw, degree = degree,
                           deriv_ord = deriv_ord, kernT = kernT, kernS = kernS,
                           parallel = parallel, cores = cores)

  theta_est_boot <- matrix(0, nrow = boot_num, ncol = length(t_eval))
  b <- 1
  while (b <= boot_num) {
    ind <- sample(n, n, replace = TRUE)
    X_boot <- X[ind,]
    Y_boot <- Y[ind]
    theta_est_boot[b, ] <- DerivEffect(Y_boot, X_boot, t_eval = t_eval, h_bar = h_bar,
                                     kernT_bar = kernT_bar, h = h, b = b, C_h = C_h,
                                     C_b = C_b, print_bw = print_bw, degree = degree,
                                     deriv_ord = deriv_ord, kernT = kernT, kernS = kernS,
                                     parallel = parallel, cores = cores)
    if (sum(is.nan(theta_est_boot[b, ])) == 0) {
      b <- b + 1
    }
  }

  # Compute the alpha% uniform confidence bands
  theta_boot_sup <- apply(abs(sweep(theta_est_boot, 2, theta_est)), 1, max)
  theta_alpha <- quantile(theta_boot_sup, alpha, na.rm = TRUE)

  # Compute the alpha% pointwise confidence intervals
  theta_boot_abs <- abs(sweep(theta_est_boot, 2, theta_est))
  theta_alpha_var <- apply(theta_boot_abs, 2, quantile, probs = alpha, na.rm = TRUE)

  # Return the results
  return(list(theta_est = theta_est, theta_est_boot = theta_est_boot,
              theta_alpha = theta_alpha, theta_alpha_var = theta_alpha_var))
}
