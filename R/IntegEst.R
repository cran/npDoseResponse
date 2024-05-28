#' The proposed integral estimator.
#'
#' This function implements our proposed integral estimator for estimating the
#' dose-response curve.
#'
#' @param Y The input n-dimensional outcome variable vector.
#' @param X The input n*(d+1) matrix. The first column of X stores the
#' treatment/exposure variables, while the other d columns are confounding variables.
#' @param t_eval The m-dimensional vector for evaluating the dose-response curve.
#' (Default: t_eval = NULL. Then, t_eval = \code{X[,1]}, which consists of the observed
#' treatment variables.)
#' @param h_bar The bandwidth parameter for the Nadaraya-Watson conditional
#' CDF estimator. (Default: h_bar = NULL. Then, the Silverman's rule of thumb
#' is applied. See Chen et al. (2016) for details.)
#' @param kernT_bar The name of the kernel function for the Nadaraya-Watson conditional
#' CDF estimator. (Default: "gaussian".)
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
#' @return The estimated dose-response curve evaluated at points \code{t_eval}.
#'
#' @author Yikun Zhang, \email{yikunzhang@@foxmail.com}
#' @references Zhang, Y., Chen, Y.-C., and Giessing, A. (2024)
#' \emph{Nonparametric Inference on Dose-Response Curves Without the Positivity Condition.}
#' \url{https://arxiv.org/abs/2405.09003}.
#'
#' @keywords curve dose-response the of estimator integral
#'
#' @examples
#' \donttest{
#'   set.seed(123)
#'   n <- 300
#'
#'   S2 <- cbind(2*runif(n) - 1, 2*runif(n) - 1)
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
#'   m_est2 = IntegEst(Y2, X2, t_eval = t_qry2, h_bar = NULL, kernT_bar = "gaussian",
#'                     h = NULL, b = NULL, C_h = 7, C_b = 3, print_bw = FALSE,
#'                     degree = 2, deriv_ord = 1, kernT = "epanechnikov",
#'                     kernS = "epanechnikov", parallel = TRUE, cores = num_workers)
#'
#'   plot(t_qry2, m_est2, type="l", col = "blue", xlab = "t", lwd=5,
#'        ylab="(Estimated) dose-response curves")
#'   lines(t_qry2, t_qry2^2 + t_qry2, col = "red", lwd=3)
#'   legend(-2, 6, legend=c("Estimated curve", "True curve"), fill = c("blue","red"))
#' }
#'
#' @export
#' @importFrom stats approx
#'

IntegEst <- function(Y, X, t_eval = NULL, h_bar = NULL, kernT_bar = "gaussian",
                     h = NULL, b = NULL, C_h = 7, C_b = 3, print_bw = TRUE,
                     degree = 2, deriv_ord = 1, kernT = "epanechnikov",
                     kernS = "epanechnikov", parallel = TRUE, cores = 6) {
  if (is.null(t_eval)) {
    t_eval <- X[, 1]
  }

  T_sort <- sort(X[, 1])
  n <- nrow(X)

  theta_est <- DerivEffect(Y, X, t_eval = T_sort, h_bar = h_bar, kernT_bar = kernT_bar,
                           h = h, b = b, C_h = C_h, C_b = C_b, print_bw = print_bw,
                           degree = degree, deriv_ord = deriv_ord, kernT = kernT,
                           kernS = kernS, parallel = parallel, cores = cores)

  T_delta <- T_sort[2:n] - T_sort[1:(n-1)]

  int_mat_up <- matrix(T_delta * (1:(n-1)) * theta_est[1:(n-1)], nrow = n-1, ncol = n)
  int_mat_up <- int_mat_up * outer(1:(n-1), 1:n, "<")

  int_mat_down <- matrix(T_delta * (n - 1:(n-1)) * theta_est[2:n], nrow = n-1, ncol = n)
  int_mat_down <- int_mat_down * outer(1:(n-1), 1:n, ">=")

  m_samp <- mean(Y) + colSums(int_mat_up - int_mat_down) / n

  m_est <- approx(T_sort, m_samp, xout = t_eval, method = "linear")$y

  return(m_est)
}
