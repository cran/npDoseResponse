#' The (partial) local polynomial regression.
#'
#' This function implements the (partial) local polynomial regression for estimating
#' the conditional mean outcome function and its partial derivatives. We use
#' higher-order local monomials for the treatment variable and first-order local
#' monomials for the confounding variables.
#'
#' @param Y The input n-dimensional outcome variable vector.
#' @param X The input n*(d+1) matrix. The first column of X stores the
#' treatment/exposure variables, while the other d columns are confounding variables.
#' @param x_eval The n*(d+1) matrix for evaluating the local polynomial regression
#' estimates. (Default: x_eval = NULL. Then, x_eval = \code{X}.)
#' @param degree Degree of local polynomials. (Default: degree = 2.)
#' @param deriv_ord The order of the estimated derivative of the conditional mean
#' outcome function. (Default: deriv_ord = 1.)
#' @param h The bandwidth parameter for the treatment/exposure variable.
#' (Default: h = NULL. Then, the rule-of-thumb bandwidth selector in Eq. (A1)
#' of Yang and Tschernig (1999) is used with additional scaling factors C_h.)
#' @param b The bandwidth vector for the confounding variables. (Default: b = NULL.
#' Then, the rule-of-thumb bandwidth selector in Eq. (A1) of Yang and Tschernig (1999)
#' is used with additional scaling factors C_b.)
#' @param C_h The scaling factor for the rule-of-thumb bandwidth parameter \code{h}.
#' @param C_b The scaling factor for the rule-of-thumb bandwidth vector \code{b}.
#' @param print_bw The indicator of whether the current bandwidth parameters
#' should be printed to the console. (Default: print_bw = TRUE.)
#' @param kernT,kernS The names of kernel functions for the treatment/exposure
#' variable and confounding variables. (Default: kernT = "epanechnikov",
#' kernS = "epanechnikov".)
#'
#' @return The estimated conditional mean outcome function or its partial
#' derivatives evaluated at points \code{x_eval}.
#'
#' @author Yikun Zhang, \email{yikunzhang@@foxmail.com}
#' @references Zhang, Y., Chen, Y.-C., and Giessing, A. (2024)
#' \emph{Nonparametric Inference on Dose-Response Curves Without the Positivity Condition.}
#' \url{https://arxiv.org/abs/2405.09003}.
#'
#' Fan, J. and Gijbels, I. (1996) \emph{Local Polynomial Modelling and its
#' Applications. Chapman & Hall/CRC.}
#' @keywords regression polynomial local (partial)
#'
#' @examples
#' library(parallel)
#' set.seed(123)
#' n <- 300

#' S2 <- cbind(2 * runif(n) - 1, 2 * runif(n) - 1)
#' Z2 <- 4 * S2[, 1] + S2[, 2]
#' E2 <- 0.2 * runif(n) - 0.1
#' T2 <- cos(pi * Z2^3) + Z2 / 4 + E2
#' Y2 <- T2^2 + T2 + 10 * Z2 + rnorm(n, mean = 0, sd = 1)
#' X2 <- cbind(T2, S2)

#' t_qry2 = seq(min(T2) + 0.01, max(T2) - 0.01, length.out = 100)
#' chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
#' if (nzchar(chk) && chk == "TRUE") {
#'   # use 2 cores in CRAN/Travis/AppVeyor
#'   num_workers <- 2L
#' } else {
#'   # use all cores in devtools::test()
#'   num_workers <- parallel::detectCores()
#' }
#' Y_est2 = LocalPolyReg(Y2, X2, x_eval = NULL, degree = 2, deriv_ord = 0,
#'                       h = NULL, b = NULL, C_h = 7, C_b = 3, print_bw = TRUE,
#'                       kernT = "epanechnikov", kernS = "epanechnikov")
#'
#' @export
#'

LocalPolyReg <- function(Y, X, x_eval = NULL, degree = 2, deriv_ord = 1, h = NULL,
                         b = NULL, C_h = 7, C_b = 3, print_bw = TRUE,
                         kernT = "epanechnikov", kernS = "epanechnikov") {
  if (is.null(x_eval)) {
    x_eval <- X
  }

  if (is.null(h) && is.null(b)) {
    bw_params <- RoTBWLocalPoly(Y, X, kernT = kernT, kernS = kernS, C_h = C_h, C_b = C_b)
    h <- bw_params$h
    b <- bw_params$b
  } else if (is.null(h)) {
    bw_params <- RoTBWLocalPoly(Y, X, kernT = kernT, kernS = kernS, C_h = C_h, C_b = C_b)
    h <- bw_params$h
  } else if (is.null(b)) {
    bw_params <- RoTBWLocalPoly(Y, X, kernT = kernT, kernS = kernS, C_h = C_h, C_b = C_b)
    b <- bw_params$b
  }

  if (print_bw) {
    cat("The current bandwidth for treatment variable in the local polynomial regression is", h, ".\n")
    cat("The current bandwidth for confounding variables in the local polynomial regression is", b, ".\n")
  }

  Y_est <- LocalPolyRegMain(Y, X, x_eval = x_eval, degree = degree, deriv_ord = deriv_ord,
                            h = h, b = b, kernT = kernT, kernS = kernS)
  return(Y_est)
}
