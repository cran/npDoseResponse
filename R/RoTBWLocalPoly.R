#' The rule-of-thumb bandwidth selector for the (partial) local polynomial regression.
#'
#' This function implements the rule-of-thumb bandwidth selector for the (partial)
#' local polynomial regression.
#'
#' @param Y The input n-dimensional outcome variable vector.
#' @param X The input n*(d+1) matrix. The first column of X stores the
#' treatment/exposure variables, while the other d columns are confounding variables.
#' @param kernT,kernS The names of kernel functions for the treatment/exposure
#' variable and confounding variables. (Default: kernT = "epanechnikov",
#' kernS = "epanechnikov".)
#' @param C_h,C_b The scaling factors for the rule-of-thumb bandwidth parameters.
#'
#' @return A list that contains two elements.
#' \item{h}{The rule-of-thumb bandwidth parameter for the treatment/exposure variable.}
#' \item{b}{The rule-of-thumb bandwidth vector for the confounding variables.}
#'
#' @author Yikun Zhang, \email{yikunzhang@@foxmail.com}
#' @references Zhang, Y., Chen, Y.-C., and Giessing, A. (2024)
#' \emph{Nonparametric Inference on Dose-Response Curves Without the Positivity Condition.}
#' \url{https://arxiv.org/abs/2405.09003}.
#'
#' Yang, L. and Tschernig, R. (1999). \emph{Multivariate Bandwidth Selection for
#' Local Linear Regression. Journal of the Royal Statistical Society Series B: Statistical Methodology,
#' 61(4), 793-815.}
#' @keywords regression polynomial local for selector bandwidth rule-of-thumb
#'
#' @export
#' @importFrom stats lm
#'

RoTBWLocalPoly <- function(Y, X, kernT = "epanechnikov", kernS = "epanechnikov",
                           C_h = 7, C_b = 3) {
  n <- nrow(X)
  d <- ncol(X) - 1

  kernel_result <- KernelRetrieval(kernT)
  kernT <- kernel_result$KernFunc
  sigmaK_sq <- kernel_result$sigmaK_sq
  K_sq <- kernel_result$K_sq

  p_coeff <- unname(lm(Y ~ poly(X[,1], 4, raw = TRUE))$coefficients)
  sec_deriv <- 12 * p_coeff[5] * X[,1] + 6 * p_coeff[4] * X[,1] + 2 * X[,1]
  C_fun <- mean(sec_deriv^2)

  T1 <- matrix(X[,1], ncol = 1)
  lhs <- as.matrix(cbind(1, X, T1^2, T1^3, T1^4))
  rcond <- .Machine$double.eps * max(dim(lhs))
  beta <- solve(t(lhs) %*% lhs, t(lhs) %*% Y, tol = rcond)

  resid <- sum((Y - lhs %*% beta)^2) * (max(X[,1]) - min(X[,1])) / (n - 5)
  sigmaK_sq <- sigmaK_sq^2

  h <- ((K_sq * resid) / (4 * n * sigmaK_sq * C_fun))^(1/5) * (n^(d / (5 * (d + 5)))) * C_h

  kernel_result <- KernelRetrieval(kernS)
  kernS <- kernel_result$KernFunc
  sigmaK_sq <- kernel_result$sigmaK_sq
  K_sq <- kernel_result$K_sq

  sec_deriv <- matrix(0, nrow = n, ncol = d)
  for (i in 2:(d + 1)) {
    p_coeff <- unname(lm(Y ~ poly(X[,i], 4, raw = TRUE))$coefficients)
    sec_deriv[,i-1] <- 12 * p_coeff[5] * X[,i] + 6 * p_coeff[4] * X[,i] + 2 * X[,i]
  }
  C_fun <- sum(diag(t(sec_deriv) %*% sec_deriv)) / n

  lhs <- as.matrix(cbind(1, X, T1^2, T1^3, T1^4))
  rcond <- .Machine$double.eps
  beta <- solve(t(lhs) %*% lhs, t(lhs) %*% Y, tol = rcond)
  S_mat = as.matrix(X[,2:(d+1)])
  resid <- sum((Y - lhs %*% beta)^2) * (apply(S_mat, 2, max) - apply(S_mat, 2, min)) / (n-5)

  sigmaK_sq <- sigmaK_sq^2
  K_sq <- K_sq^d
  b <- ((K_sq * d * resid) / (4 * n * sigmaK_sq * C_fun))^(1 / (d + 5)) * C_b

  return(list(h = h, b = b))
}
