#' The main function of the (partial) local polynomial regression.
#'
#' This function implements the main part of the (partial) local polynomial
#' regression for estimating the conditional mean outcome function and
#' its partial derivatives.
#'
#' @param Y The input n-dimensional outcome variable vector.
#' @param X The input n*(d+1) matrix. The first column of X stores the
#' treatment/exposure variables, while the other d columns are confounding variables.
#' @param x_eval The n*(d+1) matrix for evaluating the local polynomial regression
#' estimates. (Default: x_eval = NULL. Then, x_eval = \code{X}.)
#' @param degree Degree of local polynomials. (Default: degree = 2.)
#' @param deriv_ord The order of the estimated derivative of the conditional mean
#' outcome function. (Default: deriv_ord = 1.)
#' @param h,b The bandwidth parameters for the treatment/exposure variable
#' and confounding variables (Default: h = NULL, b = NULL. Then, the rule-of-thumb
#' bandwidth selector in Eq. (A1) of Yang and Tschernig (1999) is used with
#' additional scaling factors C_h and C_b, respectively.)
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
#' @keywords regression polynomial local (partial) the of function main
#'
#' @export
#'

LocalPolyRegMain <- function(Y, X, x_eval = NULL, degree = 2, deriv_ord = 1,
                             h = NULL, b = NULL, kernT = "epanechnikov",
                             kernS = "epanechnikov") {

  # Retrieve the kernel functions
  kernT_info <- KernelRetrieval(kernT)
  kernT_func <- kernT_info$KernFunc

  kernS_info <- KernelRetrieval(kernS)
  kernS_func <- kernS_info$KernFunc

  n <- nrow(X)  # Number of data points
  d <- ncol(X) - 1
  if (is.null(x_eval)) {
    x_eval <- X
  }
  Y_est <- numeric(nrow(x_eval))

  for (i in 1:nrow(x_eval)) {
    s_cur = as.numeric(x_eval[i,-1])
    weights1 <- kernT_func((X[,1] - x_eval[i,1]) / h) * apply(kernS_func(sweep(sweep(as.matrix(X[,-1]), 2, s_cur, "-"), 2, b, "/")), 1, prod)
    inds <- which(abs(weights1) > 1e-26)
    if (length(inds) == 0) {
      Y_est[i] = 0
      next
    }
    X_dat <- matrix(0, nrow = n, ncol = degree + 1 + d)

    for (p in 0:degree) {
      X_dat[, p + 1] <- (X[,1] - x_eval[i,1])^p
    }
    X_dat[, (degree + 2):(degree + 1 + d)] <- sweep(as.matrix(X[,-1]), 2, s_cur, "-")

    weight_sqrt <- sqrt(weights1[inds])
    if (length(weight_sqrt) == 1) {
      lhs <- weight_sqrt * X_dat[inds,]
      rhs <- weight_sqrt * Y[inds]
      beta <- (rhs * lhs) / sum(lhs^2)
    } else {
      lhs <- diag(weight_sqrt) %*% X_dat[inds,]
      if (det(t(lhs) %*% lhs) < 1e-26) {
        Y_est[i] = 0
        next
      }
      rhs <- weight_sqrt * Y[inds]
      beta <- solve(t(lhs) %*% lhs, t(lhs) %*% rhs, tol = 1e-36)
    }
    Y_est[i] <- factorial(deriv_ord) * beta[deriv_ord + 1]
  }

  return(Y_est)
}
