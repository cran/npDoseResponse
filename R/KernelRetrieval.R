#' The helper function for retrieving a kernel function and its associated statistics.
#'
#' This function helps retrieve the commonly used kernel function, its second moment,
#' and its variance based on the name.
#'
#' @param name The lower-case full name of the kernel function.
#'
#' @return A list that contains three elements.
#' \item{KernFunc}{The interested kernel function.}
#' \item{sigmaK_sq}{The second moment of the kernel function.}
#' \item{K_sq}{The variance of the kernel function.}
#'
#' @author Yikun Zhang, \email{yikunzhang@@foxmail.com}
#' @keywords statistics associated its and function kernel retrieve
#'
#' @examples
#' \donttest{
#'   kernel_result <- KernelRetrieval("epanechnikov")
#'   kernT <- kernel_result$KernFunc
#'   sigmaK_sq <- kernel_result$sigmaK_sq
#'   K_sq <- kernel_result$K_sq
#' }
#'
#' @export
#'

KernelRetrieval <- function(name) {
  if (name == "rectangular") {
    rectangular <- function(t) {
      ind <- (abs(t) <= 1)
      res <- abs(0.5 * ind)
      return(res)
    }
    return(list(KernFunc = rectangular, sigmaK_sq = 1/3, K_sq = 1/2))
  }

  if (name == "triangular") {
    triangular <- function(t) {
      ind <- (abs(t) <= 1)
      res <- abs((1 - abs(t)) * ind)
      return(res)
    }
    return(list(KernFunc = triangular, sigmaK_sq = 1/6, K_sq = 2/3))
  }

  if (name == "epanechnikov") {
    epanechnikov <- function(t) {
      ind <- (abs(t) <= 1)
      res <- abs(0.75 * (1 - t^2) * ind)
      return(res)
    }
    return(list(KernFunc = epanechnikov, sigmaK_sq = 1/5, K_sq = 3/5))
  }

  if (name == "biweight") {
    biweight <- function(t) {
      ind <- (abs(t) <= 1)
      res <- abs(((15/16) * (1 - t^2)^2) * ind)
      return(res)
    }
    return(list(KernFunc = biweight, sigmaK_sq = 1/7, K_sq = 5/7))
  }

  if (name == "triweight") {
    triweight <- function(t) {
      ind <- (abs(t) <= 1)
      res <- abs((35/32) * (1 - t^2)^3 * ind)
      return(res)
    }
    return(list(KernFunc = triweight, sigmaK_sq = 1/9, K_sq = 350/429))
  }

  if (name == "tricube") {
    tricube <- function(t) {
      ind <- (abs(t) <= 1)
      res <- abs((70/81) * (1 - abs(t)^3)^3 * ind)
      return(res)
    }
    return(list(KernFunc = tricube, sigmaK_sq = 35/243, K_sq = 175/247))
  }

  if (name == "gaussian") {
    gaussian <- function(t) {
      res <- (1 / sqrt(2 * pi)) * exp(-0.5 * t^2)
      return(res)
    }
    return(list(KernFunc = gaussian, sigmaK_sq = 1, K_sq = 1/(2*sqrt(pi))))
  }

  if (name == "bigaussian") {
    bigaussian <- function(t) {
      res <- (2 / sqrt(pi)) * (t^2) * exp(-t^2)
      return(res)
    }
    return(list(KernFunc = bigaussian, sigmaK_sq = 3/2, K_sq = 3*sqrt(2/pi)/8))
  }

  if (name == "cosine") {
    cosine <- function(t) {
      ind <- (abs(t) <= 1)
      res <- abs((pi/4) * cos(pi * t / 2) * ind)
      return(res)
    }
    return(list(KernFunc = cosine, sigmaK_sq = 1 - 8/(pi^2), K_sq = pi^2/16))
  }

  if (name == "logistic") {
    logistic <- function(t) {
      res <- 1 / (exp(t) + 2 + exp(-t))
      return(res)
    }
    return(list(KernFunc = logistic, sigmaK_sq = pi^2/3, K_sq = 1/6))
  }

  if (name == "sigmoid") {
    sigmoid <- function(t) {
      res <- (2 / pi) / (exp(t) + exp(-t))
      return(res)
    }
    return(list(KernFunc = sigmoid, sigmaK_sq = pi^2/4, K_sq = 2/(pi^2)))
  }

  if (name == "silverman") {
    silverman <- function(t) {
      res <- 0.5 * exp(-abs(t) / sqrt(2)) * sin(abs(t) / sqrt(2) + pi / 4)
      return(res)
    }
    return(list(KernFunc = silverman, sigmaK_sq = 0, K_sq = 3*sqrt(2)/16))
  }

  return(NULL)
}
