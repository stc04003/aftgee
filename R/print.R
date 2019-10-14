#' @export
print.aftgee <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("Coefficients:\n")
  print(x$coef.res)
  cat("\n Initial Estimator:\n")
  print(x$coef.init)
}

#' @export
print.aftsrr <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n Coefficients:\n")
  print(coef(x))
}
