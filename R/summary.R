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

#' @export
summary.aftgee <- function(object,...){
  z <- object
  if (class(z) != "aftgee"){
    stop("Most be aftgee class")
  }
  ans <- z["call"]
  TAB.ini <- NULL
  ## aftgee part
  est.gee <- z$coef.res
  if (is.null(z$var.res)) {se.gee <- rep(NaN, length(est.gee))
  } else {se.gee <- sqrt(diag(z$var.res))}
  est.temp.gee <- ifelse(se.gee == "NaN", "NaN", est.gee)
  z.val.gee <- as.numeric(est.temp.gee)/as.numeric(se.gee)
  TAB <- cbind(Estimate = round(est.gee, 3),
               StdErr = round(se.gee, 3),
               z.value = round(z.val.gee, 3),
               p.value = round(2 * pnorm(-abs(z.val.gee)), 3))
  rownames(TAB) <- names(z$coef.res)
  ## binit part
  est.ini <- z$coef.init
  res <- list(call = object$call, coefficients=TAB, binit = z$binit,
              iniEst = z$iniEst, est.ini = z$coef.init)
  class(res) <- "summary.aftgee"
  res
}

#' @export
summary.aftsrr <- function(object,...){
  z <- object
  if (class(z) != "aftsrr"){
    stop("Most be aftsrr class")
  }
  ans <- z["call"]
  var.meth <- z$var.meth[z$var.meth %in%
                         c("NULL", "bootstrap", "MB", "ZLCF", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB", "js")]
  se.count <- length(var.meth)
  se.name <- match(var.meth, names(z$covmat))
  est.srr <- z$beta
  p <- length(z$beta)
  TAB.srr <- NULL
  for (i in 1:se.count) {
      se.srr <- NA
      if (z$B != 0 & z$var.meth[1] != "NULL") {
          se.srr <- sqrt(diag(z$covmat[[se.name[i]]]))
      }
      z.val.srr <- as.numeric(est.srr)/as.numeric(se.srr)
      temp.srr <- cbind(Estimate = round(est.srr, 3), StdErr = round(se.srr, 3),
                        z.value = round(z.val.srr, 3), p.value = round(2 * pnorm(-abs(z.val.srr)), 3))
      rownames(temp.srr) <- z$vari.name
      TAB.srr <- append(TAB.srr, list(temp.srr))
  }
  res <- list(call = object$call, coefficients = TAB.srr, var.name = var.meth)
  class(res) <- "summary.aftsrr"
  res
}

#' @export
print.summary.aftgee <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("AFTGEE Estimator")
  cat("\n")
  printCoefmat(as.matrix(x$coefficients), P.values = TRUE, has.Pvalue = TRUE)
}

#' @export
print.summary.aftsrr <- function(x, ...){
  se.count <- length(x$var.name)
  cat("Call:\n")
  print(x$call)
  for (i in 1:se.count){
      cat("\n")
      cat("Variance Estimator:", as.character(x$var.name[i]))
      cat("\n")
      printCoefmat(as.data.frame(x$coefficients[i]), P.values = TRUE, has.Pvalue = TRUE)
  }
}

#' @export
coef.aftsrr <- function(object, ...){
  z <- object
  if (class(z) != "aftsrr"){
    stop("Most be aftsrr class")
  }
  ans <- z["call"]
  out <- z$beta
  names(out) <- z$vari.name
  out
}

#' @export
residuals.aftsrr <- function(object, ...){
  z <- object
  if (class(z) != "aftsrr"){
    stop("Most be aftsrr class")
  }
  ans <- z["call"]
  out <- log(z$y[,1]) - z$x %*% z$beta
  out
}

#' @export
vcov.aftsrr <- function(object, ...){
  z <- object
  if (class(z) != "aftsrr"){
    stop("Most be aftsrr class")
  }
  ans <- z["call"]
  var.meth <- z$var.meth[z$var.meth %in% c("MB", "ZLCF", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB", "js")]
  se.count <- length(var.meth)
  se.name <- match(var.meth, names(z$covmat))
  p <- length(z$beta)
  TAB.srr <- NULL
  out <- list(NULL)
  out[se.count + 1] <- NULL
  names(out) <- z$var.meth
  for (i in 1:se.count) {
      se.srr <- z$covmat[[se.name[i]]]
      rownames(se.srr) <- z$vari.name
      colnames(se.srr) <- z$vari.name
      out[[i]] <- se.srr
  }
  out
}

#' @export
coef.aftgee <- function(object, ...){
  z <- object
  if (class(z) != "aftgee"){
    stop("Most be aftgee class")
  }
  ans <- z["call"]
  out <- z$coef.res
  out
}

#' @export
vcov.aftgee <- function(object, ...){
  z <- object
  if (class(z) != "aftgee"){
    stop("Most be aftgee class")
  }
  ans <- z["call"]
  out <- z$var.res
  out
}

#' @export
residuals.aftgee <- function(object, ...){
  z <- object
  if (class(z) != "aftgee"){
    stop("Most be aftgee class")
  }
  ans <- z["call"]
  out <- log(z$y) - z$x %*% z$coef.res
  out
}

#' @export
predict.aftsrr <- function(object, newdata = NULL, se.fit = FALSE, type = "lp", ...){
  z <- object
  out <- NULL
  if (is.null(newdata)) {
      out$fit <- as.numeric(z$x %*% z$beta)
      if (type == "response") {
          out$fit <- as.numeric(exp(out$fit))
      }
  }
  if (!is.null(newdata)) {
      n <- as.matrix(newdata, ncol = length(z$beta))
      out$fit <- as.numeric(n %*% z$beta)
      if (type == "response") {
          out$fit <- as.numeric(exp(out$fit))
      }
  }
  if (se.fit == TRUE) {
      var.meth <- z$var.meth[z$var.meth %in% c("MB", "ZLCF", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB", "js")]
      se.count <- length(var.meth)
      se.name <- match(var.meth, names(z$covmat))
      p <- length(z$beta)
      TAB.srr <- NULL
      var <- list(NULL)
      var[se.count + 1] <- NULL
          names(var) <- z$var.meth
      for (i in 1:se.count) {
          se.srr <- z$covmat[[se.name[i]]]
          ## rownames(se.srr) <- z$vari.name
          ## colnames(se.srr) <- z$vari.name
          if (is.null(newdata)) {
              var[[i]] <- as.numeric(sqrt(diag(z$x %*% se.srr %*% t(z$x))))
          }
          if (!is.null(newdata)) {
              var[[i]] <- as.numeric(sqrt(diag(n %*% se.srr %*% t(n))))
          }
      }
      out$se.fit <- var
      if (type == "response") {
          out$se.fit <- lapply(out$se.fit, function(x) out$fit * x)
      }
  }
  out
}

#' @export
predict.aftgee <- function(object, newdata = NULL, se.fit = FALSE, ...){
    z <- object
    out <- NULL
    if (class(z) != "aftgee"){
        stop("Most be aftgee class")
    }
    ans <- z["call"]
    if (is.null(newdata)) {
        out$fit <- z$x %*% z$coef.res
    }
    if (!is.null(newdata)) {
        n <- as.matrix(newdata, ncol = length(z$coef.res))
        if (z$intercept == TRUE & ncol(n) < length(z$coef.res)) {
            n <- cbind(1, n)
        }
        out$fit <- n %*% z$coef.res
        if (se.fit == TRUE) {
            out$se.fit <- sqrt(diag(n %*% z$var.res %*% t(n)))
        }
    }
    out
}
