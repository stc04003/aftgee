#' @export
#' @noRd
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
#' @noRd
residuals.aftsrr <- function(object, ...){
  z <- object
  if (class(z) != "aftsrr"){
    stop("Most be aftsrr class")
  }
  ans <- z["call"]
  out <- log(z$y[,1]) - as.matrix(z$x) %*% z$beta
  out
}

#' @export
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
residuals.aftgee <- function(object, ...){
  z <- object
  if (class(z) != "aftgee"){
    stop("Most be aftgee class")
  }
  if ("(Intercept)" %in% attr(object$coef.res, "names")) z$x <- as.matrix(cbind(1, z$x))
  else z$x <- as.matrix(z$x)
  ans <- z["call"]
  out <- log(z$y) - z$x %*% z$coef.res
  out
}

#' @export
#' @noRd
predict.aftsrr <- function(object, newdata = NULL, se.fit = FALSE, type = "lp", ...){
  z <- object
  out <- NULL
  z$x <- as.matrix(z$x)
  if (is.null(newdata)) {
      out$fit <- as.numeric(z$x %*% z$beta)
      if (type == "response") {
          out$fit <- as.numeric(exp(out$fit))
      }
  }
  if (!is.null(newdata)) {
      newdata0 <- as.matrix(newdata, ncol = length(z$beta))
      out$fit <- as.numeric(newdata0 %*% z$beta)
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
          if (is.null(newdata)) {
              var[[i]] <- as.numeric(sqrt(diag(z$x %*% se.srr %*% t(z$x))))
          }
          if (!is.null(newdata)) {
              var[[i]] <- as.numeric(sqrt(diag(newdata0 %*% se.srr %*% t(newdata0))))
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
#' @noRd
predict.aftgee <- function(object, newdata = NULL, se.fit = FALSE, ...){
    z <- object
    out <- NULL
    if ("(Intercept)" %in% attr(object$coef.res, "names")) z$x <- as.matrix(cbind(1, z$x))
    else z$x <- as.matrix(z$x)
    if (class(z) != "aftgee"){
        stop("Most be aftgee class")
    }
    ans <- z["call"]
    if (is.null(newdata)) {
        out$fit <- exp(drop(z$x %*% z$coef.res)) # originally on log scale
        if (se.fit == TRUE) {
            out$se.fit <- sqrt(diag(z$x %*% z$var.res %*% t(z$x)))
        }
    }
    if (!is.null(newdata)) {
        newdata0 <- as.matrix(newdata, ncol = length(z$coef.res))
        if (z$intercept == TRUE & ncol(newdata0) < length(z$coef.res)) {
            newdata0 <- cbind(1, newdata0)
        }
        out$fit <- exp(drop(newdata0 %*% z$coef.res))
        if (se.fit == TRUE) {
            out$se.fit <- sqrt(diag(newdata0 %*% z$var.res %*% t(newdata0)))
        }
    }
    out
}
