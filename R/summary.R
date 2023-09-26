#' @export
summary.aftgee <- function(object,...){
  z <- object
  if (!is.aftgee(z)) stop("Most be aftgee class")
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
  if (!is.aftsrr(z)) stop("Most be aftsrr class")
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

#' @noRd
formatPerc <- function (probs, digits) 
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")

#' @importFrom stats qnorm setNames
#' @export
confint.aftsrr <- function(object, parm, level = 0.95, ...) {
  cf <- coef(object)
  pnames <- names(cf)
  object$covmat <- object$covmat[!is.na(object$covmat)]
  if (!length(object$covmat))
    stop("Missing covariance-variance estimate.")
  vnames <- names(object$covmat)
  ses <- lapply(object$covmat, function(e) {
    ses <- sqrt(diag(e))
    names(ses) <- pnames
    return(ses)
  })
  if (is.matrix(cf)) 
    cf <- setNames(as.vector(cf), pnames)
  if (missing(parm)) 
    parm <- pnames
  else if (is.numeric(parm)) 
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a) 
  fac <- qnorm(a)
  ciList <- lapply(ses, function(e) {
    pct <- formatPerc(a, 3)
    ci <- array(NA_real_, dim = c(length(parm), 2L), dimnames = list(parm, pct))
    ci[] <- cf[parm] + e[parm] %o% fac
    ci
  })
  for (i in 1:length(ciList)) {
    cat("\n")
    cat("Variance Estimator:", as.character(names(ciList)[i]))
    cat("\n")
    print(ciList[[i]])
  }
  invisible(ciList)  
}

#' @export
confint.aftgee <- function(object, parm, level = 0.95, ...) {
  cf <- coef(object)
  pnames <- names(cf)
  ses <- sqrt(diag(vcov(object)))
  names(ses) <- pnames
  if (is.matrix(cf)) 
    cf <- setNames(as.vector(cf), pnames)
  if (missing(parm)) 
    parm <- pnames
  else if (is.numeric(parm)) 
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a) 
  fac <- qnorm(a)
  pct <- formatPerc(a, 3)
  ci <- array(NA_real_, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ci[] <- cf[parm] + ses[parm] %o% fac
  ci
}

is.aftgee <- function(x) inherits(x, "aftgee")
is.aftsrr <- function(x) inherits(x, "aftsrr")
