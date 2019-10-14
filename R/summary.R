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

