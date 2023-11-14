#' Least-Squares Approach for Accelerated Failure Time with Generalized Estimating Equation
#'
#' Fits a semiparametric accelerated failure time (AFT) model with least-squares approach.
#' Generalized estimating equation is generalized to multivariate AFT modeling to account
#' for multivariate dependence through working correlation structures to improve efficiency.
#'
#' @param formula  a formula expression, of the form \code{response ~ predictors}.
#'     The \code{response} is a \code{Surv} object with right censoring.
#'     In the case of no censoring, \code{aftgee} will return an ordinary
#'     least estimate when \code{corstr = "independence"}.
#'     See the documentation of \code{lm}, \code{coxph} and \code{formula} for details.
#' @param data an optional data.frame in which to interpret the variables occurring
#'     in the \code{formula}.
#' @param subset an optional vector specifying a subset of observations
#'     to be used in the fitting process.
#' @param id an optional vector used to identify the clusters.
#'     If missing, then each individual row of \code{data} is presumed to
#'     represent a distinct subject.
#'     The length of \code{id} should be the same as the number of observations.
#' @param contrasts an optional list.
#' @param weights an optional vector of observation weights.
#' @param margin a \code{sformula} vector; default at 1.
#' @param corstr a character string specifying the correlation structure.
#'     The following are permitted:
#'     \itemize{
#'     \item \code{independence}
#'     \item \code{exchangeable}
#'     \item \code{ar1}
#'     \item \code{unstructured}
#'     \item \code{userdefined}
#'     \item \code{fixed}
#' }
#' @param binit an optional vector can be either a numeric vector or a character string
#'     specifying the initial slope estimator.
#'     \itemize{
#'     \item When \code{binit} is a vector, its length should be the same as the
#' dimension of covariates.
#'     \item When \code{binit} is a character string, it should be either \code{lm} for simple linear
#' regression, or \code{srrgehan} for smoothed Gehan weight estimator.
#' } The default value is "srrgehan".
#' @param B a numeric value specifies the resampling number.
#'     When B = 0, only the beta estimate will be displayed.
#' @param control controls maxiter and tolerance.
#'
#' @return An object of class "\code{aftgee}" representing the fit.
#' The \code{aftgee} object is a list containing at least the following components:
#' \describe{
#'   \item{coefficients}{a vector of initial value and a vector of point estimates}
#'   \item{coef.res}{a vector of point estimates}
#'   \item{var.res}{estimated covariance matrix}
#'   \item{coef.init}{a vector of initial value}
#'   \item{var.init.mat}{estimated initial covariance matrix}
#'   \item{binit}{a character string specifying the initial estimator.}
#'   \item{conv}{An integer code indicating type of convergence after GEE
#'   iteration. 0 indicates successful convergence; 1 indicates that the
#'   iteration limit \code{maxit} has been reached}
#'   \item{ini.conv}{An integer code indicating type of convergence for
#'   initial value. 0 indicates successful convergence; 1 indicates that the
#'   iteration limit \code{maxit} has been reached}
#'   \item{conv.step}{An integer code indicating the step until convergence}
#' }
#'
#' @references Chiou, S., Kim, J. and Yan, J. (2014) Marginal Semiparametric Multivariate
#' Accelerated Failure Time Model with Generalized Estimating Equation.
#' \emph{Lifetime Data Analysis}, \bold{20}(4): 599--618.
#' @references Jin, Z. and Lin, D. Y. and Ying, Z. (2006)
#' On Least-squares Regression with Censored Data. \emph{Biometrika}, \bold{90}, 341--353.
#'
#' @importFrom geepack geese geese.fit geese.control
#' @importFrom utils capture.output
#' @export
#' @keywords aftgee
#'
#' @example inst/examples/ex_aftgee.R
aftgee <- function(formula, data, subset, id = NULL, contrasts = NULL,
                   weights = NULL, margin = NULL, 
                   corstr = c("independence", "exchangeable", "ar1",
                              "unstructured", "userdefined", "fixed"),
                   binit = "srrgehan", B = 100,
                   control = aftgee.control()
                   ) {
  corstr <- match.arg(corstr)
  scall <- match.call()
  mnames <- c("", "formula", "data", "weights", "margin", "subset", "na.action", "id")
  cnames <- names(scall)
  cnames <- cnames[match(mnames, cnames, 0)]
  mcall <- scall[cnames]
  ##  if (is.null(mcall$id)) mcall$id <- as.name("id")
  mcall[[1]] <- as.name("model.frame")
  m <- eval(mcall, parent.frame())
  id <- model.extract(m, id)
  mterms <- attr(m, "terms")
  weights <- model.extract(m, weights) 
  obj <- unclass(m[,1]) 
  if (!is.Surv(m[[1]]) || ncol(obj) > 2)
    stop("aftgee only supports Surv object with right censoring.", call. = FALSE)
  if (is.null(id)) id <- 1:nrow(obj)
  if (is.null(weights)) weights <- rep(1, nrow(obj))
  margin <- model.extract(m, margin)
  if (is.null(margin)) margin <- rep(1, nrow(obj))
  formula[[2]] <- NULL
  ## Create DF; the first 2 columns are from Surv with time and status
  ## time, status, id, weights, margin, x1, x2, ...
  if (formula == ~1) {
    stop("No covariates are detected.")
    ## DF <- cbind(obj, zero = 0)
  } else {
    DF <- cbind(obj, id, weights, margin, model.matrix(mterms, m, contrasts))
    yint <- (sum(colnames(DF) == "(Intercept)") > 0)
    if (sum(colnames(DF) == "(Intercept)") > 0)
      DF <- DF[,-which(colnames(DF) == "(Intercept)")]
  }
  DF <- as.data.frame(DF)
  if (any(DF$time < 0))
    stop("Failure time must be positive.")
  if (any(DF$time == 0))
    DF$time <- ifelse(DF$time == 0, min(min2(DF$time) / 2, 1e-5), DF$time)
  out <- NULL
  if (sum(DF$status) == nrow(DF)) {
    cat("There is no censoring in the data, ordinary least squares approach is fitted via geese.\n")
    out$geese <- geese(as.formula(paste("log(time) ~ ", formula)[2]),
                       weights = weights, id = id, data = DF, corstr = corstr)
    out$coef.init <- out$coef.res <- as.numeric(out$geese$beta)
    out$coefficients <- cbind(out$coef.init, out$coef.res)
    out$var.res <- out$geese$vbeta
  }
  else {
    out <- aftgee.fit(DF = DF, corstr = corstr, B = B,
                      binit = binit, control = control, yint = yint)
  }
  out$data <- list(y = DF$time, d = DF$status, w = DF$weights, id = DF$id,
                   x = as.matrix(DF[,-(1:5)]))
  rownames(out$coefficients) <- names(out$coef.res) <-
    names(out$coef.init) <- colnames(model.matrix(mterms, m, contrasts))
  ## out$intercept <- (sum(x[,1]) == nrow(x))
  colnames(out$coefficients) <- c("binit", "AFTGEE")
  out$call <- scall
  class(out) <- "aftgee"
  out
}

aftgee.fit <- function(DF, corstr="independence",
                       B = 100, binit = "lm", yint = TRUE,
                       control = aftgee.control()) {
  x <- as.matrix(DF[,-(1:5), drop = FALSE])
  id <- DF$id
  n <- length(unique(id))
  rm <- NULL
  rmName <- NULL
  firstBeta <- firstSd <- firstSdMat <- firstconvergence <- NULL
  clsize <- unlist(lapply(split(id, id), length))
  N <- sum(clsize)
  if (is.numeric(binit)) {
    if (length(binit) != ncol(x) + yint * 1)
      stop("binit value length does not match with numbers of covariates.", call. = FALSE)
    firstBeta <- binit      
  }
  if (!(is.numeric(binit))) {
    if (!(binit %in% c("lm", "srrgehan"))) 
      stop("Invalid binit value method", call. = FALSE)
  }
  if (!(is.numeric(binit))) {
    if (binit == "lm") {
      if (yint) linfit <- summary(lm(log(DF$time) ~ x, subset = DF$status > 0))
      else linfit <- summary(lm(log(DF$time) ~ x - 1, subset = DF$status > 0))
      first <- list(beta = linfit$coef[,1], sd = linfit$coef[,2])
      firstBeta <- first$beta
      firstSd <- first$sd
      firstconvergence <- first$convergence
    }
    if (binit == "srrgehan") {
      engine.control <- control[names(control) %in% names(attr(getClass("gehan.is"), "slots"))]
      engine <- do.call("new", c(list(Class = "gehan.is"), engine.control))
      if (engine@b0 == 0) {
        ## engine@b0 <- as.numeric(coef(lm(log(DF$time) ~ x - 1, subset = DF$status > 0)))
        engine@b0 <- as.numeric(coef(lm(log(DF$time) ~ x))[-1])
      }
      engine@sigma0 <- diag(length(engine@b0))
      first <- rankFit.gehan.is(DF[,-5], engine, NULL)
      firstBeta <- first$beta
      firstSdMat <- NA
      firstconvergence <- first$conv
      if (yint) firstBeta <- c(mean(log(DF$time) - x %*% firstBeta), firstBeta) 
    }
  }
  if (yint) x <- cbind(1, x)
  binitValue <- list(beta = firstBeta, sd = firstSd, sdMat = firstSdMat)
  result <- aftgee.est(log(DF$time), x, DF$status, binitValue$beta, id, corstr,
                       rep(1, nrow(DF)), DF$margin, DF$weights, control)
  ## variance estimation
  ## Force trace to be false
  control$trace <- FALSE
  zsamp <- bsamp <- NULL
  if (B > 0) {
    if (!control$parallel) {
      bsamp <- matrix(NA, nrow = B, ncol = length(result$beta))
      bini <- result$beta
      if (control$seIni) {
        for (i in 1:B) {
          Z <- as.vector(rep(rexp(n, 1), time = clsize))
          zsamp[[i]] <- Z
          DF0 <- DF
          DF0$weights <- Z * DF0$weights
          bini <- rankFit.gehan.is(DF0[,-5], engine, NULL)$beta
          if (yint) bini <- c(mean(log(DF0$time) - as.matrix(DF0[,-(1:5)]) %*% bini), bini)
          bsamp[i,] <- aftgee.est(log(DF$time), x, DF$status, bini, id, corstr, Z,
                                  DF$margin, DF$weights, control)$beta
        }
      } else {
        bsamp <- resampling_No_Margin(y = log(DF$time), X = x, D = DF$status, b0 = bini,
                                      ## b0 = binitValue$beta,
                                      nt = tabulate(id), w = DF$weights, corstr = corstr, B = B,
                                      tol = control$reltol, maxit = control$maxit)
      }
      ## resampling is less stable when sample size is small
      bsamp <- bsamp[rowSums(bsamp * bsamp) < 5 * sum(bini * bini),]      
      vhat <- var(bsamp)
    } else {
      cl <- makeCluster(control$parCl)
      clusterExport(cl = cl,
                    varlist = c("n", "clsize", "DF", "control", "x", "id", "corstr", "result"),
                    envir = environment())
      bsamp <- parSapply(cl, 1:B, function(z) {
        Z <- as.vector(rep(rexp(n, 1), time = clsize))
        if (control$seIni) {
          DF0 <- DF
          DF0$weights <- Z * DF0$weights
          bini <- rankFit.gehan.is(DF0[,-5], engine, NULL)$beta
          if (yint) bini <- c(mean(log(DF0$time) - as.matrix(DF0[,-(1:5)]) %*% bini), bini)
        } else {
          bini <- result$beta
        }
        aftgee.est(log(DF$time), x, DF$status, bini, id, corstr, Z,
                   DF$margin, DF$weights, control)$beta
      })
      stopCluster(cl)
      vhat <- var(t(bsamp))
    } ## end parallel 
  }
  if (B == 0) vhat <- NULL
  ini.beta <- c(binitValue$beta)
  ini.sd <- c(binitValue$sd)
  ini.sdMat <- c(binitValue$sdMat)
  fit <- list(coefficients = cbind(ini.beta, result$beta),
              coef.res = as.numeric(result$beta),
              var.res = vhat,
              varMargin = result$gamma,
              coef.trace = result$histBeta,
              bj.trace = result$histYhat,
              alpha = result$alpha,
              coef.init = ini.beta,
              sd.init = ini.sd,
              var.init.mat = ini.sdMat,
              binit = binit,
              conv = result$convergence,
              ini.conv = firstconvergence,
              beta.resmp = bsamp,
              zsamp = zsamp,
              conv.step = result$convStep)
  class(fit) <- "aftgee.fit"
  fit
}

## aftgee.se; bootstrap or resampling, true value or aftsrr as initial value

#' Auxiliary for Controlling AFTGEE Fitting
#'
#' Auxiliary function as user interface for \code{aftgee} and \code{aftsrr} fitting.
#'
#' When \code{trace} is TRUE, output for each iteration is printed to the screen.
#' 
#' @param maxiter max number of iteration.
#' @param reltol relative error tolerance.
#' @param trace a binary variable, determine whether to display output for each iteration.
#' @param seIni a logical value indicating whether a new rank-based initial value is computed
#' for each resampling sample in variance estimation.
#' @param parallel an logical value indicating whether parallel computing is used for resampling and bootstrap.
#' @param parCl an integer value indicating the number of CPU cores used when \code{parallel = TRUE}.
#' @param gp.pwr an numerical value indicating the GP parameter when \code{rankWeights = GP}.
#' The default value is half the CPU cores on the current host.
#'
#' @export
#' @return A list with the arguments as components.
#' @seealso \code{\link{aftgee}}
aftgee.control <- function(maxiter = 50, reltol = 0.001, trace = FALSE,
                           seIni = FALSE, parallel = FALSE,
                           parCl = parallel::detectCores() / 2, gp.pwr = -999) {
  list(maxiter = maxiter, reltol = reltol, trace = trace, seIni = seIni,
       parallel = parallel, parCl = parCl, gp.pwr = gp.pwr)
}


#' This is a wraper for aftgee.est(); this check for margin and prep data
#' @noRd
aftgee.est <- function(y, x, delta, beta, id,
                       corstr = "independence",
                       Z = rep(1, length(y)),
                       margin = rep(1, length(id)),
                       weights = rep(1, length(y)),
                       control = aftgee.control()) {
  if (length(unique(margin)) == 1L) {
    if (!(corstr %in% c("independence", "exchangeable", "ar1")))
      stop("Invalid corstr when no margins are specified.")
    fit <- est_No_Margin(y = y, X = x, D = delta, b0 = beta, nt = tabulate(id),
                         w = Z * weights, corstr = corstr,
                         tol = control$reltol, maxit = control$maxiter)
    fit$hist_b <- lapply(fit$hist_b[!sapply(fit$hist_b, is.null)], drop)    
    return(list(beta = fit$b, alpha = fit$alpha, histBeta = fit$hist_b, iniBeta = beta,
                convergence = 1 * (fit$iter < control$maxiter),
                convStep = fit$iter))
  }
  est.margin(y, x, delta, beta, id, corstr, Z, margin, weights, control)
}    

#' This is now only called when margina = 1
#' @noRd
est.margin <- function(y, x, delta, beta, id,
                       corstr = "independence",
                       Z = rep(1, length(y)),
                       margin = rep(1, length(id)),
                       weights = rep(1, length(y)),
                       control = aftgee.control()) {
  iniBeta <- beta
  xmat <- as.matrix(x) 
  nobs <- length(y)
  xmatZ <- sqrt(Z * weights) * xmat
  ## If x does not have an intercept term, e.g., a column of 1s, then center xmat and xmatZ
  if (!any(apply(xmat, 2, function(e) length(unique(e))) == 1)) {
    xmat <- (diag(nrow(x)) - 1 / nrow(x)) %*% xmat
    xmatZ <- (diag(nrow(x)) - 1 / nrow(x)) %*% xmatZ
  }
  histBeta <- txts <- NULL
  for (i in 1:control$maxiter) {
    betaprev <- beta
    eres <- eres2 <- er1 <- er2 <- NULL
    e <- y - xmat %*% beta
    for (m in unique(margin)) {
      temp <- eResC2(e[margin == m], delta[margin == m], Z[margin == m])
      temp[[2]] <- ifelse(delta[margin == m] == 1, e[margin == m]^2, temp[[2]])
      eres2[m] <- mean(temp[[2]], na.rm = TRUE)
      dum <- cumsum(ifelse(margin == m, 1, 0))
      er1[[m]] <- temp[[1]][ifelse(margin == m, dum, NA)]
    }
    er1 <- c(do.call(rbind, er1))
    er1 <- er1[!is.na(er1)]
    ## sapply(1:length(y), function(i) er1[[margin[i]]][i])
    yhat <- delta * y + (1 - delta) * (er1 + xmat %*% beta)
    yhatZ <- sqrt(Z * weights) * yhat
    er2 <- as.matrix(eres2[margin])
    if (control$trace)
      txts <- capture.output(
        geefit <- geese.fit(xmatZ, yhatZ, id, zsca = er2, scale.fix = TRUE,
                            corstr = corstr, control = geese.control(trace = TRUE)))
    else         
      geefit <- geese.fit(xmatZ, yhatZ, id, zsca = er2, scale.fix = TRUE, corstr = corstr)
    beta <- geefit$beta
    if (control$trace) {
      cat("\n beta:", as.numeric(beta), "\n")
      txts <- txts[grepl("beta =", txts)]
      txts <- sub("beta = ", "", txts)
      tmp <- lapply(txts, function(x) as.numeric(unlist(strsplit(x, " "))))
      tmp$last <- as.numeric(beta)
      names(tmp) <- paste0("Iter.", 1:length(tmp))
      histBeta[[i]] <- tmp
    }
    convStep <- i
    if (max(abs(beta - betaprev) / abs(beta)) <= control$reltol) break
  } 
  ## Fitting independence structure at convergence for QIC calculation
  ## Move these to QIC calls
  ## if (length(unique(margin)) == 1L) geefit0 <- geese.fit(xmatZ, yhatZ, id, corstr = "independence")
  ## else geefit0 <- geese.fit(xmatZ, yhatZ, id, zsca = er2, scale.fix = TRUE, corstr = "independence")
  beta <- geefit$beta
  alpha <- geefit$alpha
  gamma <- geefit$gamma ## eres2
  convergence <- ifelse(i == control$maxiter, 1, 0)
  if (control$trace)
    names(histBeta) <- paste0("BJ", 1:length(histBeta))
  out <- list(beta = beta, alpha = alpha, ## gamma = gamma,
              histBeta = histBeta, iniBeta = iniBeta,
              ## iteration info.
              convergence = convergence, convStep = convStep)
  return(out)
}

min2 <- function(x) min(x[x != min(x)])
