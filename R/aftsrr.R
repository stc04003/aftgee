##############################################################################
## Function calls for different point estimators
##############################################################################
## is: with induced smoothing; X - \Phi() / \Phi()
## ns: solve non-smoothing equation directly; X - \phi() / \phi()
## mns: solve non-smoothing equation approximatly; (1 / \Phi()) * non-smooth Gehan
## mis: solve non-smoothing equation approximatly; (1 / \Phi()) * Smooth Gehan

rankFit.gehan.is <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    if (is.null(gw)) gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    gehan.obj <- function(b) {
        .C("gehan_s_obj", as.double(b), as.double(Y), as.double(X), as.double(delta),
           as.integer(clsize), as.double(engine@sigma0), as.integer(length(clsize)),
           as.integer(p), as.integer(n), as.double(W), as.double(gw),
           double(1), PACKAGE = "aftgee")[[12]]
    }
    gehan.est <- function(b) {
        .C("gehan_s_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
           as.integer(clsize), as.double(engine@sigma0), as.integer(length(clsize)),
           as.integer(p), as.integer(n), as.double(W), as.double(gw),
           double(p), PACKAGE = "aftgee")[[12]]
    }
    if (engine@solver %in% c("BBsolve", "dfsane")) {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- tryCatch(
                do.call(engine@solver,list(par = engine@b0, fn = gehan.est,
                                           quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))),
                error = function(e)
                    do.call(engine@solver,list(par = double(p), fn = gehan.est,
                                               quiet = TRUE, control = list(tol = engine@tol, trace = FALSE)))))
        end.time <- Sys.time()
    }
    if (engine@solver == "BBoptim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver, list(par = engine@b0, fn = gehan.obj,
                                               quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    if (engine@solver == "optim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver,
                           list(par = engine@b0, fn = gehan.obj,
                                control = list(abstol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    list(beta = fit$par, conv = fit$convergence, pe.time = end.time - start.time,
         value = c(fit$residual, fit$value), iter = 1)
}

rankFit.gehan.ns <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    if (is.null(gw)) gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    gehan.obj <- function(b) {
        ans <- .C("gehan_ns_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
                  as.integer(clsize), as.integer(length(clsize)), as.integer(p),
                  as.integer(n), as.double(W), as.double(gw),
                  double(p), PACKAGE = "aftgee")[[11]]
        return(sum(ans^2))
    }
    gehan.est <- function(b) {
        .C("gehan_ns_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
           as.integer(clsize), as.integer(length(clsize)), as.integer(p),
           as.integer(n), as.double(W), as.double(gw),
           double(p), PACKAGE = "aftgee")[[11]]
    }
    if (engine@solver %in% c("BBsolve", "dfsane")) {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- tryCatch(do.call(engine@solver, list(par = engine@b0, fn = gehan.est, quiet = TRUE,
                                                        control = list(tol = engine@tol, trace = FALSE))),
                            error = function(e)
                                do.call(engine@solver, list(par = double(p), fn = gehan.est, quiet = TRUE,
                                                            control = list(tol = engine@tol, trace = FALSE)))))
        end.time <- Sys.time()
    }
    if (engine@solver == "BBoptim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver, list(par = engine@b0, fn = gehan.obj, quiet = TRUE,
                                               control = list(tol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    if (engine@solver == "optim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver,
                           list(par = engine@b0, fn = gehan.obj, control = list(trace = FALSE))))
        end.time <- Sys.time()
    }
    list(beta = fit$par, conv = fit$convergence, pe.time = end.time - start.time,
         value = c(fit$residual, fit$value), iter = 1)
}

rankFit.logrank.is <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    if (is.null(gw)) gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    log.obj <- function(b) {
        ans <- .C("log_s_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
                  as.integer(clsize), as.double(engine@sigma0), as.integer(length(clsize)),
                  as.integer(p), as.integer(n), as.double(W), as.double(gw),
                  double(p), PACKAGE = "aftgee")[[12]]
        return(sum(ans^2))
    }
    log.est <- function(b) {
        .C("log_s_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
           as.integer(clsize), as.double(engine@sigma0), as.integer(length(clsize)),
           as.integer(p), as.integer(n), as.double(W), as.double(gw),
           double(p), PACKAGE = "aftgee")[[12]]
    }
    if (engine@solver %in% c("BBsolve", "dfsane")) {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- tryCatch(
                do.call(engine@solver, list(par = engine@b0, fn = log.est,
                                            quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))),
                error = function(e)
                    do.call(engine@solver, list(par = double(p), fn = log.est,
                                                quiet = TRUE, control = list(tol = engine@tol, trace = FALSE)))))
        end.time <- Sys.time()
    }
    if (engine@solver == "BBoptim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver, list(par = engine@b0, fn = log.obj,
                                               quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    if (engine@solver == "optim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver,
                           list(par = engine@b0, fn = log.obj,
                                control = list(abstol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    list(beta = fit$par, conv = fit$convergence, pe.time = end.time - start.time,
         value = c(fit$residual, fit$value), iter = 1)
}

rankFit.logrank.ns <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    if (is.null(gw)) gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    log.obj <- function(b) {
        ans <- .C("log_ns_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
                  as.integer(clsize), as.integer(length(clsize)),
                  as.integer(p), as.integer(n), as.double(W), as.double(gw),
                  double(p), PACKAGE = "aftgee")[[11]]
        return(sum(ans^2))
    }
    log.est <- function(b) {
        .C("log_ns_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
           as.integer(clsize), as.integer(length(clsize)),
           as.integer(p), as.integer(n), as.double(W), as.double(gw),
           double(p), PACKAGE = "aftgee")[[11]]
    }
    if (engine@solver %in% c("BBsolve", "dfsane")) {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- tryCatch(
                do.call(engine@solver, list(par = engine@b0, fn = log.est,
                                            quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))),
                error = function(e)
                    do.call(engine@solver, list(par = double(p), fn = log.est,
                                                quiet = TRUE, control = list(tol = engine@tol, trace = FALSE)))))
        end.time <- Sys.time()
    }
    if (engine@solver == "BBoptim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver, list(par = engine@b0, fn = log.obj,
                                               quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    if (engine@solver == "optim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver,
                           list(par = engine@b0, fn = log.obj,
                                control = list(abstol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    list(beta = fit$par, conv = fit$convergence, pe.time = end.time - start.time,
         value = c(fit$residual, fit$value), iter = 1)
}

rankFit.logrank.mns <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    b1 <- rankFit.gehan.ns(DF, engine, stdErr)$beta
    iter <- 1
    start.time <- Sys.time()
    for (i in 1:engine@maxIter) {
        gw <- 1 / .C("gehan_ns_wt", as.double(b1), as.double(Y), as.double(X), as.integer(clsize),
                     as.integer(length(clsize)), as.integer(p), as.integer(n), as.double(W),
                     double(n), PACKAGE = "aftgee")[[9]]
        gw <- ifelse(gw == Inf, 0, gw)
        b2 <- rankFit.gehan.ns(DF, engine, stdErr, gw)$beta
        engine@b0 <- b2
        iter <- iter + 1
        if (engine@trace) print(b2)
        if (max(abs(b1 - b2)) < engine@tol) {
            conv <- 0
            break
        }
        b1 <- b2
        conv <- 1
    }
    end.time <- Sys.time()
    list(beta = b2, conv = conv, pe.time = end.time - start.time,
         value = abs(b1 - b2), iter = iter)
}

rankFit.logrank.mis <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    b1 <- rankFit.gehan.is(DF, engine, stdErr)$beta
    iter <- 1
    start.time <- Sys.time()
    for (i in 1:engine@maxIter) {
        gw <- 1 / .C("gehan_s_wt", as.double(b1), as.double(Y), as.double(X), as.integer(clsize),
                     as.double(engine@sigma0), as.integer(length(clsize)), 
                     as.integer(p), as.integer(n), as.double(W),
                     double(n), PACKAGE = "aftgee")[[10]]
        gw <- ifelse(gw == Inf, 0, gw)
        b2 <- rankFit.gehan.is(DF, engine, stdErr, gw)$beta
        engine@b0 <- b2
        iter <- iter + 1
        if (engine@trace) print(b2)
        if (max(abs(b1 - b2)) < engine@tol) {
            conv <- 0
            break
        }
        b1 <- b2
        conv <- 1
    }
    end.time <- Sys.time()
    list(beta = b2, conv = conv, pe.time = end.time - start.time,
         value = abs(b1 - b2), iter = iter)
}

rankFit.pw.is <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    if (is.null(gw)) gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    pw.obj <- function(b) {
        er <- Y - X %*% b
        tmp <- survfit(Surv(er, delta) ~ 1, weights = W)
        gw <- approx(tmp$time, tmp$surv, er, "constant", yleft = 1, yright = min(tmp$surv))$y
        ans <- .C("log_s_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
                  as.integer(clsize), as.double(engine@sigma0), as.integer(length(clsize)),
                  as.integer(p), as.integer(n), as.double(W), as.double(gw),
                  double(p), PACKAGE = "aftgee")[[12]]
        return(sum(ans^2))
    }
    pw.est <- function(b) {
        er <- Y - X %*% b
        tmp <- survfit(Surv(er, delta) ~ 1, weights = W)
        gw <- approx(tmp$time, tmp$surv, er, "constant", yleft = 1, yright = min(tmp$surv))$y
        .C("log_s_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
           as.integer(clsize), as.double(engine@sigma0), as.integer(length(clsize)),
           as.integer(p), as.integer(n), as.double(W), as.double(gw),
           double(p), PACKAGE = "aftgee")[[12]]
    }
    if (engine@solver %in% c("BBsolve", "dfsane")) {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- tryCatch(
                do.call(engine@solver, list(par = engine@b0, fn = pw.est,
                                            quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))),
                error = function(e)
                do.call(engine@solver, list(par = double(p), fn = pw.est,
                                            quiet = TRUE, control = list(tol = engine@tol, trace = FALSE)))))
        end.time <- Sys.time()
    }
    if (engine@solver == "BBoptim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver, list(par = engine@b0, fn = pw.obj,
                                               quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    if (engine@solver == "optim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver,
                           list(par = engine@b0, fn = pw.obj, control = list(abstol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    list(beta = fit$par, conv = fit$convergence, pe.time = end.time - start.time,
         value = c(fit$residual, fit$value), iter = 1)
}

rankFit.pw.ns <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    if (is.null(gw)) gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    pw.obj <- function(b) {
        er <- Y - X %*% b
        tmp <- survfit(Surv(er, delta) ~ 1, weights = W)
        gw <- approx(tmp$time, tmp$surv, er, "constant", yleft = 1, yright = min(tmp$surv))$y
        ans <- .C("log_ns_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
                  as.integer(clsize), as.integer(length(clsize)),
                  as.integer(p), as.integer(n), as.double(W), as.double(gw),
                  double(p), PACKAGE = "aftgee")[[11]]
        return(sum(ans^2))
    }
    pw.est <- function(b) {
        er <- Y - X %*% b
        tmp <- survfit(Surv(er, delta) ~ 1, weights = W)
        gw <- approx(tmp$time, tmp$surv, er, "constant", yleft = 1, yright = min(tmp$surv))$y
        .C("log_ns_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
           as.integer(clsize), as.integer(length(clsize)),
           as.integer(p), as.integer(n), as.double(W), as.double(gw),
           double(p), PACKAGE = "aftgee")[[11]]
    }
      if (engine@solver %in% c("BBsolve", "dfsane")) {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- tryCatch(
                do.call(engine@solver, list(par = engine@b0, fn = pw.est,
                                            quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))),
                error = function(e)
                    do.call(engine@solver, list(par = double(p), fn = pw.est,
                                                quiet = TRUE, control = list(tol = engine@tol, trace = FALSE)))))
                end.time <- Sys.time()
    }
    if (engine@solver == "BBoptim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver, list(par = engine@b0, fn = pw.obj,
                                               quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    if (engine@solver == "optim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver,
                           list(par = engine@b0, fn = pw.obj, control = list(abstol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    list(beta = fit$par, conv = fit$convergence, pe.time = end.time - start.time,
         value = c(fit$residual, fit$value), iter = 1)
}

rankFit.pw.mns <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    b1 <- rankFit.gehan.ns(DF, engine, stdErr)$beta
    iter <- 1
    start.time <- Sys.time()
    for (i in 1:engine@maxIter) {
        er <- Y - X %*% b1
        tmp <- survfit(Surv(er, delta) ~ 1, weights = W)
        s0 <- approx(tmp$time, tmp$surv, er, "constant", yleft = 1, yright = min(tmp$surv))$y
        gw <- s0 / .C("gehan_ns_wt", as.double(b1), as.double(Y), as.double(X), as.integer(clsize),
                      as.integer(length(clsize)), as.integer(p), as.integer(n), as.double(W),
                      double(n), PACKAGE = "aftgee")[[9]]
        gw <- ifelse(gw == Inf, 0, gw)
        b2 <- rankFit.gehan.ns(DF, engine, stdErr, gw)$beta
        iter <- iter + 1
        if (engine@trace) print(b2)
        if (max(abs(b1 - b2)) < engine@tol) {
            conv <- 0
            break
        }
        b1 <- b2
        conv <- 1
    }
    end.time <- Sys.time()
    list(beta = b2, conv = conv, pe.time = end.time - start.time,
         value = abs(b1 - b2), iter = iter)
}

rankFit.pw.mis <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    b1 <- rankFit.gehan.is(DF, engine, stdErr)$beta
    iter <- 1
    start.time <- Sys.time()
    for (i in 1:engine@maxIter) {
        er <- Y - X %*% b1
        tmp <- survfit(Surv(er, delta) ~ 1, weights = W)
        s0 <- approx(tmp$time, tmp$surv, er, "constant", yleft = 1, yright = min(tmp$surv))$y
        gw <- s0 / .C("gehan_s_wt", as.double(b1), as.double(Y), as.double(X), as.integer(clsize),
                      as.double(engine@sigma0), as.integer(length(clsize)),
                      as.integer(p), as.integer(n), as.double(W),
                      double(n), PACKAGE = "aftgee")[[10]]
        gw <- ifelse(gw == Inf, 0, gw)
        b2 <- rankFit.gehan.is(DF, engine, stdErr, gw)$beta
        iter <- iter + 1
        if (engine@trace) print(b2)
        if (max(abs(b1 - b2)) < engine@tol) {
            conv <- 0
            break
        }
        b1 <- b2
        conv <- 1
    }
    end.time <- Sys.time()
    list(beta = b2, conv = conv, pe.time = end.time - start.time,
         value = abs(b1 - b2), iter = iter)
}

rankFit.gp.is <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    if (is.null(gw)) gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    pwr <- ifelse(is.null(engine@gp.pwr), 1 / ncol(X), engine@gp.pwr)
    if (pwr < 0) stop("Invalid GP class power.", call. = FALSE)
    gp.obj <- function(b) {
        er <- Y - X %*% b
        tmp <- survfit(Surv(er, delta) ~ 1, weights = W)
        gw <- approx(tmp$time, tmp$surv, er, "constant", yleft = 1, yright = min(tmp$surv))$y
        gw <- gw^pwr
        gw <- ifelse(gw == Inf, 0, gw)
        ans <- .C("log_s_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
                  as.integer(clsize), as.double(engine@sigma0), as.integer(length(clsize)),
                  as.integer(p), as.integer(n), as.double(W), as.double(gw),
                  double(p), PACKAGE = "aftgee")[[12]]
        return(sum(ans^2))
    }
    gp.est <- function(b) {
        er <- Y - X %*% b
        tmp <- survfit(Surv(er, delta) ~ 1, weights = W)
        gw <- approx(tmp$time, tmp$surv, er, "constant", yleft = 1, yright = min(tmp$surv))$y
        gw <- gw^pwr
        gw <- ifelse(gw == Inf, 0, gw)
        .C("log_s_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
           as.integer(clsize), as.double(engine@sigma0), as.integer(length(clsize)),
           as.integer(p), as.integer(n), as.double(W), as.double(gw),
           double(p), PACKAGE = "aftgee")[[12]]
    }
      if (engine@solver %in% c("BBsolve", "dfsane")) {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- tryCatch(
                do.call(engine@solver, list(par = engine@b0, fn = gp.est,
                                            quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))),
                error = function(e)
                    do.call(engine@solver, list(par = double(p), fn = gp.est,
                                                quiet = TRUE, control = list(tol = engine@tol, trace = FALSE)))))
        end.time <- Sys.time()
    }
    if (engine@solver == "BBoptim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver, list(par = engine@b0, fn = gp.obj,
                                               quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    if (engine@solver == "optim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver,
                           list(par = engine@b0, fn = gp.obj, control = list(abstol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    list(beta = fit$par, conv = fit$convergence, pe.time = end.time - start.time,
         value = c(fit$residual, fit$value), iter = 1)
}

rankFit.gp.ns <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    if (is.null(gw)) gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    pwr <- ifelse(is.null(engine@gp.pwr), 1 / ncol(X), engine@gp.pwr)
    if (pwr < 0) stop("Invalid GP class power.", call. = FALSE)
    gp.obj <- function(b) {
        er <- Y - X %*% b
        tmp <- survfit(Surv(er, delta) ~ 1, weights = W)
        gw <- approx(tmp$time, tmp$surv, er, "constant", yleft = 1, yright = min(tmp$surv))$y
        gw <- gw^pwr
        gw <- ifelse(gw == Inf, 0, gw)
        ans <- .C("log_ns_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
                  as.integer(clsize), as.integer(length(clsize)),
                  as.integer(p), as.integer(n), as.double(W), as.double(gw),
                  double(p), PACKAGE = "aftgee")[[11]]
        return(sum(ans^2))
    }
    gp.est <- function(b) {
        er <- Y - X %*% b
        tmp <- survfit(Surv(er, delta) ~ 1, weights = W)
        gw <- approx(tmp$time, tmp$surv, er, "constant", yleft = 1, yright = min(tmp$surv))$y
        gw <- gw^pwr
        gw <- ifelse(gw == Inf, 0, gw)
        .C("log_ns_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
           as.integer(clsize), as.integer(length(clsize)),
           as.integer(p), as.integer(n), as.double(W), as.double(gw),
           double(p), PACKAGE = "aftgee")[[11]]
    }
      if (engine@solver %in% c("BBsolve", "dfsane")) {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- tryCatch(
                do.call(engine@solver, list(par = engine@b0, fn = gp.est,
                                            quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))),
                error = function(e)
                do.call(engine@solver, list(par = double(p), fn = gp.est,
                                            quiet = TRUE, control = list(tol = engine@tol, trace = FALSE)))))
        end.time <- Sys.time()
    }
    if (engine@solver == "BBoptim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver, list(par = engine@b0, fn = gp.obj,
                                               quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    if (engine@solver == "optim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver,
                           list(par = engine@b0, fn = gp.obj, control = list(abstol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    list(beta = fit$par, conv = fit$convergence, pe.time = end.time - start.time,
         value = c(fit$residual, fit$value), iter = 1)
}

rankFit.gp.mns <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    b1 <- rankFit.gehan.ns(DF, engine, stdErr)$beta
    pwr <- ifelse(is.null(engine@gp.pwr), 1 / ncol(X), engine@gp.pwr)
    if (pwr < 0) stop("Invalid GP class power.", call. = FALSE)
    iter <- 1
    start.time <- Sys.time()
    for (i in 1:engine@maxIter) {
        er <- Y - X %*% b1
        tmp <- survfit(Surv(er, delta) ~ 1, weights = W)
        s0 <- approx(tmp$time, tmp$surv, er, "constant", yleft = 1, yright = min(tmp$surv))$y
        gw <- s0 / .C("gehan_ns_wt", as.double(b1), as.double(Y), as.double(X), as.integer(clsize),
                      as.integer(length(clsize)), as.integer(p), as.integer(n), as.double(W),
                      double(n), PACKAGE = "aftgee")[[9]]
        gw <- gw^pwr
        gw <- ifelse(gw == Inf, 0, gw)
        b2 <- rankFit.gehan.ns(DF, engine, stdErr, gw)$beta
        iter <- iter + 1
        if (engine@trace) print(b2)
        if (max(abs(b1 - b2)) < engine@tol) {
            conv <- 0
            break
        }
        b1 <- b2
        conv <- 1
    }
    end.time <- Sys.time()
    list(beta = b2, conv = conv, pe.time = end.time - start.time,
         value = abs(b1 - b2), iter = iter)
}

rankFit.gp.mis <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    pwr <- ifelse(is.null(engine@gp.pwr), 1 / ncol(X), engine@gp.pwr)
    if (pwr < 0) stop("Invalid GP class power.", call. = FALSE)
    b1 <- rankFit.gehan.is(DF, engine, stdErr)$beta
    iter <- 1
    start.time <- Sys.time()
    for (i in 1:engine@maxIter) {
        er <- Y - X %*% b1
        tmp <- survfit(Surv(er, delta) ~ 1, weights = W)
        s0 <- approx(tmp$time, tmp$surv, er, "constant", yleft = 1, yright = min(tmp$surv))$y
        gw <- s0 / .C("gehan_s_wt", as.double(b1), as.double(Y), as.double(X), as.integer(clsize),
                      as.double(engine@sigma0), as.integer(length(clsize)),
                      as.integer(p), as.integer(n), as.double(W),
                      double(n), PACKAGE = "aftgee")[[10]]
        gw <- gw^pwr
        gw <- ifelse(gw == Inf, 0, gw)
        b2 <- rankFit.gehan.is(DF, engine, stdErr, gw)$beta
        iter <- iter + 1
        if (engine@trace) print(b2)
        if (max(abs(b1 - b2)) < engine@tol) {
            conv <- 0
            break
        }
        b1 <- b2
        conv <- 1
    }
    end.time <- Sys.time()
    list(beta = b2, conv = conv, pe.time = end.time - start.time,
         value = abs(b1 - b2), iter = iter)
}

rankFit.user.is <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    if (is.null(gw)) gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    user.obj <- function(b) {
        gw <- engine@userRk
        ans <- .C("log_s_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
                  as.integer(clsize), as.double(engine@sigma0), as.integer(length(clsize)),
                  as.integer(p), as.integer(n), as.double(W), as.double(gw),
                  double(p), PACKAGE = "aftgee")[[12]]
        return(sum(ans^2))
    }
    user.est <- function(b) {
        gw <- engine@userRk
        .C("log_s_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
           as.integer(clsize), as.double(engine@sigma0), as.integer(length(clsize)),
           as.integer(p), as.integer(n), as.double(W), as.double(gw),
           double(p), PACKAGE = "aftgee")[[12]]
    }
    if (engine@solver %in% c("BBsolve", "dfsane")) {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- tryCatch(
                do.call(engine@solver, list(par = engine@b0, fn = user.est,
                                            quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))),
                error = function(e)
                do.call(engine@solver, list(par = double(p), fn = user.est,
                                            quiet = TRUE, control = list(tol = engine@tol, trace = FALSE)))))
                end.time <- Sys.time()
    }
    if (engine@solver == "BBoptim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver, list(par = engine@b0, fn = user.obj,
                                               quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    if (engine@solver == "optim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver,
                           list(par = engine@b0, fn = user.obj,
                                control = list(abstol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    list(beta = fit$par, conv = fit$convergence, pe.time = end.time - start.time,
         value = c(fit$residual, fit$value), iter = 1)
}

rankFit.user.ns <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    if (is.null(gw)) gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    user.obj <- function(b) {
        gw <- engine@userRk
        ans <- .C("log_ns_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
                  as.integer(clsize), as.integer(length(clsize)),
                  as.integer(p), as.integer(n), as.double(W), as.double(gw),
                  double(p), PACKAGE = "aftgee")[[11]]
        return(sum(ans^2))
    }
    user.est <- function(b) {
        gw <- engine@userRk
        .C("log_ns_est", as.double(b), as.double(Y), as.double(X), as.double(delta),
           as.integer(clsize), as.integer(length(clsize)),
           as.integer(p), as.integer(n), as.double(W), as.double(gw),
           double(p), PACKAGE = "aftgee")[[11]]
    }
      if (engine@solver %in% c("BBsolve", "dfsane")) {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- tryCatch(
                do.call(engine@solver, list(par = engine@b0, fn = user.est,
                                            quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))),
                error = function(e)
                do.call(engine@solver, list(par = double(p), fn = user.est,
                                            quiet = TRUE, control = list(tol = engine@tol, trace = FALSE)))))
                end.time <- Sys.time()
    }
    if (engine@solver == "BBoptim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver, list(par = engine@b0, fn = user.obj,
                                               quiet = TRUE, control = list(tol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    if (engine@solver == "optim") {
        start.time <- Sys.time()
        suppressWarnings(
            fit <- do.call(engine@solver,
                           list(par = engine@b0, fn = user.obj,
                                control = list(abstol = engine@tol, trace = FALSE))))
        end.time <- Sys.time()
    }
    list(beta = fit$par, conv = fit$convergence, pe.time = end.time - start.time,
         value = c(fit$residual, fit$value), iter = 1)  
}

rankFit.user.mns <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    b1 <- rankFit.gehan.ns(DF, engine, stdErr)$beta
    iter <- 1
    start.time <- Sys.time()
    for (i in 1:engine@maxIter) {
        gw <- engine@userRk
        b2 <- rankFit.gehan.ns(DF, engine, stdErr, gw)$beta
        iter <- iter + 1
        if (engine@trace) print(b2)
        if (max(abs(b1 - b2)) < engine@tol) {
            conv <- 0
            break
        }
        b1 <- b2
        conv <- 1
    }
    end.time <- Sys.time()
    list(beta = b2, conv = conv, pe.time = end.time - start.time,
         value = abs(b1 - b2), iter = iter)
}

rankFit.user.mis <- function(DF, engine, stdErr, gw = NULL) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    gw <- rep(1, n)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    b1 <- rankFit.gehan.is(DF, engine, stdErr)$beta
    iter <- 1
    start.time <- Sys.time()
    for (i in 1:engine@maxIter) {
        gw <- engine@userRk
        b2 <- rankFit.gehan.is(DF, engine, stdErr, gw)$beta
        iter <- iter + 1
        if (engine@trace) print(b2)
        if (max(abs(b1 - b2)) < engine@tol) {
            conv <- 0
            break
        }
        b1 <- b2
        conv <- 1
    }
    end.time <- Sys.time()
    list(beta = b2, conv = conv, pe.time = end.time - start.time,
         value = abs(b1 - b2), iter = iter)
}


##############################################################################
## function calls for different variance estimations
##############################################################################

rankFit.Engine.bootstrap <- function(DF, engine, stdErr, gw) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    gw <- rep(1, n)
    p <- ncol(X)
    clsz <- as.numeric(unlist(lapply(split(id, id), length)))
    fit <- rankFit(DF = DF, engine = engine, stdErr = NULL, gw = NULL)
    engine@b0 <- fit$beta
    B <- stdErr@B
    betaVar <- matrix(0, B, p)
    uID <- unique(DF$id)
    if (stdErr@parallel) {
        cl <- makeCluster(stdErr@parCl)
        clusterExport(cl = cl,
                      varlist = c("DF", "engine", "stdErr", "id"),
                      envir = environment())
        out <- parSapply(cl, 1:B, function(x) {
            sampled.id <- sample(unique(id), n, TRUE)
            ind <- unlist(sapply(sampled.id, function(x) which(id == x)))
            DF2 <- DF[ind,]
            DF2$id <- rep(1:n, clsz[sampled.id])
            rankFit(DF = DF2, engine = engine, stdErr = NULL, gw = NULL)$beta})
        stopCluster(cl)
        betaVar <- t(out)
        return(c(fit, list(betaVar = varOut(betaVar))))
    }
    for (i in 1:B) {
        sampled.id <- sample(unique(id), n, TRUE)
        ind <- unlist(sapply(sampled.id, function(x) which(id == x)))
        DF2 <- DF[ind,]
        DF2$id <- rep(1:n, clsz[sampled.id])
        betaVar[i,] <- rankFit(DF = DF2, engine = engine, stdErr = NULL, gw = NULL)$beta
    }
    return(c(fit, list(betaVar = varOut(betaVar))))
}

rankFit.Engine.mb <- function(DF, engine, stdErr, gw) {
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    gw <- rep(1, n)
    p <- ncol(X)
    clsz <- as.numeric(unlist(lapply(split(id, id), length)))
    fit <- rankFit(DF = DF, engine = engine, stdErr = NULL, gw = NULL)
    engine@b0 <- fit$beta
    B <- stdErr@B
    betaVar <- matrix(0, B, p)
    uID <- unique(DF$id)
    if (stdErr@parallel) {
        cl <- makeCluster(stdErr@parCl)
        clusterExport(cl = cl,
                      varlist=c("DF", "engine", "stdErr"),
                      envir = environment())
        out <- unlist(parLapply(cl, 1:B, function(x) {
            DF2 <- DF
            Z <- rep(rexp(length(clsz)), clsz)
            DF2$weights <- DF2$weights * Z
            rankFit(DF = DF2, engine = engine, stdErr = NULL, gw = NULL)$beta}))
        stopCluster(cl)
        betaVar <- t(matrix(out, p))
        return(c(fit, list(betaVar = varOut(betaVar))))
    }
    for (i in 1:B) {
        DF2 <- DF
        Z <- rep(rexp(length(clsz)), clsz)
        DF2$weights <- DF2$weights * Z
        betaVar[i,] <- rankFit(DF = DF2, engine = engine, stdErr = NULL, gw = NULL)$beta
    }
    return(c(fit, list(betaVar = varOut(betaVar))))
}

## rankFit.logrank.zl <- function(DF, engine, stdErr, gw) {
##     id <- DF$id
##     Y <- log(DF$time)
##     delta <- DF$status
##     X <- as.matrix(DF[,-(1:4)])
##     W <- DF$weights
##     n <- nrow(DF)
##     gw <- rep(1, n)
##     p <- ncol(X)
##     clsz <- as.numeric(unlist(lapply(split(id, id), length)))
##     fit <- rankFit(DF = DF, engine = engine, stdErr = NULL, gw = NULL)
##     engine@b0 <- fit$beta
##     B <- stdErr@B
##     Z1 <- replicate(B, rep(rexp(length(clsize)), clsize))
##     Z2 <- replicate(B, rnorm(p))
##     An <- matrix(0, p, p)
##     ## Estimate A
##     if (grep("ns", class(engine)) > 0) {
##         log.est.A <- function(Z) {
##             .C("log_ns_est", as.double(fit$beta + Z / sqrt(n)),
##                as.double(Y), as.double(X), as.double(delta),
##                as.integer(clsize), as.integer(length(clsize)),
##                as.integer(p), as.integer(n), as.double(W), as.double(gw),
##                out = double(p), PACKAGE = "aftgee")$out
##         }
##         Un <- t(apply(Z2, 2, log.est.A)) 
##     }
##     if (grep("is", class(engine)) > 0) {
##         log.est.A <- function(Z) {
##             .C("log_s_est", as.double(fit$beta + Z / sqrt(n)),
##                as.double(Y), as.double(X), as.double(delta),
##                as.integer(clsize), as.double(engine@sigma0), as.integer(length(clsize)),
##                as.integer(p), as.integer(n), as.double(W), as.double(gw),
##                out = double(p), PACKAGE = "aftgee")$out
##         }
##         Un <- t(apply(Z2, 2, log.est.A)) 
##     }
##     for (i in 1:p) {
##         An[i,] <- lm(Un[,i] ~ I(t(fit$beta + Z2 / sqrt(n))))$coef
##     }
##     ## Estimate V
##     if (grep("mb", class(stdErr)) > 0) {
##         if (grep("ns", class(engine)) > 0) {
##             log.est.V <- function(Z) {
##                 .C("log_ns_est", as.double(fit$beta),
##                    as.double(Y), as.double(X), as.double(delta),
##                    as.integer(clsize), as.integer(length(clsize)),
##                    as.integer(p), as.integer(n), as.double(W * Z), as.double(gw),
##                out = double(p), PACKAGE = "aftgee")$out
##             }
##             Vn <- varOut(t(apply(Z1, 2, log.est.V)))
##         }
##         if (grep("is", class(engine)) > 0) {
##             log.est.V <- function(Z) {
##                 .C("log_s_est", as.double(fit$beta),
##                    as.double(Y), as.double(X), as.double(delta),
##                    as.integer(clsize), as.double(engine@sigma0), as.integer(length(clsize)),
##                    as.integer(p), as.integer(n), as.double(W * Z), as.double(gw),
##                    out = double(p), PACKAGE = "aftgee")$out
##             }
##             Vn <- varOut(t(apply(Z2, 2, log.est.A)))
##         }
##     }
##     if (grep("cf", class(stdErr)) > 0) {
##         Vn <- viClo(fit$beta, Y, delta, X, id, weights, B, "logrank")$vi * n
##     }
##     if (qr(An)$rank != p) {
##         covmat <- ginv(An) %*% vi %*% t(ginv(An))
##         An.msg <- "An is singular"
##         An.inv <- 0
##     }
##     if (qr(An)$rank == p) {
##         covmat <- solve(An) %*% vi %*% solve(An)
##         An.msg <- "An is nonsingular"
##         An.inv <- 1
##     }
##     covmat <- covmat  / n
##     return(matrix(as.numeric(covmat), p))
## }

##############################################################################
# Class Definition
##############################################################################

## Point estimator
setClass("Engine",
         representation(tol = "numeric", b0 = "numeric", sigma0 = "matrix",
                        userRk = "numeric", maxIter = "numeric",
                        solver = "character", trace = "logical"),
         prototype(tol = 1e-3, b0 = 0, sigma0 = matrix(0), userRk = 0,
                   maxIter = 50, solver = "dfsane", trace = FALSE),
         contains = "VIRTUAL")

setClass("gehan.is", contains = "Engine")
setClass("gehan.ns", contains = "Engine")
setClass("gehan.mis", contains = "Engine")
setClass("gehan.mns", contains = "Engine")

setClass("logrank.is", contains = "Engine")
setClass("logrank.ns", contains = "Engine")
setClass("logrank.mis", contains = "Engine")
setClass("logrank.mns", contains = "Engine")

setClass("pw.is", contains = "Engine")
setClass("pw.ns", contains = "Engine")
setClass("pw.mis", contains = "Engine")
setClass("pw.mns", contains = "Engine")

setClass("user.is", contains = "Engine")
setClass("user.ns", contains = "Engine")
setClass("user.mis", contains = "Engine")
setClass("user.mns", contains = "Engine")

setClass("gp.ns",
         representation(gp.pwr = "numeric"),
         prototype(gp.pwr = NULL),
         contains = "Engine")
setClass("gp.is",
         representation(gp.pwr = "numeric"),
         prototype(gp.pwr = NULL),
         contains = "Engine")
setClass("gp.mis",
         representation(gp.pwr = "numeric"),
         prototype(gp.pwr = NULL),
         contains = "Engine")
setClass("gp.mns",
         representation(gp.pwr = "numeric"),
         prototype(gp.pwr = NULL),
         contains = "Engine")

## Variance
setClass("stdErr",
         representation(tol = "numeric", B = "numeric", parallel = "logical", parCl = "numeric"),
         prototype(tol = 1e-3, B = 100, parallel = FALSE, parCl = parallel::detectCores() / 2),
         contains = "VIRTUAL")

setClass("bootstrap", contains = "stdErr")
setClass("MB", contains = "stdErr")
setClass("ZLCF", contains = "stdErr")
setClass("ZLMB", contains = "stdErr")
setClass("sHCF", contains = "stdErr")
setClass("sHMB", contains = "stdErr")
setClass("ISCF", contains = "stdErr")
setClass("ISMB", contains = "stdErr")

##############################################################################
# Method Dispatch
##############################################################################
## Generic functins
setGeneric("rankFit", function(DF, engine, stdErr, gw) standardGeneric("rankFit"))
## setGeneric("geeFit", function(DF, engine, stdErr) standardGeneric("geeFit"))

## Point estimator
setMethod("rankFit", signature(engine = "gehan.is", stdErr = "NULL"), rankFit.gehan.is)
setMethod("rankFit", signature(engine = "gehan.ns", stdErr = "NULL"), rankFit.gehan.ns)
setMethod("rankFit", signature(engine = "gehan.mns", stdErr = "NULL"), rankFit.gehan.ns)
setMethod("rankFit", signature(engine = "gehan.mis", stdErr = "NULL"), rankFit.gehan.is)

setMethod("rankFit", signature(engine = "logrank.is", stdErr = "NULL"), rankFit.logrank.is)
setMethod("rankFit", signature(engine = "logrank.ns", stdErr = "NULL"), rankFit.logrank.ns)
setMethod("rankFit", signature(engine = "logrank.mns", stdErr = "NULL"), rankFit.logrank.mns)
setMethod("rankFit", signature(engine = "logrank.mis", stdErr = "NULL"), rankFit.logrank.mis)

setMethod("rankFit", signature(engine = "pw.is", stdErr = "NULL"), rankFit.pw.is)
setMethod("rankFit", signature(engine = "pw.ns", stdErr = "NULL"), rankFit.pw.ns)
setMethod("rankFit", signature(engine = "pw.mns", stdErr = "NULL"), rankFit.pw.mns)
setMethod("rankFit", signature(engine = "pw.mis", stdErr = "NULL"), rankFit.pw.mis)

setMethod("rankFit", signature(engine = "gp.is", stdErr = "NULL"), rankFit.gp.is)
setMethod("rankFit", signature(engine = "gp.ns", stdErr = "NULL"), rankFit.gp.ns)
setMethod("rankFit", signature(engine = "gp.mns", stdErr = "NULL"), rankFit.gp.mns)
setMethod("rankFit", signature(engine = "gp.mis", stdErr = "NULL"), rankFit.gp.mis)

setMethod("rankFit", signature(engine = "user.is", stdErr = "NULL"), rankFit.user.is)
setMethod("rankFit", signature(engine = "user.ns", stdErr = "NULL"), rankFit.user.ns)
setMethod("rankFit", signature(engine = "user.mns", stdErr = "NULL"), rankFit.user.mns)
setMethod("rankFit", signature(engine = "user.mis", stdErr = "NULL"), rankFit.user.mis)

## Variance
setMethod("rankFit", signature(engine = "Engine", stdErr = "bootstrap"), rankFit.Engine.bootstrap)
setMethod("rankFit", signature(engine = "Engine", stdErr = "MB"), rankFit.Engine.mb)
setMethod("rankFit", signature(engine = "Engine", stdErr = "ZLCF"), rankFit.Engine.bootstrap)
setMethod("rankFit", signature(engine = "Engine", stdErr = "ZLMB"), rankFit.Engine.bootstrap)
setMethod("rankFit", signature(engine = "Engine", stdErr = "sHCF"), rankFit.Engine.bootstrap)
setMethod("rankFit", signature(engine = "Engine", stdErr = "sHMB"), rankFit.Engine.bootstrap)
setMethod("rankFit", signature(engine = "Engine", stdErr = "ISCF"), rankFit.Engine.bootstrap)
setMethod("rankFit", signature(engine = "Engine", stdErr = "ISMB"), rankFit.Engine.bootstrap)

##############################################################################
## User's Main Function
##############################################################################

#' Accelerated Failure Time with Smooth Rank Regression
#'
#' Fits a semiparametric accelerated failure time (AFT) model with rank-based approach.
#' General weights, additional sampling weights and fast sandwich variance estimations
#' are also incorporated.
#' Estimating equations are solved with Barzilar-Borwein spectral method implemented as
#' \code{BBsolve} in package \pkg{BB}.
#'
#' @param formula a formula expression, of the form \code{response ~ predictors}.
#'     The \code{response} is a \code{Surv} object object with right censoring.
#'     See the documentation of \code{lm}, \code{coxph} and \code{formula} for details.
#' @param data an optional data frame in which to interpret the variables
#'     occurring in the \code{formula}.
#' @param subset an optional vector specifying a subset of observations
#'     to be used in the fitting process.
#' @param id an optional vector used to identify the clusters.
#'     If missing, then each individual row of \code{data} is presumed to
#'     represent a distinct subject.
#'     The length of \code{id} should be the same as the number of observation.
#' @param contrasts an optional list.
#' @param weights an optional vector of observation weights.
#' @param B a numeric value specifies the resampling number.
#'     When B = 0, only the beta estimate will be displayed.
#' @param rankWeights a character string specifying the type of general weights.
#'     The following are permitted:
#' \describe{
#'   \item{\code{logrank}}{logrank weight}
#'   \item{\code{gehan}}{Gehan's weight}
#'   \item{\code{PW}}{Prentice-Wilcoxon weight}
#'   \item{\code{GP}}{GP class weight}
#' }
#' @param eqType a character string specifying the type of the
#'     estimating equation used to obtain the regression parameters.
#'     The following are permitted:
#' \describe{
#'   \item{\code{nonsm}}{Regression parameters are estimated by directly solving the nonsmooth
#' estimating equations.}
#'   \item{\code{sm}}{Regression parameters are estimated by directly solving the
#' induced-smoothing estimating equations.}
#'   \item{\code{monosm}}{Regression parameters are estimated by iterating the
#'   monotonic smoothed estimating equations. This is typical when
#'   \code{rankWeights = "PW"} and \code{rankWeights = "GP"}.}
#' }
#' @param se a character string specifying the estimating method for the variance-covariance matrix.
#'   The following are permitted:
#' \describe{
#'   \item{\code{bootstrap}}{nonparametric bootstrap,}
#'   \item{\code{MB}}{multiplier resampling.}
#'   \item{\code{ZLCF}}{Zeng and Lin's approach with closed form Si.}
#'   \item{\code{ZLMB}}{Zeng and Lin's approach with empirical Si.}
#'   \item{\code{sHCF}}{Huang's approach with closed form Si.}
#'   \item{\code{sHMB}}{Huang's approach with empirical Si.}
#'   \item{\code{ISCF}}{Johnson and Strawderman's sandwich variance estimates with closed form Si.}
#'   \item{\code{ISMB}}{Johnson and Strawderman's sandwich variance estimates with empirical Si.}
#'   \item{\code{js}}{Johnson and Strawderman's iterating approach.}
#' }
#' @param control controls equation solver, maxiter, tolerance, and resampling variance estimation.
#' The available equation solvers are \code{BBsolve} and \code{dfsane} of the \pkg{BB} package.
#' Instead of searching for the zero crossing, options including \code{BBoptim} and \code{optim}
#' will return solution from maximizing the corresponding objective function.
#'
#' @export
#'
#' @return \code{aftsrr} returns an object of class "\code{aftsrr}" representing the fit.
#' An object of class "\code{aftsrr}" is a list containing at least the following components:
#' \describe{
#'   \item{beta}{A vector of beta estimates}
#'   \item{covmat}{A list of covariance estimates}
#'   \item{convergence}{An integer code indicating type of convergence.}
#'   \describe{
#'     \item{0}{indicates successful convergence.}
#'     \item{1}{indicates that the iteration limit \code{maxit} has been reached.}
#'     \item{2}{indicates failure due to stagnation.}
#'     \item{3}{indicates error in function evaluation.}
#'     \item{4}{is failure due to exceeding 100 step length reductions in line-search.}
#'     \item{5}{indicates lack of improvement in objective function.}
#'   }
#'   \item{bhist}{When \code{variance = "MB"}, \code{bhist} gives the bootstrap samples.}
#' }
#'
#' @references Chiou, S., Kang, S. and Yan, J. (2014)
#' Fast Accelerated Failure Time Modeling for Case-Cohort Data. \emph{Statistics and Computing},
#' \bold{24}(4): 559--568.
#' @references Chiou, S., Kang, S. and Yan, J. (2014)
#' Fitting Accelerated Failure Time Model in Routine Survival Analysis with {R} Package \pkg{Aftgee}.
#' \emph{Journal of Statistical Software}, \bold{61}(11): 1--23.
#' @references Huang, Y. (2002) Calibration Regression of Censored Lifetime Medical Cost.
#' \emph{Journal of American Statistical Association}, \bold{97}, 318--327.
#' @references Johnson, L. M. and Strawderman, R. L. (2009)
#' Induced Smoothing for the Semiparametric Accelerated Failure Time Model:
#' Asymptotic and Extensions to Clustered Data. \emph{Biometrika}, \bold{96}, 577 -- 590.
#' @references Varadhan, R. and Gilbert, P. (2009)
#' BB: An R Package for Solving a Large System of Nonlinear Equations and
#' for Optimizing a High-Dimensional Nonlinear Objective Function.
#' \emph{Journal of Statistical Software}, \bold{32}(4): 1--26
#' @references Zeng, D. and Lin, D. Y. (2008)
#' Efficient Resampling Methods for Nonsmooth Estimating Functions.
#' \emph{Biostatistics}, \bold{9}, 355--363
#'
#' @export
#' @examples
#' ## kidney data
#' library(survival)
#' data(kidney)
#' foo <- aftsrr(Surv(time, status) ~ age + sex, id = id,
#'                 data = kidney, se = c("ISMB", "ZLMB"), B = 10)
#' foo
#'
#' ## nwtco data
#' library(survival)
#' data(nwtco)
#' subinx <- sample(1:nrow(nwtco), 668, replace = FALSE)
#' nwtco$subcohort <- 0
#' nwtco$subcohort[subinx] <- 1
#' pn <- table(nwtco$subcohort)[[2]] / sum(table(nwtco$subcohort))
#' nwtco$hi <- nwtco$rel + ( 1 - nwtco$rel) * nwtco$subcohort / pn
#' nwtco$age12 <- nwtco$age / 12
#' nwtco$study <- nwtco$study - 3
#' nwtco$histol = nwtco$histol - 1
#' sub <- nwtco[subinx,]
#' fit <- aftsrr(Surv(edrel, rel) ~ histol + age12 + study, id = seqno,
#'        weights = hi, data = sub, B = 10, se = c("ISMB", "ZLMB"),
#'        subset = stage == 4)
#' summary(fit)
aftsrr <- function(formula, data, subset, id = NULL, contrasts = NULL, 
                   weights = NULL, B = 100, 
                   rankWeights = c("gehan", "logrank", "pw", "gp", "userdefined"),
                   eqType = c("is", "ns", "mns", "mis"),
                   se = c("NULL", "bootstrap", "MB", "ZLCF", "ZLMB",
                          "sHCF", "sHMB", "ISCF", "ISMB"),
                   control = list()) {
    rkWeights <- match.arg(rankWeights)
    eqType <- match.arg(eqType)
    se0 <- c("NULL", "bootstrap", "MB", "ZLCF", "ZLMB", "sHCF", "sHMB", "ISCF", "ISMB")
    if (is.null(se)) se <- "NULL"
    se <- se0[pmatch(se, se0)]
    if (length(se) == length(se0)) se <- "NULL"
    if (is.na(se)[1])
        stop("Invalid variance estimator.", call. = FALSE)
    scall <- match.call()
    mnames <- c("", "formula", "data", "weights", "subset", "na.nation", "id")
    cnames <- names(scall)
    cnames <- cnames[match(mnames, cnames, 0)]
    mcall <- scall[cnames]
    mcall[[1]] <- as.name("model.frame")
    m <- eval(mcall, parent.frame())
    id <- model.extract(m, id)
    mterms <- attr(m, "terms")
    weights <- model.extract(m, weights) 
    ## if (missing(data)) obj <- eval(formula[[2]], parent.frame())
    ## if (!missing(data)) obj <- eval(formula[[2]], data)
    obj <- unclass(m[,1]) 
    if (class(m[[1]]) != "Surv" || ncol(obj) > 2)
        stop("aftsrr only supports Surv object with right censoring.", call. = FALSE)
    if (is.null(id)) id <- 1:nrow(obj)
    if (is.null(weights)) weights <- rep(1, nrow(obj))
    formula[[2]] <- NULL
    ## Create DF; the first 2 columns are from Surv with time and status
    if (formula == ~1) DF <- cbind(obj, zero = 0)
    else {
        DF <- cbind(obj, id, weights, model.matrix(mterms, m, contrasts))
        ## remove intercept
        ## if (sum(colnames(DF) == "(Intercept)") > 0)
        ## print("Rank-based estimation assumes no-intercept model.")
        if (sum(colnames(DF) == "(Intercept)") > 0)
            DF <- DF[,-which(colnames(DF) == "(Intercept)")]
    }
    DF <- as.data.frame(DF)
    ## setup engine
    method <- paste(rkWeights, ".", eqType, sep = "")
    engine.control <- control[names(control) %in% names(attr(getClass(method), "slots"))]
    engine <- do.call("new", c(list(Class = method), engine.control))
    if (engine@b0 == 0) {
        lm.formula <- paste("log(time)", paste(formula, collapse = ""))
        engine@b0 <- as.numeric(coef(lm(lm.formula, data = DF[,-(2:4)])))[-1]
    }
    if (length(engine@b0) != ncol(DF) - 4) 
        stop ("Initial value length does not match with the numbers of covariates", call. = FALSE)
    if (engine@sigma0 == 0) engine@sigma0 <- diag(length(engine@b0))
    if (rkWeights == "userdefined" & length(engine@userRk) != nrow(DF))
        stop("Invalid userdefined rank weight values.", call. = FALSE)
    ## #####################################################################################
    ## Easy patch for now (3/15/2018);
    stdErr.control <- control[names(control) %in% names(attr(getClass(se[1]), "slots"))]
    stdErr <- do.call("new", c(list(Class = se[1]), stdErr.control))
    if (se[1] != "NULL") stdErr@B <- B
    id <- DF$id
    Y <- log(DF$time)
    delta <- DF$status
    X <- as.matrix(DF[,-(1:4)])
    W <- DF$weights
    n <- nrow(DF)
    p <- ncol(X)
    clsize <- as.numeric(unlist(lapply(split(id, id), length)))
    if (sum(se %in% "NULL")) {
        fit <- rankFit(DF = DF, engine = engine, stdErr = NULL, gw = NULL)
        covmat <- matrix(NA, p, p)
    } else {
        fit <- rankFit(DF = DF, engine = engine, stdErr = NULL, gw = NULL)
        ZLMB.An.inv <- ZLCF.An.inv <- ISMB.An.inv <- ISCF.An.inv <- js.An.inv <- 1
        vBoot <- vMB <- vZLCF <- vZLMB <- vsHCF <- vsHMB <- vISCF <- vISMB <- bstep <- NaN
        if (sum(se %in% "bootstrap") > 0) {        
            fit <- rankFit(DF = DF, engine = engine, stdErr = stdErr, gw = NULL)
            vBoot <- fit$betaVar
        }
        if (sum(se %in% "MB") > 0) {
            fit <- rankFit(DF = DF, engine = engine, stdErr = stdErr, gw = NULL)
            vMB <- fit$betaVar
        }
        if (sum(se %in% c("ISMB", "ISCF", "ZLMB", "ZLCF", "sZLCF", "sHMB", "sHCF")) > 0) {
            gw <- getGw(Y = Y, X = X, beta = fit$beta, N = nrow(DF), delta = delta,
                        clsize = clsize, sigma = engine@sigma0, weights = W,
                        rankWeights = rkWeights)
        }
        if (sum(se %in% "ZLCF") > 0) {
            ## ZLCF is only for gehan weight
            vZLCF <- zlFun(beta = fit$beta, Y = Y, delta = delta, X = X,
                           id = id, weights = W, B = B, vClose = TRUE,
                           rankWeights = rkWeights, gw = rep(1, nrow(X)),
                           method = method, stratify = TRUE, sigma = engine@sigma0)
            ZLCF.An.inv <- vZLCF$An.inv
            vZLCF <- vZLCF$covmat
        }        
        if (sum(se %in% "ZLMB") > 0) {
            if (B > 0) {
                vZLMB <- zlFun(beta = fit$beta, Y = Y, delta = delta, X = X,
                               id = id, weights = W, B = B, vClose = FALSE,
                               rankWeights = rkWeights, gw = gw,
                               method = method, stratify = TRUE, sigma = engine@sigma0)
                ZLMB.An.inv <- vZLMB$An.inv
                vZLMB <- vZLMB$covmat
            }
        }
        if (sum(se %in% "ISCF") > 0) {
            if (B > 0) {
                ## ISCF only for gehan
                vISCF <- isFun(beta = fit$beta, Y = Y, delta = delta, X = X, id = id,
                               weights = W, sigma = engine@sigma0, B = B, vClose = TRUE,
                               gw = gw, rankWeights = rkWeights, stratify = TRUE)
                ISCF.An.inv <- vISCF$An.inv
                vISCF <- vISCF$covmat
            }
        }
        if (sum(se %in% "ISMB") > 0) {
            if (B > 0) {
                vISMB <- isFun(beta = fit$beta, Y = Y, delta = delta, X = X, id = id,
                               weights = W, sigma = engine@sigma0, B = B, vClose = FALSE,
                               gw = gw, rankWeights = rkWeights)
                ISMB.An.inv <- vISMB$An.nv
                vISMB <- vISMB$covmat
            }
        }
        if (sum(se %in% "sHCF") > 0) {
            if (B > 0) {
                vsHCF <- huangFun(beta = fit$beta, Y = Y, delta = delta, X = X, id = id,
                                  weights = W, sigma = engine@sigma0, B = B, vClose = TRUE,
                                  gw = gw, rankWeights = rkWeights, stratify = TRUE)$covmat
            }
        }
        if (sum(se %in% "sHMB") > 0) {
            if (B > 0) {
                vsHMB <- huangFun(beta = fit$beta, Y = Y, delta = delta, X = X, id = id,
                                  weights = W, sigma = engine@sigma0, B = B, vClose = FALSE,
                                  gw = gw, rankWeights = rkWeights, stratify = TRUE)$covmat
            }
        }
        covmat <- list(bootstrap = vBoot, 
                       MB = vMB, ZLCF = vZLCF, ZLMB = vZLMB,
                       sHCF = vsHCF, sHMB = vsHMB,
                       ISCF = vISCF, ISMB = vISMB)        
    }
    out <- list(beta = fit$beta, covmat = covmat, convergence = fit$conv, bstep = fit$iter,
                var.meth = se, bhist = NULL)
    class(out) <- "aftsrr"
    out$call <- scall
    out$vari.name <- colnames(X)
    out$B <- B
    out$binit <- engine@b0
    return(out)
}

##############################################################################
## Background functions
##############################################################################

varOut <- function(dat, na.rm = TRUE) {
    dat[which(dat %in% boxplot(dat, plot = FALSE)$out)] <- NA
    dat <- dat[complete.cases(dat),]
    var(dat, na.rm = na.rm)
}

##############################################################################

abargehanfun <- function(beta, Y, X, delta, clsize, sigma, weights, gw = rep(1, nrow(X))) {
    p <- ncol(X)
    N <- nrow(X)
    n <- length(clsize)
    a <- vector("double", p * p)
    matrix(.C("abargehanfunC", as.double(beta), as.double(Y), as.double(X), as.double(delta),
              as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p),
              as.integer(N), as.double(weights), as.double(gw),
              as.double(a), PACKAGE = "aftgee")[[12]], nrow = p)
}


abarlogfun <- function(beta, Y, X, delta, clsize, sigma, weights, pw = rep(1, nrow(X))) {
    p <- ncol(X)
    N <- nrow(X)
    n <- length(clsize)
    a <- vector("double", p * p)
    matrix(.C("abarlogfunC", as.double(beta), as.double(Y), as.double(X), as.double(delta),
              as.integer(clsize), as.double(pw), as.double(sigma), as.integer(n),
              as.integer(p), as.integer(N), as.double(weights),
              as.double(a), PACKAGE = "aftgee")[[12]], nrow = p)
}

abarpwfun <- function(beta, Y, X, delta, clsize, sigma, weights, pw) {
    p <- ncol(X)
    N <- nrow(X)
    n <- length(clsize)
    a <- vector("double", p * p)
    pt1 <- matrix(.C("abarpwfunC", as.double(beta), as.double(Y), as.double(X),
                     as.double(delta), as.integer(clsize), as.double(pw$fhat),
                     as.double(sigma), as.integer(n), as.integer(p), as.integer(N),
                     as.double(weights), as.double(a), PACKAGE = "aftgee")[[11]], p)
    ## pt1 <- uilogFun(beta, Y, X, delta, clsize, sigma, n, Z = rep(1, nrow(X)), weights, smooth = TRUE, constant = 0, s = 0, pw = pw$fhat)
    pt2 <- abarlogfun(beta, Y, X, delta, clsize, sigma, weights, pw$Shat)
    ## rep(1, p) %o% pt1 + pt2
    ## diag(pt1) + pt2
    pt1 + pt2
}

omegaFun <- function(beta, Y, X, delta, clsize, weights) {
  p <- ncol(X)
  N <- nrow(X)
  n <- length(clsize)
  omega <- vector("double", p * p)
  matrix(.C("omegafun", as.double(beta), as.double(Y), as.double(X),
            as.double(delta), as.integer(clsize), as.integer(n), as.integer(p),
            as.integer(N), as.double(weights),
            as.double(omega), PACKAGE = "aftgee")[[10]], nrow = p)
}

getRankName <- function(rankWeights, method) {
    rktemp <- rankWeights
    if (method == "nonsm") {
        if (rktemp == "logrank") {rankWeights <- "nslogrank"}
        if (rktemp == "PW") {rankWeights <- "nsPW"}
        if (rktemp == "GP") {rankWeights <- "nsGP"}
    }
    if (method == "monosm") {
        if (rktemp == "PW") {
            rankWeights <- "mPW"
        }
        if (rktemp == "GP") {
            rankWeights <- "mGP"
        }
        if (rktemp == "logrank") {rankWeights <- "mlogrank"}
    }
    if (method == "sm") {
        if (rktemp == "gehan") {rankWeights <- "gehan"}
        if (rktemp == "logrank") {rankWeights <- "logrank"}
    }
    out <- list(rankWeights = rankWeights)
}

getSuv <- function(Y, X, beta, N, delta, weights) {
    en <- Y - X %*% beta
    Shat <- fhat <- NULL
    dummy <- 1:N
    ord <- order(en)
    ei <- en[ord]
    weightsi <- weights[ord]
    deltai <- delta[ord]
    repeats <- table(ei)
    ## Shati <- survfit(Surv(ei, deltai) ~ 1, weights = weightsi)$surv
    Shati <- exp(-1 * basehaz(coxph(Surv(ei, deltai)~1, weights = weightsi))$hazard)
    ## fhati <- diff(c(Shati[1], Shati, Shati[length(Shati)]), 2) / diff(c(ei[1], ei, ei[length(ei)]), 2)
    Shati <- rep(Shati, repeats)
    Shatlast <- rev(Shati)[1]
    ## Shatlast <- Shati[1]
    ## fhati <- rep(fhati, repeats)
    Shat[dummy[ord]] <- Shati
    ## fhat[dummy[ord]] <- fhati
    ## fhat <- ifelse(fhat < -1, -1, fhat)
    ## mu <- mean(subset(fhat, fhat > -1))
    ## fhat <- ifelse(fhat < -1, mu, fhat)
    ## list(Shat = Shat, fhat = fhat, Shatlast = Shatlast)
    list(Shat = Shat, Shatlast = Shatlast)
}

getWi <- function(Y, X, beta, N, delta, weights) {
    en <- Y - X %*% beta
    shat <- fhat <- NULL
    dummy <- 1:N
    ord <- order(en)
    ei <- en[ord]
    weightsi <- weights[ord]
    deltai <- delta[ord]
    repeats <- table(ei)
    Shati <- survfit(Surv(ei, deltai) ~ 1, weights = weightsi)$surv
    Shati <- rep(Shati, repeats)
    ## shat <- getSuv(Y, X, beta, N, delta, weights)$Shat
    what <- NULL
    whati <- cumsum(Shati)
    shat[dummy[ord]] <- Shati
    what[dummy[ord]] <- whati
    list(what = what, shat = shat)
}
   
getGehan <- function(Y, X, beta, N, delta, clsize, sigma, weights, smooth = FALSE) {
    p <- ncol(X)
    N <- nrow(X)
    n <- length(clsize)
    a <- vector("double", N)
    if (smooth == TRUE) {
        out <- matrix(.C("gehan_s_wt", as.double(beta), as.double(Y), as.double(X), as.integer(clsize),
                         as.double(sigma), as.integer(n), as.integer(p), as.integer(N), as.double(weights),
                         as.double(a), PACKAGE = "aftgee")[[10]], ncol = 1)
    }
    if (smooth == FALSE) {
        out <- matrix(.C("gehan_ns_wt", as.double(beta), as.double(Y), as.double(X), as.integer(clsize),
                         as.integer(n), as.integer(p), as.integer(N), as.double(weights),
                         as.double(a), PACKAGE = "aftgee")[[9]], ncol = 1) 
    }

    out
}

getPw <- function(Y, X, beta, N, delta, weights, rankWeights) {
    if (rankWeights == "logrank") {
        pw <- rep(1, nrow(X))
    }
    if (rankWeights == "PW") {
        pw <- getSuv(Y, X, beta, N, delta, weights)$Shat 
    }
    if (rankWeights == "GP") {
        pw <- getSuv(Y, X, beta, N, delta, weights)$Shat ^ (1 / ncol(X))
    }
    if (rankWeights == "eGP") {
        pw <- getSuv(Y, X, beta, N, delta, weights)
        pw <- ((pw$Shat - pw$Shatlast) / (1 - pw$Shatlast))^ (1 / ncol(X))
    }
    pw
}

getGw <- function(Y, X, beta, N, delta, clsize, sigma, weights, rankWeights, pw = NULL) {
    de <- getGehan(Y, X, beta, N, delta, clsize, sigma, weights)
    if (is.numeric(pw)) {
        gw <- pw / de
    }
    if (!is.numeric(pw)) {
        if (rankWeights == "gehan") {
            gw <- rep(1, nrow(X))
        }
        if (rankWeights == "logrank") {
            gw <- 1 / de
        }
        if (rankWeights == "PW") {
            ne <- getSuv(Y, X, beta, N, delta, weights)$Shat
            gw <- ne / de
        }
        if (rankWeights == "GP") {
            ne <- getSuv(Y, X, beta, N, delta, weights)$Shat ^ (1 / ncol(X))
            gw <- ne / de
        }
    }
    gw
}

getSmoothSuv <- function(Y, X, beta, N, delta, weights) {
    en <- Y - X %*% beta
    ik <- rep(1:N, each=N)
    jl <- rep(1:N, N)
    Shat <- fhat <- NULL
    dummy <- 1:N
    ord <- order(en)
    rij <- sqrt(diag((X[ord,] - X[dummy,]) %*% t((X[ord,] - X[dummy,]))))
    ei <- en[ord]
    weightsi <- weights[ord]
    deltai <- delta[ord]
    repeats <- table(ei)
    ## Shati <- survfit(Surv(ei, deltai) ~ 1, weights = weightsi)$surv
    ## assume no ties
    di <- rev(cumsum(rev(rep(1, N))))
    Z <- pnorm((en - ei) / rij)
    Z <- ifelse(is.na(Z) == T, 0, Z)
    hazi <- cumsum(deltai * Z / di)
    Shati <- exp(-1 * hazi)
    fhati <- diff(c(Shati[1], Shati, Shati[length(Shati)]), 2) / diff(c(ei[1], ei, ei[length(ei)]), 2)
    Shati <- rep(Shati, repeats)
    fhati <- rep(fhati, repeats)
    Shat[dummy[ord]] <- Shati
    fhat[dummy[ord]] <- fhati
    ## fhat <- ifelse(fhat < -1, -1, fhat)
    mu <- mean(subset(fhat, fhat > -1))
    fhat <- ifelse(fhat < -1, mu, fhat)
    list(Shat = Shat, fhat = fhat)
}

viEmp <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), B = 500,
                  mb = TRUE, zbeta = FALSE, smooth = TRUE,
                  rankWeights = "gehan", gw = gw,
                  sigma = diag(ncol(X)), gpweight = 1){
  p <- ncol(X)
  clsize <- unlist(lapply(split(id, id), length))
  ## n <- sum(unlist(lapply(split(weights, id), unique)))
  ## N <- sum(weights)
  n <- length(clsize)
  N <- length(weights)
  UnV <- zmat <- matrix(0, ncol = B, nrow = p)
  for (i in 1:B) {
    if ( mb == TRUE) {
      Z <- rep(rexp(length(clsize)), clsize)
    }
    if ( mb != TRUE) {
      Z <- rep(1, N)
    }
    if (zbeta == TRUE) {
      zb <- rnorm(p)
      newbeta <- beta + n ^ (-0.5) * zb
      zmat[,i] <- zb
    }
    if (zbeta != TRUE) {
      newbeta <- beta
    }
    sn <- vector("double", p)
    ## if (rankWeights %in% c("nsPW", "nsGP")) {
    ##     smooth <- FALSE
    ## }
    gpweight <- ifelse(rankWeights %in% c("nsGP", "GP"), gpweight, 1)
    if (smooth == TRUE) {
        if (rankWeights == "gehan") {
            UnV[,i] <- as.vector(.C("gehan_s_est", as.double(newbeta), as.double(Y),
                                    as.double(X), as.double(delta),
                                    as.integer(clsize), as.double(sigma),
                                    as.integer(length(clsize)), as.integer(p),
                                    as.integer(sum(clsize)),
                                    as.double(weights * Z), as.double(gw),
                                    as.double(sn), PACKAGE = "aftgee")[[12]]) # / n
        }
        if (rankWeights %in% c("nslogrank", "logrank")) {
            UnV[,i] <- uilogFun(beta = newbeta, Y = Y, X = X, delta = delta,
                                clsize = clsize, sigma = sigma, n = length(clsize),
                                Z = Z, weights = weights, constant = 0,
                                pw = rep(1, sum(clsize)),
                                rankWeights = "logrank", rkmethod = "sm") # / n
        }
        if (rankWeights %in% c("mPW", "mGP", "mlogrank", "userdefined")) {
            if (rankWeights == "mGP") {
                gw <- getGw(Y = Y, X = X, beta = beta, N = n, delta = delta,
                             clsize = clsize, sigma = sigma, weights = weights,
                             rankWeights = "GP")
            }
            if (rankWeights == "mlogrank") {
                gw <- getGw(Y = Y, X = X, beta = beta, N = n, delta = delta,
                             clsize = clsize, sigma = sigma, weights = weights,
                             rankWeights = "logrank")
            }
            if (rankWeights == "mPW") {
                gw <- getGw(Y = Y, X = X, beta = beta, N = n, delta = delta,
                             clsize = clsize, sigma = sigma, weights = weights,
                             rankWeights = "PW")
            }
            UnV[,i] <- as.vector(.C("gehan_s_est", as.double(newbeta), as.double(Y),
                                    as.double(X), as.double(delta),
                                    as.integer(clsize), as.double(sigma),
                                    as.integer(length(clsize)), as.integer(p),
                                    as.integer(sum(clsize)),
                                    as.double(weights * Z), as.double(gw),
                                    as.double(sn), PACKAGE = "aftgee")[[12]]) # / n
        }
        if (rankWeights %in% c("PW", "GP", "eGP")) {
            pw <- getPw(Y = Y, X = X, beta = newbeta, N = nrow(X),
                         delta = delta, weights = weights, rankWeights)
            UnV[,i] <- uilogFun(beta = newbeta, Y = Y, X = X, delta = delta,
                                clsize = clsize, sigma = sigma,
                                n = nrow(X), Z = Z, pw =pw,
                                weights = weights, constant = 0, rkmethod = "sm",
                                rankWeights = rankWeights) # / n
        }
    }
    if (smooth == FALSE) {
        if (rankWeights == "gehan") {
            UnV[,i] <- as.vector(.C("gehan_ns_est", as.double(newbeta), as.double(Y),
                                    as.double(X), as.double(delta),
                                    as.integer(clsize),
                                    as.integer(length(clsize)), as.integer(p),
                                    as.integer(sum(clsize)), 
                                    as.double(weights * Z), as.double(gw),
                                    as.double(sn), PACKAGE = "aftgee")[[11]]) # / n ## n and N?
        }
        if (rankWeights %in% c("logrank", "nslogrank")) {
            UnV[,i] <- uilogFun(beta = newbeta, Y = Y, X = X, delta = delta,
                                clsize = clsize, sigma = sigma,
                                n = length(clsize), Z = Z, weights = weights, constant = 0,
                                pw = rep(1, sum(clsize)),
                                rankWeights = "logrank", rkmethod = "nonsm") # / n
        }
        if (rankWeights %in% c("nsPW", "nsGP", "mPW", "mGP", "PW", "GP", "eGP")) {
            UnV[,i] <- uilogFun(beta = newbeta, Y = Y, X = X, delta = delta,
                                clsize = clsize, sigma = sigma,
                                n = length(clsize), Z = Z, weights = weights,
                                constant = 0, rkmethod = "nonsm",
                                pw = rep(1, sum(clsize)),
                                rankWeights = rankWeights) # / n
        }
    }
}
  vi <- var(t(UnV))
  list(vi = vi, zmat = zmat, UnV = UnV)
}


getSi <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)),
                   rankWeights = "gehan") {
    p <- ncol(X)
    clsize <- unlist(lapply(split(id, id), length))
    N <- sum(clsize)
    en <- Y - X %*% beta
    ik <- rep(1:N, each=N)
    jl <- rep(1:N, N)
    Shat <- NULL
    dummy <- 1:N
    ord <- order(en)
    ei <- en[ord]
    weightsi <- weights[ord]
    deltai <- delta[ord]
    repeats <- table(ei)
    Shati <- survfit(Surv(ei, deltai) ~ 1, weights = weightsi)$surv
    Shati <- rep(Shati, repeats)
    Shat[dummy[ord]] <- Shati
    xdif <- X[ik,] - X[jl,]
    edif <- en[ik] - en[jl]
    ind <- ifelse(edif <= 0, 1, 0)
    minEn <- ifelse(edif <= 0, ik, jl)
    si <- s <- NULL
    if (rankWeights == "gehan") {
        s <- weights[jl] * delta[ik] * ind * xdif + weights[jl] * log(Shat[minEn]) * xdif
        if (length(which(s == Inf)) > 0) {
            s <- ifelse(s == Inf, 0, s)
        }
        if (length(which(s == -Inf)) > 0) {
            s <- ifelse(s == -Inf, 0, s)
        }
        if (sum(is.na(s) > 0)) {
            s <- ifelse(is.na(s) == TRUE, 0, s)
        }
        s <- rowsum(s, ik)
        }
    if (rankWeights == "logrank") {
        haz <- -1 * log(Shat)
        haz <- ifelse(haz == Inf, 0, haz)
        haz <- ifelse(haz == -Inf, 0, haz)
        haz <- ifelse(is.na(haz) == TRUE, 0, haz)
        gamma1 <- rowsum(ind * X[jl,] * weights[jl], ik)
        gamma0 <- as.numeric(rowsum(ind * weights[jl], ik))
        si1i <-  si1 <- s2 <- matrix(0, nrow = N, ncol = p)
        si1 <- (X[1:N,] - gamma1 / gamma0)
        si1i <- si1[ord,]
        si1dif <- rbind(rep(0,p), diff(si1i))
        s2s <- apply(si1dif * haz[ord], 2, cumsum)
        s2[dummy[ord],] <- s2s
        s <- delta[1:N] * si1 - s2
    }
    s
}

getVi <- function(s, id, delta, weights, n) {
    clweights <- as.numeric(unlist(lapply(split(weights, id), unique)))
    s1 <- rowsum(s, group = id)
    si1 <- lapply(split(s1, 1:nrow(s1)), function(x) x %o% x)
    s11 <- mapply("*", si1, clweights, SIMPLIFY = FALSE)
    v1 <- Reduce("+", s11) / n
    v2 <- apply(s1 * clweights , 2, sum) / n
    v2 <- v2 %o% v2
    if (length(unique(weights)) == 1) {
        p <-  unique(weights) - 1
        vi <- p * (v1 - v2)
    }
    if (length(unique(weights)) > 1) {
        cweights <- unique(weights)
        cweights <- cweights[cweights != 1 & cweights != 0]
        ## vi <- (cweights - 1)* (v1 - v2)
        vi <- v1 - v2
    }
    list(vi = vi, v1 = v1, v2 = v2)
}


viClo <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), B = 500,
                  gw = gw, rankWeights = "gehan", stratify = TRUE) {
    s <- v1 <- vi <- v2i <- NULL
    n <- sum(unlist(lapply(split(weights, id), unique)))
    p <- ncol(X)
    clsize <- unlist(lapply(split(id, id), length))
    clid <- unlist(lapply(clsize, function(x) 1:x))
    stra <- match(weights, unique(weights))
    dim <- unique(clsize)
    s <- getSi(beta = beta, Y = Y, delta = delta, X = X, id = id,
                weights = weights, rankWeights = rankWeights)
    if (stratify) {
        v1 <- getVi(s, id, delta, weights, n)$v1
        vi <- v1
        v2 <- matrix(0, ncol = p, nrow = p)
        if (length(unique(stra)) > 1) {
            for (i in 1:length(unique(stra))) {
                ns <- sum(unlist(lapply(split(weights[stra == i], id[stra == i]), unique)))
                ## weights at cluster level
                v2i <- getVi(s = s[stra == i, ], id = id[stra == i], delta = delta[stra == i],
                              weights = weights[stra == i],
                              n = ns
                              )$vi
                strPr <- ns / n
                v2 <- v2 + v2i * strPr
            }
            vi <- v1 + v2
        }
    }
    if (!(stratify)) {
        v1 <- getVi(s, id, delta, weights, n)$v1 # / sum(clsize)
        vi <- v1
        if (length(unique(stra)) > 1) {
            v2 <- getVi(s, id, delta, weights * (1 - delta), n)$vi # / sum(clsize)
            vi <- v1 + v2
        }
    }
    list(vi = as.matrix(vi), s = s)
}


huangFun <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), B = 500,
                     vClose = TRUE, rankWeights = "gehan", gw = gw,
                     sigma = diag(ncol(X)), stratify = TRUE) {
  p <- ncol(X)
  clsize <- unlist(lapply(split(id, id), length))
  n <- sum(unlist(lapply(split(weights, id), unique)))
  N <- sum(weights)
  ## n <- length(clsize)
  ## N <- sum(clsize)
  betaTemp <- NULL
  UnMatV <- NULL
  if (vClose == TRUE) {
      vi <- viClo(beta, Y, delta, X, id, weights, B, rankWeights, stratify = TRUE)$vi * n
  }
  if (vClose != TRUE) {
      vi <- viEmp(beta, Y, delta, X, id, weights, B, mb = TRUE, zbeta = FALSE,
                  smooth = FALSE, rankWeights = rankWeights, gw = gw)$vi
  }
  qq <- chol(vi)
  qq <- t(qq)
  newBeta <- NULL
  Z <- rep(1, sum(clsize)) #
  newBeta <- matrix(0, ncol = p, nrow = p)
  for ( i in 1:p) {
      if (rankWeights == "gehan") {
          bb <- BBsolve(beta, uiFun, Y = Y, X = X, delta = delta, clsize = clsize,
                        sigma = sigma, n = length(clsize), Z = Z, gw = gw,
                        weights = weights, smooth = TRUE, constant = qq[,i] * n ^ -0.5,
                        quiet = TRUE)
          newBeta[, i] <- bb$par
      }
      if (rankWeights == "logrank") {
          bb <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize,
                        sigma = sigma, n = length(clsize), Z = Z,
                        weights = weights, smooth = TRUE, constant = qq[,i] * n ^ -0.5,
                        gpweight = 0, pw = rep(1, sum(clsize)),
                        rankWeights = rankWeights, quiet = TRUE)
          newBeta[, i] <- bb$par
      }
      if (rankWeights %in% c("PW", "Prentice-Wilcoxon", "GP")) {
          gpweight <- ifelse(rankWeights %in% c("GP", "nsGP"), gpweight, 1)
          bb <- BBsolve(beta, uilogFun, Y = Y, X = X, delta = delta, clsize = clsize,
                        sigma = sigma, n = length(clsize), Z = Z,
                        weights = weights, smooth = TRUE, constant = qq[,i] * n ^ -0.5,
                        gpweight = gpweight, pw = rep(1, sum(clsize)),
                        rankWeights = rankWeights, quiet = TRUE)
          newBeta[, i] <- bb$par
      }
  }
  dd <- newBeta - beta
  covmat <- n * t(dd) %*% dd
  list(covmat = covmat)
}

zlFun <- function(beta, Y, X, delta, id, weights = rep(1, nrow(X)),
                  B = 500, vClose = FALSE, rankWeights = "gehan", method = "sm",
                  gw = gw, stratify = TRUE,
                  sigma = diag(ncol(X))) {
    gpweight <- ifelse(rankWeights == "GP", 1/ncol(X), 1)
    rankWeights <- getRankName(rankWeights, method)$rankWeights
    p <- ncol(X)
    clsize <- unlist(lapply(split(id, id), length))
    n <- sum(unlist(lapply(split(weights, id), unique)))
    smooth <- ifelse(method != "nonsm", TRUE, FALSE)
    UnMat <- zmat <- ahat <- unTime <- NULL
    An.inv <- 1
    
    UnV <- viEmp(beta, Y, delta, X, id, weights, B, mb = FALSE, zbeta = TRUE, smooth = smooth,
                 gw = gw, rankWeights = rankWeights, sigma = sigma,
                 gpweight = gpweight)
    
    zmat <- UnV$zmat
    UnV <- UnV$UnV
    if (vClose == TRUE) {
        vi <- viClo(beta, Y, delta, X, id, weights, B, rankWeights, stratify = TRUE)$vi * n
    }
    if (vClose != TRUE) {
        vi <- viEmp(beta, Y, delta, X, id, weights, B, mb = TRUE, zbeta = FALSE, smooth = smooth,
                    gw = gw, rankWeights = rankWeights, gpweight = gpweight)$vi ## smooth
    }
    An <- matrix(0, nrow = p, ncol = p)
    for (i in 1:p) {
        An[i,] <- lm(UnV[i,] ~ matrix(t(zmat), ncol = p) - 1)$coef
    }
    if (qr(An)$rank != p) {
        covmat <- ginv(An) %*% vi %*% t(ginv(An))
        An.msg <- "An is singular"
        An.inv <- 0
    }
    if (qr(An)$rank == p) {
        covmat <- solve(An) %*% vi %*% solve(An)
        An.msg <- "An is nonsingular"
        An.inv <- 1
  }
    covmat <- covmat  / n
    covmat <- matrix(as.numeric(covmat), p)
    list(covmat = covmat, vi = vi, An.msg = An.msg, An.inv = An.inv)
}


isFun <- function(beta, Y, delta, X, id, weights = rep(1, nrow(X)), sigma, B = 500,
                  vClose = FALSE, gw = rep(1, nrow(X)), rankWeights = "gehan",
                  omega = FALSE, stratify = TRUE) {
    p <- ncol(X)
    clsize <- unlist(lapply(split(id, id), length))
    n <- sum(unlist(lapply(split(weights, id), unique)))
    ## n <- length(clsize)
    UnMat <- zmat <- ahat <- NULL
    An.inv <- 1
    if (omega == TRUE && rankWeights == "gehan") {
        vi <- omegaFun(beta, Y, X, delta, clsize, weights) * n## / n ^ 2
    }
    if (omega != TRUE && vClose == TRUE && rankWeights == "gehan") {
        vi <- viClo(beta, Y, delta, X, id, weights, B, gw = gw, stratify = TRUE)$vi * n
    }
    if (omega != TRUE && vClose == TRUE && rankWeights == "logrank") {
        vi <- viClo(beta, Y, delta, X, id, weights, B, "logrank", stratify = TRUE)$vi * n
    }
    if (omega != TRUE && vClose != TRUE) {
        vi <- viEmp(beta, Y, delta, X, id, weights, B, mb = TRUE, zbeta = FALSE,
                    gw = gw, rankWeights = rankWeights)$vi
    }
    if (rankWeights %in% c("gehan", "mPW", "mGP", "mlogrank")) {
        An <- abargehanfun(beta, Y, X, delta, clsize, sigma, weights, gw)
    }
    if (rankWeights == "logrank") {
        An <- abarlogfun(beta, Y, X, delta, clsize, sigma, weights)
    }
    if (rankWeights %in% c("PW", "GP", "Prentice-Wilcoxon")) {
        pw <- getSmoothSuv(Y, X, beta, n, delta, weights)
        An <- abarpwfun(beta, Y, X, delta, clsize, sigma, weights, pw)
    }

    if (qr(An)$rank != p) {
        covmat <- ginv(An) %*% vi %*% ginv(An)
        An.msg <- "An is singular"
        An.inv <- 0
    }
    if (qr(An)$rank == p) {
        covmat <- solve(An) %*% vi %*% solve(An)
        An.msg <- "An is nonsingular"
        An.inv <- 1
    }
    covmat <- matrix(as.numeric(covmat), p)
    list(covmat = covmat, vi = vi, An.msg = An.msg, An.inv = An.inv)
}

uilogFun <- function(beta, Y, X, delta, clsize, sigma, n, Z,
                     weights, constant = 0, pw = rep(1, nrow(X)), rankWeights, rkmethod) {
    N <- sum(clsize)
    p <- ncol(X)
    sn <- vector("double", p)
    ans <- numeric(p)
    n <- length(clsize)
    if (rkmethod != "nonsm") {
        ans <- .C("log_s_est", as.double(beta), as.double(Y), as.double(X),
                  as.double(delta), as.integer(clsize),
                  as.double(sigma), as.integer(n), as.integer(p),
                  as.integer(N), as.double(Z * weights),
                  as.double(pw), as.double(sn), PACKAGE = "aftgee")[[12]]
    }
    if (rkmethod == "nonsm") {
        pw <- getPw(Y = Y, X = X, beta = beta, N = N, delta = delta,
                    weights = weights, rankWeights)
        ans <- .C("log_ns_est", as.double(beta), as.double(Y), as.double(X),
                  as.double(delta), as.integer(clsize),
                  as.integer(n), as.integer(p),
                  as.integer(N), as.double(Z * weights),
                  as.double(pw), as.double(sn), PACKAGE = "aftgee")[[11]]
    }
    ans - constant
}

uiFun <- function(beta, Y, X, delta, clsize, sigma, n, Z, weights, smooth = TRUE, gw, constant = 0) {
  N <- nrow(X)
  p <- ncol(X)
  ans <- numeric(p)
  sn <- vector("double", p)
  ans <- .C("gehan_s_est", as.double(beta), as.double(Y), as.double(X), as.double(delta),
            as.integer(clsize), as.double(sigma), as.integer(n), as.integer(p),
            as.integer(N), as.double(Z * weights), as.double(gw),
            as.double(sn), PACKAGE = "aftgee")[[12]]
  ans <- ans - constant
}

