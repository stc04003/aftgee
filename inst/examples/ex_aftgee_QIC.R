datgen <- function(n = 100, tau = 0.3, cen = 100, dim = 2) {
    id <- rep(1:n, rep(dim, n))
    x1 <- rbinom(dim * n, 1, 0.5)
    x2 <- rnorm(dim * n)
    e <- c(t(exp(MASS::mvrnorm(n = n, mu = rep(0, dim), Sigma = tau + (1 - tau) * diag(dim)))))
    T <- exp(2 + x1 + x2 + e)
    cstime <- runif(n, 0, cen)
    delta <- (T < cstime) * 1
    Y <- pmin(T, cstime)
    out <- data.frame(T = T, Y = Y, delta = delta, x1 = x1, x2 = x2, id = rep(1:n, each = dim))
    out
}

set.seed(1)
mydata <- datgen(n = 50, dim = 2)
fit1 <- aftgee(Surv(Y, delta) ~ x1 + x2, data = mydata, id = id, corstr = "ind", B = 0)
fit2 <- aftgee(Surv(Y, delta) ~ x1 + x2, data = mydata, id = id, corstr = "ex", B = 0)

QIC(fit1)
QIC(fit2)
