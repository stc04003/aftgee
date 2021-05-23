## Simulate data from an AFT model with possible depended response
datgen <- function(n = 100, tau = 0.3, dim = 2) {
    x1 <- rbinom(dim * n, 1, 0.5)
    x2 <- rnorm(dim * n)
    e <- c(t(exp(MASS::mvrnorm(n = n, mu = rep(0, dim), Sigma = tau + (1 - tau) * diag(dim)))))
    tt <- exp(2 + x1 + x2 + e)
    cen <- runif(n, 0, 100)
    data.frame(Time = pmin(tt, cen), status = 1 * (tt < cen),
               x1 = x1, x2 = x2, id = rep(1:n, each = dim))
}
set.seed(1); dat <- datgen(n = 50, dim = 2)
fm <- Surv(Time, status) ~ x1 + x2
summary(aftgee(fm, data = dat, id = id, corstr = "ind", B = 8))
summary(aftgee(fm, data = dat, id = id, corstr = "ex", B = 8))
