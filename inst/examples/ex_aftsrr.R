## Simulate data from an AFT model
datgen <- function(n = 100) {
    x1 <- rbinom(n, 1, 0.5)
    x2 <- rnorm(n)
    e <- rnorm(n)
    tt <- exp(2 + x1 + x2 + e)
    cen <- runif(n, 0, 100)
    data.frame(Time = pmin(tt, cen), status = 1 * (tt < cen),
               x1 = x1, x2 = x2, id = 1:n)
}
set.seed(1); dat <- datgen(n = 50)
summary(aftsrr(Surv(Time, status) ~ x1 + x2, data = dat, se = c("ISMB", "ZLMB"), B = 10))

## Data set with sampling weights
data(nwtco, package = "survival")
subinx <- sample(1:nrow(nwtco), 668, replace = FALSE)
nwtco$subcohort <- 0
nwtco$subcohort[subinx] <- 1
pn <- table(nwtco$subcohort)[[2]] / sum(table(nwtco$subcohort))
nwtco$hi <- nwtco$rel + ( 1 - nwtco$rel) * nwtco$subcohort / pn
nwtco$age12 <- nwtco$age / 12
nwtco$study <- nwtco$study - 3
nwtco$histol = nwtco$histol - 1
sub <- nwtco[subinx,]
fit <- aftsrr(Surv(edrel, rel) ~ histol + age12 + study, id = seqno,
              weights = hi, data = sub, B = 10, se = c("ISMB", "ZLMB"),
              subset = stage == 4)
summary(fit)
