## R codes in JSS aftgee; http://www.jstatsoft.org/v61/i11/}
## Codes with pkg(lss) and pkg(rms) are ignored

library(aftgee)
library(survival)

datgen <- function(n = 500, tau = .7) {
    x1 <- rbinom(n, 1, 0.5)
    x2 <- rnorm(n)
    e <- rweibull(n, 1, 3)
    T <- exp(2 + x1 + x2 + e)
    cstime <- runif(n, 0, tau)
    delta <- (T < cstime) * 1
    Y <- pmin(T, cstime)
    out <- data.frame(T = T, Y = Y, delta = delta, x1 = x1, x2 = x2)
}

set.seed(1)
mydata <- datgen()

## V1.0.0
## system.time(rk.srrMB <- aftsrr(Surv(Y, delta) ~ x1 + x2, data = mydata, variance = "MB"))
## system.time(rk.srrISMB <- aftsrr(Surv(Y, delta) ~ x1 + x2, data = mydata, variance = "ISMB"))

## Github version
system.time(rk.srrMB <- aftsrr(Surv(Y, delta) ~ x1 + x2, data = mydata, se = "MB"))
system.time(rk.srrISMB <- aftsrr(Surv(Y, delta) ~ x1 + x2, data = mydata, se = "ISMB"))

system.time(ls.sur <- survreg(Surv(Y, delta) ~ x1 + x2, data = mydata, dist = "lognormal"))
system.time(ls.gee <- aftgee(Surv(Y, delta) ~ x1 + x2, data = mydata))

##################################################################################
## National Wilms' tumor study
##################################################################################
data("nwtco", package = "survival")
nwtco$age <- nwtco$age/12
head(nwtco, 5)
set.seed(1)
## V1.0.0
## system.time(fit.IS <- aftsrr(Surv(edrel, rel) ~ histol + age,
##                              data = nwtco, variance = c("ISCF", "ISMB")))

system.time(fit.IS <- aftsrr(Surv(edrel, rel) ~ histol + age, data = nwtco, se = c("ISCF", "ISMB")))
summary(fit.IS)

table(nwtco$in.subcohort, nwtco$rel)
nwtco$in.casecohort <- (nwtco$in.subcohort | nwtco$rel == 1)
nwtco$hi <- 0
nwtco$hi <- ifelse(nwtco$in.casecohort & nwtco$rel == 1, 1, nwtco$hi)
nwtco$hi <- ifelse(nwtco$in.casecohort & nwtco$rel == 0, 5.93, nwtco$hi)
table(nwtco$hi)

## V1.0.0
## system.time(fit.gh <- aftsrr(Surv(edrel, rel) ~ histol + age,
##                              weights = hi, data = nwtco, variance = "ZLMB",
##                              subset = in.casecohort))
## system.time(fit.lk <- aftsrr(Surv(edrel, rel) ~ histol + age,
##                              weights = hi, data = nwtco, variance = "ZLMB", rankWeights = "logrank",
##                              subset = in.casecohort))
## system.time(fit.pw <- aftsrr(Surv(edrel, rel) ~ histol + age,
##                              weights = hi, data = nwtco, variance = "ZLMB", rankWeights = "PW",
##                              method = "monosm", subset = in.casecohort))

system.time(fit.gh <- aftsrr(Surv(edrel, rel) ~ histol + age,
                             weights = hi, data = nwtco, se = "ZLMB", subset = in.casecohort))
system.time(fit.lk <- aftsrr(Surv(edrel, rel) ~ histol + age,
                             weights = hi, data = nwtco, se = "ZLMB", rankWeights = "logrank",
                             subset = in.casecohort))
system.time(fit.pw <- aftsrr(Surv(edrel, rel) ~ histol + age,
                             weights = hi, data = nwtco, se = "ZLMB", rankWeights = "PW",
                             eqType = "mis", subset = in.casecohort))
summary(fit.gh)
summary(fit.lk)
summary(fit.pw) ## different results b/c different tol


summary(aftsrr(Surv(edrel, rel) ~ histol + age, weights = hi, data = nwtco,
               rankWeights = "GP", eqType = "mis", subset = in.casecohort))

summary(aftsrr(Surv(edrel, rel) ~ histol + age, weights = hi, data = nwtco,
               rankWeights = "PW", eqType = "mis", subset = in.casecohort))

##################################################################################
## Kidney catheter data
##################################################################################

data("kidney", package = "survival")
set.seed(123)
fit.ind <- aftgee(Surv(time, status) ~ age + sex, id = id, data = kidney)
fit.ex <- aftgee(Surv(time, status) ~ age + sex, id = id, data = kidney, corstr = "ex")
summary(fit.ind)
summary(fit.ex)
