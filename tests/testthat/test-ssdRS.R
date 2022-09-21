context("ssdRS")

library(BayesRepDesign)

## function to compute required standard error for significance one-sided
sregionfunSig <- function(sr, alpha = 0.025) {
    successRegion(intervals = cbind(stats::qnorm(p = 1- alpha)*sr, Inf))
}

## design prior testing grid
tos <- c(2, 3, 4)
sos <- c(0.5, 1, 2)
taus <- c(0, 0.25)
mus <- c(0, 0.5)
sps <- c(Inf, 10)

## power testing grid
pows <- c(0.6, 0.7, 0.8)

## SSD numerically and analytically
testGrid <- expand.grid(to = tos, so = sos, tau = taus, mu = mus, sp = sps,
                        power = pows)
resultsDF <- do.call("rbind", lapply(X = seq(1, nrow(testGrid)), FUN = function(i) {
    to <- testGrid$to[i]
    so <- testGrid$so[i]
    tau <- testGrid$tau[i]
    mu <- testGrid$mu[i]
    sp <- testGrid$sp[i]
    power <- testGrid$power[i]
    dp <- designPrior(to = to, so = so, mu = mu, sp = sp, tau = tau)
    ## compute numerically with ssdRS
    ssd <- ssd(sregionfun = sregionfunSig, dprior = dp, power = power)
    ## compute analytically with formula
    ## TODO implement formula
    sr_analyt <- NaN
    ## return everything
    out <- data.frame(to = to, so = so, tau = tau, mu = mu, sp = sp,
                      power = power, power_ssdRS = ssd$power,
                      sr_ssdRS = ssd$sr, sr_analyt = sr_analyt)
    return(out)
}))

test_that("ssdRS: back-computed power is correct", {
    ## check that the same power
    expect_equal(resultsDF$power, resultsDF$power_ssdRS, tolerance = 0.001)
})

## test_that("Checking that probability of replication success correctly computed", {
##     ## check that the same standard error
##     expect_equal(resultsDF$sr_ssdRS, resultsDF$sr_analyt)
## })
