context("ssd")

library(BayesRepDesign)

## testing grid
tos <- c(-3, -2, 2, 3)
sos <- c(0.5, 1, 1.1)
taus <- c(0, 0.25)
mus <- c(0, 0.5)
sps <- c(Inf, 10)
alpha <- 0.025
pows <- c(0.4, 0.6, 0.7, 0.8)
testGrid <- expand.grid(to = tos, so = sos, tau = taus, mu = mus, sp = sps,
                        power = pows)

## SSD numerically and analytically
suppressWarnings({
    resultsDF <- do.call("rbind", lapply(X = seq(1, nrow(testGrid)), FUN = function(i) {
        to <- testGrid$to[i]
        so <- testGrid$so[i]
        tau <- testGrid$tau[i]
        mu <- testGrid$mu[i]
        sp <- testGrid$sp[i]
        power <- testGrid$power[i]
        dp <- designPrior(to = to, so = so, mu = mu, sp = sp, tau = tau)
        ## compute numerically with ssd
        if (sign(to) >= 0) {
            ## function to compute required standard error for significance one-sided
            sregionfunSig <- function(sr) {
                successRegion(intervals = cbind(stats::qnorm(p = 1- alpha)*sr, Inf))
            }
        } else {
            ## function to compute required standard error for significance one-sided
            sregionfunSig <- function(sr) {
                successRegion(intervals = cbind(-Inf, stats::qnorm(p = alpha)*sr))
            }
        }
        res_ssd <- ssd(sregionfun = sregionfunSig, dprior = dp, power = power,
                       searchInt = c(.Machine$double.eps, 4))

        ## compute analytically with formulae
        res_ssdSig <- ssdSig(level = alpha, dprior = dp, power = power)
        ## meta-analysis
        res_ssdMeta <- ssdMeta(level = alpha^2, dprior = dp, power = power,
                               searchInt = c(0, 10))
        ## equivalence
        res_ssdEqu <- ssdEqu(level = alpha*4, dprior = dp, power = power,
                             margin = 5, searchInt = c(0, 2))
        ## sceptical p-value
        res_ssdPs <- ssdPs(level = alpha*6, dprior = dp, power = power)
        ## replication BF
        res_ssdBFr <- ssdBFr(level = alpha*4, dprior = dp, power = power)
        ## sceptical BF
        res_ssdBFs <- ssdBFs(level = alpha*4, dprior = dp, power = power)

        ## return everything
        out <- data.frame(to = to, so = so, tau = tau, mu = mu, sp = sp,
                          power = power, power_ssd = res_ssd$powerRecomputed,
                          sr_ssd = res_ssd$sr, power_ssdSig = res_ssdSig$powerRecomputed,
                          sr_ssdSig = res_ssdSig$sr,
                          power_ssdPs = res_ssdPs$powerRecomputed,
                          sr_ssdPs = res_ssdPs$sr,
                          sr_ssdMeta = res_ssdMeta$sr,
                          power_ssdMeta = res_ssdMeta$powerRecomputed,
                          sr_ssdEqu = res_ssdEqu$sr,
                          power_ssdEqu = res_ssdEqu$powerRecomputed,
                          sr_ssdBFr = res_ssdBFr$sr,
                          power_ssdBFr = res_ssdBFs$powerRecomputed,
                          sr_ssdBFs = res_ssdBFs$sr,
                          power_ssdBFs = res_ssdBFs$powerRecomputed
                          )
        return(out)
    }))
})

## tests
testthat::test_that("ssd: back-computed power = specified power", {
    ## check that the same power
    testthat::expect_equal(resultsDF$power_ssd, resultsDF$power,
                           tolerance = 0.001)
})

testthat::test_that("ssdSig: back-computed power = specified power", {
    ## check that the same power
    testthat::expect_equal(resultsDF$power_ssdSig, resultsDF$power,
                           tolerance = 0.001)
})

testthat::test_that("ssd and ssdSig: back-computed power is the same", {
    ## check that the same power
    testthat::expect_equal(resultsDF$power_ssd, resultsDF$power_ssdSig,
                           tolerance = 0.001)
})

testthat::test_that("ssd and ssdSig: replication standard error is the same", {
    ## check that the same standard error
    testthat::expect_equal(resultsDF$sr_ssd, resultsDF$sr_ssdSig,
                           tolerance = 0.0001)
})

testthat::test_that("ssdPs: back-computed power = specified power", {
    ## check that the same power
    testthat::expect_equal(resultsDF$power_ssdPs, resultsDF$power,
                           tolerance = 0.001)
})

testthat::test_that("ssdEqu: back-computed power = specified power", {
    ## check that the same power
    nonNA <- !is.na(resultsDF$power_ssdEqu)
    testthat::expect_equal(resultsDF$power_ssdEqu[nonNA], resultsDF$power[nonNA],
                           tolerance = 0.001)
})

testthat::test_that("ssdMeta: back-computed power = specified power", {
    ## check that the same power
    nonNA <- !is.na(resultsDF$power_ssdMeta)
    testthat::expect_equal(resultsDF$power_ssdMeta[nonNA], resultsDF$power[nonNA],
                           tolerance = 0.001)
})

testthat::test_that("ssdBFr: back-computed power = specified power", {
    ## check that the same power
    nonNA <- !is.na(resultsDF$power_ssdBFr)
    testthat::expect_equal(resultsDF$power_ssdBFr[nonNA], resultsDF$power[nonNA],
                           tolerance = 0.001)
})

testthat::test_that("ssdBFs: back-computed power = specified power", {
    ## check that the same power
    nonNA <- !is.na(resultsDF$power_ssdBFs)
    testthat::expect_equal(resultsDF$power_ssdBFs[nonNA], resultsDF$power[nonNA],
                           tolerance = 0.001)
})
