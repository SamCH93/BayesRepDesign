library(BayesRepDesign)

## design prior testing grid
tos <- c(-3, 2, 0, 2, 3)
sos <- c(0.5, 1, 1.1)
taus <- c(0, 0.25)
mus <- c(0, 0.5)
sps <- c(Inf, 10, 0.5)
alpha <- 0.025
srs <- c(0.1, 0.4, 1, 2)
testGrid <- expand.grid(to = tos, so = sos, sr = srs, tau = taus, mu = mus,
                        sp = sps)

## compute probability of replication success
resultsDF <- do.call("rbind", lapply(X = seq(1, nrow(testGrid)), FUN = function(i) {
    to <- testGrid$to[i]
    so <- testGrid$so[i]
    sr <- testGrid$sr[i]
    tau <- testGrid$tau[i]
    mu <- testGrid$mu[i]
    sp <- testGrid$sp[i]
    dp <- designPrior(to = to, so = so, mu = mu, sp = sp, tau = tau)

    ## generic function with significance success region
    if (sign(to) >= 0) {
        sreg <- successRegion(intervals = cbind(stats::qnorm(p = 1- alpha)*sr, Inf))
    } else {
        sreg <- successRegion(intervals = cbind(-Inf, stats::qnorm(p = alpha)*sr))
    }
    p_pors <- pors(sregion = sreg, dprior = dp, sr = sr)

    ## significance
    p_porsSig <- porsSig(level = alpha, dprior = dp, sr = sr)

    ## TOST
    p_porsTOST <- porsTOST(level = alpha, dprior = dp, margin = 1, sr = sr)

    ## meta-analysis
    p_porsMeta <- porsMeta(level = alpha^2, dprior = dp, sr = sr)

    ## equivalence
    p_porsEqu <- porsEqu(level = 1 - 2*alpha, dprior = dp, margin = 1, sr = sr)

    ## sceptical p-value
    p_porsPs <- porsPs(level = alpha*6, dprior = dp, sr = sr)

    ## standard BF
    p_porsBF01 <- porsBF01(level = alpha*4, dprior = dp, sr = sr)

    ## replication BF
    p_porsBFr <- porsBFr(level = alpha*2, dprior = dp, sr = sr)

    ## sceptical BF
    p_porsBFs <- porsBFs(level = alpha*4, dprior = dp, sr = sr)

    ## return everything
    out <- data.frame(to, so, tau,  mu, sp, sr,
                      p_pors,
                      p_porsSig,
                      p_porsTOST,
                      p_porsMeta,
                      p_porsEqu,
                      p_porsPs,
                      p_porsBF01,
                      p_porsBFr,
                      p_porsBFs)
    return(out)
}))

## porsDF <- resultsDF
## save(object = porsDF, file = "porsDF.rda",  version = 2)

load("porsDF.rda")

expect_equal(porsDF$p_pors, resultsDF$p_pors, tolerance = 1e-05,
             info = "pors")

expect_equal(porsDF$p_porsSig, resultsDF$p_porsSig, tolerance = 1e-05,
             info = "porsSig")

expect_equal(porsDF$p_porsTOST, resultsDF$p_porsTOST, tolerance = 1e-05,
             info = "porsTOST")

expect_equal(porsDF$p_porsMeta, resultsDF$p_porsMeta, tolerance = 1e-05,
             info = "porsMeta")

expect_equal(porsDF$p_porsEqu, resultsDF$p_porsEqu, tolerance = 1e-05,
             info = "porsEqu")

expect_equal(porsDF$p_porsPs, resultsDF$p_porsPs, tolerance = 1e-05,
             info = "porsPs")

expect_equal(porsDF$p_porsBF01, resultsDF$p_porsBF01, tolerance = 1e-05,
             info = "porsBF01")

expect_equal(porsDF$p_porsBFr, resultsDF$p_porsBFr, tolerance = 1e-05,
             info = "porsBFr")

expect_equal(porsDF$p_porsBFs, resultsDF$p_porsBFs, tolerance = 1e-05,
             info = "porsBFs")
