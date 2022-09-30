#' @title Compute probability of replication success
#'
#' @description This function computes the probabiliy of replication success
#'     based on a success region for the replication effect estimate, a
#'     design prior, and a replication standard error.
#'
#'
#' @param sregion Success region for replication effect estimate
#' @param dprior Design prior object
#' @param sr Standard error of Replication effect estimate
#'
#' @return The probability of replication success
#'
#' @references
#'
#' Pawel, S., Consonni, G., and Held, L. (2022). Bayesian approaches to
#' designing replication studies. arXiv preprint.
#' \doi{10.48550/arXiv.XXXX.XXXXX}
#'
#' @author Samuel Pawel
#'
#' @examples
#' dprior <- designPrior(to = 1.1, so = 1)
#' sregion <- successRegion(intervals = cbind(1.96, Inf))
#' pors(sregion = sregion, dprior = dprior, sr = c(1, 0.9))
#'
#' @export
pors <- function(sregion, dprior, sr) {
    ## input checks
    stopifnot(
        class(sregion) == "successRegion",

        class(dprior) == "designPrior",

        length(sr) > 0,
        is.numeric(sr),
        all(is.finite(sr)),
        all(0 <= sr)
    )

    ps <- vapply(X = sr, FUN = function(sr1) {
        ## compute parameters of predictive distribution of replication effect
        ## estimate
        dpmean <- dprior$dpMean
        dpsd <- sqrt(dprior$dpVar + dprior$tau^2 + sr1^2)

        ## compute probability of replication success
        p <- sum(stats::pnorm(q = sregion[,2], mean = dpmean, sd = dpsd) -
                 stats::pnorm(q = sregion[,1], mean = dpmean, sd = dpsd))
    }, FUN.VALUE = 1)
    return(ps)
}
