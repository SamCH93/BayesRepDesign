#' @title Compute probability of replication success
#'
#' @description This function computes the probabiliy of replication success
#'     based on a success region for the replication effect estimate, a design
#'     prior, and a replication standard error. If the specified number of sites
#'     is larger than 1, the supplied success region has to be formulated in
#'     terms of the meta-analytic replication effect estimate across sites.
#'
#' @param sregion Success region for replication effect estimate
#' @param dprior Design prior object
#' @param sr Standard error of replication effect estimate
#' @param nsites Number of sites, defaults to \code{1}. The sites are assumed to
#'     have the same standard error sr
#'
#' @return The probability of replication success
#'
#' @references
#'
#' Pawel, S., Consonni, G., and Held, L. (2022). Bayesian approaches to
#' designing replication studies. arXiv preprint.
#' \doi{10.48550/arXiv.2211.02552}
#'
#' @author Samuel Pawel
#'
#' @examples
#' dprior <- designPrior(to = 1.1, so = 1)
#' sregion <- successRegion(intervals = cbind(1.96, Inf))
#' pors(sregion = sregion, dprior = dprior, sr = 1)
#'
#' @export
pors <- function(sregion, dprior, sr, nsites = 1) {
    ## input checks
    stopifnot(
        class(sregion) == "successRegion",

        class(dprior) == "designPrior",

        length(sr) > 0,
        is.numeric(sr),
        all(is.finite(sr)),
        all(0 <= sr),

        length(nsites) == 1,
        is.numeric(nsites),
        is.finite(nsites),
        nsites > 0
    )

    ps <- vapply(X = sr, FUN = function(sr1) {
        ## compute parameters of predictive distribution of (average)
        ## replication effect estimate
        predmean <- dprior$dpMean
        predsd <- sqrt(dprior$dpVar + (dprior$tau^2 + sr1^2)/nsites)

        ## compute probability of replication success
        p <- sum(stats::pnorm(q = sregion[,2], mean = predmean, sd = predsd) -
                 stats::pnorm(q = sregion[,1], mean = predmean, sd = predsd))
    }, FUN.VALUE = 1)
    return(ps)
}
