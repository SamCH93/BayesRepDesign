#' @title Sample size determination for replication success based on
#'     significance
#'
#' @description This function computes the standard error required to achieve
#'     replication success with a certain probability and based on statistical
#'     significance of the replication effect estimate.
#'
#' @param level Significance level for the replication effect estimate
#'     (one-sided and in the same direction as the original effect estimate)
#' @param dprior Design prior object
#' @param power Desired probability of replication success
#'
#' @return Returns an object of class \code{"ssdRS"}. See \code{\link{ssd}} for
#'     details.
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
#' ## specify design prior
#' to1 <- 2
#' so1 <- 0.5
#' dprior <- designPrior(to = to1, so = so1, tau = 0.1)
#' ssdSig(level = 0.025, dprior = dprior, power = 0.9)
#'
#' @export

ssdSig <- function(level, dprior, power) {
    ## input checks
    stopifnot(
        length(level) == 1,
        is.numeric(level),
        is.finite(level),
        0 < level, level < 1,

        class(dprior) == "designPrior",

        length(power) == 1,
        is.numeric(power),
        is.finite(power),
        0 < power, power < 1,
        level < power
    )

    ## extracting design prior parameters
    tau <- dprior$tau
    dpmean <- dprior$dpMean
    dpvar <- dprior$dpVar
    so <- dprior$so
    to <- dprior$to

    ## computing standard normal quantiles for power calculation
    if (sign(to) > 0) {
        za <- stats::qnorm(p = 1 - level)
    } else {
        za <- stats::qnorm(p = level)
    }
    zb <- stats::qnorm(p = power)

    ## computing bound of probability of replication success
    limP <- porsSig(level = level, dprior = dprior, sr = 0)
    if (power > limP) {
        warning(paste0("Power not achievable with specified design prior (at most ",
                       round(limP, 3), ")"))
        sr <- NaN
        outPow <- NaN
    } else {
        ## computing replication standard error sr analytically
        sr <- (dpmean*za - zb*sqrt((za^2 - zb^2)*(tau^2 + dpvar) +
                               dpmean^2))/((za^2 - zb^2))

        ## computing probability of replication success
        outPow <- porsSig(level = level, dprior = dprior, sr = sr)
    }

    ## create output object
    out <- list("designPrior" = dprior, "power" = power,
                "powerRecomputed" = outPow, "sr" = sr,
                "c" = so^2/sr^2,
                type = paste("replication p-value <=", signif(level, 3),
                             "(exact computation)"))
    class(out) <- "ssdRS"
    return(out)
}


#' @title Probability of replication success based on significance
#'
#' @description This function computes the probability to achieve replication
#'     success on statistical significance of the replication effect estimate.
#'
#' @param level Significance level for p-value of the replication effect
#'     estimate (one-sided and in the same direction as the original effect
#'     estimate)
#' @param dprior Design prior object
#' @param sr Replication standard error
#'
#' @return The probability to achieve replication success
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
#' ## specify design prior
#' to1 <- 2
#' so1 <- 1
#' dprior <- designPrior(to = to1, so = so1, tau = 0.1)
#' porsSig(level = 0.025, dprior = dprior, sr = c(0.5, 0.3))
#'
#' @export

porsSig <- function(level, dprior, sr) {
    ## input checks
    stopifnot(
        length(level) == 1,
        is.numeric(level),
        is.finite(level),
        0 < level, level < 1,

        class(dprior) == "designPrior",

        length(sr) > 0,
        is.numeric(sr),
        all(is.finite(sr)),
        all(0 <= sr)
    )

    ps <- vapply(X = sr, FUN = function(sr1) {
        ## compute probability of replication success
        sregion <- successRegionSig(sr = sr1, to = dprior$to, tau = 0,
                                    nsites = 1, level = level)
        p <- pors(sregion = sregion, dprior = dprior, sr = sr1)
        return(p)
    }, FUN.VALUE = 1)
    return(ps)
}


#' @title Success region based on significance
#'
#' @description This function returns the success region for the (meta-analytic)
#'     replication effect estimate to achieve significance
#'
#' @param sr Replication standard error
#' @param to Original effect estimate
#' @param tau Heterogeneity standard deviation used in the calculation of the
#'     meta-analytic replication effect estimate and its standard error.
#'     Defaults to \code{0} (fixed effects analysis)
#' @param nsites nsites Number of sites, defaults to \code{1}. The effect
#'     estimates from all sites are assumed to have the same standard error
#'     \code{sr}
#' @param level Significance level for p-value of the (average) replication
#'     effect estimate (one-sided and in the same direction as the original
#'     effect estimate)
#'
#' @return An object of class \code{"successRegion"}. See
#'     \code{\link{successRegion}} for details.
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
#' successRegionSig(sr = 0.05, to = 0.2, tau = 0.01, nsites = 3, level = 0.025)
#'
#' @export

successRegionSig <- function(sr, to, tau = 0, nsites = 1, level) {
    ## input checks
    stopifnot(
        length(sr) == 1,
        is.numeric(sr),
        is.finite(sr),
        0 <= sr,

        length(to) == 1,
        is.numeric(sr),
        is.finite(sr),

        length(tau) == 1,
        is.numeric(tau),
        is.finite(tau),
        0 <= tau,

        length(nsites) == 1,
        is.numeric(nsites),
        is.finite(nsites),
        nsites > 0,

        length(level) == 1,
        is.numeric(level),
        is.finite(level),
        0 < level, level < 1
    )

    ## compute standard error of weighted average
    srMA <- 1/sqrt(nsites/(sr^2 + tau^2))

    ## success region depends on direction of original estimate
    if (sign(to) >= 0) {
        sregion <- successRegion(intervals = cbind(stats::qnorm(p = 1 - level)*srMA, Inf))
    } else {
        sregion <- successRegion(cbind(-Inf, stats::qnorm(p = level)*srMA))
    }
    return(sregion)
}
