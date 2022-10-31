#' @title Sample size determination for replication success based on
#'     dual criterion
#'
#' @description This function computes the standard error required to achieve
#'     replication success with a certain probability and based on the dual
#'     criterion (significance and relevance).
#'
#' @param level Significance level (one-sided) for the replication effect
#'     estimate
#' @param relevance Minimally relevant effect size
#' @param dprior Design prior object
#' @param power Desired probability of replication success
#'
#' @return An ssdRS object containing the design prior, the recomputed power,
#'     the required replication standard error sr, and the relative sample size
#'     c = so^2/sr^2 = nr/no
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
#' ## specify design prior
#' to1 <- 2
#' so1 <- 0.5
#' dprior <- designPrior(to = to1, so = so1, tau = 0.1)
#' ssdDual(level = 0.025, dprior = dprior, power = 0.9)
#'
#' @export

ssdDual <- function(level, relevance, dprior, power) {
    ## ## input checks
    ## stopifnot(
    ##     length(level) == 1,
    ##     is.numeric(level),
    ##     is.finite(level),
    ##     0 < level, level < 1,

    ##     class(dprior) == "designPrior",

    ##     length(power) == 1,
    ##     is.numeric(power),
    ##     is.finite(power),
    ##     0 < power, power < 1,
    ##     level < power
    ## )

    ## ## extracting design prior parameters
    ## tau <- dprior$tau
    ## dpmean <- dprior$dpMean
    ## dpvar <- dprior$dpVar
    ## so <- dprior$so
    ## to <- dprior$to

    ## ## computing standard normal quantiles for power calculation
    ## if (sign(to) > 0) {
    ##     za <- stats::qnorm(p = 1 - level)
    ## } else {
    ##     za <- stats::qnorm(p = level)
    ## }
    ## zb <- stats::qnorm(p = power)

    ## ## computing bound of probability of replication success
    ## limP <- porsSig(level = level, dprior = dprior, sr = 0)
    ## if (power > limP) {
    ##     warning(paste0("Power not achievable with specified design prior (at most ",
    ##                    round(limP, 3), ")"))
    ##     sr <- NaN
    ##     outPow <- NaN
    ## } else {
    ##     ## computing replication standard error sr analytically
    ##     sr <- (dpmean*za - zb*sqrt((za^2 - zb^2)*(tau^2 + dpvar) +
    ##                            dpmean^2))/((za^2 - zb^2))

    ##     ## computing probability of replication success
    ##     outPow <- porsSig(level = level, dprior = dprior, sr = sr)
    ## }

    ## ## create output object
    ## out <- list("designPrior" = dprior, "power" = power,
    ##             "powerRecomputed" = outPow, "sr" = sr,
    ##             "c" = so^2/sr^2,
    ##             type = paste("replication p-value <=", signif(level, 3),
    ##                          "(exact computation)"))
    ## class(out) <- "ssdRS"
    ## return(out)
}


#' @title Probability of replication success based on dual criterion
#'
#' @description This function computes the probability to achieve replication
#'     success based on the dual criterion (significance and relevance).
#'
#' @param level Significance level for p-value of the replication effect
#'     estimate (one-sided in the same direction as the original effect
#'     estimate)
#' @param relevance Minimally relevant effect size
#' @param dprior Design prior object
#' @param sr Replication standard error
#'
#' @return The probability to achieve replication success
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
#' ## specify design prior
#' to1 <- 2
#' so1 <- 1
#' dprior <- designPrior(to = to1, so = so1, tau = 0.1)
#' porsDual(level = 0.025, dprior = dprior, sr = c(0.5, 0.3))
#'
#' @export

porsDual <- function(level, relevance, dprior, sr) {
    ## ## input checks
    ## stopifnot(
    ##     length(level) == 1,
    ##     is.numeric(level),
    ##     is.finite(level),
    ##     0 < level, level < 1,

    ##     class(dprior) == "designPrior",

    ##     length(sr) > 0,
    ##     is.numeric(sr),
    ##     all(is.finite(sr)),
    ##     all(0 <= sr)
    ## )

    ## ps <- vapply(X = sr, FUN = function(sr1) {
    ##     ## compute probability of replication success
    ##     sregion <- successRegionSig(sr = sr1, to = dprior$to, tau = 0,
    ##                                 nsites = 1, level = level)
    ##     p <- pors(sregion = sregion, dprior = dprior, sr = sr1)
    ##     return(p)
    ## }, FUN.VALUE = 1)
    ## return(ps)
}


#' @title Success region based on dual criteterion
#'
#' @description This function returns the success region for the replication
#'     effect estimate to achieve replication success based on the dual
#'     criterion (statistical significance and effect size relevance)
#'
#' @param sr Replication standard error
#' @param to Original effect estimate
#' @param relevance Minimally relevant effect size
#' @param level Significance level for p-value of the (average) replication
#'     effect estimate (one-sided in the same direction as the original effect
#'     estimate)
#'
#' @return A successRegion object
#'
#' @references
#'
#' Pawel, S., Consonni, G., and Held, L. (2022). Bayesian approaches to
#' designing replication studies. arXiv preprint.
#' \doi{10.48550/arXiv.XXXX.XXXXX}
#'
#' Rosenkranz, G. (2021). Replicability of studies following a dual-criterion
#' design. Statistics in Medicine, 40(18), 4068-4076. \doi{10.1002/sim.9014}
#'
#' @author Samuel Pawel
#'
#' @examples
#' successRegionDual(sr = 0.05, to = 0.2, relevance = 0.1, level = 0.025)
#'
#' @export

successRegionDual <- function(sr, to, relevance, level) {
    ## ## input checks
    ## stopifnot(
    ##     length(sr) == 1,
    ##     is.numeric(sr),
    ##     is.finite(sr),
    ##     0 <= sr,

    ##     length(to) == 1,
    ##     is.numeric(to),
    ##     is.finite(to),

    ##     length(relevance) == 1,
    ##     is.numeric(relevance),
    ##     is.finite(relevance),
    ##     sign(relevance) == sign(to),
    ##     abs(to) >= abs(relevance),

    ##     length(level) == 1,
    ##     is.numeric(level),
    ##     is.finite(level),
    ##     0 < level, level < 1
    ## )

    ## ## success region depends on direction of original estimate
    ## if (sign(to) >= 0) {
    ##     lower <- pmax(stats::qnorm(p = 1 - level)*srMA, relevance)
    ##     sregion <- successRegion(intervals = cbind(lower, Inf))
    ## } else {
    ##     upper <- pmin(stats::qnorm(p = level)*srMA, relevance)
    ##     sregion <- successRegion(cbind(-Inf, upper)
    ## }
    ## return(sregion)
}
