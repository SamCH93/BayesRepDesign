#' @title Sample size determination for replication success based on
#'     significance
#'
#' @description This function computes the standard error required to achieve
#'     replication success with a certain probability and based on statistical
#'     significance of the replication effect estimate.
#'
#' @details The probability of replication success .....
#'
#' @md
#'
#' @param level Significance level (one-sided) for the replication effect
#'     estimate
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
#' so1 <- 1
#' dprior <- designPrior(to = to1, so = so1, tau = 0.1)
#' ssdSig(level = 0.025, dprior = dprior, power = 0.8)
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

    ## computing standard normal quantiles for power calculation
    za <- stats::qnorm(p = 1 - level)
    zb <- stats::qnorm(p = 1 - power)

    ## extracting design prior parameters
    tau <- dprior$tau
    dpmean <- dprior$dpMean
    dpvar <- dprior$dpVar
    so <- dprior$so

    ## computing bound of probability of replication success
    limP <- porsSig(level = level, dprior = dprior, sr = 0)
    if (power > limP) {
        warning(paste0("Power not achievable with specified design prior (at most ",
                       round(limP, 3), ")"))
        sr <- NaN
        outPow <- NaN
    } else {
        ## computing replication standard error sr analytically
        sr <- (dpmean*za + zb*sqrt((za^2 - zb^2)*(tau^2 + dpvar) +
                               dpmean^2))/((za^2 - zb^2))

        ## computing probability of replication success
        outPow <- porsSig(level = level, dprior = dprior, sr = sr)
    }

    ## create output object
    out <- list("designPrior" = dprior, "power" = power,
                "powerRecomputed" = outPow, "sr" = sr,
                "c" = so^2/sr^2)
    class(out) <- "ssdRS"
    return(out)
}


#' @title Probability of replication success based on significance
#'
#' @description This function computes the probability to achieve replication
#'     success on statistical significance of the replication effect estimate.
#'
#' @param level Significance level for p-value of the replication effect
#'     estimate (one-sided in the same direction as the original effect
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
#' \doi{10.48550/arXiv.XXXX.XXXXX}
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## specify design prior
#' to1 <- 2
#' so1 <- 1
#' dprior <- designPrior(to = to1, so = so1, tau = 0.1)
#' porsSig(level = 0.025, dprior = dprior, sr = 0.5)
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

        length(sr) == 1,
        is.numeric(sr),
        is.finite(sr),
        0 <= sr
    )

    ## success region depends on direction of original estimate
    if (sign(dprior$to) >= 0) {
        sregion <- successRegion(intervals = cbind(stats::qnorm(p = 1 - level)*sr, Inf))
    } else {
        sregion <- successRegion(cbind(-Inf, stats::qnorm(p = level)*sr))
    }

    ## compute probability of replication success
    pors(sregion = sregion, dprior = dprior, sr = sr)
}
