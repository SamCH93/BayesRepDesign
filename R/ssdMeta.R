#' @title Sample size determination for replication success based on
#'     meta-analytic significance
#'
#' @description This function computes the standard error required to achieve
#'     replication success with a certain probability and based on statistical
#'     significance of the fixed-effects meta-analytic effect estimate.
#'
#' @details The probability of replication success .....
#'
#' @md
#'
#' @param level Significance level (one-sided) for the replication effect
#'     estimate
#' @param dprior Design prior object
#' @param power Desired probability of replication success
#' @param searchInt Interval for numerical search over replication standard
#'     errors
#'
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
#' dprior <- designPrior(to = to1, so = so1, tau = 0.1, sp = 2)
#' ssdMeta(level = 0.025^2, dprior = dprior, power = 0.9)
#'
#' @export

ssdMeta <- function(level, dprior, power, searchInt = c(0, 2)) {
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
        level < power,

        length(searchInt) == 2,
        is.numeric(searchInt),
        all(is.finite(searchInt)),
        0 <= searchInt[1], searchInt[1] < searchInt[2]
    )

    ## computing bound of probability of replication success
    limP <- porsMeta(level = level, dprior = dprior, sr = 0)
    if (power > limP) {
        warning(paste0("Power not achievable with specified design prior (at most ",
                       round(limP, 3), ")"))
        sr <- NaN
        outPow <- NaN
    } else {
        ## computing replication standard error sr
        rootFun1 <- function(sr) {
            porsMeta(level = level, dprior = dprior, sr = sr) - power
        }
        rootFun <- Vectorize(FUN = rootFun1)
        res <- try(stats::uniroot(f = rootFun, interval = searchInt)$root)
        if (inherits(res, "try-error")) {
            sr <- NaN
            outPow <- NaN
            warning("Numerical problems, try adjusting searchInt")
        } else {
            sr <- res
            ## computing probability of replication success
            outPow <- porsMeta(level = level, dprior = dprior, sr = sr)
        }

    }

    ## create output object
    out <- list("designPrior" = dprior, "power" = power,
                "powerRecomputed" = outPow, "sr" = sr,
                "c" = dprior$so^2/sr^2)
    class(out) <- "ssdRS"
    return(out)
}


#' @title Probability of replication success based on meta-analytic significance
#'
#' @description This function computes the probability to achieve replication
#'     success on statistical significance of the (fixed effects) meta-analytic
#'     effect estimate
#'
#' @param level Significance level for p-value of the meta-analytic effect
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
#' porsMeta(level = 0.025^2, dprior = dprior, sr = 0.2)
#'
#' @export

porsMeta <- function(level, dprior, sr) {
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
    so <- dprior$so
    to <- dprior$to
    if (sign(dprior$to) >= 0) {
        lowerLim <- stats::qnorm(p = 1 - level)*sr*sqrt(1 + sr^2/so^2) -
            to*sr^2/so^2
        sregion <- successRegion(intervals = cbind(lowerLim, Inf))
    } else {
        upperLim <- stats::qnorm(p = level)*sr*sqrt(1 + sr^2/so^2) -
            to*sr^2/so^2
        sregion <- successRegion(cbind(-Inf, upperLim))
    }

    ## compute probability of replication success
    pors(sregion = sregion, dprior = dprior, sr = sr)
}


## ## checking some stuff
## pmeta <- function(to, tr, so, sr) {
##     sm <- 1/sqrt(1/so^2 + 1/sr^2)
##     tm <- sm^2*(to/so^2 + tr/sr^2)
##     return(stats::pnorm(q = tm/sm, lower = FALSE))
## }
## to <- 0.5
## so <- 0.8
## sr <- 1.1
## sm <- 1/sqrt(1/so^2 + 1/sr^2)
## za <- stats::qnorm(p = 0.975)
## tr <- sr^2*(za/sm - to/so^2)
## tr <- sr*za*sqrt(1 + sr^2/so^2) - to*sr^2/so^2
## pmeta(to = to, tr = tr, so = so, sr = sr)