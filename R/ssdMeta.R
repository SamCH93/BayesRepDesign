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

    ## computing standard normal quantiles for power calculation
    za <- stats::qnorm(p = 1 - level)

    ## extracting design prior parameters
    tau <- dprior$tau
    dpmean <- dprior$dpMean
    dpvar <- dprior$dpVar
    to <- dprior$to
    so <- dprior$so

    ## probability of replication success for meta-analysis
    ## TODO implement for negative effect estimates
    prosMeta <- function(sr) {
        stats::pnorm(q = (dpmean - sr*za*sqrt(1 + sr^2/so^2) + to*sr^2/so^2)/
                         sqrt(sr^2 + tau^2 + dpvar))
    }

    ## computing bound of probability of replication success
    limP <- prosMeta(sr = 0)
    if (power > limP) {
        warning(paste0("Power not achievable with specified design prior (at most ",
                       round(limP, 3), ")"))
        sr <- NaN
        outPow <- NaN
    } else {
        ## computing replication standard error sr
        rootFun <- function(sr) {
            prosMeta(sr = sr) - power
        }
        res <- try(stats::uniroot(f = rootFun, interval = searchInt)$root)
        if (inherits(res, "try-error")) {
            sr <- NaN
            outPow <- NaN
            warning("Numerical problems, try adjusting searchInt")
        } else {
            sr <- res
            ## computing probability of replication success
            outPow <- prosMeta(sr = sr)
        }

    }

    ## create output object
    out <- list("designPrior" = dprior, "power" = power,
                "powerRecomputed" = outPow, "sr" = sr,
                "c" = so^2/sr^2)
    class(out) <- "ssdRS"
    return(out)
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
