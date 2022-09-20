#' @title Sample size determination for replication success based on
#'     effect size equivalence
#'
#' @description This function computes the standard error required to achieve
#'     replication success with a certain probability and based on establishing
#'     effect size equivalence of original and replication effect size
#'
#' @details The probability of replication success .....
#'
#' @md
#'
#' @param level 1 - confidence level of confidence interval for
#' @param dprior Design prior object
#' @param power Desired probability of replication success
#' @param margin The equivalence margin for the equivalence region [-margin,
#'     margin] of the difference in effect size of original and replication
#'     study
#' @param searchInt Interval for numerical search over replication standard
#'     errors
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
#' to1 <- 0.2
#' so1 <- 0.05
#' dprior <- designPrior(to = to1, so = so1, tau = 0.05)
#' ssdEqu(level = 0.05, dprior = dprior, power = 0.7, margin = 0.2)
#'
#' @export

ssdEqu <- function(level, dprior, power, margin, searchInt = c(0, 2)) {
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

        length(margin) == 1,
        is.numeric(margin),
        is.finite(margin),
        margin > 0,

        length(searchInt) == 2,
        is.numeric(searchInt),
        all(is.finite(searchInt)),
        0 <= searchInt[1], searchInt[1] < searchInt[2]
    )

    ## computing standard normal quantiles for power calculation
    za <- stats::qnorm(p = 1 - level/2)

    ## extracting design prior parameters
    tau <- dprior$tau
    dpmean <- dprior$dpMean
    dpvar <- dprior$dpVar
    to <- dprior$to
    so <- dprior$so

    ## probability of replication success for equivalence
    prosEqu <- function(sr) {
        sdiff <- sqrt(so^2 + sr^2)
        stats::pnorm(q = (to + margin - za*sdiff - dpmean)/
                         sqrt(sr^2 + tau^2 + dpvar)) -
            stats::pnorm(q = (to - margin + za*sdiff - dpmean)/
                             sqrt(sr^2 + tau^2 + dpvar))
    }

    ## computing bound of margin
    marginLim <- za*so
    if (margin <= marginLim) {
         warning(paste0("Equivalence not achievable with specified margin (at least ",
                           round(marginLim, 3), ")"))
    } else {

        ## computing bound of probability of replication success
        limP <- prosEqu(sr = 0)
        if (power > limP) {
            warning(paste0("Power not achievable with specified design prior (at most ",
                           round(limP, 3), ")"))
            sr <- NaN
            outPow <- NaN
        } else {
            ## computing replication standard error sr
            rootFun <- function(sr) {
                prosEqu(sr = sr) - power
            }
            res <- try(stats::uniroot(f = rootFun, interval = searchInt)$root)
            if (inherits(res, "try-error")) {
                sr <- NaN
                outPow <- NaN
                warning("Numerical problems, try adjusting searchInt")
            } else {
                sr <- res
                ## computing probability of replication success
                outPow <- prosEqu(sr = sr)
            }
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
## ciDiff <- function(to, tr, so, sr, alpha = 0.05) {
##     za <- stats::qnorm(p = 1 - alpha)
##     sdiff <- sqrt(so^2 + sr^2)
##     ci <- tr - to + c(-1, 1)*sdiff*za
##     return(ci)
## }
## to <- 0.2
## so <- 0.04
## sr <- 0.06
## sdiff <- sqrt(so^2 + sr^2)
## delta <- 0.02
## za <- stats::qnorm(p = 1 - 0.05)
## tr <- to - delta + za*sdiff
## ciDiff(to = to, tr = tr, so = so, sr = sr)
