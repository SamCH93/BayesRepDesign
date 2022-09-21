#' @title Sample size determination for replication success based on
#'     effect size equivalence
#'
#' @description This function computes the standard error required to achieve
#'     replication success with a certain probability and based on effect size
#'     equivalence of original and replication effect size
#'
#' @details The probability of replication success .....
#'
#' @md
#'
#' @param level 1 - confidence level of confidence interval for effect size
#'     difference
#' @param dprior Design prior object
#' @param power Desired probability of replication success
#' @param margin The equivalence margin > 0 for the equivalence region around
#'     zero of the difference in effect size of original and replication study
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
#' ssdEqu(level = 0.1, dprior = dprior, power = 0.8, margin = 0.2)
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

    ## computing bound of margin
    za <- stats::qnorm(p = 1 - level/2)
    so <- dprior$so
    marginLim <- za*so
    if (margin <= marginLim) {
         warning(paste0("Equivalence not achievable with specified margin (at least ",
                        round(marginLim, 3), ")"))
         sr <- NaN
         outPow <- NaN
    } else {

        ## computing bound of probability of replication success
        limP <- porsEqu(level = level, dprior = dprior, margin = margin, sr = 0)
        if (power > limP) {
            warning(paste0("Power not achievable with specified design prior (at most ",
                           round(limP, 3), ")"))
            sr <- NaN
            outPow <- NaN
        } else {
            ## computing replication standard error sr
            rootFun1 <- function(sr) {
                porsEqu(level = level, dprior = dprior, margin = margin,
                        sr = sr) - power
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
                outPow <- porsEqu(level = level, dprior = dprior, margin = margin, sr = sr)
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


#' @title Probability of replication success based on effect size equivalence
#'
#' @description This function computes the probability to achieve replication
#'     success on equivalence of original and replication effect size.
#'
#' @param level 1 - confidence level of confidence interval for effect size
#'     difference
#' @param dprior Design prior object
#' @param margin The equivalence margin > 0 for the equivalence region around
#'     zero of the difference in effect size of original and replication study
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
#' so1 <- 0.05
#' dprior <- designPrior(to = to1, so = so1, tau = 0.1)
#' porsEqu(level = 0.1, dprior = dprior, margin = 0.3, sr = 0.05)
#'
#' @export

porsEqu <- function(level, dprior, margin, sr) {
    ## input checks
    stopifnot(
        length(level) == 1,
        is.numeric(level),
        is.finite(level),
        0 < level, level < 1,

        class(dprior) == "designPrior",

        length(margin) == 1,
        is.numeric(margin),
        is.finite(margin),
        margin > 0,

        length(sr) == 1,
        is.numeric(sr),
        is.finite(sr),
        0 <= sr
    )

    ## compute probability of replication success
    to <- dprior$to
    so <- dprior$so
    sdiff <- sqrt(so^2 + sr^2)
    za <- stats::qnorm(p = 1 - level)
    if (margin <= za*sdiff) {
        p <- 0
    } else {
        sregion <- successRegion(intervals = cbind(to - margin + za*sdiff,
                                                   to + margin - za*sdiff))
        p <- pors(sregion = sregion, dprior = dprior, sr = sr)
    }
    return(p)
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