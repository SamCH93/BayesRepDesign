#' @title Sample size determination for replication success based on
#'     effect size equivalence
#'
#' @description This function computes the standard error required to achieve
#'     replication success with a certain probability and based on effect size
#'     equivalence of original and replication effect size. Effect size
#'     equivalence is defined by the confidence interval for the difference
#'     between the original and replication effect sizes falling within an
#'     equivalence region around zero defined by the specified margin.
#'
#' @param level 1 - confidence level of confidence interval for effect size
#'     difference
#' @param dprior Design prior object
#' @param power Desired probability of replication success
#' @param margin The equivalence margin > 0 for the symmetric equivalence region
#'     around zero
#' @param searchInt Interval for numerical search over replication standard
#'     errors
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
#' Anderson, S. F. and Maxwell, S. E. (2016). There's more than one way to
#' conduct a replication study: Beyond statistical significance. Psychological
#' Methods, 21(1), 1-12. \doi{10.1037/met0000051}
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
            rootFun <- function(sr) {
                porsEqu(level = level, dprior = dprior, margin = margin,
                        sr = sr) - power
            }
            res <- try(stats::uniroot(f = rootFun, interval = searchInt)$root, silent = TRUE)
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
                "c" = so^2/sr^2,
                type = paste("equivalence with confidence level =", signif(1 - level, 3),
                             "and margin =", signif(margin, 3), "(numerical computation)"))
    class(out) <- "ssdRS"
    return(out)
}


#' @title Probability of replication success based on effect size equivalence
#'
#' @description This function computes the probability to achieve replication
#'     success on equivalence of original and replication effect size. Effect
#'     size equivalence is defined by the confidence interval for the difference
#'     between the original and replication effect sizes falling within an
#'     equivalence region around zero defined by the specified margin.
#'
#' @param level 1 - confidence level of confidence interval for effect size
#'     difference
#' @param dprior Design prior object
#' @param margin The equivalence margin > 0 for the symmetric equivalence region
#'     around zero
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
#' Anderson, S. F. and Maxwell, S. E. (2016). There's more than one way to
#' conduct a replication study: Beyond statistical significance. Psychological
#' Methods, 21(1), 1-12. \doi{10.1037/met0000051}
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## specify design prior
#' to1 <- 2
#' so1 <- 0.05
#' dprior <- designPrior(to = to1, so = so1, tau = 0.1)
#' porsEqu(level = 0.1, dprior = dprior, margin = 0.3, sr = c(0.05, 0.03))
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

        length(sr) > 0,
        is.numeric(sr),
        all(is.finite(sr)),
        all(0 <= sr)
    )

    ps <- vapply(X = sr, FUN = function(sr1) {
        ## compute probability of replication success
        to <- dprior$to
        so <- dprior$so
        sdiff <- sqrt(so^2 + sr1^2)
        za <- stats::qnorm(p = 1 - level/2)
        if (margin <= za*sdiff) {
            p <- 0
        } else {
            sregion <- successRegion(intervals = cbind(to - margin + za*sdiff,
                                                       to + margin - za*sdiff))
            p <- pors(sregion = sregion, dprior = dprior, sr = sr1)
        }
        return(p)}, FUN.VALUE = 1)
    return(ps)
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
