#' @title Sample size determination for replication success based on
#'     TOST equivalence
#'
#' @description This function computes the standard error required to achieve
#'     replication success with a certain probability and based on establishing
#'     the absence of a practically relevant effect size with the Two One-Sided
#'     Tests (TOST) procedure in the replication study.
#'
#' @param level Significance level for the TOST p-value
#' @param dprior Design prior object
#' @param power Desired probability of replication success
#' @param margin The equivalence margin > 0 for the equivalence region around
#'     zero that defines a region of practically irrelevant effect sizes
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
#' to1 <- 0.05
#' so1 <- 0.05
#' dprior <- designPrior(to = to1, so = so1, tau = 0.05)
#' ssdTOST(level = 0.05, dprior = dprior, power = 0.9, margin = 0.3)
#'
#' @export

ssdTOST <- function(level, dprior, power, margin, searchInt = c(0, 2)) {
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

    ## computing bound of probability of replication success
    limP <- porsTOST(level = level, dprior = dprior, margin = margin, sr = 0)
    if (power > limP) {
        warning(paste0("Power not achievable with specified design prior (at most ",
                       round(limP, 3), ")"))
        sr <- NaN
        outPow <- NaN
    } else {
        ## computing replication standard error sr
        rootFun <- function(sr) {
            porsTOST(level = level, dprior = dprior, margin = margin,
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
            outPow <- porsTOST(level = level, dprior = dprior, margin = margin, sr = sr)
        }
    }

    ## create output object
    out <- list("designPrior" = dprior, "power" = power,
                "powerRecomputed" = outPow, "sr" = sr,
                "c" = dprior$so^2/sr^2,
                type = paste("TOST equivalence with level =", signif(level, 3),
                             "and margin =", signif(margin, 3), "(numerical computation)"))
    class(out) <- "ssdRS"
    return(out)
}

#' @title Probability of replication success based on TOST equivalence
#'
#' @description This function computes the probability to achieve replication
#'     success based on establishing the absence of a practically relevant
#'     effect size with the Two One-Sided Tests (TOST) procedure in
#'     the replication study.
#'
#' @param level Significance level for the TOST p-value
#' @param dprior Design prior object
#' @param margin The equivalence margin > 0 for the equivalence region around
#'     zero that defines a region of practically irrelevant effect sizes
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
#' porsTOST(level = 0.1, dprior = dprior, margin = 0.3, sr = c(0.05, 0.03))
#'
#' @export

porsTOST <- function(level, dprior, margin, sr) {
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

    ## compute probability of replication success
    za <- stats::qnorm(p = 1 - level/2)
    ps <- vapply(X = sr, FUN = function(sr1) {
        ints <- cbind(-margin + za*sr1, margin - za*sr1)
        if (ints[1] >= ints[2]) {
            p <- 0
        } else {
            sregion <- successRegion(intervals = ints)
            p <- pors(sregion = sregion, dprior = dprior, sr = sr1)
        }
        return(p)}, FUN.VALUE = 1)
    return(ps)
}
## TODO add tests
