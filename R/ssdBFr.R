#' @title Sample size determination for replication success based on
#'     replication Bayes factor
#'
#' @description This function computes the standard error required to achieve
#'     replication success with a certain probability and based on the
#'     replication Bayes factor under normality
#'
#' @param level Bayes factor level below which replication success is achieved
#' @param dprior Design prior object
#' @param power Desired probability of replication success
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
#' dprior <- designPrior(to = to1, so = so1, tau = 0.03)
#' ssdBFr(level = 1/10, dprior = dprior, power = 0.8)
#'
#' @export

ssdBFr <- function(level, dprior, power,
                   searchInt = c(.Machine$double.eps^0.5, 2)) {
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
    limP <- porsBFr(level = level, dprior = dprior, sr = .Machine$double.eps)
    if (power > limP) {
        warning(paste0("Power not achievable with specified design prior (at most ",
                       round(limP, 3), ")"))
        sr <- NaN
        outPow <- NaN
    } else {
        ## computing replication standard error sr
        rootFun <- function(sr) {
            porsBFr(level = level, dprior = dprior, sr = sr) - power
        }
        res <- try(stats::uniroot(f = rootFun, interval = searchInt)$root)
        if (inherits(res, "try-error")) {
            sr <- NaN
            outPow <- NaN
            warning("Numerical problems, try adjusting searchInt")
        } else {
            sr <- res
            ## computing probability of replication success
            outPow <- porsBFr(level = level, dprior = dprior, sr = sr)
        }
    }

    ## create output object
    out <- list("designPrior" = dprior, "power" = power,
                "powerRecomputed" = outPow, "sr" = sr,
                "c" = dprior$so^2/sr^2,
                type = paste("replication Bayes factor <=", signif(level, 3),
                             "(numerical computation)"))
    class(out) <- "ssdRS"
    return(out)
}


#' @title Probability of replication success based on replication Bayes factor
#'
#' @description This function computes the probability to achieve replication
#'     success based on the replication Bayes factor
#'
#' @param level Bayes factor level below which replication success is achieved
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
#' to1 <- 0.2
#' so1 <- 0.05
#' dprior <- designPrior(to = to1, so = so1, tau = 0.03)
#' porsBFr(level = 1/10, dprior = dprior, sr = c(0.05, 0.04))
#'
#' @export

porsBFr <- function(level, dprior, sr) {
    ## input checks
    stopifnot(
        length(level) == 1,
        is.numeric(level),
        is.finite(level),
        0 < level,

        class(dprior) == "designPrior",

        length(sr) > 0,
        is.numeric(sr),
        all(is.finite(sr)),
        all(0 <= sr)
    )

    ps <- vapply(X = sr, FUN = function(sr1) {
        ## compute probability of replication success
        to <- dprior$to
        so <- dprior$so
        c <- so^2/sr1^2
        A <- sr1^2*(1 +1/c)*(to^2/so^2 - 2*log(level) + log(1 + c))
        sregion <- successRegion(intervals = rbind(c(-Inf, -sqrt(A) - to/c),
                                                   c(sqrt(A) - to/c, Inf)))
        p <- pors(sregion = sregion, dprior = dprior, sr = sr1)
        return(p)
    }, FUN.VALUE = 1)
    return(ps)
}


## ## ## checking some stuff
## ## BFr <- function(to, tr, so, sr) {
## ##     stats::dnorm(x = tr, sd = sr) /
## ##         stats::dnorm(x = tr, mean = to, sd = sqrt(so^2 + sr^2))
## ## }
## to <- 0.2
## so <- 0.04
## sr <- 0.06
## zo <- to/so
## c <- so^2/sr^2
## gamma <- 1/10
## A <- sr^2*(1 + sr^2/so^2)*(to^2/so^2 - 2*log(gamma) + log(1 + so^2/sr^2))
## tr <- sqrt(A) - to*sr^2/so^2
## BFr(to, tr, so, sr)
## sRegionBFr <- function(sr) {
##     A <- sr^2*(1 + sr^2/so^2)*(to^2/so^2 - 2*log(gamma) + log(1 + so^2/sr^2))
##     successRegion(intervals = rbind(c(-Inf, -sqrt(A) - to*sr^2/so^2),
##                                     c(sqrt(A) - to*sr^2/so^2, Inf)))
## }
## sRegionBFr(0.05)
## ssdRS(sregionfun = sRegionBFr, dprior = designPrior(to = to, so = so),
##       power = 0.8)
