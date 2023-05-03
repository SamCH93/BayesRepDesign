#' @title Sample size determination for replication success based on
#'     replication Bayes factor
#'
#' @description This function computes the standard error required to achieve
#'     replication success with a certain probability and based on the
#'     replication Bayes factor under normality. The replication Bayes factor is
#'     assumed to be oriented so that values below one indicate replication
#'     success, whereas values above one indicate evidence for the null
#'     hypothesis.
#'
#' @param level Bayes factor level below which replication success is achieved
#' @param dprior Design prior object
#' @param power Desired probability of replication success
#' @param searchInt Interval for numerical search over replication standard
#'     errors
#' @param paradox Should the probability of replication success be computed
#'     allowing for the replication paradox (replication success when the effect
#'     estimates from original and replication study have a different sign)?
#'     Defaults to \code{TRUE}
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
#' Verhagen, J. and Wagenmakers, E. J. (2014). Bayesian tests to quantify the result
#' of a replication attempt. Journal of Experimental Psychology: General,
#' 145:1457-1475. \doi{10.1037/a0036731}
#'
#' Ly, A., Etz, A., Marsman, M., and Wagenmakers, E.-J. (2018). Replication Bayes
#' factors from evidence updating. Behavior Research Methods, 51(6), 2498-2508.
#' \doi{10.3758/s13428-018-1092-x}
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
                   searchInt = c(.Machine$double.eps^0.5, 2),
                   paradox = TRUE) {
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
        0 <= searchInt[1], searchInt[1] < searchInt[2],

        length(paradox) == 1,
        is.logical(paradox),
        !is.na(paradox)
    )

    ## computing bound of probability of replication success
    limP <- porsBFr(level = level, dprior = dprior, sr = .Machine$double.eps,
                    paradox = paradox)
    if (power > limP) {
        warning(paste0("Power not achievable with specified design prior (at most ",
                       round(limP, 3), ")"))
        sr <- NaN
        outPow <- NaN
    } else {
        ## computing replication standard error sr
        rootFun <- function(sr) {
            porsBFr(level = level, dprior = dprior, sr = sr,
                    paradox = paradox) - power
        }
        res <- try(stats::uniroot(f = rootFun, interval = searchInt)$root, silent = TRUE)
        if (inherits(res, "try-error")) {
            sr <- NaN
            outPow <- NaN
            warning("Numerical problems, try adjusting searchInt")
        } else {
            sr <- res
            ## computing probability of replication success
            outPow <- porsBFr(level = level, dprior = dprior, sr = sr,
                              paradox = paradox)
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
#'     success based on the replication Bayes factor. The replication Bayes
#'     factor is assumed to be oriented so that values below one indicate
#'     replication success, whereas values above one indicate evidence for the
#'     null hypothesis.
#'
#' @param level Bayes factor level below which replication success is achieved
#' @param dprior Design prior object
#' @param sr Replication standard error
#' @param paradox Should the probability of replication success be computed
#'     allowing for the replication paradox (replication success when the effect
#'     estimates from original and replication study have a different sign)?
#'     Defaults to \code{TRUE}
#'
#' @return The probability to achieve replication success
#'
#' @references
#'
#' Pawel, S., Consonni, G., and Held, L. (2022). Bayesian approaches to
#' designing replication studies. arXiv preprint.
#' \doi{10.48550/arXiv.2211.02552}
#'
#' Verhagen, J. and Wagenmakers, E. J. (2014). Bayesian tests to quantify the result
#' of a replication attempt. Journal of Experimental Psychology: General,
#' 145:1457-1475. \doi{10.1037/a0036731}
#'
#' Ly, A., Etz, A., Marsman, M., & Wagenmakers, E.-J. (2018). Replication Bayes
#' factors from evidence updating. Behavior Research Methods, 51(6), 2498-2508.
#' \doi{10.3758/s13428-018-1092-x}
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

porsBFr <- function(level, dprior, sr, paradox = TRUE) {
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
        all(0 <= sr),

        length(paradox) == 1,
        is.logical(paradox),
        !is.na(paradox)
    )

    ps <- vapply(X = sr, FUN = function(sr1) {
        ## compute probability of replication success
        to <- dprior$to
        so <- dprior$so
        c <- so^2/sr1^2
        A <- sr1^2*(1 +1/c)*(to^2/so^2 - 2*log(level) + log(1 + c))
        if (paradox) {
            ## success region that allows for replication success with
            ## replication paradox
            sregion <- successRegion(intervals = rbind(c(-Inf, -sqrt(A) - to/c),
                                                       c(sqrt(A) - to/c, Inf)))
        } else {
            ## success region depends on sign direction of original estiamte
            if (sign(to) == 1) {
                sregion <- successRegion(intervals = cbind(sqrt(A) - to/c, Inf))
            } else {
                 sregion <- successRegion(intervals = cbind(-Inf, -sqrt(A) - to/c))
            }
        }
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
