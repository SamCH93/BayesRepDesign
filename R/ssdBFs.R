#' @title Sample size determination for replication success based on
#'     the sceptical Bayes factor
#'
#' @description This function computes the standard error required to achieve
#'     replication success with a certain probability and based on the sceptical
#'     Bayes factor. The sceptical Bayes factor is assumed to be oriented so
#'     that values below one indicate replication success.
#'
#' @param level Threshold for the sceptical Bayes factor below which replication
#'     success is achieved
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
#' Pawel, S. and Held, L. (2020). The sceptical Bayes factor for the assessement
#' of replication success. Journal of the Royal Statistical Society: Series B
#' (Statistical Methodology), 84(3), 879-911. \doi{10.1111/rssb.12491}
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## specify design prior
#' to1 <- 0.2
#' so1 <- 0.05
#' dprior <- designPrior(to = to1, so = so1, tau = 0.03)
#' ssdBFs(level = 1/10, dprior = dprior, power = 0.9)
#'
#' @export

ssdBFs <- function(level, dprior, power,
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
    limP <- porsBFs(level = level, dprior = dprior, sr = .Machine$double.eps^0.5,
                    paradox = paradox)
    if (power > limP) {
        warning(paste0("Power not achievable with specified design prior (at most ",
                       round(limP, 3), ")"))
        sr <- NaN
        outPow <- NaN
    } else {
       ## computing replication standard error sr
        rootFun <- function(sr) {
            porsBFs(level = level, dprior = dprior, sr = sr,
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
            outPow <- porsBFs(level = level, dprior = dprior, sr = sr,
                              paradox = paradox)
        }
    }

    ## create output object
    out <- list("designPrior" = dprior, "power" = power,
                "powerRecomputed" = outPow, "sr" = sr,
                "c" = dprior$so^2/sr^2,
                type = paste("sceptical Bayes factor <=", signif(level, 3),
                             "(numerical computation)"))
    class(out) <- "ssdRS"
    return(out)
}


#' @title Probability of replication success based on the sceptical Bayes factor
#'
#' @description This function computes the probability to achieve replication
#'     success based on the sceptical Bayes factor. The sceptical Bayes factor
#'     is assumed to be oriented so that values below one indicate replication
#'     success.
#'
#' @param level Threshold for the sceptical Bayes factor below which replication
#'     success is achieved
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
#' Pawel, S. and Held, L. (2020). The sceptical Bayes factor for the assessement
#' of replication success. Journal of the Royal Statistical Society: Series B
#' (Statistical Methodology), 84(3), 879-911. \doi{10.1111/rssb.12491}
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## specify design prior
#' to1 <- 0.2
#' so1 <- 0.05
#' dprior <- designPrior(to = to1, so = so1)
#' porsBFs(level = 1/3, dprior = dprior, sr = 0.05)
#'
#' @export

porsBFs <- function(level, dprior, sr, paradox = TRUE) {
    ## input checks
    stopifnot(
        length(level) == 1,
        is.numeric(level),
        is.finite(level),
        0 < level, level < 1,

        class(dprior) == "designPrior",

        length(sr) > 0,
        is.numeric(sr),
        all(is.finite(sr)),
        all(0 <= sr)
    )

    ps <- vapply(X = sr, FUN = function(sr1) {
        ## compute success region
        so <- dprior$so
        to <- dprior$to
        zo <- to/so
        q <- lamW::lambertWm1(x = -zo^2/level^2*exp(-zo^2))
        s <- -zo^2/q - 1
        if (is.nan(s) | s < 0) {
            p <- 0
        } else {
            c <- so^2/sr1^2
            A <- (-2*log(level) + 2*log(1 + c) - 2*log(1 + s*c) +
                  zo^2/(1 - s))*(1/c + s)*(sr1^2 + so^2)/(1 - s)
            M <- to*(1/c + s)/(1 - s)
            if (s < 1) {
                intsBothsides <- rbind(c(-Inf, -sqrt(A) - M), c(sqrt(A) - M, Inf))
                ## replication paradox can occur in this situation
                if (paradox) {
                    ints <- intsBothsides
                } else {
                    if (sign(to) > 0) {
                        ints <- intsBothsides[2,,drop = FALSE]
                    } else {
                        ints <- intsBothsides[1,,drop = FALSE]
                    }
                }
            } else if (isTRUE(all.equal(s, 1, tolerance = 0.0001))) {
                X <- 2*log(level)*(so^2 + sr^2)/to
                if (sign(to) > 0) {
                    ints <- cbind(to - X, Inf)
                } else {
                    ints <- cbind(-Inf, to + X)
                }
            } else {
                ints <- cbind(-sqrt(A) - M, sqrt(A) - M)
            }
            sregion <- successRegion(intervals = ints)
            p <- pors(sregion = sregion, dprior = dprior, sr = sr1)
        }
        return(p)
    }, FUN.VALUE = 1)
    return(ps)
}


## ## some checks
## so <- 1.5
## sr <- 0.8
## s <- c(0.8, 1.5)
## gamma <- 1/10
## to <- 2
## tr <- 1.5
## ## should be the same
## tr^2/(sr^2 + s*so^2) - (tr - to)^2/(so^2 + sr^2)
## so^2*(1 - s)/((sr^2 + s*so^2)*(sr^2 + so^2))*(tr + to*(sr^2 + s*so^2)/so^2/(1 - s))^2 +
##     to^2/so^2/(s - 1)
