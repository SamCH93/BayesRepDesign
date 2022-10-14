#' @title Sample size determination for replication success based on
#'     the sceptical Bayes factor
#'
#' @description This function computes the standard error required to achieve
#'     replication success with a certain probability and based on the
#'     Bayes factor
#'
#' @param level Threshold for the sceptical Bayes factor below which
#'     replication success is achieved
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
#' ssdBFs(level = 1/10, dprior = dprior, power = 0.95)
#'
#' @export

ssdBFs <- function(level, dprior, power,
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
    limP <- porsBFs(level = level, dprior = dprior, sr = .Machine$double.eps^0.5)
    if (power > limP) {
        warning(paste0("Power not achievable with specified design prior (at most ",
                       round(limP, 3), ")"))
        sr <- NaN
        outPow <- NaN
    } else {
       ## computing replication standard error sr
        rootFun <- function(sr) {
            porsBFs(level = level, dprior = dprior, sr = sr) - power
        }
        res <- try(stats::uniroot(f = rootFun, interval = searchInt)$root)
        if (inherits(res, "try-error")) {
            sr <- NaN
            outPow <- NaN
            warning("Numerical problems, try adjusting searchInt")
        } else {
            sr <- res
            ## computing probability of replication success
            outPow <- porsBFs(level = level, dprior = dprior, sr = sr)
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
#'     success based on the Bayes factor
#'
#' @param level Threshold for the sceptical Bayes factor below which replication
#'     success is achieved
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
#' dprior <- designPrior(to = to1, so = so1)
#' porsBFs(level = 1/38.7, dprior = dprior, sr = 0.05)
#'
#' @export

porsBFs <- function(level, dprior, sr) {
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
        s <- BayesRep::vss(x = zo, gamma = level) # relative prior variance
        if (is.nan(s)) {
            p <- 0
        } else {
            c <- so^2/sr1^2
            A <- (-2*log(level) + 2*log(1 + c) - 2*log(1 + s*c) +
                  zo^2/(1 - s))*(1/c + s)*(sr1^2 + so^2)/(1 - s)
            M <- to*(1/c + s)/(1 - s)
            if (s < 1) {
                ints <- rbind(c(-Inf, -sqrt(A) - M),
                    c(sqrt(A) - M, Inf))
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
