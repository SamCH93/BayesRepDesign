#' @title Sample size determination for replication success based on
#'     Bayes factor
#'
#' @description This function computes the standard error required to achieve
#'     replication success with a certain probability and based on the Bayes
#'     factor under normality. The Bayes factor is oriented so that values above
#'     one indicate evidence for the null hypothesis of the effect size being
#'     zero, whereas values below one indicate evidence for the hypothesis of
#'     the effect size being non-zero (with normal prior assigned to it).
#'
#' @param level Bayes factor level below which replication success is achieved
#' @param dprior Design prior object
#' @param power Desired probability of replication success
#' @param priormean Mean of the normal prior under the alternative. Defaults to
#'     \code{0}
#' @param priorvar Variance of the normal prior under the alternative. Defaults
#'     to \code{1}
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
#' @author Samuel Pawel
#'
#' @examples
#' ## specify design prior
#' to1 <- 0.2
#' so1 <- 0.05
#' dprior <- designPrior(to = to1, so = so1, tau = 0.03)
#' ssdBF01(level = 1/10, dprior = dprior, power = 0.8)
#'
#' @export

ssdBF01 <- function(level, dprior, power, priormean = 0, priorvar = 1,
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

        length(priormean) == 1,
        is.numeric(priormean),
        is.finite(priormean),

        length(priorvar) == 1,
        is.numeric(priorvar),
        is.finite(priorvar),
        0 < priorvar,

        length(searchInt) == 2,
        is.numeric(searchInt),
        all(is.finite(searchInt)),
        0 <= searchInt[1], searchInt[1] < searchInt[2]
    )

    ## computing bound of probability of replication success
    limP <- porsBF01(level = level, dprior = dprior, sr = .Machine$double.eps,
                     priormean = priormean, priorvar = priorvar)
    if (power > limP) {
        warning(paste0("Power not achievable with specified design prior (at most ",
                       round(limP, 3), ")"))
        sr <- NaN
        outPow <- NaN
    } else {
        ## computing replication standard error sr
        rootFun <- function(sr) {
            porsBF01(level = level, dprior = dprior, sr = sr,
                     priormean = priormean, priorvar = priorvar) - power
        }
        res <- try(stats::uniroot(f = rootFun, interval = searchInt)$root, silent = TRUE)
        if (inherits(res, "try-error")) {
            sr <- NaN
            outPow <- NaN
            warning("Numerical problems, try adjusting searchInt")
        } else {
            sr <- res
            ## computing probability of replication success
            outPow <- porsBF01(level = level, dprior = dprior, sr = sr,
                               priormean = priormean, priorvar = priorvar)
        }
    }

    ## create output object
    out <- list("designPrior" = dprior, "power" = power,
                "powerRecomputed" = outPow, "sr" = sr,
                "c" = dprior$so^2/sr^2,
                type = paste("Bayes factor (in favor of H0) <=", signif(level, 3),
                             "(numerical computation)"))
    class(out) <- "ssdRS"
    return(out)
}


#' @title Probability of replication success based on Bayes factor
#'
#' @description This function computes the probability to achieve replication
#'     success based on a Bayes factor. The Bayes factor is oriented so that
#'     values above one indicate evidence for the null hypothesis of the effect
#'     size being zero, whereas values below one indicate evidence for the
#'     hypothesis of the effect size being non-zero (with normal prior assigned
#'     to it).
#'
#' @param level Bayes factor level below which replication success is achieved
#' @param dprior Design prior object
#' @param sr Replication standard error
#' @param priormean Mean of the normal prior under the alternative. Defaults to
#'     \code{0}
#' @param priorvar Variance of the normal prior under the alternative. Defaults
#'     to \code{1}
#'
#' @return The probability to achieve replication success
#'
#' @references
#'
#' Pawel, S., Consonni, G., and Held, L. (2022). Bayesian approaches to
#' designing replication studies. arXiv preprint.
#' \doi{10.48550/arXiv.2211.02552}
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## specify design prior
#' to1 <- 2
#' so1 <- 0.05
#' dprior <- designPrior(to = to1, so = so1, tau = 0.03)
#' porsBF01(level = 1/10, dprior = dprior, sr = c(0.05, 0.04))
#'
#' @export

porsBF01 <- function(level, dprior, sr, priormean = 0, priorvar = 1) {
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

        length(priormean) == 1,
        is.numeric(priormean),
        is.finite(priormean),

        length(priorvar) == 1,
        is.numeric(priorvar),
        is.finite(priorvar),
        0 < priorvar
    )

    ps <- vapply(X = sr, FUN = function(sr1) {
        ## compute probability of replication success
        g <- priorvar/sr1^2
        A <- sr1^2*(1 + 1/g)*(priormean^2/priorvar - 2*log(level) + log(1 + g))
        ## success region depends on direction of prior mean
        sregion <- successRegion(intervals = rbind(c(-Inf, -sqrt(A) - priormean/g),
                                                       c(sqrt(A) - priormean/g, Inf)))

        p <- pors(sregion = sregion, dprior = dprior, sr = sr1)
        return(p)
    }, FUN.VALUE = 1)
    return(ps)
}


## ## checking some stuff
## BF01a <- function(tr, sr, m, v) {
##     stats::dnorm(x = tr, mean = 0, sd = sr) /
##         stats::dnorm(x = tr, mean = m, sd = sqrt(v + sr^2))
## }
## BF01b <- function(tr, sr, m, v) {
##     sqrt(1 + v/sr^2)*exp(-0.5*((tr + m*sr^2/v)^2*v/sr^2/(sr^2 + v) - m^2/v))
## }
## tr <- 0.1
## sr <- 0.05
## m <- -0.2
## v <- 0.3
## BF01a(tr = tr, sr = sr, m = m, v = v)
## BF01b(tr = tr, sr = sr, m = m, v = v)
