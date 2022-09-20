#' @title Sample size determination for replication success
#'
#' @description This function computes the standard error required to achieve
#'     replication success with a certain probability and based on a certain
#'     type of success region.
#'
#' @details The probability of replication success .....
#'
#' @md
#'
#' @param sregionfun Function that returns the success region for replication
#'     effect estimate as a function of the replication standard error
#' @param dprior Design prior object
#' @param power Desired probability of replication success
#' @param searchinterval Search interval for standard errors
#' @param ... other arguments for uniroot
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
#' dprior <- designPrior(to = to1, so = so1)
#'
#' ## compute required standard error for significance at one-sided 2.5%
#' sregionfunSig <- function(sr, alpha = 0.025) {
#'     successRegion(intervals = cbind(stats::qnorm(p = 1- alpha)*sr, Inf))
#' }
#' ssdRS(sregionfun = sregionfunSig, dprior = dprior, power = 0.8)
#'
#' ## compute required standard error for replication BF = 1/10
#' ## TODO fix
#' # sregionfunBFr <- function(sr, level = 1/10, to = to1, so = so1) {
#' #   zo <- to/so
#' #   c <- so^2/sr^2
#' #   m <- sqrt((zo^2 + log(1 + c) - 2*log(level))*(1 + 1/c))
#' #   successRegion(intervals = rbind(c(-Inf, sr*(-sqrt(m) - zo/sqrt(c))),
#' #                                   c(sr*(sqrt(m) - zo/sqrt(c)), Inf)))
#' #
#' # sdRS(sregionfun = sregionfunBFr, dprior = dprior, power = 0.6)
#'
#' @export

ssdRS <- function(sregionfun, dprior, power,
                  searchinterval = c(.Machine$double.eps^0.5, 100),
                  ...) {
    ## input checks
    stopifnot(
        is.function(sregionfun),

        class(dprior) == "designPrior",

        length(power) == 1,
        is.numeric(power),
        is.finite(power),
        0 < power, power < 1,

        length(searchinterval) == 2,
        is.numeric(searchinterval),
        all(is.finite(searchinterval)),
        searchinterval[2] > searchinterval[1]
    )

    ## numerical search for log replication standard error such that probability
    ## of replication success = power
    rootFun1 <- function(logsr) {
        sregion <- sregionfun(exp(logsr))
        probRS(sregion = sregion, dprior = dprior, sr = exp(logsr)) - power
    }
    rootFun <- Vectorize(rootFun1)
    res <- try(stats::uniroot(f = rootFun, interval = log(searchinterval),
                              ... = ...)$root)
    if (inherits(res, "try-error")) {
        sr <- NaN
        outPow <- NaN
    } else {
        sr <- exp(res)
        outPow <- rootFun(log(sr)) + power
    }

    ## create output object
    out <- list("designPrior" = dprior, "power" = power,
                "powerRecomputed" = outPow, "sr" = sr,
                "c" = dprior$so^2/sr^2)
    class(out) <- "ssdRS"
    return(out)
}

#' Print method for ssdRS object
#' @method print ssdRS
#' @param x An ssdRS object
#' @param ... Other arguments
#' @examples
#' ## specify design prior
#' to1 <- 2
#' so1 <- 1
#' dprior <- designPrior(to = to1, so = so1)
#'
#' ## compute required standard error for significance at one-sided 2.5%
#' sregionfunSig <- function(sr, alpha = 0.025) {
#'     successRegion(intervals = cbind(stats::qnorm(p = 1- alpha)*sr, Inf))
#' }
#' ssd1 <- ssdRS(sregionfun = sregionfunSig, dprior = dprior, power = 0.8)
#' print(ssd1)
#' @export
print.ssdRS <- function(x, ...) {
    print(x$designPrior)
    cat("\npower =", signif(x$power, 2), "(specified)")
    cat("\npower =", signif(x$powerRecomputed, 2), "(recomputed with sr)")
    cat("\nsr =", signif(x$sr, 2))
    cat("\nc =", signif(x$c, 2))
    cat("\n\n")
    invisible(x)
}
