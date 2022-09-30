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
#' @param searchInt Search interval for standard errors
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
#' ssd(sregionfun = sregionfunSig, dprior = dprior, power = 0.8)
#'
#' @export

ssd <- function(sregionfun, dprior, power,
                searchInt = c(.Machine$double.eps^0.5, 4),
                ...) {
    ## input checks
    stopifnot(
        is.function(sregionfun),

        class(dprior) == "designPrior",

        length(power) == 1,
        is.numeric(power),
        is.finite(power),
        0 < power, power < 1,

        length(searchInt) == 2,
        is.numeric(searchInt),
        all(is.finite(searchInt)),
        0 <= searchInt[1], searchInt[1] < searchInt[2]
    )

    ## check whether specified power achievable
    sregionLim <- sregionfun(.Machine$double.eps)
    limP <- pors(sregion = sregionLim, dprior = dprior, sr = .Machine$double.eps)
    if (power > limP) {
        warning(paste0("Power not achievable with specified design prior (at most ",
                       round(limP, 3), ")"))
        sr <- NaN
        outPow <- NaN
    } else {
        ## numerical search for log replication standard error such that probability
        ## of replication success = power
        rootFun <- function(logsr) {
            sregion <- sregionfun(exp(logsr))
            pors(sregion = sregion, dprior = dprior, sr = exp(logsr)) - power
        }
        res <- try(stats::uniroot(f = rootFun, interval = log(searchInt),
                                  ... = ...)$root)
        if (inherits(res, "try-error")) {
            sr <- NaN
            outPow <- NaN
            warning("Numerical problems, try adjusting searchInt")
        } else {
            sr <- exp(res)
            outPow <- rootFun(log(sr)) + power
        }
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
#' ssd1 <- ssd(sregionfun = sregionfunSig, dprior = dprior, power = 0.8)
#' print(ssd1)
#' @export
print.ssdRS <- function(x, ...) {
    print(x$designPrior)
    cat("\npower =", signif(x$power, 2), "(specified)")
    cat("\npower =", signif(x$powerRecomputed, 2), "(recomputed with sr)")
    cat("\nsr =", signif(x$sr, 2))
    cat("\nc = so^2/sr^2 ~= nr/no =", signif(x$c, 2))
    cat("\n\n")
    invisible(x)
}
