#' @title Sample size determination for replication success
#'
#' @description This function computes the standard error of the replication
#'     effect estimate required to achieve replication success with a certain
#'     probability and based on a certain type of success region.
#'
#' @param sregionfun Function that returns the success region for replication
#'     effect estimate as a function of the replication standard error
#' @param dprior Design prior object
#' @param power Desired probability of replication success
#' @param nsites Number of sites. Defaults to \code{1}. The sites are assumed to
#'     have the same sample size
#' @param searchInt Search interval for standard errors
#' @param ... Other arguments passed to \code{uniroot}
#'
#' @return Returns an object of class \code{"ssdRS"} which is a list containing:
#' \tabular{ll}{
#'    \code{designPrior} \tab The specified \code{"designPrior"} object \cr
#'    \tab \cr
#'    \code{power} \tab The specified power \cr
#'    \tab \cr
#'    \code{powerRecomputed} \tab The recomputed power \cr
#'    \tab \cr
#'    \code{sr} \tab The required replication standard error \cr
#'    \tab \cr
#'    \code{c} \tab The required relative sample size \code{c = nr/no}
#'    (assuming \code{so = unitSD/no} and \code{sr = unitSD/nr}) \cr
#' }
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

ssd <- function(sregionfun, dprior, power, nsites = 1,
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

        length(nsites) == 1,
        is.numeric(nsites),
        is.finite(nsites),
        nsites > 0,

        length(searchInt) == 2,
        is.numeric(searchInt),
        all(is.finite(searchInt)),
        0 <= searchInt[1], searchInt[1] < searchInt[2]
    )

    ## check whether specified power achievable
    sregionLim <- sregionfun(.Machine$double.eps)
    limP <- pors(sregion = sregionLim, dprior = dprior, sr = .Machine$double.eps,
                 nsites = nsites)
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
            pors(sregion = sregion, dprior = dprior, sr = exp(logsr),
                 nsites = nsites) - power
        }
        res <- try(stats::uniroot(f = rootFun, interval = log(searchInt),
                                  ... = ...)$root, silent = TRUE)
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
                "c" = dprior$so^2/sr^2,
                type = "method agnostic success region (numerical computation)")
    class(out) <- "ssdRS"
    return(out)
}

#' Print method for class \code{"ssdRS"}
#' @method print ssdRS
#'
#' @param x Object of class \code{"ssdRS"}
#' @param ... Other arguments (for consistency with the generic)
#'
#' @return Prints text summary in the console and invisibly returns the
#'     \code{"ssdRS"} object
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
#' ssd1 <- ssd(sregionfun = sregionfunSig, dprior = dprior, power = 0.8)
#' print(ssd1)
#' @export
print.ssdRS <- function(x, ...) {
    ## cat("========================================================================\n")
    cat("       Bayesian sample size calculation for replication studies\n")
    ## cat("       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
    cat("       ========================================================\n\n")
    cat("success criterion and computation\n")
    cat("------------------------------------------------------------------------")
    cat("\n ", x$type, "\n\n")
    print(x$designPrior)
    cat("\nprobability of replication success\n")
    cat("------------------------------------------------------------------------")
    cat("\n  PoRS =", signif(x$power, 2), ": specified")
    cat("\n  PoRS =", signif(x$powerRecomputed, 2), ": recomputed with sr\n")
    cat("\nrequired sample size\n")
    cat("------------------------------------------------------------------------")
    cat("\n  sr =", signif(x$sr, 2), ": required standard error of replication effect estimate")
    cat("\n  c = so^2/sr^2 ~= nr/no =", signif(x$c, 2), ": required relative variance / sample size")
    cat("\n")
    ## cat("========================================================================\n")
    invisible(x)
}
