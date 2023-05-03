#' @title Sample size determination for replication success based on
#'     the sceptical p-value
#'
#' @description This function computes the standard error required to achieve
#'     replication success with a certain probability and based on the sceptical
#'     p-value.
#'
#' @details The sceptical p-value is assumed to be uncalibrated as in Held
#'     (2020). The package ReplicationSuccess allows for sample size and power
#'     calculations with the recalibrated sceptical p-value
#'     (\url{https://CRAN.R-project.org/package=ReplicationSuccess}).
#'
#' @param level Threshold for the (one-sided) sceptical p-value below which
#'     replication success is achieved
#' @param dprior Design prior object
#' @param power Desired probability of replication success
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
#' Held, L. (2020). A new standard for the analysis and design of replication
#' studies (with discussion). Journal of the Royal Statistical Society: Series A
#' (Statistics in Society), 183(2), 431-448. \doi{10.1111/rssa.12493}
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## specify design prior
#' to1 <- 0.2
#' so1 <- 0.05
#' dprior <- designPrior(to = to1, so = so1, tau = 0.03)
#' ssdPs(level = 0.05, dprior = dprior, power = 0.9)
#'
#' @export

ssdPs <- function(level, dprior, power) {
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
        level < power
    )

    ## computing bound of probability of replication success
    limP <- porsPs(level = level, dprior = dprior, sr = 0)
    if (power > limP) {
        warning(paste0("Power not achievable with specified design prior (at most ",
                       round(limP, 3), ")"))
        sr <- NaN
        outPow <- NaN
    } else {
        ## computing replication standard error sr
        dpmean <- dprior$dpMean
        dpvar <- dprior$dpVar
        tau <- dprior$tau
        to <- dprior$to
        so <- dprior$so
        zo <- to/so
        if (sign(to) > 0) {
            za <- stats::qnorm(p = 1 - level)
        } else {
            za <- stats::qnorm(p = level)
        }
        zb <- stats::qnorm(p = power)
        A <- dpvar + tau^2 - so^2/((zo/za)^2 - 1)
        x <- (za*dpmean - zb*sqrt(dpmean^2 + (za^2 - zb^2)*A))/(za^2 - zb^2)
        sr <- sqrt(x^2 - so^2/((zo/za)^2 - 1))
        outPow <- porsPs(level = level, dprior = dprior, sr = sr)
        ## pow <- porsPs(level = level, dprior = dprior, sr = na.omit(srs))
        ## powequal <- abs(pow - power) <= 0.0001
        ## if (any(powequal)) {
        ##     sr <- srs[powequal]
        ##     outPow <- pow[powequal]
        ## } else {
        ##     sr <- NaN
        ##     outPow <- NaN
        ## }
    }

    ## create output object
    out <- list("designPrior" = dprior, "power" = power,
                "powerRecomputed" = outPow, "sr" = sr,
                "c" = dprior$so^2/sr^2,
                type = paste("sceptical p-value <=", signif(level, 3),
                             "(exact computation)"))
    class(out) <- "ssdRS"
    return(out)
}


#' @title Probability of replication success based on the sceptical p-value
#'
#' @description This function computes the probability to achieve replication
#'     success based on the sceptical p-value.
#'
#' @details The sceptical p-value is assumed to be uncalibrated as in Held
#'     (2020). The package ReplicationSuccess allows for sample size and power
#'     calculations with the recalibrated sceptical p-value
#'     (\url{https://CRAN.R-project.org/package=ReplicationSuccess}).
#'
#' @param level Threshold for the (one-sided) sceptical p-value below which
#'     replication success is achieved
#' @param dprior Design prior object
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
#' Held, L. (2020). A new standard for the analysis and design of replication
#' studies (with discussion). Journal of the Royal Statistical Society: Series A
#' (Statistics in Society), 183(2), 431-448. \doi{10.1111/rssa.12493}
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## specify design prior
#' to1 <- 0.2
#' so1 <- 0.05
#' dprior <- designPrior(to = to1, so = so1)
#' porsPs(level = 0.025, dprior = dprior, sr = c(0.05, 0.01))
#'
#' @export

porsPs <- function(level, dprior, sr) {
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
        ## success region depends on the direction of original study
        to <- dprior$to
        so <- dprior$so
        zo <- to/so
        za <- stats::qnorm(p = 1 - level)
        if (za > abs(zo)) {
            p <- 0
        } else {
            if (sign(to) >= 0) {
                int <- cbind(za*sqrt(sr1^2 + so^2/((zo/za)^2 - 1)), Inf)
            } else {
                zaNeg <- -za
                int <- cbind(-Inf, zaNeg*sqrt(sr1^2 + so^2/((zo/zaNeg)^2 - 1)))
            }
            sregion <- successRegion(intervals = int)
            p <- pors(sregion = sregion, dprior = dprior, sr = sr1)
        }
        return(p)
    }, FUN.VALUE = 1)
    return(ps)
}
