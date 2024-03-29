#' @title Sample size related to standard error and unit standard deviation
#'
#' @description This function computes the sample size related to a specified
#'     standard error \eqn{\sigma}{\code{se}} and unit standard deviation
#'     \code{unitSD}, which is the standard deviation of one effective unit (one
#'     measurement, one pair of measurements, one event, etc.). The relationship
#'     \eqn{\sigma = \code{unitSD}/\sqrt{n}}{\code{se} = \code{unitSD}/sqrt(n)} is
#'     assumed. The unit standard deviation depends on the parameter type and
#'     the assumptions underlying the standard error calculation. The default is
#'     \code{unitSD = 2} which is, under some assumptions, a reasonable
#'     approximation to the unit standard deviation for standardized mean
#'     differences and log odds/hazard/rate ratios, see Section 2.4 in
#'     Spiegelhalter et al. (2004).
#'
#' @param se Standard error
#' @param unitSD Unit standard deviation. Defaults to \code{2}
#'
#' @return The sample size corresponding to the specified standard error and
#'     unit standard deviation
#'
#' @references
#'
#' Pawel, S., Consonni, G., and Held, L. (2023). Bayesian approaches to
#' designing replication studies. Psychological Methods.
#' \doi{10.1037/met0000604}
#'
#' Spiegelhalter, D.J., Abrams, K.R., Myles, J.P. (2004). Bayesian approaches to
#' clinical trials and health care evaluation. Wiley.
#' \doi{10.1002/0470092602}
#'
#'
#' @author Samuel Pawel
#'
#' @examples
#' smd1 <- 0.3
#' so1 <- 0.05
#' dprior <- designPrior(to = smd1, so = so1)
#' ssd1 <- ssdSig(level = 0.025, dprior = dprior, power = 0.8)
#' se2n(se = ssd1$sr, unitSD = 2) # required n
#'
#' @export

se2n <- function(se, unitSD = 2) {
    ## input checks
    stopifnot(
        length(se) > 0,
        is.numeric(se),
        all(is.finite(se)),
        all(0 <= se),

        length(unitSD) == 1,
        is.numeric(unitSD),
        is.finite(unitSD),
        0 < unitSD
    )
    ## TODO implement more precise conversions for certain effect size types?
    n <- ceiling(unitSD^2/se^2)
    return(n)
}
