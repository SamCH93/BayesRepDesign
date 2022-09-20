.cdfBFr_ <- function(level, zo, c, g = 0,
                       designPrior = c("predictive", "conditional", "EB", "H0",
                                       "H1"),
                       h = 0, mu = 0) {
    ## input checks
    stopifnot(length(level) == 1,
              is.numeric(level),
              is.finite(level),
              0 < level,

              length(zo) == 1,
              is.numeric(zo),
              is.finite(zo),

              length(c) == 1,
              is.numeric(c),
              is.finite(c),
              0 <= c,

              length(g) == 1,
              is.numeric(g),
              is.finite(g),
              0 <= g,

              !is.null(designPrior))
    designPrior <- match.arg(designPrior)

    stopifnot(length(h) == 1,
              is.numeric(h),
              is.finite(h),
              0 <= h,

              length(mu) == 1,
              is.numeric(mu),
              is.finite(mu))

    ## determine parameters of predictive distribution of zr
    if (designPrior == "H0") {
        mu <- 0
        sigma2 <- 1
    } else if (designPrior == "H1") {
        mu <- mu*sqrt(c)
        sigma2 <- 1
    } else if (designPrior == "predictive") {
        mu <- zo*sqrt(c)
        sigma2 <- 1 + c*(1 + 2*h)
    } else if (designPrior == "EB") {
        s <- pmax(1 - (1 + h)/zo^2, 0)
        mu <- s*zo*sqrt(c)
        sigma2 <- s*c*(1 + h) + 1 + c*h
    } else { ## conditional design prior
        mu <- zo*sqrt(c)
        sigma2 <- 1
    }


    ## compute probability that BFr < level
    k <- log(level)
    if (isTRUE(all.equal(g, 1))) {
        D <- (k + 0.5*zo^2/(1/c + 1))*(1 + c)/(zo*sqrt(c))
        pow <- stats::pnorm(q = sign(zo)*(D - mu)/sqrt(sigma2), mean = 0, sd = 1)
    } else {
        A <- log((1 + c)/(1 + g*c)) - 2*k + zo^2/(1 - g)
        B <- (1 - g)/(1 + c*g)/(1 + 1/c)
        M <- zo*sqrt(c)*(1/c + g)/(g - 1)
        lambda <- (mu - M)^2/sigma2
        prob <- stats::pchisq(q = A/B/sigma2, df = 1, ncp = lambda,
                              lower.tail = FALSE)

        ## when g > 1, the inequality flips and we need to take the complement
        if (g > 1) pow <- 1 - prob
        else pow <- prob
    }
    return(pow)
}

#' @title Conditional power of replication Bayes factor
#'
#' @description Computes conditional probability that the Bayes factor from
#'     \link{BFr} is smaller than specified level
#'
#' @param level Bayes factor level
#' @param zo \eqn{z}{z}-value from original study, \eqn{z_o =
#'     \hat{\theta}_o/\sigma_o}{zo = hat(theta)_o/sigma_o}, i.e. original effect
#'     estimate divided by standard error
#' @param c Relative variance \eqn{c = \sigma^2_o/\sigma^2_r}{c =
#'     sigma^2_o/sigma^2_r}, i.e the variance of the original effect estimate
#'     relative to the variance of the replication effect estimate
#' @param g Relative variance of \eqn{\mathrm{N}(0, g \cdot \sigma_o^2)}{N(0,
#'     g*sigma_o^2)} prior for the effect size under \eqn{H_S}{HS}. By default
#'     set zero, so \eqn{H_S}{HS} is a point-null hypothesis and the Bayes
#'     factor is the replication Bayes factor under normality.
#' @param H0 Logical vector indicating whether probability should be computed
#'     under \eqn{H_0}{H0} or under \eqn{H_1}{H1}
#' @param designPrior Desing prior for the effect size, can be set to
#'     "conditional" (point-mass at the original effect estimate), "predictive"
#'     (posterior-predictive distribution of the effect size conditional on
#'     original effect estimate and flat prior), or "EB" (posterior-predictive
#'     distribution of the effect size conditional on original effect estimate
#'     and g-prior with g estimated by empirical Bayes). Defaults to
#'     "predictive".
#' @param h The relative between-study heterogeneity, i.e. the ratio of the
#'     effect size heterogeneity variance to the variance of the original effect
#'     estimate. Default is 0 (no heterogeneity). Is only taken into account
#'     when \code{designPrior} = "predictive" or \code{designPrior} = "EB".
#' @param mu The true effect size divided by the standard error from the
#'     original study. Is only taken into account when \code{designPrior} = "H1".
#'
#' @return Probability of Bayes factor being smaller the specified level
#'
#' @references
#' Weiss, R. (1997). Bayesian sample size calculations for hypothesis
#' testing. Journal of the Royal Statistical Society: Series D (The Statistician),
#' 46(2), 185–191. \doi{10.1111/1467-9884.00075}
#'
#' Pawel, S. and Held, L. (2020). The sceptical Bayes factor for the assessement
#' of replication success. Preprint. <https://arxiv.org/abs/2009.01520>.
#'
#' @seealso \code{\link{ssBFr}}
#'
#' @author Samuel Pawel
#'
#' @examples
#' cdfBFr(level = 1/10, zo = 3, c = 1, g = 0)
#'
#' ## example from page 188, 189 in Weiss (1997):
#' tau1 <- 1
#' sigma <- 5
#' mu1 <- 2
#' n <- c(20, 40, 60, 65, 66, 67)
#' zo <- mu1/tau1
#' c <- tau1^2/(sigma^2/n)
#' g <- 0 # point-null
#' cdfBFr(level = 1, zo = zo, c = c, g = g, designPrior = "H0")
#' cdfBFr(level = 1, zo = zo, c = c, g = g, designPrior = "predictive")
#'
#'
cdfBFr <- Vectorize(.cdfBFr_)
