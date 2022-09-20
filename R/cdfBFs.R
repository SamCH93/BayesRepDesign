.cdfBFs_ <- function(level, zo, c,
                       designPrior = c("predictive", "conditional", "EB", "H0",
                                       "H1"),
                       h = 0, mu = 0) {

    ## input checks
    stopifnot(length(level) == 1,
              is.numeric(level),
              is.finite(level),
              0 < level, level <= 1,

              length(zo) == 1,
              is.numeric(zo),
              is.finite(zo),

              length(c) == 1,
              is.numeric(c),
              is.finite(c),
              0 < c,

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

    ## compute sufficiently sceptical prior variance for specified level
    g <- vss(x = zo, gamma = level)
    if (is.nan(g)) {
        pow <- 0
    } else {
        ## compute probability that RS at level
        ## TODO: find solution for numerical problems with non-central chi-squared
        ## calculations when g is a bit above or below 1, currently set numerical
        ## tolerance quite large
        if (isTRUE(all.equal(g, 1, tolerance = 1e-2))) {
            D <- (-log(2) + zo^2*(0.5 + 1/(1/c + 1)))/(2*zo*sqrt(c))*(1 + c)
            ## D <- 0.5*(zo*sqrt(c) + (1/c + 1)*sqrt(c)*(0.5*zo - log(2)/zo))
            pow <- stats::pnorm(q = sign(zo)*(mu - D)/sqrt(sigma2), mean = 0, sd = 1)
            ## if (zo > 0) pow <- stats::pnorm(q = D, mean = mu, sd = sqrt(sigma2),
            ##                                 lower.tail = FALSE)
            ## else pow <- stats::pnorm(q = D, mean = mu, sd = sqrt(sigma2))
        } else {
            A <- log((1/c + 1)/(1/c + g)/(1 + g)) + zo^2/(1/g + 1) + zo^2/(1 -g)
            B <- c*(1 - g)/(1 + c*g)/(1 + c)
            M <- zo*(1 + c*g)/sqrt(c)/(g - 1)
            lambda <- (mu - M)^2/sigma2
            res <- stats::pchisq(q = A/B/sigma2, df = 1, ncp = lambda,
                                 lower.tail = FALSE)
            ## when g > 1, the inequality flips and we need to take the complement...
            if (g < 1) pow <- res
            else pow <- 1 - res
        }
    }
    return(pow)
}

#' @title Conditional Power of sceptical Bayes factor
#' 
#' @description Computes power for replication success at a specified level
#'     conditional on result from original study.
#' 
#' @param level Desired level of replication success
#' @param zo \eqn{z}{z}-value from original study, \eqn{z_o =
#'     \hat{\theta}_o/\sigma_o}{zo = hat(theta)_o/sigma_o}, i.e. original effect
#'     estimate divided by standard error
#' @param c Relative variance \eqn{c = \sigma^2_o/\sigma^2_r}{c =
#'     sigma^2_o/sigma^2_r}, i.e the variance of the original effect estimate
#'     relative to the variance of the replication effect estimate
#' @param designPrior Desing prior for the effect size, can be set to
#'     "conditional" (point-mass at the original effect estimate), "predictive"
#'     (posterior-predictive distribution of the effect size conditional on
#'     original effect estimate and flat prior), "EB" (posterior-predictive
#'     distribution of the effect size conditional on original effect estimate
#'     and g-prior with g estimated by empirical Bayes), or "H0" (null
#'     hypothesis that effect size is exactly zero). Defaults to "predictive".
#' @param h The relative between-study heterogeneity, i.e. the ratio of the
#'     effect size heterogeneity variance to the variance of the original effect
#'     estimate. Default is 0 (no heterogeneity). Is only taken into account
#'     when \code{designPrior} = "predictive" or \code{designPrior} = "EB".
#' @param mu The true effect size divided by the standard error from the
#'     original study. Is only taken into account when \code{designPrior} = "H1".
#' 
#' @return Power to achieve replication success at the specified level
#' 
#' @references
#' Pawel, S. and Held, L. (2020). The sceptical Bayes factor for the assessement
#' of replication success. Preprint. <https://arxiv.org/abs/2009.01520>.
#'
#' 
#' @seealso \code{\link{ssBFs}}
#' 
#' @author Samuel Pawel
#' 
#' @examples
#' cdfBFs(level = 1/3, zo = 3, c = 2, designPrior = "conditional")
#' 
#'
cdfBFs <- Vectorize(.cdfBFs_)

## ## recreate plot from paper
## threshBFs <- levelBFs(type = "T1E", T1Elevel = 0.05*0.025)
## cseq <- exp(seq(log(1/4), log(4), length.out = 100))
## zoseq <- seq(1.75, 3, 0.25)
## dPrs <- c("conditional", "predictive")
## applyGrid <- expand.grid(dPr = dPrs, zo = zoseq)

## plotList <- lapply(X = seq(1, nrow(applyGrid)), FUN = function(i) {
##   data.frame(c = cseq, zo = applyGrid$zo[i], dPr = applyGrid$dPr[i],
##              power = cdfBFs(level = threshBFs, zo = applyGrid$zo[i],
##                               c = cseq, designPrior = applyGrid$dPr[i]))
## })
## plotDF <- do.call("rbind", plotList)

## library(ggplot2)
## ggplot(data = plotDF, aes(x = c, y = power)) +
##   geom_line(aes(linetype = dPr)) +
##   facet_wrap(~ zo, labeller = label_bquote(italic(z)[o] == .(zo))) +
##   guides(linetype = guide_legend(title = "")) + 
##   scale_y_continuous(labels = scales::percent, minor_breaks = NULL,
##                      limits = c(0, 1)) +
##   scale_x_continuous(name = bquote(italic(c)), minor_breaks = NULL) +
##   theme_bw() +
##   theme(text = element_text(family = "serif"), legend.position = "top")

## check that correct behaviour around g = 1
## zo <- 2.5
## c <- 1.5
## level <- BFo(zo = zo, g = c(0.99, 1, 1.01))
## cdfBFs(level = level, zo = zo, c = c)
