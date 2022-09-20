.ssBFr_ <- function(level = NA, t1e = NA, power, zo, h = 0, g = 0,
                    designPrior = c("predictive", "conditional", "EB"),
                    searchFac = 1) {
    ## input checks
    stopifnot(length(level) == 1,
              length(t1e) == 1,
              xor(is.na(t1e), is.na(level)),
              xor(is.numeric(t1e), is.numeric(level)),

              length(power) == 1,
              is.numeric(power),
              is.finite(power),
              0 < power, power < 1,

              length(zo) == 1,
              is.numeric(zo),
              is.finite(zo),

              length(h) == 1,
              is.numeric(h),
              is.finite(h),
              0 <= h,

              length(g) == 1,
              is.numeric(g),
              is.finite(g),
              0 <= g,

              !is.null(designPrior))
    designPrior <- match.arg(designPrior)

    stopifnot(length(searchFac) == 1,
              is.finite(searchFac),
              0 < searchFac)

    ## absolut cut-off approach
    if (!is.na(level)) {
        stopifnot(is.finite(level),
                  0 < level, level < 1)

        ## function to search with uniroot over relative variance c
        rootFun <- function(logc) {
            suppressWarnings({
                res <- cdfBFr(level = level, zo = zo, c = exp(logc), h = h, g = g,
                                designPrior = designPrior) - power
            })
            return(res)
        }
        res <- try(stats::uniroot(f = rootFun, interval = c(-5, 5)*searchFac)$root,
                   silent = TRUE)
        if (class(res) == "try-error") {
            c <- NaN
            t1e <- NaN
        } else {
            c <- exp(res)
        }

    ## relative cut-off approach
    } else {
        stopifnot(is.finite(t1e),
                  0 < t1e, t1e < power)

        ## function to search with uniroot over level for fixed t1e and c
        levelFun <- function(loglevel, c) {
            suppressWarnings({
                res <- cdfBFr(level = exp(loglevel), zo = zo, c = c, h = h, g = g,
                                designPrior = "H0") - t1e
            })
            return(res)
        }
        ## function to search with uniroot over relative variance c
        rootFun <- function(logc) {
            suppressWarnings({
                loglevel <- try(stats::uniroot(f = levelFun,
                                               interval = c(-5, 0)*searchFac,
                                               c = exp(logc))$root, silent = TRUE)
                if (class(loglevel) == "try-error") res <- NaN
                else {
                    res <- cdfBFr(level = exp(loglevel), zo = zo, c = exp(logc),
                                    h = h, g = g, designPrior = designPrior) - power
                }
            })
            return(res)
        }
        ## HACK evaluate root function for some logc to find limits
        ## for numerical search which are not NaN
        logcTest <- seq(-3, 3, length.out = 20)*searchFac
        rootTest <- sapply(X = logcTest, rootFun)
        low <- logcTest[utils::tail(which(!is.na(rootTest) & rootTest < 0), n = 1)]
        up <- logcTest[utils::head(which(!is.na(rootTest) & rootTest > 0), n = 1)]
        res <- try(stats::uniroot(f = rootFun, interval = c(low, up))$root,
                   silent = TRUE)
        if (class(res) == "try-error") {
            c <- NaN
        } else {
            c <- exp(res)
            level <- exp(stats::uniroot(f = levelFun, interval = c(-5, 0)*searchFac,
                                        c = c)$root)
        }
    }

    ## return relative variance, actual error rates, and inputs
    if (!is.na(c)) {
        suppressWarnings({
            power <- cdfBFr(level = level, zo = zo, c = c, h = h, g = g,
                              designPrior = designPrior)
            t1e <- cdfBFr(level = level, zo = zo, c = c, h = h, g = g,
                            designPrior = "H0")
        })
    }
    out <- list(c = c, level = level, t1e = t1e, power = power, zo = zo, h = h,
                g = g, designPrior = designPrior)
    return(out)
}

## stack in a matrix to transpose afterwards
.ssBFr__ <- Vectorize(FUN = .ssBFr_, SIMPLIFY = TRUE)

#' @title Sample size determination for Bayes factor analysis of replication
#' study
#' 
#' @description Computes relative variance \eqn{c = \sigma^2_o/\sigma^2_r}{c =
#'     sigma^2_o/sigma^2_r} to achieve desired power for replication success
#'     with the Bayes factor contrasting the (composite) sceptical prior to the
#'     advocacy prior (\link{BFr}). Set \code{g = 0} for the replication Bayes
#'     factor.
#'
#' For some types of effect sizes (e.g. mean differences) the relative variance
#' can directly be interpreted as the relative sample size
#' of the replication study, i.e. \eqn{c = n_r/n_o}{c = n_r/n_o}.
#' 
#' Allows to either use a *relative cut-off* (similar to Weiss, 1997) or an
#' *absolute cut-off* (similar to De Santis, 2004). In the relative cut-off
#' approach, the level to threshold the Bayes factor is determined based on the
#' specified error probabilities.
#' 
#' @param level Desired level of replication success, only specify for absolute
#'     cut-off approach
#' @param power Desired power for replication success
#' @param t1e Desired type-I error rate, only specify for relative cut-off
#'     approach
#' @param zo \eqn{z}{z}-value from original study, \eqn{z_o =
#'     \hat{\theta}_o/\sigma_o}{zo = hat(theta)_o/sigma_o}, i.e. original effect
#'     estimate divided by standard error
#' @param h The relative between-study heterogeneity, i.e. the ratio of the
#'     effect size heterogeneity variance to the variance of the original effect
#'     estimate. Default is 0 (no heterogeneity). Is only taken into account
#'     when \code{designPrior} = "predictive" or \code{designPrior} = "EB".
#' @param g Relative variance of \eqn{\mathrm{N}(0, g \cdot \sigma_o^2)}{N(0,
#'     g*sigma_o^2)} prior for the effect size under \eqn{H_S}{HS}. Set to zero
#'     for a point-null hypothesis (i.e. to obtain the replication Bayes
#'     factor). Defaults to zero.
#' @param designPrior Desing prior for the effect size, can be set to
#'     "conditional" (point-mass at the original effect estimate), "predictive"
#'     (posterior-predictive distribution of the effect size conditional on
#'     original effect estimate and flat prior), or "EB" (posterior-predictive
#'     distribution of the effect size conditional on original effect estimate
#'     and g-prior with g estimated by empirical Bayes). Defaults to
#'     "predictive".
#' @param searchFac A positive number specifying how much the search interal for
#'     the numerical search of log(c) (and log(level) if the relative cut-off
#'     approach is used) should be scaled. Defaults to 1.
#' 
#' @return A list containing the relative variance \eqn{c =
#'     \sigma^2_o/\sigma^2_r}{c = sigma^2_o/sigma^2_r} to achieve replication
#'     success at specified level. Also contains the actual type I error rate,
#'     the power, the level of replication success, and the other inputs.
#' 
#' @references Weiss, R. (1997). Bayesian sample size calculations for
#'     hypothesis testing. Journal of the Royal Statistical Society: Series D
#'     (The Statistician), 46(2), 185–191. \doi{10.1111/1467-9884.00075}
#'  
#' De Santis, F. (2004). Statistical evidence and sample size determination for
#' Bayesian hypothesis testing. Journal of Statistical Planning and Inference,
#' 124(1), 121–144. \doi{10.1016/s0378-3758(03)00198-8}
#'
#' Pawel, S. and Held, L. (2020). The sceptical Bayes factor for the assessement
#' of replication success. Preprint. <https://arxiv.org/abs/2009.01520>.
#'
#' Verhagen, J. and Wagenmakers, E. J. (2014). Bayesian tests to quantify the result
#' of a replication attempt. Journal of Experimental Psychology: General,
#' 145:1457-1475. \doi{10.1037/a0036731}
#'
#' Ly, A., Etz, A., Marsman, M., & Wagenmakers, E.-J. (2018). Replication Bayes
#' factors from evidence updating. Behavior Research Methods, 51(6), 2498-2508.
#' \doi{10.3758/s13428-018-1092-x}
#' 
#' @author Samuel Pawel
#' 
#' @examples
#' ## sample size calculation for replication BF (g = 0)
#' ssBFr(level = 1/3, power = 0.8, zo = 3, g = 0) ## absolute cut-off
#' ssBFr(t1e = 0.05*0.025, power = 0.8, zo = 3, g = 0) ## relative cut-off
#' 
#' @seealso \code{\link{cdfBFr}}
#'   
#'
ssBFr <- function(level = NA, t1e = NA, power, zo, h = 0, g = 0,
                  designPrior = "predictive", searchFac = 1) {
    t(.ssBFr__(level = level, t1e = t1e, power = power, zo = zo, g = g, h = h,
               designPrior = designPrior, searchFac = searchFac))
}
