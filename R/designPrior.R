#' @title Design prior for effect size
#' 
#' @description Creates a design prior for the effect size which can then be
#'     used for power and sample size calculations of a replication study. The
#'     design prior is obtained from updating an initial prior for the effect
#'     size by the data from the original study. A normal-normal hierarchical
#'     model is assumed, see Pawel et al. (2022) for details.
#'
#' @param to Effect estimate from original study
#' @param so Standard error of effect estimate from original study
#' @param mu The initial prior mean. Defaults to \code{0}
#' @param sp The initial prior standard deviation. Defaults to \code{Inf} (an
#'     improper uniform prior)
#' @param tau The initial prior heterogeneity standard deviation. Defaults to
#'     \code{0} (no heterogeneity)
#' @param g The relative initial prior variance \code{g} =
#'     \code{sp^2}/(\code{tau^2} + \code{so^2}) (alternative parametrization of
#'     prior standard deviation \code{sp})
#' @param h The relative initial prior heterogeneity variance \code{h} =
#'     \code{tau^2}/\code{so^2} (alternative parametrization of prior
#'     heterogeneity standard deviation \code{tau})
#' @param type Shortcut for special parameter combinations. The available
#'     options are \code{NA}, \code{"conditional"}, \code{"predictive"}, and
#'     \code{"EB"} (see details). Defaults to \code{NA}
#'
#' @details The \code{"conditional"} design prior corresponds to a point mass at
#'     the original effect estimate, i.e., assuming that the true effect size is
#'     equal to the original effect estimate. The \code{"predictive"} design
#'     prior is obtained from updating a uniform initial prior by the likelihood
#'     of the original data. The \code{"EB"} design prior is obtained by
#'     empirical Bayes estimation of the variance of the normal prior and
#'     induces adaptive shrinkage that depends on the p-value of the original
#'     effect estimate.
#'
#' @return
#'
#' Returns an object of class \code{"designPrior"} which is a list containing:
#'
#' \tabular{ll}{
#'    \code{dpMean} \tab The computed mean of the design prior \cr
#'    \tab \cr
#'    \code{dpVar} \tab The computed variance of the design prior \cr
#'    \tab \cr
#'    \code{to} \tab The specified original effect estimate \cr
#'    \tab \cr
#'    \code{so} \tab The specified original standard error \cr
#'    \tab \cr
#'    \code{mu} \tab The specified mean of the initial prior \cr
#'    \tab \cr
#'    \code{sp} \tab The specified standard deviation of the initial prior \cr
#'    \tab \cr
#'    \code{tau} \tab The specified heterogeneity variance \cr
#' }
#'
#' @author Samuel Pawel
#'
#' @references
#'
#' Pawel, S., Consonni, G., and Held, L. (2022). Bayesian approaches to
#' designing replication studies. arXiv preprint.
#' \doi{10.48550/arXiv.2211.02552}
#'
#' @seealso \code{\link{pors}}, \code{\link{ssd}}
#' 
#' @examples
#' designPrior(to = 1.1, so = 1)
#' @export
designPrior <- function(to, so, mu = 0, sp = Inf, tau = 0,
                        g = sp^2/(tau^2 + so^2), h = tau^2/so^2,
                        type = c(NA, "conditional", "predictive", "EB")) {


    ## input checks
    stopifnot(
        length(to) == 1,
        is.numeric(to),
        is.finite(to),

        length(so) == 1,
        is.numeric(so),
        is.finite(so),
        0 < so,

        length(mu) == 1,
        is.numeric(mu),
        is.finite(mu),

        length(sp) == 1,
        is.numeric(sp),
        !is.na(sp), !is.nan(sp),
        0 <= sp,

        length(tau) == 1,
        is.numeric(tau),
        ## is.finite(tau),
        0 <= tau,

        length(g) == 1,
        is.numeric(g),
        !is.na(g), !is.nan(g),
        0 <= g,

        length(h) == 1,
        is.numeric(h),
        is.finite(h),
        0 <= h,

        !is.null(type)
    )
    type <- match.arg(type)

    ## recompute absolute parameters based on relative ones
    tau <- sqrt(h)*so
    sp <- sqrt(g*(so^2 + tau^2))

    ## shortcuts
    if (!is.na(type)) {
        if (type == "conditional") {
            mu <- to
            tau <- 0
            sp <- 0
        }
        if (type == "predictive") {
            sp <- Inf
        }
        if (type == "EB") {
            sp <- sqrt(pmax((to - mu)^2 - so^2 - tau^2, 0))
            g <- sp^2/(so^2 + tau^2)
        }
        g <- sp^2/(so^2 + tau^2)
    }

    ## compute mean and variance of design prior by standard Bayesian updating
    m <- to/(1 + 1/g) + mu/(1 + g)
    v <- so^2*(1 + h)/(1 + 1/g)
    dp <- list(dpMean = m, dpVar = v,
               to = to, so = so, mu = mu, sp = sp, tau = tau)
    class(dp) <- "designPrior"
    return(dp)
}


#' Print method for class \code{"designPrior"}
#' @method print designPrior
#'
#' @param x Object of class \code{"designPrior"}
#' @param ... Other arguments (for consistency with the generic)
#'
#' @return Prints text summary in the console and invisibly returns the
#'     \code{"designPrior"} object
#'
#' @author Samuel Pawel
#'
#' @examples
#' dp <- designPrior(to = 0.5, so = 0.05, sp = 0.2, tau = 0.1)
#' print(dp)
#' @export

print.designPrior <- function(x, ...) {
    cat("original data and initial prior for effect size\n")
    cat("------------------------------------------------------------------------")
    cat("\n  to =", signif(x$to, 2), ": original effect estimate")
    cat("\n  so =", signif(x$so, 2), ": standard error of original effect estimate")
    cat("\n  tau =", signif(x$tau, 2),
        ": assumed heterogeneity standard deviation")
    cat("\n  N(mean = ", signif(x$mu, 2), ", sd = ", signif(x$sp, 2), ") ",
        ": initial normal prior", sep = "")

    cat("\n\ndesign prior for effect size\n")
    cat("------------------------------------------------------------------------")
    cat("\n  N(mean = ", signif(x$dpMean, 2), ", sd = ",
        signif(sqrt(x$dpVar), 2), ") ",
        ": normal design prior", sep = "")
    cat("\n")
    invisible(x)
}

#' Density method for class \code{"designPrior"}
#' @method density designPrior
#'
#' @param x Object of class \code{"designPrior"}
#' @param ... Other arguments passed to \code{stats::dnorm}
#'
#' @return Returns the density function of the design prior
#'
#' @author Samuel Pawel
#'
#' @examples
#' dp <- designPrior(to = 2.3123, so = 0.1, mu = 1.1, tau = 0.2)
#' f <- density(dp)
#' tseq <- seq(1, 3.5, 0.01)
#' plot(tseq, f(theta = tseq), type = "l", xlab = "theta", ylab = "Design prior density")
#' @importFrom stats density
#' @export
density.designPrior <- function(x, ...) {
    ## return design prior density function for the overall effect size (theta)
    densFun <- function(theta) {
        d <- stats::dnorm(x = theta, mean = x$dpMean, sd = sqrt(x$dpVar), ...)
        return(d)
    }
    return(densFun)
}

#' Plot method for class \code{"designPrior"}
#' @method plot designPrior
#'
#' @param x Object of class \code{"designPrior"}
#' @param ... Other arguments passed to \code{plot}
#'
#' @return Plots the density of the design prior
#'
#' @author Samuel Pawel
#'
#' @examples
#' dp <- designPrior(to = 2.3123, so = 0.1, mu = 1.1, tau = 0.2)
#' plot(dp)
#' plot(dp, xlim = c(0, 5), length.out = 500)
#' @export
plot.designPrior <- function(x, ...) {
    if (!(methods::hasArg(xlim))) {
        xlim <- x$dpMean + c(-4, 4)*sqrt(x$dpVar)
    } else {
        xlim <- list(...)$xlim
    }
    if (!(methods::hasArg(length.out))) {
        length.out <- 1000
    } else {
        length.out <- list(...)$length.out
    }
    xseq <- seq(from = xlim[1], to = xlim[2], length.out = length.out)
    densFun <- density(x)
    plot(x = xseq, y = densFun(xseq), type = "l", xlab = "Effect size",
         ylab = "Design prior density", las = 1, ...)
}
