#' @title Success region for replication effect estimate
#'
#' @description Creates a success region object which can then be used for
#'     computing the probability of replication success with \code{\link{pors}}.
#'
#' @param intervals A 2xN matrix containing N disjoint intervals, the first
#'     column containing the lower and the second column containing the upper
#'     limits
#'
#' @return Returns an object of class \code{"successRegion"} which is a matrix
#'     containing the success intervals sorted in ascending order
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
#' successRegion(intervals = rbind(c(1.96, Inf), c(-Inf, -1.96)))
#' successRegion(intervals = cbind(1.96, Inf))
#'
#' @export
successRegion <- function(intervals) {
    ## input checks
    stopifnot(
        is.matrix(intervals),
        !any(intervals[,1,drop = FALSE] > intervals[,2, drop = FALSE])
    )

    ## check whether intervals overlap
    intervalsSorted <- intervals[order(intervals[,2], decreasing = TRUE),,
                                 drop = FALSE]
    if (nrow(intervals) > 1) {
        for (i in seq(2, nrow(intervalsSorted))) {
            if (intervalsSorted[i - 1, 1] < intervalsSorted[i, 2]) {
                stop("intervals must be disjoint")
            }
        }
    }
    class(intervalsSorted) <- "successRegion"
    return(intervalsSorted)
}

#' Print method for class \code{"successRegion"}
#' @method print successRegion
#'
#' @param x Object of class \code{"successRegion"}
#' @param ... Other arguments
#'
#' @return Prints text summary in the console and invisibly returns the
#'     \code{"successRegion"} object
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## success region for two-sided significance test
#' successRegion(intervals = rbind(c(1.96, Inf), c(-Inf, -1.96)))
#' ## success region for one-sided significance test
#' successRegion(intervals = rbind(c(1.96, Inf)))
#' @export
print.successRegion <- function(x, ...) {
    parensMat <- matrix(nrow = nrow(x), ncol = ncol(x))
    parensMat[,1] <- ifelse(is.finite(x[,1]), "[", "(")
    parensMat[,2] <- ifelse(is.finite(x[,2]), "]", ")")
    intChar <- rev(paste0(paste0(parensMat[,1], x[,1], ", "),
                          paste0(x[,2], parensMat[,2])))
    cat("Success region for replication effect estimate\n")
    cat(paste0("  ", intChar), sep = "  and")
    cat("\n")
    invisible(x)
}
