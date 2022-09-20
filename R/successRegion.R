#' @title Success region for replication effect estimate
#'
#' @description Creates a success region object
#'
#' @param intervals A 2xN matrix containing N disjoint intervals, the first
#'     entry containing the lower and the second entry containing the upper
#'     limit
#'
#' @return A successRegion object
#'
#' @author Samuel Pawel
#'
#' @references
#'
#' Pawel, S., Consonni, G., and Held, L. (2022). Bayesian approaches to
#' designing replication studies. arXiv preprint.
#' \doi{10.48550/arXiv.XXXX.XXXXX}
#'
#' @seealso \code{\link{probRS}}
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
