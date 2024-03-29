% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pors.R
\name{pors}
\alias{pors}
\title{Compute probability of replication success}
\usage{
pors(sregion, dprior, sr, nsites = 1)
}
\arguments{
\item{sregion}{Success region for replication effect estimate}

\item{dprior}{Design prior object}

\item{sr}{Standard error of replication effect estimate}

\item{nsites}{Number of sites, defaults to \code{1}. The sites are assumed to
have the same standard error sr}
}
\value{
The probability of replication success
}
\description{
This function computes the probabiliy of replication success
    based on a success region for the replication effect estimate, a design
    prior, and a replication standard error. If the specified number of sites
    is larger than 1, the supplied success region has to be formulated in
    terms of the meta-analytic replication effect estimate across sites.
}
\examples{
dprior <- designPrior(to = 1.1, so = 1)
sregion <- successRegion(intervals = cbind(1.96, Inf))
pors(sregion = sregion, dprior = dprior, sr = 1)

}
\references{
Pawel, S., Consonni, G., and Held, L. (2023). Bayesian approaches to
designing replication studies. Psychological Methods.
\doi{10.1037/met0000604}
}
\author{
Samuel Pawel
}
