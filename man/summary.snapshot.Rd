% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.snapshot}
\alias{summary.snapshot}
\title{Summarize content of Gadget snapshot}
\usage{
\method{summary}{snapshot}(object, ...)
}
\arguments{
\item{object}{an object of class 'snapshot', which contains a Gadget simulation snapshot. Such objects are created, for example, using \code{\link{readsnapshot}} or \code{\link{snapshot}}. Gadget data type, which must contain at least one sublist PartType# (with #=0,1,...). Each sublist PartType# represents one type of particle (e.g. gas, stars, dark matter) and must contain at least the particles coordinates in an N-by-3 matrix \code{Coordinates}. Other optional elements of PartType# are:\cr}

\item{...}{unused arguments, used for consistency with base function \code{\link[base]{summary}}.}
}
\value{
None.
}
\description{
Summarizes the content of 'snapshot' class objects.
}
\examples{
filename = system.file('galaxy1.hdf5', package='gadgetry')
sn = readsnapshot(filename)
summary(sn)

}
\author{
Danail Obreschkow
}
