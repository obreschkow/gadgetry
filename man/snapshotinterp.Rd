% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snapshotinterp.R
\name{snapshotinterp}
\alias{snapshotinterp}
\title{Interpolate between Gadget snapshots}
\usage{
snapshotinterp(sn0, sn1, ti, fraction = FALSE, afield = NULL, dt = NULL)
}
\arguments{
\item{sn0}{gadget snapshot}

\item{sn1}{optional second gadget snapshot}

\item{ti}{single number giving the time of the interpolated/extrapolated output positions and velocities}

\item{fraction}{logical flag (default \code{FALSE}). If \code{TRUE}, the number \code{ti} is interpreted as fractional time between the two snapshots, ranging linearly between 0 and 1.}

\item{afield}{optional acceleration field for leapfrog integration. This must be a function of a n-by-d array representing n position vectors, which returns an n-by-d array with the n acceleration vectors in the same length and time units as x and v.}

\item{dt}{optional time step used for leapfrog integration; only used if \code{afield} given. The default is \code{dt=abs(t1-t0)/20}.}
}
\value{
snapshot object with interpolated particle properties.
}
\description{
Interpolate/extrapolate between adjacent Gadget snapshots, building on the interpolation routine \code{\link{particleinterp}}. Note that only positions and velocities are interpolated, while all other properties are set equal to the nearest snapshot in time. The ordering of the particles also follows that of the nearest snapshot, and if some particles change their species between the two snapshots, the interpolated snapshot uses the species of the nearest snapshot.
}
\author{
Danail Obreschkow
}