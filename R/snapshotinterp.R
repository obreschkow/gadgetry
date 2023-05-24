#' Interpolate between Gadget snapshots
#'
#' @description Interpolate/extrapolate between adjacent Gadget snapshots, building on the interpolation routine \code{\link{particleinterp}}. The ordering of the particles follows that of the nearest snapshot, and if some particles change their species between the two snapshots, the interpolated snapshot uses the species of the nearest snapshot.
#'
#' @param sn0 first gadget snapshot
#' @param sn1 second gadget snapshot
#' @param f interpolation time: 0 corresponds to time t0, 1 to time t1; should be within [0,1]
#' @param t0 time of first snapshot, if not given it is extracted from the header of sn0
#' @param t1 time of second snapshot, if not given it is extracted from the header of sn1
#' @param usevelocities logical flag. If TRUE, velocities are used to improve the interpolation of positions.
#' @param v0factor factor to be applied to velocities of sn0 to convert to length/time units of the snapshot coordinates and the t0/t1/ti times.
#' @param v1factor factor to be applied to velocities of sn1 to convert to length/time units of the snapshot coordinates and the t0/t1/ti times.
#' @param afield optional acceleration field for leapfrog integration. This must be a function of a n-by-3 array representing n position vectors, which returns an n-by-3 array with the n acceleration vectors in the same length and time units as the snapshot coordinates and t0/t1/ti times.
#' @param dt optional time step used for leapfrog integration; only used if \code{afield} given. The default is \code{dt=abs(t1-t0)/20}.
#'
#' @return snapshot object with interpolated particle properties. The interpolation time \code{ti} is returned in 'Header$Time'.
#'
#' @author Danail Obreschkow
#'
#' @export

snapshotinterp = function(sn0, sn1, f=0.5, t0=NULL, t1=NULL,
                          usevelocities=TRUE, v0factor=1, v1factor=v0factor,
                          afield=NULL, dt=NULL) {

  if (f==0) {

    sni = sn0

  } else if (f==1) {

    sni = sn1

  } else {

    # get times
    if (is.null(t0)) {
      if (is.null(sn0$Header$Time)) {
        stop('sn0$Header$Time missing; provide value of t0')
      } else {
        t0 = sn0$Header$Time
      }
    }
    if (is.null(t1)) {
      if (is.null(sn1$Header$Time)) {
        stop('sn1$Header$Time missing; provide value of t1')
      } else {
        t1 = sn1$Header$Time
      }
    }

    # make master snapshot (sn[[master]]) and slave snapshot (sn[[slave]]), such that the interpolated sn is closer to the master
    sn = list(sn0,sn1)
    sn0 = sn1 = NULL
    master = ifelse(f<0.5,1,2)
    slave = 3-master

    # copy master snapshot to interpolated snapshot
    sni = sn[[master]]

    # sort particles of slave in order of master
    idmaster = allpart(sn[[master]],'ParticleIDs')
    idslave = allpart(sn[[slave]],'ParticleIDs')
    indexslave = order(idslave)[idmaster]
    if (any(idmaster!=idslave[indexslave])) stop('the two snapshots do not contain the same particle IDs')

    # extract coordinates and velocities into simple arrays
    if (master==1) {
      x0 = allpart(sn[[1]],'Coordinates')
      x1 = allpart(sn[[2]],'Coordinates')[indexslave,]
      v0 = allpart(sn[[1]],'Velocities')*v0factor
      v1 = allpart(sn[[2]],'Velocities')[indexslave,]*v1factor
    } else {
      x0 = allpart(sn[[1]],'Coordinates')[indexslave,]
      x1 = allpart(sn[[2]],'Coordinates')
      v0 = allpart(sn[[1]],'Velocities')[indexslave,]*v0factor
      v1 = allpart(sn[[2]],'Velocities')*v1factor
    }
    if (length(unique(c(dim(x0)[1],dim(x1)[1],dim(v0)[1],dim(v1)[1])))!=1) stop('Missing coordinates/velocities')

    # interpolate positions & velocities
    ti = t0*(1-f)+t1*f
    sni$Header$Time = ti
    if (usevelocities) {
      p = particleinterp(x0=x0,x1=x1,t0=t0,t1=t1,ti=ti,v0=v0,v1=v1,afield=afield,dt=dt)
    } else {
      p = particleinterp(x0=x0,x1=x1,t0=t0,t1=t1,ti=ti)
    }

    # scale velocities back to original units
    vifactor = (1-f)*v0factor+f*v1factor
    p$v = p$v/vifactor

    # replace positions and velocities for interpolated values
    n = 0
    for (species in seq(0,10)) {
      field = sprintf('PartType%d',species)
      if (!is.null(sni[[field]])) {
        ntype = dim(sni[[field]]$Coordinates)[1]
        if (ntype>0 & !is.null(sni$Header$NumPart_ThisFile)) {
          if (species<length(sni$Header$NumPart_ThisFile)) {
            if (ntype!=sni$Header$NumPart_ThisFile[species+1]) {
              stop('Header$NumPart_ThisFile does not seem to match the particle data')
            }
          } else {
            stop('Header$NumPart_ThisFile has not enough elements')
          }
        }
        sni[[field]]$Coordinates = p$x[n+seq(ntype),]
        sni[[field]]$Velocities = p$v[n+seq(ntype),]
        n = n+ntype
      }
    }
    if (n!=dim(p$x)[1]) stop('unidentified interpolation indexing error')

  }

  class(sni) = "snapshot"
  return(sni)

}
