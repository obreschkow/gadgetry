#' Interpolate between Gadget snapshots
#'
#' @description Interpolate/extrapolate between adjacent Gadget snapshots, building on the interpolation routine \code{\link{particleinterp}}. Note that only positions and velocities are interpolated, while all other properties are set equal to the nearest snapshot in time. The ordering of the particles also follows that of the nearest snapshot, and if some particles change their species between the two snapshots, the interpolated snapshot uses the species of the nearest snapshot.
#'
#' @param sn0 gadget snapshot
#' @param sn1 optional second gadget snapshot
#' @param ti single number giving the time of the interpolated/extrapolated output positions and velocities
#' @param fraction logical flag (default \code{FALSE}). If \code{TRUE}, the number \code{ti} is interpreted as fractional time between the two snapshots, ranging linearly between 0 and 1.
#' @param afield optional acceleration field for leapfrog integration. This must be a function of a n-by-d array representing n position vectors, which returns an n-by-d array with the n acceleration vectors in the same length and time units as x and v.
#' @param dt optional time step used for leapfrog integration; only used if \code{afield} given. The default is \code{dt=abs(t1-t0)/20}.
#'
#' @return snapshot object with interpolated particle properties.
#'
#' @author Danail Obreschkow
#'
#' @export

snapshotinterp = function(sn0,sn1,ti,fraction=FALSE,afield=NULL,dt=NULL) {

  # input checks - to make sure that all the parameters are given

  # fill in sn1 (extrapolation only)
  if (is.null(sn1)) sn1=sn0

  # interpolate time
  t0 = sn0$Header$Time
  t1 = sn1$Header$Time
  if (fraction) ti = t0+(t1-t0)*ti

  if (ti==t0) {

    sni = sn0

  } else if (ti==t1) {

    sni = sn1

  } else {

    # copy properties from nearest snapshot
    f = (ti-t0)/(t1-t0)
    if (f<0.5) {
      master = 0
      slave = 1
      sni = sn0
    } else {
      master = 1
      slave = 0
      sni = sn1
    }

    # sort particles of other snapshot in same order
    idmaster = allpart(sni,'ParticleIDs')
    if (slave==0) {
      idslave = allpart(sn0,'ParticleIDs')
    } else {
      idslave = allpart(sn1,'ParticleIDs')
    }

    k = sort.int(sort.int(idmaster,index.return=TRUE)$ix,index.return=TRUE)$ix
    indexslave = sort.int(idslave,index.return=TRUE)$ix[k]
    if (any(idmaster!=idslave[indexslave])) stop('the two snapshots do not contain the same particle IDs')

    # extract ordered positions and velocities
    if (master==0) {
      x0 = allpart(sn0,'Coordinates')
      x1 = allpart(sn1,'Coordinates')[indexslave,]
      v0 = allpart(sn0,'Velocities')
      v1 = allpart(sn1,'Velocities')[indexslave,]
    } else {
      x0 = allpart(sn0,'Coordinates')[indexslave,]
      x1 = allpart(sn1,'Coordinates')
      v0 = allpart(sn0,'Velocities')[indexslave,]
      v1 = allpart(sn1,'Velocities')
    }

    # interpolate positions & velocities
    p = particleinterp(x0=x0,x1=x1,t0=t0,t1=t1,ti=ti,v0=v0,v1=v1,afield=afield,dt=dt)

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

