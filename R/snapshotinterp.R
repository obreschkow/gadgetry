#' Interpolate between Gadget snapshots
#'
#' @importFrom pointr ptr
#'
#' @description Interpolates between adjacent Gadget snapshots. Coordinates and velocities are interpolated using the low-level routine \code{\link{particleinterp}}. All other particle properties are interpolated linearly, if and only if the particles exist in both snapshots and are of the same type in both snapshots. The ordering and the types of the particles in the interpolated snapshot follows that of the nearest snapshot.
#'
#' @param sn0 first gadget snapshot; it must have at least one particle type, and each particle type must have at least the properties \code{ParticleIDs}, \code{Coordinates} and \code{Velocities} (if \code{usevelocities=TRUE})
#' @param sn1 second gadget snapshot; should have same format as sn0 and must contain exactly the same particle IDs
#' @param f fractional interpolation time: 0 corresponds to time t0, 1 to time t1; is automatically truncated to the interval [0,1]
#' @param t0 time of first snapshot, if not given it is extracted from the header of sn0
#' @param t1 time of second snapshot, if not given it is extracted from the header of sn1
#' @param usevelocities logical flag. If TRUE, velocities are used to improve the interpolation of positions.
#' @param v0factor factor to be applied to velocities of sn0 to convert to length units of the snapshot coordinates and the time units of t0 and t1.
#' @param v1factor factor to be applied to velocities of sn1 to convert to length units of the snapshot coordinates and the time units of t0 and t1.
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

  if (f<=0) {

    pointr::ptr('snmaster','sn0')

  } else if (f>=1) {

    pointr::ptr('snmaster','sn1')

  } else {

    # determine which snapshot is considered the master that
    # defines the ordering of particles and their species
    master = ifelse(f<0.5,0,1)
    if (F) snmaster = snslave = NULL# avoids warning in package check
    pointr::ptr('snmaster',sprintf('sn%d',master))
    pointr::ptr('snslave',sprintf('sn%d',1-master))

    # get start and end times
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
    if (t0>t1) stop('Time t0 must not be larger than t1.')

    # interpolate time
    ti = t0*(1-f)+t1*f
    snmaster$Header$Time = ti

    # extract coordinates and velocities into simple arrays
    x0 = allpart(sn0,'Coordinates')
    x1 = allpart(sn1,'Coordinates')
    v0 = allpart(sn0,'Velocities')*v0factor
    v1 = allpart(sn1,'Velocities')*v1factor

    # sort particles of slave in order of master
    idmaster = allpart(snmaster,'ParticleIDs')
    idslave = allpart(snslave,'ParticleIDs')
    if (length(idmaster)!=length(idslave)) stop('both snapshots must have the same number of particles in current implementation')
    indexslave = order(idslave)[idmaster]
    if (any(idmaster!=idslave[indexslave])) stop('the two snapshots do not contain the same particle IDs')
    if (master==0) {
      x1 = x1[indexslave,]
      v1 = v1[indexslave,]
    } else {
      x0 = x0[indexslave,]
      v0 = v0[indexslave,]
    }

    # interpolate positions & velocities
    if (usevelocities) {
      p = particleinterp(x0=x0,x1=x1,t0=t0,t1=t1,ti=ti,v0=v0,v1=v1,afield=afield,dt=dt)
    } else {
      p = particleinterp(x0=x0,x1=x1,t0=t0,t1=t1,ti=ti)
    }

    # scale velocities back to original units
    vifactor = (1-f)*v0factor+f*v1factor
    p$v = p$v/vifactor

    # replace positions and velocities for interpolated values
    # and interpolate all other properties
    n = 0
    for (species in seq(0,10)) {
      field = sprintf('PartType%d',species)
      if (!is.null(snmaster[[field]])) {
        imaster = snmaster[[field]]$ParticleIDs
        ntype = length(imaster)
        if (ntype>0) {

          # check number of particles against header
          if (!is.null(snmaster$Header$NumPart_ThisFile)) {
            if (species<length(snmaster$Header$NumPart_ThisFile)) {
              if (ntype!=snmaster$Header$NumPart_ThisFile[species+1]) {
                stop('Header$NumPart_ThisFile does not seem to match the particle data')
              }
            } else {
              stop('Header$NumPart_ThisFile has not enough elements')
            }
          }

          # check particle IDs
          if (any(imaster!=idmaster[n+seq(ntype)])) stop('unidentified particle ID error')

          # copy interpolated positions and velocities
          snmaster[[field]]$Coordinates = p$x[n+seq(ntype),]
          snmaster[[field]]$Velocities = p$v[n+seq(ntype),]

          # interpolate all other properties
          if (!is.null(snslave[[field]])) {
            islave = snslave[[field]]$ParticleIDs
            if (length(islave)>0) {
              g = ifelse(f<0.5,f,1-f)
              if (length(islave)!=length(imaster) || !all(islave==imaster)) {
                imaster = c(5,34,9,1,71,32,6)
                islave = c(34,71,6,5,178,2,1)
                im = which(imaster%in%islave)
                im = im[order(imaster[im])]
                is = which(islave%in%imaster)
                is = is[order(islave[is])]
                # => islave[is] = imaster[im] is a monotonically increasing vector of all matched particleIDs
              } else {
                im = is = 1:length(imaster)
              }
              if (length(imaster)>0) {
                fieldnames = names(snmaster[[field]])
                fieldnames = fieldnames[!fieldnames%in%c('ParticleIDs','Coordinates','Velocities')]
                for (fieldname in fieldnames) {
                  if (!is.null(snslave[[field]][[fieldname]])) {
                    if (length(dim(snmaster[[field]][[fieldname]]))==2) {
                      snmaster[[field]][[fieldname]][im,] = (1-g)*snmaster[[field]][[fieldname]][im,]+g*snslave[[field]][[fieldname]][is,]
                    } else {
                      snmaster[[field]][[fieldname]][im] = (1-g)*snmaster[[field]][[fieldname]][im]+g*snslave[[field]][[fieldname]][is]
                    }
                  }
                }
              }
            }
          }

          n = n+ntype

        }
      }
    }
    if (n!=dim(p$x)[1]) stop('unidentified interpolation indexing error')

  }

  class(snmaster) = "snapshot"
  return(snmaster)

}
