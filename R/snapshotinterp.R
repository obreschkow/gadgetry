#' Interpolate between Gadget snapshots
#'
#' @description Interpolate/extrapolate between adjacent Gadget snapshots, building on the interpolation routine \code{\link{particleinterp}}. The ordering of the particles follows that of the nearest snapshot, and if some particles change their species between the two snapshots, the interpolated snapshot uses the species of the nearest snapshot.
#'
#' @param sn0 first gadget snapshot
#' @param sn1 second gadget snapshot
#' @param t0 time of first snapshot
#' @param t1 time of second snapshot
#' @param ti time of interpolated output snapshot. Values outside the range [t0,t1] lead to an extrapolation.
#' @param usevelocities logical flag. If TRUE, velocities are used to improve the interpolation of positions.
#' @param v0factor factor to be applied to velocities of sn0 to convert to length/time units of the snapshot coordinates and the t0/t1/ti times.
#' @param v1factor factor to be applied to velocities of sn1 to convert to length/time units of the snapshot coordinates and the t0/t1/ti times. If not provided, this will be set equal to \code{v1factor}.
#' @param afield optional acceleration field for leapfrog integration. This must be a function of a n-by-3 array representing n position vectors, which returns an n-by-3 array with the n acceleration vectors in the same length and time units as the snapshot coordinates and t0/t1/ti times.
#' @param dt optional time step used for leapfrog integration; only used if \code{afield} given. The default is \code{dt=abs(t1-t0)/20}.
#'
#' @return snapshot object with interpolated particle properties. The interpolation time \code{ti} is returned in 'Header$Time'.
#'
#' @author Danail Obreschkow
#'
#' @export

snapshotinterp = function(sn0,sn1,t0=0,t1=1,ti=0.5,
                          usevelocities=TRUE,v0factor=1,v1factor=NULL,
                          afield=NULL,dt=NULL) {

  if (is.null(v1factor)) v1factor=v0factor

  if (ti==t0) {

    sni = sn0

  } else if (ti==t1) {

    sni = sn1

  } else {

    # copy properties from nearest snapshot
    f = max(0,min(1,(ti-t0)/(t1-t0)))
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

    # determine extra fields
    fieldnames = c()
    for (species in seq(0,10)) {
      field = sprintf('PartType%d',species)
      fieldnames = c(fieldnames,names(sn0[[field]]),names(sn1[[field]]))
    }
    fieldnames = unique(fieldnames)
    fieldnames = fieldnames[!fieldnames%in%c('Coordinates','Velocities','ParticleIDs')]

    # complete missing masses and internal energies
    if ('Masses'%in%fieldnames) {
      for (species in seq(0,10)) {
        field = sprintf('PartType%d',species)
        if (!is.null(sni[[field]])) {
          if (is.null(sn0[[field]]$Masses)) sn0[[field]]$Masses=rep(NA,dim(sn0[[field]]$Coordinates)[1])
          if (is.null(sn1[[field]]$Masses)) sn1[[field]]$Masses=rep(NA,dim(sn1[[field]]$Coordinates)[1])
        }
      }
    }
    if ('InternalEnergy'%in%fieldnames) {
      for (species in seq(0,10)) {
        field = sprintf('PartType%d',species)
        if (!is.null(sni[[field]])) {
          if (is.null(sn0[[field]]$InternalEnergy)) sn0[[field]]$InternalEnergy=rep(NA,dim(sn0[[field]]$Coordinates)[1])
          if (is.null(sn1[[field]]$InternalEnergy)) sn1[[field]]$InternalEnergy=rep(NA,dim(sn1[[field]]$Coordinates)[1])
        }
      }
    }

    # extract ordered positions and velocities
    if (master==0) {
      x0 = allpart(sn0,'Coordinates')
      x1 = allpart(sn1,'Coordinates')[indexslave,]
      v0 = allpart(sn0,'Velocities')*v0factor
      v1 = allpart(sn1,'Velocities')[indexslave,]*v1factor
      if ('Masses'%in%fieldnames) {
        m0 = allpart(sn0,'Masses')
        m1 = allpart(sn1,'Masses')[indexslave]
      }
      if ('InternalEnergy'%in%fieldnames) {
        u0 = allpart(sn0,'InternalEnergy')
        u1 = allpart(sn1,'InternalEnergy')[indexslave]
      }
    } else {
      x0 = allpart(sn0,'Coordinates')[indexslave,]
      x1 = allpart(sn1,'Coordinates')
      v0 = allpart(sn0,'Velocities')[indexslave,]*v0factor
      v1 = allpart(sn1,'Velocities')*v1factor
      if ('Masses'%in%fieldnames) {
        m0 = allpart(sn0,'Masses')[indexslave]
        m1 = allpart(sn1,'Masses')
      }
      if ('InternalEnergy'%in%fieldnames) {
        u0 = allpart(sn0,'InternalEnergy')[indexslave]
        u1 = allpart(sn1,'InternalEnergy')
      }
    }
    if ('Masses'%in%fieldnames) {
      m0[is.na(m0)] = m1[is.na(m0)]
      m1[is.na(m1)] = m0[is.na(m1)]
    }
    if ('InternalEnergy'%in%fieldnames) {
      u0[is.na(u0)] = u1[is.na(u0)]
      u1[is.na(u1)] = u0[is.na(u1)]
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
        if ('Masses'%in%fieldnames) {
          sni[[field]]$Masses = (1-f)*m0[n+seq(ntype)]+f*m1[n+seq(ntype)]
        }
        if ('InternalEnergy'%in%fieldnames) {
          sni[[field]]$InternalEnergy = exp((1-f)*log(u0[n+seq(ntype)])+f*log(u1[n+seq(ntype)]))
        }
        n = n+ntype
      }
    }
    if (n!=dim(p$x)[1]) stop('unidentified interpolation indexing error')

  }

  sni$Header$Time = ti

  class(sni) = "snapshot"
  return(sni)

}
