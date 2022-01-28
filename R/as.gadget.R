#' Create GADGET object from particle data
#'
#' @description Converts raw particle data into a \code{gadget} object.
#'
#' @param x list of N-by-3 matrices specifying the coordinates of particles of various types
#' @param v optional list of N-by-3 matrices with particle velocities
#' @param id optional list of N-vectors with particle IDs
#' @param m optional list of N-vectors with particle masses
#' @param types optional vector of particle indices. The number of elements in this vector must be the same as the number of matrices in the list \code{x}.
#'
#' @return Returns an object of class \code{gadget}, which is a structured list that closely resembles the HDF5 format of Gadget (see \url{https://wwwmpa.mpa-garching.mpg.de/gadget4/}).
#'
#' @author Danail Obreschkow
#'
#' @export
as.gadget = function(x, v=NULL, m=NULL, id=NULL, types=NULL) {

  ntypes = length(x)

  if (is.null(types)) {
    types=seq(0,ntypes-1)
  } else {
    if (length(types)!=ntypes) stop('length of "types" must match the number of matrices in "x".')
    if (min(types)<0) stop('all particle types must be 0,1,2,3,4,5')
    if (max(types)>5) stop('all particle types must be 0,1,2,3,4,5')
  }

  dat = list(Header=list(NumPart_ThisFile=rep(0,6)))

  dat = list()
  for (i in seq(ntypes)) {
    field = sprintf('PartType%d',types[i])
    dat[[field]] = list(Coordinates=x[[i]])
    dat$Header$NumPart_ThisFile[types[i]+1] = dim(x[[i]])[1]
  }

  if (!is.null(v)) {
    if (length(v)!=ntypes) stop('"v" and "x" must contain the same number of elements.')
    for (i in seq(ntypes)) {
      field = sprintf('PartType%d',types[i])
      dat[[field]]$Velocities = v[[i]]
    }
  }

  if (!is.null(id)) {
    if (length(id)!=ntypes) stop('"id" and "x" must contain the same number of elements.')
    for (i in seq(ntypes)) {
      field = sprintf('PartType%d',types[i])
      dat[[field]]$ParticleIDs = m[[i]]
    }
  }

  if (!is.null(m)) {
    if (length(m)!=ntypes) stop('"m" and "x" must contain the same number of elements.')
    for (i in seq(ntypes)) {
      field = sprintf('PartType%d',types[i])
      dat[[field]]$Masses = m[[i]]
    }
  }

  class(dat) = 'gadget'

  return(dat)

}
