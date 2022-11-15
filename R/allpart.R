#' Concatenate identical properties of different Gadget particle species
#'
#' @description Concatenate identical properties of different particle species in a Gadget snapshot
#'
#' @param snapshot an object of class 'snapshot' (or an analogously structured list), which contains a Gadget simulation snapshot. Such objects are created, for example, using \code{\link{readsnapshot}} or \code{\link{snapshot}}. Gadget data type, which must contain at least one sublist PartType# (with #=0,1,...). Each sublist PartType# represents one type of particle (e.g. gas, stars, dark matter) and must contain at least the particles coordinates in an N-by-3 matrix \code{Coordinates}. Other optional elements of PartType# are:\cr
#' @param field character specifying particle property to concatenate.
#' @param species vector of particle species to be included.
#'
#' @return Returns an N-vector or an N-by-3 matrix representing the property "field" of all the particle species given in \code{species} (0-indexed).
#'
#' @author Danail Obreschkow
#'
#' @export

allpart = function(snapshot, field = 'Coordinates', species = seq(0,5)) {

  # determine required space
  n.valid.species = 0
  for (i in species) {
    group = sprintf('PartType%d',i)
    if (!is.null(snapshot[[group]])) {
      if (!is.null(snapshot[[group]][[field]])) {
        n.valid.species = n.valid.species+1
        dm = dim(snapshot[[group]][[field]])
        d = length(dm)
      } else {
        stop(sprintf('Field %s does not exist for species %d.\n',field,i))
      }
    }
  }
  if (n.valid.species==0) stop(sprintf('Field %s not found in considered species.\n',field))

  # check data size
  if (d==2) {
    if (dm[2]!=3) stop('unknown data format')
  } else if (d!=0 & d!=1) {
    stop('unknown data format')
  }

  # write data into output array
  x = c()
  for (i in species) {
    group = sprintf('PartType%d',i)
    if (!is.null(snapshot[[group]])) {
      if (!is.null(snapshot[[group]][[field]])) {
        if (d==2) {
          x = rbind(x,snapshot[[group]][[field]])
        } else {
          x = c(x,snapshot[[group]][[field]])
        }
      } else {
        stop(sprintf('Field %s does not exist for species %d.\n',field,i))
      }
    }
  }

  return(x)
}
