#' Summarize content of Gadget snapshot
#'
#' @importFrom cooltools runif3
#' @importFrom stats median
#' @importFrom methods is
#'
#' @description Summarizes the content of 'snapshot' class objects.
#'
#' @param object an object of class 'snapshot', which contains a Gadget simulation snapshot. Such objects are created, for example, using \code{\link{readsnapshot}} or \code{\link{snapshot}}. Gadget data type, which must contain at least one sublist PartType# (with #=0,1,...). Each sublist PartType# represents one type of particle (e.g. gas, stars, dark matter) and must contain at least the particles coordinates in an N-by-3 matrix \code{Coordinates}. Other optional elements of PartType# are:\cr
#' @param ... unused arguments, used for consistency with base function \code{\link[base]{summary}}.
#'
#' @return None.
#'
#' @author Danail Obreschkow
#'
#' @examples
#' filename = system.file('galaxy1.hdf5', package='gadgetry')
#' sn = readsnapshot(filename)
#' summary(sn)
#'
#' @method summary snapshot
#' @export
summary.snapshot = function(object, ...) {

  ntot = 0
  for (type in seq(0,10)) {
    field = sprintf('PartType%d',type)
    if (!is.null(object[[field]])) {
      if (is.null(object[[field]]$Coordinates)) stop('Coordinates missing for some particles.')
      if (dim(object[[field]]$Coordinates)[2]!=3) stop('Incorrect array size for some particle coordinates.')
      n = dim(object[[field]]$Coordinates)[1]
      ntot = ntot+n
      cat(sprintf('Particle Type %d:\n',type))
      cat(sprintf('  Number of particles = %d\n',n))
      cat(sprintf('  Available properties = %s\n',paste(names(object[[field]]),collapse=', ')))
      for (j in seq(names(object[[field]]))) {
        subfield = names(object[[field]])[j]
        x = object[[field]][[subfield]]
        if (all(x==round(x))) {
          cat(sprintf('  Range of %s = [%.0f,%.0f]\n',subfield,min(x),max(x)))
        } else {
          cat(sprintf('  Range of %s = [%.5e,%.5e]\n',subfield,min(x),max(x)))
        }
      }
      x = object[[field]]$Coordinates
      x0 = apply(x,2,mean)
      cat(sprintf('  Geometric center = (%.3e, %.3e, %.3e)\n',x0[1],x0[2],x0[3]))
      r2 = rowSums(sweep(x,2,x0)^2)
      cat(sprintf('  Distance from centre = %.3e (median), %.3e (max)\n',sqrt(stats::median(r2)),sqrt(max(r2))))
    }
  }
  cat(sprintf('Total num particles = %d\n',ntot))

}
