#' Write Gadget-4 parameter file
#'
#' @importFrom utils write.table
#'
#' @description Write parameter file for the astrophysical N-body+SPH simulation code Gadget-4. Parameters that are not specified, will be taken from \code{.gadgetry.env$default.parameters}.
#'
#' @param file parameter filename
#'
#' @param param list of parameter with values
#'
#' @param use.defaults logical flag. If TRUE, default parameters as listed in \code{.gadgetry.env$default.parameters} are used, unless otherwise specified in \code{param}.
#'
#' @param suppress.checks logical flag. If TRUE, input parameter checks are suppressed.
#'
#' @return Returns the full parameter list.
#'
#' @author Danail Obreschkow
#'
#' @export

writeparam = function(file, param=NULL, use.defaults=TRUE, suppress.checks=FALSE) {

  if (use.defaults) {
    out = .gadgetry.env$default.parameters
    if (!is.null(param)) {
      name = names(param)
      for (i in seq_along(param)) {
        out[name[i]] = param[[i]]
      }
    }
  } else {
    out = param
  }

  line = {}
  line[1] = '% Gadget-4 parameterfile written with the nbody R-package'
  for (i in seq_along(out)) {
    line[i+1] = sprintf('%s %s',names(out)[i],out[[i]])
  }

  if (!suppress.checks) {
    if (out$ComovingIntegrationOn==0 & out$OmegaLambda>0) stop('If ComovingIntegrationOn=0, then OmegaLambda=0 is normally needed.')
  }


  utils::write.table(line, file=file, col.names = FALSE, row.names = FALSE, quote = FALSE)

  invisible(out)

}
