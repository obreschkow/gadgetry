#' Write Gadget-4 parameter file
#'
#' @importFrom utils write.table
#'
#' @description Write parameter file for the astrophysical N-body+SPH simulation code Gadget-4.
#'
#' @param file parameter filename
#'
#' @param param list of parameter with values
#'
#' @return Returns the full parameter list.
#'
#' @author Danail Obreschkow
#'
#' @export

writeparam = function(file, param=NULL) {

  out = .gadgetry.env$default.parameters

  if (!is.null(param)) {
    name = names(param)
    for (i in seq_along(param)) {
      out[name[i]] = param[[i]]
    }
  }

  line = {}
  line[1] = '% Gadget-4 parameterfile written with the nbody R-package'
  for (i in seq_along(out)) {
    line[i+1] = sprintf('%s %s',names(out)[i],out[[i]])
  }
  utils::write.table(line, file=file, col.names = FALSE, row.names = FALSE, quote = FALSE)

  invisible(out)

}
