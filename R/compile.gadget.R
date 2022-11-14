#' Compile Gadget simulation
#'
#' @importFrom cooltools tick tock
#'
#' @description Complies the Gadget code with custom parameter file.
#'
#' @param source path of the source code. This path must contain the Makefile ana src directory.
#' @param DIR path containing the standard configuration file (\code{Config.sh}) and executable (\code{Gadget4}, or similar) to be written.
#' @param CONFIG custom path+filename of configuration file.
#' @param EXEC custom path+filename of the executable.
#' @param n.cores number of cores to compile the code
#' @param clean logical flag indicating whether previous compilations should first be cleaned, if they exist. Make sure to specify the right directory (\code{DIR}).
#' @param verbose logical flag indicating whether to display the standard console output
#' @param measure.time logical flag indicating whether to measure the total compilation time
#'
#' @return Returns the console command used to compile Gadget.
#'
#' @author Danail Obreschkow
#'
#' @export

compile.gadget = function(source, DIR=NULL, CONFIG=NULL, EXEC=NULL, n.cores=8, clean=FALSE, verbose=FALSE, measure.time=TRUE) {

  if (!file.exists(source)) stop(sprintf('Path does not exist: %s',source))

  if (is.null(DIR)) {
    if (is.null(CONFIG)) stop('If DIR is not specified, CONFIG must be given.')
    if (is.null(EXEC)) stop('If DIR is not specified, EXEC must be given.')
  } else {
    if (!is.null(CONFIG)) stop('If DIR is specified, CONFIG must not be given.')
    if (!is.null(EXEC)) stop('If DIR is specified, EXEC must not be given.')
  }

  if (measure.time) cooltools::tick('Compile Gadget\n')

  if (clean) {
    if (is.null(DIR)) stop('DIR must be specified to clean previous builds')
    cmd = sprintf('make -C %s clean DIR=%s',source,DIR)
    error.code = system(cmd,intern=FALSE,ignore.stderr=FALSE,ignore.stdout=!verbose)
    if (error.code>0) stop('Could not clean previous builds.')
  }

  if (is.null(DIR)) {
    cmd = sprintf('make -C %s -j %d CONFIG=%s EXEC=%s',source,n.cores,CONFIG,EXEC)
  } else {
    cmd = sprintf('make -C %s -j %d DIR=%s',source,n.cores,DIR)
  }
  error.code = system(cmd,intern=FALSE,ignore.stderr=FALSE,ignore.stdout=!verbose)
  if (error.code>0) stop(sprintf('Could not compile the Gadget code using command:\n%s',cmd))

  if (measure.time) cooltools::tock('Compilation completed',fmt='%s (%.2fs).\n')

  invisible(cmd)
}
