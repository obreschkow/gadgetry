#' Run Gadget simulation
#'
#' @importFrom cooltools tick tock
#'
#' @description Executes a pre-compiled Gadget code with custom parameter file.
#'
#' @param file.gadget path+filename of compiled Gadget-4 code
#' @param file.param path+filename of parameter file
#' @param n.cores number of cores to be used with MPI
#' @param verbose logical flag indicating whether to display the standard console output
#' @param measure.time logical flag indicating whether to measure the total simulation time
#'
#' @return Returns the console command used to call Gadget.
#'
#' @author Danail Obreschkow
#'
#' @export

run.gadget = function(file.gadget, file.param, n.cores=8, verbose=FALSE, measure.time=TRUE) {

  if (!file.exists(file.gadget)) stop(sprintf('File does not exist: %s',file.gadget))
  if (!file.exists(file.param)) stop(sprintf('File does not exist: %s',file.param))

  if (measure.time) cooltools::tick('Running Gadget')

  cmd = sprintf('mpirun -np %d %s %s',n.cores,file.gadget,file.param)
  error.code = system(cmd,intern=FALSE,ignore.stderr=TRUE,ignore.stdout=!verbose)
  if (error.code>0) stop('Could not complete the Gadget run.')

  if (measure.time) cooltools::tock()

  invisible(cmd)
}
