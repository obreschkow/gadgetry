#' Write simulation snapshot in Gadget-2 binary format
#'
#' @description Write astrophysical N-body+SPH snapshot in Gadget-2 binary format, for use as initial conditions with Gadget Version 1/2/3/4 (SnapFormat=1).
#'
#' @param part is a list of initial conditions. It contains the items:\cr
#' \code{x} = N-by-3 matrix specifying the initial position in Cartesian coordinates of the N particles\cr
#' \code{v} = optional N-by-3 matrix specifying the initial velocities. If not given, velocities are set to zero.\cr
#' \code{id} = optional N-by-3 matrix specifying unique particle IDs. If not given, the particles will be numbered sequentially from 1 to N.\cr
#' \code{m} = optional vector giving the masses of the particles, whose mass is not already specified in \code{head$Massarr}. This is only used for variable particles masses. Fixed masses are best specified with the header argument \code{Massarr}.\cr
#' \code{u} = optional vector giving the internal energies of the SPH particles, i.e. the particles of type 0 (in C-indexing).\cr\cr
#'
#' @param head is a list of header arguments. For full details, please refer to the Gadget-2 manual (\url{https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf}, page 33). The most important items are:\cr
#' \code{Npart} = 6-vector giving the number of particles of each type. The first element corresponds to gas particles; the other five are treated as collision-less.\cr
#' \code{NMassarr} = optional 6-vector listing a fixed mass for each particle type. If set to 0 for a type which is present, individual particle masses must be provided in the vector part$m.\cr\cr
#'
#' @param file filename.
#'
#' @return Returns a list containing the particle data. The format of the list depends on the input format.
#'
#' @seealso \code{\link{readsnapshot}}
#'
#' @author Danail Obreschkow
#'
#' @export

writesnapshot = function (part, head, file) {

  # checks
  if (is.null(head$Npart)) stop('head$Npart must be provided')
  if (length(head$Npart)!=6) stop('head$Npart must be a vector of 6 integers')
  n = sum(head$Npart)
  n.sph = head$Npart[1]
  if (is.null(head$Massarr)) head$Massarr = rep(0,6)
  if (is.null(head$Nall)) head$Nall = head$Npart
  if (is.null(head$NallHW)) head$NallHW = rep(0,6)
  if (is.null(part$x)) stop('the matrix part$x must be provided.')
  if (length(part$x)!=3*n) stop('part$x must be a N-by-3 matrix.')
  if (dim(part$x)[1]!=n) stop('part$x must be a N-by-3 matrix.')
  if (is.null(part$v)) {
    part$v = array(0,c(n,3))
  } else {
    if (length(part$v)!=3*n) stop('part$v must be a N-by-3 matrix.')
    if (dim(part$v)[1]!=n) stop('part$v must be a N-by-3 matrix.')
  }
  if (is.null(part$id)) {
    part$id = seq(n)
  } else {
    if (length(part$id)!=n) stop('part$id msut be a N-vector.')
  }
  n.masses.needed = sum(head$Npart[head$Massarr==0])
  if (n.masses.needed==0) {
    if (!is.null(part$m)) stop('part$m must not be given if all masses are already provided in head$Massarr.')
  } else {
    if (is.null(part$m)) stop('part$m must be given if not all masses are given via head$Massarr.')
    if (length(part$m)!=n.masses.needed) stop('length of part$m inconsistent with head$Npart and head$Massarr.')
  }
  if (n.sph>0) {
    if (is.null(part$u)) stop('if gas is present (head$Npart[1]>0), part$u must be provided')
    if (length(part$u)!=n.sph) stop('part$u must be a vector of length head$Npart[1]')
  } else {
    if (!is.null(part$u)) stop('if no gas is present (head$Npart[1]=0), part$u must not be provided')
  }

  # initialize
  block.terminator = function(val) writeBin(as.integer(val), data, size=4)
  writeif = function(val, size) {
    if (is.null(val)) val = 0
    if (size==4) {
      val = as.integer(val)
    } else {
      val = as.numeric(val)
    }
    writeBin(val, data, size = size)
  }

  # open file
  data = file(file, "wb")

  # write header
  block.terminator(256)
  writeBin(as.integer(head$Npart), data, size = 4) # 24 bytes
  writeBin(as.numeric(head$Massarr), data, size = 8) # 48 bytes
  writeif(head$Time, size = 8)
  writeif(head$z, size = 8)
  writeif(head$FlagSfr, size = 4)
  writeif(head$FlagFeedback, size = 4)
  writeBin(as.integer(head$Nall), data, size = 4) # 24 bytes
  writeif(head$FlagCooling, size = 4)
  writeif(head$NumFiles, size = 4)
  writeif(head$BoxSize, size = 8)
  writeif(head$OmegaM, size = 8)
  writeif(head$OmegaL, size = 8)
  writeif(head$h, size = 8)
  writeif(head$FlagAge, size = 4)
  writeif(head$FlagMetals, size = 4)
  writeBin(as.integer(head$NallHW), data, size = 4) # 24 bytes
  writeif(head$flag_entr_ics, size = 4)
  writeBin(as.integer(rep(0, 15)), data, size = 4)
  block.terminator(256)

  # write positions
  block.terminator(12*n)
  writeBin(as.numeric(t(part$x)), data, size = 4)
  block.terminator(12*n)

  # write velocities
  block.terminator(12*n)
  writeBin(as.numeric(t(part$v)), data, size = 4)
  block.terminator(12*n)

  # write IDs
  block.terminator(4*n)
  writeBin(as.integer(part$i), data, size = 4)
  block.terminator(4*n)

  # write masses
  if (!is.null(part$m)) {
    block.terminator(4*length(part$m))
    writeBin(as.numeric(part$m), data, size = 4)
    block.terminator(4*length(part$m))
  }

  # write internal energies
  if (!is.null(part$u)) {
    block.terminator(4*length(part$u))
    writeBin(as.numeric(part$u), data, size = 4)
    block.terminator(4*length(part$u))
  }

  # close file
  close(data)
}
