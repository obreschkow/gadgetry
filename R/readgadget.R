#' Read Gadget simulation data
#'
#' @importFrom rhdf5 h5ls h5read h5readAttributes
#'
#' @description Reads astrophysical N-body data output by the Gadget code (versions 1/2/3/4, see \url{https://wwwmpa.mpa-garching.mpg.de/gadget4/}). Binary and HDF5 formats can be read.
#'
#' @param file filename of snapshot to load. If the snapshot is split into several subvolumes, an asterix symbol (*) can be used in the filename at the place of the subvolume index. In this case, all the subvolumes are loaded and concatenated into a single snapshot with adjusted header information.
#' @param type character specifying the data format. Must be either of: \code{bin} for binary format, \code{hdf} for HDF5 format or \code{auto} to automatically determine the format from file extension.
#'
#' @return Returns an object of class \code{gadget}, which is a structured list that closely resembles the HDF5 format of Gadget (see \url{https://wwwmpa.mpa-garching.mpg.de/gadget4/}). HDF5 names are used in the Header, even when reading binary files (for name conversions, see Table 4 of \url{https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf}).
#'
#' @seealso \code{\link{writegadget}}
#'
#' @author Danail Obreschkow
#'
#' @export

readgadget = function(file, type='auto') {

  fn = gsub('\\*','0',file)

  # check file
  if (!file.exists(fn)) stop(paste0('File not found: ',fn))
  if (file.access(fn,4)!=0) stop(paste0('No permission to read: ',fn))

  # handle type
  if (type=='binary') type = 'bin'
  if (type=='hdf5' | type=='h5') type = 'hdf'
  if (type=='auto') {

    # get file extension
    ext = strsplit(basename(fn), split="\\.")[[1]]
    ext = ext[length(ext)]

    # check file extension
    if (length(ext)>0) {
      if (ext%in%c('hdf','hdf5','h5')) {
        type = 'hdf'
      } else {
        type = 'bin'
      }
    } else {
      type = 'bin'
    }

  }

  if (grepl('\\*',file)) {

    # read all subvolumes
    subvol = list()
    nfiles = 0
    while (file.exists(fn)) {
      nfiles = nfiles+1
      subvol[[nfiles]] = .readgadget.single(fn, type)
      fn = gsub('\\*',sprintf('%d',nfiles),file)
    }

    # concatenate subvolumes
    dat = subvol[[1]]
    if (nfiles>1) {
      for (i in seq(2,nfiles)) {
        dat$Header$NumPart_ThisFile = dat$Header$NumPart_ThisFile+subvol[[i]]$Header$NumPart_ThisFile
        for (k in seq(0,5)) {
          field = sprintf('PartType%d',k)
          if (!is.null(dat[[field]])) {
            dat[[field]]$Coordinates = rbind(dat[[field]]$Coordinates,subvol[[i]][[field]]$Coordinates)
            dat[[field]]$Velocities = rbind(dat[[field]]$Velocities,subvol[[i]][[field]]$Velocities)
            dat[[field]]$ParticleIDs = c(dat[[field]]$ParticleIDs,subvol[[i]][[field]]$ParticleIDs)
            dat[[field]]$Masses = c(dat[[field]]$Masses,subvol[[i]][[field]]$Masses)
            dat[[field]]$InternalEnergy = c(dat[[field]]$InternalEnergy,subvol[[i]][[field]]$InternalEnergy)
          }
        }
      }
    }

    # check number of particles
    if (any(dat$Header$NumPart_ThisFile!=dat$Header$NumPart_Total)) stop('Particle number check failed.')

  } else {

    # read single file
    dat = .readgadget.single(fn, type)

  }

  class(dat) = 'gadget'

  return(dat)

}

.readgadget.single = function(file, type) {

  # load file
  if (type=='bin') {

    # initialize
    block.size = NA
    block.start = function(target=NA) {
      if (is.na(block.size)) {
        block.size <<- readBin(data, 'integer', size=4)
        if (target>2^23-1) target=NA
        if (!is.na(target)) {
          if (target!=block.size) {
            close(data)
            stop('block opening error')
          }
        }
      } else {
        close(data)
        stop('block.start must be called after block.end')
      }
    }
    block.end = function() {
      if (is.na(block.size)) {
        close(data)
        stop('block.end must be called after block.start')
      } else {
        check = readBin(data, 'integer', size=4)
        if (check!=block.size) {
          close(data)
          stop('block termination error')
        }
        block.size <<- NA
      }
    }

    # open file
    data = file(file, "rb")
    dat = list()

    # read header
    block.start(256)
    dat$Header = list()
    dat$Header$NumPart_ThisFile = readBin(data, 'integer', 6, size = 4)
    dat$Header$MassTable = readBin(data, 'numeric', 6, size = 8)
    dat$Header$Time = readBin(data, 'numeric', 1, size = 8)
    dat$Header$Redshift = readBin(data, 'numeric', 1, size = 8)
    dat$Header$Flag_Sfr = readBin(data, 'integer', 1, size = 4)
    dat$Header$Flag_Feedback = readBin(data, 'integer', 1, size = 4)
    dat$Header$NumPart_Total = readBin(data, 'integer', 6, size = 4)
    dat$Header$Flag_Cooling = readBin(data, 'integer', 1, size = 4)
    dat$Header$NumFilesPerSnapshot = readBin(data, 'integer', 1, size = 4)
    dat$Header$BoxSize = readBin(data, 'numeric', 1, size = 8)
    dat$Header$OmegaM = readBin(data, 'numeric', 1, size = 8)
    dat$Header$OmegaL = readBin(data, 'numeric', 1, size = 8)
    dat$Header$h = readBin(data, 'numeric', 1, size = 8)
    dat$Header$Flag_StellarAge = readBin(data, 'integer', 1, size = 4)
    dat$Header$Flag_Metals = readBin(data, 'integer', 1, size = 4)
    dat$Header$NumPart_Total_HW = readBin(data, 'integer', 6, size = 4)
    dat$Header$Flag_Entropy_ICs = readBin(data, 'integer', 1, size = 4)
    dat$Header$unused = readBin(data, 'integer', 15, size = 4)
    block.end()

    # read positions
    n = sum(dat$Header$Npart)
    block.start(12*n)
    for (i in seq(0,5)) {
      k = dat$Header$Npart[i+1]
      if (k>0) {
        field = sprintf('PartType%d',i)
        dat[[field]] = list()
        dat[[field]]$Coordinates = t(array(readBin(data, 'numeric', 3*k, size = 4),c(3,k)))
      }
    }
    block.end()

    # read velocities
    block.start(12*n)
    for (i in seq(0,5)) {
      k = dat$Header$Npart[i+1]
      if (k>0) {
        field = sprintf('PartType%d',i)
        dat[[field]]$Velocities = t(array(readBin(data, 'numeric', 3*k, size = 4),c(3,k)))
      }
    }
    block.end()

    # read IDs
    n = sum(dat$Header$Npart)
    block.start(4*n)
    for (i in seq(0,5)) {
      k = dat$Header$Npart[i+1]
      if (k>0) {
        field = sprintf('PartType%d',i)
        dat[[field]]$ParticleIDs = readBin(data, 'integer', k, size = 4)
      }
    }
    block.end()

    # read masses
    nmass = dat$Header$Npart[dat$Header$Massarr==0]
    if (sum(nmass)>0) {
      block.start(4*sum(nmass))
      for (i in seq(0,5)) {
        if (nmass[i+1]>0) {
          field = sprintf('PartType%d',i)
          dat[[field]]$Masses = readBin(data, 'numeric', nmass[i+1], size = 4)
        }
      }
      block.end()
    }

    # read internal energies
    ngas = dat$Header$Npart[1]
    if (ngas >0) {
      block.start(4*ngas)
      dat$PartType0$InternalEnergy = readBin(data, 'numeric', ngas, size = 4)
      block.end()
    }

    # close file
    close(data)

  } else if (type=='hdf') {

    groups = rhdf5::h5ls(file,recursive=FALSE)$name
    if ('ParticleType'%in%substring(groups,1,12)) {
      base = 'ParticleType'
    } else if ('PartType'%in%substring(groups,1,8)) {
      base = 'PartType'
    } else {
      stop('Failed to find particles groups "PartType#" or "ParticleType#" in the HDF file.')
    }

    dat = list()
    if (length(groups)<1) stop('something wrong with HDF5 file')
    if ('Config'%in%groups) dat$Config = rhdf5::h5readAttributes(file,'/Config')
    if ('Header'%in%groups) dat$Header = rhdf5::h5readAttributes(file,'/Header')
    if ('Parameters'%in%groups) dat$Parameters = rhdf5::h5readAttributes(file,'/Parameters')
    for (i in seq(0,5)) {
      field.file = sprintf('%s%d',base,i)
      field = sprintf('PartType%d',i)
      if (field%in%groups) {
        raw = rhdf5::h5read(file,field.file)
        if (!is.null(raw$Coordinates)) raw$Coordinates = t(raw$Coordinates)
        if (!is.null(raw$Velocities)) raw$Velocities = t(raw$Velocities)
        if (!is.null(raw$Accelerations)) raw$Accelerations = t(raw$Accelerations)
        dat[[field]] = raw
      }
    }

  } else {

    stop('Unknown "type".')

  }

  class(dat) = 'snapshot'

  return(dat)

}
