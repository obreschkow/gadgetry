#' Read Gadget simulation data
#'
#' @importFrom bit64 integer64
#'
#' @description Reads astrophysical N-body data output by the Gadget code (versions 1/2/3/4, see \url{https://wwwmpa.mpa-garching.mpg.de/gadget4/}). Binary and HDF5 formats can be read.
#'
#' @param file filename of snapshot to load. If the snapshot is split into several subvolumes, an asterix symbol (*) can be used in the filename at the place of the subvolume index. In this case, all the subvolumes are loaded and concatenated into a single snapshot with adjusted header information.
#' @param type character specifying the data format. Must be either of: \code{bin} for binary format, \code{hdf} for HDF5 format or \code{auto} to automatically determine the format from file extension.
#' @param bit64 logical flag specifying if 64-bit integer should be used, for example to handle large IDs. Only used with HDF5 files.
#'
#' @return Returns an object of class \code{snapshot}, which is a structured list that closely resembles the HDF5 format of Gadget (see \url{https://wwwmpa.mpa-garching.mpg.de/gadget4/}). HDF5 names are used in the Header, even when reading binary files (for name conversions, see Table 4 of \url{https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf}).
#'
#' @seealso \code{\link{writesnapshot}}
#'
#' @author Danail Obreschkow
#'
#' @export

readsnapshot = function(file, type='auto', bit64=FALSE) {

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
      subvol[[nfiles]] = .readsnapshot.single(fn, type, bit64)
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
    dat = .readsnapshot.single(fn, type, bit64)

  }

  class(dat) = 'snapshot'

  return(dat)

}

.readsnapshot.single = function(file, type, bit64) {

  # load file
  if (type=='bin') {

    # initialize
    block.start = function(target.block.size=NA) {
      block.size = readBin(data, 'integer', size=4)
      if (length(block.size)==0) {
        close(data)
        stop('End of file reached.')
      }
      #if (!is.na(target.block.size) && target.block.size>2^31-1) target.block.size=NA
      if (!is.na(target.block.size)) {
        if (target.block.size!=block.size) {
          cat(sprintf('WARNING: Expected block size (%.0f) differs from block size number in file (%.0f).\n',target.block.size,block.size))
        }
      }
      return(block.size)
    }
    block.end = function(block.size) {
      check = readBin(data, 'integer', size=4)
      if (check!=block.size) {
        cat('WARNING: Block termination error.\n')
      }
    }

    # open file
    data = file(file, "rb")
    dat = list()

    # read header
    block.size = block.start(256)
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
    block.end(block.size)

    # read positions
    n = sum(dat$Header$NumPart_ThisFile)
    block.size = block.start(12*n)
    for (i in seq(0,5)) {
      k = dat$Header$NumPart_ThisFile[i+1]
      if (k>0) {
        field = sprintf('PartType%d',i)
        dat[[field]] = list()
        dat[[field]]$Coordinates = t(array(readBin(data, 'numeric', 3*k, size = 4),c(3,k)))
      }
    }
    block.end(block.size)

    # read velocities
    block.size = block.start(12*n)
    for (i in seq(0,5)) {
      k = dat$Header$NumPart_ThisFile[i+1]
      if (k>0) {
        field = sprintf('PartType%d',i)
        dat[[field]]$Velocities = t(array(readBin(data, 'numeric', 3*k, size = 4),c(3,k)))
      }
    }
    block.end(block.size)

    # read IDs
    block.size = block.start(4*n)
    for (i in seq(0,5)) {
      k = dat$Header$NumPart_ThisFile[i+1]
      if (k>0) {
        field = sprintf('PartType%d',i)
        dat[[field]]$ParticleIDs = readBin(data, 'integer', k, size = 4)
      }
    }
    block.end(block.size)

    # read masses
    nmass = dat$Header$NumPart_ThisFile[dat$Header$MassTable==0]
    if (sum(nmass)>0) {
      block.size = block.start(4*sum(nmass))
      for (i in seq(0,5)) {
        if (nmass[i+1]>0) {
          field = sprintf('PartType%d',i)
          dat[[field]]$Masses = readBin(data, 'numeric', nmass[i+1], size = 4)
        }
      }
      block.end(block.size)
    }

    # read internal energies
    ngas = dat$Header$NumPart_ThisFile[1]
    if (ngas >0) {
      block.size = block.start(4*ngas)
      dat$PartType0$InternalEnergy = readBin(data, 'numeric', ngas, size = 4)
      block.end(block.size)
    }

    # close file
    close(data)

  } else if (type=='hdf') {

    if (!requireNamespace("rhdf5", quietly=TRUE)) {
      stop('Package rhdf5 is needed to load HDF5 data.')
    }

    groups = rhdf5::h5ls(file,recursive=FALSE)$name
    if ('ParticleType'%in%substring(groups,1,12)) {
      base = 'ParticleType'
    } else if ('PartType'%in%substring(groups,1,8)) {
      base = 'PartType'
    } else {
      stop('Failed to find particles groups "PartType#" or "ParticleType#" in the HDF file.')
    }

    vectormatrix = function(x) {
      if (length(dim(x))==2 & dim(x)[1]==3) {
        x = t(x)
      } else {
        if (length(x)%%3==0) {
          x = t(array(x,c(3,length(x)/3)))
        } else {
          stop('unknown hdf structure')
        }
      }
      return(x)
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
        if (bit64) {
          raw = rhdf5::h5read(file,field.file,bit64conversion='bit64')
        } else {
          raw = rhdf5::h5read(file,field.file)
        }
        if (!is.null(raw$Coordinates)) raw$Coordinates = vectormatrix(raw$Coordinates)
        if (!is.null(raw$Velocities)) raw$Velocities = vectormatrix(raw$Velocities)
        if (!is.null(raw$Accelerations)) raw$Accelerations = vectormatrix(raw$Accelerations)
        dat[[field]] = raw
      }
    }

  } else {

    stop('Unknown "type".')

  }

  class(dat) = 'snapshot'

  return(dat)

}
