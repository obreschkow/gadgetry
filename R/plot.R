#' Visualise one or several 3D point sets
#'
#' @importFrom cooltools lim griddata kde2 quadrupole rotation3 vectornorm nplot rasterflip lightness
#' @importFrom png writePNG
#' @importFrom grDevices pdf dev.off col2rgb rainbow
#' @importFrom graphics axis lines par rasterImage rect text arrows
#'
#' @description Produces a raster image visualizing particle positions in 3D N-body/SPH data.
#'
#' @param x object of class 'snapshot'. It must contain at least one sublist PartType# (with #=0,1,...). Each sublist PartType# represents one type of particles (e.g. gas, dark matter, stars) and must contain at least the particles coordinates in an N-by-3 matrix \code{Coordinates}. Other optional elements of PartType# are:\cr\cr
#' \code{Masses} = N-vector with particle masses, used as luminosity weights\cr\cr
#' \code{Density} or \code{Densities} = N-vector with local densities. These are only used to determine the smoothing if kde=4 and if Masses are provided.\cr\cr
#' \code{col} = color used to display this particle type. This can be a single color or a vector of colors. If a single color is provided, a range of brightness of that color is produced automatically.\cr\cr
#' \code{value} = optional N-vector of values, one for each particle, which define the position on the color scale. If not given, particles are colored according to their projected density.\cr\cr
#' \code{valrange} = optional 2-vector specifying the values corresponding to the limits of the color scale. Only used if the N-vector \code{value} is given.\cr\cr
#' \code{density.scaling} = logical flag (default TRUE), specifying whether the brightness is scaled with the particle density if \code{value} is given.\cr\cr
#' \code{lum} = overall scaling factor for the color scale.\cr\cr
#' \code{gamma} = gamma parameter setting the non-linear conversion between density and brightness/color. The larger the value, the higher the dynamic range.\cr\cr
#' \code{kde} = integer, setting the density estimation method. 0: Simple gridding, 1: Gaussian blur, 2: Fast approximate kernel density estimation relying on 2D densities, 3: More accurate 3D-kernel-density estimation, 4: Adaptive smoothing based on densities of particles (requires Density and Masses to be provided for the particle type, otherwise the 3rd algorithm will be used as default).\cr\cr
#' \code{smoothing} = linear smoothing factor (default 1) to adjust the automatically determined size of the (adaptive or fixed) smoothing kernels.\cr\cr
#' \code{ref.density} = optional number specifying a fixed reference density (number of particles per area, in units of 1/(length unit)^2) used for luminosity normalization; fixing this value allows to suppress brightness fluctuations between snapshots (e.g. in movies).\cr
#' @param types vector specifying the particle types to be displayed. These numbers must correspond to the # in the sublists PartType#. If not given, all particle species are shown.
#' @param center optional 3-vector specifying the coordinate at the center of the plot. The default is the geometric center (= center of mass, if all particle masses are equal).
#' @param rotation either an integer (1-6), a 3-vector or a 3-by-3 matrix, specifying a rotation of the 3D particle positions. In case of an integer: 1=(x,y)-plane, 2=(y,z)-plane, 3=(x,z)-plane, 4=(qmax,qmin)-plane, 5=(qmax,qmid)-plane, 6=(qmid,qmin)-plane, where qmax/qmid/qmin are the eigenvectors of the particle-quadrupole, associated with the maximum/middle/minimum eigenvalues, respectively. If \code{rotation} is a vector, its direction specifies the rotation axis and its norm the rotation angle in radians in the positive geometric sense.
#' @param rot.center optional 3-vector specifying the center of rotation. If not given, this is assumed to be equal to the centre of the plot.
#' @param width optional horizontal range of the image in simulation units. The default corresponds to the full range of particle positions. If a field-of-view is specified by the parameter \code{fov}, the width of the view is superseded by the latter and the value of \code{width} is used to normalised the smoothing lengths.
#' @param fov optional field-of-view in degrees. If \code{NULL}, an orthogonal projection is applied. If specified, a stereographic projection is applied, such that the width of the image shows exactly the number of degrees specified by \code{fov}.
#' @param depth optional depth of the slice/field to be shown. In the case of an orthogonal projection (fov=NULL) all particles that lie within -depth/2 to +depth/2 are shown. In the case of a stereographic projection, all particles within a distance \code{depth} from the observer (placed at \code{center}) are shown.
#' @param taper logical flag, which is only used if a finite \code{depth} is specified. If \code{taper} is TRUE, the particles are faded out towards the edges of the depth range, using a cubic spline kernel (Monaghan 1992), whose characteristic width corresponds to \code{depth}.
#' @param aspect aspect ratio (= width/height).
#' @param npixels number of grid cells per side in the output image. If the image is not square npixels is interpreted as the geometric mean between the horizontal and the vertical number of pixels, such that the total number of pixels remains about npixels^2.
#' @param hdcolors logical flag. If TRUE, input color scales are slightly smoothed to avoid spurious slight color jumps steps due to 8-bit color representation.
#' @param lum default value used in \code{dat$PartType#} (see above).
#' @param gamma default value used in \code{dat$PartType#} (see above).
#' @param kde default value used in \code{dat$PartType#} (see above).
#' @param smoothing default value used in \code{dat$PartType#} (see above).
#' @param shadows overall contrast value.
#' @param sample.fraction fraction of particles to be used. If 1, all particles are used, if <1 a random subsample is drawn.
#' @param screen logical flag specifying whether the image is displayed on the screen.
#' @param pngfile optional png-filename to save the image as raster image.
#' @param pdffile optional pdf-filename to save the image as pdf-file.
#' @param title Text to be added to the figure.
#' @param title.origin optional 2-vector specifying the position of the title
#' @param arrows logical flag, specifying if axis-arrows are drawn
#' @param arrow.origin optional 2-vector specifying the bottom left corner of the arrows
#' @param arrow.length optional number specifying the length of the arrows
#' @param arrow.lwd line width of arrows
#' @param scale logical flag, specifying if a length scale is shown
#' @param scale.origin optional 2-vector specifying the right end of the length scale
#' @param scale.length optional number specifying the length of the length scale (is always rounded to one significant digit)
#' @param scale.lwd line width of length scale
#' @param length.unit character string of length unit (e.g. "m" or "kpc")
#' @param xlab label on horizontal arrow (only shown if \code{arrows=TRUE}).
#' @param ylab label on vertical arrow (only shown if \code{arrows=TRUE}).
#' @param cex text size
#' @param text.offset optional extra offset for the arrow and scale labels
#' @param text.col color of text, arrows and scale
#' @param ... additional parameters for \code{\link[graphics]{plot}}
#'
#' @return Returns a structured list of parameters and arrays useful to analyse and reproduce the figure.
#'
#' @seealso \code{\link{plot4}}
#'
#' @author Danail Obreschkow
#'
#' @examples
#' ## EXAMPLE 1: COLIBRE GALAXY COLOURED BY GAS IONIZATION FRACTION
#' # load a test galaxy from the COLIBRE simulation (L25m5)
#' filename = system.file('galaxy2.hdf5', package='gadgetry')
#' sn = readsnapshot(filename)
#'
#' # set gas rendering parameters, such that it is colored by the
#' # hydrogen ionization fraction (HII/H)
#' sn$PartType0$gamma = .8
#' sn$PartType0$value = sn$PartType0$SpeciesFractions[,3] # hydrogen ionization fraction
#' sn$PartType0$col = colorRampPalette(c('#1117ff','#ff0533'))(256)
#' sn$PartType0$smoothing = .7
#'
#' # set star rendering parameters
#' sn$PartType4$col = 'white'
#' sn$PartType4$lum = 0.7
#' sn$PartType4$gamma = .8
#'
#' # generate plot
#' par(mar=c(0,0,0,0))
#' plot(sn, center=0, width=0.07, rotation=c(0,1,0.3), npixels=500, arrows=FALSE)
#'
#' ## EXAMPLE 2: COLIBRE GALAXY COLOURED BY GAS TEMPERATURE
#' # load a test galaxy from the COLIBRE simulation (L25m5)
#' filename = system.file('galaxy2.hdf5', package='gadgetry')
#' sn = readsnapshot(filename)
#'
#' # set gas rendering parameters, such that it is colored by temperature
#' sn$PartType0$gamma = .8
#' sn$PartType0$value = log10(sn$PartType0$Temperatures)
#' sn$PartType0$valrange = c(2.5,4.2)
#' sn$PartType0$col = rev(rainbow(100,end=2/3))
#'
#' # generate plot, showing only the gas
#' par(mar=c(0,0,0,0))
#' plot(sn, center=0, width=0.07, rotation=c(0,1.6,0.55), types=0, npixels=500, arrows=FALSE)
#'
#' ## EXAMPLE 3: LOW-RES TEST GALAXY WITH DARK MATTER
#' # load galaxy
#' filename = system.file('galaxy1.hdf5', package='gadgetry')
#' sn = readsnapshot(filename)
#'
#' # rendering parameters of dark matter
#' sn$PartType1$col = '#0511ff'
#' sn$PartType1$smoothing = 2
#'
#' # rendering parameters of stars
#' sn$PartType2$col = '#ffff44'
#' sn$PartType2$smoothing = 0.7
#' sn$PartType2$lum = 0.7
#' sn$PartType2$gamma = 0.8
#'
#' # plot snapshot
#' out = plot(sn, length.unit='kpc', width=500, kde=3, rotation=c(0,0,0.5), arrows=FALSE)
#'
#' # add legend
#' col = cooltools::lightness(allpart(out, 'col'), 0.65)
#' legend(out$header$xlim[1],out$header$ylim[2],c('Dark matter','Stars'),
#'        col=col,pch=16,text.col=col,bty='n')
#'
#' @method plot snapshot
#' @export
plot.snapshot = function(x, center=NULL, rotation=1, rot.center=NULL, width=NULL, fov=NULL, depth=NULL, taper=FALSE,
                         aspect=1, npixels=300, kde=4, smoothing=1,
                         types=NULL, sample.fraction=1,
                         lum=1, gamma=1, shadows=1, hdcolors=TRUE,
                         screen=TRUE, pngfile=NULL, pdffile=NULL,
                         title=NULL, title.origin = NULL,
                         arrows = TRUE, arrow.origin = NULL, arrow.length = NULL, arrow.lwd = 1.5,
                         scale = TRUE, scale.origin = NULL, scale.length = NULL, scale.lwd = 1.5, length.unit = '',
                         xlab = NULL, ylab = NULL, cex=1, text.offset = 0, text.col = 'white', ...) {

  snapshot = x; x = NULL # needed, because the first argument of the generic plot function is "x"

  out = list()

  # determine particle types to be considered
  all.types = c()
  for (type in seq(0,20)) {
    field = sprintf('PartType%d',type)
    if (!is.null(snapshot[[field]]$Coordinates)) {
      all.types=c(all.types,type)
    }
  }
  if (is.null(types)) {
    types = all.types
  } else {
    types = sort(unique(types[types%in%all.types]))
  }

  # handle fov
  if (!is.null(fov)) {
    if (fov>120) stop('fov cannot be larger than 120 degrees.')
    if (fov<=0) stop('fov must be larger than 0.')
  }

  # check and pre-process input arguments ######################################
  if (min(types)<0) stop('particle types must be non-negative integers')

  for (type in types) {

    field = sprintf('PartType%d',type)
    out[[field]] = list()

    # handle brightness parameters
    if (is.null(snapshot[[field]]$lum)) snapshot[[field]]$lum=lum
    if (is.null(snapshot[[field]]$gamma)) snapshot[[field]]$gamma=gamma

    # specify smoothing kernel
    if (is.null(snapshot[[field]]$smoothing)) snapshot[[field]]$smoothing = smoothing
    if (is.null(snapshot[[field]]$kde)) snapshot[[field]]$kde = kde

    # convert "Densities" (used in SWIFT) to "Density" (used in Gadget)
    names(snapshot[[field]])[names(snapshot[[field]])=="Densities"]="Density"

    # handle values
    snapshot[[field]]$color.by.property = !is.null(snapshot[[field]]$value)
    if (snapshot[[field]]$color.by.property) {
      if (is.null(snapshot[[field]]$valrange)) snapshot[[field]]$valrange=range(snapshot[[field]]$value)
    } else {
      snapshot[[field]]$valrange = c(0,1)
    }
    out[[field]]$valrange = snapshot[[field]]$valrange

    # determine particle color
    if (is.null(snapshot[[field]][['col']])) {
      if (snapshot[[field]]$color.by.property) {
        snapshot[[field]][['col']] = grDevices::rainbow(256,end=5/6)
      } else {
        if (type<=5) {
          snapshot[[field]][['col']] = c('#ff0010', '#1515ff', 'white', 'yellow', 'orange', 'purple')[type+1]
        } else {
          snapshot[[field]][['col']] = 'white'
        }
      }
    }
    out[[field]][['col']] = snapshot[[field]][['col']]

    if (length(snapshot[[field]]$col)==1) {
      col = grDevices::col2rgb(snapshot[[field]]$col)/255
      snapshot[[field]]$colrgb = cbind(col%*%seq(0,1,length=5000),c(0,0,0))
    } else {
      snapshot[[field]]$colrgb = grDevices::col2rgb(snapshot[[field]]$col)/255
      if (hdcolors) {
        # smooth color scale to overcome 8-bit representation
        ncol = dim(snapshot[[field]]$colrgb)[2]
        if (ncol>3) {
          for (d in seq(3)) {
            snapshot[[field]]$colrgb[d,] = cooltools::smoothfun(seq(ncol),c(snapshot[[field]]$colrgb[d,]),df=max(2,min(20,floor(sqrt(ncol)))))(seq(ncol))
          }
        }
      }
      snapshot[[field]]$colrgb = cbind(cooltools::lim(snapshot[[field]]$colrgb),c(0,0,0))
    }

  }
  # end input handling #########################################################

  if (is.null(center) | is.null(width)) {

    # extract all particles positions
    x = allpart(snapshot,'Coordinates',species = all.types)

    # determine if the data is a cubic box
    L = apply(apply(x,2,range),2,diff)
    dx = apply(apply(x,2,range),1,diff)
    Lmax = max(L)
    is.box = all(L/Lmax>1-1e-3) & all(abs(dx/Lmax)<1e-3)
  }

  # determine geometric center
  if (is.null(center)) {
    if (is.box) {
      center = rep((max(x)+min(x))/2,3)
    } else {
      center = colSums(x)/dim(x)[1] # geometric centre
    }
  } else {
    if (length(center)==1) center=rep(center,3)
    if (length(center)==2) center=c(center,0)
  }

  # determine plotting limits
  if (is.null(width)) {
    if (is.box) {
      width = Lmax
    } else {
      width = 2*sqrt(2)*max(cooltools::vectornorm(t(t(x)-center)))/sqrt(1+1/aspect^2)
    }
  }
  height = width/aspect
  xlim = c(-1,1)*width/2
  ylim = c(-1,1)*height/2
  mean.length = sqrt(width*height)

  # determine arrows, scale, titles
  if (is.null(arrow.origin)) arrow.origin = c(xlim[1],ylim[1])+0.05*mean.length
  if (is.null(arrow.length)) arrow.length = 0.1*mean.length
  if (is.null(scale.origin)) scale.origin = c(xlim[2]-0.05*mean.length,ylim[1]+0.05*mean.length)
  if (is.null(scale.length)) scale.length = 0.1*mean.length
  scale.length = signif(scale.length,1)
  if (is.null(title.origin)) title.origin = c(xlim[1]+0.06*mean.length,ylim[2]-0.06*mean.length)
  dtext = 0.01*mean.length*(1+text.offset)

  # make rotation matrix
  if (length(rotation)==3) {
    rot = t(cooltools::rotation3(rotation))
    if (is.null(xlab)) xlab = expression(e[1])
    if (is.null(ylab)) ylab = expression(e[2])
  } else if (length(rotation)==1) {
    if (rotation>3) {
      if (is.null(x)) x = allpart(snapshot)
      e = eigen(cooltools::quadrupole(t(t(x)-center)))$vectors
    }
    if (rotation==1) {
      rot = diag(3)
      if (is.null(xlab)) xlab = 'x'
      if (is.null(ylab)) ylab = 'y'
    } else if (rotation==2) {
      rot = rbind(c(0,1,0),c(0,0,1),c(1,0,0))
      if (is.null(xlab)) xlab = 'y'
      if (is.null(ylab)) ylab = 'z'
    } else if (rotation==3) {
      rot = rbind(c(1,0,0),c(0,0,1),c(0,1,0))
      if (is.null(xlab)) xlab = 'x'
      if (is.null(ylab)) ylab = 'z'
    } else if (rotation==4) {
      rot = t(e[,c(1,3,2)])
      if (is.null(xlab)) xlab = expression(lambda ['max'])
      if (is.null(ylab)) ylab = expression('  '*lambda ['min'])
    } else if (rotation==5) {
      rot = t(e[,c(1,2,3)])
      if (is.null(xlab)) xlab = expression(lambda ['max'])
      if (is.null(ylab)) ylab = expression('  '*lambda ['mid'])
    } else if (rotation==6) {
      rot = t(e[,c(2,3,1)])
      if (is.null(xlab)) xlab = expression(lambda ['mid'])
      if (is.null(ylab)) ylab = expression('  '*lambda ['min'])
    } else {
      stop('rotation must be an integer 1,...,6, a real 3-vector, or a 3-by-3 matrix.')
    }
  } else if (length(rotation)==9) {
    rot = rotation
  } else {
    stop('rotation must be an integer 1,...,6, a real 3-vector, or a 3-by-3 matrix.')
  }

  rm(x) # to free up memory
  gc()

  # prepare grid
  nx = round(npixels*width/mean.length)
  ny = round(npixels*height/mean.length)
  if (max(nx,ny)>8000) stop('Not more than 8000 pixels allowed on each side.')
  dx = width/nx # pixel size in simulation length units

  # put particles on a grid
  for (type in types) {

    field = sprintf('PartType%d',type)

    # get positions
    x = snapshot[[field]]$Coordinates
    if (!is.array(x)) x=rbind(x)

    # subsampling
    if (!is.null(sample.fraction)) {
      if (sample.fraction<1) {
        nsub = max(1,round(dim(x)[1]*sample.fraction)) # number of particles to select
        sel = sample(dim(x)[1],nsub)
        x = x[sel,,drop = FALSE]
        if (snapshot[[field]]$color.by.property) snapshot[[field]]$value=snapshot[[field]]$value[sel]
        if (!is.null(snapshot[[field]]$Density)) snapshot[[field]]$Density=snapshot[[field]]$Density[sel]
        if (!is.null(snapshot[[field]]$Masses)) snapshot[[field]]$Masses=snapshot[[field]]$Masses[sel]
      }
    }

    # compute total number of particles
    out[[field]]$n.tot = dim(x)[1]

    # translate particles to custom center and rotate
    if (is.null(rot.center)) {
      x = t(rot%*%(t(x)-center)) # probably the fastest way to do this
    } else {
      x = t(rot%*%(t(x)-rot.center)+rot.center-center)
    }

    # stereographic projection
    if (!is.null(fov)) {
      sel = which(x[,3]>0)
      if (snapshot[[field]]$color.by.property) snapshot[[field]]$value=snapshot[[field]]$value[sel]
      if (!is.null(snapshot[[field]]$Density)) snapshot[[field]]$Density=snapshot[[field]]$Density[sel]
      if (!is.null(snapshot[[field]]$Masses)) snapshot[[field]]$Masses=snapshot[[field]]$Masses[sel]
      stretch = width/2/tan(fov/180*pi/2)
      x = x[sel,, drop = FALSE]
      distance = vectornorm(x)
      x = x[,1:2]/x[,3]*stretch
    }

    # sub-select slice thickness
    weight = NULL
    kernel = function(x,h) {
      q = x/h
      ifelse(q<1,2/(3*h)*(1-3*q^2/2*(1-abs(q)/2)),pmax(0,1/(6*h)*(2-q)^3))
    }
    if (!is.null(depth)) {
      if (is.null(fov)) {
        if (taper) {
          sel = which(abs(x[,3])<=depth)
          weight = kernel(abs(x[,3]),depth/2)
        } else {
          sel = which(abs(x[,3])<=depth/2)
        }
      } else {
        if (taper) {
          sel = which(distance<=2*depth)
          distance = distance[sel]
          weight = kernel(distance,depth)
        } else {
          sel = which(distance<=depth)
        }
      }
      x = x[sel,, drop = FALSE]
      if (snapshot[[field]]$color.by.property) snapshot[[field]]$value=snapshot[[field]]$value[sel]
      if (!is.null(snapshot[[field]]$Density)) snapshot[[field]]$Density=snapshot[[field]]$Density[sel]
      if (!is.null(snapshot[[field]]$Masses)) snapshot[[field]]$Masses=snapshot[[field]]$Masses[sel]
    }

    # weights
    if (!is.null(snapshot[[field]]$Masses)) {
      if (is.null(weight)) {
        weight = snapshot[[field]]$Masses
      } else {
        weight = weight*snapshot[[field]]$Masses
      }
    }

    #  raster particle field
    sigma = NULL
    if (snapshot[[field]]$kde==4) {
      if (is.null(snapshot[[field]]$Density) | is.null(snapshot[[field]]$Masses)) {
        snapshot[[field]]$kde = 3
      } else {
        sigma = (snapshot[[field]]$Masses/snapshot[[field]]$Density)^(1/3)
      }
    }
    out[[field]]$density = cooltools::kde2(x, w=weight, xlim=xlim, ylim=ylim, n=nx, smoothing=snapshot[[field]]$smoothing,
                                sigma=sigma, algorithm=snapshot[[field]]$kde)$field
    out[[field]]$density[out[[field]]$density<0] = 0
    out[[field]]$n.eff = sum(out[[field]]$density)*dx^2

    # raster particle values (optional)
    if (snapshot[[field]]$color.by.property) {
      if (is.null(weight)) {
        w = as.vector(snapshot[[field]]$value)
      } else {
        w = as.vector(snapshot[[field]]$value)*weight
      }
      out[[field]]$value = cooltools::kde2(x, w=w, sigma=sigma, xlim=xlim, ylim=ylim, n=nx, smoothing=snapshot[[field]]$smoothing, algorithm=snapshot[[field]]$kde)$field/out[[field]]$density
      out[[field]]$value[!is.finite(out[[field]]$value)] = 0
      out[[field]]$value[out[[field]]$value<0] = 0
    }

  }

  # deterime mean density across all types
  mass = 0
  for (type in types) {
    field = sprintf('PartType%d',type)
    mass = mass+out[[field]]$n.eff
  }
  mean.density = mass/dx^2/nx/ny/length(types)

  # turn density and value matrices into RGB layers
  nlayers = length(out)
  img4 = array(dim=c(nx,ny,3,nlayers))
  layer = 0

  for (type in types) {

    field = sprintf('PartType%d',type)
    layer = layer+1

    # convert density to brightness
    filling.factor = sum(out[[field]]$density>mean(out[[field]]$density)*0.1)/npixels^2
    if (is.null(snapshot[[field]]$ref.density)) {
      linear.scaling = filling.factor/mean.density
    } else {
      linear.scaling = 1/snapshot[[field]]$ref.density
    }
    brightness = 0.1*(snapshot[[field]]$lum*linear.scaling*out[[field]]$density)^snapshot[[field]]$gamma

    # if no values provided use density as values
    if (snapshot[[field]]$color.by.property) {
      val = out[[field]]$value
    } else {
      val = brightness
    }

    # turn value into color
    nvalcol = dim(snapshot[[field]]$colrgb)[2]-1
    for (k in seq(3)) {
      normalised.value = (as.vector(val)-snapshot[[field]]$valrange[1])/diff(snapshot[[field]]$valrange)
      index = round(pmax(0,pmin(1,normalised.value))*(nvalcol-1)+1)
      index[is.na(index)] = nvalcol+1
      img4[,,k,layer] = snapshot[[field]]$colrgb[k,index]*(1+pmax(0,normalised.value-1))
    }

    # adjust brightness as a function of density if hue represents a property
    if (snapshot[[field]]$color.by.property) {
      density.scaling = TRUE
      if (!is.null(snapshot[[field]]$density.scaling)) density.scaling=snapshot[[field]]$density.scaling
      if (density.scaling) {
        for (k in seq(3)) {
          img4[,,k,layer] = img4[,,k,layer]*brightness
        }
      }
    }
  }

  # combine layers
  if (layer!=nlayers) stop('wrong number of layers')

  # sum color layers
  img = array(0,dim=c(nx,ny,3))
  for (layer in seq(nlayers)) img=img+img4[,,,layer]

  # apply non-linear brightness scale to each RGB channel individually
  img = atan(img)/pi*2
  f = 10^max(0,2.2*shadows)
  img = log10(f*img+1)/log10(f+1)
  img = cooltools::lim(img) # to get rid of small numerical outliers

  # save raster image as png
  if (!is.null(pngfile)) {
    png::writePNG(cooltools::rasterflip(img),pngfile)
  }

  # show on screen and save as pdf
  for (mode in seq(2)) {

    make = FALSE

    if (mode==1 & screen) {
      make = TRUE
    }

    if (mode==2 & !is.null(pdffile)) {
      make = TRUE
      par(mar=rep(0,4))
      grDevices::pdf(pdffile,width=7*width/mean.length,height=7*height/mean.length)
      par(mar=rep(0,4))
    }

    if (make) {

      # initialize plot
      cooltools::nplot(xlim=xlim, ylim=ylim, asp=1, ...)

      # plot raster
      graphics::rasterImage(cooltools::rasterflip(img),xlim[1],ylim[1],xlim[2],ylim[2])

      # arrows
      if (arrows) {
        graphics::arrows(arrow.origin[1],arrow.origin[2],arrow.origin[1]+arrow.length,arrow.origin[2],col=text.col,length = 0.1,angle=20,lwd=arrow.lwd)
        graphics::arrows(arrow.origin[1],arrow.origin[2],arrow.origin[1],arrow.origin[2]+arrow.length,col=text.col,length = 0.1,angle=20,lwd=arrow.lwd)
        if (!is.null(xlab)) {
          graphics::text(arrow.origin[1]+arrow.length+dtext,arrow.origin[2],xlab,pos=4,col=text.col,cex=cex)
        }
        if (!is.null(ylab)) {
          graphics::text(arrow.origin[1],arrow.origin[2]+arrow.length+dtext,ylab,pos=3,col=text.col,cex=cex)
        }
      }

      # length scale
      if (scale) {
        graphics::text(scale.origin[1]-scale.length/2,scale.origin[2]+dtext,sprintf('%s %s',signif(scale.length,1),length.unit),pos=3,col=text.col,cex=cex)
        graphics::lines(scale.origin[1]-c(0,scale.length),rep(scale.origin[2],2),col=text.col,lwd=scale.lwd)
      }

      # title
      if (!is.null(title)) graphics::text(title.origin[1], title.origin[2], title, pos=4, col=text.col, offset=-0.4)

      if (mode==2) {grDevices::dev.off()}

    }

  }

  # return results
  out$rgb = img
  out$header = list(xlim=xlim, ylim=ylim, center=center, rotationmatrix=rot)
  invisible(out)

}
