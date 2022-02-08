#' Visualise one or several 3D point sets
#'
#' @importFrom cooltools lim griddata2 kde2 quadrupole rotation3 vectornorm nplot rasterflip lightness
#' @importFrom png writePNG
#' @importFrom grDevices pdf dev.off col2rgb rainbow
#' @importFrom graphics axis lines par rasterImage rect text arrows
#'
#' @description Produces a raster image visualizing particle positions in 3D N-body/SPH data.
#'
#' @param x Gadget data type, which must contain at least one sublist PartType# (with #=0,1,...). Each sublist PartType# represents one type of particle (e.g. gas, stars, dark matter) and must contain at least the particles coordinates in an N-by-3 matrix \code{Coordinates}. Other optional elements of PartType# are:\cr
#' \code{col} = color used to display this particle type. This can be a single color or a vector of colors. If a single color is provided, a range of brightness is produced automatically following the \code{color.model}. If a vector is provided, it is interpreted as a color scale applied to the particles.\cr
#' \code{color.model} = character string that must be either \code{hsv} or \code{hsl}. In the former case, the color \code{col} is expanded into a color scale ranging from black to that particular color. In the latter case, the color scale ranges from black to white passing through the color \code{col}.\cr
#' \code{value} = optional N-vector of values, one for each particle, which define the position on the color scale. If not given, particles are colored according to their projected density.\cr
#' \code{valrange} = optional 2-vector specifying the values corresponding to the limits of the color scale. Only used if the N-vector \code{value} is given.\cr
#' \code{density.scaling} = logical flag (default TRUE), specifying whether the brightness is scaled with the particle density if \code{value} is given.\cr
#' \code{lum} = overall scaling factor for the color scale.\cr
#' \code{hdr} = parameter specifying the size of the dynamic range to be represented by the color scale. The larger the value, the higher the dynamic range.\cr
#' \code{kde} = logical flag (default TRUE), whether the particles are smoothed using an adaptive kernel density estimator.\cr
#' \code{smoothing} = smoothing scale in simulation units. If not given, a value is set automatically to one percent of the geometric mean of width and height.\cr\cr
#' @param types vector specifying the particle types to be displayed. These numbers must correspond to the # in the sublists PartType#. If not given, all particle species are shown.
#' @param combine.fun function used to combined the image layers corresponding to different particle species. Typically, this is set to \code{sum} or \code{mean}.
#' @param center optional 3-vector specifying the coordinate at the center of the plot. The default is the geometric center (= center of mass, if all particle masses are equal).
#' @param rotation either an integer (1-6), a 3-vector or a 3-by-3 matrix, specifying a rotation of the 3D particle positions. In case of an integer: 1=(x,y)-plane, 2=(y,z)-plane, 3=(x,z)-plane, 4=(qmax,qmin)-plane, 5=(qmax,qmid)-plane, 6=(qmid,qmin)-plane, where qmax/qmid/qmin are the eigenvectors of the particle-quadrupole, associated with the maximum/middle/minimum eigenvalues, respectively. If \code{rotation} is a vector, its direction specifies the rotation axis and its norm the rotation angle in radians in the positive geometric sense.
#' @param width optional horizontal range of the image in simulation units. The default corresponds to the full range of particle positions.
#' @param thickness optional thickness of the slice to be plotted in simulation units. If not given, all particles are shown.
#' @param aspect aspect ratio (= width/height).
#' @param ngrid number of grid cells per side in the output image. If the image is not square ngrid is interpreted as the geometric mean between the horizontal and the vertical number of pixels, such that the total number of pixels remains about ngrid^2.
#' @param fix.luminosity logical flag specifying whether the brightness should scale with the number of particles in the field of view. Set this to \code{FALSE} to avoid luminosity fluctuations in movies.
#' @param color.model default value used in \code{dat$PartType#} (see above).
#' @param lum default value used in \code{dat$PartType#} (see above).
#' @param hdr default value used in \code{dat$PartType#} (see above).
#' @param kde default value used in \code{dat$PartType#} (see above).
#' @param smoothing default value used in \code{dat$PartType#} (see above).
#' @param sample.fraction fraction of particles to be used. Ff 1, all particles are used, if <1 a random subsample is drawn.
#' @param screen logical flag specifying whether the images is displayed on the screen.
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
#' @author Danail Obreschkow
#'
#' @examples
#' ## Example of 2x1e4 particles
#'
#' # Produce a density plot with different colors for the two particle species
#' dat = as.snapshot(list(cooltools::runif3(1e4), array(rnorm(3e4),c(1e4,3))))
#' oldpar = par(mar=c(1,1,1,1))
#' plot(dat, width=5)
#' par(oldpar)
#'
#' # Only plot particles of species "1" and vary hue with distance from centre
#' xyradius = as.vector(cooltools::vectornorm(dat$PartType1$Coordinates[,1:2]))
#' dat$PartType1$value = xyradius
#' dat$PartType1$valrange = c(0,2)
#' oldpar = par(mar=c(1,1,1,6))
#' out = plot(dat, width=5, types=1)
#' cooltools::colorbar(2.7,-2.5,3,2.5,
#'                     col=out$PartType1$col,clim=out$PartType1$valrange,
#'                     text='xy-distance from center')
#' par(oldpar)
#'
#' @method plot snapshot
#' @export
plot.snapshot = function(x, center=NULL, rotation=1, thickness=NULL, width=NULL, aspect=1, ngrid=300, kde=TRUE, smoothing=NULL, types=NULL,
                       sample.fraction=1, lum=1, color.model='hsv', hdr=1, fix.luminosity=FALSE, screen=TRUE, pngfile=NULL, pdffile=NULL,
                       combine.fun = sum, title=NULL, title.origin = NULL,
                       arrows = TRUE, arrow.origin = NULL, arrow.length = NULL, arrow.lwd = 1.5,
                       scale = TRUE, scale.origin = NULL, scale.length = NULL, scale.lwd = 1.5, length.unit = '',
                       xlab = NULL, ylab = NULL, cex=1, text.offset = 0, text.col = 'white', ...) {

  dat = x; x = NULL

  out = list()

  # determine particle types to be considered
  if (is.null(types)) {
    types = c()
    for (type in seq(0,20)) {
      field = sprintf('PartType%d',type)
      if (!is.null(dat[[field]])) types=c(types,type)
    }
  }

  # check and pre-process input arguments ######################################
  if (min(types)<0) stop('particle types must be non-negative integers')

  for (type in types) {

    field = sprintf('PartType%d',type)
    out[[field]] = list()

    # handle brightness parameters
    if (is.null(dat[[field]]$lum)) dat[[field]]$lum=lum
    if (is.null(dat[[field]]$hdr)) dat[[field]]$hdr=hdr

    # handle values
    dat[[field]]$color.by.property = !is.null(dat[[field]]$value)
    if (dat[[field]]$color.by.property) {
      if (is.null(dat[[field]]$valrange)) dat[[field]]$valrange=range(dat[[field]]$value)
    } else {
      dat[[field]]$valrange = c(0,1)
    }
    out[[field]]$valrange = dat[[field]]$valrange

    # color model
    if (is.null(dat[[field]]$color.model)) {
      dat[[field]]$color.model = color.model
    }
    if (!dat[[field]]$color.model%in%c('hsl','hsv')) {
      stop('if specified, the color.model must be equal to "hsl" or "hsv".')
    }

    # determine particle color
    if (is.null(dat[[field]]$col)) {
      if (dat[[field]]$color.by.property) {
        dat[[field]]$col = grDevices::rainbow(256,end=5/6)
      } else {
        if (type<=5) {
          dat[[field]]$col = c('#ff0010', '#0515ff', 'green', 'orange', 'yellow', 'purple')[type+1]
        } else {
          dat[[field]]$col = 'white'
        }
      }
    }
    if (length(dat[[field]]$col)==1) {
      if (dat[[field]]$color.model=='hsv') {
        dat[[field]]$col = cooltools::lightness(dat[[field]]$col,seq(0,0.5,length=1000))
      } else {
        dat[[field]]$col = cooltools::lightness(dat[[field]]$col,seq(0,1,length=1000))
      }
    }
    out[[field]]$col = dat[[field]]$col
    dat[[field]]$colrgb = cbind(grDevices::col2rgb(dat[[field]]$col)/255,c(0,0,0))

  }
  # end input handling #########################################################

  x = allpart(dat,'Coordinates')

  # determine geometric center
  if (is.null(center)) {
    center = apply(x,2,mean) # geometric centre
  } else {
    if (length(center)==1) center=rep(center,3)
    if (length(center)==2) center=c(center,0)
  }

  # determine plotting limits
  if (is.null(width)) width = 2*sqrt(2)*max(cooltools::vectornorm(x))/sqrt(1+1/aspect^2)
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
    if (rotation>3) e = eigen(cooltools::quadrupole(x))$vectors
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
    rot = matrix(rotation,3,3)
  } else {
    stop('rotation must be an integer 1,...,6, a real 3-vector, or a 3-by-3 matrix.')
  }

  x = NULL # to free up memory

  # prepare grid
  nx = round(ngrid*width/mean.length)
  ny = round(ngrid*height/mean.length)
  if (max(nx,ny)>8000) stop('Not more than 8000 pixels allowed on each side.')
  dx = width/nx # pixel size in simulation length units

  # determine default smoothing
  if (is.null(smoothing)) smoothing = mean.length*0.01 # [length units of sim]

  # put particles on a grid
  for (type in types) {

    field = sprintf('PartType%d',type)

    # specify smoothing kernel
    if (is.null(dat[[field]]$smoothing)) dat[[field]]$smoothing = smoothing
    if (is.null(dat[[field]]$kde)) {
      kde = kde
    } else {
      kde = dat[[field]]$kde
    }
    out[[field]]$smoothing = dat[[field]]$smoothing
    out[[field]]$kde = dat[[field]]$kde

    # get positions
    x = dat[[field]]$Coordinates

    # compute total number of particles
    out[[field]]$n.tot = dim(x)[1]

    # translate particles to custom center
    x = sweep(x, 2, center)

    # rotate coordinates (active rotation)
    x = t(rot%*%t(x))

    # default selection
    sel = seq(dim(x)[1])

    # sub-select slice thickness
    if (!is.null(thickness)) {
      sel = sel[abs(x[,3])<=thickness/2]
    }

    # subsampling
    if (!is.null(sample.fraction)) {
      if (sample.fraction<1) {
        nsub = max(1,round(length(sel)*sample.fraction)) # number of particles to select
        sel = sample(sel,nsub)
      }
    }

    #  raster particle data
    if (dat[[field]]$smoothing==0) {
      g = cooltools::griddata2(x[sel,1], x[sel,2], w=as.vector(dat[[field]]$value[sel]), xlim=xlim, ylim=ylim, n=c(nx,ny))
      out[[field]]$density = g$n
      if (dat[[field]]$color.by.property) out[[field]]$value = g$m/out[[field]]$density
    } else {
      if (kde) {
        g = cooltools::kde2(x[sel,1], x[sel,2], xlim=xlim, ylim=ylim, n=c(nx,ny), s=dat[[field]]$smoothing/8/dx, sd.max=dat[[field]]$smoothing*2/dx, cpp=TRUE)
        out[[field]]$density = g$d
        if (dat[[field]]$color.by.property) {
          g = cooltools::kde2(x[sel,1], x[sel,2], w=as.vector(dat[[field]]$value[sel]), xlim=xlim, ylim=ylim, n=c(nx,ny), s=dat[[field]]$smoothing/8/dx, sd.max=dat[[field]]$smoothing*2/dx, cpp=TRUE)
          out[[field]]$value = g$d/out[[field]]$density
        }
      } else {
        if (!requireNamespace("EBImage", quietly=TRUE)) {
          stop('Package EBImage is needed in function plot.gadget if kde=FALSE. Consider using kde=TRUE.')
        }
        g = cooltools::griddata2(x[sel,1], x[sel,2], w=as.vector(dat[[field]]$value[sel]), xlim=xlim, ylim=ylim, n=c(nx,ny))
        out[[field]]$density = EBImage::gblur(g$n, dat[[field]]$smoothing/dx)
        if (dat[[field]]$color.by.property) out[[field]]$value = EBImage::gblur(g$m, dat[[field]]$smoothing/dx)/out[[field]]$density
      }
    }
    out[[field]]$n.eff = sum(out[[field]]$density)
  }

  # turn density and value matrices into RGB layers
  nlayers = sum(dat$Header$NumPart_ThisFile[types+1]>0)
  img = array(0,c(nx,ny,3,nlayers))
  layer = 0

  for (type in types) {
    if (dat$Header$NumPart_ThisFile[type+1]>0) {

      field = sprintf('PartType%d',type)
      layer = layer+1

      # convert density to brightness
      # 1) linearly rescale the density field as a function of grid cells, particle numbers and custom 'lum' parameter
      scaling = (0.5*dat[[field]]$lum)*ngrid^2/max(1,ifelse(fix.luminosity,out[[field]]$n.tot,out[[field]]$n.eff)) # linear luminosity scaling factor
      brightness = scaling*out[[field]]$density
      # 2) remap density field from the domain [0,infinity) to [0,1] using a generalised Sigmoid function
      sigmoid = function(x,a=1) {
        f = 1 # value of x where sigmoid(x)=x/2
        b = ((2/(2-f))^a-1)/(a*f)
        1-1/(x*a*b+1)^(1/a)
      }
      brightness = sigmoid(brightness,dat[[field]]$hdr)
      # 3) ensure that the density field is strictly contained in [0,1], by cropping minute outlines caused by floating-point errors
      brightness = cooltools::lim(brightness)

      # if no values provided use density as values
      if (dat[[field]]$color.by.property) {
        val = out[[field]]$value
      } else {
        val = brightness
      }

      # turn value into color
      nvalcol = dim(dat[[field]]$colrgb)[2]-1
      for (k in seq(3)) {
        index = round(pmax(0,pmin(1,(as.vector(val)-dat[[field]]$valrange[1])/diff(dat[[field]]$valrange)))*(nvalcol-1)+1)
        index[is.na(index)] = nvalcol+1
        img[,,k,layer] = dat[[field]]$colrgb[k,index]
      }

      # adjust brightness as a function of density if hue represents a property
      if (dat[[field]]$color.by.property) {
        density.scaling = TRUE
        if (!is.null(dat[[field]]$density.scaling)) density.scaling=dat[[field]]$density.scaling
        if (density.scaling) {
          for (k in seq(3)) {
            img[,,k,layer] = img[,,k,layer]*brightness
          }
        }
      }

    }
  }

  # combine layers
  if (layer!=nlayers) stop('wrong number of layers')
  img = cooltools::lim(apply(img,1:3,combine.fun))

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
      grDevices::pdf(pdffile,width=7*width/mean.length,height=7*height/mean.length)
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
  out$header = list(xlim=xlim, ylim=ylim, center=center, rotationmatrix = rot)
  invisible(out)

}
