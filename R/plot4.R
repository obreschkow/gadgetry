#' Display a simulation snapshot in four different projections
#'
#' @importFrom cooltools nplot
#' @importFrom grDevices pdf dev.off graphics.off
#' @importFrom graphics lines par text
#'
#' @description Produces a 2-by-2 panel figure, where each panel is made of an \code{\link{plot.snapshot}} image in a different rotation.
#'
#' @param x Gadget data type, which must contain at least one sublist PartType# (with #=0,1,...). For details, see \code{\link{plot.snapshot}}.
#' @param rotations 4-vector or 4-element list specifying four different values/vectors of the argument \code{rotation} in \code{\link{plot.snapshot}}. These different rotations are shown in the bottom left, bottom right, top left, and top right panel, respectively, in this order.
#' @param screen logical flag specifying whether the images is displayed on the screen.
#' @param pdffile optional pdf-filename to save the image as pdf-file.
#' @param title Text to be added to the figure.
#' @param scale logical flag, specifying if a length scale is shown
#' @param width optional horizontal range of the image in the length units of \code{x}. The default corresponds to the full range of particle positions.
#' @param ... additional parameters for \code{\link{plot.snapshot}}.
#'
#' @return None
#'
#' @seealso \code{\link{plot.snapshot}}
#'
#' @export

plot4 = function(x, rotations=c(2,4,1,3), screen = TRUE, pdffile = NULL, title = NULL, scale = TRUE, width = NULL, ...) {

  oldpar = par(no.readonly = TRUE)
  on.exit(par(oldpar))

  for (mode in seq(2)) {

    make = FALSE
    if (mode==1 & screen) {
      make = TRUE
    }
    if (mode==2 & !is.null(pdffile)) {
      make = TRUE
      grDevices::pdf(pdffile,width=7,height=7)
    }

    if (make) {

      if (mode==1) grDevices::graphics.off()
      par(pty='s', mar=c(0,0,0,0))
      cooltools::nplot(xlim=c(0,1), ylim=c(0,1))
      p = par()$plt

      ix = iy = 0

      for (i in seq(4)) {

        par(fig=c(p[1]+(p[2]-p[1])*ix/2,p[1]+(p[2]-p[1])*(ix/2+0.5),p[3]+(p[4]-p[3])*iy/2,p[3]+(p[4]-p[3])*(iy/2+0.5)),
            new=TRUE, mar=c(0,0,0,0) )

        if (i==1) {
          xlim = plot.snapshot(x, rotation = rotations[[i]], pngfile = NULL, pdffile = NULL, screen = TRUE,
                      scale=FALSE, width=width, ...)$xlim
          if (is.null(width)) {width=diff(xlim)}
        } else {
          plot.snapshot(x, rotation = rotations[[i]], pngfile = NULL, pdffile = NULL, screen = TRUE,
                      scale=ifelse(i==2,scale,FALSE), width=width, ...)
        }
        ix = (ix+1)%%2
        if (ix==0) iy = iy+1

      }

      # lines between panels
      par(fig=p, new=TRUE, mar=c(0,0,0,0))
      cooltools::nplot(xlim=c(0,1), ylim=c(0,1))
      graphics::lines(c(1,1)/2,c(0,1),col='grey')
      graphics::lines(c(0,1),c(1,1)/2,col='grey')

      # title
      if (!is.null(title)) graphics::text(0.06, 1.94, title, pos=4, col='white', offset=-0.4)

      if (mode==2) grDevices::dev.off()

    }

  }

}
