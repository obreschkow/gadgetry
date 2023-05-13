#' Standard spherical mass distributions in astrophysics
#'
#' @importFrom graphics hist curve abline
#' @importFrom stats runif
#'
#' @description Generates functions associated with spherically symmetric halo models. Currently only Hernquist model implemented. Extension planned. Everything is in units dimensionless units where M=R=G=1, where M is the total halo mass, R is the characteristic radius (or truncation radius if no characteristic radius exists) and G is the gravitational constant.
#'
#' @param model name of model. Currently only "Hernquist".
#' @param rmax truncation radius; profiles are normalized, such that this radius contains the full mass
#' @param plots logical flag indicating whether to produce a series of test plots
#' @param rplot maximum plotting radius; only used if plots=TRUE.
#'
#' @return list of functions for the selected model
#'
#' @author Danail Obreschkow
#'
#' @examples
#'
#' f = halomodel(plots = TRUE, rplot = 5, rmax=3)
#'
#' @export

halomodel = function(model='hernquist', rmax=NULL, plots=FALSE, rplot=10) {

  f = list(model = tolower(model), notes=c())

  if (f$model=='hernquist') {

    # truncation radius
    if (is.null(rmax)) {
      f$rmax = Inf # default truncation radius
    } else {
      f$rmax = rmax
    }

    # 3D Density dM/dV
    # Must be normalised such that integral(4*pi*r^2*rho*dr,0,rmax)=1
    normalization = ifelse(is.infinite(f$rmax),1,1/(1-(2*f$rmax+1)/(f$rmax+1)^2))
    f$rho = function(r) ifelse(r<=f$rmax,normalization/(2*pi*r*(1+r)^3),0)

    # 1D Probability density dP/dr=4*pi*r^2*rho (since M=1)
    # is normalised such that integral(density)
    f$pdf = function(r) ifelse(r<=f$rmax,normalization*2*r/(1+r)^3,0)

    # Enclosed mass M(<r)=integral(pdf,0,r)
    f$mass = function(r) ifelse(r<=f$rmax,normalization*(1-(2*r+1)/(r+1)^2),1)

    # Isotropic velocity dispersion (for r<rmax)
    variance = function(r) r*(1+r)^3*log((1+r)/r)-(r*(25+52*r+42*r^2+12*r^3))/(12*(1+r)) # for rmax=Inf, runs into floating-point issues at r>1e2, seriously r>1e3
    f$sigma = function(r) sqrt(normalization*ifelse(r<=200,variance(r),exp(-1.6152-log(r))))
    f$notes = c(f$notes,'isotropic velocity dispersion sigma only valid for r<rmax')

    # Potential phi(r)=integral(pdf(r)*mass(<r)/r,0,r)+const, (since G=1)
    # with const defined such that phi vanishes at r=Inf
    const = ifelse(is.infinite(f$rmax),-1,-1/f$rmax-f$rmax/(1+f$rmax))
    f$potential = function(r) ifelse(r<=f$rmax,r/(1+r)+const,-1/r)

    # Quantile function for 1D PDF
    f$quantile = function(p) (p+sqrt(normalization*p))/(normalization-p)

  } else {

    stop('unknown model')

  }

  # Derive a few more quantities using general equations for any spherical mass distr.:

  # Random radius generator
  f$random = function(n) f$quantile(stats::runif(n))

  # Circular velocity v=sqrt(M(<r)/r)
  f$vcirc = function(r) sqrt(f$mass(r)/r)

  # Local escape velocity [equation true for any spherical mass distribution]
  f$vesc = function(r) sqrt(-2*f$potential(r))

  if (plots) {

    mycurve = function(x,...) {
      graphics::curve(x,xlab='Radius r',n=1e3,...)
      if (rplot>f$rmax) graphics::abline(v=f$rmax,lty=3)
    }

    fct = function(x) ifelse(x<=f$rmax,f$rho(x),NA)
    mycurve(fct,0.1,rplot,ylab='f$rho',log='xy')
    mycurve(f$mass,0,rplot,ylab='f$mass')
    mycurve(f$potential,0,rplot,ylab='f$potential')
    mycurve(f$vcirc,0,rplot,ylab='f$vcirc (f$sigma)')
    mycurve(f$sigma,0,rplot,add=TRUE,lty=2)
    n = 1e6
    r = f$random(n)
    r = r[r<rplot]
    graphics::hist(r,breaks=seq(0,rplot,length=101),xlim=c(0,rplot),main='Histogram of 1e6 random radii')
    fct2 = function(x) f$pdf(x)*n*rplot/100
    mycurve(fct2,add=TRUE,col='red',lwd=1.5)

    invisible(f)

  } else {

    return(f)

  }

}
