##' Generate standard spherical mass distributions in astrophysics
#'
#' @importFrom graphics hist curve abline
#' @importFrom stats runif integrate
#'
#' @description Generates potential-density pairs and auxiliary functions associated with spherically symmetric mass distributions commonly used in astrophysics. The units of all outputs are the same as those of \code{rmax}, \code{rscale}, \code{mtot} and \code{gravity}.
#'
#' @param model name of model. Currently available models: "uniform" (spherical top-hat), "isothermal" (singular isothermal sphere with logarithmic potential), "hernquist" (Hernquist 1990), .
#' @param rmax optional truncation radius. If not given, a natural model-dependent value is assigned, which is infinite for models where the mass converges in this case.
#' @param rscale scale radius for models, which admit a characteristic scale in addition to a possible truncation radius. For all other models, this parameter is ignored. \code{rmax} and \code{rscale} should be provided in the same units.
#' @param mtot total mass enclosed by the truncation radius
#' @param gravity gravitational constant
#' @param plots logical flag indicating whether to produce a series of test plots
#' @param rplot maximum plotting radius in the same units as \code{rmax} and \code{rscale}; only used if plots=TRUE.
#'
#' @return list of functions for the selected model
#'
#' @author Danail Obreschkow
#'
#' @examples
#'
#' halomodel('hernquist', plots = TRUE, rplot = 5, rmax=3)
#'
#' @export

halomodel = function(model='hernquist', rmax=NULL, rscale=1, mtot=1, gravity=1, plots=FALSE, rplot=2) {

  # Step 0
  # Some input checks and initialisations

  # some basic input checks
  if (rscale<=0) stop('rscale must be positive')
  if (mtot<=0) stop('mtot must be positive')
  if (gravity<=0) stop('gravity must be positive')
  if (rplot<=0) stop('rplot must be positive')

  # make list of units for later use in Step 3
  unit = list(mass = mtot,
              gravity = gravity)

  # Step 1
  # The following code generates potential-density pairs and auxiliary functions in a
  # non-dimensional way. That is, all function arguments and values are in units where the
  # gravitational constant, the total halo mass (to the truncation radius), and a characteristic radius
  # are all equal to one. What this characteristic radius is (e.g. rmax or rscale), can depend on the
  # model and is set via the variable unit$length.

  f = list(model = tolower(model))

  if (f$model=='hernquist') {

    # define the length unit (unit$length) for the non-dimensional formulation of this model,
    # as well as the truncation radius (f$rmax) and, if applicable, scale radius (f$rscale),
    # both in units of unit$length
    unit$length = rscale
    f$rscale = 1
    f$rmax = ifelse(is.null(rmax), Inf, rmax/unit$length)
    if (f$rmax<=0) stop('rmax must be larger than 0')

    # 3D Density dM/dV
    # Must be normalised such that integral(4*pi*r^2*rho*dr,0,rmax)=1
    normalization = ifelse(is.infinite(f$rmax),1,1/(1-(2*f$rmax+1)/(f$rmax+1)^2))
    f$rho = function(r) ifelse(r<=f$rmax,normalization/(2*pi*r*(1+r)^3),0)

    # 1D Probability density dM/dr=dP/dr=4*pi*r^2*rho (since M=1)
    f$pdf = function(r) ifelse(r<=f$rmax,normalization*2*r/(1+r)^3,0)

    # Enclosed mass M(<r)=integral(pdf,0,r)
    f$mass = function(r) ifelse(r<=f$rmax & !is.infinite(r),normalization*(1-(2*r+1)/(r+1)^2),1)

    # Circular velocity v=sqrt(M(<r)/r)
    f$vcirc = function(r) ifelse(r<=f$rmax,sqrt(normalization*(1/r-(2+1/r)/(r+1)^2)),sqrt(1/r))

    # Potential phi(r)=integral(pdf/r,0,r)+const
    const = ifelse(is.infinite(f$rmax),-1,-(normalization*f$rmax^2+f$rmax+1)/(f$rmax*(f$rmax+1)))
    f$potential = function(r) ifelse(r<=f$rmax & !is.infinite(r),normalization*r/(1+r)+const,-1/r)

    # Isotropic velocity dispersion from Jeans equation (for r<=f$rmax; NA otherwise)
    # d(rho*var)/dr = -dphi/dr; sigma = sqrt(var)
    variance = function(r) r*(1+r)^3*log((1+r)/r)-(r*(25+52*r+42*r^2+12*r^3))/(12*(1+r)) # for rmax=Inf, runs into floating-point issues at r>1e2, seriously r>1e3
    f$sigma = function(r) ifelse(r<=f$rmax,sqrt(normalization*ifelse(r<=200,variance(r),exp(-1.6152-log(r)))),NA)

    # Quantile function for 1D PDF
    f$quantile = function(p) (p+sqrt(normalization*p))/(normalization-p)

  } else if (f$model=='uniform') {

    # define the length unit (unit$length) for the non-dimensional formulation of this model,
    # as well as the truncation radius (f$rmax) and, if applicable, scale radius (f$rscale),
    # both in units of unit$length
    f$rmax = 1
    unit$length = ifelse(is.null(rmax), 1, rmax)
    if (unit$length<=0) stop('rmax must be larger than 0')
    if (is.infinite(unit$length)) stop('rmax must be finite')

    # 3D Density dM/dV
    f$rho = function(r) ifelse(r<=1,3/(4*pi),0)

    # 1D Probability density dM/dr=dP/dr=4*pi*r^2*rho (since M=1)
    f$pdf = function(r) ifelse(r<=1,3*r^2,0)

    # Enclosed mass M(<r)=integral(pdf,0,r)
    f$mass = function(r) ifelse(r<=1,r^3,1)

    # Circular velocity v=sqrt(M(<r)/r)
    f$vcirc = function(r) ifelse(r<=1,r,sqrt(1/r))

    # Potential phi(r)=integral(M(<r)/r^2,0,r)+const
    f$potential = function(r) ifelse(r<=1,r^2/2-3/2,-1/r)

    # Isotropic velocity dispersion from Jeans equation (for r<=f$rmax; NA otherwise)
    # d(rho*var)/dr = -dphi/dr; sigma = sqrt(var)
    f$sigma = function(r) ifelse(r<=1,sqrt((1-pmin(r^2,1))/2),NA) # pmin needed despite ifelse to avoid wired warnings around r=1 in plotting

    # Quantile function for 1D PDF
    f$quantile = function(p) p^(1/3)

  } else if (f$model=='isothermal') {

    # define the length unit (unit$length) for the non-dimensional formulation of this model,
    # as well as the truncation radius (f$rmax) and, if applicable, scale radius (f$rscale),
    # both in units of unit$length
    f$rmax = 1
    unit$length = ifelse(is.null(rmax), 1, rmax)
    if (unit$length<=0) stop('rmax must be larger than 0')
    if (is.infinite(unit$length)) stop('rmax must be finite')

    # 3D Density dM/dV
    f$rho = function(r) ifelse(r<=1,1/(4*pi*r^2),0)

    # 1D Probability density dM/dr=dP/dr=4*pi*r^2*rho (since M=1)
    f$pdf = function(r) ifelse(r<=1,1,0)

    # Enclosed mass M(<r)=integral(pdf,0,r)
    f$mass = function(r) ifelse(r<=1,r,1)

    # Circular velocity v=sqrt(M(<r)/r)
    f$vcirc = function(r) ifelse(r<=1,1,sqrt(1/r))

    # Potential phi(r)=integral(M(<r)/r^2,0,r)+const
    f$potential = function(r) ifelse(r<=1,log(r)-1,-1/r)

    # Isotropic velocity dispersion from Jeans equation (for r<=f$rmax; NA otherwise)
    # d(rho*f)/dr = -rho*dphi/dr; sigma=sqrt(f)
    f$sigma = function(r) ifelse(r<=1,1/sqrt(2),NA)

    # Quantile function for 1D PDF
    f$quantile = function(p) p

  } else {

    stop('unknown model')

  }

  # Step 2
  # Post-processing and checks of non-dimensional mass models

  # check if all required properties have been assigned
  desired.names = c('model','rmax','rho','pdf','mass','vcirc','sigma','potential','quantile')
  actual.names = names(f)
  if ('rscale'%in%actual.names) actual.names=actual.names[-which(actual.names=='rscale')]
  if (!all(desired.names%in%actual.names)) stop('not all desired properties exist')
  if (!all(actual.names%in%desired.names)) stop('more than the desired properties exist')

  # Derive a few more quantities in dimensionless units using general equations for any
  # spherical mass distribution; also run several checks

  # check mass distribution
  if (f$mass(0)!=0) stop('mass distribution error')

  # check quantile function
  cooltools::is.equal(f$quantile(0),0,stoptext='quantile q(0)')
  cooltools::is.equal(f$quantile(1),f$rmax,stoptext='quantile q(1)')

  # Random radius generator
  f$random = function(n) f$quantile(stats::runif(n))

  # Local escape velocity [equation true for any spherical mass distribution]
  f$vesc = function(r) sqrt(-2*f$potential(r))

  # Radial acceleration a=-M(<r)/r^2
  f$acc = function(r) ifelse(r==0,0,f$mass(r)/r^2)

  # Total potential energy, while checking consistency of potential-density pair
  f1 = function(r) -f$mass(r)*f$pdf(r)/r
  f2 = function(r) -(f$mass(r)/r)^2/2
  f3 = function(r) f$potential(r)*f$pdf(r)/2
  u1 = stats::integrate(f1,0,f$rmax)$value
  u2 = stats::integrate(f2,0,f$rmax)$value-1/2/f$rmax # the second term is the integral from rmax to infinity
  u3 = stats::integrate(f3,0,f$rmax)$value
  cooltools::is.equal(u1,u2,u3,stoptext='potential energy')
  f$epot = u3

  # Total kinetic energy for pure isotropic dispersion
  y = function(r) 3*f$sigma(r)^2/2*f$pdf(r) # factor three because of three dimensions
  f$ekin = stats::integrate(y,0,f$rmax)$value

  # Step 3
  # Turn non-dimensional mass models in to dimensional ones, by scaling them with the custom values of
  # mass, radius and gravitational constant

  # Complete set of units
  unit$velocity = sqrt(unit$gravity*unit$mass/unit$length)
  unit$acceleration = unit$gravity*unit$mass/unit$length^2
  unit$potential = unit$gravity*unit$mass/unit$length
  unit$density = unit$mass/unit$length^3
  unit$energy = unit$gravity*unit$mass^2/unit$length

  # Scale all quantities
  h = list(model = f$model,
           rho = function(r) unit$density*f$rho(r/unit$length),
           mass = function(r) unit$mass*f$mass(r/unit$length),
           potential = function(r) unit$potential*f$potential(r/unit$length),
           sigma = function(r) unit$velocity*f$sigma(r/unit$length),
           vcirc = function(r) unit$velocity*f$vcirc(r/unit$length),
           vesc = function(r) unit$velocity*f$vesc(r/unit$length),
           acc = function(r) unit$acceleration*f$acc(r/unit$length),
           epot = unit$energy*f$epot,
           ekin = unit$energy*f$ekin,
           pdf = function(r) f$pdf(r/unit$length)/unit$length,
           quantile = function(p) f$quantile(p)*unit$length,
           random = function(n) unit$length*f$random(n),
           rmax = f$rmax*unit$length)
  if (!is.null(f$rscale)) h$rscale = f$rscale*unit$length
  #print(names(f))
  #print(names(h))
  if (!all(names(h)%in%names(f))) stop('h in f problem')
  if (!all(names(f)%in%names(h))) stop('f in h problem')
  h$rscale = rscale
  h$mtot = mtot
  h$gravity = gravity

  # additional checks on dimensional quantities, just for extra safety
  cooltools::is.equal(integrate(function(r) 4*pi*r^2*h$rho(r),0,h$rmax)$value,h$mass(h$rmax),h$mass(2*h$rmax),h$mtot,stoptext='mass')
  cooltools::is.equal(h$potential(Inf),0,stoptext='potential')
  cooltools::is.equal(h$vesc(h$rmax)/sqrt(2),h$vcirc(h$rmax),sqrt(h$gravity*h$mtot/h$rmax),stoptext='velocity')
  cooltools::is.equal(h$quantile(0),0,stoptext='quantile q(0)')
  cooltools::is.equal(h$quantile(1),h$rmax,stoptext='quantile q(1)')

  if (plots) {

    mycurve = function(x,...) {
      graphics::curve(x,xlab='Radius r',n=1e3,...)
      if (rplot>h$rmax) graphics::abline(v=h$rmax,lty=3)
    }

    fct = function(x) ifelse(x<=h$rmax,h$rho(x),NA)
    mycurve(fct,rplot*1e-3,rplot,ylab='3D-density',log='xy')
    mycurve(h$mass,0,rplot,ylab='Enclosed mass')
    mycurve(h$potential,0,rplot,ylab='Potential')
    mycurve(h$vcirc,0,rplot,ylab='Circular velocity (1D-dispersion)')
    mycurve(h$sigma,0,rplot,add=TRUE,lty=2)
    n = 1e6
    r = h$random(n)
    r = r[r<rplot]
    graphics::hist(r,breaks=seq(0,rplot,length=101),xlim=c(0,rplot),main=sprintf('Histogram of 1e6 random radii in a %s profile',h$model))
    fct2 = function(x) h$pdf(x)*n*rplot/100
    mycurve(fct2,add=TRUE,col='red',lwd=1.5)

    invisible(h)

  } else {

    return(h)

  }

}
