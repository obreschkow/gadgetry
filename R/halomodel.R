##' Generate standard spherical mass distributions in astrophysics
#'
#' @importFrom graphics hist curve abline plot
#' @importFrom stats runif integrate optim
#' @importFrom cooltools midseq bindata is.equal
#'
#' @description Generates potential-density pairs and auxiliary functions associated with spherically symmetric mass distributions commonly used in astrophysics. The units of all outputs are the same as those of \code{rmax}, \code{rscale}, \code{mtot} and \code{gravity}.
#'
#' @param model name of model. Currently available models: "uniform" (spherical top-hat), "isothermal" (singular isothermal sphere with logarithmic potential), "hernquist" (Hernquist 1990), .
#' @param rmax optional truncation radius. If not given, a natural model-dependent value is assigned, which is infinite for models where the mass converges in this case.
#' @param rscale scale radius for models, which admit a characteristic scale in addition to a possible truncation radius. For all other models, this parameter is ignored. \code{rmax} and \code{rscale} should be provided in the same units.
#' @param mtot total mass enclosed by the truncation radius
#' @param gravity gravitational constant
#'
#' @return list of functions for the selected model
#'
#' @author Danail Obreschkow
#'
#' @examples
#'
#' # sample n points from Hernquist distribution
#' n = 1e6
#' h = halomodel('hernquist')
#' set = h$rradvel(n)
#'
#' # plot radial density of sample
#' xlim = c(0.01,100)
#' bin = cooltools::bindata(set$r,bins=10^seq(-2,2,by=0.05),method='custom')
#' bin$dV = 4*pi/3*(bin$xright^3-bin$xleft^3)
#' bin$rho = bin$count/bin$dV/n
#' graphics::plot(bin$xmedian,bin$rho,type='p',pch=16,cex=0.7,col='blue',
#'                xlim=xlim,ylim=c(1e-9,2e1),xaxs='i',yaxs='i',log='xy',xlab='Radius',ylab='Density')
#'
#' # compare radial density profile to analytic profile
#' graphics::curve(h$rho(x),xlim[1],xlim[2],add=TRUE)
#' graphics::abline(v=1,lty=3)
#'
#' # plot velocity distribution in a small radius interval
#' rlim = c(0.95,1.05) # interval centered around the scale radius r=1
#' xlim = c(0,h$vesc(rlim[1]))
#' sel = set$r>=rlim[1] & set$r<=rlim[2]
#' graphics::hist(set$v[sel],breaks=seq(0,h$vesc(rlim[1]),length=50),col='lightblue',border=NA,
#'                xaxs='i',yaxs='i',freq=FALSE,ylim=c(0,2.2),main=NA,
#'                xlab='Velocity norm',ylab='Probablity density')
#'
#' # compare against non-Maxwellian analytical solution
#' fvel = function(v) {
#'   g = function(v) h$dradvel(mean(rlim),v)
#'   gnorm = integrate(g,0,h$vesc(rlim[1]))$value
#'   return(g(v)/gnorm)
#' }
#' graphics::curve(fvel(x),0,1.5,1e3,add=TRUE)
#' x = cooltools::midseq(0,sqrt(2),1e4)
#' graphics::abline(v=sum(x*fvel(x))/sum(fvel(x))) # mean
#'
#' # compare against Jeans approximation
#' sigma = h$sigma(mean(rlim))
#' dmaxwell = function(x) sqrt(2/pi)*x^2/sigma^3*exp(-x^2/2/sigma^2)
#' graphics::curve(dmaxwell(x),add=TRUE,lty=2)
#' graphics::abline(v=sigma*sqrt(8/pi),lty=2) # mean of Maxwell distribution
#'
#' @export

halomodel = function(model='hernquist', rmax=NULL, rscale=1, mtot=1, gravity=1) {

  # Step 0
  # Some input checks and initialisations

  # some basic input checks
  if (rscale<=0) stop('rscale must be positive')
  if (mtot<=0) stop('mtot must be positive')
  if (gravity<=0) stop('gravity must be positive')

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
    f$dradius = function(r) ifelse(r<=f$rmax,normalization*2*r/(1+r)^3,0)

    # Enclosed mass M(<r)=integral(f$dradius,0,r)
    f$mass = function(r) ifelse(is.infinite(r),1,pmax(0,pmin(1,normalization*r^2/(1+r)^2)))

    # Circular velocity v=sqrt(M(<r)/r)
    f$vcirc = function(r) ifelse(r<=f$rmax,sqrt(normalization*(1/r-(2+1/r)/(r+1)^2)),sqrt(1/r))

    # Potential phi(r)=integral(M(<r)/r^2,0,r)+const, such that phi(Inf)=0
    const = ifelse(is.infinite(f$rmax),-1,-(normalization*f$rmax^2+f$rmax+1)/(f$rmax*(f$rmax+1)))
    f$potential = function(r) ifelse(r<=f$rmax & !is.infinite(r),normalization*r/(1+r)+const,-1/r)

    # Isotropic velocity dispersion from Jeans equation (for r<=f$rmax; NA otherwise)
    # d(rho*s)/dr = -rho*dphi/dr; sigma=sqrt(s)
    # NB: choose the integration constant, such that s(r) does not diverge anywhere
    variance = function(r) r*(1+r)^3*log((1+r)/r)-(r*(25+52*r+42*r^2+12*r^3))/(12*(1+r)) # for rmax=Inf, runs into floating-point issues at r>1e2, seriously r>1e3
    f$sigma = function(r) ifelse(r<=f$rmax,sqrt(normalization*ifelse(r<=200,variance(r),exp(-1.6152-log(r)))),NA)

    # Quantile function for 1D PDF of radius
    f$qradius = function(p) (p+sqrt(normalization*p))/(normalization-p)

    # Optional mass density per unit of 6D phase space, as a function of energy, for isotropic velocity distribution. This can be used to construct more accurate,
    # non-Maxwellian velocity distributions, not accounted for by the Jeans equation
    if (is.infinite(f$rmax)) {
      f$f6D = function(x) {
        # DF from eq. (17), Hernquist 1990
        x[x>0 | x<(-1)] = 0
        return((3*asin(sqrt(abs(-x)))+sqrt(-x-x^2)*(1+2*x)*(8*x^2+8*x-3))/(8*sqrt(2)*pi^3*(1+x)^(5/2)))
      }
    }

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
    f$dradius = function(r) ifelse(r<=1,3*r^2,0)

    # Enclosed mass M(<r)=integral(f$dradius,0,r)
    f$mass = function(r) pmin(1,pmax(0,r)^3)

    # Circular velocity v=sqrt(M(<r)/r)
    f$vcirc = function(r) ifelse(r<=1,r,sqrt(1/r))

    # Potential phi(r)=integral(M(<r)/r^2,0,r)+const, such that phi(Inf)=0
    f$potential = function(r) ifelse(r<=1,r^2/2-3/2,-1/r)

    # Isotropic velocity dispersion from Jeans equation (for r<=f$rmax; NA otherwise)
    # d(rho*s)/dr = -rho*dphi/dr; sigma=sqrt(s)
    # NB: choose the integration constant, such that s(r) does not diverge anywhere
    f$sigma = function(r) ifelse(r<=1,sqrt((1-pmin(r^2,1))/2),NA) # pmin needed despite ifelse to avoid wired warnings around r=1 in plotting

    # Quantile function for 1D PDF of radius
    f$qradius = function(p) p^(1/3)

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
    f$dradius = function(r) pmax(0,pmin(1,ceiling(r)))

    # Enclosed mass M(<r)=integral(f$dradius,0,r)
    f$mass = function(r) pmax(0,pmin(1,r))

    # Circular velocity v=sqrt(M(<r)/r)
    f$vcirc = function(r) ifelse(r<=1,1,sqrt(1/r))

    # Potential phi(r)=integral(M(<r)/r^2,0,r)+const, such that phi(Inf)=0
    f$potential = function(r) ifelse(r<=1,log(r)-1,-1/r)

    # Isotropic velocity dispersion from Jeans equation (for r<=f$rmax; NA otherwise)
    # d(rho*s)/dr = -rho*dphi/dr; sigma=sqrt(s)
    # NB: choose the integration constant, such that s(r) does not diverge anywhere
    f$sigma = function(r) ifelse(r<=1,1/sqrt(2),NA)

    # Quantile function for 1D PDF of radius
    f$qradius = function(p) p

  } else {

    stop('unknown model')

  }

  # Step 2
  # Post-processing and checks of non-dimensional mass models

  # check if all required properties have been assigned
  required.names = c('model','rmax','rho','mass','vcirc','sigma','potential','dradius','qradius')
  optional.names = c('rscale','f6D')
  allowed.names = c(required.names,optional.names)
  actual.names = names(f)
  if ('rscale'%in%actual.names) actual.names=actual.names[-which(actual.names=='rscale')]
  if (!all(required.names%in%actual.names)) stop('not all required properties exist')
  if (!all(actual.names%in%allowed.names)) stop('more than the allowed properties exist')

  # Derive a few more quantities in dimensionless units using general equations for any
  # spherical mass distribution; also run several checks

  # check mass distribution
  if (f$mass(0)!=0) stop('mass distribution error')

  # check quantile function
  cooltools::is.equal(f$qradius(0),0,stoptext='quantile q(0)')
  cooltools::is.equal(f$qradius(1),f$rmax,stoptext='quantile q(1)')

  # Complete dpqr functions for radius
  f$pradius = f$mass
  f$rradius = function(n) f$qradius(stats::runif(n))

  # Local escape velocity [equation true for any spherical mass distribution]
  f$vesc = function(r) sqrt(-2*f$potential(r))

  # Radial acceleration a=-M(<r)/r^2
  f$acc = function(r) ifelse(r==0,0,f$mass(r)/r^2)

  # Total potential energy, while checking consistency of potential-density pair
  f1 = function(r) -f$mass(r)*f$dradius(r)/r
  f2 = function(r) -(f$mass(r)/r)^2/2
  f3 = function(r) f$potential(r)*f$dradius(r)/2
  u1 = stats::integrate(f1,0,f$rmax)$value
  u2 = stats::integrate(f2,0,f$rmax)$value-1/2/f$rmax # the second term is the integral from rmax to infinity
  u3 = stats::integrate(f3,0,f$rmax)$value
  cooltools::is.equal(u1,u2,u3,stoptext='potential energy')
  f$epot = u3

  # Total kinetic energy for pure isotropic dispersion
  y = function(r) 3*f$sigma(r)^2/2*f$dradius(r) # factor three because of three dimensions
  f$ekin = stats::integrate(y,0,f$rmax)$value

  # Generate 2D probability density from phase space density as a function of energy
  if (!is.null(f$f6D)) {

    # recast energy distribution function in to a 2D probability density,
    # as a function of radius and velocity norm
    f$dradvel = function(r,v) f$f6D(v^2/2+f$potential(r))*r^2*v^2*(4*pi)^2

    # determine maximum
    fct = function(x) f$dradvel(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)])
    dradvel.max = stats::optim(c(1,0.1),fct,control=list(fnscale=-1,reltol=1e-12))$value

    # 2D random number generator
    f$rradvel = function(n) {
      r = f$rradius(n)
      vmax = f$vesc(r)
      v = rep(0,n)
      sel = seq(n)
      iteration = 0
      while (n>0) {
        v[sel] = stats::runif(n,max=vmax[sel])
        sel = sel[f$dradvel(r[sel],v[sel])<=stats::runif(n,max=dradvel.max)]
        n = length(sel)
        iteration = iteration+1
        if (iteration%%10==0) dradvel.max=dradvel.max/2
      }
      return(data.frame(r=r,v=v))
    }
  }

  # Step 3
  # Turn non-dimensional mass models in to dimensional ones, by scaling them with the custom values of
  # mass, radius and gravitational constant

  # Complete set of units
  unit$velocity = sqrt(unit$gravity*unit$mass/unit$length)
  unit$acceleration = unit$gravity*unit$mass/unit$length^2
  unit$potential = unit$gravity*unit$mass/unit$length
  unit$density = unit$mass/unit$length^3
  unit$energy = unit$gravity*unit$mass^2/unit$length

  # Scale and finalize all quantities
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
           dradius = function(x) f$dradius(x/unit$length)/unit$length,
           pradius = function(q) f$pradius(q/unit$length),
           qradius = function(p) {
             if (any(p<0) || any(p>1)) stop('p must lie inside [0,1]')
             return(f$qradius(p)*unit$length)
           },
           rradius = function(n) unit$length*f$rradius(n),
           rmax = f$rmax*unit$length)
  if (!is.null(f$rscale)) h$rscale = f$rscale*unit$length
  if (!is.null(f$f6D)) {
    h$f6D = function(x) f$f6D(x/unit$energy)/unit$length^3/unit$velocity^3
    h$dradvel = function(r,v) f$dradvel(r/unit$length,v/unit$velocity)/unit$length
    h$rradvel = function(n) {
      out = f$rradvel(n)
      return(data.frame(r=out$r*unit$length,v=out$v*unit$velocity))
    }
  }

  # check if all functions are present
  if (!all(names(f)%in%names(h))) {
    print(names(f)[!names(f)%in%names(h)])
    stop('f in h problem')
  }
  if (!all(names(h)%in%names(f))) {
    print(names(h)[!names(h)%in%names(f)])
    stop('h in f problem')
  }

  h$rscale = rscale
  h$mtot = mtot
  h$gravity = gravity

  # additional checks on dimensional quantities, just for extra safety
  cooltools::is.equal(integrate(function(r) 4*pi*r^2*h$rho(r),0,h$rmax)$value,h$mass(h$rmax),h$mass(2*h$rmax),h$mtot,stoptext='mass')
  cooltools::is.equal(h$potential(Inf),0,stoptext='potential')
  cooltools::is.equal(h$vesc(h$rmax)/sqrt(2),h$vcirc(h$rmax),sqrt(h$gravity*h$mtot/h$rmax),stoptext='velocity')
  cooltools::is.equal(h$qradius(0),0,stoptext='quantile q(0)')
  cooltools::is.equal(h$qradius(1),h$rmax,stoptext='quantile q(1)')

  return(h)

}
