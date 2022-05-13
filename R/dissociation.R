#' Moment difference of two point sets
#'
#' @importFrom cooltools sphericalharmonics runif3
#'
#' @description Evaluates the multipole moment differences of two 3D point sets; designed to study large-scale dissociation between dark matter and gas in simulations.
#'
#' @param x n-by-3 matrix containing the first set of points (e.g. dark matter). It can also be a n-by-2 matrix, if the line-of-sight coordinate is missing (as in observations). In this case, a correction is applied to estimate the 3D dissociation indices.
#' @param y m-by-3 matrix containing the second set of points (e.g. gas). It can also be a m-by-2 matrix, if the line-of-sight coordinate is missing (as in observations). In this case, a correction is applied to estimate the 3D dissociation indices.
#'
#' @return list of dissociation indices and moments
#'
#' @author Danail Obreschkow
#'
#' @examples
#' # Generate a mock bullet-cluster and output its main dissociation index
#' x.gas = cooltools::runif3(4e3)
#' x.cdm = t(cbind(t(cooltools::runif3(5e3))-c(1,0,0),t(cooltools::runif3(5e3))+c(1,0,0)))
#' cat(sprintf('Primary dissociation index = %.2f\n',dissociation(x.cdm,x.gas)$D2))
#' \dontrun{
#' sn = snapshot(list(x.gas,x.cdm))
#' plot(sn)
#' }
#'
#' @export
dissociation = function(x,y) {

  lmax = 4

  # ensure that input is always a matrix
  x = as.matrix(x)
  y = as.matrix(y)

  # check if coordinates are 2D or 3D
  d = dim(x)[2]
  if (dim(y)[2]!=d) stop('x and y must have the same dimension')
  if (d==2) {
    x = cbind(x,0)
    y = cbind(y,0)
  }

  # number index
  nx = dim(x)[1]
  ny = dim(y)[1]
  numberindex = ifelse(nx+ny>0, (nx-ny)/(nx+ny), 0)

  D.index = P.index = array(0,lmax+1)

  if (nx>0 & ny>0) {

    # center particles to mid-point between CM of mass of DM and gas
    m = c(rep(1/nx,nx),rep(1/ny,ny))
    cg = colSums(rbind(x,y)*m)/sum(m) # geometric centre
    x = t(t(x)-cg)
    y = t(t(y)-cg)

    # radius vector
    rx = sqrt(x[,1]^2+x[,2]^2+x[,3]^2)
    ry = sqrt(y[,1]^2+y[,2]^2+y[,3]^2)
    rxmean = mean(rx)
    rymean = mean(ry)

    # moment analysis
    for (l in seq(0,lmax)) {
      m = seq(-l,l)
      fx = fy = rep(NA,2*l+1)
      for (i in seq(2*l+1)) {
        fx[i] = mean(cooltools::sphericalharmonics(l,m[i],x)*rx)
        fy[i] = mean(cooltools::sphericalharmonics(l,m[i],y)*ry)
      }
      if (d==2) {
        fx[l+1] = 0
        fy[l+1] = 0
      }
      prefactor = sqrt(4*pi/(2*l+1))
      kx = sqrt(sum(abs(fx)^2))
      ky = sqrt(sum(abs(fy)^2))
      dk = sqrt(sum(abs(fx-fy)^2))
      P.index[l+1] = prefactor*dk/(rxmean+rymean) # moment of difference
      D.index[l+1] = prefactor*(kx-ky)/max(rxmean,rymean) # difference of moments
    }
  }

  # return list
  return(list(numberindex = numberindex,
              D0 = D.index[1], # monopole difference
              D1 = D.index[2], # dipole difference
              D2 = D.index[3], # quadrupole difference
              D3 = D.index[4], # octupole difference
              D4 = D.index[5], # hexadecapole difference
              P0 = P.index[1], # monopole difference
              P1 = P.index[2], # dipole difference
              P2 = P.index[3], # quadrupole difference
              P3 = P.index[4], # octupole difference
              P4 = P.index[5]))# hexadecapole difference


}
