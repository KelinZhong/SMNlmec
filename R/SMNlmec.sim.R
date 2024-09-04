


#' @title Generating Censored UNC, DEC, CAR errors with Mixed Effects, for normal, student's-t and slash distribution.
#' @import lmec
#' @import MASS
#' @import tmvtnorm
#' @import mvtnorm
#' @import mnormt
#' @param m Number of individuals.
#' @param x Design matrix of the fixed effects of order \code{N x p}, corresponding to vector of fixed effects.
#' @param z Design matrix of the random effects of order\code{N x d}, corresponding to vector of random effects.
#' @param tt Vector \code{1 x N} with the time the measurements were made, where \code{N} is the total number of measurements for all individuals.
#' @param nj Vector \code{1 x m} with the number of measurements of each individual, where \code{m} is the total number of individuals.
#' @param beta Vector of values fixed effects.
#' @param sigma2 Values of the scalar of the variance matrix.
#' @param D Variance matrix of the random effects of order \code{d x d}.
#' @param phi Vector of parameter in the \code{DEC} and \code{CAR} structure. NULL for \code{UNC}, c(phi_1,phi_2) for \code{DEC} and c(phi_1,1) for \code{CAR}.
#' @param struc Structure for the simulated data. Available options are \code{UNC}, \code{DEC} and \code{CAR}.
#' @param typeModel Distribution of the simulated data. Available options are \code{Normal}, \code{Student} and \code{Slash}.
#' @param p.cens Percentage of censored measurements in the responses. The default value is 0.1.
#' @param n.cens Number of censored measurements in the responses. The default value is NULL.
#' @param cens_type The direction of cesoring. Available options are \code{left} and \code{right}.
#' @param nu_set degrees of freedom for student's-t or slash simulated data. The default value is NULL.
#' @return return list:
#' \item{cc}{Vector of censoring indicators.}
#' \item{y_cc}{Vector of responses censoring.}
#' \item{y_origin}{Vector of original complete responses.}
#' @export



SMNlmec.sim <- function(m, x, z, tt, nj, beta, sigma2, D, phi,
                        struc = "UNC", typeModel = "Normal",
                        p.cens = 0.1, n.cens = NULL, cens_type = "right",
                        nu_set) {

  sim_obj <- MMsimu.mod(m=m,x=x,z=z,tt=tt,nj=nj,beta=beta,sigmae=sigmae2,D=D,phi= phi,
                        struc="DEC",typeModel="Normal",percCensu=p.cens,
                        nivel.Censu=NULL,cens.type=cens_type,nu=NULL)

  return(sim_obj)

}













