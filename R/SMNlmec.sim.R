


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
#' @examples
#'
#' p.cens <- 0.1
#' m <- 50
#' D <- matrix(c(0.049,0.001,0.001,0.002),2,2)
#' sigma2_set <- 0.15
#' beta <- c(-2.83,-0.18)
#' nu <- 2
#' phi <- c(0.6,2)
#' nj <- rep(6,m)
#' tt <- rep(1:6,length(nj))
#' X1 <- rep(1,sum(nj))
#' X2 <- tt
#' x <- as.matrix(cbind(X1,X2))
#' Z1 <- rep(1,sum(nj))
#' Z2 <- tt
#' z <- as.matrix(cbind(Z1,Z2))
#'
#' Slash_DEC_sim <- SMNlmec.sim(m = m,x = x,z = z,tt = tt,nj = nj,beta = beta,
#'                       sigma2 = sigma2_set,D = D,phi= phi,struc ="DEC",
#'                       typeModel="Slash",p.cens = p.cens,n.cens = NULL,
#'                       cens_type="right",nu_set=nu)
#'
#' head(Slash_DEC_sim$cc)
#' head(Slash_DEC_sim$y_cc)
#'
#' y_com <- as.numeric(Slash_DEC_sim$y_cc)
#' rho_com <- as.numeric(Slash_DEC_sim$cc)
#' ycen <- y_com[rho_com == 1]
#' l_set <- length(beta)
#' q1_set <- dim(D)[1]
#' cens_N <- sum(rho_com)
#' N_com <- length(y_com)
#' ind_set <- numeric()
#' log_nj <- 0
#' for(i in 1:m){
#'   for(j in 1:nj[i]){
#'     ind_set[log_nj + j] <- i
#'   }
#'   log_nj <- log_nj + nj[i]
#' }
#' cens_nj <- numeric()
#' log_nj <- 0
#' cens_count <- 0
#'
#' for(i in 1:length(nj)){
#'   for(j in 1:nj[i]){
#'     if(rho_com[log_nj + j]==1) {
#'       cens_count <- cens_count + 1
#'     }
#'   }
#'   log_nj <- log_nj + nj[i]
#'   cens_nj[i] <- cens_count
#'   cens_count <- 0
#' }
#'
#' Slash_DEC_est <- SMNlmec.est(N_complete = N_com, N_cens = cens_N,
#'                              N_individuals = m, beta_length = l_set,
#'                              D_dim = q1_set, x_set = x, z_set = z, tt = tt,
#'                              y_complete = y_com, y_censor = ycen,
#'                              censor_vector = rho_com, nj_vector = nj,
#'                              censor_nj_vector = cens_nj,
#'                              y_ind = ind_set, dist = "Slash",
#'                              struc = "DEC", direction = "right",
#'                              thin_num = 1, chains_num = 3, iter_num = 1500,
#'                              burn_percen = 0.1, seed_set = 9955,
#'                              adapt_delta_set = 0.8)
#'
#' SMNlmec.summary(Slash_DEC_est)
#' @export



SMNlmec.sim <- function(m, x, z, tt, nj, beta, sigma2, D, phi,
                        struc = "UNC", typeModel = "Normal",
                        p.cens = 0.1, n.cens = NULL, cens_type = "right",
                        nu_set = NULL) {

  sim_obj <- MMsimu.mod(m=m,x=x,z=z,tt=tt,nj=nj,beta=beta,sigmae=sigma2,
                        D=D,phi= phi,struc=struc,typeModel=typeModel,
                        percCensu=p.cens,nivel.Censu=n.cens,
                        cens.type=cens_type,nu=nu_set)

  return(sim_obj)

}













