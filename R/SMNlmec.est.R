




#' @title Bayesian Censored Mixed-Effects Models with Damped Exponential Correlation Structures for Scale Mixture of Normal distributions error
#' @import rstan
#' @import StanHeaders
#' @import MASS
#' @import tmvtnorm
#' @import mvtnorm
#' @import mnormt
#' @description This function fits left, right censored mixed-effects linear model, with scale mixture of normal distribution errors, using the Stan. It returns estimates, standard errors and LPML, AIC, BIC and DIC.
#' @param N_complete Total number of observations \code{N} in the response vector.
#' @param N_cens Total number of censored observations \code{N_cens} in response vector.
#' @param N_individuals Total number of subjects \code{n} in the data.
#' @param beta_length Value of p in the the \code{1 x p} coefficient vector.
#' @param D_dim Dimension of the \code{d x d} variance matrix of the random error.
#' @param x_set Design matrix of the fixed effects of order \code{N x p}.
#' @param z_set Design matrix of the random effects of order \code{N x d}.
#' @param tt Vector \code{1 x N} with the time the measurements were made, where \code{N} is the total number of measurements for all individuals. Default it's considered regular times.
#' @param y_complete Vector \code{1 x N} of the complete responses.
#' @param y_censor Vector \code{1 x N_cens} of the censored responses in the complete responses.
#' @param censor_vector Vector \code{1 x N} of the indicator vector of censored responses.
#' @param nj_vector Vector \code{1 x n} of the number of responses for each individual.
#' @param censor_nj_vector Vector \code{1 x n} of the number of censored responses for each individual.
#' @param y_ind Vector \code{1 x N} to show which subject does a specific response belong to.
#' @param dist Distribution of the random effects and random error. Available options are \code{Normal}, \code{Student} and \code{Slash}.
#' @param struc Structure of the correlation structure. Available options are \code{UNC}, \code{DEC}, \code{CAR}.
#' @param direction Direction of censoring type. Available options are \code{left} and \code{right}.
#' @param thin_num A positive integer specifying the period for saving samples. The default is 5. See more details in rstan::stan().
#' @param chains_num A positive integer specifying the number of chains generating by rstan::stan(). The default is 3.
#' @param iter_num A positive integer specifying the number of iterations for each chain (including warmup). The default is 5000.
#' @param burn_percen A percentage of the warm-up iterations in each chain the Stan. The default is 0.2.
#' @param seed_set A random seed. The default is NULL.
#' @param adapt_delta_set A parameter to control the sampler's behavior. The default is 0.8. See rstan::stan() for more details.
#' @return Return a S4 class SMNlmecfit object. Using function \code{SMNlmec.summary()} to obtain the estimation of parameters and model selection criteria. The SMNlmecfit include:
#' \item{stan_object}{A stanfit object from rstan::stan().}
#' \item{model_criteria}{A data frame includes LPML, DIC, EAIC, EBIC.}
#' \item{dist_set}{The setting of distribution of the stan model.}
#' \item{struc_set}{The setting of correlation structure of the stan model.}
#' @references Bayesian analysis of censored linear mixed-effects models for  heavy-tailed  irregularly observed repeated measures (Under review)
#' @examples
#'
#' \donttest{
#' require(rstan)
#' require(StanHeaders)
#' require(MASS)
#' require(tmvtnorm)
#' require(mvtnorm)
#' require(mnormt)
#'
#' data("UTIdata_sub")
#' data1 <- UTIdata_sub
#' subjects <- unique(data1$Patid)
#' cluster <- c(match(data1$Patid,subjects))
#' m <- length(subjects)
#' N <-length(cluster)
#' y1 <- c(log10(data1$RNA))
#' cc <- (data1$RNAcens==1)+0
#' x <- cbind((data1$Fup==0)+0, (data1$Fup==1)+0, (data1$Fup==3)+0, (data1$Fup==6)+0, (data1$Fup==9)+0, (data1$Fup==12)+0, (data1$Fup==18)+0, (data1$Fup==24)+0)
#' tem <- data1$Fup
#' nj<- rep(0,m)
#' for (j in 1:m) {
#'   nj[j] <- sum(cluster==j)
#' }
#'
#' z <- matrix(rep(1, length(y1)), ncol=1)
#'
#' ####### extract inputs for the SMNlmec.est()
#'
#' y_com <- as.numeric(y1)
#' rho_com <- as.numeric(cc)
#' ycen <- y_com[rho_com==1]
#'
#' l_set <- dim(x)[2]
#' q1_set <- 1
#'
#' cens_N <- sum(rho_com)
#' N_com <- length(y_com)
#' m_com <- m
#'
#' ind_set <- numeric()
#' log_nj <- 0
#' for(i in 1:length(nj)){
#'   for(j in 1:nj[i]){
#'     ind_set[log_nj + j] <- i
#'   }
#'   log_nj <- log_nj + nj[i]
#' }
#'
#' cens_nj <- numeric()
#' log_nj <- 0
#' cens_count <- 0
#'
#' for(i in 1:length(nj)){
#'   for(j in 1:nj[i]){
#'     if(rho_com[log_nj + j]==1) {
#'       cens_count = cens_count + 1
#'     }
#'   }
#'   log_nj <- log_nj + nj[i]
#'   cens_nj[i] <- cens_count
#'   cens_count <- 0
#' }
#'
#' UTI_T_DEC <- SMNlmec.est(N_complete = N_com, N_cens = cens_N,
#'                              N_individuals = m_com, beta_length = l_set,
#'                              D_dim = q1_set, x_set = x, z_set = z, tt = tem,
#'                              y_complete = y_com, y_censor = ycen,
#'                              censor_vector = rho_com, nj_vector = nj,
#'                              censor_nj_vector = cens_nj,
#'                              y_ind = ind_set, dist = "Student",
#'                              struc = "DEC", direction = "left",
#'                              thin_num = 1, chains_num = 3, iter_num = 3000,
#'                              burn_percen = 0.2, seed_set = 9955,
#'                              adapt_delta_set = 0.8)
#'
#' SMNlmec.summary(UTI_T_DEC)
#' }
#'
#' @export


SMNlmec.est <- function(N_complete, N_cens, N_individuals,
                        beta_length, D_dim, x_set, z_set, tt,
                        y_complete, y_censor, censor_vector, nj_vector,
                        censor_nj_vector, y_ind, dist = "Normal",
                        struc = "UNC", direction = "left",
                        thin_num = 5, chains_num = 3, iter_num = 5000,
                        burn_percen = 0.2, seed_set = NULL,
                        adapt_delta_set = 0.8) {

  N_com <- N_complete
  cens_N <- N_cens
  obs_N <- N_com - cens_N
  m <- N_individuals
  l_set <- beta_length
  rho_com <- censor_vector
  ind_set <- y_ind
  q1_set <- D_dim
  x <- x_set
  z <- z_set
  y_com <- y_complete
  ycen <- y_censor
  nj <- nj_vector
  cens_nj <- censor_nj_vector



  if(dist == "Normal") {
    if(struc == "UNC")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-UNC-N-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens=cens_N, N_obs= obs_N,
                                     n = m, l = l_set,yobs=yobs,ycen=ycen,x=x,
                                     ind=ind_set,rho=rho_com),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                              tt = tem, espac=5,
                              stanobj = stan_obj, distr="Normal",
                              depstr = "UNC", cens.type="left",LI=NULL, LS=NULL)

      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-UNC-N-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens=cens_N, N_obs= obs_N,
                                     n = m, l = l_set,yobs=yobs,ycen=ycen,x=x,
                                     ind=ind_set,rho=rho_com),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tem, espac=5,
                                     stanobj = stan_obj, distr="Normal",
                                     depstr = "UNC", cens.type="right",LI=NULL, LS=NULL)

      }
    }

    if(struc == "DEC")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-DEC-N-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                              data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                          n=m, l=l_set, q1 = q1_set,
                                          x=x, z=z, timevar = tt,
                                          y_complete = y_com, rho = rho_com,
                                          njvec =nj, cens_nj = cens_nj,
                                          ind = ind_set),
                              thin = thin_num, chains = chains_num,
                              iter = iter_num, warmup = iter_num*burn_percen,
                              seed = seed_set, control = list(adapt_delta=adapt_delta_set))


        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                              tt = tt, espac=5,
                              stanobj = stan_obj, distr="Normal",
                              depstr = "DEC", cens.type="left",LI=NULL, LS=NULL)
      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-DEC-N-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                     n=m, l=l_set, q1 = q1_set,
                                     x=x, z=z, timevar = tt,
                                     y_complete = y_com, rho = rho_com,
                                     njvec =nj, cens_nj = cens_nj,
                                     ind = ind_set),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))


        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = stan_obj, distr="Normal",
                                     depstr = "DEC", cens.type="right",LI=NULL, LS=NULL)
      }
    }

    if(struc == "CAR")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-CAR-N-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                     n=m, l=l_set, q1 = q1_set,
                                     x=x, z=z, timevar = tt,
                                     y_complete = y_com, rho = rho_com,
                                     njvec =nj, cens_nj = cens_nj,
                                     ind = ind_set),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))


        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = fit_stanCAR_N, distr="Normal",
                                     depstr = "CAR", cens.type="left",LI=NULL, LS=NULL)
      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-CAR-N-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                     n=m, l=l_set, q1 = q1_set,
                                     x=x, z=z, timevar = tt,
                                     y_complete = y_com, rho = rho_com,
                                     njvec =nj, cens_nj = cens_nj,
                                     ind = ind_set),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))


        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = fit_stanCAR_N, distr="Normal",
                                     depstr = "CAR", cens.type="right",LI=NULL, LS=NULL)

      }
    }

  }

  if(dist == "Student") {
    if(struc == "UNC")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-UNC-t-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens=cens_N, N_obs= obs_N,
                                     n = m, l = l_set,yobs=yobs,ycen=ycen,x=x,
                                     ind=ind_set,rho=rho_com),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                              tt = tem, espac=5,
                              stanobj = stan_obj, distr="Student",
                              depstr = "UNC", cens.type="left",LI=NULL, LS=NULL)
      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-UNC-t-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens=cens_N, N_obs= obs_N,
                                     n = m, l = l_set,yobs=yobs,ycen=ycen,x=x,
                                     ind=ind_set,rho=rho_com),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tem, espac=5,
                                     stanobj = stan_obj, distr="Student",
                                     depstr = "UNC", cens.type="right",LI=NULL, LS=NULL)
      }
    }

    if(struc == "DEC")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-DEC-t-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                              data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                          n=m, l=l_set, q1 = q1_set,
                                          x=x, z=z, timevar = tt,
                                          y_complete = y_com, rho = rho_com,
                                          njvec =nj, cens_nj = cens_nj,
                                          ind = ind_set),
                              thin = thin_num, chains = chains_num,
                              iter = iter_num, warmup = iter_num*burn_percen,
                              seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                              tt = tt, espac=5,
                              stanobj = stan_obj, distr="Student",
                              depstr = "DEC", cens.type="left",LI=NULL, LS=NULL)
      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-DEC-t-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                     n=m, l=l_set, q1 = q1_set,
                                     x=x, z=z, timevar = tt,
                                     y_complete = y_com, rho = rho_com,
                                     njvec =nj, cens_nj = cens_nj,
                                     ind = ind_set),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = stan_obj, distr="Student",
                                     depstr = "DEC", cens.type="right",LI=NULL, LS=NULL)
      }
    }

    if(struc == "CAR")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-CAR-t-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                              data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                          n=m, l=l_set, q1 = q1_set,
                                          x=x, z=z, timevar = tt,
                                          y_complete = y_com, rho = rho_com,
                                          njvec =nj, cens_nj = cens_nj,
                                          ind = ind_set),
                              thin = thin_num, chains = chains_num,
                              iter = iter_num, warmup = iter_num*burn_percen,
                              seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                              tt = tt, espac=5,
                              stanobj = stan_obj, distr="Student",
                              depstr = "CAR", cens.type="left",LI=NULL, LS=NULL)
      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-CAR-t-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                     n=m, l=l_set, q1 = q1_set,
                                     x=x, z=z, timevar = tt,
                                     y_complete = y_com, rho = rho_com,
                                     njvec =nj, cens_nj = cens_nj,
                                     ind = ind_set),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = stan_obj, distr="Student",
                                     depstr = "CAR", cens.type="right",LI=NULL, LS=NULL)
      }
    }

  }

  if(dist == "Slash") {
    if(struc == "UNC")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-UNC-SL-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                               data = list(N_complete = N_com, N_cens=cens_N, N_obs= obs_N,
                                           n = m, l = l_set,yobs=yobs,ycen=ycen,x=x,
                                           ind=ind_set,rho=rho_com),
                               thin = thin_num, chains = chains_num,
                               iter = iter_num, warmup = iter_num*burn_percen,
                               seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                               tt = tem, espac=5,
                               stanobj = stan_obj, distr="Slash",
                               depstr = "UNC", cens.type="left",LI=NULL, LS=NULL)

      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-UNC-SL-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens=cens_N, N_obs= obs_N,
                                     n = m, l = l_set,yobs=yobs,ycen=ycen,x=x,
                                     ind=ind_set,rho=rho_com),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tem, espac=5,
                                     stanobj = stan_obj, distr="Slash",
                                     depstr = "UNC", cens.type="right",LI=NULL, LS=NULL)

      }
    }

    if(struc == "DEC")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-DEC-SL-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                               data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                           n=m, l=l_set, q1 = q1_set,
                                           x=x, z=z, timevar = tt,
                                           y_complete = y_com, rho = rho_com,
                                           njvec =nj, cens_nj = cens_nj,
                                           ind = ind_set),
                               thin = thin_num, chains = chains_num,
                               iter = iter_num, warmup = iter_num*burn_percen,
                               seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                               tt = tt, espac=5,
                               stanobj = stan_obj, distr="Slash",
                               depstr = "DEC", cens.type="left",LI=NULL, LS=NULL)
      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-DEC-SL-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                     n=m, l=l_set, q1 = q1_set,
                                     x=x, z=z, timevar = tt,
                                     y_complete = y_com, rho = rho_com,
                                     njvec =nj, cens_nj = cens_nj,
                                     ind = ind_set),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = stan_obj, distr="Slash",
                                     depstr = "DEC", cens.type="right",LI=NULL, LS=NULL)
      }
    }

    if(struc == "CAR")   {
      if(direction == "left") {

        stan_file <- system.file("stan", "lmm-CAR-SL-censored-left.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                               data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                           n=m, l=l_set, q1 = q1_set,
                                           x=x, z=z, timevar = tt,
                                           y_complete = y_com, rho = rho_com,
                                           njvec =nj, cens_nj = cens_nj,
                                           ind = ind_set),
                               thin = thin_num, chains = chains_num,
                               iter = iter_num, warmup = iter_num*burn_percen,
                               seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                               tt = tt, espac=5,
                               stanobj = stan_obj, distr="Slash",
                               depstr = "CAR", cens.type="left",LI=NULL, LS=NULL)
      }

      if(direction == "right") {

        stan_file <- system.file("stan", "lmm-CAR-SL-censored-right.stan", package = "SMNlmec")

        stan_obj <- rstan::stan(file= stan_file,
                         data = list(N_complete = N_com, N_cens = cens_N, ycen=ycen,
                                     n=m, l=l_set, q1 = q1_set,
                                     x=x, z=z, timevar = tt,
                                     y_complete = y_com, rho = rho_com,
                                     njvec =nj, cens_nj = cens_nj,
                                     ind = ind_set),
                         thin = thin_num, chains = chains_num,
                         iter = iter_num, warmup = iter_num*burn_percen,
                         seed = seed_set, control = list(adapt_delta=adapt_delta_set))

        SMNlmec_criteria <- criteria(cc = rho_com, nj = nj, y = y_com, x = x, z = z,
                                     tt = tt, espac=5,
                                     stanobj = stan_obj, distr="Slash",
                                     depstr = "CAR", cens.type="right",LI=NULL, LS=NULL)
      }
    }
  }

  SMNlmec.est_object <- SMNlmecfit.creator(stan_obj, SMNlmec_criteria, dist, struc)

  return(SMNlmec.est_object)

}






