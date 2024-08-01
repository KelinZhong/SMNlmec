
require(lmec)
require(rstan)
require(StanHeaders)
require(MASS)
require(tmvtnorm)
require(mvtnorm)
require(mnormt)


#### just one arg for data set
#### consider list object for data set, x; z; d_column; nj_vector; y_cens;

#### set default value of some args


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
  cens_nj <-



  if(dist = "Normal") {
    if(struc = "UNC")   {
      if(direction = "left") {
        stan_obj <- rstan::stan(file='lmm-UNC-N-censored-left.stan',
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

        #### return(stan_obj[stanfit(from rstan)], criteria[list/matrix])
        #### return(list(stan_obj, criteria))

      }

      if(direction = "right") {
        stan_obj <- rstan::stan(file='lmm-UNC-N-censored-right.stan',
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

    if(struc = "DEC")   {
      if(direction = "left") {
        stan_obj <- rstan::stan(file='lmm-DEC-N-censored-left.stan',
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

      if(direction = "right") {
        stan_obj <- rstan::stan(file='lmm-DEC-N-censored-right.stan',
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

    if(struc = "CAR")   {
      if(direction = "left") {
        stan_obj <- rstan::stan(file='lmm-CAR-N-censored-left.stan',
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

      if(direction = "right") {
        stan_obj <- rstan::stan(file='lmm-CAR-N-censored-right.stan',
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

  if(dist = "Student") {
    if(struc = "UNC")   {
      if(direction = "left") {
        stan_obj <- rstan::stan(file='lmm-UNC-t-censored-left.stan',
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

      if(direction = "right") {
        stan_obj <- rstan::stan(file='lmm-UNC-t-censored-right.stan',
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

    if(struc = "DEC")   {
      if(direction = "left") {
        stan_obj <- rstan::stan(file='lmm-DEC-t-censored-left.stan',
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

      if(direction = "right") {
        stan_obj <- rstan::stan(file='lmm-DEC-t-censored-right.stan',
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

    if(struc = "CAR")   {
      if(direction = "left") {
        stan_obj <- rstan::stan(file='lmm-CAR-t-censored-UTI.stan',
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

      if(direction = "right") {
        stan_obj <- rstan::stan(file='lmm-CAR-t-censored-right.stan',
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

  if(dist = "Slash") {
    if(struc = "UNC")   {
      if(direction = "left") {
        stan_obj <- rstan::stan(file='lmm-UNC-SL-censored-left.stan',
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

      if(direction = "right") {
        stan_obj <- rstan::stan(file='lmm-UNC-SL-censored-right.stan',
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

    if(struc = "DEC")   {
      if(direction = "left") {
        stan_obj <- rstan::stan(file='lmm-DEC-SL-censored-left.stan',
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

      if(direction = "right") {
        stan_obj <- rstan::stan(file='lmm-DEC-SL-censored-right.stan',
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

    if(struc = "CAR")   {
      if(direction = "left") {
        stan_obj <- rstan::stan(file='lmm-CAR-SL-censored-left.stan',
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

      if(direction = "right") {
        stan_obj <- rstan::stan(file='lmm-CAR-SL-censored-right.stan',
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

}






