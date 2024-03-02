#' @export
cond_loglik_one_it <-
  function(Z_obs,
           Z_hold,
           obs_coord,
           hold_coord,
           mcmc_obj,
           mcmc_it,
           return_log = T){

    # Do search and replace on this
    samples <- mcmc_obj$samples



    RandomFields::RFoptions(pch = "")
    RandomFields::RFoptions(spConform = FALSE)
    # Constants
    nloc_obs <- nrow(Z_obs)
    nrep_obs <- ncol(Z_obs)
    nloc_hold <- nrow(Z_hold)
    nrep_hold <- ncol(Z_hold)
    nloc_all <- nloc_obs + nloc_hold
    all_coord <- rbind(obs_coord, hold_coord)
    # Transform obs to normal scale
    Y_obs <- qnorm(stablemix::pevdM(Z_obs,
                                    loc = matrix(samples$gev_loc[mcmc_it,], nrow = nloc_obs, ncol = nrep_obs),
                                    scale= matrix(samples$gev_scale[mcmc_it,], nrow = nloc_obs, ncol = nrep_obs),
                                    shape = matrix(samples$gev_shape[mcmc_it], nrow = nloc_obs, ncol = nrep_obs)))

    # Impute missing with post-pred draws
    Y_obs[mcmc_obj$miss_ind_mat] <- mcmc_obj$post_pred_Y[mcmc_it, ]
    if(any(is.infinite(Y_obs))){
      stop("Observations are Inf after transformation")
    }
    # Models for GEV parameter GPs
    gev_loc_mod <- RandomFields::RMmatern(
      var = samples$gev_loc_var[mcmc_it],
      scale = samples$gev_loc_scale[mcmc_it],
      nu = 1/2,
      notinvnu = TRUE
    )
    gev_scale_mod <- RandomFields::RMmatern(
      var = samples$gev_scale_var[mcmc_it],
      scale = samples$gev_scale_scale[mcmc_it],
      nu = 1/2,
      notinvnu = TRUE
    )
    gp_mod <- RandomFields::RMmatern(
      var = 1,
      scale = samples$gp_scale[mcmc_it],
      nu = samples$gp_smooth[mcmc_it],
      notinvnu = TRUE
    )
    # Make conditional draws of GEV GPs at holdout locs
    gev_loc_hold <- RandomFields::RFsimulate(
      n = 1,
      model = gev_loc_mod,
      x = hold_coord[,1],
      y = hold_coord[,2],
      data = data.frame(x = obs_coord[,1],
                        y = obs_coord[,2],
                        G = samples$gev_loc[mcmc_it, ]))
    gev_scale_hold <- exp(RandomFields::RFsimulate(
      n = 1,
      model = gev_scale_mod,
      x = hold_coord[,1],
      y = hold_coord[,2],
      data = data.frame(x = obs_coord[,1],
                        y = obs_coord[,2],
                        G = log(samples$gev_scale[mcmc_it, ]))))

    # Transform holdout GEV data to GP scale
    Y_hold <- qnorm(stablemix::pevdM(Z_hold,
                                     loc = matrix(gev_loc_hold, nrow = nloc_hold, ncol = nrep_hold),
                                     scale= matrix(gev_scale_hold, nrow = nloc_hold, ncol = nrep_hold),
                                     shape = matrix(samples$gev_shape[mcmc_it], nrow = nloc_hold, ncol = nrep_hold)))

    # Conditional MVN distributions -------------------------------------------
    Zcov <- RFcovmatrix(gp_mod, all_coord)
    Zcov_obs <- Zcov[1:nloc_obs, 1:nloc_obs]
    Zcov_hold <- Zcov[(nloc_obs + 1):nloc_all, (nloc_obs + 1):nloc_all]
    Zcov_hold_obs <- Zcov[(nloc_obs + 1):nloc_all, 1:nloc_obs]

    C12_i22 <- Zcov_hold_obs%*%solve(Zcov_obs)
    cond_cov <- Zcov_hold - C12_i22%*%t(Zcov_hold_obs)

    # Calculate likelihoods
    ll <- 0
    for(i in 1:nrep_hold){
      cond_mn <- C12_i22%*%Y_obs[,i]
      is_missing <- is.na(Y_hold[,i])
      cond_cov_sub <- (cond_cov[!is_missing, !is_missing] + t(cond_cov[!is_missing, !is_missing]))/2
      ll <- ll + mvtnorm::dmvnorm(Y_hold[!is_missing, i], mean = cond_mn[!is_missing],
                       sigma = cond_cov_sub, log = T)
    }
    if(return_log){
      return(ll)
    }
    return(exp(ll))
  }

#'@export
cond_loglik <- function(Z_obs,
                        Z_hold,
                        obs_coord,
                        hold_coord,
                        mcmc_obj,
                        mcmc_it,
                        burnin = 0,
                        thin_int = 1,
                        return_log = T){
  nmcmc <- length(mcmc_obj$samples$gp_scale)
  sub_seq <- seq(burnin + 1, nmcmc, by = thin_int)
  ll <- rep(NA_real_, length(sub_seq))
  cnt <- 0
  for(i in sub_seq){
    cnt <- cnt + 1
    ll[cnt] <- cond_loglik_one_it(
                       Z_obs,
                       Z_hold,
                       obs_coord,
                       hold_coord,
                       mcmc_obj,
                       mcmc_it = i,
                       return_log)
  }
  return(ll)
}


