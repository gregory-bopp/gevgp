#' @export
mcmc_gevgp <-
  function(Z,
           coord,
           init,
           tune_var = list(gp_scale = 0.1,
                           gp_smooth = 0.1,
                           gev_loc = 0.1,
                           gev_scale = 0.1,
                           gev_shape = 0.1,
                           gev_loc_var = 0.1,
                           gev_loc_scale = 0.1,
                           gev_scale_var = 0.1,
                           gev_scale_scale = 0.1
                           ),
           const = list(gp_scale_range = NULL,
                        gp_smooth_range = c(1/2,3/2),
                        gev_loc_var_range = c(0, Inf),
                        gev_loc_scale_range = NULL,
                        gev_scale_var_range =  c(0, Inf),
                        gev_scale_scale_range = NULL),
           clusters = NULL,
           niter = 100,
           thin_int = 1,
           opt_rt = 0.4,
           n_core = NULL,
           parallelize = TRUE,
           progress_bar = FALSE,
           save_on_error = TRUE,
           save_ppred = TRUE) {
    if(save_on_error){
      options(error = quote(dump.frames("gevgp_mcmc_dump", to.file = TRUE, include.GlobalEnv = TRUE)))
      # to load, run load("intensity_dump.rda") and debugger(intensity_dump)
    }
    RFoptions(spConform=FALSE)
    # Overwrite any defaults in lists that have been passed -------------------
    # Get the default arguments in this function and set them
    form_arg <- formals(mcmc_gevgp)
    form_arg$const <- call_to_list(form_arg$const)
    form_arg$tune_var <- call_to_list(form_arg$tune_var)
    const <- update_list(const, form_arg$const)
    tune_var <- update_list(tune_var, form_arg$tune_var)
    coord_dist <- dist(coord)
    default_coord_range <- c(min(coord_dist[coord_dist>0])/2, 4*max(coord_dist))
    if(is.null(const$gp_scale_range)){
      const$gp_scale_range <- default_coord_range
    }
    if(is.null(const$gev_loc_scale_range)){
      const$gev_loc_scale_range <- default_coord_range
    }
    if(is.null(const$gev_scale_scale_range)){
      const$gev_scale_scale_range <- default_coord_range
    }

    # Define constants---------------------------------------------------------
    nrep <- ncol(Z)
    nloc <- nrow(Z)
    nsamples <- floor(niter / thin_int)

    # Setup parallelization ----------------------------------------------------
    if (parallelize) {
      if (is.null(n_core)) {
        n_core = parallel::detectCores()
      }
      chunk_size <- ceiling(nrep / n_core)
      par_opt <- list(chunkSize = chunk_size)
    }


    # Define Clusters ---------------------------------------------------------
    if(is.null(clusters)){
      km <- kmeans(coord, floor(nrow(coord)/20))
      clusters <- list(gev_loc = km$cluster, gev_scale = km$cluster, gev_shape  = km$cluster)
    }
    nclus <- lapply(clusters, function(x)length(unique(x)))

    # Matrices to store mcmc samples--------------------------------------------
    samples <- list(
      gp_scale =  rep(NA_real_, nsamples),
      gp_smooth = rep(NA_real_, nsamples),
      gev_loc = matrix(NA_real_, nrow = nsamples, ncol = nloc),
      gev_scale = matrix(NA_real_, nrow = nsamples, ncol = nloc),
      gev_shape = matrix(NA_real_, nrow = nsamples),
      gev_loc_var = rep(NA, nsamples),
      gev_loc_scale = rep(NA, nsamples),
      gev_scale_var = rep(NA, nsamples),
      gev_scale_scale = rep(NA, nsamples)
    )

    # Automatic Tuning of Proposal variance-------------------------------------
    c0 <- 10
    c1 <- 0.8
    tune_k <- 3
    win_len <- min(niter, 50)

    # Proposal Variances--------------------------------------------------------
    prop_var <- list(
      gp_scale = c(tune_var$gp_scale, rep(NA_real_, win_len - 1)),
      gp_smooth = c(tune_var$gp_smooth, rep(NA_real_, win_len - 1)),
      gev_loc = matrix(NA_real_, nrow = win_len, ncol = nclus[["gev_loc"]]),
      gev_scale = matrix(NA_real_, nrow = win_len, ncol = nclus[["gev_scale"]]),
      gev_shape = c(tune_var$gev_loc, rep(NA_real_, win_len - 1)),
      gev_loc_var = c(tune_var$gev_loc_var, rep(NA_real_, win_len - 1)),
      gev_loc_scale = c(tune_var$gev_loc_scale, rep(NA_real_, win_len - 1)),
      gev_scale_var= c(tune_var$gev_scale_var, rep(NA_real_, win_len - 1)),
      gev_scale_scale = c(tune_var$gev_scale_scale, rep(NA_real_, win_len - 1))
    )
    prop_var$gev_loc[1,] <- tune_var$gev_loc
    prop_var$gev_scale[1,] <- tune_var$gev_scale

    acpt <- list(
      gp_scale = rep(NA_real_, win_len),
      gp_smooth = rep(NA_real_, win_len),
      gev_loc = matrix(NA_real_, nrow = win_len, ncol = nclus[["gev_loc"]]),
      gev_scale = matrix(NA_real_, nrow = win_len, ncol = nclus[["gev_scale"]]),
      gev_shape = rep(NA_real_, win_len),
      gev_loc_var = rep(NA_real_, win_len),
      gev_loc_scale = rep(NA_real_, win_len),
      gev_scale_var = rep(NA_real_, win_len),
      gev_scale_scale = rep(NA_real_, win_len)
    )

    # Get initial values--------------------------------------------
    gp_scale_cur <- init$gp_scale
    gp_smooth_cur <- init$gp_smooth
    gev_loc_cur <- init$gev_loc
    gev_scale_cur <- init$gev_scale
    gev_shape_cur <- init$gev_shape
    gev_loc_var_cur <- init$gev_loc_var
    gev_loc_scale_cur <- init$gev_loc_scale
    gev_scale_var_cur <- init$gev_scale_var
    gev_scale_scale_cur <- init$gev_scale_scale

    # Get indices of missing values -------------------------------------------
    miss_ind_mat <- is.na(Z)

    if(save_ppred){
      post_pred_Y <- matrix(NA_real_, nrow = nsamples, ncol = sum(miss_ind_mat))
    }


    # Transform from GEV to GP scale ------------------------------------------
    Y_cur <- qnorm(stablemix::pevdM(Z, loc = matrix(gev_loc_cur, nrow = nloc, ncol = nrep),
                                   scale = matrix(gev_scale_cur, nrow = nloc, ncol = nrep),
                                   shape = matrix(gev_shape_cur, nrow = nloc, ncol = nrep)))

    ##################################
    ######## Begin MCMC loop #########
    ##################################
    smp_loc <- 1
    if (progress_bar) {
      pb <- txtProgressBar(min = 0,
                           max = nsamples,
                           style = 3)
    }
    for (i in 1:niter) {
      if (progress_bar) {
        setTxtProgressBar(pb, smp_loc)
      }
      gamma1 <- c0 / (i + tune_k) ^ (c1)            # Automatic tuning constant
      pos <- (i - 1) %% win_len + 1
      nxt_pos <- i %% win_len + 1

      # Conditionally simulate missing values from posterior predictive dist.
      Y_cur <- cond_sim_y(
        Y = Y_cur,
        miss_ind_mat = miss_ind_mat,
        coord = coord,
        gp_scale = gp_scale_cur,
        gp_smooth = gp_smooth_cur
      )
      Z_cur <- stablemix::qevdM(pnorm(Y_cur),
                        loc = matrix(gev_loc_cur, nrow = nloc, ncol = nrep),
                        scale = matrix(gev_scale_cur, nrow = nloc, ncol = nrep),
                        shape = matrix(gev_shape_cur, nrow = nloc, ncol = nrep))

      ll_cur <- llik_gp(Y_cur, coord, 1, gp_scale_cur, gp_smooth_cur)

      # Update gp_scale ---------------------------------------------------------
      updates <- update_gp_par(ll_cur = ll_cur, Y = Y_cur, obs_coord = coord,
                               which_par = "gp_scale", gp_var = 1,
                               gp_scale = gp_scale_cur, gp_smooth = gp_smooth_cur,
                               prop_var$gp_scale[pos], par_range = const$gp_scale_range)
      ll_cur <- updates$ll_cur
      gp_scale_cur <- updates$cur_par
      acpt$gp_scale[pos] <- updates$acpt
      prop_var$gp_scale[nxt_pos] <-
        update_var(
          prop_var$gp_scale[pos],
          acpt_rt = mean(acpt$gp_scale, na.rm = TRUE),
          opt_rt = opt_rt,
          gamma1 = gamma1
        )

      # Update gp_smooth --------------------------------------------------------
      updates <- update_gp_par(ll_cur = ll_cur, Y = Y_cur, obs_coord = coord,
                               which_par = "gp_smooth", gp_var = 1,
                               gp_scale = gp_scale_cur, gp_smooth = gp_smooth_cur,
                               prop_var$gp_smooth[pos], par_range = const$gp_smooth_range)
      ll_cur <- updates$ll_cur
      gp_smooth_cur <- updates$cur_par
      acpt$gp_smooth[pos] <- updates$acpt
      prop_var$gp_smooth[nxt_pos] <-
        update_var(
          prop_var$gp_smooth[pos],
          acpt_rt = mean(acpt$gp_smooth, na.rm = TRUE),
          opt_rt = opt_rt,
          gamma1 = gamma1
        )

      # Update GEV parameters ---------------------------------------------------
      # Location
      ll_cur <- llik_gp(gev_loc_cur, obs_coord = coord, gp_var = gev_loc_var_cur,
                        gp_scale = gev_loc_scale_cur, gp_smooth = 1/2) +
        llik_gp(Y_cur, obs_coord = coord, gp_var = 1,
                gp_scale = gp_scale_cur, gp_smooth = gp_smooth_cur)  +
        ljacob(Z_cur, Y_cur,
               gev_loc = matrix(gev_loc_cur, nrow = nloc, ncol = nrep),
               gev_scale = matrix(gev_scale_cur, nrow = nloc, ncol = nrep),
               gev_shape = matrix(gev_shape_cur, nrow = nloc, ncol = nrep))
      gev_loc_prop <- gev_loc_cur
      for(j in 1:nclus[["gev_loc"]]){
        sub <- (clusters[["gev_loc"]] == j)
        cmod_cur <- RandomFields::RMexp(var = gev_loc_var_cur, scale = gev_loc_scale_cur)
        Cmat <- RandomFields::RFcovmatrix(cmod_cur, coord)
        Csub <- Cmat[sub, sub]
        Lt <- chol(Csub)
        white <- gev_loc_cur[sub]%*%solve(Lt)
        white_shifted <- white + rnorm(sum(sub), mean = 0, sd = sqrt(prop_var$gev_loc[pos, j]))
        gev_loc_prop[sub] <- white_shifted%*%Lt
        Y_prop <- qnorm(stablemix::pevdM(Z_cur, loc = matrix(gev_loc_prop, nrow = nloc, ncol = nrep),
                                                      scale = matrix(gev_scale_cur, nrow = nloc, ncol = nrep),
                                                      shape = matrix(gev_shape_cur, nrow = nloc, ncol = nrep)))
        if(any(abs(Y_prop) == Inf)){
          ll_prop <- -Inf
        }
        else{
        ll_prop <- llik_gp(gev_loc_prop, obs_coord = coord, gp_var = gev_loc_var_cur,
                          gp_scale = gev_loc_scale_cur,gp_smooth = 1/2) +
                  llik_gp(Y_prop, obs_coord = coord, gp_var = 1,
                          gp_scale = gp_scale_cur, gp_smooth = gp_smooth_cur) +
                  ljacob(Z_cur, Y_prop,
                         gev_loc = matrix(gev_loc_prop, nrow = nloc, ncol = nrep),
                         gev_scale = matrix(gev_scale_cur, nrow = nloc, ncol = nrep),
                         gev_shape = matrix(gev_shape_cur, nrow = nloc, ncol = nrep))
        }
        if(log(runif(1)) < ll_prop - ll_cur){
          acpt$gev_loc[pos, j] <- 1
          gev_loc_cur[sub] <- gev_loc_prop[sub]
          ll_cur <- ll_prop
          Y_cur <- Y_prop
        }
        else{
          acpt$gev_loc[pos, j] <- 0
          gev_loc_prop[sub]  <- gev_loc_cur[sub]
        }
        prop_var$gev_loc[nxt_pos, j ] <-
          update_var(
            prop_var$gev_loc[pos, j],
            acpt_rt = mean(acpt$gev_loc, na.rm = TRUE),
            opt_rt = opt_rt,
            gamma1 = gamma1
          )
      }

      # Log-scale
      gev_scale_prop <- gev_scale_cur
      ll_cur <- llik_gp(log(gev_scale_cur), obs_coord = coord, gp_var = gev_scale_var_cur,
                        gp_scale = gev_scale_scale_cur,gp_smooth = 1/2) +
                llik_gp(Y_cur, obs_coord = coord, gp_var = 1,
                        gp_scale = gp_scale_cur, gp_smooth = gp_smooth_cur) +
                ljacob(Z_cur, Y_cur,
                       gev_loc = matrix(gev_loc_cur, nrow = nloc, ncol = nrep),
                       gev_scale = matrix(gev_scale_cur, nrow = nloc, ncol = nrep),
                       gev_shape = matrix(gev_shape_cur, nrow = nloc, ncol = nrep))
      for(j in 1:nclus[["gev_scale"]]){
        sub <- (clusters[["gev_scale"]] == j)
        cmod_cur <- RandomFields::RMexp(var = gev_scale_var_cur, scale = gev_scale_scale_cur)
        Cmat <- RandomFields::RFcovmatrix(cmod_cur, coord)
        Csub <- Cmat[sub, sub]
        Lt <- chol(Csub)
        white <- log(gev_scale_cur[sub])%*%solve(Lt)
        white_shifted <- white + rnorm(sum(sub), mean = 0, sd = sqrt(prop_var$gev_scale[pos, j]))
        gev_scale_prop[sub] <- exp(white_shifted%*%Lt)
        Y_prop <- qnorm(stablemix::pevdM(Z_cur, loc = matrix(gev_loc_cur, nrow = nloc, ncol = nrep),
                                         scale = matrix(gev_scale_prop, nrow = nloc, ncol = nrep),
                                         shape = matrix(gev_shape_cur, nrow = nloc, ncol = nrep)))
        if(any(abs(Y_prop) == Inf)){
          ll_prop <- -Inf
        }
        else{
          ll_prop <- llik_gp(log(gev_scale_prop), obs_coord = coord, gp_var = gev_scale_var_cur,
                             gp_scale = gev_scale_scale_cur,gp_smooth = 1/2) +
                    llik_gp(Y_prop, obs_coord = coord, gp_var = 1,
                            gp_scale = gp_scale_cur, gp_smooth = gp_smooth_cur) +
                    ljacob(Z_cur, Y_prop,
                           gev_loc = matrix(gev_loc_cur, nrow = nloc, ncol = nrep),
                           gev_scale = matrix(gev_scale_prop, nrow = nloc, ncol = nrep),
                           gev_shape = matrix(gev_shape_cur, nrow = nloc, ncol = nrep))
        }
        if(log(runif(1)) < ll_prop - ll_cur){
          acpt$gev_scale[pos, j] <- 1
          gev_scale_cur[sub] <- gev_scale_prop[sub]
          ll_cur <- ll_prop
          Y_cur <- Y_prop
        }
        else{
          acpt$gev_scale[pos, j] <- 0
          gev_scale_prop[sub]  <- gev_scale_cur[sub]
        }
        prop_var$gev_scale[nxt_pos, j ] <-
          update_var(
            prop_var$gev_scale[pos, j],
            acpt_rt = mean(acpt$gev_scale, na.rm = TRUE),
            opt_rt = opt_rt,
            gamma1 = gamma1
          )
      }
      # Shape
      ll_cur <- llik_gp(Y_cur, obs_coord = coord, gp_var = 1,
                gp_scale = gp_scale_cur, gp_smooth = gp_smooth_cur) +
                ljacob(Z_cur, Y_cur,
                       gev_loc = matrix(gev_loc_cur, nrow = nloc, ncol = nrep),
                       gev_scale = matrix(gev_scale_cur, nrow = nloc, ncol = nrep),
                       gev_shape = matrix(gev_shape_cur, nrow = nloc, ncol = nrep))
        while(TRUE){
          gev_shape_prop <- rnorm(1, gev_shape_cur, sd = sqrt(prop_var$gev_shape[pos]))
          if((gev_shape_prop > -1) & (gev_shape_prop < 1)){ break}
        }
        Y_prop <- qnorm(stablemix::pevdM(Z_cur, loc = matrix(gev_loc_cur, nrow = nloc, ncol = nrep),
                                         scale = matrix(gev_scale_cur, nrow = nloc, ncol = nrep),
                                         shape = matrix(gev_shape_prop, nrow = nloc, ncol = nrep)))
        if(any(abs(Y_prop) == Inf)){
          ll_prop <- -Inf
        }
        else{
          ll_prop <- llik_gp(Y_prop, obs_coord = coord, gp_var = 1,
                            gp_scale = gp_scale_cur, gp_smooth = gp_smooth_cur) +
                    ljacob(Z_cur, Y_prop,
                           gev_loc = matrix(gev_loc_cur, nrow = nloc, ncol = nrep),
                           gev_scale = matrix(gev_scale_cur, nrow = nloc, ncol = nrep),
                           gev_shape = matrix(gev_shape_prop, nrow = nloc, ncol = nrep))
        }
        if(log(runif(1)) < ll_prop - ll_cur){
          acpt$gev_shape[pos] <- 1
          gev_shape_cur <- gev_shape_prop
          ll_cur <- ll_prop
          Y_cur <- Y_prop
        }
        else{
          acpt$gev_shape[pos] <- 0
        }
        prop_var$gev_shape[nxt_pos] <-
          update_var(
            prop_var$gev_shape[pos],
            acpt_rt = mean(acpt$gev_shape, na.rm = TRUE),
            opt_rt = opt_rt,
            gamma1 = gamma1
          )

      # Update GEV Location Covariance parameters ------------------------------
      ll_cur <- llik_gp(gev_loc_cur, coord, gev_loc_var_cur, gev_loc_scale_cur, 1/2)
      # GEV Loc Cov Var
      updates <- update_gp_par(ll_cur = ll_cur, Y = gev_loc_cur, obs_coord = coord,
                               which_par = "gp_var", gp_var = gev_loc_var_cur,
                               gp_scale = gev_loc_scale_cur, gp_smooth = 1/2,
                               prop_var$gev_loc_var[pos],
                               par_range = const$gev_loc_var_range)
      ll_cur <- updates$ll_cur
      gev_loc_var_cur <- updates$cur_par
      acpt$gev_loc_var[pos] <- updates$acpt
      prop_var$gev_loc_var[nxt_pos] <-
        update_var(
          prop_var$gev_loc_var[pos],
          acpt_rt = mean(acpt$gev_loc_var, na.rm = TRUE),
          opt_rt = opt_rt,
          gamma1 = gamma1
        )
      # GEV Loc Cov Scale
      updates <- update_gp_par(ll_cur = ll_cur, Y = gev_loc_cur, obs_coord = coord,
                               which_par = "gp_scale", gp_var = gev_loc_var_cur,
                               gp_scale = gev_loc_scale_cur, gp_smooth = 1/2,
                               prop_var$gev_loc_scale[pos],
                               par_range = const$gev_loc_scale_range)
      ll_cur <- updates$ll_cur
      gev_loc_scale_cur <- updates$cur_par
      acpt$gev_loc_scale[pos] <- updates$acpt
      prop_var$gev_loc_scale[nxt_pos] <-
        update_var(
          prop_var$gev_loc_scale[pos],
          acpt_rt = mean(acpt$gev_loc_scale, na.rm = TRUE),
          opt_rt = opt_rt,
          gamma1 = gamma1
        )

      # Update GEV Scale Covariance parameters ---------------------------------
      ll_cur <- llik_gp(log(gev_scale_cur), coord, gev_scale_var_cur,
                        gev_scale_scale_cur, 1/2)
      # GEV Scale Cov Var
      updates <- update_gp_par(ll_cur = ll_cur, Y = log(gev_scale_cur), obs_coord = coord,
                               which_par = "gp_var", gp_var = gev_scale_var_cur,
                               gp_scale = gev_scale_scale_cur, gp_smooth = 1/2,
                               prop_var$gev_scale_var[pos],
                               par_range = const$gev_scale_var_range)
      ll_cur <- updates$ll_cur
      gev_scale_var_cur <- updates$cur_par
      acpt$gev_scale_var[pos] <- updates$acpt
      prop_var$gev_scale_var[nxt_pos] <-
        update_var(
          prop_var$gev_scale_var[pos],
          acpt_rt = mean(acpt$gev_scale_var, na.rm = TRUE),
          opt_rt = opt_rt,
          gamma1 = gamma1
        )
      # GEV Scale Cov Scale
      updates <- update_gp_par(ll_cur = ll_cur, Y = log(gev_scale_cur), obs_coord = coord,
                               which_par = "gp_scale", gp_var = gev_scale_var_cur,
                               gp_scale = gev_scale_scale_cur, gp_smooth = 1/2,
                               prop_var$gev_scale_scale[pos],
                               par_range = const$gev_scale_scale_range)

      ll_cur <- updates$ll_cur
      gev_scale_scale_cur <- updates$cur_par
      acpt$gev_scale_scale[pos] <- updates$acpt
      prop_var$gev_scale_scale[nxt_pos] <-
        update_var(
          prop_var$gev_scale_scale[pos],
          acpt_rt = mean(acpt$gev_scale_scale, na.rm = TRUE),
          opt_rt = opt_rt,
          gamma1 = gamma1
        )

      # Save samples -----------------------------------------------------------
      if (i %% thin_int == 0) {
        samples$gp_scale[smp_loc] <- gp_scale_cur
        samples$gp_smooth[smp_loc] <- gp_smooth_cur
        samples$gev_loc[smp_loc,] <- gev_loc_cur
        samples$gev_scale[smp_loc,] <- gev_scale_cur
        samples$gev_shape[smp_loc] <- gev_shape_cur
        samples$gev_loc_var[smp_loc] <- gev_loc_var_cur
        samples$gev_loc_scale[smp_loc] <- gev_loc_scale_cur
        samples$gev_scale_var[smp_loc] <- gev_scale_var_cur
        samples$gev_scale_scale[smp_loc] <- gev_scale_scale_cur
        if(save_ppred){
          post_pred_Y[smp_loc,] <- c(Y_cur[miss_ind_mat])
        }
        smp_loc <- smp_loc  + 1
      }
    }
    acpt_rts  <- list(
      gp_scale = mean(acpt$gp_scale, na.rm = TRUE),
      gp_smooth = mean(acpt$gp_smooth, na.rm = TRUE),
      gev_loc <- colMeans(acpt$gev_loc, na.rm = T),
      gev_scale <- colMeans(acpt$gev_scale, na.rm = T),
      gev_shape <-mean(acpt$gev_shape, na.rm = TRUE),
      gev_loc_var <-mean(acpt$gev_loc_var, na.rm = TRUE),
      gev_loc_scale <- mean(acpt$gev_loc_scale, na.rm = TRUE),
      gev_scale_var <- mean(acpt$gev_scale_var, na.rm = TRUE),
      gev_scale_scale <- mean(acpt$gev_scale_scale, na.rm = TRUE)
    )

    if(!save_ppred){
      post_pred_Y <- NULL
    }

    return(
      list(
        samples = samples,
        prop_var = prop_var,
        acpt_rts = acpt_rts,
        const = const,
        post_pred_Y = post_pred_Y,
        miss_ind_mat = miss_ind_mat
      )
    )
  }
