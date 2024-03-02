#' @export
mcmc_gp <-
  function(Y,
           coord,
           init,
           tune_var = list(gp_var = 0.1,
                           gp_scale = 0.1,
                           gp_smooth = 0.1),
           const = list(gp_var_range = c(0, Inf),
                        gp_scale_range = NULL,
                        gp_smooth_range = c(1/2,3/2)),
           niter = 100,
           thin_int = 1,
           opt_rt = 0.4,
           n_core = NULL,
           parallelize = TRUE,
           progress_bar = FALSE,
           save_on_error = TRUE,
           save_ppred = TRUE) {
    if(save_on_error){
      options(error = quote(dump.frames("gp_mcmc_dump", to.file = TRUE, include.GlobalEnv = TRUE)))
      # to load, run load("intensity_dump.rda") and debugger(intensity_dump)
    }
    # Overwrite any defaults in lists that have been passed -------------------
    # Get the default arguments in this function and set them
    form_arg <- formals(mcmc_gp)
    form_arg$const <- call_to_list(form_arg$const)
    form_arg$tune_var <- call_to_list(form_arg$tune_var)
    const <- update_list(const, form_arg$const)
    tune_var <- update_list(tune_var, form_arg$tune_var)
    if(is.null(const$gp_scale_range)){
      coord_dist <- dist(coord)
      const$gp_scale_range <- c(min(coord_dist[coord_dist>0])/2, 4*max(coord_dist))
    }

    # Define constants---------------------------------------------------------
    nrep <- ncol(Y)
    nsamples <- floor(niter / thin_int)

    # Setup parallelization ----------------------------------------------------
    if (parallelize) {
      if (is.null(n_core)) {
        n_core = parallel::detectCores()
      }
      chunk_size <- ceiling(nrep / n_core)
      par_opt <- list(chunkSize = chunk_size)
    }

    # Matrices to store mcmc samples--------------------------------------------
    samples <- list(
      gp_var = rep(NA_real_, nsamples),
      gp_scale =  rep(NA_real_, nsamples),
      gp_smooth = rep(NA_real_, nsamples)
    )

    # Automatic Tuning of Proposal variance-------------------------------------
    c0 <- 10
    c1 <- 0.8
    tune_k <- 3
    win_len <- min(niter, 50)

    # Proposal Variances--------------------------------------------------------
    prop_var <- list(
      gp_var = c(tune_var$gp_var, rep(NA_real_, win_len - 1)),
      gp_scale = c(tune_var$gp_scale, rep(NA_real_, win_len - 1)),
      gp_smooth = c(tune_var$gp_smooth, rep(NA_real_, win_len - 1))
    )

    acpt <- list(
      gp_var = rep(NA_real_,  win_len),
      gp_scale = rep(NA_real_, win_len),
      gp_smooth = rep(NA_real_, win_len)
    )

    # Get initial values--------------------------------------------
    gp_var_cur <- init$gp_var
    gp_scale_cur <- init$gp_scale
    gp_smooth_cur <- init$gp_smooth

    # Get indices of missing values -------------------------------------------
    miss_ind_mat <- is.na(Y)

    if(save_ppred){
      post_pred_Y <- matrix(NA_real_, nrow = nsamples, ncol = sum(miss_ind_mat))
    }

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
      Y <- cond_sim_y(
        Y = Y,
        miss_ind_mat = miss_ind_mat,
        coord = coord,
        gp_var = gp_var_cur,
        gp_scale = gp_scale_cur,
        gp_smooth = gp_smooth_cur
      )

      ll_cur <- llik_gp(Y, coord, gp_var_cur, gp_scale_cur, gp_smooth_cur)

      # Update gp_var -----------------------------------------------------------
      updates <- update_gp_par(ll_cur = ll_cur, Y = Y, obs_coord = coord,
                    which_par = "gp_var", gp_var = gp_var_cur,
                    gp_scale = gp_scale_cur, gp_smooth = gp_smooth_cur,
                    prop_var$gp_var, par_range = const$gp_var_range)
      ll_cur <- updates$ll_cur
      gp_var_cur <- updates$cur_par
      acpt$gp_var[pos] <- updates$acpt
      prop_var$gp_var[nxt_pos] <-
        update_var(
          prop_var$gp_var[pos],
          acpt_rt = mean(acpt$gp_var, na.rm = TRUE),
          opt_rt = opt_rt,
          gamma1 = gamma1
        )

      # Update gp_scale ---------------------------------------------------------
      updates <- update_gp_par(ll_cur = ll_cur, Y = Y, obs_coord = coord,
                               which_par = "gp_scale", gp_var = gp_var_cur,
                               gp_scale = gp_scale_cur, gp_smooth = gp_smooth_cur,
                               prop_var$gp_scale, par_range = const$gp_scale_range)
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
      updates <- update_gp_par(ll_cur = ll_cur, Y = Y, obs_coord = coord,
                               which_par = "gp_smooth", gp_var = gp_var_cur,
                               gp_scale = gp_scale_cur, gp_smooth = gp_smooth_cur,
                               prop_var$gp_smooth, par_range = const$gp_smooth_range)
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

      # Save samples -----------------------------------------------------------
      if (i %% thin_int == 0) {
        samples$gp_var[smp_loc] <- gp_var_cur
        samples$gp_scale[smp_loc] <- gp_scale_cur
        samples$gp_smooth[smp_loc] <- gp_smooth_cur
        if(save_ppred){
          post_pred_Y[smp_loc,] <- c(Y[miss_ind_mat])
        }
        smp_loc <- smp_loc  + 1
      }
    }
    acpt_rts  <- list(
      gp_var = mean(acpt$gp_var, na.rm = TRUE),
      gp_scale = mean(acpt$gp_scale, na.rm = TRUE),
      gp_smooth = mean(acpt$gp_smooth, na.rm = TRUE)
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
