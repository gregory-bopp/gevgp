
#' @title Automatic proposal variance tuning update
#' @export
update_var <- function(cur_var, acpt_rt, opt_rt, gamma1) {
  exp(log(cur_var) + gamma1 * (acpt_rt - opt_rt))
}

#' @export
call_to_list <- function(list_or_call){
  if(class(list_or_call) == "call"){
    return(eval(list_or_call))
  }
}

#' Add elements of sup_list to sub_list
#' @param sup_list high-priority list (will overwrite sub_list values)
#' @param sub_list low-priority list (will be overwritten by sup_list values)
#'
#' @return combine sub_list and sup_list. Any element that is in both sup_list
#' and sub_list will be overwritten by sup_list
#' @export
#'
#' @examples
#' const <- list(x = 1)
#' my_const <- list(y = 3)
#' update_list(my_const, const)
update_list <- function(sup_list, sub_list){
  add_names <- names(sup_list)
  if(!is.null(add_names)){
    for(name in add_names){
      sub_list[[name]] <- sup_list[[name]]
    }
  }
  return(sub_list)
}

#' @export
cond_sim_y <-
  function(Y,
           miss_ind_mat,
           coord,
           gp_var,
           gp_scale,
           gp_smooth) {
    nrep <- ncol(Y)
    RandomFields::RFoptions(pch = "")
    RandomFields::RFoptions(spConform = FALSE)
    model <- RandomFields::RMmatern(
      var = gp_var,
      scale = gp_scale,
      nu = gp_smooth,
      notinvnu = TRUE
    )
    for (i in 1:nrep) {
      miss_i <- miss_ind_mat[, i]
      if (any(miss_i)) {
        Y[miss_i, i] <- RandomFields::RFsimulate(
          n = 1,
          model = model,
          x = coord[miss_i, 1],
          y = coord[miss_i, 2],
          data = data.frame(x = coord[!miss_i, 1],
                            y = coord[!miss_i, 2],
                            G = Y[!miss_i, i])
        )
      }
    }
    return(Y)
  }

#' @export
llik_gp <- function(Y,
                obs_coord,
                gp_var,
                gp_scale,
                gp_smooth) {
  RandomFields::RFoptions(spConform = FALSE)
  model <-
    RandomFields::RMmatern(
      var = gp_var,
      scale = gp_scale,
      nu = gp_smooth,
      notinvnu = TRUE
    )
  return(RandomFields::RFlikelihood(model, x = obs_coord, data = Y)$loglikelihood)
}

#' @export
update_gp_par <-
  function(ll_cur,
           Y,
           obs_coord,
           which_par = c("gp_var", "gp_scale", "gp_smooth"),
           gp_var,
           gp_scale,
           gp_smooth,
           prop_var,
           par_range = c(0, Inf)) {
    cur_par <- eval(parse(text = which_par))
    prop_par <- rnorm(1, cur_par, sd = sqrt(prop_var))
    acpt <- 0
    if ((prop_par > par_range[1]) & (prop_par < par_range[2])) {
      if (which_par == "gp_var") {
        ll_prop <- llik_gp(Y, obs_coord, prop_par, gp_scale, gp_smooth)
      }
      if (which_par == "gp_scale") {
        ll_prop <- llik_gp(Y, obs_coord, gp_var, prop_par, gp_smooth)
      }
      if (which_par == "gp_smooth") {
        ll_prop <- llik_gp(Y, obs_coord, gp_var, gp_scale, prop_par)
      }
      if (log(runif(1)) < ll_prop - ll_cur) {
        acpt <- 1
        cur_par <- prop_par
        ll_cur <- ll_prop
      }
    }
    return(list(acpt = acpt,
                cur_par = cur_par,
                ll_cur = ll_cur
                ))
  }


