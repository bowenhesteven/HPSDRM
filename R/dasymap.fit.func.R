################################################################################# 
## R code for "Incorporating spatial autocorrelation in dasymetric             ## 
## analysis: A hierarchical poisson spatial disaggregation regression model"   ##
## Authors: Bowen He, Jonathan M. Gilligan, Janey V. Camp                      ## 
#################################################################################

#########################################
####Dasymetric Mapping Fit Functions#####
#########################################

dasymap_make_model_object <- function(data, 
                                      priors = NULL, 
                                      family = 'poisson', 
                                      link = 'log', 
                                      pcfield = FALSE,
                                      gausfield = TRUE,
                                      car = TRUE,
                                      carmatrix,
                                      iid = FALSE,
                                      silent = TRUE) {
  
  
  # Check that binomial model has sample_size values supplied
  if(family == 'binomial') {
    if(sum(is.na(data$polygon_data$N)) != 0) {
      stop("There are NAs in the sample sizes. These must be supplied for a binomial likelihood")
    }
  }
  
  if(family == 'gaussian') {
    family_id = 0
  } else if(family == 'binomial') {
    family_id = 1
  } else if(family == 'poisson') {
    family_id = 2
  } else {
    stop(paste(family, "is not a valid likelihood"))
  }
  
  if(link == 'logit') {
    link_id = 0
  } else if(link == 'log') {
    link_id = 1
  } else if(link == 'identity') {
    link_id = 2
  } else {
    stop(paste(link, "is not a valid link function"))
  }
  
  if(family == 'gaussian' & iid) {
    warning('You are using both a gaussian likeihood and an iid effect. Using both of these is redundant as they are 
            having the same effect on the model. Consider setting iid = FALSE.')
  }
  
  if(is.null(data$mesh)) {
    stop('Your data object must contain an INLA mesh.')
  }
  
  nu = 1
  # Sort out mesh bits
  spde <- (INLA::inla.spde2.matern(data$mesh, alpha = nu + 1)$param.inla)[c("M0", "M1", "M2")]	
  Apix <- INLA::inla.mesh.project(data$mesh, loc = data$coordsForFit)$A
  n_s <- nrow(spde$M0)
  
  cov_matrix <- as.matrix(data$covariate_data[, -c(1:2)])
  
  # If we have exactly one column we don't have to transpose. Sure this 
  #   this could be cleaner but I don't know how.
  if(ncol(cov_matrix) == 1){
    cov_matrix <- as.matrix(apply(cov_matrix, 1, as.numeric))
  } else {
    cov_matrix <- t(apply(cov_matrix, 1, as.numeric))
  }
  
  # Construct sensible default field hyperpriors
  limits <- sp::bbox(data$polygon_shapefile)
  hypontenuse <- sqrt((limits[1,2] - limits[1,1])^2 + (limits[2,2] - limits[2,1])^2)
  prior_rho <- hypontenuse/3
  
  prior_sigma <- sd(data$polygon_data$response/mean(data$polygon_data$response))
  
  # Default priors if they are not specified
  default_priors <- list(priormean_intercept = -4.0,
                         priorsd_intercept = 2.0,
                         priormean_slope = 0.0,
                         priorsd_slope = 0.5,
                         prior_rho_min = prior_rho,
                         prior_rho_prob = 0.01,
                         prior_sigma_max = prior_sigma,
                         prior_sigma_prob = 0.01,
                         prior_logrho_gaussian_mean = 0,
                         prior_logrho_gaussian_sd = 1,
                         prior_logsigma_gaussian_mean = 0,
                         prior_logsigma_gaussian_sd = 1,
                         prior_logit_car_rho_mean = 1,
                         prior_logit_car_rho_sd = 0.5,
                         prior_car_gamma_shape = 0.3,
                         prior_car_gamma_scale = 1,
                         prior_iideffect_sd_max = 0.1,
                         prior_iideffect_sd_prob = 0.01)
  
  # Replace with any specified priors
  if(!is.null(priors)) {
    final_priors <- default_priors
    prior_names <- names(priors)
    # Check all input priors are named correctly
    if(sum(!(prior_names %in% names(final_priors))) != 0) {
      stop(paste(prior_names[!(prior_names %in% names(final_priors))], 'is not the name of a prior'))
    }
    # Check priors are not specified multiple times
    if(sum(duplicated(names(priors))) != 0) {
      message(paste(names(priors)[duplicated(names(priors))],'are specified multiple times. Will only take first value'))
      priors <- priors[!duplicated(names(priors))]
    }
    # Replace default value with new prior value
    for(i in 1:length(priors)) {
      prior_to_replace <- prior_names[i]
      final_priors[[prior_to_replace]] <- priors[[prior_to_replace]]
    }
  } else {
    final_priors <- default_priors
  }
  
parameters <- list(intercept = 2,
                     slope = rep(1, ncol(cov_matrix)),
                     log_tau_gaussian = 8,
                     iideffect = rep(0, nrow(data$polygon_data)),
                     iideffect_log_tau = 1,
                     log_sigma = 0,
                     log_rho = -0.91,
                     logit_car_rho = 0.85,
                     log_car_prec = -1.1,
                     car_u = rep(0, nrow(data$polygon_data)),
                     nodemean = rep(0, n_s))
  
input_data <- list(x = cov_matrix,
                   aggregation_values = data$aggregation_pixels,
                   Apixel = Apix,
                   spde = spde,
                   mymatrix = carmatrix,
                   startendindex = data$startendindex,
                   polygon_response_data = data$polygon_data$response,
                   response_sample_size = data$polygon_data$N,
                   family = family_id,
                   link = link_id,
                   nu = nu,
                   pcfield = as.integer(pcfield),
                   gausfield = as.integer(gausfield),
                   car = as.integer(car),
                   iid = as.integer(iid))
  
  input_data <- c(input_data, final_priors)
  
  tmb_map <- list()
  if(!pcfield & !gausfield) {
    tmb_map <- c(tmb_map, list(log_sigma = as.factor(NA),
                               log_rho = as.factor(NA),
                               nodemean = factor(rep(NA, n_s))))
  }
  
  if(!iid) {
    tmb_map <- c(tmb_map, list(iideffect_log_tau = as.factor(NA),
                               iideffect = factor(rep(NA, nrow(data$polygon_data)))))
  }
  
  if(!car) {
    tmb_map <- c(tmb_map, list( logit_car_rho = as.factor(NA),
                                log_car_prec = as.factor(NA),
                                car_u = factor(rep(NA, nrow(data$polygon_data)))))
  }
  
  if(family_id != 0) { # if not gaussian do not need a dispersion in likelihood
    tmb_map <- c(tmb_map, list(log_tau_gaussian = as.factor(NA)))
  }
  
  random_effects <- c()
  if(pcfield) {
    random_effects <- c(random_effects, 'nodemean')
  }
  if(gausfield) {
    random_effects <- c(random_effects, 'nodemean')
  }
  if(car) {
    random_effects <- c(random_effects, 'car_u')
  }
  if(iid) {
    random_effects <- c(random_effects, 'iideffect')
  }
  
  obj <- TMB::MakeADFun(
    data = input_data, 
    parameters = parameters,
    map = tmb_map,
    random = random_effects,
    silent = silent,
    DLL = "dasymap")
  
  return(obj)
}

dasymap_fit_func <- function(data, 
                             priors,
                             family, 
                             link,
                             iterations,
                             pcfield,
                             gausfield,
                             car,
                             carmatrix,
                             iid, 
                             silent){

obj <- dasymap_make_model_object(data = data,
                                 priors = priors,
                                 family = family,
                                 link = link,
                                 pcfield = pcfield,
                                 gausfield = gausfield,
                                 car = car,
                                 carmatrix = carmatrix,
                                 iid = iid,
                                 silent = silent)

message('Fitting model. This may be slow.')
opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = iterations, trace = 0))
sd_out <- TMB::sdreport(obj, getJointPrecision = TRUE)
model_output <- list(obj = obj,
                     opt = opt,
                     sd_out = sd_out,
                     data = data,
                     model_setup = list(family = family, link = link, pcfield = pcfield, gausfield = gausfield, car = car, iid = iid))
return(model_output)
}
