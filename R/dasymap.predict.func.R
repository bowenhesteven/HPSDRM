################################################################################# 
## R code for "Incorporating spatial autocorrelation in dasymetric             ## 
## analysis: A hierarchical poisson spatial disaggregation regression model"   ##
## Authors: Bowen He, Jonathan M. Gilligan, Janey V. Camp                      ## 
#################################################################################

#########################################
#Dasymetric Mapping Prediction Functions#
#########################################

predict.dasymap_model <- function(object, newdata = NULL, predict_iid = FALSE, predict_car = FALSE, N = 100, CI = 0.95, ...) {
  
  mean_prediction <- dasymap.predict_dasymap_model(object, newdata = newdata, predict_iid, predict_car)
  
  uncertainty_prediction <- dasymap.predict_uncertainty(object, newdata = newdata, predict_iid, predict_car, N, CI)
  
  prediction <- list(mean_prediction = mean_prediction,
                     uncertainty_prediction = uncertainty_prediction)
  #class(prediction) <- c('disag_prediction', 'list')
  return(prediction)
}
dasymap.predict_dasymap_model <- function(model_output, newdata = NULL, predict_iid = FALSE, predict_car = FALSE) {
  
  objects_for_prediction <- dasymap.setup_objects(model_output, newdata = newdata, predict_iid, predict_car)
  
  pars <- model_output$obj$env$last.par.best
  pars <- split(pars, names(pars))
  
  prediction <- dasymap.predict_single_raster(pars, 
                                      objects_for_prediction,
                                      link_function = model_output$model_setup$link) 
  return(prediction)
}
dasymap.predict_uncertainty <- function(model_output, newdata = NULL, predict_iid = FALSE, predict_car = FALSE, N = 100, CI = 0.95) {
  
  objects_for_prediction <- dasymap.setup_objects(model_output, newdata = newdata, predict_iid, predict_car)
  
  parameters <- model_output$obj$env$last.par.best
  
  # If we have either of the random effects, we have the jointPrecision matrix.
  #   but if we have neither, we don't get that matrix and should use the
  #   covariance matrix instead
  if(model_output$model_setup$iid | model_output$model_setup$pcfield | model_output$model_setup$car | model_output$model_setup$gausfield){
    ch <- Matrix::Cholesky(model_output$sd_out$jointPrecision)
    par_draws <- sparseMVN::rmvn.sparse(N, parameters, ch, prec = TRUE)
  } else {
    covariance_matrix <- Matrix::Matrix(model_output$sd_out$cov.fixed, sparse = TRUE)
    ch <- Matrix::Cholesky(covariance_matrix)
    par_draws <- sparseMVN::rmvn.sparse(N, parameters, ch, prec = FALSE)
  }
  
  predictions <- list()
  
  for(r in seq_len(N)) {
    
    p <- split(par_draws[r, ], names(parameters))
    
    prediction_result <- dasymap.predict_single_raster(p, 
                                               objects_for_prediction,
                                               link_function = model_output$model_setup$link) 
    
    predictions[[r]] <- prediction_result$prediction
  }
  
  predictions <- raster::stack(predictions)
  
  probs <- c((1 - CI) / 2, 1 - (1 - CI) / 2)
  predictions_ci <- raster::calc(predictions, function(x) stats::quantile(x, probs = probs, na.rm = TRUE))
  names(predictions_ci) <- c('lower CI', 'upper CI')
  
  uncertainty <- list(realisations = predictions,
                      predictions_ci = predictions_ci)
  return(uncertainty)
}
dasymap.getCoords <- function(data) {
  
  points_raster <- data$covariate_rasters[[1]]
  points_raster[is.na(points_raster)] <- -9999
  raster_pts <- raster::rasterToPoints(points_raster, spatial = TRUE)
  coords <- raster_pts@coords
  
  return(coords)
}
dasymap.getAmatrix <- function(mesh, coords) {
  
  spde <- (INLA::inla.spde2.matern(mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
  n_s <- nrow(spde$M0)						
  
  Amatrix <- INLA::inla.mesh.project(mesh, loc = as.matrix(coords))$A
  
  return(Amatrix)
}
# check and sort out new raster data.
dasymap.check_newdata <- function(newdata, model_output){
  if(is.null(newdata)) return(NULL)
  if(!is.null(newdata)){
    if(!(inherits(newdata, c('RasterStack', 'RasterBrick', 'RasterLayer')))){
      stop('newdata should be NULL or a RasterStack or a RasterBrick')
    } 
    if(!all(names(model_output$data$covariate_rasters) %in% names(newdata))){
      stop('All covariates used to fit the model must be in newdata')
    }
    # Take just the covariates we need and in the right order
    newdata <- newdata[[names(model_output$data$covariate_rasters)]]
  }
  return(newdata)
}
# Function to setup covariates, field and iid , and car objects for prediction
dasymap.setup_objects <- function(model_output, newdata = NULL, predict_iid = FALSE, predict_car = FALSE) {
  
  newdata <- dasymap.check_newdata(newdata, model_output)
  
  # Pull out original data
  data <- model_output$data
  
  # Decide which covariates to predict over
  if(is.null(newdata)){
    covariates <- data$covariate_rasters
  } else {
    covariates <- newdata
  }
  
  data$covariate_rasters <- covariates
  
  # If there is no iid effect in the model, it cannot be predicted
  if(!model_output$model_setup$iid) {
    predict_iid <- FALSE
  }
  
  # If there is no car effect in the model, it cannot be predicted
  if(!model_output$model_setup$car) {
    predict_car <- FALSE
  }
  
  if(model_output$model_setup$pcfield) {
    if(is.null(newdata)) {
      coords <- data$coordsForPrediction
    } else {
      coords <- dasymap.getCoords(data)
    }
    Amatrix <- dasymap.getAmatrix(data$mesh, coords)
    pcfield_objects <- list(coords = coords, Amatrix = Amatrix)
  } else {
    pcfield_objects <- NULL
  }
  
  if(model_output$model_setup$gausfield) {
    if(is.null(newdata)) {
      coords <- data$coordsForPrediction
    } else {
      coords <- dasymap.getCoords(data)
    }
    Amatrix <- dasymap.getAmatrix(data$mesh, coords)
    gausfield_objects <- list(coords = coords, Amatrix = Amatrix)
  } else {
    gausfield_objects <- NULL
  }
  
  if(predict_iid) {
    tmp_shp <- model_output$data$polygon_shapefile
    tmp_shp@data <- data.frame(area_id = factor(model_output$data$polygon_data$area_id))
    shapefile_raster <- raster::rasterize(tmp_shp, 
                                          model_output$data$covariate_rasters, 
                                          field = 'area_id')
    shapefile_ids <- raster::unique(shapefile_raster)
    iid_objects <- list(shapefile_raster = shapefile_raster, shapefile_ids = shapefile_ids)
  } else {
    iid_objects <- NULL
  }
  
  if(predict_car) {
    tmp_shp <- model_output$data$polygon_shapefile
    tmp_shp@data <- data.frame(area_id = factor(model_output$data$polygon_data$area_id))
    shapefile_raster <- raster::rasterize(tmp_shp, 
                                          model_output$data$covariate_rasters, 
                                          field = 'area_id')
    shapefile_ids <- raster::unique(shapefile_raster)
    car_objects <- list(shapefile_raster = shapefile_raster, shapefile_ids = shapefile_ids)
  } else{
    car_objects <- NULL
  }
  
  return(list(covariates = covariates,
              pcfield_objects = pcfield_objects,
              gausfield_objects = gausfield_objects,
              iid_objects = iid_objects, 
              car_objects = car_objects))
}
# Function to take model parameters and predict a single raster
dasymap.predict_single_raster <- function(model_parameters, objects, link_function) {
  
  # Create linear predictor
  covs_by_betas <- list()
  for(i in seq_len(raster::nlayers(objects$covariates))){
    covs_by_betas[[i]] <- model_parameters$slope[i] * objects$covariates[[i]]
  }
  
  cov_by_betas <- raster::stack(covs_by_betas)
  if(raster::nlayers(cov_by_betas) > 1){
    sum_cov_by_betas <- sum(cov_by_betas)
  } else { 
    # With only 1 covariate, there's nothing to sum. Do this to avoid warnings.
    sum_cov_by_betas <- cov_by_betas
  }
  cov_contribution <- sum_cov_by_betas + model_parameters$intercept
  
  linear_pred <- cov_contribution  
  
  if(!is.null(objects$pcfield_objects)){
    # Extract field values
    field <- (objects$pcfield_objects$Amatrix %*% model_parameters$nodemean)[, 1]
    field_ras <- raster::rasterFromXYZ(cbind(objects$pcfield_objects$coords, field))
    linear_pred <- linear_pred + field_ras
  } 
  else if(!is.null(objects$gausfield_objects)){
    # Extract field values
    field <- (objects$gausfield_objects$Amatrix %*% model_parameters$nodemean)[, 1]
    field_ras <- raster::rasterFromXYZ(cbind(objects$gausfield_objects$coords, field))
    linear_pred <- linear_pred + field_ras
  } 
  else {
    field_ras <- NULL
  }
  
  if(!is.null(objects$iid_objects)) {
    iid_ras <- objects$iid_objects$shapefile_raster
    iideffect_sd <- 1/sqrt(exp(model_parameters$iideffect_log_tau))
    for(i in seq_along(model_parameters$iideffect)) {
      iid_ras@data@values[which(objects$iid_objects$shapefile_raster@data@values == objects$iid_objects$shapefile_ids[i])] <- 
        model_parameters$iideffect[i]
      na_pixels <- which(is.na(iid_ras@data@values))
      na_iid_values <- stats::rnorm(length(na_pixels), 0, iideffect_sd)
      iid_ras@data@values[na_pixels] <- na_iid_values
    }
    if(raster::extent(iid_ras) != raster::extent(linear_pred)) {
      # Extent of prediction space is different to the original model. Keep any overlapping iid values but predict to the new extent
      raster_new_extent <- linear_pred
      raster_new_extent@data@values <- NA
      iid_ras <- raster::merge(iid_ras, raster_new_extent, ext = raster::extent(raster_new_extent))
      missing_pixels <- which(is.na(iid_ras@data@values))
      missing_iid_values <- stats::rnorm(length(missing_pixels), 0, iideffect_sd)
      iid_ras@data@values[missing_pixels] <- missing_iid_values
    }
    linear_pred <- linear_pred + iid_ras
  } else {
    iid_ras <- NULL
  }
  
  if(!is.null(objects$car_objects)) {
    car_ras <- objects$car_objects$shapefile_raster
    #car_prec <- exp(log_car_prec)
    #car_rho <- invlogit(logit_car_rho)
    for(i in seq_along(model_parameters$car_u)) {
      car_ras@data@values[which(objects$car_objects$shapefile_raster@data@values == objects$car_objects$shapefile_ids[i])] <- 
        model_parameters$car_u[i]
      #na_pixels <- which(is.na(car_ras@data@values))
      #na_iid_values <- stats::rnorm(length(na_pixels), 0, iideffect_sd)
      #iid_ras@data@values[na_pixels] <- na_iid_values
    }
    if(raster::extent(car_ras) != raster::extent(linear_pred)) {
      # Extent of prediction space is different to the original model. Keep any overlapping car values but predict to the new extent?
      raster_new_extent <- linear_pred
      raster_new_extent@data@values <- NA
      car_ras <- raster::merge(car_ras, raster_new_extent, ext = raster::extent(raster_new_extent))
      #missing_pixels <- which(is.na(iid_ras@data@values))
      #missing_iid_values <- stats::rnorm(length(missing_pixels), 0, iideffect_sd)
      #iid_ras@data@values[missing_pixels] <- missing_iid_values
    }
    linear_pred <- linear_pred + car_ras
  } else {
    car_ras <- NULL
  }
  
  if(link_function == 'logit') {
    prediction_ras <- 1 / (1 + exp(-1 * linear_pred))
  } else if(link_function == 'log') {
    prediction_ras <- exp(linear_pred)
  } else if(link_function == 'identity') {
    prediction_ras <- linear_pred
  }
  
  predictions <- list(prediction = prediction_ras, 
                      field = field_ras,
                      iid = iid_ras,
                      car = car_ras,
                      covariates = cov_contribution)
  
  return(predictions)
}

