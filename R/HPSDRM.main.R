################################################################################ 
## R code for "" #Incorporating spatial autocorrelation in dasymetric         ## 
## analysis: A hierarchical poisson spatial disaggregation regression model   ##
## Authors: Bowen He, Jonathan M. Gilligan, Janey V. Camp                     ## 
################################################################################

################################################################################
##### Hierarchical Poisson Spatial Disaggregation Regression Model(HPSDRM) #####
################################################################################

# In this R code file, the decennial 2020 population in the 
# Davidson County, Nashville is spatially interpolated from the 
# census tract scale to the 30m*30m*5*5 = 22500 m^2 resolution 
# grids using the proposed Hierarchical Poisson Spatial Disaggregation 
# Regression Model (HPSDRM).

# The NLCD land cover data is used as ancillary information 
# and the selection of the land cover covariates in the model 
# as the predictors is based on the results of regularized 
# regression algorithms: (1) DOS; (2) DLI; (3) DMI.

# Based on the regularized analysis results, 3 land cover 
# covariates are selected as the significant land cover 
# predictors in the proposed Hierarchical Poisson Spatial 
# Regression Model (HPSDRM):
# (1) DOS: developed open space;
# (2) DLI: developed low intensity; 
# (3) DMI: developed medium intensity.

# Load required libraries
library('sf')
library('tidyverse')
library('dplyr')
library('INLA') 
library('raster') 
library('sp') 
library('spdep')
library('spatialreg')
library('rgdal')
library('parallel')
library('Metrics')
library('Matrix')
library('exactextractr')
library('ggplot2')
library('tidycensus')
library('rprojroot')
library('units')
library('fasterize')
library('disaggregation')
library('TMB')
library('tidycensus')
library('LaplacesDemon')
proj_root <- find_root(is_rstudio_project)

#compile(paste0(proj_root,"/src/dasymap.cpp"))
dyn.load(dynlib(paste0(proj_root,"/src/dasymap")))

source(paste0(proj_root,"/R/dasymap.fit.func.R"))
source(paste0(proj_root,"/R/dasymap.predict.func.R"))
source(paste0(proj_root,"/R/dasymap.scale.func.R"))
source(paste0(proj_root,"/R/plot.prior.post.func.R"))

#### Prediction raster
# cov.dos: developed open space covariate layer for prediction
# cov.dli: developed low intensity covariate layer for prediction
# cov.dmi: developed medium intensity covariate layer for prediction


# shape: Census Tract Polygon Responses SpatialPolygons Shapes file
# mymatrix: Leroux prior precision matrix
# priors: list of parameters associated with the priors for the HPSDRM
# dis_data: prepared data lists for the HPSDRM
# fitted_model: list of fitted results for the HPSDRM
# prediction_model: list of predicted results for the HPSDRM
# scale_func:a scaling function to fulfill the pycnophylactic property of the posterior mean prediction of the HPSDRM

# Define the CRS of the study region
CRS <- "+proj=utm +zone=16 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

# download in the decennial 2020 population data
pop.tract <- get_decennial(geography = "tract", sumfile = "pl", state = 'TN', 
                           county = 'Davidson', variables = "P1_001N", 
                           year = 2020, geometry = TRUE, cache_table = TRUE) %>%
  mutate(area = st_area(.)) %>%
  st_transform(crs = CRS)

# download the Decennial 2020 blocks data of the Davidson County, Nashville
pop.block <- get_decennial(geography = "block", sumfile = "pl", state = 'TN', 
                           county = 'Davidson', variables = "P1_001N", 
                           year = 2020, geometry = TRUE, cache_table = TRUE) %>%
  mutate(area = st_area(.)) %>%
  filter(value > 0)  %>%
  st_transform(crs = CRS)

# load in the NLCD 2019 data
nlcd_raster <- raster::raster(file.path(proj_root, "data/src/NLCD-2019/nlcd-2019.img")) %>%
  projectRaster(crs = CRS, method="ngb")

# Census Tract Polygon Responses
shapes <- as(pop.tract %>% dplyr::select(GEOID, value) %>% st_transform(CRS), 'Spatial')

# Raster to points
nlcd_raster_pts <- as.data.frame(raster::rasterToPoints(nlcd_raster)) 

# covariate layer for prediction
# Developed Open Space (dos)
cov.dos <- nlcd_raster_pts %>%
  mutate(nlcd.2019 = ifelse(nlcd.2019 == 21, 1, 0)) %>%
  rasterFromXYZ() %>%
  raster::aggregate(fact = 5, fun = sum) 

# Developed Low Intensity (dli)
cov.dli <- nlcd_raster_pts %>%
  mutate(nlcd.2019 = ifelse(nlcd.2019 == 22, 1, 0)) %>%
  rasterFromXYZ() %>%
  raster::aggregate(fact = 5, fun = sum) 

# Developed Medium Intensity (dmi)
cov.dmi <- nlcd_raster_pts %>%
  mutate(nlcd.2019 = ifelse(nlcd.2019 == 23, 1, 0)) %>%
  rasterFromXYZ() %>%
  raster::aggregate(fact = 5, fun = sum) 

# DOS, DLI, and DMI as Fixed Effects Components
covariate_stack <- stack(cov.dos, cov.dli,  cov.dmi)
names(covariate_stack) <- c('dos', 'dli', 'dmi')

# obtain the polygon response neighborhood matrix
area.nb <- poly2nb(shapes)
area.nb_B <- nb2listw(area.nb, style="B", zero.policy=TRUE)

sym.matrix <- as(area.nb_B, "symmetricMatrix")
one.matrix <- rep(1, nrow(shapes))
iden.matrix <- Diagonal(x = as.numeric(sym.matrix%*%one.matrix))
mymatrix <- iden.matrix - sym.matrix

# prepare data
dis_data <- prepare_data(polygon_shapefile = shapes,
                         covariate_rasters = covariate_stack,
                         mesh.args = list(max.edge = c(1, 3),
                                          cut = 0.02,
                                          offset = c(3, 5)),
                         id_var = 'GEOID',
                         response_var = 'value',
                         ncores = 8)
#print('start model fitting')

# set up the priors for the HPSDRM
priors <- list(priormean_intercept = 0,
               priorsd_intercept = 2,
               priormean_slope = rep(0, nlayers(covariate_stack)),
               priorsd_slope = 1,
               prior_rho_min = 1.5,
               prior_rho_prob = 0.01,
               prior_sigma_max = 0.25,
               prior_sigma_prob = 0.01,
               prior_logit_car_rho_mean = logit(0.5),
               prior_logit_car_rho_sd = 15,
               prior_car_gamma_shape = 1,
               prior_car_gamma_scale = 2)

fitted_model <- dasymap_fit_func(data = dis_data,
                                 priors = priors,
                                 family = 'poisson',
                                 link = 'log',
                                 iterations = 2000,
                                 pcfield = TRUE,
                                 gausfield = FALSE,
                                 car = TRUE,
                                 carmatrix = mymatrix,
                                 iid = FALSE,
                                 silent = TRUE)

prediction_model <- predict.dasymap_model(fitted_model, predict_iid = FALSE, predict_car = TRUE)

model.result <- list(dis_data = dis_data,
                     fitted_model= fitted_model,
                     prediction_model = prediction_model)

saveRDS(model.result, file = (paste0(proj_root, "/data/gen/sensitivity_analysis/HPSDRM.res.logitlambda.mean0.sd15.prec.shape1.scale8.rds")))
#scaling function to fulfill the pycnophylactic property
pred.raster <- model.result$prediction_model$mean_prediction$prediction

# plot the scaled dasymetric mapping graph
hpsdrm_dasy_pop_raster <- scale_func(pop.tract, pred.raster)
plot(hpsdrm_dasy_pop_raster, main = 'HPSDRM Grids Population')

#plot dis.data----------------------------------------------------------------------------
#plot Polygon Response Data
polygon.ras<-rasterize(model.result$dis_data$polygon_shapefile, model.result$prediction_model$mean_prediction$prediction, field = 'value')
plot(polygon.ras, main='Polygon Response Data')

# plot covariate rasters
covariate_rasters_plot <- raster::reclassify(model.result$dis_data$covariate_rasters, cbind(0 , NA))
plot(covariate_rasters_plot)

# plot INLA mesh for Matérn spatial field
plot_mesh <- function(x, main = '', col = 'blue', lwd = 0.3, linecol = 'darkgrey', size = 1.2) {

  mesh <- x
  # extract point data
  d <- data.frame(x = mesh$loc[, 1], y = mesh$loc[, 2], type = 'evertices')
  levels(d$type) <- c('evertices', 'adata')
  d[mesh$idx$loc, 'type'] <- 'adata'
  # extract lines data.
  # mesh$graph$tv column 1, 2, 3 are points in triangles.
  # Therefore need 1 to 2, 2 to 3 and 3 to 1.
  idx = rbind(mesh$graph$tv[, 1:2, drop = FALSE],
              mesh$graph$tv[, 2:3, drop = FALSE],
              mesh$graph$tv[, c(3, 1), drop = FALSE])
  segments <- data.frame(mesh$loc[idx[, 1], 1:2], mesh$loc[idx[, 2], 1:2], type = 'bsegments')

  innerouter <- data.frame(mesh$loc[mesh$segm$bnd$idx[, 1], 1:2],
                           mesh$loc[mesh$segm$bnd$idx[, 2], 1:2],
                           type = 'cbinding', stringsAsFactors = FALSE)
  if(nrow(mesh$segm$int$idx) > 0){
    innerouter <- rbind(innerouter,
                        data.frame(mesh$loc[mesh$segm$int$idx[, 1], 1:2],
                                   mesh$loc[mesh$segm$int$idx[, 2], 1:2],
                                   type = 'dinternal'))
  } else {
    #innerouter <- rbind(innerouter,
    #                    NA)
    #innerouter[nrow(innerouter), 5] <- 'dinternal'
    innerouter$type = factor(innerouter$type, levels = c('dinternal', 'cbinding'))
  }


  names(segments) <- c('x1', 'y1', 'x2', 'y2', 'type')
  names(innerouter) <- c('x1', 'y1', 'x2', 'y2', 'type')

  segments <- rbind(segments, innerouter)


  p <- ggplot2::ggplot(data = d,
                       ggplot2::aes_string('x', 'y',
                                           colour = 'type',
                                           size = 'type')) +
    ggplot2::geom_segment(data = segments,
                          ggplot2::aes_string(x = 'x1', y = 'y1', xend = 'x2', yend = 'y2')) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = 'none')
  #stroke
  p <- p +
    ggplot2::scale_colour_manual(values = c(col, linecol, 'black', 'black', 'black'), drop = FALSE) +
    ggplot2::scale_size_manual(values = c(size, lwd, 1.3, 1.3, 0), drop = FALSE) +
    ggtitle(main) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 50),
          axis.title.x = element_text(size=35),
          axis.title.y = element_text(size=35),
          axis.text.y = element_text(size = 35),
          axis.text.x = element_text(size = 35))

  return(invisible(p))
}
mesh.plot<-plot_mesh(model.result$dis_data$mesh, main ='INLA mesh for Matérn spatial field', lwd = 0.3, size = 0.5, col = 'blue')
mesh.plot

# plot pred.data------------------------------------------------------------------------------
# covaraite contribution
covariate_contribution_raster <- raster::mask(model.result$prediction_model$mean_prediction$covariates, model.result$dis_data$polygon_shapefile, inverse = FALSE)
plot(covariate_contribution_raster, main = 'Land Cover Covariates Contributions', col = brewer.pal(11, "PRGn"),
     zlim = c(-max(abs(covariate_contribution_raster@data@values)[!is.na(abs(covariate_contribution_raster@data@values))]),
              max(abs(covariate_contribution_raster@data@values)[!is.na(abs(covariate_contribution_raster@data@values))]))
)

# CAR effect raster plot
pred.car.raster <- model.result$prediction_model$mean_prediction$car
plot(pred.car.raster, main = 'Predicted Mean Leroux CAR Random Effects', col = brewer.pal(11, "RdBu"),
     zlim = c(-max(abs(pred.car.raster@data@values)[!is.na(abs(pred.car.raster@data@values))]),
              max(abs(pred.car.raster@data@values)[!is.na(abs(pred.car.raster@data@values))])))

# Matérn spatial field raster plot
pred.field.raster <- raster::mask(model.result$prediction_model$mean_prediction$field, model.result$dis_data$polygon_shapefile, inverse = FALSE)
plot(pred.field.raster, main = 'Predicted Mean Matérn Spatial Random Field Effects', col = brewer.pal(11, "PuOr"),
     zlim = c(-max(abs(pred.field.raster@data@values)[!is.na(abs(pred.field.raster@data@values))]),
              max(abs(pred.field.raster@data@values))[!is.na(abs(pred.field.raster@data@values))]))

# plot histogram of posterior parameter draws
ch <- Matrix::Cholesky(model.result$fitted_model$sd_out$jointPrecision)
par_draws <- sparseMVN::rmvn.sparse(100, model.result$fitted_model$obj$env$last.par.best, ch, prec = TRUE) %>%
  as.data.frame()
colnames(par_draws) <- names(model.result$fitted_model$obj$env$last.par.best)

# plot the fitted posterior parameters results
posteriors <- as.data.frame(summary(model.result$fitted_model$sd_out, select = 'fixed'))
posteriors <- dplyr::mutate(posteriors, name = c('Intercept','DOS','DLI','DMI','Log\nSigma','Log\nRho','Logit\nLambda','Log\nPrec'))
names(posteriors) <- c('mean', 'sd', 'parameter')

fixedeffects <- ggplot() +
  geom_errorbar(posteriors, mapping = aes(x = parameter, ymin = mean - sd, ymax = mean + sd), width = 0.4, color = "blue", size = 0.4) +
  geom_point(posteriors, mapping = aes(x = parameter, y = mean)) +
  ggtitle("Fixed effects") +
  theme_bw() +
  theme(plot.title = element_text(size=20, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size =15),
        axis.text.x = element_text(size =15))

fixedeffects

# plot/summary the in-sample performance
report<-model.result$fitted_model$obj$report()
observed_data = report$polygon_response_data
predicted_data = report$reportprediction_cases

insample.data.df <- data.frame(obs = observed_data, pred = predicted_data)

obspred <- ggplot(insample.data.df, aes(x = obs, y = pred)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = 'blue') +
  ggtitle('In sample performance:\nTract Population') +
  theme_bw() +
  theme(plot.title = element_text(size=20, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size =15),
        axis.text.x = element_text(size=15))

obspred
insample.metrics <- dplyr::summarise(insample.data.df,
                            RMSE = sqrt(mean((pred - obs) ^ 2)),
                            MAE = mean(abs(pred - obs)),
                            pearson = cor(pred, obs, method = 'pearson'),
                            spearman = cor(pred, obs, method = 'spearman'),
                            log_pearson = cor(log1p(pred), log1p(obs), method = 'pearson'))

# plot parameter samples from posterior
# intercept
plot_norm_prior_post(data = par_draws$intercept,
                     mean = priors$priormean_intercept,
                     sd = priors$priorsd_intercept,
                     par = 'Intercept',
                     par_val = model.result$fitted_model$opt$par[[1]],
                     limits = c(-3.5 * priors$priorsd_intercept + priors$priormean_intercept,
                                3.5 * priors$priorsd_intercept + priors$priormean_intercept), 
                     bins = 60)

# DOS
plot_norm_prior_post(data = par_draws[,2],
                     mean = priors$priormean_slope[[1]],
                     sd = priors$priorsd_slope,
                     par = "DOS",
                     par_val = model.result$fitted_model$opt$par[[2]],
                     limits = c(-3.5 * priors$priorsd_slope + priors$priormean_slope[[1]],
                                3.5 * priors$priorsd_slope + priors$priormean_slope[[1]]),
                     bins = 80)

# DLI
plot_norm_prior_post(data = par_draws[,3],
                     mean = priors$priormean_slope[[2]],
                     sd = priors$priorsd_slope,
                     par = "DLI",
                     par_val = model.result$fitted_model$opt$par[[3]],
                     limits = c(-3.5 * priors$priorsd_slope + priors$priormean_slope[[2]],
                                3.5 * priors$priorsd_slope + priors$priormean_slope[[2]]),
                     bins = 80)


# DMI
plot_norm_prior_post(data = par_draws[,4],
                     mean = priors$priormean_slope[[3]],
                     sd = priors$priorsd_slope,
                     par = "DMI",
                     par_val = model.result$fitted_model$opt$par[[4]],
                     limits = c(-3.5 * priors$priorsd_slope + priors$priormean_slope[[3]],
                                3.5 * priors$priorsd_slope + priors$priormean_slope[[3]]),
                     bins = 80)

# Log Sigma
log.sigma.samp <- par_draws$log_sigma
hist(log.sigma.samp, xlab="Estimated Log-Sigma",  main="Histogram of Estimated Log-Sigma Value", col = 'blue')
abline(v = model.result$fitted_model$opt$par[[5]], col="red", lwd=3, lty=2)

# Log Rho
log.rho.samp <- par_draws$log_rho
hist(log.rho.samp, xlab="Estimated Log-Rho",  main="Histogram of Estimated Log-Rho Value", col = 'blue')
abline(v = model.result$fitted_model$opt$par[[6]], col="red", lwd=3, lty=2)

# Lambda
plot_norm_prior_post(data = par_draws$logit_car_rho,
                     mean = priors$prior_logit_car_rho_mean,
                     sd = priors$prior_logit_car_rho_sd,
                     par = "Logit-Lambda",
                     par_val = model.result$fitted_model$opt$par[[7]],
                     limits = c(-3.5 * priors$prior_logit_car_rho_sd + priors$prior_logit_car_rho_mean,
                                3.5 * priors$prior_logit_car_rho_sd + priors$prior_logit_car_rho_mean),
                     bins = 60)

# Precision
plot_gamma_prior_post(data = exp(par_draws$log_car_prec),
                     shape = priors$prior_car_gamma_shape,
                     scale = priors$prior_car_gamma_scale,
                     par = "Precision",
                     par_val = exp(model.result$fitted_model$opt$par[[8]]),
                     limits = c(0,8),
                     bins = 60)

# plot uncertainty prediction--------------------------------------------------------------------------------------------------
probs <- c((1 - 0.95) / 2, 1 - (1 - 0.95) / 2)

# The prediction CI is in a temp file cannot be read back in
# resample it here
predictions_ci <- raster::calc(model.result$prediction_model$uncertainty_prediction$realisations, function(x) stats::quantile(x, probs = probs, na.rm = TRUE))
names(predictions_ci) <- c('lower CI', 'upper CI')

# plot the upper and lower CI
spplot(predictions_ci, col.regions = rev(terrain.colors(25)), main = 'Uncertainty Predictions')

# plot 4 realizations of spatially disaggregated population by HPSDRM
reaz1 <- model.result$prediction_model$uncertainty_prediction$realisations$layer.1
reaz2 <- model.result$prediction_model$uncertainty_prediction$realisations$layer.2
reaz3 <- model.result$prediction_model$uncertainty_prediction$realisations$layer.3
reaz4 <- model.result$prediction_model$uncertainty_prediction$realisations$layer.4
plot(reaz1, main = 'Realization 1')
plot(reaz2, main = 'Realization 2')
plot(reaz3, main = 'Realization 3')
plot(reaz4, main = 'Realization 4')

#---------Aggregate the spatially interpolated grids population to 2020 Decennial blocks-------------------
# aggregate the spatially interpolated grids population to 2020 Decennial blocks
hpsdrm_dasy_pop_block <- exactextractr::exact_extract(hpsdrm_dasy_pop_raster, pop.block, fun = 'sum')

#-------Evaluate the accuracy/performance of the HPSDRM spatially interpolated grids population-----------------
r_square <- function(sample, predict){
   rss <- sum((sample - predict)^2)
   tss <- sum((sample - mean(sample))^2)
   r2 <- 1 - rss/tss
}

# R square
HPSDRM.r2 <- r_square(pop.block$value, hpsdrm_dasy_pop_block)

# RMSE
HPSDRM.rmse <- Metrics::rmse(pop.block$value, hpsdrm_dasy_pop_block)

# MAE
HPSDRM.mae <- Metrics::mae(pop.block$value, hpsdrm_dasy_pop_block)

# dataframe of true and pred block population
hpsdrm.block.df <- pop.block %>%
  rename(true = value) %>%
  cbind(pred = hpsdrm_dasy_pop_block)

# plot out-of-sample-block performance
hpsdrm.block.obspred <- ggplot(hpsdrm.block.df, aes(x = true, y = pred)) +
  geom_point(alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = 'blue') +
  ggtitle('HPSDRM out sample performance:\nBlock Population') +
  theme_bw() +
  theme(plot.title = element_text(size=20, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size =15),
        axis.text.x = element_text(size=15))

hpsdrm.block.obspred

hpsdrm.outsample.metrics <- dplyr::summarise(hpsdrm.block.df %>% st_drop_geometry(),
                                     RMSE = sqrt(mean((pred - true) ^ 2)),
                                     MAE = mean(abs(pred - true)),
                                     pearson = cor(pred, true, method = 'pearson'),
                                     spearman = cor(pred, true, method = 'spearman'),
                                     log_pearson = cor(log1p(pred), log1p(true), method = 'pearson'))

HPSDRM.res <- list(model.result = model.result,
                   HPSDRM.r2 = HPSDRM.r2,
                   HPSDRM.rmse = HPSDRM.rmse,
                   HPSDRM.mae = HPSDRM.mae,
                   hpsdrm.block.df = hpsdrm.block.df
)

saveRDS(HPSDRM.res, file = (paste0(proj_root, "/data/gen/HPSDRM.res.rds")))