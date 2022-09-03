################################################################################# 
## R code for "Incorporating spatial autocorrelation in dasymetric             ## 
## analysis: A hierarchical poisson spatial disaggregation regression model"   ##
## Authors: Bowen He, Jonathan M. Gilligan, Janey V. Camp                      ## 
#################################################################################

######################################################
########## FIELD-CAR HPSDRM SIMULATION STUDY #########
######################################################

# Conduct (hyper)priors' parameters sensitivity analysis on the 
# integrated field-car hpsdrm simulation model using the same 100
# simulated realizations.

# load required library
library(TMB)
library(disaggregation)
library(tidyverse)
library(raster)
library(sf)
library(sp)
library(INLA)
library(rprojroot)
library(scales)
library(vctrs)
library(tidycensus)
library(fasterize)
library(LaplacesDemon)
library(spdep)
library(spatialreg)
proj_root <- find_root(is_rstudio_project)

#compile(paste0(proj_root,"/src/dasymap.cpp"))
dyn.load(dynlib(paste0(proj_root,"/src/dasymap")))

source(paste0(proj_root, "/R/dasymap.fit.func.R"))
source(paste0(proj_root, "/R/dasymap.predict.func.R"))
source(paste0(proj_root, "/R/dasymap.scale.func.R"))


field.car.simu.res <- readRDS(file.path(proj_root, "data/gen/sensitivity_analysis/field.car.simu.res.sd1.rds"))

n <- 150
beta0 <- 2
beta1 <- 1
beta2 <- 0.5
sigma2x <- 1
kappa <- 7
rho <- sqrt(8)/kappa
nu <- 1
log_sigma <- log(sqrt(sigma2x))
log_rho <- log(rho)

# pixel aggregation factor
fct <- 5

# car parameter
car_rho <- 0.7
car_var <- 3
car_prec <- 1/car_var
logit_car_rho <- logit(car_rho)
log_car_prec <- log(car_prec)

#prior_gamma_scale <- (-car_prec + sqrt(car_prec^2 + 4))/2
#prior_gamma_shape <- 1 + car_prec/prior_gamma_scale

area.raster <- raster(ncol = n/fct, nrow = n/fct, xmn=0, xmx=1, ymn=0, ymx=1)
area <- area.raster %>%
  as('SpatialPolygonsDataFrame') %>%
  st_as_sf() %>%
  mutate(ID = row_number())

area.nb <- poly2nb(area)
area.nb_B <- nb2listw(area.nb, style="B", zero.policy=TRUE)

sym.matrix <- as(area.nb_B, "symmetricMatrix")
one.matrix <- rep(1, nrow(area))
a <- Diagonal(x = as.numeric(sym.matrix%*%one.matrix))
mymatrix <- a - sym.matrix

# set up the priors using the selected parameters values in the experiment
priors <- list(priormean_intercept = beta0,
               priorsd_intercept = 1,
               priormean_slope = c(beta1, beta2),
               priorsd_slope = 1,
               prior_logrho_gaussian_mean = log_rho,
               prior_logrho_gaussian_sd = 1,
               prior_logsigma_gaussian_mean = log_sigma,
               prior_logsigma_gaussian_sd = 1,
               prior_rho_min = 1.3,
               prior_sigma_max = 1,
               prior_logit_car_rho_mean = logit_car_rho,
               prior_logit_car_rho_sd = 1,
               prior_car_gamma_shape = car_prec^2,
               prior_car_gamma_scale = 1/car_prec)

i <- 0
field.car.simu.res.sd1.gamma <- list()
n_realizations <- 100
while (i < n_realizations) {
  print(i)
  simu.true.raster <- field.car.simu.res[[i+1]]$simu.true.raster
  cov.raster <- field.car.simu.res[[i+1]]$cov.raster
  matern_field.raster <- field.car.simu.res[[i+1]]$matern_field.raster
  car_ras <- field.car.simu.res[[i+1]]$car_ras
  dis_data_simu <- field.car.simu.res[[i+1]]$dis_data_simu
  
  fitted_model_simu <- dasymap_fit_func(data = dis_data_simu,
                                        priors = priors,
                                        family = 'poisson',
                                        link = 'log',
                                        iterations = 1000,
                                        pcfield = FALSE,
                                        gausfield = TRUE,
                                        car = TRUE,
                                        carmatrix = mymatrix,
                                        iid = FALSE,
                                        silent = TRUE)
  
  pred_model_simu <- predict.dasymap_model(fitted_model_simu, predict_iid = FALSE, predict_car = TRUE)
  
  field.car.simu.output <- list(simu.true.raster = simu.true.raster,
                                cov.raster = cov.raster,
                                matern_field.raster = matern_field.raster,
                                car_ras =  car_ras,
                                dis_data_simu = dis_data_simu,
                                fitted_model_simu = fitted_model_simu,
                                pred_model_simu = pred_model_simu) 
  i = i + 1
  field.car.simu.res.sd1.gamma[[i]] <- field.car.simu.output
}
saveRDS(field.car.simu.res.sd1.gamma, file = (paste0(proj_root, "/data/gen/sensitivity_analysis/field.car.simu.res.sd1.gamma.rds")))
