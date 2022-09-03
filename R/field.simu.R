################################################################################# 
## R code for "Incorporating spatial autocorrelation in dasymetric             ## 
## analysis: A hierarchical poisson spatial disaggregation regression model"   ##
## Authors: Bowen He, Jonathan M. Gilligan, Janey V. Camp                      ## 
#################################################################################

#####################################################
########### FIELD PSDRM SIMULATION STUDY ############
#####################################################

# Poisson Spatial Disaggregation Regression Model with the only Matérn FIELD 
# covariance gaussian latent field model as the latent component.
# To model the population count data, the Poisson model is applied. 
# log(lambda) = intercept + covariate1 + covariate2 + GP (Gaussian random process/field model)
# GP ~ Matérn spatial random field
# population ~ Poisson(lambda)
# In this case, the linear predictor is considered as the logarithm of lambda value (mean population counts)
# and could be used to simulate the population counts in the latent finer grid.

# This 100-simu experiment finds that even the single level PSDRM spatial model exceeds
# the traditional Areal Weighting model in terms of interpolation accuracy. 

# Load required library
library('TMB')
library('disaggregation')
library('MASS')
library('INLA')
library('dplyr')
library('tidyverse')
library('sp')
library('sf')
library('rgdal')
library('spdep')
library('maptools')
library('Matrix')
library('LaplacesDemon')
library('ggplot2')
library('raster')
library('gtools')
library('parallel')
library('rprojroot')
library('exactextractr')
library('spatialreg')
library('fasterize')
library('RColorBrewer')
proj_root <- find_root(is_rstudio_project)

#compile(paste0(proj_root,"/src/dasymap.cpp"))
dyn.load(dynlib(paste0(proj_root,"/src/dasymap")))

source(paste0(proj_root, "/R/dasymap.fit.func.R"))
source(paste0(proj_root, "/R/dasymap.predict.func.R"))
source(paste0(proj_root, "/R/dasymap.scale.func.R"))

# Set up the simulation study parameter space----------------------------------
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

grid.x <- n
grid.y <- n

grid <- expand.grid(x = seq(0, 1, length.out = grid.x), y = seq(0, 1, length.out = grid.y))

# dmat <- dist(grid)
# mcor <- as.matrix(2^(1-nu)*(kappa*dmat)^nu*besselK(dmat*kappa,nu)/gamma(nu))
# 
# diag(mcor) <- 1
# mcov <- sigma2x*mcor
# 
# L <- chol(mcov)

L <- readRDS(file.path(proj_root, "data/gen/L.rds"))

# obtain the areal value
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
               prior_sigma_max = 1)

i <- 0
field.simu.res <- list()
n_realizations <- 100

# 100-realizations experiment----------------------------------------------
while (i < n_realizations) {
  
print(i)
matern_field <- data.frame(x=grid$x, y = grid$y, value = drop(rnorm(n^2)%*%L))
matern_field.raster <- rasterFromXYZ(matern_field, crs = "+proj=longlat +datum=WGS84 +no_defs")
  
cov1.df <- data.frame(x=grid$x, y = grid$y, value = rnorm(n = n^2))
cov1.raster <- rasterFromXYZ(cov1.df, crs = "+proj=longlat +datum=WGS84 +no_defs")
cov2.df <- data.frame(x=grid$x, y = grid$y, value = rnorm(n = n^2))
cov2.raster <- rasterFromXYZ(cov2.df, crs = "+proj=longlat +datum=WGS84 +no_defs")

# synthetic linear predictor
cov.raster <- beta0 + beta1*cov1.raster + beta2*cov2.raster
eta.raster <- beta0 + beta1*cov1.raster + beta2*cov2.raster + matern_field.raster
simu.true.raster <- eta.raster 
simu.true.raster@data@values <- rpois(length(exp(eta.raster@data@values)), lambda = exp(eta.raster@data@values))

polygon.response.raster <- raster::aggregate(simu.true.raster, fact = fct, fun = 'sum') 
polygon.response <- polygon.response.raster %>%
  as('SpatialPolygonsDataFrame') %>%
  st_as_sf() %>%
  mutate(ID = row_number()) %>%
  rename('value' = 'layer')

shapes <- as(polygon.response %>% dplyr::select(ID, value), 'Spatial')
covariate_stack <- stack(cov1.raster, cov2.raster)
names(covariate_stack) <- c('cov1', 'cov2')

# prepare data
dis_data_simu <- disaggregation::prepare_data(polygon_shapefile = shapes,
                              covariate_rasters = covariate_stack,
                              mesh.args = list(max.edge = c(0.05, 0.1),
                                               cut = 0.01,
                                               offset = c(0.03, 0.05)),
                              id_var = 'ID',
                              response_var = 'value',
                              na.action = TRUE,
                              ncores = 10)
# fit the spatial model
fitted_model_simu <- dasymap_fit_func(data = dis_data_simu,
                                      priors = priors,
                                      family = 'poisson',
                                      link = 'log',
                                      iterations = 1000,
                                      pcfield = FALSE,
                                      gausfield = TRUE,
                                      car = FALSE,
                                      carmatrix = mymatrix,
                                      iid = FALSE,
                                      silent = TRUE)

# Model prediction section
pred_model_simu <- predict.dasymap_model(fitted_model_simu, predict_iid = FALSE, predict_car = FALSE)

simu.res.output <- list(simu.true.raster = simu.true.raster,
                        cov.raster = cov.raster,
                        matern_field.raster = matern_field.raster,
                        dis_data_simu = dis_data_simu,
                        fitted_model_simu = fitted_model_simu,
                        pred_model_simu = pred_model_simu)
i = i + 1
field.simu.res[[i]] <- simu.res.output
}
saveRDS(field.simu.res, file = (paste0(proj_root, "/data/gen/field.simu.res.gaus.sd1.rds")))
#-------------------------Experiment Results Evaluation--------------------
#Evaluation field simulation results
field.eva.output <- list()
for (i in 1:length(field.simu.res)){
 pred.raster <- field.simu.res[[i]]$pred_model_simu$mean_prediction$prediction

# best/optimal parameter
intercept <- field.simu.res[[i]]$fitted_model_simu$opt$par[['intercept']]
slope <- field.simu.res[[i]]$fitted_model_simu$opt$par[['slope']]
slope.1 <- field.simu.res[[i]]$fitted_model_simu$opt$par[[3]]
log_sigma_best <- field.simu.res[[i]]$fitted_model_simu$opt$par[['log_sigma']]
log_rho_best <- field.simu.res[[i]]$fitted_model_simu$opt$par[['log_rho']]
sigma_best <- exp(log_sigma_best)
kappa_best <- sqrt(8)/exp(log_rho_best)

# parameter sd
sd.df <- as.data.frame(summary(field.simu.res[[i]]$fitted_model_simu$sd_out))
intercept.sd <-sd.df['intercept',]$`Std. Error`
slope.sd <-sd.df['slope',]$`Std. Error`
slope.1.sd <-sd.df['slope.1',]$`Std. Error`
log_sigma_sd <- sd.df['log_sigma',]$`Std. Error`
log_rho_sd <- sd.df['log_rho',]$`Std. Error`

# Bayesian posterior information
ch <- Matrix::Cholesky(field.simu.res[[i]]$fitted_model_simu$sd_out$jointPrecision)
par_draws <- as.data.frame(sparseMVN::rmvn.sparse(1000, field.simu.res[[i]]$fitted_model_simu$obj$env$last.par.best, ch, prec = TRUE))
colnames(par_draws) <- names(field.simu.res[[i]]$fitted_model_simu$obj$env$last.par.best)

# true parameter quantile
intercept_q <- ecdf(par_draws$intercept)(beta0)
slope_q <- ecdf(par_draws[,2])(beta1)
slope1_q <- ecdf(par_draws[,3])(beta2)
log_sigma_q <- ecdf(par_draws$log_sigma)(log_sigma)
log_rho_q <- ecdf(par_draws$log_rho)(log_rho)


output <- list(pred.raster = pred.raster,
               intercept = intercept,
               intercept.sd = intercept.sd,
               intercept_q = intercept_q,
               slope = slope,
               slope.sd = slope.sd,
               slope_q = slope_q,
               slope.1 = slope.1,
               slope.1.sd = slope.1.sd,
               slope1_q = slope1_q,
               log_sigma_best = log_sigma_best,
               log_sigma_sd = log_sigma_sd,
               log_sigma_q = log_sigma_q,
               log_rho_best = log_rho_best,
               log_rho_sd = log_rho_sd,
               log_rho_q = log_rho_q,
               par_draws = par_draws)

field.eva.output[[i]] <- output
}

# Fixed effects parameters assessment----------
# intercept--------------------------------------
intercept_vtr <- c()
for (i in 1:length(field.eva.output)) {
  intercept_vtr[i]  <- field.eva.output[[i]]$intercept
}
hist(intercept_vtr, xlab="Estimated Intercept",  main="Histogram of Estimated Intercept Value", col = 'blue')
abline(v = beta0, col="red", lwd=3, lty=2)

# Find the corresponding quantile of the true intercept
intercept_q_vtr <- c()
for (i in 1:length(field.eva.output)) {
  intercept_q_vtr[i]  <- field.eva.output[[i]]$intercept_q
}

hist(intercept_q_vtr, xlab="True Intercept Quantile", main="Histogram of True Intercept Quantiles")

# intercept sd
intercept_sd_vtr <- c()
for (i in 1:length(field.eva.output)) {
  intercept_sd_vtr[i] <- field.eva.output[[i]]$intercept.sd
}

intercept_df <- data.frame(mean = intercept_vtr, sd = intercept_sd_vtr, quantile = intercept_q_vtr)

ggplot(intercept_df, aes(mean, sd)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = beta0, linetype="dashed",
             color = "red", size=1.0)+
  labs(title="Estimated Intercept Mean vs. SD", x = "Estimated Mean", y = "Estimated SD")  +
  theme_bw() +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=10, face = 'bold'),
        legend.text = element_text(size=9),
        axis.title.x = element_text(size=13, face="bold" ),
        axis.title.y = element_text(size=13, face="bold"),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))

ggplot(intercept_df, aes(mean, quantile)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = beta0, linetype="dashed",
             color = "red", size=1.0)+
  labs(title="Estimated Intercept Mean vs. True Quantile", x = "Estimated Mean", y = "Estimated True Quantile")  +
  theme_bw() +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=10, face = 'bold'),
        legend.text = element_text(size=9),
        axis.title.x = element_text(size=13, face="bold" ),
        axis.title.y = element_text(size=13, face="bold"),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))


# slope-----------------------------------------------------------------------------------------------------
slope_vtr <- c()
for (i in 1:length(field.eva.output)) {
  slope_vtr[i]  <- field.eva.output[[i]]$slope
}
hist(slope_vtr, xlab="Estimated Slope.1",  main="Histogram of Estimated Slope.1 Value", col = 'blue')
abline(v = beta1, col="red", lwd=3, lty=2)

# Find the corresponding quantile of the true intercept
slope_q_vtr <- c()
for (i in 1:length(field.eva.output)) {
  slope_q_vtr[i]  <- field.eva.output[[i]]$slope_q
}

hist(slope_q_vtr, xlab="True Slope.1 Quantile", main="Histogram of True Slope.1 Quantiles")

# slope sd
slope_sd_vtr <- c()
for (i in 1:length(field.eva.output)) {
  slope_sd_vtr[i] <- field.eva.output[[i]]$slope.sd
}

slope_df <- data.frame(mean = slope_vtr, sd = slope_sd_vtr, quantile = slope_q_vtr)

ggplot(slope_df, aes(mean, sd)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = beta1, linetype="dashed",
             color = "red", size=1.0)+
  labs(title="Estimated Slope.1 Mean vs. SD", x = "Estimated Mean", y = "Estimated SD")  +
  theme_bw() +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=10, face = 'bold'),
        legend.text = element_text(size=9),
        axis.title.x = element_text(size=13, face="bold" ),
        axis.title.y = element_text(size=13, face="bold"),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))

ggplot(slope_df, aes(mean, quantile)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = beta1, linetype="dashed",
             color = "red", size=1.0)+
  labs(title="Estimated Slope.1 Mean vs. True Quantile", x = "Estimated Mean", y = "Estimated True Quantile")  +
  theme_bw() +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=10, face = 'bold'),
        legend.text = element_text(size=9),
        axis.title.x = element_text(size=13, face="bold" ),
        axis.title.y = element_text(size=13, face="bold"),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))

# slope1-----------------------------------------------------------------------------------------------------
slope1_vtr <- c()
for (i in 1:length(field.eva.output)) {
  slope1_vtr[i]  <- field.eva.output[[i]]$slope.1
}
hist(slope1_vtr, xlab="Estimated Slope.2",  main="Histogram of Estimated Slope.2 Value", col = 'blue')
abline(v = beta2, col="red", lwd=3, lty=2)

# Find the corresponding quantile of the true intercept
slope1_q_vtr <- c()
for (i in 1:length(field.eva.output)) {
  slope1_q_vtr[i]  <- field.eva.output[[i]]$slope1_q
}

hist(slope1_q_vtr, xlab="True Slope.2 Quantile", main="Histogram of True Slope.2 Quantiles")

# slope sd
slope1_sd_vtr <- c()
for (i in 1:length(field.eva.output)) {
  slope1_sd_vtr[i] <- field.eva.output[[i]]$slope.1.sd
}

slope1_df <- data.frame(mean = slope1_vtr, sd = slope1_sd_vtr, quantile = slope1_q_vtr)

ggplot(slope1_df, aes(mean, sd)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = beta2, linetype="dashed",
             color = "red", size=1.0)+
  labs(title="Estimated Slope.2 Mean vs. SD", x = "Estimated Mean", y = "Estimated SD")  +
  theme_bw() +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=10, face = 'bold'),
        legend.text = element_text(size=9),
        axis.title.x = element_text(size=13, face="bold" ),
        axis.title.y = element_text(size=13, face="bold"),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))

ggplot(slope1_df, aes(mean, quantile)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = beta2, linetype="dashed",
             color = "red", size=1.0)+
  labs(title="Estimated Slope.2 Mean vs. True Quantile", x = "Estimated Mean", y = "Estimated True Quantile")  +
  theme_bw() +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=10, face = 'bold'),
        legend.text = element_text(size=9),
        axis.title.x = element_text(size=13, face="bold" ),
        axis.title.y = element_text(size=13, face="bold"),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))

# Matérn field (hyper)parameters estimation assessment----------------------------------------------------------------
# logsigma-----------------------
logsigma_vtr <- c()
for (i in 1:length(field.eva.output)) {
  logsigma_vtr[i]  <- field.eva.output[[i]]$log_sigma_best
}
hist(logsigma_vtr, xlab="Estimated Log-Sigma",  main="Histogram of Estimated Log-Sigma Value", col = 'blue')
abline(v = log_sigma, col="red", lwd=3, lty=2)

# Find the corresponding quantile of the true intercept
logsigma_q_vtr <- c()
for (i in 1:length(field.eva.output)) {
  logsigma_q_vtr[i]  <- field.eva.output[[i]]$log_sigma_q
}

hist(logsigma_q_vtr, xlab="True Log-Sigma Quantile", main="Histogram of True Log-Sigma Quantiles")

# slope sd
logsigma_sd_vtr <- c()
for (i in 1:length(field.eva.output)) {
  logsigma_sd_vtr[i] <- field.eva.output[[i]]$log_sigma_sd
}

logsigma_df <- data.frame(mean = logsigma_vtr, sd = logsigma_sd_vtr, quantile = logsigma_q_vtr)

ggplot(logsigma_df, aes(mean, sd)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = log_sigma, linetype="dashed",
             color = "red", size=1.0)+
  labs(title="Estimated Log-Sigma Mean vs. SD", x = "Estimated Mean", y = "Estimated SD")  +
  theme_bw() +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=10, face = 'bold'),
        legend.text = element_text(size=9),
        axis.title.x = element_text(size=13, face="bold" ),
        axis.title.y = element_text(size=13, face="bold"),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))

ggplot(logsigma_df, aes(mean, quantile)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = log_sigma, linetype="dashed",
             color = "red", size=1.0)+
  labs(title="Estimated Log-Sigma Mean vs. True Quantile", x = "Estimated Mean", y = "Estimated True Quantile")  +
  theme_bw() +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=10, face = 'bold'),
        legend.text = element_text(size=9),
        axis.title.x = element_text(size=13, face="bold" ),
        axis.title.y = element_text(size=13, face="bold"),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))

# logrho----------------------------------------------------
logrho_vtr <- c()
for (i in 1:length(field.eva.output)) {
  logrho_vtr[i]  <- field.eva.output[[i]]$log_rho_best
}
hist(logrho_vtr, xlab="Estimated Log-Rho",  main="Histogram of Estimated Log-Rho Value", col = 'blue')
abline(v = log_rho, col="red", lwd=3, lty=2)

# Find the corresponding quantile of the true intercept
logrho_q_vtr <- c()
for (i in 1:length(field.eva.output)) {
  logrho_q_vtr[i]  <- field.eva.output[[i]]$log_rho_q
}

hist(logrho_q_vtr, xlab="True Log-Rho Quantile", main="Histogram of True Log-Rho Quantiles")

# slope sd
logrho_sd_vtr <- c()
for (i in 1:length(field.eva.output)) {
  logrho_sd_vtr[i] <- field.eva.output[[i]]$log_rho_sd
}

logrho_df <- data.frame(mean = logrho_vtr, sd = logrho_sd_vtr, quantile = logrho_q_vtr)

ggplot(logrho_df, aes(mean, sd)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = log_rho, linetype="dashed",
             color = "red", size=1.0)+
  labs(title="Estimated Log-Rho Mean vs. SD", x = "Estimated Mean", y = "Estimated SD")  +
  theme_bw() +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=10, face = 'bold'),
        legend.text = element_text(size=9),
        axis.title.x = element_text(size=13, face="bold" ),
        axis.title.y = element_text(size=13, face="bold"),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))

ggplot(logrho_df, aes(mean, quantile)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = log_rho, linetype="dashed",
             color = "red", size=1.0)+
  labs(title="Estimated Log-Rho Mean vs. True Quantile", x = "Estimated Mean", y = "Estimated True Quantile")  +
  theme_bw() +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=10, face = 'bold'),
        legend.text = element_text(size=9),
        axis.title.x = element_text(size=13, face="bold" ),
        axis.title.y = element_text(size=13, face="bold"),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))

# Simulation Visualization Example-----------------------------------------------------
# vis1
# underlying true simulation grids value
simu.true1 <- field.simu.res[[3]]$simu.true.raster
plot(simu.true1, main = 'Underlying Simulated True Grids Value',
     zlim = c(0, max(simu.true1@data@values)))

# Simulated Matérn spatial field
simu.field.raster1 <- field.simu.res[[1]]$matern_field.raster
plot(simu.field.raster1, main = 'Simulated Matérn Spatial Random Field', col = brewer.pal(11, "PuOr"),
     zlim = c(-max(abs(simu.field.raster1@data@values)), max(abs(simu.field.raster1@data@values))))

# Spatial Disaggregated Grids Value
pred.raster1 <- field.simu.res[[3]]$pred_model_simu$mean_prediction$prediction
field.simu.vis1<-scale_func_simu(n=5, simulation.true.raster = simu.true1, pred.raster = pred.raster1)
plot(field.simu.vis1$pred.scaled.raster, main = "Spatially Disaggregated Grids Value",
     zlim = c(0, max(simu.true1@data@values)))

# Predicted Matérn spatial field
pred.field.raster1 <- field.simu.res[[1]]$pred_model_simu$mean_prediction$field
plot(pred.field.raster1, main = 'Predicted Matérn Spatial Random Field', col = brewer.pal(11, "PuOr"),
     zlim = c(-max(abs(simu.field.raster1@data@values)), max(abs(simu.field.raster1@data@values))))

# Polygon Response raster plot
polygon.ras1<-rasterize(field.simu.res[[3]]$dis_data_simu$polygon_shapefile, field.simu.res[[3]]$simu.true.raster, field = 'value')
plot(polygon.ras1, main='Polygon Response Data',
     zlim = c(0, max(polygon.ras1@data@values)))

# Covariances rasters plot
# cov1.ras1 <- plot(cov1.raster, main = 'cov1', col = brewer.pal(11, "PuOr"), zlim = c(min(cov2.raster@data@values), max(cov2.raster@data@values)))
# cov2.ras1 <- plot(cov2.raster, main = 'cov2', col = brewer.pal(11, "PuOr"), zlim = c(min(cov2.raster@data@values), max(cov2.raster@data@values)))
# Model residuals error evaluation------------------------------------
field.simu.res.df <- field.simu.res[[1]]$simu.true.raster %>%
  rasterToPoints() %>%
  as.data.frame() %>%
  dplyr::select('x', 'y')

field.simu.res.vct <- list()
i <- 0
for (i in 1:length(field.simu.res)){
  field.simu <- field.simu.res[[i]]$simu.true.raster %>%
    raster::rasterToPoints() %>%
    as.data.frame()

  field.simu.pred <- scale_func_simu(n=10, simulation.true.raster = field.simu.res[[i]]$simu.true.raster,
                               pred.raster = field.simu.res[[i]]$pred_model_simu$mean_prediction$prediction) %$%
    pred.scaled.raster %>%
    raster::rasterToPoints() %>%
    as.data.frame()

  field.simu.res.df <- cbind(field.simu.res.df,
                            tibble(!!paste0('simu',i) := field.simu$layer),
                            tibble(!!paste0('pred',i) := field.simu.pred$layer),
                            tibble(!!paste0('resi.error',i) := field.simu.pred$layer - field.simu$layer))
  field.simu.res.vct[[i]] <- list(pred = field.simu.pred$layer,
                                 simu = field.simu$layer)
}

# Mean Raw Count (MRC) Evaluation
field.simu.raw_count_diff <- field.simu.res.df %>%
  as_tibble %>%
  dplyr::select(x, y, matches("resi\\.error[0-9]+")) %>%
  rowwise() %>%
  transmute(x, y, mean_raw_count_diff = mean(c_across(-c(x,y))))

field.simu.raw.count.diff.sf <- st_as_sf(field.simu.raw_count_diff, coords = c("x", "y"))
field.simu.raw.count.diff.raster <- rasterize(field.simu.raw.count.diff.sf, field.simu.res[[1]]$simu.true.raster, field = 'mean_raw_count_diff', fun='last')

plot(field.simu.raw.count.diff.raster,
     main = "Mean Prediction Error (Raw Count Difference)", col = brewer.pal(11, "RdBu"),
     zlim = c(-max(abs(field.simu.raw.count.diff.raster@data@values)), max(abs(field.simu.raw.count.diff.raster@data@values))))

# RMSE Evaluation
field.simu.rmse <- c()
for (i in 1:length(field.simu.res.vct)) {
  field.simu.rmse[i] <- Metrics::rmse(field.simu.res.vct[[i]]$simu,  field.simu.res.vct[[i]]$pred)
}

hist(field.simu.rmse, main = 'Histogram of RMSE', xlab = 'RMSE')

# RMSE Grid Evaluation
field.simu.rmse.grid <- field.simu.res.df %>%
  as_tibble %>%
  dplyr::select(x, y, matches("resi\\.error[0-9]+")) %>%
  rowwise() %>%
  transmute(x, y, rmse = sqrt(mean(c_across(-c(x,y))^2)))

field.simu.rmse.grid.sf <- st_as_sf(field.simu.rmse.grid, coords = c("x", "y"))
field.simu.rmse.grid.raster <- rasterize(field.simu.rmse.grid.sf, field.simu.res[[1]]$simu.true.raster, field = 'rmse', fun='last')

plot(field.simu.rmse.grid.raster,
     main = "Spatial RMSE")

# R-square Evaluation
# R2 comparison
r_square_calc <- function(sample, predict){
  rss <- sum((sample - predict)^2)
  tss <- sum((sample - mean(sample))^2)
  r2 <- 1 - rss/tss
}

field.simu.r.square <- c()
for (i in 1:length(field.simu.res.vct)) {
  field.simu.r.square[i] <-r_square_calc(field.simu.res.vct[[i]]$simu, field.simu.res.vct[[i]]$pred)
}

hist(field.simu.r.square, main = expression(bold("Histogram of R"^"2")), xlab = expression("R"^2), breaks = 10)

# (Example)Realization Interpolation accuracy Comparison with AWM---------------------------------------
## Absolute value difference between Areal Weighting model and HPSRM
field.simu.vis1$sf$abs.diff <- abs(field.simu.vis1$sf$aw.resi.error) - abs(field.simu.vis1$sf$scaled.resi.error)
field.simu.vis1$sf$sqrd.diff <- (field.simu.vis1$sf$aw.resi.error)^2 - (field.simu.vis1$sf$scaled.resi.error)^2

field.simu.vis1.abs.diff <- fasterize(field.simu.vis1$sf, field.simu.vis1$pred.scaled.raster, field = 'abs.diff', fun = 'last')
# %>%
#   as.data.frame(xy=TRUE) %>%
#   rename('abs.diff' = 'layer')

plot(field.simu.vis1.abs.diff, main = "Absolute Difference of Residual Error \nbetween the AWM and the HPSRM")

# ggplot() +
#   geom_raster(data = field.simu.vis1.abs.diff, aes(x=x, y=y, fill = abs.diff)) +
#   scale_fill_gradient2(limits=c(-max(abs(field.simu.vis1.abs.diff$abs.diff)), max(abs(field.simu.vis1.abs.diff$abs.diff)))) +
#   labs(title="Absolute Difference of Residual Error of AWM compared to HPSRM")  +
#   theme_bw() +
#   theme(plot.title = element_text(size = 30, hjust = 0.5, face = 'bold'),
#         legend.title = element_text(size = 25, face = 'bold'),
#         legend.text = element_text(size = 25),
#         axis.title.x = element_text(size = 25, face="bold" ),
#         axis.title.y = element_text(size = 25, face="bold"),
#         axis.text.y = element_text(size = 25),
#         axis.text.x = element_text(size = 25))

# Squared value difference between Areal Weighting model and HPSRM
field.simu.vis1.sqrd.diff <- fasterize(field.simu.vis1$sf, field.simu.vis1$pred.scaled.raster, field = 'sqrd.diff', fun = 'last')
# %>%
#   as.data.frame(xy=TRUE) %>%
#   rename('squared.diff' = 'layer')

plot(field.simu.vis1.sqrd.diff, main = "Squared Difference of Residual Error \nbetween the AWM and the HPSRM")

ggplot() +
  geom_raster(data = field.simu.vis1.sqrd.diff, aes(x=x, y=y,fill = squared.diff)) +
  scale_fill_gradient2(limits=c(-max(abs(field.simu.vis1.sqrd.diff$squared.diff)), max(abs(field.simu.vis1.sqrd.diff$squared.diff)))) +
  labs(title="Squared Difference of Residual Error of AWM compared to HPSRM")  +
  theme_bw() +
  theme(plot.title = element_text(size = 25, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size = 20, face = 'bold'),
        legend.text = element_text(size = 20),
        axis.title.x = element_text(size = 20, face="bold"),
        axis.title.y = element_text(size = 20, face="bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20))