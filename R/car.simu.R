################################################################################ 
## R code for "Incorporating spatial autocorrelation in dasymetric            ## 
## analysis: A hierarchical poisson spatial disaggregation regression model"  ##
## Authors: Bowen He, Jonathan M. Gilligan, Janey V. Camp                     ## 
################################################################################

################################################################################
#########################  CAR SIMULATION STUDY ################################
################################################################################

# Poisson Spatial Disaggregation Regression Model with the only 
# Leroux CAR prior as the latent components.
# To model the population count data, the Poisson model is applied. 
# log(lambda) = CAR (Leroux CAR prior)
# CAR ~ Leroux CAR prior
# population ~ Poisson(lambda)
# In this case, the linear predictor is considered as the logarithm of lambda value
# (mean population counts) and could be used to simulate the CAR area random effects.

# laod required libraries
library('TMB')
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
proj_root <- find_root(is_rstudio_project)

#compile(paste0(proj_root,"/src/car.cpp"))
dyn.load(dynlib(paste0(proj_root,"/src/car")))

source(paste0(proj_root,"/R/plot.prior.post.func.R"))
# set up parameters and hyperparameters---------------------------
rho <- 0.7
var <- 3
prec <- 1/var
logit_rho <- logit(rho)
log_prec <- log(prec)

initial_parameter <- list(
  logit_rho = logit_rho,
  log_prec = log_prec
)

# set up hyperpriors parameters
prior_logit_rho_mean <- logit_rho
prior_logit_rho_sd <- 1
prior_gamma_shape <- prec^2
prior_gamma_scale <- 1/prec

fct <- 5
raster <- raster(ncol=30, nrow=30, xmn=0, xmx=1, ymn=0, ymx=1)
raster.grid <- raster(ncol=30*fct, nrow=30*fct, xmn=0, xmx=1, ymn=0, ymx=1)
area <- raster %>%
  as('SpatialPolygonsDataFrame') 

area.nb <- poly2nb(area)
area.nb_B <- nb2listw(area.nb, style="B", zero.policy=TRUE)

sym.matrix <- as(area.nb_B, "symmetricMatrix")
one.matrix <- rep(1, nrow(area))
a <- Diagonal(x = as.numeric(sym.matrix%*%one.matrix))
mymatrix <- a - sym.matrix
Q <- prec*(rho*(mymatrix) + (1-rho)*diag(nrow(area)))

n_realizations <- 100  

# 100 realization experiments
i <- 0
car.simu.res <- list()
while (i < n_realizations) {

print(i)
simu_car <- inla.qsample(1, Q) 
lambda <- exp(simu_car)
simu <- rpois(n = length(lambda), lambda = lambda)

car_simu <- st_as_sf(area) %>%
  dplyr::select(-'layer') %>%
  mutate(simu = simu) 

input_data<-list(
  response_data = car_simu$simu,
  mymatrix = mymatrix,
  prior_logit_rho_mean = prior_logit_rho_mean,
  prior_logit_rho_sd = prior_logit_rho_sd,
  prior_gamma_shape = prior_gamma_shape,
  prior_gamma_scale = prior_gamma_scale
  )

parameters <- list(
  logit_rho = initial_parameter$logit_rho,
  log_prec = initial_parameter$log_prec,
  u = rep(0, nrow(car_simu))
)

obj <- TMB::MakeADFun(
  data = input_data, 
  parameters = parameters,
  random = c('u'),
  silent = TRUE,
  DLL = "car")

opt <- stats::nlminb(obj$par, obj$fn, control = list(iter.max = 1000, trace = 0))

# Calc uncertainty using the fixed hessian from above.
sd_out <- TMB::sdreport(obj, getJointPrecision = TRUE)

pred <- obj$env$last.par.best[3:902]
car_simu$pred <- exp(pred)

# best/optimal parameter
logit_rho_best <- opt$par[[1]]
log_prec_best <- opt$par[[2]]

# parameter sd
sd.df <- as.data.frame(summary(sd_out))
logit_rho_sd <- sd.df['logit_rho',]$`Std. Error`
log_prec_sd <- sd.df['log_prec',]$`Std. Error`
car_simu$sd <- sd.df[3:902,]

# Bayesian posterior information 
ch <- Matrix::Cholesky(sd_out$jointPrecision)
par_draws <- sparseMVN::rmvn.sparse(1000, obj$env$last.par.best, ch, prec = TRUE)

# true parameter quantile
logit_rho_q <- ecdf(par_draws[,1])(logit_rho)
log_prec_q <- ecdf(par_draws[,2])(log_prec)

output <- list(car_simu = car_simu,
               logit_rho_best = logit_rho_best,
               log_prec_best = log_prec_best,
               logit_rho_sd = logit_rho_sd,
               log_prec_sd = log_prec_sd,
               logit_rho_q = logit_rho_q,
               log_prec_q = log_prec_q,
               sd.df = sd.df,
               par_draws = par_draws)
i = i + 1
car.simu.res[[i]] <- output
}
saveRDS(car.simu.res, file = (paste0(proj_root, "/data/gen/car.simu.res.rds")))
# CAR model estimation assessment
# precision estimation-----------------------------------------------
log_prec_best_vtr <- c()
for (i in 1:length(car.simu.res)) {
  log_prec_best_vtr[i]  <- car.simu.res[[i]]$log_prec_best
}

plot_gamma_prior_post_simu(data = exp(log_prec_best_vtr), 
                     shape = prior_gamma_shape, 
                     scale = prior_gamma_scale, 
                     par = 'Precision', 
                     par_val = prec, 
                     limits = c(0, 1)) 
#hist(log_prec_best_vtr, xlab="Estimated Log-Precision Value", main="Histogram of Estimated Log-Precision Value", col = 'blue')
#abline(v = log_prec, col="red", lwd=3, lty=2)

# Find the corresponding quantile of the true log-precision
log_prec_q_vtr <- c()
for (i in 1:length(car.simu.res)) {
  log_prec_q_vtr[i]  <- car.simu.res[[i]]$log_prec_q
}

hist(log_prec_q_vtr, xlab="True Log-Precision Quantile", main="Histogram of True Log-Precision Quantiles")

# log-precision sd
log_prec_sd_vtr <- c()
for (i in 1:length(car.simu.res)) {
  log_prec_sd_vtr[i] <- car.simu.res[[i]]$log_prec_sd
}

log_prec_df <- data.frame(mean = log_prec_best_vtr, sd = log_prec_sd_vtr, quantile = log_prec_q_vtr)

ggplot(log_prec_df, aes(mean, sd)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = log_prec, linetype="dashed",
             color = "red", size=1.0)+
  labs(title="Estimated Log-Precision Mean vs. SD", x = "Estimated Mean", y = "Estimated SD")  +
  theme_bw() +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=10, face = 'bold'),
        legend.text = element_text(size=9),
        axis.title.x = element_text(size=13, face="bold" ),
        axis.title.y = element_text(size=13, face="bold"),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))

ggplot(log_prec_df, aes(mean, quantile)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = log_prec, linetype="dashed",
             color = "red", size=1.0)+
  labs(title="Estimated Log-Precision Mean vs. \nTrue Quantile", x = "Estimated Mean", y = "Estimated True Quantile")  +
  theme_bw() +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=10, face = 'bold'),
        legend.text = element_text(size=9),
        axis.title.x = element_text(size=13, face="bold" ),
        axis.title.y = element_text(size=13, face="bold"),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))

# rho estimation-----------------------------------------------
# logit rho
logit_rho_best_vtr <- c()
for (i in 1:length(car.simu.res)) {
  logit_rho_best_vtr[i]  <- car.simu.res[[i]]$logit_rho_best
}

plot_norm_prior_post_simu(data = logit_rho_best_vtr, 
                      mean = prior_logit_rho_mean, 
                      sd = prior_logit_rho_sd, 
                      par = 'Logit Lambda', 
                      par_val = logit_rho, 
                      limits = c(-3.5*prior_logit_rho_sd + prior_logit_rho_mean, 
                                 3.5*prior_logit_rho_sd + prior_logit_rho_mean))

#hist(logit_rho_best_vtr, xlab="Estimated Logit-Lambda Value", main="Histogram of Estimated Logit-Lambda Value", col = 'blue')
#abline(v = logit_rho, col="red", lwd=3, lty=2)

# Find the corresponding quantile of the true logit rho
logit_rho_q_vtr <- c()
for (i in 1:length(car.simu.res)) {
  logit_rho_q_vtr[i]  <- car.simu.res[[i]]$logit_rho_q
}

hist(logit_rho_q_vtr, xlab="True Logit-Lambda Quantile", main="Histogram of True Logit-Lambda Quantiles")

# logit-lambda sd
logit_rho_sd_vtr <- c()
for (i in 1:length(car.simu.res)) {
  logit_rho_sd_vtr[i] <- car.simu.res[[i]]$logit_rho_sd
}

logit_rho_df <- data.frame(mean = logit_rho_best_vtr, sd = logit_rho_sd_vtr, quantile = logit_rho_q_vtr)

ggplot(logit_rho_df, aes(mean, sd)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = logit_rho, linetype="dashed",
             color = "red", size=1.0)+
  labs(title="Estimated Logit-Lambda Mean vs. SD", x = "Estimated Mean", y = "Estimated SD")  +
  theme_bw() +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=10, face = 'bold'),
        legend.text = element_text(size=9),
        axis.title.x = element_text(size=13, face="bold" ),
        axis.title.y = element_text(size=13, face="bold"),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))

ggplot(logit_rho_df, aes(mean, quantile)) +
  geom_point() +
  geom_smooth() +
  geom_vline(xintercept = logit_rho, linetype="dashed",
             color = "red", size=1.0)+
  labs(title="Estimated Logit-Lambda Mean vs. \nTrue Quantile", x = "Estimated Mean", y = "Estimated True Quantile")  +
  theme_bw() +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size=10, face = 'bold'),
        legend.text = element_text(size=9),
        axis.title.x = element_text(size=13, face="bold" ),
        axis.title.y = element_text(size=13, face="bold"),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))

# CAR model interpolation assessment
# 3 realizations visualizations
#vis #1--------------------------------------------------
car.vis1 <- car.simu.res[[1]]$car_simu
car.vis1.simu.ras <- rasterize(car.vis1 %>% dplyr::select('simu'), raster.grid, field = 'simu')
plot(car.vis1.simu.ras, main='Simulated CAR Area Random Effects', col = brewer.pal(7, "PuRd"), zlim = c(0, max(car.vis1.simu.ras@data@values)))

car.vis1.pred.ras <- rasterize(car.vis1 %>% dplyr::select('pred'), raster.grid, field = 'pred')
plot(car.vis1.pred.ras, main='Predicted CAR Area Random Effects', col = brewer.pal(7, "PuRd"), zlim = c(0, max(car.vis1.simu.ras@data@values)))

car.vis1.sd.ras <- rasterize(car.vis1$sd %>% dplyr::select(`Std. Error`) %>% bind_cols(car.vis1 %>% dplyr::select('simu','pred')) %>% st_as_sf() %>%
                               rename('std.error' = `Std. Error`), raster.grid, field ='std.error')
plot(car.vis1.sd.ras, main = 'Predicted CAR Area Random Effects \nStandard Error', 
     col = brewer.pal(7, "Greys"), 
     zlim = c(0, max(car.vis1.sd.ras@data@values)))

# True simulated car random effect
# ggplot(car.vis1) +
#   geom_sf(aes(fill = simu), color = NA) +
#   scale_fill_distiller('value', palette = "RdBu", limits = c(0, max(abs(car.vis1$simu)))) +
#   theme_bw() +
#   ggtitle("Simulated CAR Random Effect") +
#   xlab("x") + ylab("y") +
#   theme(plot.title = element_text(size=20),
#         legend.title = element_text(size=15),
#         legend.text = element_text(size=15),
#         axis.title.x = element_text(size=15),
#         axis.title.y = element_text(size=15),
#         axis.text.y = element_text(size = 15),
#         axis.text.x = element_text(size = 15))


# pred car random effect
# ggplot(car.vis1) +
#   geom_sf(aes(fill = pred), color = NA) +
#   scale_fill_distiller('pred', palette = "RdBu", limits = c(0, max(abs(car.vis1$simu)))) +
#   theme_bw() +
#   ggtitle("Predicted CAR Random Effect") +
#   xlab("x") + ylab("y") +
#   theme(plot.title = element_text(size=20),
#         legend.title = element_text(size=15),
#         legend.text = element_text(size=15),
#         axis.title.x = element_text(size=15),
#         axis.title.y = element_text(size=15),
#         axis.text.y = element_text(size = 15),
#         axis.text.x = element_text(size = 15))

# pred car random effect sd
# ggplot(car.vis1) +
#   geom_sf(aes(fill = car.vis1$sd$`Std. Error`), color = NA) +
#   scale_fill_distiller('sd', palette = "RdBu") +
#   theme_bw() +
#   ggtitle("Predicted CAR Standard Deviation") +
#   xlab("x") + ylab("y") +
#   theme(plot.title = element_text(size=20),
#         legend.title = element_text(size=15),
#         legend.text = element_text(size=15),
#         axis.title.x = element_text(size=15),
#         axis.title.y = element_text(size=15),
#         axis.text.y = element_text(size = 15),
#         axis.text.x = element_text(size = 15))

#vis #2-----------------------------------------------------------
car.vis2 <- car.simu.res[[2]]$car_simu
car.vis2.simu.ras <- rasterize(car.vis2 %>% dplyr::select('simu'), raster.grid, field = 'simu')
plot(car.vis2.simu.ras, main='Simulated CAR Area Random Effects', col = brewer.pal(7, "PuRd"), zlim = c(0, max(car.vis2.simu.ras@data@values)))

car.vis2.pred.ras <- rasterize(car.vis2 %>% dplyr::select('pred'), raster.grid, field = 'pred')
plot(car.vis2.pred.ras, main='Predicted CAR Area Random Effects', col = brewer.pal(7, "PuRd"), zlim = c(0, max(car.vis2.simu.ras@data@values)))

car.vis2.sd.ras <- rasterize(car.vis2$sd %>% dplyr::select(`Std. Error`) %>% bind_cols(car.vis2 %>% dplyr::select('simu','pred')) %>% st_as_sf() %>%
                               rename('std.error' = `Std. Error`), raster.grid, field ='std.error')
plot(car.vis1.sd.ras, main = 'Predicted CAR Area Random Effects \nStandard Error', 
     col = brewer.pal(7, "Greys"), 
     zlim = c(0, max(car.vis2.sd.ras@data@values)))

# True simulated car random effect

# ggplot(car.vis2) +
#   geom_sf(aes(fill = simu), color = NA) +
#   scale_fill_distiller('value', palette = "RdBu", limits = c(0, max(abs(car.vis2$simu)))) +
#   theme_bw() +
#   ggtitle("Simulated CAR Random Effect") +
#   xlab("x") + ylab("y") +
#   theme(plot.title = element_text(size=20),
#         legend.title = element_text(size=15),
#         legend.text = element_text(size=15),
#         axis.title.x = element_text(size=15),
#         axis.title.y = element_text(size=15),
#         axis.text.y = element_text(size = 15),
#         axis.text.x = element_text(size = 15))


# pred car random effect
# ggplot(car.vis2) +
#   geom_sf(aes(fill = pred), color = NA) +
#   scale_fill_distiller('pred', palette = "RdBu", limits = c(0, max(abs(car.vis2$simu)))) +
#   theme_bw() +
#   ggtitle("Predicted CAR Random Effect") +
#   xlab("x") + ylab("y") +
#   theme(plot.title = element_text(size=20),
#         legend.title = element_text(size=15),
#         legend.text = element_text(size=15),
#         axis.title.x = element_text(size=15),
#         axis.title.y = element_text(size=15),
#         axis.text.y = element_text(size = 15),
#         axis.text.x = element_text(size = 15))

# pred car random effect sd
# ggplot(car.vis2) +
#   geom_sf(aes(fill = car.vis2$sd$`Std. Error`), color = NA) +
#   scale_fill_distiller('sd', palette = "RdBu") +
#   theme_bw() +
#   ggtitle("Predicted CAR Standard Deviation") +
#   xlab("x") + ylab("y") +
#   theme(plot.title = element_text(size=20),
#         legend.title = element_text(size=15),
#         legend.text = element_text(size=15),
#         axis.title.x = element_text(size=15),
#         axis.title.y = element_text(size=15),
#         axis.text.y = element_text(size = 15),
#         axis.text.x = element_text(size = 15))

#vis #3---------------------------------------------
car.vis3 <- car.simu.res[[3]]$car_simu
car.vis3.simu.ras <- rasterize(car.vis3 %>% dplyr::select('simu'), raster.grid, field = 'simu')
plot(car.vis3.simu.ras, main='Simulated CAR Area Random Effects', col = brewer.pal(7, "PuRd"), zlim = c(0, max(car.vis3.simu.ras@data@values)))

car.vis3.pred.ras <- rasterize(car.vis3 %>% dplyr::select('pred'), raster.grid, field = 'pred')
plot(car.vis3.pred.ras, main='Predicted CAR Area Random Effects', col = brewer.pal(7, "PuRd"), zlim = c(0, max(car.vis3.simu.ras@data@values)))

car.vis3.sd.ras <- rasterize(car.vis3$sd %>% dplyr::select(`Std. Error`) %>% bind_cols(car.vis3 %>% dplyr::select('simu','pred')) %>% st_as_sf() %>%
                               rename('std.error' = `Std. Error`), raster.grid, field ='std.error')
plot(car.vis3.sd.ras, main = 'Predicted CAR Area Random Effects \nStandard Error', 
     col = brewer.pal(7, "Greys"), 
     zlim = c(0, max(car.vis3.sd.ras@data@values)))


# True simulated car random effect

# ggplot(car.vis3) +
#   geom_sf(aes(fill = simu), color = NA) +
#   scale_fill_distiller('value', palette = "RdBu", limits = c(0, max(abs(car.vis3$simu)))) +
#   theme_bw() +
#   ggtitle("Simulated CAR Random Effect") +
#   xlab("x") + ylab("y") +
#   theme(plot.title = element_text(size=20),
#         legend.title = element_text(size=15),
#         legend.text = element_text(size=15),
#         axis.title.x = element_text(size=15),
#         axis.title.y = element_text(size=15),
#         axis.text.y = element_text(size = 15),
#         axis.text.x = element_text(size = 15))


# pred car random effect
# ggplot(car.vis3) +
#   geom_sf(aes(fill = pred), color = NA) +
#   scale_fill_distiller('pred', palette = "RdBu", limits = c(0, max(abs(car.vis3$simu)))) +
#   theme_bw() +
#   ggtitle("Predicted CAR Random Effect") +
#   xlab("x") + ylab("y") +
#   theme(plot.title = element_text(size=20),
#         legend.title = element_text(size=15),
#         legend.text = element_text(size=15),
#         axis.title.x = element_text(size=15),
#         axis.title.y = element_text(size=15),
#         axis.text.y = element_text(size = 15),
#         axis.text.x = element_text(size = 15))

# pred car random effect sd
# ggplot(car.vis3) +
#   geom_sf(aes(fill = car.vis3$sd$`Std. Error`), color = NA) +
#   scale_fill_distiller('sd', palette = "RdBu") +
#   theme_bw() +
#   ggtitle("Predicted CAR Standard Deviation") +
#   xlab("x") + ylab("y") +
#   theme(plot.title = element_text(size=20),
#         legend.title = element_text(size=15),
#         legend.text = element_text(size=15),
#         axis.title.x = element_text(size=15),
#         axis.title.y = element_text(size=15),
#         axis.text.y = element_text(size = 15),
#         axis.text.x = element_text(size = 15))

