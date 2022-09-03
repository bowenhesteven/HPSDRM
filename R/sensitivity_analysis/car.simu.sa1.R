################################################################################# 
## R code for "Incorporating spatial autocorrelation in dasymetric             ## 
## analysis: A hierarchical poisson spatial disaggregation regression model"   ##
## Authors: Bowen He, Jonathan M. Gilligan, Janey V. Camp                      ## 
#################################################################################

################################################################################
#########################  CAR SIMULATION STUDY ################################
################################################################################

# Conduct (hyper)priors' parameters sensitivity analysis on the 
# Leroux CAR simulation model using the same 100
# simulated realizations.

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

car.simu.res <- readRDS(file.path(proj_root, "data/gen/sensitivity_analysis/car.simu.res.gaus.sd1.rds"))
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
prior_gamma_shape <- prec
prior_gamma_scale <- 1
#prior_gamma_scale <- (-prec + sqrt(prec^2 + 4))/2
#prior_gamma_shape <- 1 + prec/prior_gamma_scale

raster <- raster(ncol=30, nrow=30, xmn=0, xmx=1, ymn=0, ymx=1)
area <- raster %>%
  as('SpatialPolygonsDataFrame') 

area.nb <- poly2nb(area)
area.nb_B <- nb2listw(area.nb, style="B", zero.policy=TRUE)

sym.matrix <- as(area.nb_B, "symmetricMatrix")
one.matrix <- rep(1, nrow(area))
a <- Diagonal(x = as.numeric(sym.matrix%*%one.matrix))
mymatrix <- a - sym.matrix

n_realizations <- 100  

# 100 realization experiments
i <- 0
car.simu.res.gaus.sd1 <- list()
while (i < n_realizations) {
  
  print(i)
  car_simu <- car.simu.res[[i+1]]$car_simu  
  
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
  car.simu.res.gaus.sd1[[i]] <- output
}
saveRDS(car.simu.res.gaus.sd1, file = (paste0(proj_root, "/data/gen/sensitivity_analysis/car.simu.res.gaus.sd1.rds")))
