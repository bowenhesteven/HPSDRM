#############################################################################
## R code for "Incorporating spatial autocorrelation in dasymetric analysis:# 
## A hierarchical spatial regression model"                                 #
## Authors: Bowen He, Jonathan M. Gilligan, Janey V. Camp                  ##  
#############################################################################

# Regularized Regression Analysis-
# Ridge, Lasso, ElasticNet Regression Analysis #

# In this R code document, the decennial 2020 census block population in the 
# Davidson County, Nashville (as the dependent variable data) and the NLCD 
# land cover data (as the independent data) is analyzed in the three regularized 
# regression algorithms: (1) Ridge, (2) Lasso and (3) ElasticNet Regression
# to identify the optimal land cover predictors/covariates combination for the proposed 
# poisson spatial regression model.

# Based on the regularized regression results, 8 land cover covariates are selected:
# (1) Developed open space
# (2) Developed low intensity
# (3) Developed medium intensity
# (4) Developed high intensity
# (5) Shrub
# (6) Grassland
# (7) Barren
# (8) Cultivated Crops

## Load required libraries
library(glmnet)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(foreach)
library(parallel)
proj_root <- find_root(is_rstudio_project)

## Prepare data set 

# Read in the cb_ct_land_cover_pop_df.Rds
cb_ct_land_cover_pop_df <- readRDS(file.path(proj_root, "data", 'src/cb_ct/cb_ct_land_cover_pop_df.Rds')) %>%
  st_drop_geometry()

# land cover predictors/covariates are from the column 6 to column 20 in the dataframe cb_ct_land_cover_pop_df
x <- model.matrix(~.-1, data=cb_ct_land_cover_pop_df[, c(6:20)]) 

# the dependent variable is the population which is the column 2 in the dataframe cb_ct_land_cover_pop_df
y <- cb_ct_land_cover_pop_df %>% 
  dplyr::select(population) %>%
  as.matrix()

# Lasso Regression Analysis---------------------------------------------------
# Run cross-validation

# The cv.glmnet function does k-fold cross-validation for glmnet. more information can be found at:
# https://www.rdocumentation.org/packages/glmnet/versions/4.1-3/topics/cv.glmnet
lasso.mod_cv <- cv.glmnet(x=x, y=y, family='gaussian')

# Combine the lambda value used in the fits and the mean cross-validated error value cvm in a tibble
lasso.lambda_cvm <- tibble(lambda = lasso.mod_cv$lambda, cvm = lasso.mod_cv$cvm)

# plot the lambda vs. cvm using ggplot function
ggplot(data=lasso.lambda_cvm, aes(x=log(lambda), y=cvm, group=1)) +
  geom_line()+
  geom_point(color = "blue", shape=18, size = 3) +
  geom_vline(xintercept=log(lasso.mod_cv$lambda)[lasso.mod_cv$cvm == min(lasso.mod_cv$cvm)],color='red') +
  ggtitle("Cross-Validation Mean Error (CVME) vs. Log(lambda)") + xlab("Log(lambda)") + ylab("CVME") + 
  theme_classic() + 
  theme(plot.title = element_text(size= 20, face="bold", hjust = 0.5), 
        axis.title.x = element_text(size= 20, face="bold" ), 
        axis.title.y = element_text(size= 20, face="bold"), 
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))

# retrieve the value of 'lambda' that gives the minimum 'cvm'
lasso.mod_cv$lambda.min

coef(lasso.mod_cv, lasso.mod_cv$lambda.min)

# plot coefficient vs. log(lambda)

# retrieve the regression coefficients for each land cover predictors/covariants for each lambda value fitted in the experiments
lasso.betas <- as.matrix(lasso.mod_cv$glmnet.fit$beta)

# retrieve the lambda values fitted in the experiments
lasso.lambdas <- lasso.mod_cv$lambda

# assign a name/index for each lambda experiment
names(lasso.lambdas) = colnames(lasso.betas)

# Assemble the dataframe that contains the regression coefficients for each land cover predictors/covariants for each lambda value 
# in the experiments and also the fitted lambda value for each experiments
lasso.betas <- as.data.frame(lasso.betas) %>% 
  tibble::rownames_to_column("variable") %>% 
  pivot_longer(-variable) %>% 
  mutate(lambda=lasso.lambdas[name]) %>%
  rename(`Land Cover` = variable)

# plot the land cover predictors coefficients vs. log(lambda) using ggplot function
  ggplot(lasso.betas %>% filter(!`Land Cover` %in% c('open_water', 'herbaceous_wetlands')), aes(x=log(lambda), y=value, color = `Land Cover`, linetype = `Land Cover`)) + 
  geom_line() +
  geom_vline(xintercept=log(lasso.mod_cv$lambda)[lasso.mod_cv$cvm == min(lasso.mod_cv$cvm)], color = 'red') +
  xlab("Log(lambda)") + 
  ylab("Coefficients") +
  ggtitle('Regression Coefficients vs. Log(Lambda)') + theme_bw() +
    theme(plot.title = element_text(size = 25, face="bold", hjust = 0.5), 
        legend.title = element_text(size = 15, face = 'bold'), 
        legend.text = element_text(size = 13), 
        axis.title.x = element_text(size = 20, face="bold" ), 
        axis.title.y = element_text(size = 20, face="bold"),
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20)) 

# Ridge Regression Analysis------------------------------------------------------
# The ridge regression analysis follows a very similar procedure like the previous lasso regression analysis, 
# except that alpha = 0 makes the difference in the cv.glmnet() function

# The alpha = 0 in the cv.glmnet() function conducts the ridge regression analysis, the default alpha=1 is the lasso regression analysis
ridge.mod_cv <- cv.glmnet(x=x, y=y, family='gaussian', alpha = 0)
ridge.lambda_cvm <- tibble(lambda = ridge.mod_cv$lambda, cvm = ridge.mod_cv$cvm)
ggplot(data=ridge.lambda_cvm, aes(x=log(lambda), y=cvm, group=1)) +
  geom_line()+
  geom_point(color = "blue", shape=18, size = 3) +
  geom_vline(xintercept=log(ridge.mod_cv$lambda)[ridge.mod_cv$cvm == min(ridge.mod_cv$cvm)], color='red') +
  ggtitle("Cross-Validation Mean Error (CVME) vs. Log(lambda)") + xlab("Log(lambda)") + ylab("CVME") + 
  theme_classic() + 
  theme(plot.title = element_text(size=20, face="bold", hjust = 0.5), 
        axis.title.x = element_text(size=20, face="bold" ), 
        axis.title.y = element_text(size=20, face="bold"), 
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))

coef(ridge.mod_cv, ridge.mod_cv$lambda.min)

# plot coefficient vs. log(lambda)

ridge.betas <- as.matrix(ridge.mod_cv$glmnet.fit$beta)
ridge.lambdas <- ridge.mod_cv$lambda
names(ridge.lambdas) = colnames(ridge.betas)

ridge.betas <- as.data.frame(ridge.betas) %>% 
  tibble::rownames_to_column("variable") %>% 
  pivot_longer(-variable) %>% 
  mutate(lambda=ridge.lambdas[name]) %>%
  rename(`Land Cover` = variable)

ggplot(ridge.betas %>% filter(!`Land Cover` %in% c('open_water', 'herbaceous_wetlands')), aes(x=log(lambda), y=value, color = `Land Cover`, linetype = `Land Cover`)) + 
  geom_line() +
  geom_vline(xintercept=log(ridge.mod_cv$lambda)[ridge.mod_cv$cvm == min(ridge.mod_cv$cvm)], color = 'red') +
  xlab("Log(lambda)") + 
  ylab("Coefficients") +
  ggtitle('Regression Coefficients vs. Log(Lambda)') + theme_bw() +
  theme(plot.title = element_text(size = 25, face="bold", hjust = 0.5), 
        legend.title = element_text(size = 15, face = 'bold'), 
        legend.text = element_text(size = 13), 
        axis.title.x = element_text(size = 20, face="bold" ), 
        axis.title.y = element_text(size = 20, face="bold"),
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20)) 

# ElasticNet Regression Analysis--------------------------------------------------------------

# Search for the best alpha value, the default alpha = 1 is the lasso, the alpha = 0 is the ridge, and when the alpha is between 0 and 1,
# the ElasticNet regression analyze the data. Now, we have two parameters to tune, the optimal alpha and the optimal lambda

# First, we test the alpha with a 0.05 interval from 0-1.
alpha <- seq(0, 1, 0.05)

# retrieve the lambda value that gives the minimum cvm for each tested alpha value set up above  
search <- foreach(i = alpha, .combine = rbind) %dopar% {
  cv <- cv.glmnet(x, y, family = "gaussian", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.min = cv$lambda.min, alpha = i)
}

# retrieve the smallest cvm and its associated lambda and alpha 
ENmodel_parameter <- search[search$cvm == min(search$cvm), ]

# build the penalty model based on the chosen lambda and alpha value calculated above
ENmodel <- glmnet(x, y, family = "gaussian", lambda = ENmodel_parameter$lambda.min, alpha = ENmodel_parameter$alpha)
coef(ENmodel)

# plot the alpha vs. lambda and CVME
search <- search %>% 
  pivot_longer(!alpha, names_to = 'variable', values_to = 'value')

# plot two subplots in one big ggplot: (1) alpha vs. CVM; (2) alpha vs. lambda.min (the optimal lambda which gives the smallest cvm) 
# so that you can clearly see how alpha can alter the choice of the optimal lambda and where the cvm reach the minimum, the associated lambda
# and the alpha is the one.
ggplot(data = search, aes(x = alpha, y = value))+
  geom_line() +  geom_vline(xintercept=ENmodel_parameter$alpha, color = 'red') +
  facet_wrap(~variable, scales = "free", labeller = as_labeller(c(`cvm` = 'CVME', `lambda.min` = 'Best Lambda'))) +
  xlab("Alpha") + 
  ylab("Value") +
  ggtitle('Elastic Net Regression Alpha Experiments Results') + theme_bw() +
  theme(plot.title = element_text(size = 25, face="bold", hjust = 0.5), 
        legend.title = element_text(size = 15, face = 'bold'), 
        legend.text = element_text(size = 13), 
        axis.title.x = element_text(size = 20, face="bold" ), 
        axis.title.y = element_text(size = 20, face="bold"),
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20))
