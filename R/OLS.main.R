###############################################################################
## R code for "Incorporating spatial autocorrelation in dasymetric           ## 
## analysis: A hierarchical poisson spatial disaggregation regression model" ##
## Authors: Bowen He, Jonathan M. Gilligan, Janey V. Camp                    ## 
###############################################################################

###############################################################################
########### Ordinary Least Squares Linear Regression Model (OLS) ##############
###############################################################################

# In this R code file, the Decennial 2020 population of the Davidson County
# Nashville is spatially interpolated from the tract scale to the 5*5*30*30m^2 
# (22500 m^2) resolution scale performed by the Ordinary Linear Regression Model (OLS).

# Then, this spatially interpolated grids population computed by the OLS (ordinary least squares)
# is aggregated to the Decennial blocks to compare with the true 2020 Decennial blocks population
# data of the Davidson County, Nashville downloaded from the 2020 Decennial census to evaluate 
# the accuracy and performance of the OLS dasymetric mapping model.

# Three performance metrics are used: (1) R2;
# (2) Root Mean Squared Error (RMSE); 
# (3) Mean Absolute Error (MAE)

# OLS (ordinary least squares) Linear Regression Model
# The grid i population density is modeled:

# Y_i = LC*Beta
# LC is the selected significant land cover covariates value, Beta is the land covariates
# coefficients. And the grids population is scaled within the tract to perform the pycnophylactic
# property.

# Based on the regularized analysis results, 3 land cover 
# covariates are selected as the significant land cover 
# predictors in the OLS linear regression model.

# 3 land cover covariates in the model as the predictors:
# (1) Developed open space
# (2) Developed low intensity
# (3) Developed medium intensity

## Load required libraries
library('sf')
library('tidyverse')
library('dplyr')
library('INLA') 
library('raster') 
library('maptools') 
library('gtools')
library('sp') 
library('spdep') 
library('parallel')
library('Metrics')
library('exactextractr')
library('ggplot2')
library('tidycensus')
library('rprojroot')
library('units')

proj_root <- find_root(is_rstudio_project)
source(paste0(proj_root,"/R/IDM.func.R"))

# Define the CRS of the study region
CRS <- "+proj=utm +zone=16 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#-------------------Read in OLS dasymetric mapping Input---------------------------------------
# read in the decennial 2020 population data
pop.tract <- get_decennial(geography = "tract", sumfile = "pl", state = 'TN', county = 'Davidson', variables = "P1_001N", year = 2020, geometry = TRUE) %>%
  mutate(area = st_area(.)) %>%
  st_transform(crs = CRS)

# download the Decennial 2020 blocks data of the Davidson County, Nashville
pop.block <- get_decennial(geography = "block", sumfile = "pl", state = 'TN', county = 'Davidson', variables = "P1_001N", year = 2020, geometry = TRUE) %>%
  mutate(area = st_area(.)) %>%
  filter(value > 0) %>%
  st_transform(crs = CRS)

# read in the 2019 NLCD land cover raster data
nlcd_raster <- raster::raster(file.path(proj_root, "data", 'src/NLCD-2019/nlcd-2019.img')) %>%
  projectRaster(crs = CRS, method="ngb")

# read in the joined NLCD land cover grids-census tract polygon dataframe 
nlcd_grids_ct_df <- read_rds(file.path(proj_root, "data", 'src/nlcd_grids_ct_df.Rds'))
nlcd_grids_ct_df <- nlcd_grids_ct_df$nlcd_grids_ct_df 

# read in the 2019 NLCD land cover points data
nlcd_pts <- read_rds(file.path(proj_root, "data", 'src/nlcd_pts.Rds'))
nlcd_pts <- nlcd_pts$nlcd_pts

# read in the joind census tract polygon-NLCD land cover summary dataframe
ct_land_cover_df <- read_rds(file.path(proj_root, "data", 'src/ct_land_cover_df.Rds')) 
ct_land_cover_df <- ct_land_cover_df$ct_land_cover_df %>%
  st_transform(crs = CRS)
ct_land_cover_pop_df <- ct_land_cover_df %>%
  left_join(pop.tract %>% 
              st_drop_geometry() %>%
              dplyr::select(GEOID, pop.tract = value), by = 'GEOID')

#---------------------OLS dasymetric mapping function/process---------------------
ols.model <- lm(pop.tract ~ 0 + developed_open_space + 
                  developed_low_intensity + 
                  developed_medium_intensity, 
                data = ct_land_cover_pop_df)

# Redistribute population based on the land cover coefficient
ols.density <- ct_land_cover_df %>%
  mutate(developed_open_space = ols.model$coefficients[[1]], 
         developed_low_intensity = ols.model$coefficients[[2]],
         developed_medium_intensity = ols.model$coefficients[[3]],
         developed_high_intensity = 0, 
         open_water = 0, 
         deciduous_forest = 0, 
         evergreen_forest = 0, 
         mixed_forest = 0, 
         grassland = 0, 
         pasture = 0, 
         herbaceous_wetlands = 0, 
         barren = 0, 
         cultivated_crops = 0, 
         shrub = 0,
         woody_wetlands = 0) %>% 
  left_join(pop.tract %>% dplyr::select(GEOID, value) %>% st_drop_geometry(), by = 'GEOID')

# Calc land cover density of grids in each tract------------------------------------------------------------------
ols_pop_density <- ols.density %>%
  st_drop_geometry() %>%
  left_join(ct_land_cover_df %>% st_drop_geometry() %>% dplyr::select(GEOID, 
                                                                      open_water_area = open_water, developed_open_space_area = developed_open_space,
                                                                      developed_low_intensity_area = developed_low_intensity, developed_medium_intensity_area = developed_medium_intensity,
                                                                      developed_high_intensity_area = developed_high_intensity, deciduous_forest_area = deciduous_forest, 
                                                                      evergreen_forest_area = evergreen_forest, mixed_forest_area = mixed_forest, 
                                                                      grassland_area = grassland, pasture_area = pasture, herbaceous_wetlands_area = herbaceous_wetlands,
                                                                      barren_area = barren, cultivated_crops_area = cultivated_crops, shrub_area = shrub,
                                                                      woody_wetlands_area = woody_wetlands), by = 'GEOID') %>%
  mutate(tot_pop = open_water * open_water_area + developed_open_space * developed_open_space_area + developed_low_intensity * developed_low_intensity_area
         + developed_medium_intensity * developed_medium_intensity_area + developed_high_intensity * developed_high_intensity_area + deciduous_forest * deciduous_forest_area
         + evergreen_forest * evergreen_forest_area + mixed_forest*mixed_forest_area + grassland* grassland_area
         + pasture*pasture_area + herbaceous_wetlands * herbaceous_wetlands_area + barren * barren_area + cultivated_crops*cultivated_crops_area + shrub*shrub_area + woody_wetlands*woody_wetlands_area)%>%
  mutate(open_water = value*open_water/tot_pop, developed_open_space = value*developed_open_space/tot_pop,
         developed_low_intensity = value*developed_low_intensity/tot_pop, developed_medium_intensity = value*developed_medium_intensity/tot_pop,
         developed_high_intensity = value*developed_high_intensity/tot_pop, deciduous_forest = value*deciduous_forest/tot_pop,
         evergreen_forest = value*evergreen_forest/tot_pop, mixed_forest = value*mixed_forest/tot_pop, grassland = value*grassland/tot_pop,
         pasture = value*pasture/tot_pop, herbaceous_wetlands = value*herbaceous_wetlands/tot_pop, barren = value*barren/tot_pop, 
         cultivated_crops  = value*cultivated_crops/tot_pop, shrub = value*shrub/tot_pop, woody_wetlands = value*woody_wetlands/tot_pop) %>%
  dplyr::select(!contains("area")) %>% 
  pivot_longer(!c(GEOID, CT_NAME, value, tot_pop), names_to = 'land_use', values_to = 'pop.density.value')

# calculate the final OLS dasymetric population raster
ols_dasy_pop_raster <- dasy_pop_raster_func(nlcd_grids_ct_df, ols_pop_density, nlcd_pts) %>%
  raster::aggregate(fact = 5, fun = sum) # function sourced in the IDM-tool-func.R doc

# plot the OLS dasymetric population raster
plot(ols_dasy_pop_raster, main = 'The OLS Mapped Grids Population')

#---------Aggregate the spatially interpolated grids population to 2020 Decennial blocks-------------------
# aggregate the spatially interpolated grids population to 2020 Decennial blocks
ols_dasy_pop_block <- exactextractr::exact_extract(ols_dasy_pop_raster, pop.block, fun = 'sum')

#-------Evaluate the accuracy/performance of the OLS spatially interpolated grids population-----------------
r_square <- function(sample, predict){
  rss <- sum((sample - predict)^2)
  tss <- sum((sample - mean(sample))^2)
  r2 <- 1 - rss/tss
}

# R square 
OLS.r2 <- r_square(pop.block$value, ols_dasy_pop_block)

# RMSE
OLS.rmse <- Metrics::rmse(pop.block$value, ols_dasy_pop_block)

# MAE
OLS.mae <- Metrics::mae(pop.block$value, ols_dasy_pop_block)

# dataframe of true and pred block population
ols.block.df <- pop.block %>%
  rename(true = value) %>%
  cbind(pred = ols_dasy_pop_block)

# plot out-of-sample-block performance
ols.block.obspred <- ggplot(ols.block.df, aes(x = true, y = pred)) + 
  geom_point(color = 'brown', alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 1, color = 'blue') + 
  ggtitle('OLS out sample performance:\nBlock Population') +
  theme_bw() +
  theme(plot.title = element_text(size=20, hjust = 0.5, face = 'bold'), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15), 
        axis.title.x = element_text(size=15), 
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size =15), 
        axis.text.x = element_text(size=15))

ols.block.obspred

ols.outsample.metrics <- dplyr::summarise(ols.block.df %>% st_drop_geometry(), 
                                             RMSE = sqrt(mean((pred - true) ^ 2)),
                                             MAE = mean(abs(pred - true)),
                                             pearson = cor(pred, true, method = 'pearson'),
                                             spearman = cor(pred, true, method = 'spearman'),
                                             log_pearson = cor(log1p(pred), log1p(true), method = 'pearson'))

OLS.res <- list(ols_dasy_pop_raster = ols_dasy_pop_raster,
                OLS.r2 = OLS.r2,
                OLS.rmse = OLS.rmse,
                OLS.mae = OLS.mae,
                ols.block.df = ols.block.df
)
saveRDS(OLS.res, file = (paste0(proj_root, "/data/gen/OLS.res.rds")))