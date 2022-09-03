################################################################################
###### R code for "Incorporating spatial autocorrelation in dasymetric    ###### 
###### analysis: A hierarchical spatial regression model"                 ######
#####  Authors: Bowen He, Jonathan M. Gilligan, Janey V. Camp             ###### 
################################################################################

################################################################################
###################### Areal Weighting Mapping (AWM) ###########################
################################################################################


# In this R code file, the Decennial 2020 population of the Davidson County
# Nashville is spatially interpolated from the tract scale to the  5*5*30*30m^2 
# (22500 m^2) resolution scale performed by the traditional Areal Weighting Mapping (AWM).


# Then, this spatially interpolated grids population interpolated by the Areal
# Weighting Mapping (AWM) is aggregated to the 2020 Decennial blocks to compare 
# with the true 2020 Decennial blocks population data of the Davidson County, Nashville 
# downloaded from the 2020 Decennial census to evaluate the accuracy and performance of 
# the traditional AWM approach.

# Three performance metrics are used: (1) R2; 
# (2) Root Mean Squared Error (RMSE); 
# (3) Mean Absolute Error (MAE)


# The Areal Weighting Mapping (AWM) model
# The assumption is that the population is homogeneously spatially distributed across the grids 
# with no variances within the tract except that the predefined uninhabited land cover (open water and
# herbaceous wetlands) have absolute no population distributed.

# In the Areal Weighting Mapping (AWM) model, all the land cover type is included in the model except 
# that no population is distributed on the (1) open water; and (2) herbaceous wetlands land cover type.


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
#-------------------AWM dasymetric mapping Input---------------------------------------
# read in the decennial 2020 population data
pop.tract <- get_decennial(geography = "tract", sumfile = "pl", state = 'TN', county = 'Davidson', variables = "P1_001N", year = 2020, geometry = TRUE) %>%
  mutate(area = st_area(.))

# download the Decennial 2020 blocks data of the Davidson County, Nashville
pop.block <- get_decennial(geography = "block", sumfile = "pl", state = 'TN', county = 'Davidson', variables = "P1_001N", year = 2020, geometry = TRUE) %>%
  mutate(area = st_area(.)) %>%
  filter(value > 0)

# read in the 2019 NLCD land cover raster data
nlcd_raster <- raster::raster(file.path(proj_root, "data", 'src/NLCD-2019/nlcd-2019.img'))%>%
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


# Redistribute population based on the Areal Weighting Mapping (AWM) assumption
awm.density <- ct_land_cover_df %>%
  left_join(pop.tract %>% dplyr::select(GEOID, value) %>% st_drop_geometry(), by = 'GEOID') %>%
  mutate(tot.grids = developed_open_space + developed_low_intensity + developed_medium_intensity + developed_high_intensity + deciduous_forest + evergreen_forest + mixed_forest 
          + grassland + pasture + barren + cultivated_crops + shrub + woody_wetlands, awm.density = value/tot.grids ) %>%
  mutate(developed_open_space = awm.density, 
         developed_low_intensity = awm.density, 
         developed_medium_intensity = awm.density, 
         developed_high_intensity = awm.density, 
         open_water = 0, 
         deciduous_forest = awm.density, 
         evergreen_forest = awm.density, 
         mixed_forest = awm.density, 
         grassland = awm.density, 
         pasture = awm.density, 
         herbaceous_wetlands = 0,
         barren = awm.density, 
         cultivated_crops = awm.density, 
         shrub = awm.density,
         woody_wetlands = awm.density) %>%
  dplyr::select(-tot.grids, -awm.density)

awm_pop_density <- awm.density %>%
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
awm_dasy_pop_raster <- dasy_pop_raster_func(nlcd_grids_ct_df, awm_pop_density, nlcd_pts) %>%
  raster::aggregate(fact = 5, fun = sum) # function sourced in the IDM-tool-func.R doc# function sourced in the IDM-tool-func.R doc

# plot the Areal Weighting Dasymetric population raster
plot(awm_dasy_pop_raster, main = 'The AWM Grids Population')

#---------Aggregate the spatially interpolated grids population to 2020 Decennial blocks-------------------
# aggregate the spatially interpolated grids population to 2020 Decennial blocks
awm_dasy_pop_block <- exactextractr::exact_extract(awm_dasy_pop_raster, pop.block, fun = 'sum')

#-------Evaluate the accuracy/performance of the IDM spatially interpolated grids population-----------------
r_square <- function(sample, predict){
  rss <- sum((sample - predict)^2)
  tss <- sum((sample - mean(sample))^2)
  r2 <- 1 - rss/tss
}

# R square 
AWM.r2 <- r_square(pop.block$value, awm_dasy_pop_block)

# RMSE
AWM.rmse <- Metrics::rmse(pop.block$value, awm_dasy_pop_block)

# MAE
AWM.mae <- Metrics::mae(pop.block$value, awm_dasy_pop_block)

# dataframe of true and pred block population
awm.block.df <- pop.block %>%
  rename(true = value) %>%
  cbind(pred = awm_dasy_pop_block)

# plot out-of-sample-block performance
awm.block.obspred <- ggplot(awm.block.df, aes(x = true, y = pred)) + 
  geom_point(color = 'green', alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 1, color = 'blue') + 
  ggtitle('AWM out sample performance:\nBlock Population') +
  theme_bw() +
  theme(plot.title = element_text(size=20, hjust = 0.5, face = 'bold'), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15), 
        axis.title.x = element_text(size=15), 
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size =15), 
        axis.text.x = element_text(size=15))

awm.block.obspred

awm.outsample.metrics <- dplyr::summarise(awm.block.df %>% st_drop_geometry(), 
                                          RMSE = sqrt(mean((pred - true) ^ 2)),
                                          MAE = mean(abs(pred - true)),
                                          pearson = cor(pred, true, method = 'pearson'),
                                          spearman = cor(pred, true, method = 'spearman'),
                                          log_pearson = cor(log1p(pred), log1p(true), method = 'pearson'))


AWM.res <- list(awm_dasy_pop_raster = awm_dasy_pop_raster,
                AWM.r2 = AWM.r2,
                AWM.rmse = AWM.rmse,
                AWM.mae = AWM.mae,
                awm.block.df = awm.block.df
)
saveRDS(AWM.res, file = (paste0(proj_root, "/data/gen/AWM.res.rds")))