########################################################################
## R code for "Incorporating spatial autocorrelation in dasymetric    ## 
## analysis: A hierarchical spatial regression model"                 ##
## Authors: Bowen He, Jonathan M. Gilligan, Janey V. Camp             ## 
########################################################################

########################################################################
############## Intelligent Dasymetric Mapping Model (IDM) ##############
########################################################################

# In this R code file, the Decennial 2020 population of the Davidson County
# Nashville is spatially interpolated from the tract scale to the 5*5*30*30m^2 
# (22500 m^2) resolution scale performed by the Intelligent Dasymetric Mapping 
# (IDM) algorithm proposed by Jeremy Mennis and Torrin Hultgren (2006).

# Then, this spatially interpolated grids population computed by the IDM is 
# aggregated to the Decennial blocks to compare with the true 2020 Decennial 
# blocks population data of the Davidson County, Nashville downloaded from the 
# 2020 Decennial census to evaluate the accuracy and performance of the Mennis
# and Hultgren (2006) IDM spatial interpolation approach.

# Three performance metrics are used: (1) R2; 
# (2) Root Mean Squared Error (RMSE); 
# (3) Mean Absolute Error (MAE)

# citation: Mennis, J., & Hultgren, T. (2006). Intelligent dasymetric mapping and 
# its application to areal interpolation. Cartography and Geographic Information Science, 
# 33(3), 179-194.

## Load required libraries
library(tidycensus)
library(sf)
library(dplyr)
library(tidyverse)
library(gt)
library(rprojroot)
library(raster)
library(rgeos)
library(units)
library(sp)
library(exactextractr)

proj_root <- find_root(is_rstudio_project)
source(paste0(proj_root,"/R/IDM.func.R"))

# Define the CRS of the study region
CRS <- "+proj=utm +zone=16 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
#-------------------Intelligent Dasymetric Mapping Input---------------------------------------
# read in the decennial 2020 population data
pop.tract <- get_decennial(geography = "tract", sumfile = "pl", state = 'TN', county = 'Davidson', variables = "P1_001N", year = 2020, geometry = TRUE) %>%
  mutate(area = st_area(.))%>%
  st_transform(crs = CRS)

# download the Decennial 2020 blocks data of the Davidson County, Nashville
pop.block <- get_decennial(geography = "block", sumfile = "pl", state = 'TN', county = 'Davidson', variables = "P1_001N", year = 2020, geometry = TRUE) %>%
  mutate(area = st_area(.)) %>%
  filter(value > 0)

# read in the decennial 2020 census tract data of the Davidson County, Nashville
census_tract <- st_read(file.path(proj_root, "data", 'src/census_tract/census_tract.shp'))

# read in the 2019 NLCD land cover raster data
nlcd_raster <- raster::raster(file.path(proj_root, "data", 'src/NLCD-2019/nlcd-2019.img'))%>%
  projectRaster(crs = CRS, method="ngb")

# set the threshold as 0.5
threshold <- 0.5

# preset the unhibitated land cover
preset_land_cover_no_people <- c('open_water','herbaceous_wetlands') 

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

#--------Intelligent Dasymetric Mapping Function----------------------
idm_dasy_pop_raster <- IDM(ct_land_cover_df, pop.tract, nlcd_raster, threshold, preset_land_cover_no_people, nlcd_pts, nlcd_grids_ct_df) %>%
  raster::aggregate(fact = 5, fun = sum) # function sourced in the IDM-tool-func.R doc

# plot the dasymetric population raster
plot(idm_dasy_pop_raster, main = 'The IDM Grids Population')

#---------Aggregate the spatially interpolated grids population to 2020 Decennial blocks-------------------
# aggregate the spatially interpolated grids population to 2020 Decennial blocks
idm_dasy_pop_block <- exactextractr::exact_extract(idm_dasy_pop_raster, pop.block, fun = 'sum')

#-------Evaluate the accuracy/performance of the IDM spatially interpolated grids population-----------------
r_square <- function(sample, predict){
  rss <- sum((sample - predict)^2)
  tss <- sum((sample - mean(sample))^2)
  r2 <- 1 - rss/tss
}

# R square 
IDM.r2 <- r_square(pop.block$value, idm_dasy_pop_block)

# RMSE
IDM.rmse <- Metrics::rmse(pop.block$value, idm_dasy_pop_block)

# MAE
IDM.mae <- Metrics::mae(pop.block$value, idm_dasy_pop_block)

# dataframe of true and pred block population
idm.block.df <- pop.block %>%
  rename(true = value) %>%
  cbind(pred = idm_dasy_pop_block)

# plot out-of-sample-block performance
idm.block.obspred <- ggplot(idm.block.df, aes(x = true, y = pred)) + 
  geom_point(color = 'red', alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 1, color = 'blue') + 
  ggtitle('IDM out sample performance:\nBlock Population') +
  theme_bw() +
  theme(plot.title = element_text(size=20, hjust = 0.5, face = 'bold'), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15), 
        axis.title.x = element_text(size=15), 
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size =15), 
        axis.text.x = element_text(size=15))

idm.block.obspred

idm.outsample.metrics <- dplyr::summarise(idm.block.df %>% st_drop_geometry(), 
                                          RMSE = sqrt(mean((pred - true) ^ 2)),
                                          MAE = mean(abs(pred - true)),
                                          pearson = cor(pred, true, method = 'pearson'),
                                          spearman = cor(pred, true, method = 'spearman'),
                                          log_pearson = cor(log1p(pred), log1p(true), method = 'pearson'))


IDM.res <- list(idm_dasy_pop_raster = idm_dasy_pop_raster,
                IDM.r2 = IDM.r2,
                IDM.rmse = IDM.rmse,
                IDM.mae = IDM.mae,
                idm.block.df = idm.block.df
)
saveRDS(IDM.res, file = (paste0(proj_root, "/data/gen/IDM.res.rds")))