###############################################################################
## R code for "Incorporating spatial autocorrelation in dasymetric analysis: ## 
## A hierarchical spatial regression model"                                  ##
## Authors: Bowen He, Jonathan M. Gilligan, Janey V. Camp                    ##  
###############################################################################

# In this script, I create R functions to perform Intelligent Dasymetric Mapping 
# (IDM) Algorithm by Jeremy Mennis and Torrin Hultgren (2006)

# citation: Mennis, J., & Hultgren, T. (2006). Intelligent dasymetric mapping and its application to areal interpolation. 
# Cartography and Geographic Information Science, 33(3), 179-194.
#---------------------------------------dasymetric mapping func---------------------------------------------------
#------function-create the joined population data and the spatial joined census tract and NLCD land cover raster data----------
pop_ct_land_cover_df <- function(ct_land_cover_df, population) { 
  df <- ct_land_cover_df %>%
  left_join(population %>% dplyr::select('GEOID', "value") %>% st_drop_geometry(), by = 'GEOID') %>%
  mutate(open_water = open_water * 900/area, developed_open_space = developed_open_space * 900/area, developed_low_intensity = developed_low_intensity * 900/area, developed_medium_intensity = developed_medium_intensity * 900/area,
         developed_high_intensity = developed_high_intensity * 900/area, barren = barren * 900/area, deciduous_forest = deciduous_forest * 900/area, evergreen_forest = evergreen_forest * 900/area,
         mixed_forest = mixed_forest * 900/area, shrub = shrub * 900/area, grassland = grassland * 900/area, pasture = pasture * 900/area, cultivated_crops = cultivated_crops * 900 /area,woody_wetlands = woody_wetlands *900/area, herbaceous_wetlands = herbaceous_wetlands * 900/area)
}
# ------------------------------------dasymetric sampling func-------------------------------------------------------
dasy_sample <- function(pop_ct_land_cover_df, threshold = 0.5){
for (i in 1:ncol(pop_ct_land_cover_df)) {
  if (!colnames(pop_ct_land_cover_df[, i])[1] %in% c('CT_NAME', "area", "geometry", "value", "GEOID")) {
    colname <- colnames(pop_ct_land_cover_df[,i])[1] 
  
    df <- pop_ct_land_cover_df %>% 
      filter(unlist(pop_ct_land_cover_df[,i] %>% st_drop_geometry()) > threshold) %>%
      dplyr::select('CT_NAME', "area", "geometry", "value", "GEOID")
 
  density <- sum(df$value)/sum(df$area) 
  assign(colname, list(density = density, df = df))
  }
}

# create sampling result
sampling_result <- list(open_water = open_water, developed_open_space = developed_open_space, developed_low_intensity = developed_low_intensity, developed_medium_intensity = developed_medium_intensity, 
                        developed_high_intensity = developed_high_intensity, barren = barren, deciduous_forest = deciduous_forest, evergreen_forest = evergreen_forest, mixed_forest = mixed_forest,
                        shrub = shrub, grassland = grassland, pasture = pasture, cultivated_crops = cultivated_crops, woody_wetlands = woody_wetlands, herbaceous_wetlands = herbaceous_wetlands)
}
# ---------------------------synthesizing the density of sampled land cover func---------------------------------------
sampled_land_cover_density <- function(pop_ct_land_cover_df, sampled_density, preset_land_cover_no_people){   
  df <- pop_ct_land_cover_df %>% 
  mutate(open_water = sampled_density$open_water$density, herbaceous_wetlands = sampled_density$herbaceous_wetlands$density, 
         developed_open_space = sampled_density$developed_open_space$density, developed_low_intensity = sampled_density$developed_low_intensity$density,
         developed_medium_intensity = sampled_density$developed_medium_intensity$density, developed_high_intensity = sampled_density$developed_high_intensity$density,
         barren = sampled_density$barren$density, deciduous_forest = sampled_density$deciduous_forest$density, evergreen_forest = sampled_density$evergreen_forest$density, 
         mixed_forest = sampled_density$mixed_forest$density, 
         shrub = sampled_density$shrub$density, 
         grassland = sampled_density$grassland$density, 
         pasture = sampled_density$pasture$density, 
         cultivated_crops = sampled_density$cultivated_crops$density, 
         woody_wetlands = sampled_density$woody_wetlands$density) %>%
  drop_units()
  
  df[, preset_land_cover_no_people] <- 0
  
  return(df)
}
#------------------redistribute population into every land cover type and calculate the densities of each land cover type for each census tract---------
redistribute_pop <- function(sampled_land_cover_density_df, ct_land_cover_df){

sampled_land_cover_density <- sampled_land_cover_density_df %>%
  dplyr::select(-'CT_NAME', -"area", -"value", -"GEOID") %>%
  st_drop_geometry()
  
ct_land_cover_df_subset <- ct_land_cover_df %>%
  dplyr::select(-'CT_NAME', -"area", -"GEOID") %>%
  st_drop_geometry() 

# calculate sampled land cover total population
tot_pop <- sampled_land_cover_density * ct_land_cover_df_subset * 900
tot_pop_sampled <- tot_pop %>%
  mutate(tot_pop_sampled = rowSums(., na.rm = TRUE)) 
tot_pop_sampled <- tot_pop_sampled[ , apply(tot_pop_sampled, 2, function(x) !any(is.na(x)))]

# calculate total unsampled land cover population
rest_pop_unsampled <- sampled_land_cover_density_df$value - tot_pop_sampled$tot_pop_sampled

# unsampled land cover area
unsampled_area <- ct_land_cover_df_subset[,!colnames(ct_land_cover_df_subset) %in% colnames(tot_pop_sampled)] 
tot_unsampled_area <- rowSums(unsampled_area, na.rm = TRUE)
  
# calculate unsampled land cover total population
tot_pop_unsampled <- unsampled_area * rest_pop_unsampled/tot_unsampled_area

# combine the two dataframes: unsampled land cover population and sampled land cover population
df1 <- cbind(tot_pop_unsampled, tot_pop_sampled) %>%
   cbind(sampled_land_cover_density_df %>% dplyr::select('CT_NAME', "area", "value", "GEOID"), rest_pop_unsampled)
 
df2 <- cbind(unsampled_area, rest_pop_unsampled)
 
ouput = list(redistribute_pop = df1, unsampled_area = df2)
}
#-------------------------------calculate the unsampled land cover density----------------------------------------------------
unsampled_land_cover_density <- function(redistribute_pop){
  
unsampled_density_filter_ct <- redistribute_pop$redistribute_pop %>%
  filter(rest_pop_unsampled > 0) 

unsampled_area_population <- unsampled_density_filter_ct[,colnames(unsampled_density_filter_ct) %in% colnames(redistribute_pop$unsampled_area)]  
 
unsampled_area_filter_ct <- redistribute_pop$unsampled_area %>%
   filter(rest_pop_unsampled > 0) 

tot_pop_unsampled <- apply(unsampled_area_population 
      %>% dplyr::select(-'rest_pop_unsampled'), 2, sum, na.rm = TRUE)

tot_pop_area <-  apply(unsampled_area_filter_ct 
                       %>% dplyr::select(-'rest_pop_unsampled'), 2, sum)

unsampled_area_density <- tot_pop_unsampled/(tot_pop_area*900)

unsampled_area_density[is.na(unsampled_area_density)] <- 0

return(unsampled_area_density)
}
# ------------------------------comprehensive land cover density dataframe----------------------------------------------
density <- function(sampled_land_cover_density_df, unsampled_area_density){
   
  unsampled_density_df <- sampled_land_cover_density_df[, apply(sampled_land_cover_density_df, 2, function(x) any(is.na(x)))] %>%
  st_drop_geometry()

  sampled_density_df <- sampled_land_cover_density_df[, apply(sampled_land_cover_density_df, 2, function(x) !any(is.na(x)))]
  
for (i in 1:ncol(unsampled_density_df)) {
  unsampled_density_df[,i] <- unsampled_area_density[i]
}

 density <- st_as_sf(cbind(unsampled_density_df, sampled_density_df))
}
#--------------calculate the dasymetric population for each land cover grid for each census tract------------------------
pop_density <- function(density, ct_land_cover_df){
pop_density <- density %>%
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
  dplyr::select(!contains("area"))

df <- pop_density %>% 
  pivot_longer(!c(GEOID, CT_NAME, value, tot_pop), names_to = 'land_use', values_to = 'pop.density.value')
}
# --------------------------------Generate the dasymetric population distribution raster------------------------------------------------------------
# create the final IDM function that creates the IDM produced dasymetric grid population raster 
dasy_pop_raster_func <- function(nlcd_grids_ct_df, pop_density, nlcd_pts){
# change the name of the land cover type and left join the population density dataframe
nlcd_grids_ct_pop_df <- nlcd_grids_ct_df %>% 
  left_join(pop_density, by = c("GEOID", "CT_NAME","land_use"))

coords <- nlcd_pts %>%
  st_coordinates()

# left join grids geometry
dasy_pop <- nlcd_grids_ct_pop_df %>%
  cbind(coords) %>%
  dplyr::select('X', 'Y', 'pop.density.value') 

dasy_pop_raster <- rasterFromXYZ(dasy_pop)
}

# --------------------------------IDM function---------------------------------
IDM <- function(ct_land_cover_df, population, nlcd_raster, threshold, preset_land_cover_no_people, nlcd_pts, nlcd_grids_ct_df){
  pop_ct_land_cover_df <- pop_ct_land_cover_df(ct_land_cover_df, population)
  sampled_density <- dasy_sample(pop_ct_land_cover_df, threshold = threshold)
  sampled_land_cover_density_df <- sampled_land_cover_density(pop_ct_land_cover_df, sampled_density, preset_land_cover_no_people)  
  redistribute_pop <- redistribute_pop(sampled_land_cover_density_df, ct_land_cover_df)
  unsampled_area_density <- unsampled_land_cover_density(redistribute_pop)
  density <- density(sampled_land_cover_density_df, unsampled_area_density)
  pop_density <- pop_density(density, ct_land_cover_df)
  dasy_pop_raster <- dasy_pop_raster_func(nlcd_grids_ct_df, pop_density, nlcd_pts)
}