################################################################################# 
## R code for "" #Incorporating spatial autocorrelation in dasymetric          ## 
## analysis: A hierarchical poisson spatial disaggregation regression model    ##
## Authors: Bowen He, Jonathan M. Gilligan, Janey V. Camp                      ## 
#################################################################################

# scale function for the HPSDRM
scale_func <- function(pop.tract, pred.raster){
  
  pred.sf <- as.data.frame(rasterToPoints(pred.raster)) %>%
    st_as_sf(coords = c("x", "y")) %>%
    st_set_crs(st_crs(pop.tract))%>%
    mutate(ID = row_number())
  
  # transform to sf object to fasterize
  pred.polygon.sf <- rasterToPolygons(pred.raster) %>%
    st_as_sf() %>%
    mutate(ID = row_number())
  
  # combine the areal observational value to each underlying grids
  factor.df <- st_join(pred.sf %>% st_centroid(), pop.tract, join = st_intersects) %>%
    st_drop_geometry() %>%
    group_by(GEOID) %>%
    summarise(pred.areal.value = sum(layer), factor = value/pred.areal.value, .groups = "drop") %>%
    unique() 
  
  main.df <- st_join(pred.sf %>% st_centroid(), pop.tract, join = st_intersects) %>%
    st_drop_geometry() %>%
    left_join(factor.df, by = 'GEOID') %>%
    mutate(scaled.pred.grid.value = factor * layer) %>%
    left_join(pred.polygon.sf %>% dplyr::select('ID'), by = 'ID') %>%
    st_as_sf()
  
  # Fasterize the HPSRM dasymetric population raster  
  dasy_pop_raster <- fasterize::fasterize(main.df %>% st_as_sf, pred.raster, field = 'scaled.pred.grid.value', fun='last')
} 
# scale function for (SPDE)-CAR simulation model
scale_func_simu <- function(n, simulation.true.raster, pred.raster){
  
  raster <- raster(ncol=n, nrow=n, xmn=0, xmx=1, ymn=0, ymx=1)
  
  area <- raster %>%
    as('SpatialPolygonsDataFrame') %>%
    st_as_sf() %>%
    mutate('layer' = row_number()) %>%
    rename('ID' = 'layer')
  
  simulation.true.sf <- simulation.true.raster %>%
    as('SpatialPolygonsDataFrame') %>%
    st_as_sf() %>%
    rename('true.grid.value'='layer') %>%
    st_centroid()
  
  pred.sf <- pred.raster %>%
    as('SpatialPolygonsDataFrame') %>%
    st_as_sf() %>%
    rename('pred.grid.value'='layer') %>%
    mutate('Grid.ID' = row_number())
  
  pred.sf$true.grid.value <- simulation.true.sf$true.grid.value
  
  factor.df <- st_join(pred.sf, area, join = st_intersects) %>%
    filter(!duplicated(Grid.ID)) %>%
    st_drop_geometry() %>%
    group_by(ID) %>%
    summarise(pred.areal.value = sum(pred.grid.value), true.area.value = sum(true.grid.value), factor = true.area.value/pred.areal.value, .groups = "drop") %>%
    unique()
  
  # combine the areal observational value to each underlying grids
  df <- st_join(pred.sf, area, join = st_intersects) %>%
    filter(!duplicated(Grid.ID)) %>%
    st_drop_geometry() %>%
    left_join(factor.df, by = 'ID') %>%
    mutate(scaled.pred.grid.value = factor * pred.grid.value,
           resi.error = pred.grid.value - true.grid.value,
           scaled.resi.error = scaled.pred.grid.value - true.grid.value ,
           aw.resi.error = true.area.value/n^2 -  true.grid.value) %>%
    left_join(pred.sf %>% dplyr::select('Grid.ID'), by = 'Grid.ID')
  
  r_square_calc <- function(sample, predict){
    rss <- sum((sample - predict)^2)
    tss <- sum((sample - mean(sample))^2)
    r2 <- 1 - rss/tss
  }
  
  areal_weighting_r2 <- r_square_calc(df$true.grid.value, df$true.area.value/n^2)
  areal_weighting_rmse <- Metrics::rmse(df$true.grid.value, df$true.area.value/n^2)
  areal_weighting_mae <- Metrics::mae(df$true.grid.value, df$true.area.value/n^2)
  
  model_r2 <- r_square_calc(df$true.grid.value, df$pred.grid.value)
  model_rmse <- Metrics::rmse(df$true.grid.value, df$pred.grid.value)
  model_mae <- Metrics::mae(df$true.grid.value, df$pred.grid.value)
  
  scaled_model_r2 <- r_square_calc(df$true.grid.value, df$scaled.pred.grid.value)
  scaled_model_rmse <- Metrics::rmse(df$true.grid.value, df$scaled.pred.grid.value)
  scaled_model_mae <- Metrics::mae(df$true.grid.value, df$scaled.pred.grid.value)
  
  sf <- df %>% st_as_sf()
  pred.scaled.raster <- fasterize(sf, pred.raster, field = 'scaled.pred.grid.value', 'last')
  
  result <- list(pred.scaled.raster = pred.scaled.raster,
                 areal_weighting_r2 = areal_weighting_r2,
                 areal_weighting_rmse = areal_weighting_rmse,
                 areal_weighting_mae = areal_weighting_mae,
                 model_r2 = model_r2,
                 model_rmse = model_rmse,
                 model_mae = model_mae,
                 scaled_model_rmse = scaled_model_rmse,
                 scaled_model_r2 = scaled_model_r2,
                 scaled_model_mae = scaled_model_mae,
                 sf = sf)
}