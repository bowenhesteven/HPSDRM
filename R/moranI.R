#############################################################################
## R code for "Incorporating spatial autocorrelation in dasymetric analysis:# 
## A hierarchical spatial regression model"                                 #
## Authors: Bowen He, Jonathan M. Gilligan, Janey V. Camp                  ##  
#############################################################################

## Moran's I Analysis ##

# In this script, I will be conducting Moran's I analysis on the 
# Decennial 2020 Davidson County, Nashville's census blocks population data

## Load required libraries
library(ggplot2)
library(sf)
library(spdep)
library(tmap)
library(tidycensus)
library(units)
proj_root <- find_root(is_rstudio_project)

# Load the shapefile
pop.tract <- get_decennial(geography = "tract",
                                variables = "P1_001N",
                                state = 'tn',
                                county = 'Davidson',
                                year = 2020,
                                cache_table = TRUE, 
                                geometry = TRUE,
                                sumfile = "pl")

pop.cbg <- get_decennial(geography = "block group",
                                variables = "P1_001N",
                                state = 'tn',
                                county = 'Davidson',
                                year = 2020,
                                cache_table = TRUE, 
                                geometry = TRUE,
                                sumfile = "pl")

pop.cb <- get_decennial(geography = "block",
                                variables = "P1_001N",
                                state = 'tn',
                                county = 'Davidson',
                                year = 2020,
                                cache_table = TRUE, 
                                geometry = TRUE,
                                sumfile = "pl")

hist(pop.cb$value, breaks = 35, main="Population Histogram", xlab="Census Block Population")

pop.cb$area <- st_area(pop.cb) %>% 
  set_units(km^2)
pop.cb$pop.density <- pop.cb$value/pop.cb$area

boxplot(pop.cb$pop.density %>% 
          as.vector(), horizontal = TRUE, main="Population Density Boxplot", xlab = 'Block Population Density (1/km^2)')

tm_shape(population_ct) + tm_fill(col="pop.density", style="quantile", n=8, palette="Greens") +
  tm_legend(outside=TRUE)

# Define neighboring polygons
nb <- poly2nb(pop.cb, queen = TRUE)

# Assign weights to the neighbors
lw <- nb2listw(nb, style="W", zero.policy=TRUE)

# Compute the weighted neighbor mean population values
pop.density.lag <- lag.listw(lw, pop.cb$pop.density %>% as.vector())
pop.density.lag

plot(pop.density.lag ~ pop.cb$pop.density %>% as.vector(), pch=16, asp=1, main="Weighted Neighbor Mean Population Density", xlab="Block Population Density (1/km^2)")
M1 <- lm(pop.density.lag ~ pop.cb$pop.density %>% as.vector())
abline(M1, col="blue", lwd = 3, lty = 2)

# Computing the Moran's I statistic
I <- moran(pop.cb$pop.density %>% as.vector(), lw, length(nb), Szero(lw), zero.policy = TRUE)[1]
I

# Performing a hypothesis test
moran.test(pop.cb$pop.density %>% as.vector(),lw, alternative="greater", zero.policy = TRUE)

MC<- moran.mc(pop.cb$pop.density %>% as.vector(), lw, nsim=999, alternative="greater", zero.policy = TRUE)

# View results (including p-value)
MC
plot(MC, xlab = NULL, main = "Density Plot of Permutation Outcomes")


set.seed(131)
pop.cb$rand1 <- sample(pop.cb$pop.density %>% as.vector(), length(pop.cb$pop.density %>% as.vector()), replace = FALSE)
pop.cb$rand2 <- sample(pop.cb$pop.density %>% as.vector(), length(pop.cb$pop.density %>% as.vector()), replace = FALSE)
pop.cb$rand3 <- sample(pop.cb$pop.density %>% as.vector(), length(pop.cb$pop.density %>% as.vector()), replace = FALSE)

tm_shape(pop.cb) + 
  tm_fill(col=c("pop.density", "rand1", "rand2", "rand3"),
  style="quantile", n=8, palette="Blues", legend.show = FALSE) +
  tm_facets(nrow=1)
