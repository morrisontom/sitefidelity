#' ---
#' title: Introduction to spExtract package
#' author: Heather McDevitt
#' date: March 2020
#' output: html_document
#' ---
 
# This script contains the code used in Intro_spExtract.html.
# See that file for more detailed explanations.


### Install & load the package ----

# To install a package from GitHub you need to use the install_github() 
# function in the devtools package.

library(devtools)
install_github("hmcdevitt/spExtract", build_vignettes = TRUE)
# you need to specifiy to also install the vignettes when installing from github

# The package can be loaded as normal from now on
library(spExtract) 

# To check if there have been updates and install them use
devtools::update_packages(packages = "spExtract")

# See list of functions in package and their documentation
help(package = "spExtract")


### Extract from a matrix ----

# load elevation raster
library(raster)
elevation <- raster::raster("data/DEM_Serengeti.tif")

# import location data
locations <- read.csv("data/locs.csv")

# turn into spatial object for plotting
library(sp)
loc_sp <- sp::SpatialPoints(coords = locations[,c("x","y")],
                            proj4string = raster::crs(elevation)) 

# plot point on elevation raster
plot(elevation)
plot(loc_sp, pch=20, add=T)

# create matrix
elev_mat <- raster::as.matrix(elevation)
# create extent object
elev_ext <- raster::extent(elevation)
# create cell resolution object
elev_res <- raster::res(elevation)

# extract elevation values at locations
locations$elev <- spExtract::extract_matrix(matrix = elev_mat, 
                                            coords = locations[,c("x","y")], 
                                            raster_extent = elev_ext, 
                                            raster_res = elev_res)

# view dataframe
head(locations, 10)


### Extract from an array ----

# to be added


### Extract NDVI values through time ----

# Generate description of which files to compile into the NDVI array
x<-Sys.time()
files_desc <- spExtract::build_ndvi_desc(folder = "data/", # folder where NDVI tiffs are stored
                                         start.date = "2018-01-01", # start date of data
                                         end.date = "2018-03-31") # end date of data

# Build array of NDVI layers for files identified in ndvi_files
ndvi_array <- suppressWarnings(spExtract::build_ndvi_array(files_desc = files_desc)) 
Sys.time()-x

x<-Sys.time()
w <- list.files(path='./data',pattern='MOD')
st <- stack(paste0('./data/',w))
q <- extract(st,loc_sp)
Sys.time()-x

# note: 
#the warning messages being suppressed here are not a problem and will be fixed inside the function soon


# Import mean NDVI raster for aNDVI calculation
mean_raster <- raster("data/NDVI_mean_2000049-2019033.tif") # Folder path for mean NDVI rasters
mean_raster_utm <- projectRaster(from = mean_raster, crs = raster::crs(elevation)) # covert to same projection as everything else
ndvi_mean_matrix <- as.matrix(mean_raster_utm) # create as matrix
ndvi_extent <- extent(mean_raster_utm) # create raster extent object
ndvi_res <- res(mean_raster_utm) # create raster cell resolution object

# note: 
# Calculating aNDVI is optionl in ndvi_extraction
# if not needed specify avg_matrix = NULL (deault option)


# Extract iNDVI, dNDVI and aNDVI at locations
ndvi_results <- spExtract::ndvi_extraction(files_desc = files_desc,
                                           array = ndvi_array,
                                           extent = ndvi_extent,
                                           resolution = ndvi_res,
                                           avg_matrix = ndvi_mean_matrix,
                                           date = as.Date(locations$timestamp),
                                           coords = locations[,c("x","y")])

# Add results to dataframe
locations <- cbind(locations, ndvi_results)

# View results
head(locations, 10)
