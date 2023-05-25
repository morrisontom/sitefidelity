##################################
##### Function to extract tempral heterogeneity info from a spatialpointsdataframe.
##### Value is the Coefficient of Variation (and SD if wanted) for each metric
##### By Jerod Merkle, 5 June 2015
##################################

##### SOURCE THE FOLLOWING CODE BEFORE RUNNING THE ABOVE CODE
tempheteroIntegNDVI <- function(XYdata=XYdata,juldate=NA,year=NA, NDVIfolder="E:/ndvi_bisc_calcs", 
                       wantSDtoo=TRUE){ 
  #checks....
  require(raster)
  require(rgdal)
  if(juldate > 365)
    stop("juldate must be less than 365")
  if(class(XYdata)!="SpatialPointsDataFrame")
    stop("XYdata must be a Spatial Points DataFrame")
  drs <- dir(NDVIfolder)
  if(is.na(is.projected(XYdata))==TRUE)
    stop("You must specify a projection for your XYdata.")
  
  #math starts here
    stk <- stack(paste(NDVIfolder, "/csumNDVI_", year, ".grd", sep="")) 
    vals <- extract(stk[[juldate]], spTransform(XYdata, CRS(projection(stk))))
    MEAN <- rowMeans(vals, na.rm=T)
    SD <- apply(vals, 1, sd, na.rm=T)
    XYdata$IntegNDVI_MEAN <- MEAN
    if(wantSDtoo==TRUE){
      XYdata$IntegNDVI_SD <- SD
    }
    XYdata$IntegNDVI_CV <- (SD/MEAN)*100
    
  return(XYdata)
} #end of function


