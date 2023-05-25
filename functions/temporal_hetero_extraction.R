##################################
##### Function to extract tempral heterogeneity info from a spatialpointsdataframe.
##### Value is the Coefficient of Variation (and SD if wanted) for each metric
##### By Jerod Merkle, 5 June 2015
##################################

##### SOURCE THE FOLLOWING CODE BEFORE RUNNING THE ABOVE CODE
temphetero <- function(XYdata=XYdata, NDVIfolder="E:/StatewideNDVI/ndvi_bisc_calcs/", dates=dates, #as POSIX object 
                       AID=AID,
                       metrics=c("maxIRG","maxIRGdate","springLength","greenupLength", "IRG","DFP"),
                       wantSDtoo=TRUE){ 
  #checks....
  require(raster)
  require(rgdal)
  require(lubridate)
  if(class(XYdata)!="SpatialPointsDataFrame")
    stop("XYdata must be a Spatial Points DataFrame")
  drs <- dir(NDVIfolder)
  if(all(c("maxIRG.grd", "springStart.grd", "springEnd.grd", "maxIRGdate.grd", "maxNDVIdate.grd","dlcparams.grd") %in% drs)==FALSE)
    print('WARNING: you do not have all phenology metrics in NDVIfolder')
  if(is.na(is.projected(XYdata))==TRUE)
    stop("You must specify a projection for your XYdata.")
  
  # math starts here
  if("IRG"%in% metrics){#gives both IRG and abs(days from peak (DFP))
    
    # make sure dates and XYdata match
    XYdata<-XYdata[order(dates), ]
    dates<-sort(dates)
    
    params<-stack(paste(NDVIfolder,"/dlcparams.grd", sep=""))
    XYdata <- spTransform(XYdata,CRSobj = crs(params))
    
    years<-unique(year(dates))
    jdays<-yday(dates)  #converts POSIX into julian day
    IRGvals<-NULL
    DFPvals<-NULL
    
    for(y in 1:length(years)){
      irgMax <- params[[paste("xmidS_", years[y], sep="")]]
      scale <- params[[paste("scalS_", years[y], sep="")]]
      dats.y <- XYdata[year(dates)==years[y], ]
      irgVals <- raster::extract(irgMax, dats.y)  #gives error messages 
      scaleVals <- raster::extract(scale, dats.y)
      
      # run IRG function
      irg <- IRG(irgVals, scaleVals, jdays[year(dates)==years[y]])
      IRGvals <- c(IRGvals, irg)
      dfp <- abs(irgVals-jdays[year(dates)==years[y]])
      DFPvals <- c(DFPvals, dfp)
    }
    
    XYdata$IRG <- IRGvals
    XYdata$DFP <- DFPvals
    
    #reorder XYdata to original order
    XYdata<-XYdata[order(XYdata$AID, ymd_hms(XYdata$date)), ]
  }
  
  
  if("maxIRG" %in% metrics){
    stk <- stack(paste(NDVIfolder, "/maxIRG.grd", sep="")) 
    vals <- extract(stk, spTransform(XYdata, CRS(projection(stk))))
    MEAN <- rowMeans(vals, na.rm=T)
    SD <- apply(vals, 1, sd, na.rm=T)
    XYdata$maxIRG_MEAN <- MEAN
    if(wantSDtoo==TRUE){
      XYdata$maxIRG_SD <- SD
    }
    XYdata$maxIRG_CV <- (SD/MEAN)*100
  }
  
  if("maxIRGdate"%in% metrics){
    stk <- stack(paste(NDVIfolder, "/maxIRGdate.grd", sep="")) 
    vals <- extract(stk, spTransform(XYdata, CRS(projection(stk))))
    MEAN <- rowMeans(vals, na.rm=T)
    SD <- apply(vals, 1, sd, na.rm=T)
    XYdata$maxIRGdate_MEAN <- MEAN
    if(wantSDtoo==TRUE){
      XYdata$maxIRGdate_SD <- SD
    }
    XYdata$maxIRGdate_CV <- (SD/MEAN)*100
  }

  if("springLength"%in% metrics){
    start <- stack(paste(NDVIfolder, "/springStart.grd", sep="")) 
    start <- extract(start, spTransform(XYdata, CRS(projection(start))))
    end <- stack(paste(NDVIfolder, "/springEnd.grd", sep="")) 
    end <- extract(end, spTransform(XYdata, CRS(projection(end))))
    
    length <- end-start
    MEAN <- rowMeans(length, na.rm=T)
    SD <- apply(length, 1, sd, na.rm=T)
    mstart <- rowMeans(start, na.rm=T)
    mend <- rowMeans(end, na.rm=T)
    XYdata$springLength_MEAN <- MEAN
    XYdata$springStart_MEAN <- mstart
    XYdata$springEnd_MEAN <- mend
    if(wantSDtoo==TRUE){
      XYdata$springLength_SD <- SD
    }
    XYdata$springLength_CV <- (SD/MEAN)*100
  }
  
  if("greenupLength"%in% metrics){
    if("springLength"%in% metrics==FALSE){
      start <- stack(paste(NDVIfolder, "/springStart.grd", sep="")) 
      start <- extract(start, spTransform(XYdata, CRS(projection(start))))
    }
    end <- stack(paste(NDVIfolder, "/maxNDVIdate.grd", sep="")) 
    end <- extract(end, spTransform(XYdata, CRS(projection(end))))
    
    length <- end-start
    MEAN <- rowMeans(length, na.rm=T)
    SD <- apply(length, 1, sd, na.rm=T)
    XYdata$greenupLength_MEAN <- MEAN
    if(wantSDtoo==TRUE){
      XYdata$greenupLength_SD <- SD
    }
    XYdata$greenupLength_CV <- (SD/MEAN)*100
  }
  
  return(XYdata)
} #end of function


