# Convert raw MODIS files to geotiffs, clip scenes, 
# Extract NDVI values at spatial points

###### Convert raw HDF files to geotiffs, clip to size of Serengeti extent

convert.hdf.geotiff <- function(folder,    # file location where hdf/hdr files are stored
                                newfolder  # file location where geotiffs are stored
){
  library(sp)
  library(lubridate)
  library(raster) 
  library(rgdal)
  library(rgeos)
  library(gdalUtils)
  library(snowfall) #parallel processing
  
  z<-list.files(path = folder,pattern=".hdf") #list all HDF files in your directory
  
  if(!missing(newfolder)) setwd(newfolder)
  
  # loop across all hdf files and extract NDVI and EVI
  for(i in 1:length(z)){
    gdalinfo(z[i])
    
    # Which subdatasets are within my hdf4 MODIS files and makes them into a list
    # sd[[1]] == NDVI
    # sd[[2]] == EVI
    sds <- get_subdatasets(z[i])
    
    # Only select first subdataset = NDVI and convert+write to a geotiff
    gdal_translate(sds[1], dst_dataset = paste0("NDVI_",substr(z[i],1,41),".tif")) #ndvi
    gdal_translate(sds[2], dst_dataset = paste0("EVI_",substr(z[i],1,41),".tif")) #evi
    # now you can read the new object as a geotif if youd like
  }
}

clip.ndvi <- function(folder,      #folder location where whole scenes are saved
                      clipfolder,  #folder location where clipped images will be saved
                      newdts,      #vector of dates in MODIS scene to be clipped
                      #Extent of clip for Serengeti 
                      xmin=33.66042, 
                      xmax=36.26042,
                      ymin=-4.031250,
                      ymax=-0.2708333
){
  
  # clip it to Serengeti scene
  latlong <- ("+proj=longlat +datum=WGS84")
  modcrs<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
  
  setwd(folder)
  z<-list.files(pattern="NDVI_") #list all HDF files in your directory
  if(!missing(newdts)) z<-z[substr(z,15,21) %in% newdts] #subset list by selected dates
  
  # extent used for Serengeti NDVI scenes
  ext<-extent(c(xmin,xmax,ymin,ymax)) #define extent of study area to be clipped
  extsp<-as(ext,'SpatialPolygons')
  crs(extsp)<-latlong
  extsp<-spTransform(extsp,CRSobj = crs(modcrs))
  
  for(i in 1:length(z)){
    mod<-raster(z[i])
    mod<-crop(mod,extsp)
    mod<-projectRaster(mod,crs=latlong) #reproject to latlong
    
    if(maxValue(mod)>1000) mod<-mod/10000000   
    
    # Write to file  
    writeRaster(mod,paste0(
      clipfolder,
      substr(z[i],1,32),".tif"),
      format="GTiff",overwrite=T)
  }
}

####### Functions to extract ndvi data

build.ndvi.description <- function(folder, # full path location of folder containing NDVI geotifs
                                   date.loc=c(15,21), # character location in filenames of NDVI containing dates 
                                   date.format="%Y%j", # format of dates
                                   start.date = NULL, end.date = NULL){ 
  library(sp)
  library(lubridate)
  library(raster) 
  library(rgdal)
  library(rgeos)
  library(snowfall)
  
  ndvi.data<-data.frame(filenames=list.files(folder,pattern = "MOD13Q1"),
                        year=NA,date=ymd("1970-01-01"))
  
  ndvi.data$filenames <- as.character(ndvi.data$filenames)
  
  ndvi.data$date<-as.Date(substr(ndvi.data$filenames, date.loc[1], date.loc[2]),
                          format=date.format,
                          tz="GMT")
  
  ndvi.data$year<-year(ndvi.data$date)
  
  ndvi.data$folder<-folder
  
  ## order by date
  ndvi.data<-ndvi.data[order(ndvi.data$date),]
  ndvi.data[,1]<-as.character(ndvi.data[,1])
  
  ## filter by start and end date if provided
  if (!is.null(start.date)){
    start.date <- as.Date(start.date)
    ndvi.data <- ndvi.data[ndvi.data$date >= start.date,]
  }
  
  if (!is.null(end.date)){
    end.date <- as.Date(end.date)
    ndvi.data <- ndvi.data[ndvi.data$date <= end.date,]
  }
  
  return(ndvi.data)
}

build.ndvi.stack <- function(ndvi.description){
  ndvi.stack <- stack()
  
  for(i in 1:nrow(ndvi.description)){
    ndvi.stack <- stack(ndvi.stack, raster(dsn='.', paste(ndvi.description$folder[i],
                                                          ndvi.description$filenames[i],
                                                          sep="/")))
  }
  
  return(ndvi.stack)
}

extract.ndvi <- function(folder,         # folder in which NDVI_mean.tif is saved
                         stack,          # Raster stack output from above build.ndvi.stack
                         stack.desc,     # Stack description from above build.ndvi.description
                         loc,            # All x-y locations projected in latlong
                         date,           # as.Dates from loc, e.g. date column from your GPS points
                         cpus=3,         # number of cpus to use to run extract 
                         extract.andvi=T # extract anomally NDVI?
                         ){
  # Warning messages
  if (length(date) > 1 && length(date) != nrow(loc))
    stop("Dates must be either one date for all points, or one date for each point")
  if (class(loc) != "SpatialPointsDataFrame")
    stop("Locations must be a SpatialPointsDataFrame with same projection as raster stack")    
  if (substr(as.character(crs(loc)),1,26) != substr(as.character(crs(stack)),1,26))
    stop("Locations and raster stack do not have same CRS projection")    
  
  # Check if extraction is across multiple dates
  if (length(date) == 1) {multi.extract <- T}else{multi.extract <- F}
  
  # ndvi.ret <- matrix(rep(rep(NA, 3), nrow(loc)), nrow=nrow(loc))
  
  # calculate mean NDVI across all years
  ave.ndvi <- raster(paste0(folder,"/NDVI_mean_2000049-2019033.tif"))

  ## Loop across each date in the dates provided

  #initiate snowfall
  sfInit(parallel=T,cpus=cpus, type="SOCK")   
  
  #specify any homemade functions and objects to be read into the multithread environment
  sfExport("ave.ndvi","date","stack.desc","stack","loc","multi.extract","extract.andvi")
  
  #specify all packages to be used in the multithread environment
  pks <- c("raster")
  for(i in 1:length(pks)){
    sfLibrary(pks[i],character.only=TRUE)
  }
  rm(pks)
  
  ### SNOWFALLL - parallel process 
  #loop through all individuals i in uaid

  results <- do.call(rbind,sfClusterApplyLB(1:length(date), function(i){  #
    
  # for (i in 1:length(date)){
    
    nidx<-which.min(abs(as.numeric(stack.desc$date-(date[i]))))
    
    ## initialize the before and after layer indices to the nearest layer
    loidx<-nidx
    hiidx<-nidx
    
    ## now calculate whether we want to use the previous or next raster for interpolation;
    ## switch is faster for performing these kind of ifelse operations, so scale and normalize
    diff<-as.numeric(stack.desc$date[nidx]) - as.numeric(date[i])      # -n, 0, n
    if (diff != 0) diff<-diff/abs(diff)
    diff <- diff + 2                                                   #  1, 2, 3
    
    ## 1 means interpolate forward, 3 means interpolate backwards; 2 means we are on an NDVI
    ## layer date, so no interpolation needed
    switch(diff, hiidx<-nidx+1, interp<-1, loidx<-nidx-1)
    
    ## get the scaling value for interpolation if needed; check if we are interpolating
    if (hiidx != loidx){
      ## calculate the linear scaling value from the previous to the next NDVI values
      interp<-(as.numeric(date[i]) - as.numeric(stack.desc$date[loidx])) /
        (as.numeric(stack.desc$date[hiidx]-stack.desc$date[loidx]))
    } else {
      
      ## we are bang on an NDVI date; check how we need to adjust where we are interpolating to:
      ## full interpolation from the previous NDVI date
      if (loidx > 1){
        loidx <- loidx - 1
        ## no interpolation from the current NDVI date
      } else {
        hiidx <- hiidx + 1
        interp <- 0
      }
    }
    
    ## check if we actually have the later layer, ignore if not
    if (hiidx <= nrow(stack.desc)){
      
      ## load both layers (previous and next)
      loras<-raster(stack, loidx)
      hiras<-raster(stack, hiidx)
      
      ## get the before and after NDVI values
      if (multi.extract){
        ndvi_start<-extract(loras,loc) 
        ndvi_end<-extract(hiras,loc)
      } else {
        ndvi_start<-extract(loras,loc[i,]) 
        ndvi_end<-extract(hiras,loc[i,])        
      }
            
      ## calculate NDVI changes, and then interpolate the estimated NDVI value
        indvi <- ndvi_start + (interp * (ndvi_end - ndvi_start)) #NDVI, interpolated
        dndvi <- (ndvi_end - ndvi_start) / as.numeric(stack.desc$date[hiidx] - stack.desc$date[loidx], units = "days") #delta ndvi

      ## calculate days different from max IRG for this time period     
        dtrange <- (loidx-30):(hiidx+30)
        fs <- extract(stack[[dtrange]],loc[i,]) # time consuming step
        
        # 
        # require(zoo)
        # Sys.time()<-x
        # fs <- extract(stack,loc[i,]) # time consuming step
        # x-Sys.time()
        
        fs.na <- c(fs, NA)
        #calcuate rolling mean of 3 values at a time
        fs.na <- rollapply(fs.na, width = 5, by = 1, FUN = mean, na.rm = TRUE, align = "left")
        
        plot(stack.desc$date[1:length(fs)],fs,type='l',lty=2,col='grey',lwd=1,xlab='',ylab='NDVI')
        lines(stack.desc$date[1:length(fs.na)],fs.na,type='l',col='red',lwd=1.7)
        critn <- 6
        # identify the set of dates where IRG is maximal
        xx <- (sapply(1:length(fs.na), function(x) {
          if(x>critn & x<(length(fs.na)-critn)){
            p0 <- fs.na[x]-fs.na[x-1]
            p1 <- sapply(1:critn,function(z) fs.na[x+z]-fs.na[x+(z-1)])
            if(p0<0 & all(p1>0)) {
              green=1
            }else{
              green=0
            }
            k0 <- fs.na[x+1]-fs.na[x]
            k1 <- sapply(1:critn,function(z) fs.na[x-(z-1)]-fs.na[x-z])
            if(k0<0 & all(k1>0)) {
              green=2
            }
          }else{
            green=0
          }
          return(green)
        }))
        
        # wx <- which(x > qnorm(0.95,mean(x),sd(x)))
        startx <- which(xx == 1)
        endx <- which(xx == 2)
        points(stack.desc$date[startx],fs.na[startx],col='darkred',pch=16)
        points(stack.desc$date[endx],fs.na[endx],col='blue',pch=16)
        mtext(paste(as.data.frame(loc[i,])[,2],as.data.frame(loc[i,])[,3]))
        stack.desc$date[wx]
        
        dfs <- which.max(sapply(2:length(fs), function(x) (fs[x]-fs[x-1]))) #1:nb
        maxirg <- median(stack.desc$date[dtrange[c(dfs,dfs+1)]])
        daysdiffirg <- as.numeric(difftime(date[i],maxirg,units = 'days'))
        
      ## set the relevant NDVI and dNDVI values
      if (multi.extract){
        ndvi.ret <- data.frame(indvi=indvi,dndvi=dndvi,daysdiffirg=daysdiffirg)
      } else {
        ndvi.ret <- data.frame(indvi=indvi,dndvi=dndvi,daysdiffirg=daysdiffirg)
      }
    }
    return(ndvi.ret)
  }))
  sfStop()
  
  if(extract.andvi){
    ave.ndvi.val<-extract(ave.ndvi,loc)
    andvi <- (results$indvi-ave.ndvi.val)/ave.ndvi.val # anomaly NDVI
    results$andvi <- andvi
  }
  return(results)
}

## EXAMPLE RUN CODE


# # Load wb GPS data
# wb <- read.csv("C:/Users/tmorri19/OneDrive - University of Glasgow/Serengeti Wildebeest Project/Tracking Database/WB_ZB_SNP_DATA_19990225-20190610_Deads_deleted.csv")
# wb<-as.data.frame(wb)
# 
# # Subset wb GPS data for dates in which 
# wb2 <- wb[as.Date(wb$Date)<min(sd$date) | as.Date(wb$Date)>max(sd$date),]
# wb3 <- wb[as.Date(wb$Date)>=min(sd$date) & as.Date(wb$Date)<=max(sd$date),]
# loc<-wb3[,names(wb3) %in% c("x","y","Date")]
# 
# # Convert locs to spatial object
# coordinates(loc) <- ~x+y
# utmproj<-("+proj=utm +south +zone=36 +init=EPSG:21036") #spatial coordinate reference system in Serengeti
# latlongproj<-("+proj=longlat +datum=WGS84")
# proj4string(loc)<-utmproj
# loc<-spTransform(loc,CRSobj = latlongproj)
# date<-as.Date(loc$Date)
 
# run NDVI, dNDVI and aNDVI extraction for all points within date range
# {w<-Sys.time()
#   xx<-extract.ndvi(folder="C:/Users/Administrator/OneDrive - University of Glasgow/Serengeti Wildebeest Project/SerengetiGIS/NDVI_mean",      # folder in which NDVI_mean.tif is saved
#                    stack=st,       # Raster stack output from above build.ndvi.stack
#                    stack.desc=sd,  # Stack description from above build.ndvi.description
#                    loc=loc,         # All x-y locations projected in latlong
#                    date=date,        # as.Date from loc, e.g. date column from your GPS points
#                    cpus=3,
#                    extract.andvi=T)
# Sys.time()-w}

# combine data with WB GPS datasets
# wb3 <- cbind(wb3,xx)
# wb2$indvi <- NA
# wb2$dndvi <- NA
# wb2$andvi <- NA 
# wb3 <- rbind(wb3,wb2)
# write.csv(wb3,"WB_ZB_SNP_DATA_19990225-20190610_Deads_deleted.csv",row.names = F)
