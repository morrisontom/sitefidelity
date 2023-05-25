

extract_Cspace <- function(dat,dat2,modfold){  

  # find min and max dates of xy dataset +/- 16 days  
  mindt <- min(as.Date(dat$date),na.rm=T)-16
  maxdt <- max(as.Date(dat$date),na.rm=T)+16
  
  ### read in NDVI layers that correspond with data
  lf <- list.files(path=modfold,pattern='.tif')  
  if(!dat2$title[1] %in% c("WB_SE","ZB_SE","WB_AM") ){
    ndvi.dates <- as.Date(substr(lf,10,17),format='%Y%j')
  }else{
    ndvi.dates <- as.Date(substr(lf,15,21),format='%Y%j')
  }
  
  # exclude amboseli wb because of missing ndvi data
  if(!dat2$title[1]=='WB_AM'){
    
    # subset to include only dates that are within time period of the data
    lf <- lf[ndvi.dates >= mindt & ndvi.dates <= maxdt ]
    ndvi.dates <- ndvi.dates[ndvi.dates >= mindt & ndvi.dates <= maxdt ]
    
    # create stack of NDVI 
    maps <- stack(paste0(modfold,'/',lf))
    
    # transform xy data if they are in different CRS
    dat <- spTransform(dat,CRSobj = crs(maps))
    
    # crop to all individuals to reduce next 2 computations
    # maps <- crop(maps,dat)
    
    # recalc if NDVI data are not scaled correctly
    # if(any(maxValue(maps,layer=1)>1000)) maps <- calc(maps,fun = function(x) {x/100000000})
    
    if(!substr(dat2$title[1],1,2) %in% c('RAND')){
      
      #initiate snowfall
      sfInit(parallel=T,cpus=3, type="SOCK")
      
      #specify any homemade functions and objects to be read into the multithread environment
      sfExport('maps', 'ndvi.dates', 'dat', 'dat2')
      
      #specify all packages to be used in the multithread environment
      pks <- c("raster",'rgdal','graphics')
      for(i in 1:length(pks)){
        sfLibrary(pks[i],character.only=TRUE)
      }
      rm(pks)
      
      ### SNOWFALLL multicore
      
      # loop through all individuals i in uaid
      df2 <- do.call(rbind,sfClusterApplyLB(1:nrow(dat), function(i){  #
        wm <- abs(ndvi.dates-as.Date(dat$date[i]))
        
        ndvi <- maps[[which.min(wm)]]
        dat.aid <- dat[dat$AID==dat$AID[i],]
        
        # crop NDVI to extent for each individual - individuals with 1 pt are buffered first
        if(nrow(dat.aid)>1) {
          ndvi <- crop(ndvi,dat.aid)
        }else{
          dbuff <- gBuffer(dat.aid[dat.aid$AID==dat.aid$AID[i],],byid = T,width = 0.01)
          ndvi <- crop(ndvi,dbuff)
          print(paste0('WARNING:',dat$AID[i],'Cpsace value is incorrect'))
        }
        
        # remove NA's, rescale to 0-1
        Vndvi = getValues(ndvi)
        
        # set to NA any water related values
        Vndvi[Vndvi<0] <- NA
        
        # only Cspace for dates with corresponding iyd val + where NDVI layer is available within 8 days 
        if(any(!is.na(dat2[i,c("iydREAL","iyd0","iyd20","iyd40","iyd60","iyd80",
                               "iyd100","iyd120","iyd140","iyd160","iyd180")])) & 
           min(wm)<=8 & any(!is.na(Vndvi))){
          
          # rescale so NDVI is between 0 and 1 - caused by a NASA issue with the donload data
          if(max(Vndvi,na.rm=T)>10000){
            Vndvi <- Vndvi/100000000
          }
          
          #------------------------------------------------------------------------------------------------
          ## define the NDVI intervals among which the NDVI data will be broken down.
          #-----------------------------------------------------------------------------------------------
          # ranges = cellStats(ndvi, "range")
          n = 100                                                      # number of classes
          
          # intervals = seq(from = min(ranges[1,]), to = max(ranges[2,]), length.out = n+1)
          # mval <- minValue(ndvi)
          # mxval <- maxValue(ndvi)
          
          mval <- min(Vndvi,na.rm=T)
          mxval <- max(Vndvi,na.rm=T)
          intervals = seq(from = mval, to = mxval, length.out = n+1)
          
          #------------------------------------------------------------------------------------------------------
          ##-----------------------------------------------------------------------------------------------------
          ### Constancy in space
          ##-----------------------------------------------------------------------------------------------------
          #------------------------------------------------------------------------------------------------------
          # This is a measure of the environment's spatial homogeneity (at a given time) 
          
          #--------------------------------------------------------------------------------
          ## 1st possible index: Colwell's constancy in space
          # also called compositional entropy  (see Vranken et al. 2015). 
          # It is equal to 1 if the environment is homogeneous, and is very small otherwise.
          
          # original code
          # eff = sapply(1:(length(intervals) - 1)
          #              , function(x) sum((Vndvi > intervals[x]) & (Vndvi <= intervals[x+1]),na.rm=TRUE))
          
          # new code - use hist instead of sapply - 12x faster      
          ww <- Vndvi[!is.na(Vndvi)]
          g <- graphics::hist(ww,plot=FALSE,breaks=intervals)
          eff <- g$counts
          
          # sum over all elements in interval
          p = eff/sum(eff)
          Cspace1 = 1 + sum(p[p > 0] * log(p[p > 0]))/log(n)
          
          #-------------------------------------------------------------------------------
          ## 2nd possible index: inverse of standard deviation
          
          Cspace2 = 1/sd(Vndvi, na.rm = TRUE)
        }else{
          Cspace1 = NA
          Cspace2 = NA
        }
        # add vals to dataframe
        anom <- data.frame(
          Cspace1_entropy = Cspace1,
          Cspace2_1sd = Cspace2)  
        
        return(anom)
      })) # end snowfall
      
      # if dataset is CA, WB or ZE just loop
    }else{  
      
      # to run snowfall need: maps, ndvi.dates, dat, dat2
      for(i in 1:nrow(dat)){
        wm <- abs(ndvi.dates-as.Date(dat$date[i]))
        
        ndvi <- maps[[which.min(wm)]]
        dat.aid <- dat[dat$AID==dat$AID[i],]
        
        # crop to reduce computations
        if(nrow(dat.aid)>1) {
          ndvi <- crop(ndvi,dat.aid)
        }else{
          dbuff <- gBuffer(dat.aid[dat.aid$AID==dat.aid$AID[i],],byid = T,width = 0.01)
          ndvi <- crop(ndvi,dbuff)
          print(paste0('WARNING:',dat$AID[i],'Cpsace value is incorrect'))
        }
        # remove NA's, rescale to 0-1
        Vndvi = getValues(ndvi)
        Vndvi[Vndvi<0] <- NA
        
        # only Cspace for dates with corresponding iyd val + where NDVI layer is available within 8 days 
        if(any(!is.na(dat2[i,c("iydREAL","iyd0","iyd20","iyd40","iyd60","iyd80",
                               "iyd100","iyd120","iyd140","iyd160","iyd180")])) & 
           min(wm)<=8 & any(!is.na(Vndvi))){
          
          if(max(Vndvi,na.rm=T)>10000){
            Vndvi <- Vndvi/100000000
          }
          #------------------------------------------------------------------------------------------------
          ## define the NDVI intervals among which the NDVI data will be broken down.
          #-----------------------------------------------------------------------------------------------
          # ranges = cellStats(ndvi, "range")
          n = 100                                                      # number of classes
          # intervals = seq(from = min(ranges[1,]), to = max(ranges[2,]), length.out = n+1)
          # mval <- minValue(ndvi)
          # mxval <- maxValue(ndvi)
          
          mval <- min(Vndvi,na.rm=T)
          mxval <- max(Vndvi,na.rm=T)
          intervals = seq(from = mval, to = mxval, length.out = n+1)
          
          #------------------------------------------------------------------------------------------------------
          ##-----------------------------------------------------------------------------------------------------
          ### Constancy in space
          ##-----------------------------------------------------------------------------------------------------
          #------------------------------------------------------------------------------------------------------
          # This is a measure of the environment's spatial homogeneity (at a given time) 
          
          #--------------------------------------------------------------------------------
          ## 1st possible index: Colwell's constancy in space
          # also called compositional entropy  (see Vranken et al. 2015). 
          # It is equal to 1 if the environment is homogeneous, and is very small otherwise.
          
          # original code
          # eff = sapply(1:(length(intervals) - 1)
          #              , function(x) sum((Vndvi > intervals[x]) & (Vndvi <= intervals[x+1]),na.rm=TRUE))
          
          # new code - use hist instead of sapply - 12x faster      
          ww <- Vndvi[!is.na(Vndvi)]
          g <- hist(ww,plot=FALSE,breaks=intervals)
          eff <- g$counts
          
          # sum over all elements in interval
          p = eff/sum(eff)
          if(i==1){
            Cspace1 = 1 + sum(p[p > 0] * log(p[p > 0]))/log(n)
          }else{
            Cspace1 = c(Cspace1,
                        1 + sum(p[p > 0] * log(p[p > 0]))/log(n)
            )
          }
          
          #-------------------------------------------------------------------------------
          ## 2nd possible index: inverse of standard deviation
          if(i==1){
            Cspace2 = 1/sd(Vndvi, na.rm = TRUE)
          }else{
            Cspace2 = c(Cspace2,
                        1/sd(Vndvi, na.rm = TRUE)
            )
          }
        }else{
          if(i==1){
            Cspace1 = NA
            Cspace2 = NA
          }else{
            Cspace1 = c(Cspace1,
                        NA)
            Cspace2 = c(Cspace2,
                        NA)
          }
        }     
      }
      
      # add vals to dataframe
      df2 <- data.frame(
        Cspace1_entropy = Cspace1,
        Cspace2_1sd = Cspace2)  
    }# end if-then WY Species
    
  }else{
    df2 <- data.frame(
      Cspace1_entropy = NA,
      Cspace2_1sd = NA)  
  }#end if title==WB_AM exclude
  return(df2)
}  
