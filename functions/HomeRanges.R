library(snowfall)
library(BBMM)
library(raster)

######## FUNCTIONS #######################################
# 1. Subset data monthly, run bbmm function, run rasterize function, create 12 surfaces
# 2. Generates a brownian bridge for each unique individual
# 3. converts bbmm (results) from list into a raster

####
# 1. Subset data monthly, run bbmm function, run rasterize function, create 12 surfaces
beale.monthly<-function(bbmm, # wb telemtry data
                        rasn, # grid on which bbmm is calculated
                        swd   # working directory
){
  
  for(j in 1:12){
    
    #generate bbmm surfaces by month
    bbmm.w<-bbmm[bbmm$SPECIES=="WB" & bbmm$month==j & bbmm$migrant==1,]
    uaid<-unique(bbmm.w$AID)
    results2<-fallsnow(bbmm.x=bbmm.w,rasn=rasn,uaid=uaid) #function above
    
    #rasterize all bbmm's for month j
    all<-rasterize.bbmm(results2,wb,uaid,j,utmproj) #this may not have correct arguments
    
    #exclude cells with very small values (less than 1 wildebeest)
    all[all<0.000001]<-NA
    
    #recalculate so that sum across all squares is 1.0
    tm<-all*(1/sum(getValues(all),na.rm=T))
    writeRaster(tm, file=paste0("bbmm-wb-month_",j,"_",Sys.Date(),".tif"),driver="GEO Tif")
  }
  return(tm)
}

# claculate MCP by individual
mcp.by.id <- function(xy){
  
  # save crs
  proj <- crs(xy)
  
  # loop through possible window frames to find one that encompasses 
  # the largest number of individuals

  x <- as.data.frame(xy) %>% 
    group_by(AID) %>% 
    filter(date <= (min(date) + (364*60*60*24)))

  x <- as.data.frame(x)
  x <- x[,c('x','y','AID')]
  coordinates(x) <- ~x+y
  proj4string(x) <- proj
  x$AID <- as.factor(x$AID)
  
  z <- mcp.area(x,percent = 100,unout = 'km2',plotit = F)  
  
  # w<-kernelUD(xy = x,h = 'href')
  # kernel.area(w,percent = 95,unout = 'km2')
  # zz <- mcp(x,percent = 100)  
  return(z)
}


#####
# 2. Generates a brownian bridge for each unique individual
fallsnow<-function(bbmm.x # wildebeest telemtery data
                   ){

  # save crs
  proj <- crs(bbmm.x)
    
  # remove all but first year of data with at least 180 days of data
  bbmm.x <- as.data.frame(bbmm.x) %>% 
    group_by(AID) %>% 
    filter(date <= (min(date) + (365*60*60*24))) %>%
    filter(length(date) > 180)
  
  if(nrow(bbmm.x)>0){
    bbmm.x <- as.data.frame(bbmm.x)
    coordinates(bbmm.x) <- ~x+y 
    proj4string(bbmm.x) <- proj
    
    # identify unique individuals
    uaid <- unique(bbmm.x$AID)
    
    #initiate snowfall
    sfInit(parallel=T,cpus=3, type="SOCK")
    
    #specify any homemade functions and objects to be read into the multithread environment
    sfExport("bbmm.x","uaid",'proj')
    
    #specify all packages to be used in the multithread environment
    pks <- c("BBMM","raster",'rgdal')
    for(i in 1:length(pks)){
      sfLibrary(pks[i],character.only=TRUE)
    }
    rm(pks)
    
    ### SNOWFALLL - 
    #loop through all individuals i in uaid
    # test<-NULL
    # results<-NULL
    results <- do.call(rbind,sfClusterApplyLB(1:length(uaid), function(i){  #
      
      # sqkmhr95 <- NULL
      # for(i in 1:length(uaid)){
      tmp<-bbmm.x[bbmm.x$AID==uaid[i],]
      tmp<-tmp[!is.na(tmp$date),]
      
      # create raster template
      ext <- extent(tmp)
      dims <- c(ext@ymax-ext@ymin,ext@xmax-ext@xmin)
      maxdim <- dims[which.max(dims)]*0.05
      
      ext@xmin <- ext@xmin - maxdim
      ext@xmax <- ext@xmax + maxdim
      ext@ymin <- ext@ymin - maxdim
      ext@ymax <- ext@ymax + maxdim
      
      # identify which dimension (ht or width) maximizes the landscape
      dims <- c(ext@ymax-ext@ymin,ext@xmax-ext@xmin)
      maxdim <- which.max(dims)
      ress <- c(NA,NA)
      ress[maxdim] <- 200
      ress[-maxdim] <- round(dims[-maxdim]/(dims[maxdim]/200))
      
      # create raster
      rasn <- raster(ext, nrows=ress[1], ncols=ress[2],vals=NULL,crs=proj)
      
      #raster management
      area.grid <- rasterToPoints(rasn)[,c(1,2)]
      
      #make sure sufficient number of points per individual to estimate monthly bbmm
      if (nrow(tmp)>=20){
        x<-tmp$x
        y<-tmp$y
        
        # difference in minutes between each point
        time.lag <- as.numeric(tmp$date[2:length(tmp$date)])/60-as.numeric(tmp$date[1:(length(tmp$date)-1)])/60
        
        maxlagx<- 60*36 #exclude any consecutive points that are greater than 36 hours (in minutes)
        test<-brownian.bridge(x,y,
                              time.lag = time.lag,
                              location.error = 5,
                              area.grid = area.grid,
                              max.lag = maxlagx)
        
        # calc home range size - 99% level
        contours = bbmm.contour(test, levels=95,plot=FALSE) 
        ww <- test$probability
        ww <- ww[ww>=contours$Z]
        ww <- ww[order(ww)]
        sqkmhr95 <- data.frame(hr = (length(ww) * res(rasn)[1] * res(rasn)[2] / (1000^2)),
                               AID = uaid[i])
        return(sqkmhr95)
      }
    }
    ))
    
  }else{
    results <- NULL
  } #end if statement
  return(results)
}

######
# 3. converts bbmm (results) from list into a raster
rasterize.bbmm<-function(results,bbmm.w,uaid,k,utmproj){
  
  all1<-NULL
  
  for(i in 1:length(results)){
    try1<-results[[i]]
    tmp <- data.frame(x=try1$x, y=try1$y, z=try1$probability)
    P<-bbmm.w[bbmm.w$AID==uaid[i] & bbmm.w$month==k & bbmm.w$migrant==1,]
    
    if(class(P)=="data.frame") {stop("BBMM must be projected")}
    
    # rasterize the bbmm's
    if(nrow(P)>=20){
      m <- SpatialPixelsDataFrame(points = tmp[c("x", "y")], data=as.data.frame(tmp$z),
                                  proj4string = CRS(utmproj))
      m <- as(m, "SpatialGridDataFrame")
      d<-as(m,"RasterLayer")
      d[d==0]<-NA
      
      if(is.null(all1)==TRUE){
        all1<-d
      }else{
        all1<-sum(all1,d,na.rm=T)
      }
    }
    print(paste(i," out of ",length(results)," - ",uaid[i],sep=""))
  }
  return(all1)  
}

####### RUN CODE ###########

# projection
# utmproj<-("+proj=utm +south +zone=36 +init=epsg:21036") #spatial coordinate reference system in Serengeti
# 
# # read grid for Brownian Bridges
# ras<-raster("NDVI_MOD13Q1.A2000049.h21v09.006.tif")
# ras<-projectRaster(ras,crs = utmproj)
# ras <- aggregate(ras, 6) #lower the resolution of grid by 2 times
# 
# # read and manage wildebeest trajectory data
# wb<-read.csv("WB_ZB_SNP_DATA_19990225-20190610_Deads_deleted.csv")
# wb$Date <- as.POSIXct(as.character(wb$Date), "%Y-%m-%d %H:%M:%S",tz="UTC")
# wb<-wb[order(wb$AID,wb$Date),]
# attributes(wb$Date)$tzone<-"Africa/Nairobi" #change from UTC to Nairobi timezone
# wb <- wb[wb$OWNER=="HOPCRAFT" & wb$SPECIES=="WB" & wb$year>2007 & wb$migrant==1,] #subset to be only Grant's wb data
# 
# # remove individuals with fewer than 2 weeks data
# dd<-aggregate((wb$Date),by=list(wb$AID),FUN=function(x) max(x)-min(x))
# dd<-dd[dd[,2]<15,1]
# wb<-wb[!wb$AID %in% dd,]
# 
# # make wb into spatial object
# coordinates(wb)<- ~x+y
# proj4string(wb)<-CRS(utmproj)
# 
# # define working directory
# swd<-"C:/Users/tmorri19/OneDrive - University of Glasgow/Documents/Glasgow/_STUDENTS/Freja/"
# 
# # run bbmm functions - this will save monthly bbmm's as a raster to your working directory (swd) 
# beale.monthly(
#   bbmm = wb,
#   rasn = ras,
#   swd)