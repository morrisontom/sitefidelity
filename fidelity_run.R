### Code to compile fidelity dataset from a set of GPS telemetry datasets
### T. Morrison 
rm(list=ls(all=TRUE))
setwd('C:/Users/Administrator/OneDrive - University of Glasgow/Documents/Side projects/Fidelity/ungulatefidelity')

####### INPUT PARAMETERS
wind.x <- seq(0,180,20)              ## 84 days orig windwo size no. days within each window for comparison !should be even number!
crit <- 365+round(min(wind.x)/2)     ## no. days required for each individual to be included in analysis
cps <- 5                             ## min number of individuals required to have a single population
nshifts <- 250                     ## number of random trajectories in resource reliability hypothesis
elbuff <- 2500                     ## min distance (in meters) in which two individuals must occur on same day to be included in edgelist function
fold <- paste0("output_",Sys.Date())

# projections
projlatlong<-("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") #lat long projection
aeaproj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

spp<-c("BS","CA","EK","MD","MS","PH","WB","ZB")

# seasonal data
dtrut<-data.frame(spp=spp,
                  rut=c(318,296,273,333,272,277,152,15),
                  parturition=c(144,151,151,169,142,158,46,15),
                  sum1=c(rep(182,length(spp)-2),106,106),       #JUL 1 - TEMPERATE; APR 1 - TROPICS
                  sum2=c(rep(246,length(spp)-2),166,166),       #AUG 31 - TEMPERATE; MAY 31 - TROPICS
                  win1=c(rep(1,length(spp)-2),213,213),     #JAN 1 - TEMPERATE; AUG 1 - TROPICS
                  win2=c(rep(60,length(spp)-2),304,304),    #FEB 28 - TEMPERATE; OCT 31 - TROPICS
                  spr1=c(95,105,88,82,110,83,305,305), # based on modelled NDVI data across all points within individuals 
                  spr2=c(171,165,171,160,170,151,334,334), # (cont...) mean by individuals, then by species 
                  window = c(72,24,88,44,84,72,86,110) ) # window sizes based on SD*2 from parameters of NSD model

# datasets
datasets <- list.files(path = './data')

# AID's to skip in NSD analysis
pp <- read.csv("./NSD/toSkip_NSD.csv")[,1]
# pp <- pp[!duplicated(pp)]
# write.csv(pp,"./NSD/toSkip_NSD.csv",row.names = F)

####### LIBRARY
libs <- c("moveHMM","foreign","lme4","maptools","graphics","stats","cluster",
          "raster","MASS","igraph","CircStats","sp","rgdal","adehabitatHR",
          "adehabitatLT","geoR","stringr","rgeos","snowfall","lubridate",
          "spatstat",'amt',
          'migrateR',  # run in R V 3
          'tidyverse') # run in R v 4

sapply(libs,library,character=TRUE)

####### FUNCTIONS
funks<-c(
  'Cspace.r',
  'ndvi.r',
  'HomeRanges.R',
  "iyd.r",
  'edgelist.r',
  'mcp_all.r',
  'ovlap.r',
  "make_clust.R",
  'temporal_hetero_extraction.r',
  'temporal_NDVIintegrated.r',
  'temporal_NDVIintegrated.r',
  "season.R",
  "IRGfunction.R",
  "rm_crit.R",
  'smp.R',
  "torus_shift.R",
  "find.encountr.R",
  'run.moveHMM.R',
  'wind.size.NSD.R')

for(i in 1:length(funks)){
  source(paste0('./functions/',funks[i])) 
}

# datasets <- datasets[c(1:4,22:25)]

# RUN DATA
runtime <- Sys.time()
for(q in 1:length(datasets)){
  
  # define dataset-specific parameters
  title<-datasets[q]
  
  # read/manage data 
  dat<-read.csv(paste0('./data/',title,'/',title,'.csv'),header=TRUE)
  print(paste("yo, this be dataset:", datasets[q],q))
  spp <- substr(title,1,2)
  wind <- c(dtrut$window[dtrut$spp==spp],wind.x)
  
  region <- 'noregion'
  if(!fold %in% list.files(path = './output/') ){ dir.create(paste0('./output/',fold))}

  # format data
  names(dat)[grep(pattern ='AID',names(dat))] <-"AID"
  if(length(names(dat)[names(dat)=="TelemDate"])>0) names(dat)[names(dat)=="TelemDate"] <- "TelemDt"
  dat$date <- as.POSIXct(strptime(as.character(dat$TelemDt), format = "%Y-%m-%d %H:%M:%S", tz ="GMT"))
  dat<-dat[,!names(dat) %in% c("SEX","SQL_Shp","Elvtn_m","Date","Time","TmLg_Al","TelemDt","Status")]
  dat<-dat[!is.na(strptime(dat$date,"%Y-%m-%d %H:%M:%S")),]
  if(length(names(dat)[names(dat)=="AEA_X"])>0) names(dat)[names(dat)=="AEA_X"] <- "x"
  if(length(names(dat)[names(dat)=="AEA_Y"])>0) names(dat)[names(dat)=="AEA_Y"] <- "y"
  if(length(names(dat)[names(dat)=="x"])==0|length(names(dat)[names(dat)=="y"])==0 
     |length(names(dat)[names(dat)=="AID"])==0|length(names(dat)[names(dat)=="date"])==0) 
    stop("rename columns in input data")  
  
  #remove any missing dates and locs
  dat <- dat[!is.na(dat$date),]
  dat <- dat[!is.na(dat$x),]
  dat <- dat[!is.na(dat$y),]
  
  # region specific stuff
  if(!title %in% c("WB_SE","ZB_SE","WB_AM","CA_AL")) {
    dat<-dat[dat$x>470000,]
    attributes(dat$date)$tz <- 'MST'
    # dat$date <- with_tz(dat$date,'MST')
    uproj <- "+proj=utm +north +zone=12 +datum=WGS84"
    modfold <- 'C:/Users/tm73z/Documents/Predictability/studysites/Wyoming/NDVI'
    if(region!='WY'){ # used here so that stack does not have to be read multi times
      pred <- stack('./predictability/Ctime_Wyoming.tif')
      region <- 'WY'
    }
  }
  if(title %in% c("WB_SE","ZB_SE")) {
    dat$date <- with_tz(dat$date, 'Africa/Nairobi')
    uproj <- "+proj=utm +south +zone=36 +datum=WGS84"
    modfold <- 'C:/Users/tm73z/OneDrive - University of Glasgow/Serengeti Wildebeest Project/SerengetiGIS/NDVI/'
    if(region!='SE'){
      pred <- stack('./predictability/Ctime_WB_SE.tif')
      region <- 'SE'
    }
  }
  if(title %in% c("WB_AM")) {
    dat$date <- with_tz(dat$date, 'Africa/Nairobi')
    uproj <- "+proj=utm +south +zone=37 +datum=WGS84"
    modfold <- 'C:/Users/tm73z/Documents/Predictability/studysites/Mara-Amboseli-Athi-Kapeitei/NDVI/'
    pred <- stack('./predictability/Ctime_Mara-Amboseli-Athi-Kapeitei.tif')
  }
  if(title %in% c("CA_AL")) {
    dat$date <- with_tz(dat$date, 'America/Anchorage')
    uproj <- "+proj=utm +north +zone=6 +datum=WGS84"
    modfold <- 'C:/Users/tm73z/Documents/Predictability/studysites/Alaska/Alaska_NDVI/'
    pred <- stack('./predictability/Ctime_Alaska.tif')
  }

  # remove any duplicated date/times
  dat <- dat[!duplicated(paste0(dat$date,dat$AID)),]
  
  # remove unwanted columns
  dat <- dat[,names(dat) %in% c('AID','x','y','date')]
  
  # Compile individuals with >crit days of data
  dat2 <- rm_crit(dat,crit)
  
  # reorder
  dat2 <- dat2[order(dat2$AID,dat2$date),]
  
  # add dataset name
  dat2$title <- title
  dat2$AID <- paste0(dat2$title,"_",dat2$AID)
    
  # subsample data to 1 fix per day
  dat2 <- smp(dat2)
  
 # x <- smp(dat2)
 # x[x$AID=='EK_AB_8',]
 # dat2[dat2$AID=='EK_AB_8',]
 
  # classify individuals into populations - unneccessary if random group ID is study ID
  # if(!is.null(dat2)){
  #   dat2el<-smp_edgelist(dat2)
  #   el<-net(dat2el,elbuff,title,getwd())
  #   names(el) <- c("from","to","freq")
  #   clsts <- make_clust(data=el, 
  #                       freq_cutoff=1, 
  #                       use_weights=TRUE, 
  #                       name=title,
  #                       plot_it=FALSE)
  #   
  #   # remove clusters w/o pop ID
  #   pp<-match(dat2$AID,clsts$AID)
  #   dat2$pop<-clsts$cluster[pp]
  #   dat2<-dat2[!is.na(dat2$pop),]
  #   
  #   # identify clusters with fewer than crit number of individuals
  #   tmp<-NULL
  #   for(i in unique(clsts$cluster)){
  #     if(table(clsts$cluster)[[i]]>=cps) tmp<-rbind(tmp,clsts[clsts$cluster==i,])
  #   }
  #   
  #   # add exclude=1 to column in data to know that it should be excluded later
  #   dat2$exclude <- 1
  #   dat2$exclude[dat2$AID %in% tmp$cluster] <- 0
  #   
  # }
  
  # add lat/long
  latlong <- dat2[,c("x","y")]
  coordinates(latlong) <- ~x+y
  proj4string(latlong) <- uproj
  latlong <- as.data.frame(spTransform(latlong,projlatlong))
  dat2$longitude <- latlong$x
  dat2$latitude <- latlong$y
  rm(latlong)
  
  # convert to spatial data
  xy <- dat2[,c('AID','date','x','y')]
  coordinates(xy) <- ~x + y
  proj4string(xy) <- uproj
  xy <- spTransform(xy,CRSobj = crs(pred))
  date <- as.Date(dat2$date)
  
  # # build ndvi datasets
  # sd<-build.ndvi.description(folder = modfold)
  # sd <- sd[sd$date > (min(date)-102) & sd$date < (max(date)+102),]
  # st<-build.ndvi.stack(sd)
  # 
  # # build features to go into ndvi calcs
  # loc <- spTransform(xy,CRSobj = crs(st))
  # 
  # # extract ndvi - replace with Heather's function?
  # # calculate date of peak IRG for each used cell?
  # xx<-extract.ndvi(folder="C:/Users/tm73z/OneDrive - University of Glasgow/Serengeti Wildebeest Project/SerengetiGIS/NDVI_mean",      # folder in which NDVI_mean.tif is saved
  #                  stack=st,       # Raster stack output from above build.ndvi.stack
  #                  stack.desc=sd,  # Stack description from above build.ndvi.description
  #                  loc=loc,         # All x-y locations projected in latlong
  #                  date=date,        # as.Date from loc, e.g. date column from your GPS points
  #                  cpus=3,
  #                  extract.andvi=T)
  # 
  # # combine with data
  # dat2 <- cbind(dat2,xx)
  # rm(loc,date,xx)
  
  # run days from peak IRG code
  # if(spp %in% c('BS','EK','MD','MS','PH')){
  # tndvi <- temphetero(XYdata=xy,
  #                     NDVIfolder="C:/Users/tm73z/OneDrive - University of Glasgow/Documents/Side projects/Fidelity/StatewideNDVI/ndvi_bisc_calcs",
  #                     dates=ymd_hms(xy$date),
  #                     AID=unique(xy$AID),
  #                     metrics=c("IRG","DFP"),
  #                     wantSDtoo=TRUE)
  # print('Done: IRG calcs')
  # }else{
  #   tndvi <- cbind(as.data.frame(xy),
  #                 data.frame(
  #                      IRG=NA, DFP=NA))
  # }
  # tndvi <- as.data.frame(tndvi)
  # tndvi <- tndvi[,names(tndvi) %in% c('IRG','DFP')]
  # dat2 <- cbind(dat2,tndvi)
                
  # # idenitfy window sizes for each species using NSD
  # wind.size.NSD(dat2)
  
  # predictability - Ctime metrics
  # pred.val <- extract(pred,xy)
  # dat2$ctime_entropy <- pred.val[,1]
  # dat2$ctime_1sd <- pred.val[,2]
  # dat2$ctime_periodicity <- pred.val[,3]
  # rm(pred,pred.val)
  
  # calculate fidelity with various window sizes
  # print('start iyd')
  # vv <- iyd(dat2,wind,uproj)

  # combine with data
  # dat2 <- cbind(dat2,vv)
  # rm(date,vv)
  # print('Done: iyd')
  
  # predictability - Cspace metrics
  # Cspace <- extract_Cspace(xy,dat2,modfold)
  # dat2 <- cbind(dat2,Cspace)
  # print('Done: Predictability metrics')
  
  # identify movement state using HMM
  # dat2 <- as.data.frame(dat2)
  # dat2$states<-hmm2state(dat2)

  # add other data
  dat3 <- data.frame(dat2,
                     dtrut[which(dtrut$spp==spp),]
  )

  # calculate homerange sizes
  xy <- spTransform(xy,CRSobj = uproj)
  if(q==1) {
    hr <- fallsnow(xy)
    # hr.x <- as.data.frame(mcp.by.id(xy))
    # hr <- data.frame(AID=names(hr.x),hr=unlist(hr.x))
    
  }else{
    zed <- fallsnow(xy)
    if(all(names(zed) == names(hr) & !is.null(zed))){
      hr <- rbind(hr,zed)
    }
    
    # hr.x <- as.data.frame(mcp.by.id(xy))
    # hr <- rbind(hr,data.frame(AID=names(hr.x),hr=unlist(hr.x)))
    # 
  }
  print('Done: homeranges')
  
  # make sure everything is ordered correctly
  dat3 <- dat3[,!names(dat3) %in% c("hour","year" ,"jul","hour_since_smphour")]
 
  # add together
  if(q==1) {
    alldat <- dat3
  }else{
    alldat <- rbind(alldat,dat3)
  }
  print(paste('done',title,q, Sys.time()))
  
  # write.csv(alldat,paste0('./output/',fold,'/','fidelity_all.csv'),row.names = FALSE)
  # write.csv(hr,paste0('./output/',fold,'/','MCPsqkm_all_.csv'),row.names = FALSE)
}
Sys.time()-runtime


write.csv(alldat,paste0('./output/output_2020-07-21/','fidelity_DFP.csv'),row.names = FALSE)

x <- read.csv("./NSD/nsdsumJune2020.csv")

dat2 %>%
  group_by(AID) %>%
  summarise(iydREAL)

### examine the effect of window size
# meanfid <- colMeans(log(alldat[,18:ncol(alldat)]),na.rm=T)
# u95 <- meanfid + (1.96*apply(alldat[,18:ncol(alldat)],2,function(x) sd(log(x),na.rm=T)))
# l95 <- meanfid - (1.96*apply(alldat[,18:ncol(alldat)],2,function(x) sd(log(x),na.rm=T)))
# plot(wind,meanfid,type='l',ylim=c(min(l95,na.rm=T),max(u95,na.rm=T)))
# lines(wind,u95,type='l',lty=2)
# lines(wind,l95,type='l',lty=2)

# predictability
predak <- stack('./predictability/Ctime_Alaska.tif')
predwy <- stack('./predictability/Ctime_Wyoming.tif')
predse <- stack('./predictability/Ctime_WB_SE.tif' )

par(mfcol=c(3,3),mai=c(0.2,0.4,0.6,0.2))
#ak
plot(predak[[1]],main='entropy')
GISTools::map.scale(x=-140, y=65, len = 9.6, ndivs = 2, subdiv = 400, units = "km")
aksd <- predak[[2]]
aksd[aksd>10] <- NA
plot(aksd,main='1/sd')
mtext('Alaska north slope',side = 3)
plot(predak[[3]],main='periodicity')


#wy
plot(predwy[[1]],main='entropy')
GISTools::map.scale(x=-106.5, y=40.9, len = 4.8, ndivs = 2, subdiv = 200, units = "km")
wysd <- predwy[[2]]
wysd[wysd>10] <- NA
plot(wysd,main='1/sd')
mtext('Wyoming',side = 3)
plot(predwy[[3]],main='periodicity')
p <- predwy[[3]]
summary(getValues(p))
#serengeti
plot(predse[[1]],main='entropy')
GISTools::map.scale(x=35, y=-3.1, len = 1.2, ndivs = 2, subdiv = 50, units = "km")
plot(predse[[2]],main='1/sd')
mtext('Serengeti',side = 3)
plot(predse[[3]],main='periodicity')
