## calculates fidelity (inter-annual distances) and relates it to behavioural states
## updated Aug 2 2019

rm(list=ls(all=TRUE))
require(lme4)
require(lubridate)
require(sp)
require(rgdal)
require(rgeos)

# packages
libs <- c("moveHMM","foreign","lme4","maptools","graphics","stats","cluster","raster","Hmisc","MASS","igraph","CircStats",
          "sp","rgdal","adehabitatHR","adehabitatLT","geoR","stringr","rgeos","snowfall","lubridate","spatstat")

lapply(libs, require, character.only = TRUE)

# load functions
swd<-"C:/Users/tmorri19/OneDrive - University of Glasgow/Documents/Side projects/Fidelity/"

setwd(paste(swd,"functions",sep=""))
funks<-c("iyd.r",'edgelist.r','smp_edgelist.r','mcp_all.r','ovlap.r',"make_clust.R",
         'temporal_hetero_extraction.r','temporal_NDVIintegrated.r','temporal_NDVIintegrated.r',
         "season.R","IRGfunction.R","rm_crit.R","torus_shift.R","find.encountr.R")
for(i in 1:length(funks)){
  source(funks[i]) 
}

setwd("C:/users/tmorri19/OneDrive - University of Glasgow/Documents/Side projects/Fidelity/functions/")
source('downsample_GPS_pts.R')

# projections
wyproj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
ca_alproj <- "+proj=utm +north +zone=6 +datum=WGS84"
wb_amproj <- "+proj=utm +south +zone=37 +datum=WGS84" #amboseli
wb_seproj <- "+proj=utm +south +zone=36 +init=epsg:21036" #serengeti
latlongcrs<-("+proj=longlat +datum=WGS84")

# compile datasets and calculate inter-annual distances 
spp<-c("BS","CA","EK","MD","MS","PH","WB","ZB")
cps<-5 # min number of individuals required to have a single population
elbuff<-2500 #min distance (in meters) in which two individuals must occur on same day to be included in edgelist function
fold<-"hmm_analysis"
wind<-84                    ## no. days within each window for comparison !should be even number!
crit<-365+round(wind/2)     ## no. days required for each individual to be included in analysis
dtrut<-data.frame(spp=spp,
                  rut=c(318,296,273,333,272,277,152,15),
                  parturition=c(144,151,151,169,142,158,46,15),
                  sum1=c(rep(182,length(spp)-2),106,106),     #JUL 1 - TEMPERATE; APR 1 - TROPICS
                  sum2=c(rep(246,length(spp)-2),166,166),   #AUG 31 - TEMPERATE; MAY 31 - TROPICS
                  win1=c(rep(1,length(spp)-2),213,213),     #JAN 1 - TEMPERATE; AUG 1 - TROPICS
                  win2=c(rep(60,length(spp)-2),304,304))    #FEB 28 - TEMPERATE; OCT 31 - TROPICS


setwd(paste(swd,"data",sep=""))
datasets <- list.files()

dat5 <- NULL
for(q in 28:length(datasets)){ 
  
  print(paste("yo, this be dataset:", datasets[q],q))
  title<-datasets[q]
  setwd(paste0(swd,'data/',title))
  
  dat<-read.csv(paste0(title,'.csv'),header=TRUE)
  
  if( !fold %in% list.files(swd) ){ dir.create(paste0(swd,fold))}
  if( !title %in% list.files(paste0(swd,fold)) ){ dir.create(paste0(swd,fold,"/",title))}
  setwd(paste0(swd,fold,"/",title))
  
  ####### FORMAT DATA
  if(length(names(dat)[names(dat)=="?..AID"])>0) names(dat)[names(dat)=="?..AID"] <-"AID"
  if(length(names(dat)[names(dat)=="ï..AID"])>0) names(dat)[names(dat)=="ï..AID"] <-"AID"
  if(length(names(dat)[names(dat)=="TelemDate"])>0) names(dat)[names(dat)=="TelemDate"] <- "TelemDt"
  
  # assign tzone and CRS to objects     
  if(!title %in% c("WB_SE","ZB_SE","WB_AM","CA_AL") & nrow(dat3)>0){
    proj <- wyproj
    tzone <- "MST"}
  if(title %in% c("WB_SE","ZB_SE") & nrow(dat3)>0){
    proj <- wb_seproj
    tzone <- "Africa/Nairobi"}
  if(title %in% c("WB_AM") & nrow(dat3)>0){
    proj <- wb_amproj
    tzone <- "Africa/Nairobi"}
  if(title %in% c("CA_AL") & nrow(dat3)>0){
    proj <- ca_alproj
    tzone <- "America/Juneau"}
  
  dat$date <- as.POSIXct(strptime(as.character(dat$TelemDt), format = "%Y-%m-%d %H:%M:%S", tz = "GMT"))
  attributes(dat$date)$tzone <- tzone
  
  dat<-dat[,names(dat) %in% c("SEX","SQL_Shp","Elvtn_m","Date","Time","TmLg_Al","TelemDt","Status")==FALSE] #c(1,5:9,12,13)]
  dat<-dat[is.na(strptime(dat$date,"%Y-%m-%d %H:%M:%S"))==FALSE,]
  if(length(names(dat)[names(dat)=="AEA_X"])>0) names(dat)[names(dat)=="AEA_X"] <- "x"
  if(length(names(dat)[names(dat)=="AEA_Y"])>0) names(dat)[names(dat)=="AEA_Y"] <- "y"
  if(length(names(dat)[names(dat)=="x"])==0|length(names(dat)[names(dat)=="y"])==0 
     |length(names(dat)[names(dat)=="AID"])==0|length(names(dat)[names(dat)=="date"])==0) 
    stop("rename columns in input data")  
  if(title!="WB_SE" & title!="ZB_SE" & title!="WB_AM" & title!="CA_AL") dat<-dat[dat$x>470000,]
  
  # remove individuals without crit amount of data
  dat3<-rm_crit(dat,crit)
  
  # Sample for 1 relocation per day
  dat3<-smp(dat3)
  dat3$AID<-as.factor(as.character(dat3$AID))
  dat3<-dat3[order(dat3$AID,dat3$date),]
  
  uAID<-unique(dat3$AID)
  
  dat3$iyd<-NA
  dat4 <- NULL
  
  # loop across all unique individuals
  for(i in 1:length(uAID)){ 
    tmp<-dat3[dat3$AID==uAID[i],]
    sd<-min(as.Date(tmp$date),na.rm=TRUE)+365+(wind/2)
    
    # loop across all unique dates for this individuals after date=sd
    for(f in nrow(tmp[as.Date(tmp$date)<sd,]):nrow(tmp)) { #(nrow(tmp[as.Date(tmp$date)<sd,])+10)){  #
      
      #define focal date
      cd<-as.Date(tmp$date[f],tz="MST") #focal date yr2
      
      # define window
      dt1<-cd-365-(wind/2)     #define wind yr1 start date
      dt2<-cd-365+(wind/2)     #define wind yr1 end date 
      
      # first check that the window size is at least 42 days
      if(length(tmp$x[as.Date(tmp$date)>=dt1 & as.Date(tmp$date)<=dt2]) > wind/2){
        
        # calculate the log-min distance btw Year 2 and Year 1    
        tmp$iyd[f] <- log(min(dist(rbind(
          cbind(tmp$x[f],tmp$y[f]), # x-y from year 2
          (cbind(tmp$x[as.Date(tmp$date)>=dt1 & as.Date(tmp$date)<=dt2],     # x from year 1
                 tmp$y[as.Date(tmp$date)>=dt1 & as.Date(tmp$date)<=dt2]))))  # y from year 1
          [1:length(tmp$x[as.Date(tmp$date)>=dt1 & as.Date(tmp$date)<=dt2])] # select elements for yr 2 vs yr 1
          ,na.rm=T))
      }# end if statement
      
      #compile

    }# end loop across dates
    
    # calculate states for each individual using a hmm model
    tmp$states<-hmm2state(tmp,title)
    dat4<-rbind(dat4,tmp)

  }# end loop across individuals
  dat4<-dat4[,names(dat4) %in% c("AID","x","y","date","hour", "year", "jul","iyd", "states")] 
  dat4$pop <- title
  dat5 <- rbind(dat5,dat4)
  
}# end loop across datasets

x<-dat5
dat5 <- dat5[!duplicated(paste(dat5$AID,dat5$date)),]
setwd(swd)
setwd("C:/Users/tmorri19/OneDrive - University of Glasgow/Documents/Side projects/Fidelity/hmm_analysis/")
dat5$species <- substr(dat5$pop,1,2)
write.csv(dat5,"hmm_iyd_alldata.csv",row.names = F)

#########################################################
setwd("C:/Users/tmorri19/OneDrive - University of Glasgow/Documents/Side projects/Fidelity/hmm_analysis/")
dat<-read.csv("hmm_iyd_alldata.csv")
dat<-dat[!is.na(dat$iyd),]
dat<-dat[dat$iyd != -Inf,]


mod1 <- lmer(iyd ~ states+species + (1|pop/AID),data=dat)
summary(mod1)

mod2 <- lmer(iyd ~ states*species + (1|pop/AID),data=dat)
summary(mod2)

AIC(mod1,mod2)

## bootstrap CI on the prediction
pred <- predict(mod2,re.form=NA)
predFun<-function(x,newdat){
  predict(x,newdat)
}
boot<-bootMer(mod2,FUN=predFun,nsim=100,re.form=NA)
CI.lo<-apply(boot$t, 2, function(x) as.numeric(quantile(x, probs=0.025, na.rm=TRUE)))
CI.hi<-apply(boot$t, 2, function(x) as.numeric(quantile(x, probs=0.975, na.rm=TRUE)))

# std.err<-apply(boot$t,2,sd)
# CI.lo <- pred - std.err*1.96
# CI.hi <- pred + std.err*1.96

ranef(mod2)[[2]]

tmp<-paste0(dat$AID,dat$pop)
tmp<-unique(tmp)

library(merTools)
newDat <- expand.grid(states=c(1,2),
                      species=unique(dat$species),
                      AID=unique(dat$AID),
                      pop=unique(dat$pop))
newDat <- newDat[paste0(newDat$AID,newDat$pop) %in% tmp,]
newDat <- newDat[substr(newDat$pop,1,2) == newDat$species,]
# newDat <- newDat[!duplicated((paste0(newDat$states,newDat$species,newDat$pop))),]

preds <- predictInterval(mod2, newdata = newDat, n.sims = 999,level = 0.5)
newDat <- cbind(newDat,preds)

write.csv(newDat,file = "predictionIntervals.csv",row.names = F)

z<-aggregate(newDat$fit, by=list(newDat$states,newDat$species),mean)
zsd<-aggregate(newDat$fit, by=list(newDat$states,newDat$species),sd)
zlen<-aggregate(newDat$fit, by=list(newDat$states,newDat$species),length)

zsd[,3]<-zsd[,3]/sqrt(zlen[,3])

z<-z[c(7,8,9,10,1,2,5,6,11,12,13,14,15,16,3,4),]
zsd <- zsd[c(7,8,9,10,1,2,5,6,11,12,13,14,15,16,3,4),]

colr<-c("blue","blue","purple","purple",
  "red","red","yellow","yellow","orange","orange",
  "green","green","black","black","lightblue","lightblue")
### PLOTTING
plot(seq(1:nrow(z)),z[,3],xaxt='n',#xlab="Species",
     ylab="Inter-year distance (log-km)",pch=15,cex=1.1,ylim=c(4,12),
     col=colr,
     xlab=paste("Mule deer","Moose","Bighorn sheep","Elk","Pronghorn","Wildebeest","Zebra","Caribou",
     sep="   "),cex.lab=0.8)
arrows(x0 = seq(1:nrow(z)),y0 = z[,3],
       x1 = seq(1:nrow(z)),y1 = z[,3]+zsd[,3]*1.96,angle = 90,length=0.07)
arrows(x0 = seq(1:nrow(z)),y0 = z[,3],
       x1 = seq(1:nrow(z)),y1 = z[,3]-zsd[,3]*1.96,angle = 90,length=0.07)
points(seq(1:nrow(z)),z[,3],pch=15,cex=1.1,col=colr)
points(seq(1:nrow(z)),z[,3],pch=0,cex=1.1)

### MCP analysis
require(adehabitatHR)
bhs<-dat[dat$species=="BS",]
upop<-unique(bhs$pop)

setwd("C:/Users/tmorri19/OneDrive - University of Glasgow/Documents/Side projects/Fidelity/")
# pro<-read.table("WY_projection_AEA.txt",header=FALSE)
pro<-read.table("WY_projection_UTM12.txt",header=FALSE)
pro<-as.character(pro[[1]])

coordinates(bhs)<-~x+y
proj4string(bhs)<-CRS(pro)
bhs$pop<-factor(bhs$pop)
upop<-unique(bhs$pop)

x<-bhs[bhs$pop=="BS_DC",]
mcp(x[names(x)=="pop"],percent=99)

mcp99<-mcp(bhs[names(bhs)=="pop"],percent=99)
mcp95<-mcp(bhs[names(bhs)=="pop"],percent=95)


library(mapview)
library(raster)
library(rgdal)
library(rgeos)

sqrt(gArea(mcp99))/1000

mapView(bhs[bhs$pop=="BS_DC",]) + mapView(mcp99)
