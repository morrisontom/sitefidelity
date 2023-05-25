library(sp)
library(rgdal)

#### PRONGHORN
setwd("E:/morrison/SiteFidelityProj/other data/")
dat<-read.csv("PH1.txt",header=FALSE)
dim(dat)
head(dat)
colnames(dat)<-c("locationID","AID","del","date","year","month","day","hour","lat","long","x","y","alt","horiz")
dat$date<-strptime(as.character(dat$date), format = "%m/%d/%Y", tz ="MST")+(dat$hour*60*60)
coordinates(dat) <- c("x","y")
crs<-CRS("+proj=utm +north +zone=12 +datum=NAD83")
proj4string(dat)<-crs
writeOGR(dat,"E:/morrison/SiteFidelityProj/maps","PH_RD3",driver="ESRI Shapefile")
names(dat$TelemDt)<-names(dat$date)
dat<-as.data.frame(dat)
dat$TelemDt<-dat$date
write.table(dat,file="MD_RD.txt")

#### ROCK SPRINGS PRONGHORN
setwd("E:/morrison/SiteFidelityProj/data/")
dat<-read.csv("DEER.txt",header=FALSE)
head(dat)
colnames(dat)<-c("locationID","AID","del","date","year","month","day","hour","lat","long","x","y","alt","horiz")
dat$date<-strptime(as.character(dat$date), format = "%m/%d/%Y", tz ="MST")+(dat$hour*60*60)
coordinates(dat) <- c("x","y")
crs<-CRS("+proj=longlat +datum=NAD83")
crs<-CRS("+proj=utm +north +zone=12 +datum=NAD83")
proj4string(dat)<-crs
writeOGR(dat,"E:/morrison/SiteFidelityProj/maps","MD_RD3",driver="ESRI Shapefile")
names(dat$TelemDt)<-names(dat$date)
dat<-as.data.frame(dat)
dat$TelemDt<-dat$date
write.table(dat,file="MD_RD.txt")



u<-as.character(unique(dat$captureLoc))
for (i in 1:length(unique(dat$captureLoc))){
  dat2<-dat[as.character(dat$captureLoc)==u[i],]
  write.table(dat2,paste("EK_",unique(dat$captureLoc)[i],".txt",sep=""),col.names=TRUE)
}


### CODY ELK

setwd("E:/morrison/SiteFidelityProj/data/EK_CY")
dat<-read.csv("EK_CY.csv")

dim(dat)
dat<-dat[!is.na(dat$Latitude),]
head(dat)
dat<-dat[order(dat$DATE2),]
check<-aggregate(dat$DATE2,by=list(dat$AID),function(x){(unique(x)[1])})

dat<-dat[dat$DATE2>(41708.04+7),] ## Remove 1-2 weeks of data

coordinates(dat)<-c("Longitude","Latitude")
proj4string(dat)<-CRS("+proj=longlat +datum=NAD83") 
dat2<-spTransform(dat, CRS("+proj=utm +zone=12 ellps=NAD83")) 
dat3<-as.data.frame(dat2)
colnames(dat3)[6]<-"x"
colnames(dat3)[7]<-"y"
head(dat3)
write.table(dat3,"EK_CY.txt",col.names=TRUE)

### AB elk
setwd("E:/morrison/SiteFidelityProj/data/EK_AB")
dat<-read.delim("AB_EK.txt")
head(dat)
write.table(dat3,"EK_AB.txt",col.names=TRUE)

### Bighorn
setwd("E:/morrison/SiteFidelityProj/")
pro<-read.table("WY_projection.txt",header=FALSE)
pro<-as.character(pro[[1]])

setwd("E:/morrison/SiteFidelityProj/data/")
vv<-c("BS_GT")#"BS_JK","BS_LP","BS_SE","BS_WB") #"BS_DC","BS_EM",
for (j in 1:length(vv)){
  setwd(paste("E:/morrison/SiteFidelityProj/data/",vv[j],sep=""))
  file<-list.files()
  dat<-read.delim(file,header=TRUE)
  if(length(names(dat)[names(dat)=="AEA_X"])>0) names(dat)[names(dat)=="AEA_X"] <- "x"
  if(length(names(dat)[names(dat)=="AEA_Y"])>0) names(dat)[names(dat)=="AEA_Y"] <- "y"
  if(length(names(dat)[names(dat)=="Easting"])>0) names(dat)[names(dat)=="Easting"] <- "x"
  if(length(names(dat)[names(dat)=="Northing"])>0) names(dat)[names(dat)=="Northing"] <- "y"
  if(length(names(dat)[names(dat)=="TelemDate"])>0) names(dat)[names(dat)=="TelemDate"] <- "TelemDt"
  # head(dat)
  coordinates(dat)<-c("x","y")
  proj4string(dat)<-CRS(pro) 
  dat2<-spTransform(dat, CRS("+proj=utm +zone=12 ellps=NAD83")) 
  dat3<-as.data.frame(dat2)
  write.table(dat3,paste(vv[j],".txt",sep=""),col.names=TRUE)
}

### Mule deer
setwd("E:/morrison/SiteFidelityProj/data/MD_PD")
dat<-readOGR('.',"pd_3d_fix_wyckoff_point_data")
dat2<-as.data.frame(dat)

if(length(names(dat2)[names(dat2)=="AEA_X"])>0) names(dat2)[names(dat2)=="AEA_X"] <- "x"
if(length(names(dat2)[names(dat2)=="AEA_Y"])>0) names(dat2)[names(dat2)=="AEA_Y"] <- "y"
if(length(names(dat2)[names(dat2)=="Easting"])>0) names(dat2)[names(dat2)=="Easting"] <- "x"
if(length(names(dat2)[names(dat2)=="Northing"])>0) names(dat2)[names(dat2)=="Northing"] <- "y"
if(length(names(dat2)[names(dat2)=="Telemdate"])>0) names(dat2)[names(dat2)=="TelemDate"] <- "TelemDt"
coordinates(dat2)<-c("x","y")
proj4string(dat2)<-CRS(pro) 
dat2<-spTransform(dat2, CRS("+proj=utm +zone=12 ellps=NAD83")) 
dat3<-as.data.frame(dat2)
write.table(dat3,"MD_PD.txt",col.names=TRUE)

### PRONGY's
setwd("E:/morrison/SiteFidelityProj/data/PH_KE")
dat<-read.delim("KE_PH.txt")
head(dat)
dat2<-dat
if(length(names(dat2)[names(dat2)=="AEA_X"])>0) names(dat2)[names(dat2)=="AEA_X"] <- "x"
if(length(names(dat2)[names(dat2)=="AEA_Y"])>0) names(dat2)[names(dat2)=="AEA_Y"] <- "y"
if(length(names(dat2)[names(dat2)=="Easting"])>0) names(dat2)[names(dat2)=="Easting"] <- "x"
if(length(names(dat2)[names(dat2)=="Northing"])>0) names(dat2)[names(dat2)=="Northing"] <- "y"
if(length(names(dat2)[names(dat2)=="Telemdate"])>0) names(dat2)[names(dat2)=="TelemDate"] <- "TelemDt"
coordinates(dat2)<-c("x","y")
proj4string(dat2)<-CRS(pro) 
dat2<-spTransform(dat2, CRS("+proj=utm +zone=12 ellps=NAD83")) 
dat3<-as.data.frame(dat2)
write.table(dat3,"PH_KE.txt",col.names=TRUE)

