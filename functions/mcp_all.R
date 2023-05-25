# 
# setwd(paste('C:/Users/tmorri19/Documents/GitHub/SiteFidelityProj/data/',title,sep=""))
# dat3<-read.csv("PlatteValleyMuleDeer_smp.csv")
# 
# migdist<-mcp_all(vv=dat2[dat2$pop==u[k],],title=title,pop=u[k]) #Generates MCPS and returns migration distances btw mean summer and winter ranges

####### CREATE SPATIAL DATA WITH WY PROJECTION AND RETURN MIGRATION DISTANCE

mcp_all<-function(vv,title,pop,swd) {
  xdir<-getwd()
  setwd(swd)
  pro<-read.table("WY_projection.txt",header=FALSE)
  if(title=="WB_SE" | title=="ZB_SE") pro<-read.table("SE_projection.txt",header=FALSE)
  if(title=="WB_AM") pro<-read.table("AM_projection.txt",header=FALSE)
  if(title=="CA_AL") pro<-read.table("AL_projection.txt",header=FALSE)
  
  setwd(xdir)
  nset<-vv[1,]
  rownames(nset)<-"0"
  pro<-as.character(pro[[1]])
  coordinates(vv)<-~x+y; coordinates(nset)<-~x+y
  proj4string(vv)<-CRS(pro); proj4string(nset)<-CRS(pro)
  vv@data$AID<-as.factor(vv@data$AID); nset@data$AID<-as.factor(nset@data$AID)
  uaid<-unique(vv$AID)

  #### make MCP for entire pop
  vv$AID<-factor(vv$AID)
  mcp99<-mcp(vv[names(vv)=="AID"],percent=99)
  writeOGR(mcp99,".",paste(title,"_",pop,"_ALL_MCP_99",sep=""),driver="ESRI Shapefile",overwrite_layer=TRUE)
  mcp95<-mcp(vv[names(vv)=="AID"],percent=95)
  writeOGR(mcp95,".",paste(title,"_",pop,"_ALL_MCP_95",sep=""),driver="ESRI Shapefile",overwrite_layer=TRUE)
  writeOGR(vv,".",paste(title,"_",pop,"_ALL_POINTS",sep=""),driver="ESRI Shapefile",overwrite_layer=TRUE)
  
  setwd(paste(swd,"/data/",title,sep=""))
  writeOGR(mcp95,".",paste(title,"_",pop,"_ALL_MCP_95",sep=""),driver="ESRI Shapefile",overwrite_layer=TRUE)
  setwd(xdir)
  
  #### generate summer and winter range MCP for each year and each animal
  uy<-unique(format(as.POSIXct(vv$date,tz="MST"),"%Y"))
  migd<-yrs<-NULL
  
  migdist <- do.call(rbind, lapply(1:length(uaid), function(j){ #
    tmp<-vv[as.character(vv$AID)==uaid[j],]
    sumyr<-max(year(tmp$date[months(tmp$date)=="August"]))
    winyr<-max(year(tmp$date[months(tmp$date)=="February"]))
    tmp<-tmp[paste(year(tmp$date),months(tmp$date),sep="")==paste(sumyr,"August",sep="") |
              paste(year(tmp$date),months(tmp$date),sep="")==paste(sumyr,"July",sep="") |
              paste(year(tmp$date),months(tmp$date),sep="")==paste(winyr,"January",sep="") |
              paste(year(tmp$date),months(tmp$date),sep="")==paste(winyr,"February",sep="") ,]

    ####### Calc migration distance
    win<-tmp[months(as.POSIXct(tmp$date,tz="MST"))=="January"|months(as.POSIXct(tmp$date,tz="MST"))=="February",]
    sum<-tmp[months(as.POSIXct(tmp$date,tz="MST"))=="July"|months(as.POSIXct(tmp$date,tz="MST"))=="August",]
  
    if(title=="WB_SE" | title=="ZB_SE" | title=="WB_AM"){
      tmp<-vv[as.character(vv$AID)==uaid[j],]
      sumyr<-max(year(tmp$date[months(tmp$date)=="May"]))
      winyr<-max(year(tmp$date[months(tmp$date)=="October"]))
      tmp<-tmp[paste(year(tmp$date),months(tmp$date),sep="")==paste(sumyr,"March",sep="") |
                 paste(year(tmp$date),months(tmp$date),sep="")==paste(sumyr,"April",sep="") |
                 paste(year(tmp$date),months(tmp$date),sep="")==paste(sumyr,"May",sep="") |
                 paste(year(tmp$date),months(tmp$date),sep="")==paste(winyr,"August",sep="") |
                 paste(year(tmp$date),months(tmp$date),sep="")==paste(winyr,"September",sep="") |
                 paste(year(tmp$date),months(tmp$date),sep="")==paste(winyr,"October",sep="") 
               ,]
      
      ####### Calc migration distance
      win<-tmp[months(as.POSIXct(tmp$date,tz="MST"))=="August"|months(as.POSIXct(tmp$date,tz="MST"))=="September"|months(as.POSIXct(tmp$date,tz="MST"))=="October",]
      sum<-tmp[months(as.POSIXct(tmp$date,tz="MST"))=="March"|months(as.POSIXct(tmp$date,tz="MST"))=="April"|months(as.POSIXct(tmp$date,tz="MST"))=="May",]
    }
    
    mwin<-c(mean(win$x,na.rm=TRUE),mean(win$y,na.rm=TRUE))  #Identify centroid of winter
    msum<-c(mean(sum$x,na.rm=TRUE),mean(sum$y,na.rm=TRUE))  #Identify centroid of summer
    migd<-c(migd,dist(rbind(mwin,msum))) ## Calc distance between win and summ ranges
    return(migd)
  }))

  migdist<-matrix(migdist,ncol=1)
  
  return(migdist)
}