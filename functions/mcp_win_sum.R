
# setwd(paste('C:/Users/tmorri19/Documents/GitHub/SiteFidelityProj/data/',title,sep=""))
# dat3<-read.csv("PlatteValleyMuleDeer_smp.csv")


####### CREATE SPATIAL DATA WITH WY PROJECTION

mcps<-function(vv,title) {
  require('SP')
  require("rgdal")
  require("adehabitatHR")
  
  setwd('C:/Users/tmorri19/Documents/GitHub/SiteFidelityProj/')
  pro<-read.table("WY_projection.txt",header=FALSE)
  pro<-as.character(pro[[1]])
  coordinates(vv)<-~x+y
  proj4string(vv)<-CRS(pro)
  vv@data$AID<-as.factor(vv@data$AID)
  uaid<-unique(vv$AID)
  #### generate summer and winter range MCP for each year and each animal
  setwd(paste('C:/Users/tmorri19/Documents/GitHub/SiteFidelityProj/data/',title,sep=""))
  uy<-unique(format(as.POSIXct(vv$date,tz="MST"),"%Y"))
  migd<-yrs<-NULL
  for(i in uy) {
    tmp<-vv[format(as.POSIXct(vv$date,tz="MST"),"%Y")==i,]
    
    win<-tmp[months(as.POSIXct(tmp$date,tz="MST"))=="January"|months(as.POSIXct(tmp$date,tz="MST"))=="February",]
    if(nrow(win)>5){
      mcpwin<-mcp(win[,2],percent=95)
      writeOGR(mcpwin,".",paste(title,"_",as.numeric(i),"_mcpwinter",sep=""),driver="ESRI Shapefile")
    }  
    sum<-tmp[months(as.POSIXct(tmp$date,tz="MST"))=="July"|months(as.POSIXct(tmp$date,tz="MST"))=="August",]
    if(nrow(sum)>5){
      mcpsum<-mcp(sum[,2],percent=95)
      writeOGR(mcpsum,".",paste(title,"_",as.numeric(i),"_mcpsummer",sep=""),driver="ESRI Shapefile")
    }
    
    for (j in 1:length(uaid)){
      tmpW<-win[win$AID==uaid[j],]
      tmpS<-sum[sum$AID==uaid[j],]
      
      mwin<-c(mean(tmpW$x),mean(tmpW$y))  
      msum<-c(mean(tmpS$x),mean(tmpS$y))
      migd<-c(migd,dist(rbind(mwin,msum)))
    }
    yrs<-c(yrs,rep(as.numeric(i),length(uaid)))
  }
  migd<-cbind(migd,rep(as.numeric(uaid),length(uy)),yrs)
}  