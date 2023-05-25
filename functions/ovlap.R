
####### CREATE SPATIAL DATA WITH WY PROJECTION AND RETURN MIGRATION DISTANCE

ovlap<-function(vv,title,pop,swd,wind) {
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
  
  overlap.index<-list()
  mean.overlap<-NULL
  for(g in 1:length(uaid)){
    tmp<-vv[vv$AID==uaid[g],]
    stdt<-min(tmp$date)+60*60*24*(365+wind/2)
    endt<-max(tmp$date)
    em<-tmp$date[tmp$date>=stdt & tmp$date<=endt]
    pval<-NULL
    
    for(h in 1:length(em)){
      emly<-c(em[h]-60*60*24*(365+wind/2),em[h]-60*60*24*(365-wind/2))
      tmp2<-tmp[tmp$date>=emly[1] & tmp$date<=emly[2],]    
      if(length(tmp2)>=(wind/2)) {
        kde<-kernelUD(tmp2)
        ras<-raster(as(kde,"SpatialPixelsDataFrame"))
        # ras[ras==Inf]<-NA
        # ras[ras==-Inf]<-NA
        # ras[ras<0.000000000001]<-NA
        ras<-ras*(1/sum(getValues(ras),na.rm=T))
        pval<-c(pval,extract(ras,tmp[tmp$date==em[h],]))
      }else{
        pval<-c(pval,NA)
      }
      # plot(ras)
      # plot(tmp[tmp$date==em[h],],add=T,col="blue")
    }
    overlap.index[[g]]<-data.frame(AID=uaid[g],pval=pval,dates=em)
    mean.overlap<-rbind(mean.overlap,data.frame(AID=uaid[g],pval=mean(pval,na.rm=T),length=length(em)))
  }
  return(mean.overlap)
}