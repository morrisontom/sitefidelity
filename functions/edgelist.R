##############################################################################################
# Generates an edgelist with the #days in which each pair of individual was within some distance
# T. Morrison
# May 22, 2015

net<-function(x,buff,title,swd){  ## x is object with $AID $date $x $y from single population 
  x<-x[is.na(x$AID)==FALSE,]
  x$date<-as.POSIXct(x$date)
  u<-unique(as.Date(x$date,tz='MST'))
  aid<-unique(x$AID)
  aid<-aid[order(aid)]
  zz<-expand.grid(aid,aid)
  zz[,3]<-0
  
  ## Remove all duplicated pairs 
  srn<-sqrt(nrow(zz))
  xx<-NULL
  for(z in 2:length(aid)){
    xx<-c(xx,seq((srn*(z-1))+1,(srn*(z-1))+(z-1),1))
  }  
  zz<-zz[xx,]
  zz<-zz[is.na(zz[,1])==FALSE,]
  if(nrow(zz)!=(length(aid)*(length(aid)-1)/2)) STOP("Something is screwy")
  zz<-zz[order(zz[,1]),]
  
  ## NEW
  setwd(swd)
  #   pro<-read.table("WY_projection.txt",header=FALSE)
  #   if(title=="WB_SE" | title=="ZB_SE") pro<-read.table("SE_projection.txt",header=FALSE)
  #   if(title=="WB_AM") pro<-read.table("AM_projection.txt",header=FALSE)
  #   if(title=="CA_AL") pro<-read.table("AL_projection.txt",header=FALSE)
  #   pro<-as.character(pro[[1]])
  #   coordinates(x)<-~x+y
  #   proj4string(x)<-CRS(pro)
  #  2842	30086
  
  ## NEW: 
  
  for(i in 1:nrow(zz)){
    add<-0
    
    ## Does #1 overlap with #2
    tmp1<-x[x$AID==zz[i,1],]
    tmp2<-x[x$AID==zz[i,2],]
    
    if(min(tmp1$date)>max(tmp2$date) | min(tmp2$date)>max(tmp1$date)){
      zz[i,3]<-0
    }else{
      for(j in 1:length(tmp1$date)){
        xx<-as.matrix(dist(rbind(c(tmp1$x[j],tmp1$y[j]),cbind(tmp2$x[abs(tmp2$date-tmp1$date[j])<=24],tmp2$y[abs(tmp2$date-tmp1$date[j])<=24]) )))
        if(dim(xx)[1]>1) {xx<-min(xx[2:nrow(xx),1],na.rm=TRUE)
        }else{xx<-buff+1}
        if(xx<=buff) add<-add+1
      }
      
      ## Does #2 overlap with #1?
      tmp1<-x[x$AID==zz[i,2],]
      tmp2<-x[x$AID==zz[i,1],]
      for(j in 1:length(tmp1$date)){
        xx<-as.matrix(dist(rbind(c(tmp1$x[j],tmp1$y[j]),cbind(tmp2$x[abs(tmp2$date-tmp1$date[j])<=24],tmp2$y[abs(tmp2$date-tmp1$date[j])<=24]))))
        if(dim(xx)[1]>1) {xx<-min(xx[2:nrow(xx),1],na.rm=TRUE)}
        else{xx<-buff+1}
        if(xx<=buff) add<-add+1
      }
      zz[i,3]<-add
      # print(paste("row",i,add))
    }
  }  
  
  #     ptm <- proc.time()
  #     fe<-find.encountr(tmp1,tmp1$date,tmp2,tmp2$date,1:nrow(tmp1),1:nrow(tmp2),time.buff=1440,spatial.buff=buff)
  #     proc.time() - ptm
  #     
  #     zz[i,3]<-nrow(fe)
  #   }
  
  ## Calc distances between each point on each unique day
  #   for(i in 1:length(u)) {  
  #     tmp<-x[as.Date(x$date,tz='MST')==u[i],]
  #     tmp<-tmp[order(tmp$AID),]
  #     
  #     dists<-as.matrix(dist(cbind(tmp$x,tmp$y)))
  #     dists[dists<buff]<-1
  #     dists[dists>buff]<-0
  #     id<-unique(tmp$AID)
  #     
  #     rownames(dists)<-colnames(dists)<-id
  #     cols<-match(zz[,1],colnames(dists))
  #     rows<-match(zz[,2],rownames(dists))
  #     g<-cbind(cols,rows)
  #     add<-dists[g]
  #     add[is.na(add)==TRUE]<-0
  #     zz[,3]<-zz[,3]+add
  #   }
  return(zz)
}
