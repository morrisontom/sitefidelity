rm_crit<-function(x,crit){
  zz<-dat2<-NULL
  for (j in 1:length(unique(x$AID)) ) {
    tmp <- x[x$AID==unique(x$AID)[j],]
    if ((max(tmp$date,na.rm=TRUE) - min(tmp$date,na.rm=TRUE) > crit) == TRUE ) {
      yy <- data.frame(AID=unique(x$AID)[j], MINDATE=as.numeric(min(tmp$date,na.rm=TRUE)),
                       MAXDATE=as.numeric(max(tmp$date,na.rm=TRUE)),TOTDATE=max(tmp$date,na.rm=TRUE) - min(tmp$date,na.rm=TRUE),
                       xMINDATE=tmp$x[tmp$date==min(tmp$date,na.rm=TRUE)][1],yMINDATE=tmp$y[tmp$date==min(tmp$date,na.rm=TRUE)][1])
      zz <- as.data.frame(rbind(zz,yy))
      tmp<-tmp[is.numeric(tmp$x),]
      tmp<-tmp[is.numeric(tmp$y),]
      dat2 <- rbind(dat2,tmp)
    }
  }
  return(dat2)
}