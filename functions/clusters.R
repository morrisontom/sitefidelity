## Classifies populations by startiung locations in year 1 using kmeans. Clusters with fewer than... 
##critical pop size (cps) individual are thrown out and loop starts again until clusters stabilize
## T. Morrison
## May 2, 2015

clst <- function(x,cps){
  t<-0
  repeat{
    ncl<-pamk(x[,5:6])[[2]]  #identify no. clusters
    cl<-kmeans(x[,5:6],ncl)  #group by clusters
    cs<-cl$cluster
    x<-x[order(cs),]
    cs<-cs[order(cs)]
    
    cdist<-NULL
    for (i in 1:ncl) {
      if(length(cs[cs==i])<cps) {
        x<-x[cs!=i,]
        cs<-cs[cs!=i]
      }
      else {
        cdist<-c(cdist,dist(rbind(cl[[2]][i,],x[cs==i,5:6]))[1:table(cs)[[i]]]/1000)
      }  
    }
    x<-x[cdist<centdist,]
    cs<-cs[cdist<centdist]
    t<-t+1
    
    if(min(table(cs))>5 & max(cdist)<centdist){ 
      break
    }
  }
  x<-cbind(x,cs)
  x<-data.frame(x)
  colnames(x)<-c("AID","START","END","DAYS","START_X","START_Y","pop")
  x$START<-as.POSIXct(x$START,origin="1970-01-01",tz="MST")
  x$END<-as.POSIXct(x$START,origin="1970-01-01",tz="MST")
  
  return(x)
}
