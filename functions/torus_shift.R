
TorusShift<-function(xy,                  #spatial data with id in xy
                     id.identifier="AID", #name of column containing animal id 
                     SpringLength,        #spring length stack
                     INDVI,               #INTREGRATED ndvi stack
                     buffer=500,          #Buffer of study area (in same units as projection)
                     nshifts=250){    #Number of times to rotate the points              
   
  crs.xy<-xy@proj4string
  xy.df<-as.data.frame(xy) 
  
  #get extent for torus
  bb<-bbox(xy)
  bb[, 1]<-bb[, 1]-buffer
  bb[, 2]<-bb[, 2]+buffer
  
  window.ppp<-owin(bb[1, ], bb[2, ])
  #rename id and time column for easy manipulation
  #deal with names...not sure if this is an issue with your data 
  n<-names(xy.df)
  id.idx<-which(n==id.identifier)    
  n[id.idx]<-"ID"
  names(xy.df)<-n

  ids<-unique(xy.df$ID)
  
  summary.df<-data.frame(ID=ids[1], iteration=1, SpringLength=0, INDVI_CV=0, INDVI=0)[0, ]
  for(i in 1:length(ids)){
    xy.i<-xy.df[xy.df$ID==ids[i], ]
    xy.ppp<-ppp(x=xy.i$x, y=xy.i$y, bb[1, ], bb[2, ])
    
    xy.shift.sp<-as.SpatialPoints.ppp(xy.ppp)
    proj4string(xy.shift.sp)<-crs.xy
    
    shift.indvi<-extract(INDVI, xy.shift.sp)
    
    shift.indvi.CV<-mean(apply(shift.indvi, MARGIN=1, FUN=sd, na.rm=T)/apply(shift.indvi, MARGIN=1, FUN=mean, na.rm=T), na.rm=T)  #NEW
    shift.indvi<-mean(apply(shift.indvi, MARGIN=1, FUN=mean, na.rm=T), na.rm=T)
    
    shift.sl<-extract(SpringLength, xy.shift.sp)
    shift.sl<-mean(apply(shift.sl, MARGIN=1, FUN=mean, na.rm=T), na.rm=T)
    s.df<-data.frame(ID=ids[i], iteration=0, SpringLength=shift.sl, INDVI_CV=shift.indvi.CV, INDVI=shift.indvi)
    summary.df<-rbind(summary.df, s.df)
    
    for(n in 1:nshifts){
      xy.shift<-rshift(xy.ppp, edge="torus")
      xy.shift.sp<-as.SpatialPoints.ppp(xy.shift)
      proj4string(xy.shift.sp)<-crs.xy
      
      shift.indvi<-extract(INDVI, xy.shift.sp)
      
      shift.indvi.CV<-mean(apply(shift.indvi, MARGIN=1, FUN=sd, na.rm=T)/apply(shift.indvi, MARGIN=1, FUN=mean, na.rm=T), na.rm=T)  #NEW
      shift.indvi<-mean(apply(shift.indvi, MARGIN=1, FUN=mean, na.rm=T), na.rm=T)
      shift.sl<-extract(SpringLength, xy.shift.sp)
      shift.sl<-mean(apply(shift.sl, MARGIN=1, FUN=mean, na.rm=T), na.rm=T)
      s.df<-data.frame(ID=ids[i], iteration=n, SpringLength=shift.sl, INDVI_CV=shift.indvi.CV, INDVI=shift.indvi)
      summary.df<-rbind(summary.df, s.df)
    }
  }
  return(summary.df)
}



