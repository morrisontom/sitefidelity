
########################################
### Calculate Weibull
########################################
calc.weibull<-function(z){
  shape<-scale<-NULL  
  for (j in 1:ncol(z)) {
    xx<-z[,j]
    xx[xx==-Inf]<-NA
    xx<-xx[is.na(xx)==FALSE]
    if(length(xx)>20){
      wd<-fitdist(xx,'weibull')
      shape<-rbind(shape,cbind(wd$estimate[1],wd$sd[1]))
      scale<-rbind(scale,cbind(wd$estimate[2],wd$sd[2]))
    }else {
      shape<-rbind(shape,cbind(NA,NA))
      scale<-rbind(scale,cbind(NA,NA))
    }  
  }
  parm<-cbind(shape,scale)
  colnames(parm)<-c("shape","shape.sd","scale","scale.sd")
  return(parm)
}

########################################
### Calculate normal
########################################
calc.normal<-function(z){
  mean1<-sd1<-NULL  
  for (j in 1:ncol(z)) {
    xx<-z[,j]
    xx[xx==-Inf]<-NA
    xx<-xx[is.na(xx)==FALSE]
    if(length(xx)>20){
      nd<-fitdistr(xx,'normal')
      mean1<-rbind(mean1,nd$estimate[[1]])
      sd1<-rbind(sd1,nd$estimate[[2]])
    }else {
      mean1<-rbind(mean1,NA)
      sd1<-rbind(sd1,NA)
    }  
  }
  parm<-cbind(mean1,sd1)
  colnames(parm)<-c("mean","sd")
  return(parm)
}


#########################################
# Calculate min IYD's by individual
#########################################

iyd<-function(x,    # df of GPS data with colnames: AID='animalID',date in POSIX, x,y in UTM)
              wind, # 
              uproj
              ) {  
  
  x <- as.data.frame(x)
  u<-unique(x$AID)
  allanim <- NULL
  for(i in 1:length(u)){
    a1<-x[x$AID==u[i],]
    
    # create new window sizes for each possible window
    windownames <- paste0('iyd',c('REAL',wind[2:length(wind)]))
    a1[,windownames] <- NA
    
    # identify max/min dates
    date2 <- as.Date(a1$date)
    startdt <- min(date2)
    enddt <- max(date2)
    
    if(as.numeric(enddt-startdt) > (365+min(wind))){
      
      for(p in 1:length(wind)){
        startfid <- startdt + wind[p] + 365
        fidcalcdt <- which(date2>=startfid)
        
        for(j in fidcalcdt){
          dt1 <- date2[j]
          focal <- a1[j,]
          
          # points around focal in year 1
          y1 <- a1[date2>=(dt1-365-(wind[p]/2)) & 
                     date2<=(dt1-365+(wind[p]/2)),]
          
          # points around focal in year 2
          y2 <- a1[date2>=(dt1-(wind[p]/2)) & 
                     date2<=(dt1+(wind[p]/2)),]
  
          # only calc IYD for individuals with at least 80% of points in y1
          if(nrow(y1) > (wind[p]*0.8)){
            
            # calculate min distance using utm - thanks Pythagorean
            a1[j,windownames[p]] <- min(sqrt((focal$y-y1$y)^2+(focal$x-y1$x)^2))
            
            # calculate home range overlap
            y1$yr <- 'yr1'
            y2$yr <- 'yr2'
          }else{
            a1[j,windownames[p]] <- NA
          }
          
          # if(nrow(y1) > 5 & nrow(y2) > 5){
          #   aover <- rbind(y1,y2)
          #   coordinates(aover) <- ~x+y
          #   proj4string(aover) <- uproj
          #   aover <- aover[,'yr']
          #   kd.yr <- kerneloverlap(aover,method = 'PHR',percent=95)
          #   a1[j,windownames.overlap[p]] <- kd.yr[2,1]
          # }else{
          #   a1[j,windownames.overlap[p]] <- NA
          # }
          
          # RSF of distance to previous points
          # coordinates(y1) <- ~x+y
          # proj4string(y1) <- uproj
          # coordinates(focal) <- ~x+y
          # proj4string(focal) <- uproj
          # 
          # xproj <- a1
          # coordinates(xproj) <- ~x+y
          # proj4string(xproj) <- uproj
          # 
          # ext <- extent(xproj)
          # cell.size <- 100
          # emptygrid <- raster::raster(
          #   nrow=(ext@xmax-ext@xmin)/cell.size, 
          #   ncol=(ext@ymax-ext@ymin)/cell.size,
          #   xmn=ext@xmin,
          #   xmx=ext@xmax,
          #   ymn=ext@ymin,
          #   ymx=ext@ymax,
          #   crs=crs(uproj),
          #   resolution=cell.size,
          #   vals=0)
          # 
          # mask by home range
          # dist.yr1 <- dist_raster(y1,r=emptygrid)
          # plot(uproj,add=T,pch=16)
        }
      }
    }
    # allanim <- rbind(allanim,a1[,c(windownames,windownames.overlap)])
    allanim <- rbind(allanim,a1[,c(windownames)])
  }
  return(allanim)
}        

dist_raster <- function(
  s1,   # spatial points/lines/polygons
  r){   # empty raster into which distances will be calcuated
  
  require(rgeos)
  
  dd = gDistance(s1, as(r,"SpatialPoints"), byid=TRUE)
  r[] <- apply(dd,1,min)
  return(r)
} #End function


##################
### Plot Min IYD's
##################
iyd.plot<-function(z,wind,title) {     #z = compall
  par(mfrow=c(2,2),mar=c(3,5,1,1), oma=c(0,0,2,1))
  cg<-c("red","blue","black","green","orange","purple","yellow","lightblue","lightred","lightgreen","darkgrey")  
  
  z[[1]]<-z[[1]][order(z[[1]])]
  z[[2]]<-z[[2]][order(z[[1]]),]
  z[[3]]<-z[[3]][order(z[[1]]),]
  alldates<-z[[1]]
  sd<-min(alldates,na.rm=TRUE)+365+floor(wind/2)
  sd<-"2014-05-01"
  mcdt<-alldates[sd<=alldates]
  ed<-max(alldates)
  z[[2]][z=="-Inf"] <- NA
  z[[2]][z=="Inf"] <- NA
  z[[3]][z=="-Inf"] <- NA
  z[[3]][z=="Inf"] <- NA
  
  ## DISTANCES
  for (j in 1:ncol(z[[2]])) {      
    ci<-z[[2]][,j]
    if(j==1) {plot(alldates,exp(ci),col="lightgrey",type="l",lwd=0.6,xlim=c(as.Date(sd,origin="1970-01-01"),ed),
                ylim=c(0,exp(max(z[[2]],na.rm=TRUE))*4),main="Distances",xlab="Date",ylab="Min IYD (m)",
                xaxp=c(min(mcdt),max(mcdt),24))
    }else{
    lines(alldates,exp(ci),type="l",lwd=0.6,col=c("lightgrey",alpha=0.5))
    }
  }

  cp<-z[[2]]  
  cp[cp==-Inf] <- NA
  cp<-(rowMeans(cp,na.rm=TRUE))
  lines(alldates,exp(cp),lwd=5,col=cg[1])  
  
  np<-z[[3]]
  np[np==-Inf] <- NA
  np<-(rowMeans(np,na.rm=TRUE))
  lines(alldates,exp(np),lwd=5,col=cg[2])  

  ## LOG DISTANCES
  for (j in 1:ncol(z[[2]])) {      
    ci<-z[[2]][,j]
    if(j==1) {plot(alldates,(ci),col="lightgrey",type="l",lwd=0.6,xlim=c(as.Date(sd,origin="1970-01-01"),ed),
                   ylim=c(0,max(z[[2]],na.rm=TRUE)*4),main="LOG(Distances)",xlab="Date",ylab="Min IYD (m)",
                   xaxp=c(min(mcdt),max(mcdt),24))
    }else{
      lines(alldates,(ci),type="l",lwd=0.6,col=c("lightgrey",alpha=0.5))
    }
  }
  
  cp<-z[[2]]
  cp[cp==-Inf]<-NA
  cp<-(rowMeans(cp,na.rm=TRUE))
  lines(alldates,(cp),lwd=2.5,col=cg[1])  
  
  np<-z[[3]]
  np[np==-Inf]<-NA
  np<-(rowMeans(np,na.rm=TRUE))
  lines(alldates,(np),lwd=2.5,col=cg[2])  
  
  ################
  ## CDF DISTANCES
  ################
  
  ## CDF
  mcall<-NULL
  tt<-0
  for(j in 1:ncol(z[[2]])){
    ci<-z[[2]][j]
    ci<-ci[is.na(ci)==FALSE]

    if(length(ci)>0){
      tt<-tt+1
      ci<-exp(ci[order(ci)])
      ci<-cbind(ci,(1:length(ci))/length(ci))  
      if(tt==1) {plot(ci,col="lightgrey",type="l",lwd=0.6,xlim=c(0,max(exp(z[[2]]),na.rm=TRUE)),
                   ylab="Proportion of IYDs",xlab="Minimum Interyear-Distances (m)")
      }else{lines(ci,col="lightgrey",type="l",lwd=0.6)}
      mcall<-c(mcall,ci[,1])    
    }
  }  
  mcall<-mcall[order(mcall)]
  mcall<-cbind(mcall,(1:length(mcall))/length(mcall))
  lines(mcall,col=c(cg[1],alpha=0.2),lwd=2.5)
  lines(c(0,max(mcall)+5),c(1,1),lty=2,col=c("grey"))
  
  ## LOG CDF
  mcall<-NULL
  tt<-0
  for(j in 1:ncol(z[[2]])){
    ci<-z[[2]][j]
    ci<-ci[is.na(ci)==FALSE]
    
    if(length(ci)>0){
      tt<-tt+1
      ci<-(ci[order(ci)])
      ci<-cbind(ci,(1:length(ci))/length(ci))
      if(tt==1) {plot(ci,col="lightgrey",type="l",lwd=0.6,xlim=c(0,max(z[[2]],na.rm=TRUE)),
                   ylab="Proportion of IYDs",xlab="Minimum Interyear-Distances (m)")
      }else{lines(ci,col="lightgrey",type="l",lwd=0.6)}
      mcall<-c(mcall,ci[,1])    
    }  
  }
  mcall<-mcall[order(mcall)]
  mcall<-cbind(mcall,(1:length(mcall))/length(mcall))
  lines(mcall,col=c(cg[1],alpha=0.2),lwd=2.5)
  lines(c(0,max(mcall)+5),c(1,1),lty=2,col=c("grey"))
  mtext(paste(title,"n =",ncol(z[[2]])), side=3, line=1, outer=TRUE, cex=1, font=2)
} # END FUNCTION


##################
### Plot log IYD's
##################
iyd.plot.log<-function(z,wind,title,irgdt,crit) { 
  cg<-c("red","blue","black","green","orange","purple","yellow","lightblue","lightred","lightgreen","darkgrey")  
  
  ## reorder all datasets by date
  z[[1]]<-z[[1]][order(z[[1]])]
  z[[2]]<-z[[2]][order(z[[1]]),]
  z[[3]]<-z[[3]][order(z[[1]]),]
  z[[4]]<-z[[4]][order(z[[1]]),]
  
  z[[2]][z=="-Inf"] <- NA
  z[[2]][z=="Inf"] <- NA
  z[[3]][z=="-Inf"] <- NA
  z[[3]][z=="Inf"] <- NA
  z[[4]][z=="-Inf"] <- NA
  z[[4]][z=="Inf"] <- NA
  
  alldates<-z[[1]]
  sd<-min(alldates,na.rm=TRUE)+365+floor(wind/2)
  mcdt<-alldates[sd<=alldates]
  
  nrows <- (apply( z[[2]], 1, function(f) length(na.omit(f))))
  alldates<-alldates[nrows>2]
  z[[1]]<-z[[1]][nrows>2]
  z[[2]]<-z[[2]][nrows>2,]
  z[[3]]<-z[[3]][nrows>2,]
  z[[4]]<-z[[4]][nrows>2,]
  ed<-max(alldates)
  yy<-c(0,50000)
  
  ## LOG DISTANCES
  cp<-z[[2]]
  cp[cp==-Inf]<-NA
  cp<-(rowMeans(cp,na.rm=TRUE))
  plot(alldates,(cp),type="l",lwd=2.5,col=cg[1],xlim=c(as.Date(sd,origin="1970-01-01"),ed),
       ylim=yy,xlab="Date",ylab="Min IYD (m)",
       xaxp=c(min(mcdt),max(mcdt),24),main=(paste(title,"n =",ncol(z[[2]]))))  
  
  np<-z[[3]]
  np[np==-Inf]<-NA
  np<-(rowMeans(np,na.rm=TRUE))
  lines(alldates,(np),lwd=2.5,col=cg[2])  
  
  for (j in 1:ncol(z[[2]])) {      
    ci<-z[[2]][,j]
    if(j==1) {lines(alldates,(ci),col="lightgrey",type="l",lwd=0.6)
    }else{
      lines(alldates,(ci),type="l",lwd=0.6,col=c("lightgrey",alpha=0.5))
    }
  }
  lines(alldates,(cp),lwd=2.5,col=cg[1])  
  lines(alldates,(np),lwd=2.5,col=cg[2])

  dtsp<-c(irgdt-15,irgdt+15)
  season.plot(z,dtsp,"spring",title,crit,"green")
  
  dtsu<-c(182,246)
  season.plot(z,dtsu,"summer",title,crit,"red")
  
  dtwi<-c(1,60)
  season.plot(z,dtwi,"winter",title,crit,"blue")
  
}  # END FUNCTION


############################
### Plot seasonal fidelity
############################
pop.plot<-function(z,i,j,title) {
 
  if(substr(title,1,2)=="EK") colr=rainbow(5)[1]
  if(substr(title,1,2)=="BS") colr=rainbow(5)[2]
  if(substr(title,1,2)=="PH") colr=rainbow(5)[3]
  if(substr(title,1,2)=="MS") colr=rainbow(5)[4]
  if(substr(title,1,2)=="MD") colr=rainbow(5)[5]
      
  z[[2]][z=="Inf"] <- NA
  nrows <- (apply( z[[2]], 1, function(f) length(na.omit(f))))
  z[[1]]<-z[[1]][nrows>2]
  z[[2]]<-z[[2]][nrows>2,]
  alldates<-z[[1]]
  alldates<-as.numeric(strftime(alldates, format = "%j", tz = tz))
  alldates<-alldates[order(alldates)]
  z[[2]]<-z[[2]][order(alldates),]
  yy<-c(0,35)
  
  ## DISTANCES
  cp<-z[[2]]
  cp<-rowMeans(cp,na.rm=TRUE)
  #cp<-log(cp)
  if(i==1){
  plot(alldates,(cp),type="l",col=(colr),xlim=c(0,365),
       ylim=yy,xlab="Date",ylab="Fidelity",
       xaxp=c(min(mcdt),max(mcdt),24),main=substr(title,1,2))
  }else{
    lines(alldates,(cp),col=(colr)  )
  }
} # END FUNCTION


############################
### Make colors transparent
############################

makeTransparent<-function(someColor, alpha=25)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
} #ENd function
