
######## Functions to aggregate data by season

season<-function(z,dts,seas,title) {
  uAID<-unique(z[[6]])
  sfid<-sobs<-NULL
  all2<-NULL
  
  for(k in 1:length(uAID)){
    istrue<-NA
    inddates<-c(as.Date(z[[7]][1,k]),z[[7]][2,k])
    alldates<-z[[1]][z[[1]]>=inddates[1] & z[[1]]<=inddates[2]]  ## ALL DATES FOR THE INDIVIUDAL
    
    sd<-min(alldates,na.rm=TRUE)+(365)                                            ## YEAR 2 sd and ed
    ed<-max(alldates,na.rm=TRUE)
    yr<-year(ed)
    sobs<-z[[2]][z[[1]]>=sd & z[[1]]<=ed,k]
    sfid<-z[[4]][z[[1]]>=sd & z[[1]]<=ed,k]
    alldates<-alldates[alldates>=sd & alldates<=ed]
    all<-data.frame(AID=NA,sfid=NA,sobs=NA,nobs=NA,season=NA,year=NA,title=NA)
    
    if(any(yday(alldates)>dts[1] & yday(alldates)<dts[2])){
      sobs<-sobs[yday(alldates)>dts[1] & yday(alldates)<dts[2]]
      sfid<-sfid[yday(alldates)>dts[1] & yday(alldates)<dts[2]]
      sobs<-sobs[!is.na(sobs)]
      sfid<-sfid[!is.na(sfid)]
      
      istrue<-unlist(lapply(sfid,function(x){length(x[is.na(x)==FALSE])}))
      istrue<-length(istrue[istrue>0])
      if(istrue==0) {
        istrue<-NULL
      }else{
        all<-data.frame(AID=z[[6]][k],sfid=mean(sfid,na.rm=TRUE),sobs=mean(sobs,na.rm=TRUE),
                        nobs=istrue,season=seas,year=yr,title=title)
      }
    }
    all2[[k]]<-all
  }
  return(all2)
} # END FUNCTION

season.annual<-function(z,seas,title) {
  #uy<-unique(as.numeric(strftime(z[[1]], format = "%Y", tz = tz)))
  
  uAID<-unique(z[[6]])
  sfid<-sobs<-NULL
  all2<-NULL
  
  for(k in 1:length(uAID)){
    istrue<-NULL
    inddates<-c(as.Date(z[[7]][1,k]),z[[7]][2,k])
    alldates<-z[[1]][z[[1]]>=inddates[1] & z[[1]]<=inddates[2]]  ## ALL DATES FOR THE INDIVIUDAL
    
    sd<-min(alldates,na.rm=TRUE)+(365)                                            ## YEAR 2 sd and ed
    ed<-max(alldates,na.rm=TRUE)
    yr<-year(ed)
    sobs<-z[[2]][z[[1]]>sd & z[[1]]<ed,k]
    sfid<-z[[4]][z[[1]]>sd & z[[1]]<ed,k]
    alldates<-alldates[alldates>=sd & alldates<=ed]
    
    all<-data.frame(AID=NA,sfid=NA,sobs=NA,nobs=NA,season=NA,year=NA,title=NA)
    if(any((alldates)>(sd) & (alldates)<(ed))){
      sobs<-sobs[(alldates)>(sd) & (alldates)<(ed)]
      sfid<-sfid[(alldates)>(sd) & (alldates)<(ed)]
      sobs<-sobs[!is.na(sobs)]
      sfid<-sfid[!is.na(sfid)]
      
      istrue<-unlist(lapply(sfid,function(x){length(x[is.na(x)==FALSE])}))
      istrue<-length(istrue[istrue>0])
      if(istrue==0) {
        istrue<-NULL
        all<-data.frame(AID=z[[6]][k],sfid=NA,sobs=mean(sobs,na.rm=TRUE),
                        nobs=NA,season=seas,year=NA,title=title)
        }else{
          all<-data.frame(AID=z[[6]][k],sfid=mean(sfid,na.rm=TRUE),sobs=mean(sobs,na.rm=TRUE),
                          nobs=istrue,season=seas,year=yr,title=title)
        }
    }

    all2[[k]]<-all
  }  
  return(all2)
}

compseason<-function(x,dt,title,rr){
  #WINTER FIDELITY - must have at least 30 days of data
  dtwi <-c(dt[match(substr(title,1,2),dt[,1]),6],dt[match(substr(title,1,2),dt[,1]),7])
  fidwin<-do.call("rbind", season(x,dtwi,"winter",title))
  
  #SPRING FIDELITY
  dtsp<-c(dt[match(substr(title,1,2),dt[,1]),8],dt[match(substr(title,1,2),dt[,1]),9])
  fidspr<-do.call("rbind", season(x,dtsp,"spring",title))
  
  #PARTURITION FIDELITY
  part <-c(dt[match(substr(title,1,2),dt[,1]),3]-rr,dt[match(substr(title,1,2),dt[,1]),3]+rr)
  fidpart<-do.call("rbind", season(x,part,"rut",title))
  
  #SUMMER FIDELITY - must have at least 30 days of data
  dtsu <-c(dt[match(substr(title,1,2),dt[,1]),4],dt[match(substr(title,1,2),dt[,1]),5])
  fidsum<-do.call("rbind", season(x,dtsu,"summer",title))
  
  #RUT FIDELITY
  rut<-c(dt[match(substr(title,1,2),dt[,1]),2]-rr,dt[match(substr(title,1,2),dt[,1]),2]+rr)
  fidrut<-do.call("rbind", season(x,rut,"rut",title))
  
  #ALL FIDELITY
  fid<-do.call("rbind", season.annual(x,"annual",title))
  
  fidelity<-cbind(fid,fidwin[,2:4],fidspr[,2:4],fidpart[,2:4],fidsum[,2:4],fidrut[,2:4])
  #fidelity<-fidelity[!fidelity$year==min(fidelity$year),]
  names(fidelity)<-c("AID","fid","lobs", "nobs", "season", "year","title","fidwin","lobswin","nobswin",
                     "fidspr","lobsspr","nobsspr","fidpart","lobspart","nobspart",
                     "fidsum","lobssum","nobssum","fidrut","lobsrut","nobsrut")
  return(fidelity)
}

############################
### Plot seasonal fidelity
############################
season.plot<-function(z,dts,seas,title,critN,colr) {
  uy<-unique(as.numeric(strftime(z[[1]], format = "%Y", tz = tz)))
  all2<-NULL
  for(k in 1:length(uy)){
    alldates2<-z[[1]][as.numeric(strftime(z[[1]], format = "%Y", tz = tz))==uy[k]]
    if(any(as.numeric(strftime(alldates2, format = "%j", tz = tz))<dts[2] & as.numeric(strftime(alldates2, format = "%j", tz = tz))>dts[1])){
      sx<-c(min(alldates2[as.numeric(strftime(alldates2, format = "%j", tz = tz))<dts[2] & 
                            as.numeric(strftime(alldates2, format = "%j", tz = tz))>dts[1]],na.rm=TRUE),
            max(alldates2[as.numeric(strftime(alldates2, format = "%j", tz = tz))<dts[2] & 
                            as.numeric(strftime(alldates2, format = "%j", tz = tz))>dts[1]],na.rm=TRUE))
      sy<-c(-2,yy[2]+10000,yy[2]+10000,-2)
      polygon(c(rep(sx[1],2),rep(sx[2],2)),sy,col=makeTransparent(colr),border=makeTransparent(colr)) ### change alpha
    }
  }
} # END FUNCTION

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