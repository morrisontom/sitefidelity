

gen.traj<-function(loc){
  loc<-as.data.frame(loc)
  traj2<-as.ltraj(xy=loc[,c("x","y")],date=loc$date,id=loc$AID)
  return(traj2)
}

### Fit normal and wrapped-cauchy dist to data for each individual
fitwc<-function(xx){
  xx$AID<-as.character(xx$AID)
  traj2<-gen.traj(xx)
  dofit<-do.call(rbind,lapply(1:length(traj2),function(p){
  
  for(p in 1:length(traj2)){    
  # for(p in 1:length(traj2)){
    lens<-traj2[[p]]$dist
    lens<-lens[!is.na(lens)]
    ang<-traj2[[p]]$rel.angle  
    ang<-ang[!is.na(ang)]
    lens[lens==0]<-0.001
    ne<-fitdistr((lens),densfun="lognormal")  #weibull

    mu0 <- circ.mean(ang)
    rho0 <- est.rho(ang)
    wc<-wrpcauchy.ml(ang, mu0, rho0)
  }
    
    return(data.frame(AID=p,mean=(ne[[1]][[1]]),sd=(ne[[1]][[2]]),mu=wc[[1]],rho=wc[[2]])) 
  }))
}

### For each point, sample 10 nearby points
gen_avail<-function(x){ 
  nsim<-10
  x<-x[!is.na(x$AID),]
  datdist<-fitwc(x) ## function above 
  dat4<-as.data.frame(x) 
  u<-unique(dat4$AID) 
  #for(k in 1:length(u)){
  rloc<-do.call(rbind,lapply(1:length(u),function(k){
    tmp<-dat4[dat4$AID==u[k],]
    #for(p in 1:length(unique(tmp$date))){
    rn<-do.call(rbind,lapply(1:length(unique(tmp$date)),function(p){
      tmp2<-tmp[tmp$date==unique(tmp$date)[p],]
      l<-exp(rnorm(nsim*10,datdist[k,2],datdist[k,3]))  #multiplied by 10 to generate enough points that fall more than 250m away
      l<-l[l>250]  # make sure the new point is not same modis pixel as before
      l<-l[1:nsim]
      theta<-runif(nsim,0,360)
      theta<-(theta*pi)/180
      newx<-tmp2$x-cos(theta)*l
      newy<-tmp2$y-sin(theta)*l

      return(data.frame(x=newx,y=newy,date=unique(tmp$date)[p]))
    }))
    return(data.frame(AID=u[k],rn))
  }))

}


