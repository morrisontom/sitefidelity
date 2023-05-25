wind<-90
title<-"FG_EK"

setwd("C:/Users/tmorri19/Documents/GitHub/SiteFidelityProj/data/",title,"/",sep=""))  
dat<-read.csv("BFH_elk_database_complete_smp.csv")
u<-unique(dat$pop)

setwd(paste("C:/Users/tmorri19/Documents/GitHub/SiteFidelityProj/data/",title,"/",title,"_run",sep=""))
source('iyd.R')

for (i in 1:length(u)){
  obs<-read.csv(paste(title,"_pop",u[i],"_obs.csv",sep=""))
  exp<-read.csv(paste(title,"_pop",u[i],"_exp.csv",sep=""))
  dates<-read.csv(paste(title,"_pop",u[i],"_dates.csv",sep=""))

  obs<-obs[,-1]
  exp<-exp[,-1]
  dates<-dates[,-1]
  dates<-as.Date(dates)

  f<-NULL
  comp<-list(dates,obs,exp)

  jpeg(paste(title,u[i],".jpg",sep=""),width = 800, height = 800)
  iyd.plot(comp,wind,title)    # plots in 2x2 plot all iyd data
  dev.off()
}