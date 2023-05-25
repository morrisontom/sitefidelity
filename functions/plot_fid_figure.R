## IYD over time (FIG3)
## updated Jul 30 2020

rm(list=ls(all=TRUE))
require(lme4)
require(lubridate)
require(sp)
require(rgdal)
require(rgeos)

### GENERATE FIDELITY PLOTS
setwd("C:/users/administrator/OneDrive - University of Glasgow/Documents/Side projects/Fidelity/ungulatefidelity/functions")
source('iyd.R')
source('fig_fid_overtime.R')

setwd("C:/users/administrator/OneDrive - University of Glasgow/Documents/Side projects/Fidelity/ungulatefidelity/output/output_2020-07-21/")
results <- read.csv("fidelity_all.csv")
results$date <- as.POSIXct(results$date,format="%Y-%m-%d %H:%M:%S")
results <- results[!is.na(results$iydREAL) & !is.na(results$date),]
results$weeks <- week(results$date)

# remove resident wb - not sure about this
setwd('C:/Users/Administrator/OneDrive - University of Glasgow/Serengeti Wildebeest Project/Tracking Database')
wb <- read.csv('WB_ZB_EL_SNP_DATA_19990225-2020-05-20_Deads_deleted.csv')
mwb <- unique(wb$AID[wb$migrant==1 & wb$SPECIES=='WB'])
mwb <- paste0('WB_SE_',mwb)
x <- results[results$spp=='WB',]
results <- results[results$spp!='WB',]
x <- x[x$AID %in% mwb,]
results <- rbind(results,x)

# add homerange
setwd("C:/Users/Administrator/OneDrive - University of Glasgow/Documents/Side projects/Fidelity/ungulatefidelity/output/output_2020-07-21/")
hr <- read.csv('HRsqkm_all_.csv')
hr$spp <- substr(hr$AID,1,2)
hr$title <- substr(hr$AID,4,5)

# some home ranges are clearly not correct, so set to NA if below log-1.5 = 4 kmsq
hr$hr[which(log(hr$hr)<1.5)] <- NA
results$hr <- hr$hr[match(results$AID,hr$AID)]
results$hr_log <- log(results$hr)

# unique species
uSPP <- unique(results$spp)

res <- results %>%
  group_by(AID,weeks) %>%
  summarise(iydREAL=mean(iydREAL,na.rm=T),
            hr_log = hr_log[1])
res$spp <- substr(res$AID,1,2)
res$title <- substr(res$AID,1,5)

# fit weekly model to iydREAL
for(j in 1:length(uSPP)){
  
  obs <- res[res$spp==uSPP[j],]
  obs <- obs[obs$weeks!=53,]
  obs$iyd_log <- log(obs$iydREAL)
  
  if(!uSPP[j] %in% c("CA","WB","ZB")){
    obsmod<-lmer(iyd_log ~ hr_log + as.factor(weeks) + (1|title/AID),data=obs)
    
  }else{
    obsmod<-lmer(iyd_log ~ hr_log + as.factor(weeks) + (1|AID),data=obs)
  }
  # summary(obsmod)
  # predicted effects
  mySumm <- function(.) { s <- sigma(.)
  c(beta =getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) }
  (t0 <- mySumm(obsmod)) 
  bs <- bootMer(obsmod,FUN = mySumm,nsim=1000)
  
  dfb <- as.data.frame(bs)
  dfb <- data.frame(mean = apply(dfb,2,mean),
                    se = apply(dfb,2,sd)) # calc se
  dfb <- dfb[1:52,]
  dfb$mean[2:nrow(dfb)] <- dfb$mean[1] + dfb$mean[2:nrow(dfb)]
  dfb$l95 <- dfb$mean - (dfb$se*1.96)
  dfb$u95 <- dfb$mean + (dfb$se*1.96)
  dfb$spp <- uSPP[j]  
  
  if(j==1) {
    fn <- dfb
  }else{
    fn <- rbind(fn,dfb)
  }
  print(j)
}
  
#### FIGURE 3
setwd("C:/users/administrator/OneDrive - University of Glasgow/Documents/Side projects/Fidelity/MS/J Anim Ecol/Submit_Jun20//")
fig3(fn,dtrut)
write.csv(fn, 'seasonal_fidelity.csv',row.names = T)

# pdf(paste("FIG3b_",Sys.Date(),".pdf",sep=""))
doh<-fig3b(fn,on,dtrut,results) 
# dev.off()



#############
## TABLE 2
#############

# fn2<-fn[fn$spp!="ZB",]
# on2<-on[on$spp!="ZB",]
on2<-on
fn2<-fn
on2$rn<-as.numeric(substr(on2$rn,17,nchar(on2$rn)))
fn2$rn<-as.numeric(substr(fn2$rn,17,nchar(fn2$rn)))

uSPP<-unique(on2$spp)
mxf<-mnf<-mxo<-mno<-NoInd<-NoPop<-wkfx<-wkfn<-wkox<-wkon<-mt<-mt2<-NULL
for(i in 1:length(uSPP)){
  res<-results[results$species==uSPP[i],]
  NoInd<-c(NoInd,nrow(res))
  NoPop<-c(NoPop,length(unique(res$pop)))
  
  tmp<-fn2$Estimate[fn2$spp==uSPP[i]][1]+fn2$Estimate[fn2$spp==uSPP[i]][2:length(fn2$Estimate[fn2$spp==uSPP[i]])]
  wk<-fn2$rn[fn2$spp==uSPP[i]][2:length(fn2$Estimate[fn2$spp==uSPP[i]])]
  wkfx<-c(wkfx,as.character(as.POSIXct(paste(2015,wk[tmp==max(tmp)],1,sep="-"),format="%Y-%U-%u")))
  wkfn<-c(wkfn,as.character(as.POSIXct(paste(2015,wk[tmp==min(tmp)],1,sep="-"),format="%Y-%U-%u")))
  mxf<-c(mxf,max(tmp))
  mnf<-c(mnf,min(tmp))
  mt<-c(mt,mean(tmp))
  
  as.Date(paste(2015,rut,1,sep="-"),format="%Y-%j")
  
  tmp2<-on2$Estimate[on2$spp==uSPP[i]][1]+on2$Estimate[on2$spp==uSPP[i]][2:length(on2$Estimate[on2$spp==uSPP[i]])]
  wk<-on2$rn[on2$spp==uSPP[i]][2:length(on2$Estimate[on2$spp==uSPP[i]])]
  wkox<-c(wkox,as.character(as.POSIXct(paste(2015,wk[tmp2==max(tmp2)],1,sep="-"),format="%Y-%U-%u")))
  wkon<-c(wkon,as.character(as.POSIXct(paste(2015,wk[tmp2==min(tmp2)],1,sep="-"),format="%Y-%U-%u")))
  mxo<-c(mxo,max(tmp2))
  mno<-c(mno,min(tmp2))
  mt2<-c(mt2,mean(tmp2))
  
}
mxo<-exp(mxo)/1000
mno<-exp(mno)/1000

tab2<-data.frame(Species=as.character(c("Odocoileus hemionus","Alces alces","Connochaetes taurinus","Antilocapra americana","Ovis canadensis","Cervus elephas","Rangifer tarandus")),
                 CommonNames=as.character(c("Mule deer","Moose","Wildebeest","Pronghorn","Bighorn sheep","Elk","Caribou")),
                 spp<-uSPP,
                 NoInd=NoInd,NoPop=NoPop,IndivSD=indsd,PopSD=popsd,
                 HomeRange=mhr,
                 MeanFid=meeny[,1],SEFid=meeny[,2],
                 MaxFid=mxf,MinFid=mnf,
                 MaxFiddt=wkfx,MinFiddt=wkfn,
                 MeanObs=exp(mt2)/1000,
                 MaxObs=mxo,MinObs=mno,
                 MaxObsdt=wkox,MinObsdt=wkon
                 )

setwd("C:/users/administrator/OneDrive - University of Glasgow/Documents/Side projects/Fidelity/results/Figures/FINAL/")
write.csv(tab2,"Table2_out.csv",row.names = FALSE)

#########################
### PLOT by species ALL


#########################



par(mfcol=c(2,1),mar=c(1,3,1,0.2))

#### MEAN 50%
of<-aggregate(obsall$mean,by=list(results$species,results$pop),FUN=function(x) mean(x,na.rm=TRUE))
of.2<-aggregate(of[,3],by=list(of[,1]),FUN=function(x) mean(x,na.rm=TRUE))
of.se<-aggregate(of[,3],by=list(of[,1]),FUN=function(x) 1.96*(sd(x,na.rm=TRUE)/sqrt(length(x[!is.na(x)]))))
of.2<-cbind(of.2,of.se[,2])
of.2[,2]<-exp(of.2[,2])/1000
of.2[,3]<-exp(of.2[,3])
of[,3]<-exp(of[,3])/1000
of.2<-of.2[order(of.2[,2]),]
cole<-rainbow(5)
cole<-cole[c(5,4,1,2,3)]
spp<-spp[match(of.2[,1],spp)]

plot(c(1:5),of.2[,2],col="black",pch='-',cex=6,ylim=c(0,3.5),xaxt="n",xlab="",ylab="",cex.axis=1.2)
for(k in 1:length(unique(of.2[,1]))){
  num<-of[of[,1]==spp[k],3]
  rnum<-runif(length(num),-0.04,0.04)
  points(rep(k,length(num))+rnum,num,col=makeTransparent(cole[k],alpha=100),pch=16,cex=1.5)
  points(rep(k,length(num))+rnum,num,cex=1.5)
}
points(c(1:5),of.2[,2],col=makeTransparent("black",alpha=100),pch='-',cex=6,ylim=c(0,45000),xaxt="n",xlab="",ylab="")


obsall<-obsall[order(match(obsall$species,spp)),]
head(obsall)

far<-aggregate(obsall$mean,by=list(obsall$pop),FUN=function(x) mean(x,na.rm=TRUE))
far$species<-substr(far[,1],1,2)

md<-aggregate(obsall$species,by=list(obsall$pop),length)
md$species<-substr(md[,1],1,2)
look<-aggregate(md[,1],by=list(md$species),length)
look<-look[order(match(look[,1],spp)),]

look$ed<-look$st<-NA
look$st[1]<-1
look$ed[1]<-look[1,2]
for(k in 2:nrow(look)){
  look$st[k]<-look$ed[(k-1)]+1
  look$ed[k]<-look[k,2]+look$ed[(k-1)]
}
look<-c(look$st[1]:look$ed[1],(look$st[2]+2):(look$ed[2]+2),(look$st[3]+4):(look$ed[3]+4),
        (look$st[4]+6):(look$ed[4]+6),(look$st[5]+8):(look$ed[5]+8))

obsall<-obsall[order(obsall)]
boxplot((obsall$mean)~obsall$pop,at=look)
length(unique(obsall$pop))

#### MEAN 95%

of<-aggregate(obsall$dd,by=list(results$species,results$pop),FUN=function(x) mean(x,na.rm=TRUE))
of.2<-aggregate(of[,3],by=list(of[,1]),FUN=function(x) mean(x,na.rm=TRUE))
of.se<-aggregate(of[,3],by=list(of[,1]),FUN=function(x) 1.96*(sd(x,na.rm=TRUE)/sqrt(length(x[!is.na(x)]))))
of.2<-cbind(of.2,of.se[,2])
of.2[,2]<-exp(of.2[,2])/1000
of.2[,3]<-exp(of.2[,3])
of[,3]<-exp(of[,3])/1000
of.2<-of.2[match(spp,of.2[,1]),]

plot(c(1:5),of.2[,2],col="black",pch='-',cex=6,ylim=c(0,50),xaxt="n",xlab="",ylab="",cex.axis=1.2)
for(k in 1:length(unique(of.2[,1]))){
  num<-of[of[,1]==spp[k],3]
  rnum<-runif(length(num),-0.04,0.04)
  points(rep(k,length(num))+rnum,num,col=makeTransparent(cole[k],alpha=100),pch=16,cex=1.5)
  points(rep(k,length(num))+rnum,num,cex=1.5)
}
points(c(1:5),of.2[,2],col=makeTransparent("black",alpha=100),pch='-',cex=6,ylim=c(0,45000),xaxt="n",xlab="",ylab="")



### Summarize datasets
of<-aggregate(results$fid,by=list(results$species,results$pop),FUN=function(x) length(x))
osd<-aggregate(results$AID,by=list(results$species,results$pop),FUN=function(x) length(x))

# #setwd("E:/morrison/SiteFidelityProj/")
# #write.csv(ff,"Data summary.csv")
# 
# res2<-cbind(mm,se[,2])
# colnames(res2)<-c("Species","Mean","SE")
# res2<-res2[(order(res2[,2])),]
# nspp<-1:length(unique(results$species))
# plot(nspp,res2[,2],cex=2,ylim=c(0.2,max(results$fid,na.rm=TRUE)),lwd=4,pch=1,xaxt="n",xlab="Species",ylab="Fidelity")
# arrows(nspp,res2[,2],nspp,res2[,2]-res2[,3],lwd=2,angle=0)
# arrows(nspp,res2[,2],nspp,res2[,2]+res2[,3],lwd=2,angle=0)
# 
# points(rep(1,length(results$fid[results$species=="MD"]))+runif(length(results$fid[results$species=="MD"]),-0.1,0.1),
#        results$fid[results$species=="MD"],col="gray50")
# points(rep(2,length(results$fid[results$species=="MS"]))+runif(length(results$fid[results$species=="MS"]),-0.1,0.1),
#        results$fid[results$species=="MS"],col="gray50")
# points(rep(3,length(results$fid[results$species=="EK"]))+runif(length(results$fid[results$species=="EK"]),-0.1,0.1),
#        results$fid[results$species=="EK"],col="gray50")
# points(rep(4,length(results$fid[results$species=="PH"]))+runif(length(results$fid[results$species=="PH"]),-0.1,0.1),
#        results$fid[results$species=="PH"],col="gray50")
# points(rep(5,length(results$fid[results$species=="BS"]))+runif(length(results$fid[results$species=="BS"]),-0.1,0.1),
#        results$fid[results$species=="BS"],col="gray50")
# points(nspp,res2[,2],cex=2,ylim=c(0,1),pch=1,xaxt="n",lwd=4)

#########################
#### PLOT BY SEASON
#########################
ss1<-ss[is.na(ss$sfid)==FALSE,]
ss1$species<-substr(ss1$title,1,2)
ss1$sfid<-1/ss1$sfid
# mm<-aggregate(ss1$sfid,by=list(ss1$pop,ss1$species,ss1$season),FUN=function(x) mean(x,na.rm=TRUE))

par(mfcol=c(2,1),mar=c(0.3,3,1.5,1))

### OBSERVED
# mm<-aggregate(ss1$sobs,by=list(ss1$species,ss1$season),FUN=function(x) mean(x,na.rm=TRUE))
# ci<-aggregate(ss1$sobs,by=list(ss1$species,ss1$season),FUN=function(x) (sd(x,na.rm=TRUE)/sqrt(length(x))))

mm<-aggregate(ss1$sobs,by=list(ss1$pop,ss1$species,ss1$season),FUN=function(x) mean(x,na.rm=TRUE))
colnames(mm)<-c("pop","spp","season","sobs")
mm$sobs<-exp(mm$sobs)
ci<-aggregate(mm$sobs,by=list(mm$spp,mm$season),FUN=function(x) (sd(x,na.rm=TRUE)/sqrt(length(x))))
mm<-aggregate(mm$sobs,by=list(mm$spp,mm$season),FUN=function(x) mean(x,na.rm=TRUE))

res3<-cbind(mm,ci[,3])
colnames(res3)<-c("species","season","mean","CI")
nspp<-c(1:length(unique(res3$species)),1:length(unique(res3$species))+8,1:length(unique(res3$species))+18)
res3<-res3[(order(res3[,2],res3[,3])),]
res3$col<-0
for(i in 1:length(unique(res3[,1]))){
  res3$col[res3[,1]==unique(res3[,1])[i]]<-i
}
cole<-rainbow(length(unique(res3[,1])))[res3$col]
ylim<-c(-250,10000)
plot(nspp,(res3[,3]),cex=2,pch=16,ylim=ylim,
     xaxt="n",xlab="Species",ylab="Fidelity",col=cole,xlim=c(0,max(nspp)+1),cex.axis=1.7) #ylim=c(1,max(res3$mean,na.rm=TRUE)+0.2)
arrows(nspp,res3[,3]-res3[,4],nspp,res3[,3],lwd=2,angle=0)
arrows(nspp,res3[,3]+res3[,4],nspp,res3[,3],lwd=2,angle=0)
points(nspp,(res3[,3]),cex=2,ylim=c(1,max(res3$mean,na.rm=TRUE)+0.2),pch=1,lwd=1,
       xaxt="n",xlab="Species",ylab="Fidelity",col="black")

lines(c(7,7),c(-10000,ylim[2]+1000),lty=2)
lines(c(17,17),c(-10000,ylim[2]+1000),lty=2)
text(nspp-0.5,res3$mean+650,res3$species,col="black",cex=1)
text(3,ylim[2],"Spring",col="black",cex=1.4)
text(12,ylim[2],"Summer",col="black",cex=1.4)
text(21,ylim[2],"Winter",col="black",cex=1.4)

### EXPECTED 

# mm<-aggregate(ss1$sfid,by=list(ss1$species,ss1$season),FUN=function(x) mean(x,na.rm=TRUE))
# ci<-aggregate(ss1$sfid,by=list(ss1$species,ss1$season),FUN=function(x) (sd(x,na.rm=TRUE)/sqrt(length(x))))

mm<-aggregate(ss1$sfid,by=list(ss1$pop,ss1$species,ss1$season),FUN=function(x) mean(x,na.rm=TRUE))
colnames(mm)<-c("pop","spp","season","fid")
ci<-aggregate(mm$fid,by=list(mm$spp,mm$season),FUN=function(x) (sd(x,na.rm=TRUE)/sqrt(length(x))))
mm<-aggregate(mm$fid,by=list(mm$spp,mm$season),FUN=function(x) mean(x,na.rm=TRUE))

res2<-cbind(mm,ci[,3])
colnames(res2)<-c("species","season","mean","CI")
res2<-res2[(order(res2[,2],res2[,3])),]
res2$col<-res3$col[match(res2[,1],res3[,1])]
cole<-rainbow(length(unique(res2[,1])))[res2$col]

plot(nspp,res2[,3],cex=2,ylim=c(1,max(res2$mean,na.rm=TRUE)+0.2),pch=16,
     xaxt="n",xlab="Species",ylab="Fidelity",col=cole,xlim=c(0,max(nspp)+1),cex.axis=1.7)
arrows(nspp,res2[,3],nspp,res2[,3]-res2[,4],lwd=2,angle=0)
arrows(nspp,res2[,3],nspp,res2[,3]+res2[,4],lwd=2,angle=0)
points(nspp,res2[,3],cex=2,ylim=c(1,max(res2$mean,na.rm=TRUE)+0.2),pch=1,lwd=1,
       xaxt="n",xlab="Species",ylab="Fidelity",col="black")
lines(c(7,7),c(0,5),lty=2)
lines(c(17,17),c(0,5),lty=2)
text(nspp+0.9,res2$mean-0.05,res2$species,col="black",cex=1)




#########################
### MIGDIST by species
#########################

results$fidwin<-1/results$fidwin

par(mfrow=c(3,2),cex.lab=1.8)
tlim<-2000000
spn<-"MD"
plot(results$migdist[results$species==spn & results$migdist<tlim],results$fidsum[results$species==spn & results$migdist<tlim],xlab="Migration distance",ylab="fidsumelity",main=spn,pch=16,col=makeTransparent(rainbow(5)[5],alpha=100),cex=1.3)
points(results$migdist[results$species==spn & results$migdist<tlim],results$fidsum[results$species==spn & results$migdist<tlim],pch=1,cex=1.3,lwd=0.1)
mod<-lm(results$fidsum[results$species==spn & results$migdist<tlim]~results$migdist[results$species==spn & results$migdist<tlim])
#abline(mod)
summary(mod)

spn<-"MS"
plot(results$migdist[results$species==spn],results$fidsum[results$species==spn],xlab="Migration distance",ylab="fidsumelity",main=spn,pch=16,col=makeTransparent(rainbow(5)[4],alpha=100),cex=1.3)
points(results$migdist[results$species==spn],results$fidsum[results$species==spn],pch=1,cex=1.3)
mod<-lm(results$fidsum[results$species==spn]~results$migdist[results$species==spn])
#abline(mod)
summary(mod)
results<-results[results$fidsum<3.7,] ### Major outlier in pronghorn data
spn<-"PH"
plot(results$migdist[results$species==spn],results$fidsum[results$species==spn],xlab="Migration distance",ylab="fidsumelity",xlim=c(0,120000),main=spn,pch=16,col=makeTransparent(rainbow(5)[3],alpha=100),cex=1.3)
points(results$migdist[results$species==spn],results$fidsum[results$species==spn],pch=1,cex=1.3)
mod<-lm(results$fidsum[results$species==spn]~results$migdist[results$species==spn])
#abline(mod)
summary(mod)
spn<-"EK"
plot(results$migdist[results$species==spn],results$fidsum[results$species==spn],xlab="Migration distance",ylab="fidsumelity",main=spn,pch=16,col=makeTransparent(rainbow(5)[1],alpha=100),cex=1.3)
points(results$migdist[results$species==spn],results$fidsum[results$species==spn],pch=1,cex=1.3)
mod<-lm(results$fidsum[results$species==spn]~results$migdist[results$species==spn])
#abline(mod)
summary(mod)
spn<-"BS"
plot(results$migdist[results$species==spn],results$fidsum[results$species==spn],xlab="Migration distance",ylab="fidsumelity",main=spn,pch=16,col=makeTransparent(rainbow(5)[2],alpha=100),cex=1.3)
points(results$migdist[results$species==spn],results$fidsum[results$species==spn],pch=1,cex=1.3)
mod<-lm(results$fidsum[results$species==spn]~results$migdist[results$species==spn])
abline(mod)
summary(mod)

####################################
### PLOT FID OVER TIME BY POPULATION
####################################

## Does not produce nice figure because some fidelity values are extremely large (due to very very high fidelity)
results<-rbind(results[results$species=="MD",],results[results$species=="MS",],results[results$species=="EK",],results[results$species=="PH",],results[results$species=="BS",])
us<-unique(results$species)
par(mfrow=c(1,5),mar=c(3,5,1,1), oma=c(0,0,2,1))
for(j in 1:length(us)){
  u<-unique(results$pop[results$species==us[j]])
  for (i in 1:length(u)){
    title<-substr(u[i],1,5)
    setwd(paste("E:/morrison/SiteFidelityProj/data/",title,sep=""))
    fid<-read.csv(paste(title,"_pop",substr(u[i],7,7),"_fid.csv",sep=""))
    dates<-read.csv(paste(title,"_pop",substr(u[i],7,7),"_dates.csv",sep=""))
    dates<-dates[,-1]
    fid<-fid[,-1]
    fid<-1/fid
    dates<-as.Date(dates)
    comp<-list(dates,fid)
    pop.plot(comp,i,j,title)
  }
}

###############################################
### PLOT by population variation by species ALL
###############################################

mm<-aggregate(popvar$popvar,by=list(popvar$species),FUN=function(x) mean(x,na.rm=TRUE))
se<-aggregate(popvar$popvar,by=list(popvar$species),FUN=function(x) sd(x,na.rm=TRUE)/sqrt(length(x)))

res2<-cbind(mm,se[,2])
colnames(res2)<-c("Species","Mean","SE")
res2<-res2[(order(res2[,2])),]

exp(res2$Mean)

nspp<-1:length(unique(popvar$species))
colr<-c(rainbow(5)[2],rainbow(5)[5],rainbow(5)[4],rainbow(5)[1],rainbow(5)[3])
plot(nspp[1:5],res2[1:5,2],cex=2,xlim=c(1,5.2),ylim=c(0,max(popvar$popvar,na.rm=TRUE)),col=colr,pch=16,lwd=4,xaxt="n",xlab="Species",ylab="Population variability")
arrows(nspp[1:5],res2[1:5,2],nspp[1:5],res2[1:5,2]-res2[1:5,3],lwd=2,angle=0)
arrows(nspp[1:5],res2[1:5,2],nspp[1:5],res2[1:5,2]+res2[1:5,3],lwd=2,angle=0)

points(rep(2,length(popvar$popvar[popvar$species=="MD"]))+runif(length(popvar$popvar[popvar$species=="MD"]),-0.01,0.01),
       popvar$popvar[popvar$species=="MD"],col="gray50")
points(rep(3,length(popvar$popvar[popvar$species=="MS"]))+runif(length(popvar$popvar[popvar$species=="MS"]),-0.01,0.01),
       popvar$popvar[popvar$species=="MS"],col="gray50")
points(rep(4,length(popvar$popvar[popvar$species=="EK"]))+runif(length(popvar$popvar[popvar$species=="EK"]),-0.01,0.01),
       popvar$popvar[popvar$species=="EK"],col="gray50")
points(rep(5,length(popvar$popvar[popvar$species=="PH"]))+runif(length(popvar$popvar[popvar$species=="PH"]),-0.01,0.01),
       popvar$popvar[popvar$species=="PH"],col="gray50")
points(rep(1,length(popvar$popvar[popvar$species=="BS"]))+runif(length(popvar$popvar[popvar$species=="BS"]),-0.01,0.01),
       popvar$popvar[popvar$species=="BS"],col="gray50")
points(nspp[1:5],res2[1:5,2],cex=2,ylim=c(0,1),pch=1,xaxt="n",lwd=1,col="black")
text(nspp[1:5]+0.2,res2[1:5,2],unique(popvar$species)[c(3,1,2,4,5)])

