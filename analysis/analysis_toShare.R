# Analysis of site fidelity data
# 4-Sept-20 
# T. Morrison <thomas.morrison@glasgow.ac.uk>

rm(list=ls(all=TRUE))
require(lme4)
require(lubridate)
require(tidyr)
require(tidyverse)
require(broom)
require(MuMIn)
require(car)

# source Fig3-4 plotting functions
source("figures3-4.R")

# read dataset
dat <- read.csv('./output/output_2020-07-21/fidelity_toshare.csv')
dat$date <- as.POSIXct(dat$date,format="%Y-%m-%d %H:%M:%S")
dat$AID <- as.character(dat$AID)

# change bad values
dat[dat$iyd<0 & !is.na(dat$iyd),] <- NA

###############################################
#### fidelity vs environmental variability
###############################################

# summarize data per individual
dat2 <- dat %>% 
  filter(!is.na(iyd)) %>%
  group_by(AID) %>%
  summarize(iyd_log = mean(iyd_log,na.rm=T),
            spp = spp[1],
            Ctime_entropy = mean(ctime_entropy,na.rm=T),
            Cspace1_entropy = mean(Cspace1_entropy,na.rm=T),
            Ctime_periodicity = mean(ctime_periodicity,na.rm=T),
            pop = studyID[1],
            age = Age[1],
            hr = hr[1],
            hr_log=hr_log[1])

setwd('C:/Users/Administrator/OneDrive - University of Glasgow/Documents/Side projects/Fidelity/MS/J Anim Ecol/Nov20/dataset')
write.csv(dat2,"dataset_sitefidelity_toShare.csv",row.names = F)

# create Table 2
tab2 <- dat2 %>%
  group_by(spp) %>%
  summarise(
    naid=length(unique(AID)),
    nsites=length(unique(pop)),
    q25=summary(hr,na.rm=T)[[2]],
    q75=summary(hr,na.rm=T)[[5]],
    med=summary(hr,na.rm=T)[[3]],
    iqr=(iqr75=summary(hr,na.rm=T)[[5]]-summary(hr,na.rm=T)[[2]]))


# make df
dat2 <- as.data.frame(dat2)

# check for correlation between variables
cor(dat2[,c('Ctime_entropy','Ctime_periodicity','Cspace1_entropy','age')],
    use="pairwise.complete.obs")

# Ctime & Cspace metric somewhat correlated

############ 1. Constancy time ############## 
#subset data
dat3Ct <- dat2[!is.na(dat2$Ctime_entropy),]

# run model with species and var of interest interacting
mCt1 <- lmer(iyd_log ~ Ctime_entropy * spp + hr_log + (1|pop) ,data=dat3Ct)

# examine out
summary(mCt1)
fixParam<-fixef(mCt1)
ranParam<-ranef(mCt1)$pop

# examine significance
coefs <- data.frame(coef(summary(mCt1)))
coefs$p.z <- 1.96 * (1 - pnorm(abs(coefs$t.value)))
coefs

## Predictions
# create newdata.frame
nd <- expand.grid(Ctime_entropy=seq(min(dat3Ct$Ctime_entropy,na.rm=T),
                                        max(dat3Ct$Ctime_entropy,na.rm=T),length.out = 20),
                  hr_log=mean(dat3Ct$hr_log,na.rm=T),
                  pop=unique(dat3Ct$pop))
nd$spp <- substr(nd$pop,1,2)

# remove pops that are not incl in model
nd <- nd[nd$pop %in% row.names(ranef(mCt1)$pop),]

# make predictions
p <- predict(mCt1,newdata=nd,re.form=NULL,type='response')
nd$iyd_log <- p

# remove pop effects
nd2 = nd %>%
  group_by(spp,Ctime_entropy) %>%
  summarise(iyd_log=mean(iyd_log),
            spp=spp[1],
            Ctime_entropy_log=Ctime_entropy[1]
  )

############ 2. Constancy space ############## 
#subset data
dat3Cp <- dat2[!is.na(dat2$Cspace1_entropy),]
dat3Cp <- dat3Cp[!is.na(dat3Cp$iyd_log),]
dat3Cp <- dat3Cp[!is.na(dat3Cp$spp),]
dat3Cp <- dat3Cp[!is.na(dat3Cp$hr_log),]

# run model
mCp1 <- lmer(iyd_log ~ Cspace1_entropy * spp + hr_log + (1|pop) ,data=dat3Cp)

# examine out
summary(mCp1)
fixParam<-fixef(mCp1)
ranParam<-ranef(mCp1)$pop
coefs <- data.frame(coef(summary(mCp1)))
coefs$p.z <- 1.96 * (1 - pnorm(abs(coefs$t.value)))
coefs

# calc conditional and marginal r2
r.squaredGLMM(mCp1)
?r.squaredGLMM

## assumptions
# 1. residuals are not linear or curvilear
plot(dat3Cp$iyd_log,resid(mCp1)) 

# 2. homoscedasticity - unclear how important this test is given lmm absorbs within pop variation
mr <- data.frame(Model.F.Res=residuals(mCp1)) #extracts the residuals and places them in a new column in our original data table
mr$Abs.Model.F.Res <- abs(mr$Model.F.Res) #creates a new column with the absolute value of the residuals
mr$Model.F.Res2 <- mr$Abs.Model.F.Res^2 #squares the absolute values of the residuals to provide the more robust estimate
mr$pop <- dat3Cp$pop
Levene.Model.F <- lm(Model.F.Res2 ~ pop, data=mr) #ANOVA of the squared residuals
anova(Levene.Model.F) #displays the results

# 2b. location scale plot - looking for horizontal line with equal spread of points
plot(mCp1)

# 3. normality of residuals
require(lattice)
qqmath(mCp1, id=0.05)

# predictions
mean.by.spp <- dat3Cp %>%
  group_by(spp) %>%
  summarise(mean(hr_log,na.rm=T))

nd <- expand.grid(Cspace1_entropy=seq(min(dat3Cp$Cspace1_entropy,na.rm=T),
                                      max(dat3Cp$Cspace1_entropy,na.rm=T),length.out = 25),
                  pop=unique(dat3Cp$pop))
                  
nd$spp <- substr(nd$pop,1,2)
nd$hr_log <- mean.by.spp$`mean(hr_log, na.rm = T)`[match(nd$spp,mean.by.spp$spp)]

# remove pops that are not incl in model
nd <- nd[nd$pop %in% row.names(ranef(mCp1)$pop),]

# make predictions
p <- predict(mCp1,newdata=nd,re.form=NULL,type='response')
nd$iyd_log <- p

# summarize predictions
nd2mCp = nd %>%
  group_by(spp,Cspace1_entropy) %>%
  summarise(iyd_log=mean(iyd_log),
            spp=spp[1],
            Cspace1_entropy=Cspace1_entropy[1]
  )

############ 3. Temporal periodicity ############## 
#subset data
dat3Cper <- dat2[!dat2$spp %in% c('CA') & !is.na(dat2$Ctime_periodicity),]
mCper1 <- lmer(iyd_log ~ Ctime_periodicity * spp + hr_log + (1|pop) ,data=dat3Cper)

# examine output
summary(mCper1)
fixParam<-fixef(mCper1)
ranParam<-ranef(mCper1)$pop
coefs <- data.frame(coef(summary(mCper1)))
coefs$p.z <- 1.96 * (1 - pnorm(abs(coefs$t.value)))
coefs

# calc conditional and random r2
r.squaredGLMM(mCper1)

## Maringal Predictions
# create newdata.frame - set wb as the ref pop
nd <- expand.grid(Ctime_periodicity=seq(min(dat3Cper$Ctime_periodicity,na.rm=T),
                                        max(dat3Cper$Ctime_periodicity,na.rm=T),length.out = 25),
                  hr_log=mean(dat3Cper$hr_log,na.rm=T),
                  pop='WB_SE')
nd$spp <- substr(nd$pop,1,2)

# remove pops that are not incl in model
nd <- nd[nd$pop %in% row.names(ranef(mCper1)$pop),]

# make predictions
p <- predict(mCper1,newdata=nd,re.form=NULL,type='response')
nd$iyd_log <- p

# remove pop effects - unneccesary if pop = WB_SE
nd2mCper = nd %>%
  group_by(Ctime_periodicity) %>%
  summarise(iyd_log=mean(iyd_log)
  )


############ 4. WIN-STAY LOSE SWITCH ############## 

# subset data - require more steps
dat3 <- dat[dat$spp %in% c('BS','EK','MD','MS','PH'),]

# subset spring Year 1
y1 <- dat3 %>% 
  group_by(AID) %>% 
  filter(date <= (min(date) + (365*60*60*24))) %>%
  filter(yday(date) >= spr1[1]) %>%
  filter(yday(date) <= spr2[1]) %>% 
  summarise(
    AID=AID[1],
    DFP = mean(DFP,na.rm=T),
    hr = mean(hr,na.rm=T),
    ctime_periodicity = mean(ctime_periodicity,na.rm=T),
    spp=spp[1])

# subset spring Year 2
spr2 <- dat3 %>% 
  group_by(AID) %>% 
  filter(date > (min(date) + (365*60*60*24)) & date <= (min(date) + (365*60*60*24*2))) %>%
  filter(yday(date) >= spr1[1]) %>%
  filter(yday(date) <= spr2[1])  %>%
  summarise(
    AID=AID[1],
    DFP = mean(DFP,na.rm=T),
    ctime_periodicity = mean(ctime_periodicity,na.rm=T),
    iyd_logspr = mean(iyd_log,na.rm=T))

# combine together
y1$iyd_logspr <- spr2$iyd_logspr[match(y1$AID,spr2$AID)]
y1$hr_log <- log(y1$hr) 
y1$pop <- substr(y1$AID,1,5)
y1$DFP_log <- log(y1$DFP)
y1 <- y1[!is.na(y1$DFP) & !is.na(y1$hr) & !is.na(y1$iyd_logspr) & (y1$iyd_logspr)>0, ]

# spring model
m1<-lmer(iyd_logspr ~ DFP_log * spp + hr_log + (1|pop),data=y1)

# examine output
summary(m1)
coefs <- data.frame(coef(summary(m1)))
coefs$p.z <- 1.96 * (1 - pnorm(abs(coefs$t.value)))
coefs

# predictions
nd <- expand.grid(DFP_log=seq(min(y1$DFP_log,na.rm=T),
                              max(y1$DFP_log,na.rm=T),length.out = 25),
                  hr_log=mean(y1$hr_log,na.rm=T),
                  pop=unique(y1$pop))
nd$spp <- substr(nd$pop,1,2)

# remove pops that are not incl in model
nd <- nd[nd$pop %in% row.names(ranef(m1)$pop),]

# make predictions
p <- predict(m1,newdata=nd,re.form=NULL,type='response')
nd$iyd_log <- p

# summarize predictions
nd2dfp = nd %>%
  group_by(spp,DFP_log) %>%
  summarise(iyd_log=mean(iyd_log),
            spp=spp[1],
            DFP_log=DFP_log[1]
  )

############ 5. Age - Experience ############## 

# subset data
dat3 <- dat2[!is.na(dat2$age),]
dat3$spp <- as.character(dat3$spp)
dat3$spp[dat3$spp=='PH'] <- 'aPH'
dat3$spp <- as.factor(dat3$spp)

# create Table 2
x2 <- dat2 %>%
  group_by(spp) %>%
  summarise(
    naid=length(unique(AID)),
    q25=summary(age,na.rm=T)[[2]],
    q75=summary(age,na.rm=T)[[5]],
    med=summary(age,na.rm=T)[[3]],
    iqr=(iqr75=summary(age,na.rm=T)[[5]]-summary(age,na.rm=T)[[2]]))


# model age
mA1 <- lmer(iyd_log ~ age * spp  + hr_log + (1|pop) ,data=dat3)

# output
summary(mA1)
fixParam<-fixef(mA1)
ranParam<-ranef(mA1)$pop
coefs <- data.frame(coef(summary(mA1)))
coefs$p.z <- 1.96 * (1 - pnorm(abs(coefs$t.value)))
coefs

# predictions
nd <- expand.grid(age=seq(min(dat3$age,na.rm=T),
                              max(dat3$age,na.rm=T),length.out = 25),
                  hr_log=mean(dat3$hr_log,na.rm=T),
                  pop=unique(dat3$pop))
nd$spp <- substr(nd$pop,1,2)
nd$spp[nd$spp=='PH'] <- 'aPH'

# remove pops that are not incl in model
nd <- nd[nd$pop %in% row.names(ranef(mA1)$pop),]

# make predictions
p <- predict(mA1,newdata=nd,re.form=NULL,type='response')
nd$iyd_log <- p

# summarize output
ndA = nd %>%
  group_by(spp,age) %>%
  summarise(iyd_log=mean(iyd_log),
            spp=spp[1],
            age=age[1]
  )


##########################
##### SUMMARIZE
##########################

## PLOT ALL 
par(mfcol=c(1,3))
fig4(nd2 = nd2mCp,var = 'Cspace1_entropy',x = dat3Cp)
fig4(nd2 = nd2mCper,var = 'Ctime_periodicity',x = dat3Cper)
fig4(nd2 = nd2dfp,var = 'DFP_log',x = y1)

# create table of coefficients
out_Ct <- tidy(mCt1)
out_Cp <- tidy(mCp1)
out_Cper <- tidy(mCper1)
out_dfp <- tidy(m1)
out_age <- tidy(mA1)

write.csv(out_Ct,paste0(response,'_Ct.csv'),row.names = F)
write.csv(out_Cp,paste0(response,'_Cp.csv'),row.names = F)
write.csv(out_Cper,paste0(response,'_Cper.csv'),row.names = F)
write.csv(out_dfp,paste0(response,'_dfp.csv'),row.names = F)
write.csv(out_age,paste0(response,'_age.csv'),row.names = F)

#######################
#### species model ####
#######################

sm2 <- lmer(iyd_log ~ spp + hr_log + (1|pop),data=dat2)
summary(sm2)

# output from model
summary(sm2)
dat2 <- dat2[!is.na(dat2$AID),]

### predictions for figure
hr_log <- dat2 %>%
  group_by(spp) %>%
  summarise(hr_log=mean(hr_log))

# new data
nd <- data.frame(pop=c('BS_LP','CA_AL','EK_S2','MD_WR',
                       'MS_SR','PH_SB','WB_SE','ZB_SE'), 
                  hr_log=hr_log$hr_log)
nd$spp <- substr(nd$pop,1,2)

# remove pops that are not incl in model
nd <- nd[nd$pop %in% row.names(ranef(sm2)$pop),]

# make predictions
p <- predict(sm2,newdata=nd,re.form=NULL,type='response')
nd$iyd_log <- p

# Use bootstrapp to get CI predictons
mySumm <- function(.) { s <- sigma(.)
c(beta =getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) }
(t0 <- mySumm(sm2)) # just three parameters
bsp <- bootMer(sm2,FUN = mySumm,nsim=1000) 

# 
nd$se <- apply(bsp$t[,1:8],2,sd)
nd$l95 <- nd$iyd_log-(nd$se*1.96)
nd$u95 <- nd$iyd_log+(nd$se*1.96)

# plot species level effects
colr <- c('red','cyan','palegreen3','blue','purple','orange','green','darkred')
colr <- colr[order(nd$iyd_log)]
nd <- nd[order(nd$iyd_log),]

# plot
par(mfcol=c(1,1))
plot(1:nrow(nd),nd$iyd_log,col=colr,xaxt='n',cex=1.9,pch=16,
     ylim=c(min(dat2$iyd_log,na.rm=T),max(dat2$iyd_log,na.rm=T)),
     # xlim=c(0,9),
     xlab='Species',ylab='Min inter-year distance (log-meters)')

# add raw points
for(i in 1:nrow(nd)){
  ff <- dat2$iyd_log[dat2$spp==nd$spp[i]]
  points(jitter(rep(i,length(ff)),amount = 0.10),ff,
         col=makeTransparent('grey70'),pch=1,cex=0.6)
}

# add CI
arrows(x0 = 1:nrow(nd),y0 = nd$iyd_log,x1 = 1:nrow(nd),y1 = nd$l95,angle=90,length = 0.05,lwd=2.6)
arrows(x0 = 1:nrow(nd),y0 = nd$iyd_log,x1 = 1:nrow(nd),y1 = nd$u95,angle=90,length = 0.05,lwd=2.6)
points(1:nrow(nd),nd$iyd_log,xaxt='n',col=(colr),cex=1.9,pch=16)
points(1:nrow(nd),nd$iyd_log,xaxt='n',cex=1.9,lwd=2.6,pch=1)

# test equal means across populations within species - 
dd <- dat2[dat2$spp=='BS',]
Anova(lm(iyd_log ~ hr_log + pop,data=dd),type=2)
dd <- dat2[dat2$spp=='EK',]
Anova(lm(iyd_log ~ hr_log + pop,data=dd),type=2)
dd <- dat2[dat2$spp=='MD',]
Anova(lm(iyd_log ~ hr_log + pop,data=dd),type=2)
dd <- dat2[dat2$spp=='MS',]
Anova(lm(iyd_log ~ hr_log + pop,data=dd),type=2)
dd <- dat2[dat2$spp=='PH',]
Anova(lm(iyd_log ~ hr_log + pop,data=dd),type=2)

#######################
#### SEASONAL model ###
#######################

# add week of year to full 
dat$weeks <- week(results$date)

# fit weekly model to iyd_log
for(j in 1:length(uSPP)){
  
  obs <- dat[dat$spp==uSPP[j],]
  obs <- obs[obs$weeks!=53,]
  obs$iyd_log <- log(obs$iyd)
  
  if(!uSPP[j] %in% c("CA","WB","ZB")){
    obsmod<-lmer(iyd_log ~ hr_log + as.factor(weeks) + (1|title/AID),data=obs)
    
  }else{
    obsmod<-lmer(iyd_log ~ hr_log + as.factor(weeks) + (1|AID),data=obs)
  }
  
  # predicted effects using Bootstrapped CI's
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

# Plot Seasonal data 3
setwd("C:/users/administrator/OneDrive - University of Glasgow/Documents/Side projects/Fidelity/MS/J Anim Ecol/Submit_Jun20//")
fig3(fn,dtrut)
