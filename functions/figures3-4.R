## Plots 3 & 4 site fidelity ms
## updated Aug 2020

# Plot transparency function
makeTransparent<-function(someColor, alpha=130)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
} #ENd function

#####################
# PLOT FIGURE 3
fig3 <- function(fn,dtrut){
  reord <- c(4,5,1,3,6,7,8,2)
  colr <- c('red','cyan','palegreen3','blue','purple','orange','green','darkred')
  splo <- c("Bighorn sheep",'Caribou','Elk','Mule deer','Moose','Pronghorn','Wildebeest','Zebra')
  uSPP<-unique(fn$spp)
  colr <- colr[reord]
  splo <- splo[reord]
  uSPP <- uSPP[reord]
  
  #plotting parms
  alp <- 40
  par(mfrow=c(2,4),mar=c(0,3,2,1),oma=c(3,4,2,4),font.axis=1.2, cex.axis=1.2)
  
  for(i in 1:length(uSPP)){
    x <- fn[fn$spp==uSPP[i],]
    seas<-dtrut[as.character(dtrut$spp)==uSPP[i],]
    mono <- as.Date(paste(2015,1:52,1,sep="-"), format="%Y-%U-%u")
    
    # set up parameters
    part <- c(seas$parturition-14,seas$parturition-14,seas$parturition+14,seas$parturition+14)
    part <- as.Date(paste(2015,part,1,sep="-"),format="%Y-%j")
    rut <- c(seas$rut-14,seas$rut-14,seas$rut+14,seas$rut+14)
    rut <- as.Date(paste(2015,rut,1,sep="-"),format="%Y-%j")
    spr <- c(rep(as.Date(paste(2015,seas$spr1,sep="-"),format="%Y-%j"),2),
             rep(as.Date(paste(2015,seas$spr2,sep="-"),format="%Y-%j"),2))
    dip <- c(as.Date("2014-12-20"),as.Date("2016-01-12"))
    xlim <- as.Date(c("2015-01-01","2015-12-31"))
    xaxis <- 'n'
    
    # species specific parms
    if(uSPP[i] %in% c('BS','MD','MS','EK')) ylim<- c(3.5,8)
    if(uSPP[i] %in% c('PH')) ylim<- c(6,9)
    if(uSPP[i] == c('PH')) xaxis <- NULL
    if(uSPP[i] %in% c('WB')) {
      ylim<- c(7,10) 
      spr <- c(rep(as.Date(paste(2015,305,sep="-"),format="%Y-%j"),2),
               rep(as.Date(paste(2015,334,sep="-"),format="%Y-%j"),2))
      spr2 <- c(rep(as.Date(paste(2015,60,sep="-"),format="%Y-%j"),2),
               rep(as.Date(paste(2015,105,sep="-"),format="%Y-%j"),2))
      xaxis <- NULL
    }
    if(uSPP[i] %in% c('ZB')) {
      ylim<- c(7,10)
      part <- c(NA,NA,NA,NA)
      rut <- c(NA,NA,NA,NA)
      xaxis <- NULL
      spr <- c(rep(as.Date(paste(2015,305,sep="-"),format="%Y-%j"),2),
               rep(as.Date(paste(2015,334,sep="-"),format="%Y-%j"),2))
      spr2 <- c(rep(as.Date(paste(2015,60,sep="-"),format="%Y-%j"),2),
                rep(as.Date(paste(2015,105,sep="-"),format="%Y-%j"),2))
    }
    if(uSPP[i] %in% c('CA')) {
      ylim<- c(9,14)
      xaxis <- NULL
    }
    # plot
    plot(mono,x$mean,type="l",
         lwd=2,xlab="",
         ylab="",xaxt=xaxis,
         yaxt="n",col="blue",
         ylim=ylim,xlim=xlim,
         cex.axis=1.2)
    lines(mono,rep(mean(x$mean),length(mono)),
          type="l",lty=2,lwd=2,
          col="grey")
    polygon(c(spr),c(ylim[1],ylim[2],ylim[2],ylim[1]),col=makeTransparent("green",alpha = alp),border=NA)
    polygon(c(part),c(ylim[1],ylim[2],ylim[2],ylim[1]),col=makeTransparent("yellow",alpha = alp),border=NA)
    polygon(c(rut),c(ylim[1],ylim[2],ylim[2],ylim[1]),col=makeTransparent("purple",alpha = alp),border=NA)
    if(uSPP[i] %in% c('WB','ZB')){ 
      polygon(c(spr2),c(ylim[1],ylim[2],ylim[2],ylim[1]),col=makeTransparent("green",alpha = alp),border=NA)
    }
    polygon(c(mono,rev(mono)),c(x$l95,rev(x$u95)),col=makeTransparent("blue",alpha = alp),border=NA)
    lines(mono,x$mean,lwd=2,col="blue")
    # lines(mono,x$l95,lwd=1,lty=2,col="blue")  
    # lines(mono,x$u95,lwd=1,lty=2,col="blue")  
    axis(2,c(seq(0,ylim[2],1)),col.axis="black",cex.axis=1.2)
    mtext(splo[i],side = 3)
  }
}


###########################################################
###########################################################

# PLOT FIGURE 4
fig4 <- function(nd2,var,x){
  
  spp <- c("BS","CA","EK","MD","MS","PH","WB","ZB")
  
  # plot 1
  if(var=='Cspace1_entropy'){
    uspp <- unique(nd2$spp)
    ltyr <- c(1,1,1,1,1,3,3,3)
    lwdr <- c(3,3,3,3,3,2,2,2)
    ci95 <- 4.26538652*1.96
    w=data.frame(x=c(0.0284,0.419),y=c(7.91,6.18))
  }
  
  # plot 2
  if(var=="Ctime_periodicity"){
    uspp<-c("BS","EK","MD","MS","PH","WB","ZB")
    nd2 <- as.data.frame(nd2)
    ss <- nd2
    ci95 <- 0.66914111*1.96 # est bootMer() nsim =1000
  }
  
  # plot 3
  if(var %in% c('DFP_log','DFP')) {
    # uspp <- 'EK'
    uspp <- unique(nd2$spp)
    x <- as.data.frame(x)
    x$iyd_log <- x$iyd_logspr
    x <- x[!is.na(x$hr),]
    ltyr <- c(3,2,3,3,3)
    lwdr <- c(2,3,2,2,2)
    ci95 <- 0.05033335*1.96 # bootstrapped SE using bootMer(nsim=1000) 
  }
  
  # plot 4 - age
  if(var=='age'){
    uspp <- unique(nd2$spp)
    ltyr <- c(1,3,1,1,1,3,3,3)
    lwdr <- c(3,2,3,3,3,2,2,2)
    ci95 <- 0.97273723*1.96
  }
  
  # colors
  colr <- c('red','cyan','palegreen3','blue','purple','orange','green','darkred')
  splo <- c("bighorn sheep",'caribou','elk','mule deer','moose','pronghorn','wildebeest','zebra')
  colr <- colr[which(spp %in% uspp)]
  splo <- splo[which(spp %in% uspp)]
  
  if(var!='Ctime_periodicity'){
    for(i in 1:length(uspp)){    
      ss <- nd2[nd2$spp==uspp[i],]
      ss <- as.data.frame(ss)
      ss <- ss[which(ss[,var]>=min(x[x$spp==uspp[i],var],na.rm=T) &
                       ss[,var]<=max(x[x$spp==uspp[i],var],na.rm=T)) ,]
      
      if(i==1) {
        plot(ss[,var],ss$iyd_log,
             xlim=c(min(x[,var],na.rm=T),
                    max(x[,var],na.rm=T)),
             ylim=c(3.5,12),
             # ylim=c(0,20),
             type='l',col=colr[i],xlab=var,
             ylab="Inter-year distance log-meters",
             lty=ltyr[i],
             lwd=lwdr[i])
        
        if(var=='Cspace1_entropy'){
          
          # lines(w$x,w$y,lwd=3)
          
          # add CI
          # polygon(x=c(min(w$x),max(w$x),
          #             max(w$x),min(w$x)),
          #         y=c(w$y[1]-ci95,w$y[nrow(w)]-ci95,
          #             w$y[nrow(w)]+ci95,w$y[1]+ci95),
          #         col=makeTransparent('grey60',alpha=20),
          #         border=makeTransparent('grey60',alpha=50))
          
          
        }
        
      }
      if(nrow(ss)>0){
        points(x[,var][x$spp==uspp[i]],
               x$iyd_log[x$spp==uspp[i]],col=makeTransparent(colr[i],alpha = 80),
               pch=16,cex=0.8)
      }
      
      # add species specific lines
      lines(ss[,var],ss$iyd_log,col=colr[i],
            lwd=lwdr[i],lty=ltyr[i])
      text(mean(ss[,var]),mean(ss$iyd_log),labels=splo[i])
      
      # if(var=='DFP' & uspp[i]=='EK'){
      #   polygon(x=c(min(ss$DFP),max(ss$DFP),max(ss$DFP),min(ss$DFP)),
      #           y=c(ss$iyd_log[1]-ci95,ss$iyd_log[nrow(ss)]-ci95,
      #               ss$iyd_log[nrow(ss)]+ci95,ss$iyd_log[1]+ci95),
      #           col=makeTransparent(colr[i],alpha=20),border=makeTransparent(colr[i],alpha=50))  
      # }
    }
  }
  
  if(var=='Ctime_periodicity'){ 
    for(i in 1:length(uspp)){   
      if(i==1)
        plot(nd2[,var],nd2$iyd_log,lwd=3,
             xlab=var,ylab="Inter-year distance log-meters",
             type='l',
             xlim=c(min(x[,var],na.rm=T),
                    max(x[,var],na.rm=T)),
             ylim=c(3.5,12))
      # add points
      points(x[,var][x$spp==uspp[i]],
             x$iyd_log[x$spp==uspp[i]],col=makeTransparent(colr[i],alpha = 80),
             pch=16,cex=0.8)
    }
    #plot again on top of points
    lines(nd2[,var],nd2$iyd_log,lwd=3)
    
    # add CI
    # polygon(x=c(min(ss$Ctime_periodicity),max(ss$Ctime_periodicity),
    #             max(ss$Ctime_periodicity),min(ss$Ctime_periodicity)),
    #         y=c(ss$iyd_log[1]-ci95,ss$iyd_log[nrow(ss)]-ci95,
    #             ss$iyd_log[nrow(ss)]+ci95,ss$iyd_log[1]+ci95),
    #         col=makeTransparent('grey60',alpha=20),border=makeTransparent('grey60',alpha=50))
    # 
  }
}
