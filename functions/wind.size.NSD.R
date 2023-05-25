# function to estimate window sizes for each species using NSD
# July 7 2020

wind.size.NSD <- function(dat2) {
  
  # truncate data to be <= 1year per individual for best fit (Spitz etal)
  
  # Year 1
  xy2 <- dat2 %>% 
    group_by(AID) %>% 
    filter(date <= (min(date) + (365*60*60*24))) %>%
    filter(length(date) >= 300)
  
  if(nrow(xy2)>0) {
    
    xy2$AID2 <- paste0(xy2$AID,"_1")
    
    # Year 2
    xy3 <- dat2 %>% 
      group_by(AID) %>% 
      filter(date <= (min(date) + (365*60*60*24*2)) & date > (min(date) + (365*60*60*24))) %>%
      filter(length(date) >= 300)
    
    # paste together
    if(nrow(xy3)>0) {
      
      xy3$AID2 <- paste0(xy3$AID,"_2")
      
      xy2 <- rbind(xy2,xy3)
      
      # skip individuals that cause convergence problems
      xy2 <- xy2[!xy2$AID2 %in% pp,]
      
      # identify start date - minumum julian date
      w <- xy2 %>% 
        group_by(AID) %>% 
        summarise(yday(min(date)))
      
      if(nrow(w)>0){
        mindt <- substr(as.character((as.Date(paste0("2010-",min(w[,2])),format="%Y-%j"))),6,10)
        
        # convert to ltraj
        elk <- as.ltraj(cbind(as.data.frame(xy2)$x,as.data.frame(xy2)$y),date=xy2$date,id = xy2$AID2)
        
        # fit models
        nsd <- mvmtClass(elk,fam='nsd',stdt = mindt)
        
        ## use the following code to trouble shoot fitting errors
        # for(i in 1:length(elk)){
        # nsd <- mvmtClass(elk[i],fam='nsd',stdt = mindt)
        # }
        # i
        # z <- NULL
        # z <- c(z,unique(xy2$AID2)[i])
        # pp <- c(pp,z)
        
        # identify which individuals are migrants, mixmigrants and dispersers 
        winc <- which(fullmvmt(nsd) & names(topmvmt(nsd)) %in% c("migrant","mixmig","disperser"))
        
        # extract start and end dates of the migration
        if(length(winc)>0){
          
          if(q==1) nsdsum <- NULL
          
          for(i in 1:length(nsd)){
            name.top <- names(topmvmt(nsd[[i]]))
            
            if(name.top %in% c("migrant","mixmig","disperser")){
              
              # extract start end of migration from top model
              x <- mvmt2dt(nsd[[i]], mod=name.top)
              
              if(!is.null(x)){
                # compile dataset  
                nsdsum <- rbind(nsdsum,
                                data.frame(
                                  spp= substr(unique(xy2$AID2)[i],1,2),           
                                  AID = unique(xy2$AID2)[i],
                                  startdate = as.Date(min(xy2$date[xy2$AID2==unique(xy2$AID2)[i]])), #as.POSIXct(tapply(xy2$date,INDEX = as.factor(xy2$AID2),FUN = function(x) min(x,na.rm=T)),origin="1970-01-01"),
                                  enddate = as.Date(max(xy2$date[xy2$AID2==unique(xy2$AID2)[i]])),  #as.POSIXct(tapply(xy2$date,INDEX = as.factor(xy2$AID2),FUN = function(x) max(x,na.rm=T)),origin="1970-01-01"),
                                  autoclassNSD = fullmvmt(nsd)[i],
                                  topNSDmod = name.top,
                                  str1 = yday(x[[1]][1,2]),
                                  end1 = yday(x[[1]][2,2]),
                                  str2 = yday(x[[1]][3,2]),
                                  end2 = yday(x[[1]][4,2])
                                ))
              }
            }
          }# end for loop i
        }
      } 
      
      write.csv(nsdsum, "./NSD/nsdsum.csv",row.names = F)  
      print('STOP: need to manually exclude rows of data with overlapping dates spr and aut migrations')
      
    } # end if-then data available to fit in both years 
  } # end if-then data available to fit in both years 
  
  # # manage
  x$AID.x <- substr(as.character(x$AID),1,(nchar(as.character(x$AID))-2))
  x$per1 <- "spring"
  x$per1[x$str1 >= 182] <- "autumn"
  x$per2 <- "autumn"
  x$per2[x$per1 == "autumn"] <- "spring"
  x$spp <- substr(x$AID,1,2)

  # calculate difference in start/end dates within individuals
  dd <- x %>%
    group_by(AID.x) %>%
    filter(autoclassNSD==TRUE & exclude!='x') %>%
    filter(length(unique(AID))>1) %>%
    summarise(diffsprstr=sd(str1),
              diffsprend=sd(end1),
              diffautstr=sd(str2),
              diffautend=sd(end2))

  dd$spp <- substr(dd$AID.x,1,2)
  aa <- data.frame(diff = c(dd$diffsprstr,dd$diffsprend), spp = c(dd$spp,dd$spp), AID=c(dd$AID.x,dd$AID.x))
  spr <- lmer(diff ~ spp + (1|AID),data=aa)

  aa <- data.frame(diff = c(dd$diffautstr,dd$diffautend), spp = c(dd$spp,dd$spp), AID=c(dd$AID.x,dd$AID.x))
  aut <- lmer(diff ~ spp + (1|AID),data=aa)

  out <- data.frame(spr = c(fixef(spr)[1],fixef(spr)[1]+fixef(spr)[2:6]),
                    aut = c(fixef(aut)[1],fixef(aut)[1]+fixef(aut)[2:6]))
  row.names(out) <- c('CA','EK','MD','MS','WB','ZB')

  
}# end wind.size.NSD function
