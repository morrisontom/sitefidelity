##################################################################
### function to subsample a set of locations per day per animal id
##################################################################
### T. Morrison 9-Nov-18
# note: that the smphours default is once per day (6am). you can add vector of multiple hours per day, 
# but some individuals/days may NOT have points near the sample hours, in which case you may get fewer samples 
# per day than expected

smp <- function(data=data, 
                sampledays = 1, #sampledays: sample frequency (days per sample) 1= sample each day, 2 = sample every other day 
                smphours = 6){  #smphour: vector of hour(s) (i.e. 1-24) for which points should be near within days
  #some checks 
  if(inherits(data, "data.frame") != TRUE) 
    stop("data is not a dataframe")
  if(all(c("date","AID") %in% names(data))==FALSE)
    stop("you need at least a date and AID column in your data")
  if(inherits(data$date, "POSIXct") != TRUE & inherits(data$date, "POSIXxt") != TRUE ) 
    stop("data$date is not POSIXct")
  
  tz <- attr(data$date,"tzone")
  data$hour<- hour(data$date)
  data$year <- year(data$date)
  data$jul <- as.numeric(strftime(data$date, format = "%j", tz = tz))
  
  check.integer <- function(x) {
    x == round(x)
  }
  
  #exclude days that are not divisible by sampledays
  data<-data[check.integer(data$jul/sampledays),]
  
  all<-NULL  
  for(i in 1:length(smphours)){
    data$hour_since_smphour<-(data$hour-smphours[i])
    
    if(i==1){
      #exclude data before smphour if available and select data point nearest to smphour
      new<-aggregate(data$hour_since_smphour,by=list(data$AID,as.Date(data$date)),function(x) {
        out<-min(abs(x))  
        return(out)
      })
    }
    if(i>1){
      
      new<-aggregate(data$hour_since_smphour,by=list(data$AID,as.Date(data$date)),function(x) {
        crit<-smphours[i]-smphours[i-1]
        
        #take hour that is closest to smphour, but exclude hours near and before previous smphour
        y<-x[abs(x)>=crit & (x+smphours[i]) > smphours[i-1]] 
        if(length(y)>0){
          out<-min(y)
        }else{           
          #if no points exist before smphour, take next closest hour after smphour - some points may be duplicated with previous smphours  
          out<-min(abs(x))    
        } 
        return(out)
      }) 
    }
    new[,4]<-paste(new[,3],new[,1],new[,2],sep="")  #hrsincedawn+AID+DATE
    data$tmp<-paste(data$hour_since_smphour,data$AID,as.Date(data$date),sep="")
    all<-rbind(all,data[data$tmp%in%new[,4],])
  }
  
  #exclude duplicated data, as days with fewer datapoints than smpperday will have duplicated data
  all<-all[!duplicated(paste(all$x,all$date)),]
  
  #reorder everything back to original
  all <- all[order(all$AID, all$date),]
  return(all[,names(all) %in% c("tmp", "hour_since_dawn")==FALSE])
}
