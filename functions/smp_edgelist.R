#################################
### function to take 1 random location per day per id per year between time 10am and 2pm
###############################################################

#sampledays is the number of days to sample from. Currenlty, only works with 1

smp_edgelist <- function(data=data, sampledays = 1){
  #some checks 
  if(inherits(data, "data.frame") != TRUE) 
    stop("data is not a dataframe")
  if(all(c("date", "AID") %in% names(data))==FALSE)
    stop("you need at least a date and AID column in your data")
  if(inherits(data$date, "POSIXct") != TRUE) 
    stop("data$date is not POSIXct")
  if(sampledays != 1)
    stop("sampledays must be 1, for now")
  
  tz <- attr(data$date,"tzone")
  year <- as.numeric(strftime(data$date, format = "%Y", tz = tz))
  jul <- as.numeric(strftime(data$date, format = "%j", tz = tz))
  data$id_day_yr <- paste(data$AID, year, jul, sep="_")
  data$rand <- sample(nrow(data),nrow(data))
  data$unique <- 1:nrow(data)
  
  data$time1<-as.numeric(strftime(data$date,format="%H",tz="MST"))  
  
  ### Change time window here
  crittime<-6
  if(nrow(data[data$time1<crittime,])==0) crittime<-8
  if(nrow(data[data$time1<crittime,])==0) crittime<-10
  if(nrow(data[data$time1<crittime,])==0) crittime<-12
  
  data<-data[data$time1<crittime,]   
  ###
  
  data <- data[order(data$rand),]
  data <- data[duplicated(data$id_day_yr)==FALSE,]
  data <- data[order(data$AID, data$date),]
return(data[,names(data) %in% c("unique", "id_day_yr", "rand")==FALSE])
}

