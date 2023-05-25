#################################
### function to take 1 random location per day per id per year
###############################################################

#sampledays is the number of days to sample from. Currenlty, only works with 1


smp <- function(data=data){
  #some checks 
  if(inherits(data, "data.frame") != TRUE) 
    stop("data is not a dataframe")
  if(all(c("date", "AID") %in% names(data))==FALSE)
    stop("you need at least a date and AID column in your data")
  if(inherits(data$date, "POSIXct") != TRUE) 
    stop("data$date is not POSIXct")
  tz <- attr(data$date,"tzone")
  
  year <- as.numeric(strftime(data$date, format = "%Y", tz = tz))
  jul <- as.numeric(strftime(data$date, format = "%j", tz = tz))
  data$unique <- 1:nrow(data)
  
  data$id_day_yr <- paste(data$AID, year, jul, sep="_") 
  
  data <- data[duplicated(data$id_day_yr)==FALSE,]
  data <- data[order(data$AID, data$date),]
  data <- data[,names(data) %in% c("unique", "id_day_yr")==FALSE]
  return(data)
}

