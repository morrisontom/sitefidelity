#####################################################
### A function to find when two GPS points from two
### different animal location datasets come within
### a specified distance from each other, given a
### time buffer. 
### Written by: Jerod A. Merkle, November 2015
#####################################################
# 

find.encountr <- function(dat1=spatialpoints(), #spatial points dataframe of your first dataset
                          time1=dat1$POSIX,     #the times (in POSIXct) associated with dat1
                          dat2=spatialpoints(), #spatial points dataframe of your second dataset
                          time2=dat2$POSIX,     #the times (in POSIXct) associated with dat2
                          uniqueid1=1:nrow(dat1),  #unique values associated with dat1, must be basically 1:nrow(dat1)
                          uniqueid2=1:nrow(dat2),  #unique values associated with dat2, must be basically 1:nrow(dat2)
                          time.buff=30,         #buffer in time, in minutes
                          spatial.buff=100)     #buffer in space, in dat1 units
                          { 
                          #cpus=5){              #how many cpus do you want to spread this over
  #prepare packages
  if(all(c("sp","snowfall") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: sp, snowfall")
  #require(sp)
  #require(snowfall)
  
  #some checks...
  if(class(dat1)!="SpatialPointsDataFrame")
    stop("dat1 must be a SpatialPointsDataFrame.")
  if(class(dat2)!="SpatialPointsDataFrame")
    stop("dat2 must be a SpatialPointsDataFrame.")
  if(is.na(is.projected(dat1))==TRUE)
    stop("dat1 must be projected.")
  if(is.na(is.projected(dat2))==TRUE)
    stop("dat2 must be projected.")
  if(proj4string(dat1)!=proj4string(dat2))
    stop("dat1 and dat2 do not have the same projection")
  if(any(class(time1)%in%"POSIXct")==FALSE)
    stop("time1 muset be a POSIXct object.")
  if(any(class(time2)%in%"POSIXct")==FALSE)
    stop("time2 muset be a POSIXct object.")
  if(min(as.numeric(time1))>max(as.numeric(time2)) & 
       max(as.numeric(time1)) < min(as.numeric(time2)))
    stop("Your times do not overlap.")
  if(length(unique(uniqueid1)) < length(uniqueid1))
    stop("Your uniqueid1 is not all unique values")
  if(length(unique(uniqueid2)) < length(uniqueid2))
    stop("Your uniqueid2 is not all unique values")
  if(length(time1) != nrow(dat1))
    stop("Your time1 is not the same length as dat1")
  if(length(time2) != nrow(dat2))
    stop("Your time2 is not the same length as dat2")
  if(anyNA(time1)==TRUE)
    stop("You have NAs in your time1.")
  if(anyNA(time2)==TRUE)
    stop("You have NAs in your time2.")
  tz <- strftime(time1[1], "%Z")
  if(tz %in% strftime(time2, "%Z")==FALSE)
    stop("Your time zones do not match  for time1 and time2.")
  
  #just grab the xy values from dat1 and dat2
  dat1 <- coordinates(dat1)
  dat2 <- coordinates(dat2)
  
  #loop over each point in dat1
#   sfInit(parallel = T, cpus = cpus)   #must change the cpus
#   sfExport("dat1","dat2","time1","time2","time.buff","spatial.buff","uniqueid1","uniqueid2")
#   result <- do.call(rbind, sfClusterApplyLB(uniqueid1, function(i){
  result <- do.call(rbind, lapply(uniqueid1, function(i){
    #grab the xy value from dat1 for this iteration
    xy <- dat1[uniqueid1==i] 
    #calculate distance from xy to all points in dat2
    dist <- sqrt((xy[1]-dat2[,1])^2 + (xy[2]-dat2[,2])^2) 
    #calculate minutes from xy to all points in dat2
    time <- abs(as.numeric(difftime(time1[uniqueid1==i],time2, units="mins"))) 
    
    #NA all dist and time values that are outside of the spatial and time buffer
    dist[dist > spatial.buff] <- NA
    time[time > time.buff] <- NA
    
    #identify the rows in dat2 where there were interactions with i in dat1
    rows <- uniqueid2[is.na(dist)==FALSE & is.na(time)==FALSE]
    
    if(length(rows)==0){
      return(data.frame(uniqueid1=NA, x1=NA, y1=NA, time1=NA,
                        uniqueid2=NA, x2=NA, y2=NA, 
                        time2=NA, dist=NA, timediff=NA))
    }else{
      return(data.frame(uniqueid1=i, x1=xy[1], y1=xy[2], time1=as.character(time1)[uniqueid1==i],
                        uniqueid2=rows, x2=dat2[rows,1], y2=dat2[rows,2], 
                        time2=as.character(time2)[rows], dist=dist[rows], timediff=time[rows]))      
    }
  }))
  #sfStop()
  result <- result[is.na(result$uniqueid1)==FALSE,]
  rownames(result) <- NULL
  result$time1 <- as.POSIXct(result$time1, format="%Y-%m-%d %H:%M:%S", tz=tz)
  result$time2 <- as.POSIXct(result$time2, format="%Y-%m-%d %H:%M:%S", tz=tz)
  if(nrow(result)==0)
    print("There were no interactions!")
  return(result)
} 



