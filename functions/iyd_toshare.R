### Calculate min IYD's by individual
# Used in Morrison et al. 2021 "Drivers of Site Fidelity in ungulates" J Anim Ecology
# T Morrison thomas.morrison@glasgow.ac.uk, 28-Dec-2020

iyd <- function(x,    # df of GPS data with colnames: AID='animalID',date in POSIX, x,y in UTM)
              wind,   # vector of possible window sizes (in days)
              ) {  
  
  # convert data to df if a spatial object
  x <- as.data.frame(x)
  
  # identify unique individuals
  u <- unique(x$AID)
  
  # loop through each unique individual and return a data frame allanim
  allanim <- NULL
  for(i in 1:length(u)){
    
    # subset to individual i
    a1 <- x[x$AID==u[i],]
    
    # create new columns for each possible window size
    windownames <- paste0('iyd',c('REAL',wind[2:length(wind)]))
    a1[,windownames] <- NA
    
    # identify max/min dates of that individual
    date2 <- as.Date(a1$date)
    startdt <- min(date2)
    enddt <- max(date2)
    
    # check that trajectory spans at least 1 year + window size
    if(as.numeric(enddt-startdt) > (365+min(wind))){
      
      # loop through each window size p
      for(p in 1:length(wind)){
        
        # start date for calculating iyd
        startfid <- startdt + wind[p] + 365
        
        # all dates for which iyd is calculated
        fidcalcdt <- which(date2 >= startfid)
        
        # loop through each date j and calculate min distance
        for(j in fidcalcdt){
          dt1 <- date2[j]
          focal <- a1[j,]
          
          # points around focal date j in year 1
          y1 <- a1[date2 >= (dt1-365-(wind[p]/2)) & 
                     date2 <= (dt1-365+(wind[p]/2)),]
          
          # points around focal date j in year 2
          y2 <- a1[date2>=(dt1-(wind[p]/2)) & 
                     date2<=(dt1+(wind[p]/2)),]
  
          # Only calc IYD for individuals with at least 80% of points in y1 - 
          # Needed if there are missing dates
          if(nrow(y1) > (wind[p]*0.8)){
            
            # calculate min distance using utm - thanks Pythagorean
            a1[j,windownames[p]] <- min(sqrt((focal$y-y1$y)^2+(focal$x-y1$x)^2))
            
          }else{
            a1[j,windownames[p]] <- NA
          }# end ifthen 80%
        }# end loop date j
      }# end loop window p
    }# end loop indiv i
    
    # collate new columns of data, which can be appended to original dataset
    allanim <- rbind(allanim,a1[,c(windownames)])
  }
  return(allanim)
}        
