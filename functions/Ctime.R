##################################################################################
# Function to calc 3 indices of temporal predictability from time series data
# Louise Riotte-Lambert (June 2020)
# <louise.riotte-lambert@glasgow.ac.uk>
# adapted by Tom Morrison
# see: Morrison et al. 2021 Drivers of site fidelity in ungulates, J Animal ecology
###################################################################################

Ctime <- function(ind,       # xy locations to be calculated (spatial pts dataframe object) 
                  maps.mask, # raster stack
                  intervals=100) # number of intervals to be used in calc of constancy (n=100)
{
  ## calculate for pixel "ind" the two possible measures of constancy of the NDVI time series
  ## and the measure of contingency on time (periodicity)
  
  serie_temp = extract(maps.mask, ind) 
  
  # if the series is not only NAs
  if(length(which(is.na(serie_temp))) < length(serie_temp))
  {
    ## measure 1: Colwell's constancy measure
    # equal to 1 if NDVI is constant in time
    # very small otherwise
    eff = sapply(1:(length(intervals) - 1)
                 , function(x) length(which((serie_temp > intervals[x]) & (serie_temp <= intervals[x+1]))))
    p = eff/sum(eff)
    Ctime1 = 1 + sum(p[p > 0] * log(p[p > 0]))/log(length(p))
    
    ## measure 2: inverse of standard deviation
    Ctime2 = 1/sd(serie_temp, na.rm=TRUE)
    
    ## periodicity: spectral entropy
    # equal to 1 if the time series is perfectly periodic with no noise
    # is very small if it is not periodic
    # it doesn't work if thre is an NA in the time series. 
    # remove the NAs that are at the beginning or the end. if NAs are in the middle, impossible to calculate
    
    indNotNa = which(!is.na(serie_temp))
    serie_temp2 = serie_temp[1:indNotNa[length(indNotNa)]]
    if(min(indNotNa) > 1){
      serie_temp3 = serie_temp2[-(1:(indNotNa[1]-1))]
    } 
    else {
      serie_temp3 = serie_temp2
    }
    
    if(length(which(is.na(serie_temp3))) == 0 & length(serie_temp3)>1) # change the above if statement from original
    {
      speca = spectrum(serie_temp3, plot=FALSE)
      relative = speca$spec/sum(speca$spec)
      Perio = 1 + sum(relative*log(relative))/log(length(relative))
    }else{Perio = NA}
  }else{
    Ctime1 = NA
    Ctime2 = NA
    Perio = NA
  }
  
  return(c(Ctime1, Ctime2, Perio))
}

