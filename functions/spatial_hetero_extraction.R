################# code to calculate spatial heterogeneity of a polygon
##### by Jerod Merkle, 3 June 2015 ###############


#### SOURCE THE FOLLOWING CODE FIRST
spathetero <- function(spdf = spdf, #this is your spatial polygons dataframe with an id and area column
                       buff=30000, #size of buffer around polygon, in meters
                       NDVIfolder="E:/StatewideNDVI/ndvi_bisc_calcs", #this is the folder where the NDVI calculation data are stored
                       year=1995,  #the year in question. Must change this to correspond with the year of polygons in spdf 
                       name=title, #name of the population/species in question
                       percentOFcell=0.5) # when cells are not completely covered by the polygon, what percent coverage is the cutt-of to keep the cell in the analysis?
{
  
  #require(sgeostat)
  
  #some checks
  if(all(c("id")%in% names(spdf))==FALSE)
    stop("your spdf must contain an id column.")
  drs <- dir(NDVIfolder)
  if(all(c("maxIRGdate.grd", paste("IRGVals_", year, ".grd", sep="")) %in% drs)==FALSE)
    stop("you do not have the following files in your NDVIfolder: maxIRGdate.grd and IRGVals_year.grd")
  if(class(spdf)!="SpatialPolygonsDataFrame")
    stop("spdf must be a Spatial polygons DataFrame")
  if(year==1995) stop("I'm not convinced you specified the right year for the spdf object")
  
  #add buffer to polygon
  spdf$area2 <- gArea(spdf, byid=TRUE)/1000000
  spdf <- gBuffer(spdf, byid=TRUE, width=buff)
  spdf$areaWbuff <- gArea(spdf, byid=TRUE)/1000000
  #swith to character so can deal with multiple types of id formats
  spdf$id <- as.character(spdf$id)
  
  #first, identify mean maxIRG date for each year in each polygon
  yr <- paste("X", year, sep="")
  yr2 <- year
  maxIRGdate <- stack(paste(NDVIfolder, "/maxIRGdate.grd", sep=""))
  maxIRGdate <- maxIRGdate[[yr]]
  
  #specify u (each polygon's id) to loop over
  u <- spdf$id
  params <- do.call(rbind, lapply(1:length(u), function(e){
    print(e)
    temp <- extract(maxIRGdate, spTransform(spdf[spdf$id==u[e],], CRS(projection(maxIRGdate))),
                    weights=T, df=T, cellnumbers=TRUE, normalizeWeights=FALSE)

    
    if(nrow(temp)==0){
      return(data.frame(shapefileID=u[e], year= yr2, area=spdf$area[e], name=name,
                        range=NA, sill=NA, nugget=NA, comments="Polygon fell outside raster"))
    }else{ 
      
      temp$x <- xFromCell(maxIRGdate,temp$cell)
      temp$y <- yFromCell(maxIRGdate,temp$cell)
      
      temp$x <- temp$x/1000
      temp$x <- temp$x-min(temp$x, na.rm=T)
      temp$y <- temp$y/1000
      temp$y <- temp$y-min(temp$y, na.rm=T)
      
#         vls <- temp[,yr]
#         vls <- (vls-min(vls, na.rm=TRUE))/(max(vls, na.rm=T)-min(vls, na.rm=TRUE))
#         plot(temp$x[is.na(temp[,yr])==FALSE], temp$y[is.na(temp[,yr])==FALSE], 
#               col=gray(level=vls[is.na(temp[,yr])==FALSE]))
#         rm(vls)
      #  this removes cells that are covered by less than 50% by the polygon
      temp <- temp[temp$weight >=percentOFcell,]
      
      #like in Teitelbaum et al. 2015 (Ecol Lett), randomly reduce the points to reduce computation time
     #temp <- temp[sample(1:nrow(temp),nrow(temp)/64),]  #subset to about 2km resolution
      if(nrow(temp) > 6000){
        temp <- temp[sample(1:nrow(temp),6000),]
      }
      
      temp[,yr] <- temp[,yr]-min(temp[,yr],na.rm=T) #scale values so that the minimum is 0
      #temp[,yr] <- (temp[,yr]-min(temp[,yr], na.rm=TRUE))/(max(temp[,yr], na.rm=T)-min(temp[,yr], na.rm=TRUE))
      #temp[,names(r)] <- temp[,names(r)]*100
      
      #develop semivariogram
      gd <- as.geodata(obj=temp, coords.col=(1:ncol(temp))[names(temp)%in%c("x","y")], 
                       data.col=(1:ncol(temp))[names(temp)==yr])
      v <- variog(gd, uvec=20)
      #plot(v)
      #fit a model to the variogram
      m <- variofit(v, ini.cov.pars=c(500,5), cov.model="exponential", messages=FALSE, fix.nugget=TRUE,
                    control=list(maxit=10000), limits=pars.limits(sigmasq=c(0,1)))
      #lines(m) 
      #round(m$cov.pars,4)
      #m$message
      
      return(data.frame(shapefileID=u[e], year= yr2, area_orig=spdf$area2[e], areaWbuff=spdf$areaWbuff[e], 
                        name=name, range=m$cov.pars[2], comments=m$message))
      
    }
  }))
  return(params)
}