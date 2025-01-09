#### RescaleSpatRast.R
#### Rudolf Cesaretti, 6/10/2022

#### "RescaleSpatRast" 
#### 
#### 
#### 
#### 

pak <- c("rgdal", "raster", "terra")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)

###############################################################
######################  RescaleSpatRast  #######################
###############################################################

RescaleSpatRast <- function(x, minx, maxx){
  minx <- x@ptr$range_min
  maxx <- x@ptr$range_max
  tmp<-(x-minx)/(maxx-minx)
  return(tmp)
}


