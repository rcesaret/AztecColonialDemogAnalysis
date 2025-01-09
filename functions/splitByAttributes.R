#### splitByAttributes.R
#### Rudolf Cesaretti, 5/31/2022

#### "splitByAttributes" splits sp spatial data objects by an 
#### attribute variable, outputting a list; option to export 
#### as separate .gpkg GIS vector files

###################################################################
#######################  splitByAttributes  #######################
###################################################################

pak <- c("rgdal", "sp", "sf", "GISTools", "lwgeom")#, "tidyverse", "tidyr", "data.table", "zoo"
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)

splitByAttributes = function(spdata,       # sp package SPDF class object
                             attr,         # attribute name to slit by
                             prefix = "",  # prefix for output object names
                             suffix = "",  # suffix for output object names
                             export = F,   # whether to save/export the output
                             dirpath = NA  # directory path to save/export output
                             ){
  
  y <-list() #temp storage list
  
  num = which( colnames(spdata@data)==paste(attr)) 
  
  ATTRS <- as.character(unique(spdata[[num]]))
  
  outnames <- paste0(prefix,ATTRS,suffix)
  
  for (i in 1:length(ATTRS)) {
    tmp <- spdata[spdata[[num]] == ATTRS[i], ]
    y[[i]] <- tmp
    
    if (export == T){
      dirpath = ifelse(is.na(dirpath), getwd(), dirpath)
      writeOGR(tmp, paste0(outnames[i],".gpkg"), paste0(outnames[i]), 
               driver = "GPKG", overwrite_layer=TRUE)
    }
    
  }
  
  names(y) <- outnames
  return(y)
  
}