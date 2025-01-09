#### SpInt_Search.R
#### Rudolf Cesaretti, 6/29/2022

#### "SpInt_Search" 
#### 
#### 
#### 
#### 

pak <- c("rgdal", "sp", "sf", "GISTools", "lwgeom", "tidyverse", "tidyr", "data.table", "zoo")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)


library(raster)
library(terra)
library(scales)
library(pracma)
library(modelsummary)
library(cowplot)


###############################################################
######################  SpInt_Search  #######################
###############################################################