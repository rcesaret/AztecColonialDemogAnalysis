#### Script1_HelperFunctions.R
#### Rudolf Cesaretti, 5/31/2022

pak <- c("rgdal", "sp", "sf", "GISTools", "lwgeom", "tidyverse", "tidyr", "data.table", "zoo")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)


###############################################################
####################### Helper Functions ######################
###############################################################

pasteuniquesep <- function(x) {
  z = paste(unique(x),sep="; ")
  return(z)
}

pasteunique <- function(x) {
  z = paste(unique(x),collapse="; ")
  return(z)
}

NAPH <- function(x) {
  #x <- NA
  y <- pmax(x)
  y <- ifelse(x == x, NA, NA)
  return(y)
}

MeanNA <- function(x) {
  y <- mean(x, na.rm=T)
  return(y)
}

MaxNA <- function(x) {
  y <- max(x, na.rm=T)
  return(y)
}

MinNA <- function(x) {
  y <- min(x, na.rm=T)
  return(y)
}

MedNA <- function(x) {
  y <- median(x, na.rm=T)
  return(y)
}

SumNA <- function(x) {
  y <- sum(x, na.rm=T)
  return(y)
}

Area_ha_Num <- function(spdf) {
  vec <- as.numeric(st_area(st_as_sf(spdf))*0.0001)
  return(vec)
}

Perim_m2_Num <- function(spdf) {
  vec <- as.numeric(st_perimeter(st_as_sf(spdf)))
  return(vec)
}

centroid_coords <- function(spdf) {
  x <- coordinates(gCentroid(spdf,byid=TRUE))
  df <- data.frame(x=c(x[,1]),y=c(x[,2]))
  return(df)
}