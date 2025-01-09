#### CatchmentPolys.R
#### Rudolf Cesaretti, 6/10/2022

#### "CatchmentPolys" 
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
######################  CatchmentPolys  #######################
###############################################################
#dissolve polys first!
#do parallel!!
#https://stackoverflow.com/questions/37738251/how-to-expand-the-polygon-to-reach-nearby-line-in-r

CatchmentPolys = function(sitepolys, 
                          border, 
                          dist_lim, 
                          method = c("centroids","borders"), 
                          plot.results = F){
  
  poly <- sitepolys
  poly = as(st_make_valid(st_as_sf(poly)), "Spatial")
  
  
  if (method == "centroids"){
    v =voronoi(poly, ext=extent(border))
    c <- crop(v, border)
    
    c=c[order(c$AggSite),]
    poly=poly[order(poly$AggSite),]
    
    plist <- splitByAttributes(spdata = poly, attr = "AggSite", suffix="_SitePoly") 
    clist <- splitByAttributes(spdata = c, attr = "AggSite", suffix="_VorPoly")
    
    vp2out <- list()
    for (j in 1:length(plist)){
      sp <- plist[[j]]
      sp = as(st_make_valid(st_as_sf(sp)), "Spatial")
      vp <- clist[[j]]
      vp = as(st_make_valid(st_as_sf(vp)), "Spatial")
      vpid <- sapply(slot(vp, "polygons"), function(x) slot(x, "ID"))
      b <- gBuffer(sp, width=dist_lim)
      vp2 <- raster::intersect(vp, b)
      vp2out[[j]] <- vp2
    }
    vp3 <- do.call(bind, vp2out)
    vp3 = as(st_make_valid(st_as_sf(vp3)), "Spatial")
    #vp3@data <- left_join(vp3@data,poly@data,by="AggSite")
    vp3@data$Catchment_ha <- area(vp3)/10000
    vp3@data$CatchmentBeyond_ha <- (area(vp3)-area(poly))/10000
  }
  
  if (method == "borders"){
    # get the coordinates of the polygons    
    p <- unique(geom(poly))
    v <- voronoi(p[, c('x', 'y')], ext=extent(border))
    v = as(st_make_valid(st_as_sf(v)), "Spatial")
    
    # assign group id to the new polygons
    v$group <- p[v$id, 1]
    
    # aggregate (dissolve) polygons by group id
    a <- aggregate(v, 'group')
    # remove areas outside of the zone
    c <- crop(a, border)
    c = as(st_make_valid(st_as_sf(c)), "Spatial")
    # add another identifier 
    c$AggSite <- poly$AggSite[c$group]
    
    c=c[order(c$AggSite),]
    poly=poly[order(poly$AggSite),]
    
    
    plist <- splitByAttributes(spdata = poly, attr = "AggSite", suffix="_SitePoly") 
    clist <- splitByAttributes(spdata = c, attr = "AggSite", suffix="_VorPoly")
    vp2out <- list()
    
    for (j in 1:length(plist)){
      sp <- plist[[j]]
      sp = as(st_make_valid(st_as_sf(sp)), "Spatial")
      vp <- clist[[j]]
      vp = as(st_make_valid(st_as_sf(vp)), "Spatial")
      vpid <- sapply(slot(vp, "polygons"), function(x) slot(x, "ID"))
      b <- gBuffer(sp, width=dist_lim)
      vp2 <- raster::intersect(vp, b)
      if (length(vp2@polygons[[1]]@Polygons) > 1){##
        #r=rast(ncols=1000, nrows=1000,crs=vp@proj4string@projargs,extent=extent(vp))
        #vp2=vect(vp, extent=extent(vp), crs=vp@proj4string@projargs)
        r=rast(ncols=750, nrows=750,crs=vp2@proj4string@projargs,extent=extent(vp2))##
        vp2=vect(vp, extent=extent(vp2), crs=vp2@proj4string@projargs)##
        rv <- terra::rasterize(vp2, r, 'group')
        vr = terra::as.polygons(rv)
        vp2 <- as(vr, "Spatial")
        vp2@data$AggSite <- vp@data$AggSite##
        vp2 <- spChFIDs(vp2, vpid)
        vp2 <- SpatialPolygonsDataFrame(vp2,vp@data)
        #vp2=sp::spTransform(vp2, CRSobj = vp@proj4string)
        vp2 = as(st_make_valid(st_as_sf(vp2)), "Spatial")
        #vp <- vp2
      }
      #b <- gBuffer(sp, width=catch_dist_lim)
      #vp2 <- raster::intersect(vp, b)
      vp2out[[j]] <- vp2
    }
    vp3 <- do.call(bind, vp2out)
    vp3 = as(st_make_valid(st_as_sf(vp3)), "Spatial")
    vp3@data <- left_join(vp3@data,poly@data,by="AggSite") %>% dplyr::select(-group)
    vp3@data$Catchment_ha <- area(vp3)/10000
    vp3@data$CatchmentBeyond_ha <- (area(vp3)-area(poly))/10000
  }
  
  if (plot.results == T){
    #plot(v, lwd=1)
    plot(c, col="green", lwd=2)
    plot(vp3, add=TRUE, lwd=2, col="blue")
    plot(poly, add=TRUE, col="red")
    #text(poly, "AggSite", cex=0.5)
  }
  return(vp3)
}
