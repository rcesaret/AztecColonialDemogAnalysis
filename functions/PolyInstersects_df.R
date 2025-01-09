#### PolyInstersects_df.R
#### Rudolf Cesaretti, 5/31/2022

#### "PolyInstersects_df" calculates which polygons overlap with 
#### each other in order to be merged by the "DissolvePoly" 
#### function

#### DEPENDS ON Script1_HelperFunctions.R

pak <- c("igraph", "rgdal", "sp", "sf", "GISTools", "lwgeom", "tidyverse", "tidyr", "data.table", "zoo")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak[2:length(pak)], library, character.only = TRUE))
rm(pak,ip)

###################################################################
######################  PolyInstersects_df  #######################
###################################################################

PolyInstersects_df = function(spdf,     # sp package SPDF class object
                              Expand=F  # Should other sherd dens + demog vars be calculated?
                              ){
  
  spdf1 = as(st_make_valid(st_as_sf(spdf)), "Spatial")
  spdf2 = as(st_make_valid(st_as_sf(spdf)), "Spatial")
  spdf1 <- spChFIDs(spdf1, spdf1$Site)
  spdf2 <- spChFIDs(spdf2, spdf2$Site)
  
  x = gIntersects(spdf1,byid=TRUE,returnDense=F)
  
  x =enframe(x)
  COLZ <- c("V1","V2","V3","V4","V5")
  
  l = as.vector(sapply(x$value, length))
  COLZ <- COLZ[1:max(l)]
  
  nm1 <- c(setdiff(names(x), 'value'), COLZ)
  x=unnest_wider(x, col = "value", names_repair = ~ nm1)
  
  colnames(x)[1] <- "Site"
  x$SiteID = as.numeric(rownames(x))
  for (i in 1:length(COLZ)) {
    
    new <- x[,c(COLZ[i])]
    x[ , ncol(x) + 1] <- new
    
    colnames(x)[ncol(x)] <- paste0(COLZ[i],"ID")
    
  }
  s = paste0(COLZ,"ID")
  for (i in 1:length(COLZ)) {
    
    new <- plyr::mapvalues(as.numeric(unlist(x[,c(s[i])])), 
                           from=x$SiteID, 
                           to=x$Site)
    x[ , ncol(x) + 1] <- new
    
    colnames(x)[ncol(x)] <- paste0(COLZ[i],"Site")
  }
  
  library(igraph)
  g = st_intersects(st_as_sf(spdf1),st_as_sf(spdf1))
  G = graph_from_adj_list(g)
  Gc = components(G)
  df <- data.frame(SiteID = c(1:nrow(x)), PairID = Gc$membership)
  df2 <- data.frame(PairID = c(1:max(Gc$membership)), NumOverlaps = Gc$csize)
  df2$Overlap <- ifelse(df2$NumOverlaps > 1, 1, 0)
  detach("package:igraph", unload = TRUE)
  
  
  df <- left_join(df,df2, by="PairID")
  x <- left_join(x,df, by="SiteID")
  
  if (Expand == T){
    xx = x[,c("Site","PairID")]
    
    x2 = x %>% filter(Overlap==1) %>% filter(!duplicated(paste0(pmax(V1, V2), pmin(V1, V2))))
    spdf1b <- subset(spdf1, Site %in% x2$V1Site)
    spdf2b <- subset(spdf2, Site %in% x2$V2Site)
    
    sf1 = st_as_sf(spdf1b)
    sf2 = st_as_sf(spdf2b)
    
    sf2=sf2 %>% setdiff(sf1)
    
    tmp <- st_intersection(sf1, sf2)
    
    y = colnames(tmp)
    zz = (ncol(tmp) - 1)/2
    y1 = paste0(y[1:zz],".1")
    y2 = paste0(gsub('.{2}$', '', y[(zz+1):(zz*2)]),".2")
    colnames(tmp)<-c(y1,y2,"geometry")
    st_geometry(tmp) <- "geometry"
    colnames(xx)<-c("Site.1","PairID")
    tmp = left_join(tmp, xx, by="Site.1")
    
    x3 = data.frame(PairID = tmp$PairID, Site1 = tmp$Site.1, Site2 = tmp$Site.2)
    
    x3$OverlapArea <- Area_ha_Num(tmp)
    x3$Site1Area = tmp$Area_ha.1
    x3$Site2Area = tmp$Area_ha.2
    x3$Site1_NoOvlp_Area = x3$Site1Area - x3$OverlapArea
    x3$Site2_NoOvlp_Area = x3$Site2Area - x3$OverlapArea
    x3$CombinedArea = x3$Site1_NoOvlp_Area + x3$Site2_NoOvlp_Area + x3$OverlapArea
    
    x3$Site1PopDens = tmp$PopDens.1
    x3$Site2PopDens = tmp$PopDens.2
    x3$OverlapPopDens <- tmp$PopDens.1 + tmp$PopDens.2 #######MEAN????
    x3$CombinedPop <- tmp$O_Population.1 + tmp$O_Population.2
    x3$CombinedAssemb <-  tmp$Tot.Assemb.1 + tmp$Tot.Assemb.2
    x3$CombinedPopDens <- x3$CombinedPop / x3$CombinedArea
    x3$CombinedSherdDens <- x3$CombinedAssemb / (x3$CombinedArea * 10000)
    
    x4 = x3[,c("PairID", "Site1","Site2", "CombinedPopDens", "CombinedSherdDens", "CombinedAssemb")]
    
    x5 = left_join(x, x4, by="PairID")
    #x5 = left_join(x5, x4, by="Site2")
    
    x5$PopDens = coalesce(x5$CombinedPopDens, spdf1@data$PopDens)
    x5$SherdDens = coalesce(x5$CombinedSherdDens, spdf1@data$SherdDens)
    x5$Tot.Assemb = coalesce(x5$CombinedAssemb, spdf1@data$Tot.Assemb)
    
    x = x5
  }
  
  return(x)
}

