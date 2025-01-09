#### DissolvePoly.R
#### Rudolf Cesaretti, 5/31/2022

#### "DissolvePoly" merges the polygons and attribute variables  
#### of any group of polygons according to a specified ID 
#### variable input 

#### DEPENDS ON Script1_HelperFunctions.R

pak <- c("rgdal", "sp", "sf", "GISTools", "lwgeom", "tidyverse", "tidyr", "data.table", "zoo")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)


##################################################################
########################  DissolvePoly  ##########################
##################################################################

DissolvePoly = function(spdf,           # sp package SPDF class object
                        id,             # ID variable by which to dissolve polys
                        varz,           # variables passed to pmap
                        funz,           # functions passed to pmap
                        orderz,         # ordering of output variables
                        remove.id = T,  # whether to remove id variable for outputs
                        rename = T,     # whether to rename the sites == paste0(SurvReg,"-",Period,"-",Number)
                        demog = F,      # whether to recalculate demographic variables
                        cont = F        # whether to recalculate Fw & Bw Continuity variables
                        ){
  #components = F, 
  idz = spdf@data %>% pull(!!sym(id)) 
  
  Dpoly <- unionSpatialPolygons(spdf, IDs = idz)
  
  df = spdf@data
  colnames(df) <- colnames(spdf@data)
  
  df2 <- list(.vars = varz, .funs = funz) %>% 
    pmap(~ df %>% group_by( !!sym(id) ) %>% summarise_at(.x, .y)) %>% 
    reduce(inner_join, by = id) %>%
    select( !!!syms(orderz) )
  
  if (rename == T){
    df2 <- df2 %>% group_by(SurvReg) %>% 
      mutate(
        Number = seq_along(SurvReg)) %>% 
      ungroup() %>% 
      mutate(
        Site = paste0(SurvReg,"-",Period,"-",Number))
  }
  
  df2 <- as.data.frame(df2)
  rownames(df2) = df2[,c(id)]
  
  outpoly <- SpatialPolygonsDataFrame(Dpoly, df2)
  coordzz <- centroid_coords(outpoly)
  outpoly@data$East  <- coordzz$x
  outpoly@data$North <- coordzz$y
  outpoly@data$Area_ha <- Area_ha_Num(outpoly)
  outpoly@data$Perim_m2 <- Perim_m2_Num(outpoly)
  
  if (demog == T){
    outpoly@data <- outpoly@data %>% rowwise() %>% mutate(
      UrbanPop = ifelse(Population > UrbanThresh, Population - UrbanThresh, 0),
      RuralPop = ifelse(Population < UrbanThresh, Population, UrbanThresh)) %>%
      ungroup() %>% mutate(
        PctUrban = UrbanPop / Population,
        PctRural = RuralPop / Population,
        PopDens = Population / Area_ha,
        SherdDens = Tot.Assemb / (Area_ha * 10000),
        UrbanScale = Population / UrbanThresh)
  }
  
  if (cont == T){
    outpoly@data <- outpoly@data %>% rowwise() %>% mutate(
      AreaBwCont = ifelse(PeriodType == "Overlap", 1, 
                          ifelse(Period == "EF", NA, BwOvlp.Area / Area_ha)), 
      PopBwCont = ifelse(PeriodType == "Overlap", 1, 
                         ifelse(Period == "LA", NA, BwOvlp.Pop / Population)),
      PopFwCont = ifelse(PeriodType == "Overlap", 1, 
                         ifelse(Period == "EF", NA, FwOvlp.Pop / Population)), 
      AreaFwCont = ifelse(PeriodType == "Overlap", 1, 
                          ifelse(Period == "LA", NA, FwOvlp.Area / Area_ha))) %>% 
      ungroup()
  }
  
  
  if (remove.id == T){
    outpoly@data <- outpoly@data %>% select(- !!sym(id) )
  }
  
  return(outpoly)  
}