#### OverlapPoly.R
#### Rudolf Cesaretti, 5/31/2022

#### "OverlapPoly" creates a new polygon file for the overlap of
#### sites in adjacent phases, and merges their attribute data 
#### together 

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

###################################################################
##########################  OverlapPoly  ##########################
###################################################################


OverlapPoly = function(sf1, # sf package class polygon object
                       sf2, # sf package class polygon object
                       per, # new overlap period name
                       ord  # ordering of output variables
                       ){
  
  tmp <- st_intersection(sf1, sf2)
  y = colnames(tmp)
  zz = (ncol(tmp) - 1)/2
  y1 = paste0(y[1:zz],".1")
  y2 = paste0(gsub('.{2}$', '', y[(zz+1):(zz*2)]),".2")
  colnames(tmp)<-c(y1,y2,"geometry")
  st_geometry(tmp) <- "geometry"
  coordzz <- centroid_coords(as(tmp, "Spatial"))
  tmp$East  <- coordzz$x
  tmp$North <- coordzz$y
  tmp$Period <- per
  tmp$CerPhase <- per
  tmp$PeriodType = "Overlap"
  tmp$SurvReg <- tmp$SurvReg.1
  tmp <- tmp[order(tmp$SurvReg),]
  tmp <- tmp %>% mutate(Number = rowid(SurvReg))
  tmp$Site = paste0(tmp$SurvReg,"-",tmp$Period,"-",tmp$Number)
  tmp$Area_ha <- as.numeric(st_area(tmp)*0.0001)
  tmp$Perim_m2 <- as.numeric(st_perimeter(tmp))
  tmp$OccSeqLoc <-  tmp$OccSeqLoc.1
  tmp$SherdDens <- tmp$SherdDens.1 + tmp$SherdDens.2
  tmp$Tot.Assemb <- tmp$SherdDens * tmp$Area_ha * 10000
  tmp$PopDens <- tmp$PopDens.1 + tmp$PopDens.2
  tmp$Population <- tmp$PopDens * tmp$Area_ha
  tmp$ComponentNum <- 1
  tmp$ComponentSites <- tmp$Site
  tmp$FwOvlp.Area = tmp$Area_ha
  tmp$FwOvlp.Pop = tmp$Population
  tmp$BwOvlp.Area = tmp$Area_ha
  tmp$BwOvlp.Pop = tmp$Population
  tmp$BwOvlp.Assemb <- tmp$SherdDens.1 * (tmp$Area_ha * 10000)
  tmp$FwOvlp.Assemb <- tmp$SherdDens.2 * (tmp$Area_ha * 10000)
  tmp$Net.Assemb = 0
  tmp$AreaBwCont <- 1
  tmp$AreaFwCont <- 1
  tmp$PopBwCont <- 1
  tmp$PopFwCont <- 1
  tmp <- tmp  %>% rowwise() %>% mutate(
    UrbanPop = ifelse(Population > UrbanThresh, Population - UrbanThresh, 0),
    RuralPop = ifelse(Population < UrbanThresh, Population, UrbanThresh)) %>% ungroup() %>% mutate(             
      UrbanScale = Population / UrbanThresh,
      PctUrban = UrbanPop / Population,
      PctRural = RuralPop / Population)
  tmp$FwOvlp.Sites <- paste(tmp$Site, tmp$Site.2, sep=" + ")
  tmp$BwOvlp.Sites <- paste(tmp$Site.1, tmp$Site, sep=" + ")
  
  NAVARS = c("PeriodLength", "PeriodNum", "PeriodBegin", "PeriodEnd", "OccSeqLoc.Sites", 
             "SubOccSeqLoc", "SubOccSeqLoc.Sites", "Occ.EF", "Occ.EF_MF", "Occ.MF", 
             "Occ.MF_LF", "Occ.LF", "Occ.LF_TF", "Occ.TF", "Occ.TF_CL", "Occ.CL", 
             "Occ.CL_ET", "Occ.ET", "Occ.ET_LTAzI", "Occ.LTAzI", "Occ.LTAzI_EA", 
             "Occ.EA", "Occ.EA_LA", "Occ.LA", "Occ.TOT", "SubOcc.EF", "SubOcc.EF_MF", 
             "SubOcc.MF", "SubOcc.MF_LF", "SubOcc.LF", "SubOcc.LF_TF", "SubOcc.TF", 
             "SubOcc.TF_CL", "SubOcc.CL", "SubOcc.CL_ET", "SubOcc.ET", "SubOcc.ET_LTAzI", 
             "SubOcc.LTAzI", "SubOcc.LTAzI_EA", "SubOcc.EA", "SubOcc.EA_LA", "SubOcc.LA", 
             "SubOcc.TOT", "Found", "FoundInit", "Abandon", "Persist", "DewarType", "OccuIntertia")
  
  for (j in 1:length(NAVARS)) {
    new <- rep(NA, nrow(tmp))
    tmp[ , ncol(tmp) + 1] <- new
    colnames(tmp)[ncol(tmp)] <- paste0(NAVARS[j])
  }
  
  PASTEVARS = c("M_Sites", "M_SiteCode", "M_SiteName", "M_FieldSite.Region", 
                "M_FieldSite.Period", "M_SurveyYearNumber", "M_Supervisor", "M_Map", "O_Elev", 
                "O_ElevMed", "O_ElevMin", "O_ElevMax", "O_EZcode", "O_EnvironmentalZone", 
                "O_Soil", "O_SoilMed", "O_SoilMin", "O_SoilMax", "O_Erosion", "O_ErosionMed",                  
                "O_ErosionMin", "O_ErosionMax", "O_ModernUse", "O_ModernSettlement",                  
                "O_Rainfall", "O_Area", "O_MoundDomestic", "O_MoundCeremonial",                    
                "O_MoundQuestionable", "O_MoundTotal", "O_MoundRecorded", "O_DMoundArea",                    
                "O_Architecture", "O_TerraceConfidence", "O_TerraceExtent", "O_Sherd",                  
                "O_SherdMed", "O_SherdMin", "O_SherdMax", "O_Rubble", "O_RubbleMed", "O_RubbleMin", 
                "O_RubbleMax", "O_Population", "O_PopMin", "O_PopMax", "O_PopMethod", "O_stcode", 
                "O_SiteType", "O_SubPeriod1", "O_SubPeriod2", "O_OccEF", "O_OccMF", "O_OccLF", 
                "O_OccTF", "O_OccCL", "O_OccEC", "O_OccMC", "O_OccLC", "O_OccET", "O_OccLT", 
                "O_OccAZ", "O_OccEA", "O_OccLA", "O_OccTot", "O_OccSeqLoc", "O_SubOc1", 
                "O_SubOc2", "O_PdDupSite", "O_Group", "O_Comments")
  
  for (j in 1:length(PASTEVARS)) {
    z <- paste0(PASTEVARS[j])
    z1 <- paste0(PASTEVARS[j],".1")
    z2 <- paste0(PASTEVARS[j],".2")
    z3 <- c(z1,z2)
    z4 <- c(z,z1,z2)
    zz = tmp %>% select(!!!syms(z3)) 
    zz = as.data.frame(st_drop_geometry(zz))
    zz$g <- 1:nrow(zz)
    zz2 = zz %>% group_by(g) %>% summarise(
      VARIA = paste(unique(!!!syms(z3)), collapse="; ")) %>% select(VARIA)
    tmp[ , ncol(tmp) + 1] <- zz2$VARIA
    colnames(tmp)[ncol(tmp)] <- paste0(PASTEVARS[j])
  }
  
  tmp2 = as(tmp, "Spatial")
  
  tmp2@data <- tmp2@data %>% select(!!!syms(ord))
  
  tmp2@data$rownames <- tmp2@data$Site
  tmp2@data <- column_to_rownames(tmp2@data, var = "rownames")
  #rownames(tmp2@data) <- tmp2@data$Site
  tmp2 <- spChFIDs(tmp2, tmp2$Site)
  return(tmp2) 
  
}


