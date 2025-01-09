#### ContinuityCalc.R
#### Rudolf Cesaretti, 5/31/2022

#### "ContinuityCalc" calculates the inputs to forward and 
#### backward continuity variables

pak <- c("rgdal", "sp", "sf", "GISTools", "lwgeom", "tidyverse", "tidyr", "data.table", "zoo")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)

######################################################################
##########################  ContinuityCalc  ##########################
######################################################################


ContinuityCalc = function(sf1, # sf package class polygon object
                          sf2, # sf package class polygon object
                          per, # period name
                          mode # either "forward" or "backward"
                          ){
  
  tmp <- st_intersection(sf1, sf2)
  y = colnames(tmp)
  zz = (ncol(tmp) - 1)/2
  y1 = paste0(y[1:zz],".1")
  y2 = paste0(gsub('.{2}$', '', y[(zz+1):(zz*2)]),".2")
  colnames(tmp)<-c(y1,y2,"geometry")
  st_geometry(tmp) <- "geometry"
  tmp$Area_ha <- as.numeric(st_area(tmp)*0.0001)
  tmp$Period <- per
  tmp$CerPhase <- per
  tmp$PeriodType = "Overlap"
  tmp$SurvReg <- tmp$SurvReg.1
  tmp <- tmp[order(tmp$SurvReg),]
  tmp <- tmp %>% mutate(Number = rowid(SurvReg))
  tmp$Site = paste0(tmp$SurvReg,"-",tmp$Period,"-",tmp$Number)
  tmp$OccSeqLoc <-  tmp$OccSeqLoc.1
  tmp$SherdDens <- tmp$SherdDens.1 + tmp$SherdDens.2
  tmp$Assemb1 <- tmp$SherdDens.1 * (tmp$Area_ha * 10000)
  tmp$Assemb2 <- tmp$SherdDens.2 * (tmp$Area_ha * 10000)
  tmp$PopDens <- tmp$PopDens.1 + tmp$PopDens.2
  tmp$Population <- tmp$PopDens * tmp$Area_ha
  
  if (mode == "forward"){
    fwc = as.data.frame(tmp$Site.1)
    colnames(fwc) = "Site"
    fwc$OvlpSite = tmp$Site
    fwc$PhaseSite1 = tmp$Site.1
    fwc$PhaseSite1_OccSeqLoc = tmp$OccSeqLoc.1
    fwc$PhaseSite1_Pop = tmp$Population.1
    fwc$PhaseSite1_Area = tmp$Area_ha.1
    fwc$OvlpSite_Pop = tmp$Population
    fwc$OvlpSite_Area = tmp$Area_ha
    fwc$OvlpSite_Area = tmp$Area_ha
    fwc$FwOvlp.Assemb = tmp$Assemb1
    fwc2 <- fwc %>% group_by(Site) %>% summarise(
      OccSeqLoc = paste(unique(PhaseSite1_OccSeqLoc),collapse="; "),
      OvlpSites = paste(unique(OvlpSite),collapse=" + "),
      FwOvlp.Sites = paste(unique(OvlpSite,PhaseSite1),collapse=" + "),
      PhaseSite1_Pop = mean(PhaseSite1_Pop),
      PhaseSite1_Area = mean(PhaseSite1_Area),
      FwOvlp.Area = sum(OvlpSite_Area),
      FwOvlp.Pop = sum(OvlpSite_Pop),
      FwOvlp.Assemb = sum(FwOvlp.Assemb))
    
    out <- fwc2
  }
  
  if (mode == "backward"){
    bwc = as.data.frame(tmp$Site.2)
    colnames(bwc) = "Site"
    bwc$OvlpSite = tmp$Site
    bwc$PhaseSite2 = tmp$Site.2
    bwc$PhaseSite2_OccSeqLoc = tmp$OccSeqLoc.2
    bwc$PhaseSite2_Pop = tmp$Population.2
    bwc$PhaseSite2_Area = tmp$Area_ha.2
    bwc$OvlpSite_Pop = tmp$Population
    bwc$OvlpSite_Area = tmp$Area_ha
    bwc$BwOvlp.Assemb = tmp$Assemb2
    bwc2 <- bwc %>% group_by(Site) %>% summarise(
      OccSeqLoc = paste(unique(PhaseSite2_OccSeqLoc),collapse="; "),
      OvlpSites = paste(unique(OvlpSite),collapse=" + "),
      BwOvlp.Sites = paste(unique(OvlpSite,PhaseSite2),collapse=" + "),
      PhaseSite2_Pop = mean(PhaseSite2_Pop),
      PhaseSite2_Area = mean(PhaseSite2_Area),
      BwOvlp.Area = sum(OvlpSite_Area),
      BwOvlp.Pop = sum(OvlpSite_Pop),
      BwOvlp.Assemb = sum(BwOvlp.Assemb))
    
    out <- bwc2
  }
  
  return(out)
}