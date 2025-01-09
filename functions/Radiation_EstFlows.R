#### Radiation_EstFlows.R
#### Rudolf Cesaretti, 6/31/2022

#### "Radiation_EstFlows" 
#### 
#### 
#### 
#### https://book.archnetworks.net/spatialinteraction

pak <- c("compiler", "msm", "snowfall", "parallel", "doParallel", "rgdal", "sp", 
         "sf", "GISTools", "lwgeom", "tidyverse", "tidyr", "data.table", "zoo", 
         "scales", "igraph", "tidygraph", "centiserve", "CINNA", "sfnetworks")
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
library(stars)
library(modelsummary)
library(cowplot)

# Multicore setup
cl <- makeCluster(detectCores())
registerDoParallel(cl)

###############################################################
##########################  Radiation_EstFlows  #########################
###############################################################



#df = Poly_List[[17]]@data
#var = "Population.s2"
#d_mat = CD.mats[[17]]
#xcoords = Poly_List[[17]]@data$East
#ycoords = Poly_List[[17]]@data$North
#crs_coords = 26914
#extended = T #T
#alpha = 1.2
#scale = "variant" #c("invariant", "variant", "none")
#commuters = "var" #c("var", "input", "none")
#Ti_input = NULL
#outputs = c("adjmat", "sflines", "sfpts_netstats")#, "sfnetwork"



Radiation_EstFlows <- function(df,
                      var,
                      d_mat, 
                      pts.sf,
                      xcoords,
                      ycoords,
                      crs_coords,
                      extended = FALSE,
                      alpha = NULL,
                      scale = c("invariant", "variant", "none"),
                      commuters = c("var", "input", "none"),
                      Ti_input = NULL,
                      outputs = c("adjmat", "igraph", "sflines", "sfpts_netstats", "netstats", "sfnetwork")
                      ){
  v <- as.numeric(df[,var])
  
  if (commuters == "var"){T_i <- v}
  if (commuters == "input"){T_i <- Ti_input}
  if (commuters == "none"){T_i <- rep(1, length(v))}
  if (scale == "invariant" | scale == "none"){scalevar <- rep(1, length(v))}
  if (scale == "variant"){scalevar <- (1 - v/sum(v))}
  
  ## create square matrix with rows and columns for every site
  #out <-matrix(0, length(pop), length(pop))
  d_mat2 <- d_mat
  diag(d_mat2) <- NA
  df1 <- df2 <- data.frame(AggSite = rownames(d_mat2), x1 = xcoords, y1 = ycoords, m = v)
  colnames(df2) <- c("AggSite","x2","y2", "n")
  Ti.df <- data.frame(AggSite = rownames(d_mat2), Ti = (T_i/scalevar))
  
  d_mat_df <- as.data.frame(d_mat2)
  s_df <- d_mat_df <- d_mat_df %>% rownames_to_column(var = "Site") %>% 
    tidyr::pivot_longer(!Site,names_to="Site2",values_to ="r_ij") %>% 
    left_join(df1, by = c("Site" = "AggSite")) %>% 
    left_join(df2, by = c("Site2" = "AggSite")) %>% 
    left_join(Ti.df, by = c("Site" = "AggSite")) %>% rowwise() %>%
    mutate(self = ifelse(is.na(r_ij), 0, 1)) %>% ungroup()
    
  s_df$sel_circle <- sel_circle <- unlist((foreach(i=1:nrow(d_mat2)) %:% foreach(j=1:ncol(d_mat2)) %dopar% (names(which(d_mat2[i,-j] <= d_mat2[i,j]))) ), recursive = F)
  s_df$s <- unlist((foreach(i=1:length(sel_circle)) %dopar% (sum(df2[df2$AggSite %in% sel_circle[[i]],4], na.rm=T)) ))
  s_df$n_s <- unlist(lapply(sel_circle, length))
  sel_circle <- list()
  s_df <- s_df %>% select(-sel_circle) %>%  mutate(
                          P = self*((m * n) / ((m + s) * (m + n + s))),
                          Tij = P*Ti)
  
  if (extended == T){
    s_df$a <- alpha
    s_df <- s_df %>% mutate(P = self*((((m + n + s)^a) - ((m + s)^a))*((m^a) + 1))/((((m + s)^a) + 1)*(((m + n + s)^a) + 1)),
                            Tij = P*Ti)
  }
  
  Tij_adjmat <- matrix(s_df$Tij, nrow=nrow(d_mat), ncol=ncol(d_mat), dimnames= list(rownames(d_mat), colnames(d_mat)))
  
  Tij_ig <- graph_from_adjacency_matrix(Tij_adjmat,mode = "directed",
                                           weighted=T,add.rownames = NULL)
  
  netstat <- net_stats2(ig=Tij_ig,ID=rownames(d_mat2),prefix="Rad_",direc=T, weights=T)
  
  rename_geom_col <- function(x, new_name){
    current_name  <-  attr(x, "sf_column")
    names(x)[names(x)==current_name]  <- new_name
    st_geometry(x) <- new_name
    return(x)
  }
  
  pts1 <- mapply(function(x,y) st_point(c(x,y)),s_df$x1,s_df$y1,SIMPLIFY = FALSE)
  pts2 <- mapply(function(x,y) st_point(c(x,y)),s_df$x2,s_df$y2,SIMPLIFY = FALSE)
  lines <- mapply(function(p1,p2) st_as_sf(st_sfc(st_linestring(c(p1,p2))), crs = crs_coords) ,pts1,pts2,SIMPLIFY = FALSE)
  
  
  #linesdf <- do.call(rbind, lines)
  linesdf <- data.table::rbindlist(lines)
  colnames(linesdf) <- "geometry"
  #linesdf <- rename_geom_col(linesdf, "geometry")
  
  colnames(df1) <- c("AggSite","x","y", paste(var))
  sfpts <- st_as_sf(df1, coords = c("x", "y"), crs = crs_coords, agr = "constant")
  sfpts <- cbind(sfpts,netstat[,-1])
  
  sdf <- s_df %>% select(-c(x1,y1,x2,y2)) %>%  rename(from = Site, to = Site2)
  
  sflines <- cbind(sdf,linesdf)
  sflines <- st_as_sf(sflines, crs = crs_coords, agr = "constant")
  
  outlist <- list()
  q=1
  if ("adjmat" %in% outputs){outlist[[q]] <- Tij_adjmat ; q=q+1}
  if ("igraph" %in% outputs){outlist[[q]] <- Tij_ig ; q=q+1}
  if ("sflines" %in% outputs){outlist[[q]] <- sflines ; q=q+1}
  if ("sfpts_netstats" %in% outputs){outlist[[q]] <- sfpts ; q=q+1}
  if ("netstats" %in% outputs){outlist[[q]] <- netstat ; q=q+1}
  
  if ("sfnetwork" %in% outputs){
    ig_sf = sfnetwork(sfpts, sflines)
     
    ig_sf <- ig_sf %>%
      activate("edges") %>%
        mutate(weight = Tij) %>%
      activate("nodes") %>%
        mutate(weight = !!sym(var))
    outlist[[q]] <- ig_sf ; q=q+1
    }
  names(outlist) <- outputs
  return(outlist)
}
