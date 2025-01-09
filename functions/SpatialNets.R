#### SpatialNets.R
#### Rudolf Cesaretti, 7/4/2022

#### "SpatialNets" 
#### 
#### 
#### 
#### 

pak <- c("rgdal", "rgeos", "sp", "sf", "GISTools", "raster", "Matrix", "gdistance", "lwgeom", "tidyverse",
         "tidyr", "stars", "dismo", "purrr", "spatialEco", "whitebox", "classInt", "deldir", "spdep",
         "igraph", "dbscan", "cowplot", "deldir", "cccd", "ggraph", "geosphere", "statnet", "intergraph",
         "ggnewscale")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)

#library(scales)
#library(pracma)
#library(modelsummary)

###############################################################
########################  SpatialNets  ########################
###############################################################


SpatialNets <- function(CD.mat, #CD.mats[[p]]/3600
                        nodes, #Pts_List[[p]]@data
                        coords, #coordinates(Pts_List[[p]])
                        Pts.sf, #st_as_sf(Pts_List[[p]])
                        Sites.sf, #st_as_sf(Poly_List[[p]])
                        Catch.sf, #st_as_sf(Catch_List[[p]])
                        limits, #CatchLims
                        bounds = extent(limits),
                        window = c(bounds@xmin, bounds@xmax, bounds@ymin, bounds@ymax),
                        ggplots = c("MultiPlot","PlotResults", "none"),
                        ggbasemap, #ggbasemap
                        mode = c("single", "multi", "multicombine"),
                        method = c("gabriel", "rnn", "delaunay", "soi", "poly", "knn", "maxdist"),
                        combine = c("weighted","unweighted"),
                        plotNodeStat = c("centr_deg_n","centr_btw_n","centr_clos_n","centr_hrmo_n",
                                         "centr_hub_n","centr_pwr_n","centr_auth_n","centr_pgrk_n",
                                         "centr_alph_n","centr_avg_n","trans_loc_n"),
                        ggPlotTitle = "Spatial Network Models SBOM",
                        knn_k = 4,
                        maxdist = 1,
                        outputs = c("net_distmat", "netstats_df", "igraph", 
                                    "igraph_wt", "lines_sf", "adjmat", "adjmatCD",
                                    "nb")){
  
  outlist <- list()
  MultiPlotList <- list()
  o <- rownames(CD.mat)
  
  if ("net_distmat" %in% outputs){net_distmat_list <- list() ; ia = 1}
  if ("netstats_df" %in% outputs){netstats_df_list <- list() ; ib = 1}
  if ("igraph" %in% outputs){igraph_list <- list() ; ic = 1}
  if ("igraph_wt" %in% outputs){igraph_wt_list <- list() ; id = 1}
  if ("lines_sf" %in% outputs){lines_sf_list <- list() ; ie = 1}
  if ("adjmat" %in% outputs){adjmat_list <- list() ; ih = 1}
  if ("adjmatCD" %in% outputs){adjmatCD_list <- list() ; ii = 1}
  if ("nb" %in% outputs){nb_list <- list() ; ij = 1}
  
  #singlePlotResults <- (method == "single" & ggplots == "PlotResults")
  #mp <- (ggplots == "MultiPlot")
  
  ggpltz <- function(ggbasemap, netstats_df, CatchLims.sf=CatchLims.sf, 
                     poly=Catch.sf, subtitle, NodeStat=plotNodeStat, 
                     lines_sf, geomet="geometry", s=c(1,5),wt=F,subtit=T){
    netstats.sf <- st_as_sf(x = netstats_df, coords = c("East", "North"),crs = st_crs(CatchLims.sf))
    
    ggplt <- ggbasemap + ggnewscale::new_scale_fill()+
      geom_sf(data = CatchLims.sf, aes(geometry = geometry), color="black", size=2, alpha = 0) +
      geom_sf(data = poly, aes(geometry = !!sym(geomet)),color="black", size=0.35, alpha = 0) +
      {if(wt == FALSE) geom_sf(data = lines_sf, aes(geometry = geometry),color="red2", size=1)}+
      #{if(wt == TRUE) geom_sf(data = lines_sf, aes(geometry = geometry, color=wt), size=1)}+
      {if(wt == TRUE) geom_segment(data = lines_sf, aes(x = X1,y = Y1,xend = X2,yend = Y2, color=wt), size=1)}+
      {if(wt == TRUE) scale_color_gradient(low = "white", high = "red")}+
      {if(wt == TRUE) guides(color="none")}+
      geom_sf(data = netstats.sf, aes(geometry = geometry, size=(!!sym(NodeStat))),
              color="black", shape=21,fill="firebrick", alpha=0.8) +
      scale_size(range = s)+
      coord_sf() + theme_void()+
      {if(subtit == TRUE) labs(subtitle = subtitle, size=NodeStat)}+
      {if(subtit == FALSE) labs(size=NodeStat)}+
      theme(plot.subtitle = element_text(hjust = 0.5, face="bold", size=14),
            plot.background = element_rect(fill = "white", colour = NA),
            legend.position=c(0.15,0.3), legend.spacing.y = unit(0.1,"line"),
            legend.text=element_text(size=8, face="bold"),
            legend.title=element_text(size=10, face="bold"),
            legend.key = element_rect(colour = "transparent", fill = "white"), 
            legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
    return(ggplt)
  }
  j=1
  
  if ("gabriel" %in% method){
    #Gabriel Graph
    gab_ig <- cccd::gg(CD.mat)
    V(gab_ig)$vertex.names <- V(gab_ig)$name <- row.names(CD.mat)
    gab_am <- as.matrix(as_adj(gab_ig,type = "both",attr = NULL,names = TRUE))
    gab_am <- gab_am[o, o]
    gab_amCD <- CD.mat*gab_am
    gab_ig_wt <- graph_from_adjacency_matrix(gab_amCD,mode = "undirected",
                                             weighted=T,add.colnames = NULL,
                                             add.rownames = NULL)
    gab_nb <- mat2listw(gab_am,style="B")
    gab_nb <- gab_nb$neighbours
    gab_lines_sf <- nb2lines(gab_nb, coords=coords, as_sf=T)
    st_crs(gab_lines_sf) <- st_crs(CatchLims.sf)
    gab_path_CDmat <- distances(gab_ig_wt,  mode = "all")
    gab_netstats = net_stats(gab_ig,o,"gabriel",coords)
    
    if ("net_distmat" %in% outputs){net_distmat_list[[ia]] <- gab_path_CDmat ; ia=ia+1}
    if ("netstats_df" %in% outputs){netstats_df_list[[ib]]<-gab_netstats ; ib=ib+1}
    if ("igraph" %in% outputs){igraph_list[[ic]]<-gab_ig ; ic=ic+1}
    if ("igraph_wt" %in% outputs){igraph_wt_list[[id]]<-gab_ig_wt ; id=id+1}
    if ("lines_sf" %in% outputs){lines_sf_list[[ie]]<-gab_lines_sf ; ie=ie+1}
    if ("adjmat" %in% outputs){adjmat_list[[ih]]<-gab_am ; ih=ih+1}
    if ("adjmatCD" %in% outputs){adjmatCD_list[[ii]]<-gab_amCD ; ii=ii+1}
    if ("nb" %in% outputs){nb_list[[ij]]<-gab_nb ; ij=ij+1}
    
    if (ggplots == "MultiPlot"){
      gg_gab <- ggpltz(ggbasemap=ggbasemap, netstats_df=gab_netstats, 
                       CatchLims.sf=CatchLims.sf,poly=Catch.sf, s=c(0.75,1.5),
                       subtit=F,#subtitle = "Gabriel Graph",
                       lines_sf = gab_lines_sf, geomet="geometry")
      MultiPlotList[[j]] <- gg_gab+theme(legend.position="none",
                                         plot.margin = unit(c(0,0,0,0), "cm"))
      j=j+1
    }
    if (mode == "single" & ggplots == "PlotResults"){
      ggout <- ggpltz(ggbasemap=ggbasemap, netstats_df=gab_netstats, 
                      CatchLims.sf=CatchLims.sf,poly=Catch.sf, 
                      subtitle = "Gabriel Graph", 
                      lines_sf = gab_lines_sf, geomet="geometry")
      ggout <- ggout  + labs(title=ggPlotTitle) + 
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
    }
  }
  
  
  if ("rnn" %in% method){
    #Relative Neighborhood Network
    rnn_ig <- cccd::gg(CD.mat)
    V(rnn_ig)$vertex.names <- V(rnn_ig)$name <- row.names(CD.mat)
    rnn_am <- as.matrix(as_adj(rnn_ig,type = "both",attr = NULL,names = TRUE))
    rnn_am <- rnn_am[o, o]
    rnn_amCD <- CD.mat*rnn_am
    rnn_ig_wt <- graph_from_adjacency_matrix(rnn_amCD,mode = "undirected",
                                             weighted=T,add.colnames = NULL,
                                             add.rownames = NULL)
    rnn_nb <- mat2listw(rnn_am,style="B")
    rnn_nb <- rnn_nb$neighbours
    rnn_lines_sf <- nb2lines(rnn_nb, coords=coords, as_sf=T)
    st_crs(rnn_lines_sf) <- st_crs(CatchLims.sf)
    rnn_path_CDmat <- distances(rnn_ig_wt,  mode = "all")
    rnn_netstats = net_stats(rnn_ig,o,"rnn",coords)
    
    if ("net_distmat" %in% outputs){net_distmat_list[[ia]] <- rnn_path_CDmat ; ia=ia+1}
    if ("netstats_df" %in% outputs){netstats_df_list[[ib]]<-rnn_netstats ; ib=ib+1}
    if ("igraph" %in% outputs){igraph_list[[ic]]<-rnn_ig ; ic=ic+1}
    if ("igraph_wt" %in% outputs){igraph_wt_list[[id]]<-rnn_ig_wt ; id=id+1}
    if ("lines_sf" %in% outputs){lines_sf_list[[ie]]<-rnn_lines_sf ; ie=ie+1}
    if ("adjmat" %in% outputs){adjmat_list[[ih]]<-rnn_am ; ih=ih+1}
    if ("adjmatCD" %in% outputs){adjmatCD_list[[ii]]<-rnn_amCD ; ii=ii+1}
    if ("nb" %in% outputs){nb_list[[ij]]<-rnn_nb ; ij=ij+1}
    
    if (ggplots == "MultiPlot"){
      gg_rnn <- ggpltz(ggbasemap=ggbasemap, netstats_df=rnn_netstats, 
                       CatchLims.sf=CatchLims.sf,poly=Catch.sf, s=c(0.75,1.5),
                       subtit=F,#subtitle = "Relative Neighborhood", 
                       lines_sf = rnn_lines_sf, geomet="geometry")
      MultiPlotList[[j]] <- gg_rnn+theme(legend.position="none",
                                         plot.margin = unit(c(0,0,0,0), "cm"))
      j=j+1
    }
    if (mode == "single" & ggplots == "PlotResults"){
      ggout <- ggpltz(ggbasemap=ggbasemap, netstats_df=rnn_netstats, 
                      CatchLims.sf=CatchLims.sf,poly=Catch.sf, 
                      subtitle = "Relative Neighborhood", 
                      lines_sf = rnn_lines_sf, geomet="geometry")
      ggout <- ggout  + labs(title=ggPlotTitle) + 
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
    }
  }
  
  if ("delaunay" %in% method){
    #Delaunay Triangulation
    box <- st_sfc(bbox_polygon(CatchLims.sf))
    vor <- st_voronoi(st_union(Pts.sf), box)
    x=st_intersection(st_cast(vor), CatchLims.sf)
    x=x %>% st_cast("MULTIPOLYGON")
    xx=st_sf(x, agr="identity", a=3)
    vor_sf = st_join(xx,Pts.sf)
    rownames(vor_sf)<-vor_sf$AggSite
    vor_sp=as(vor_sf,"Spatial")
    
    del_nb <- poly2nb(vor_sp)
    del_lines_sf <- nb2lines(del_nb, coords=vor_sp, as_sf=T)
    st_crs(del_lines_sf) <- st_crs(CatchLims.sf)
    del_am <- nb2mat(del_nb, style="B")
    colnames(del_am) <- rownames(del_am)
    del_am <- del_am[o, o]
    del_amCD <- CD.mat*del_am
    del_ig_wt <- graph_from_adjacency_matrix(del_amCD,mode = "undirected",
                                             weighted=T,add.colnames = NULL,
                                             add.rownames = NULL)
    del_ig <- graph_from_adjacency_matrix(del_am,mode = "undirected", 
                                          add.colnames = NA,add.rownames = NULL)
    del_path_CDmat <- distances(del_ig_wt,  mode = "all")
    del_netstats = net_stats(del_ig,o,"delaunay",coords)
    
    if ("net_distmat" %in% outputs){net_distmat_list[[ia]] <- del_path_CDmat ; ia=ia+1}
    if ("netstats_df" %in% outputs){netstats_df_list[[ib]]<-del_netstats ; ib=ib+1}
    if ("igraph" %in% outputs){igraph_list[[ic]]<-del_ig ; ic=ic+1}
    if ("igraph_wt" %in% outputs){igraph_wt_list[[id]]<-del_ig_wt ; id=id+1}
    if ("lines_sf" %in% outputs){lines_sf_list[[ie]]<-del_lines_sf ; ie=ie+1}
    if ("adjmat" %in% outputs){adjmat_list[[ih]]<-del_am ; ih=ih+1}
    if ("adjmatCD" %in% outputs){adjmatCD_list[[ii]]<-del_amCD ; ii=ii+1}
    if ("nb" %in% outputs){nb_list[[ij]]<-del_nb ; ij=ij+1}
    
    if (ggplots == "MultiPlot"){
      gg_del <- ggpltz(ggbasemap=ggbasemap, netstats_df=del_netstats, 
                       CatchLims.sf=CatchLims.sf,poly=Catch.sf, s=c(0.75,1.5),
                       subtit=F,#subtitle = "Delaunay Triangulation", 
                       lines_sf = del_lines_sf, geomet="geometry")
      MultiPlotList[[j]] <- gg_del+theme(legend.position="none",
                                         plot.margin = unit(c(0,0,0,0), "cm"))
      j=j+1
    }
    if (mode == "single" & ggplots == "PlotResults"){
      ggout <- ggpltz(ggbasemap=ggbasemap, netstats_df=del_netstats, 
                      CatchLims.sf=CatchLims.sf,poly=Catch.sf, 
                      subtitle = "Delaunay Triangulation", 
                      lines_sf = del_lines_sf, geomet="geometry")
      ggout <- ggout  + labs(title=ggPlotTitle) + 
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
    }
  }
  
  if ("soi" %in% method){
    #Sphere of Influence Graph (SOI)
    x <- tri2nb(coords, row.names= nodes$AggSite)
    soi_nb <- graph2nb(soi.graph(x, coords, quadsegs=50), row.names=nodes$AggSite)
    soi_lines_sf <- nb2lines(soi_nb, coords=coords, as_sf=T)
    st_crs(soi_lines_sf) <- st_crs(CatchLims.sf)
    soi_am <- nb2mat(soi_nb, style="B")
    colnames(soi_am) <- rownames(soi_am)
    soi_am <- soi_am[o, o]
    soi_amCD <- soi_am*CD.mat
    soi_ig_wt <- graph_from_adjacency_matrix(soi_amCD,mode = "undirected",
                                             weighted=T,add.colnames = NULL,
                                             add.rownames = NULL)
    soi_ig <- graph_from_adjacency_matrix(soi_am,mode = "undirected", 
                                          add.colnames = NA,add.rownames = NULL)
    soi_path_CDmat <- distances(soi_ig_wt,  mode = "all")
    soi_netstats = net_stats(soi_ig,o,"soi",coords)
    
    if ("net_distmat" %in% outputs){net_distmat_list[[ia]] <- soi_path_CDmat ; ia=ia+1}
    if ("netstats_df" %in% outputs){netstats_df_list[[ib]]<-soi_netstats ; ib=ib+1}
    if ("igraph" %in% outputs){igraph_list[[ic]]<-soi_ig ; ic=ic+1}
    if ("igraph_wt" %in% outputs){igraph_wt_list[[id]]<-soi_ig_wt ; id=id+1}
    if ("lines_sf" %in% outputs){lines_sf_list[[ie]]<-soi_lines_sf ; ie=ie+1}
    if ("adjmat" %in% outputs){adjmat_list[[ih]]<-soi_am ; ih=ih+1}
    if ("adjmatCD" %in% outputs){adjmatCD_list[[ii]]<-soi_amCD ; ii=ii+1}
    if ("nb" %in% outputs){nb_list[[ij]]<-soi_nb ; ij=ij+1}
    
    if (ggplots == "MultiPlot"){
      gg_soi <- ggpltz(ggbasemap=ggbasemap, netstats_df=soi_netstats, 
                       CatchLims.sf=CatchLims.sf,poly=Catch.sf, s=c(0.75,1.5),
                       subtit=F,#subtitle = "Sphere of Influence", 
                       lines_sf = soi_lines_sf, geomet="geometry")
      MultiPlotList[[j]] <- gg_soi+theme(legend.position="none",
                                         plot.margin = unit(c(0,0,0,0), "cm"))
      j=j+1
    }
    if (mode == "single" & ggplots == "PlotResults"){
      ggout <- ggpltz(ggbasemap=ggbasemap, netstats_df=soi_netstats, 
                      CatchLims.sf=CatchLims.sf,poly=Catch.sf, 
                      subtitle = "Sphere of Influence", 
                      lines_sf = soi_lines_sf, geomet="geometry")
      ggout <- ggout  + labs(title=ggPlotTitle) + 
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
    }
  }
  
  
  if ("poly" %in% method){
    # Catchment Neighbors
    catch_nb <- poly2nb(Catch_List[[p]], row.names= Catch_List[[p]]$AggSite)
    catch_lines_sf <- nb2lines(catch_nb, coords=Catch_List[[p]], as_sf=T)
    st_crs(catch_lines_sf) <- st_crs(CatchLims.sf)
    catch_am <- nb2mat(catch_nb, style="B", zero.policy	=T)
    colnames(catch_am) <- rownames(catch_am)
    catch_am <- catch_am[o, o]
    catch_amCD <- CD.mat*catch_am
    catch_ig_wt <- graph_from_adjacency_matrix(catch_amCD,mode = "undirected",
                                               weighted=T,add.colnames = NULL,
                                               add.rownames = NULL)
    catch_ig <- graph_from_adjacency_matrix(catch_am,mode = "undirected", 
                                            add.colnames = NULL,add.rownames = NULL)
    catch_path_CDmat <- distances(catch_ig_wt,  mode = "all")
    catch_netstats = net_stats(catch_ig,o,"catch",coords)
    
    if ("net_distmat" %in% outputs){net_distmat_list[[ia]] <- catch_path_CDmat ; ia=ia+1}
    if ("netstats_df" %in% outputs){netstats_df_list[[ib]]<-catch_netstats ; ib=ib+1}
    if ("igraph" %in% outputs){igraph_list[[ic]]<-catch_ig ; ic=ic+1}
    if ("igraph_wt" %in% outputs){igraph_wt_list[[id]]<-catch_ig_wt ; id=id+1}
    if ("lines_sf" %in% outputs){lines_sf_list[[ie]]<-catch_lines_sf ; ie=ie+1}
    if ("adjmat" %in% outputs){adjmat_list[[ih]]<-catch_am ; ih=ih+1}
    if ("adjmatCD" %in% outputs){adjmatCD_list[[ii]]<-catch_amCD ; ii=ii+1}
    if ("nb" %in% outputs){nb_list[[ij]]<-catch_nb ; ij=ij+1}
    
    if (ggplots == "MultiPlot"){
      gg_catch <- ggpltz(ggbasemap=ggbasemap, netstats_df=catch_netstats, 
                         CatchLims.sf=CatchLims.sf,poly=Catch.sf, s=c(0.75,1.5),
                         subtit=F,#subtitle = "Catchment Poly Net", 
                         lines_sf = catch_lines_sf, geomet="geometry")
      MultiPlotList[[j]] <- gg_catch+theme(legend.position="none",
                                           plot.margin = unit(c(0,0,0,0), "cm"))
      j=j+1
    }
    if (mode == "single" & ggplots == "PlotResults"){
      ggout <- ggpltz(ggbasemap=ggbasemap, netstats_df=catch_netstats, 
                      CatchLims.sf=CatchLims.sf,poly=Catch.sf, 
                      subtitle = "Catchment Poly Net", 
                      lines_sf = catch_lines_sf, geomet="geometry")
      ggout <- ggout  + labs(title=ggPlotTitle) + 
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
    }
  }
  
  
  if ("knn" %in% method){
    #K Nearest Neighbors (KNN)
    knn_ig <- nng(CD.mat, k = knn_k)
    V(knn_ig)$vertex.names <- V(knn_ig)$name <- row.names(CD.mat)
    knn_am <- as.matrix(as_adj(knn_ig,type = "both",attr = NULL,names = TRUE))
    knn_am <- knn_am[o, o]
    knn_amCD <- CD.mat*knn_am
    knn_ig_wt <- graph_from_adjacency_matrix(knn_amCD,mode = "max",
                                             weighted=T,add.colnames = NULL,
                                             add.rownames = NULL)
    #edgelist <- get.edgelist(knn_ig)
    #x=as.data.frame(as_edgelist(knn_ig))
    #edges <- as.data.frame(matrix(NA, nrow(edgelist), 4))
    #colnames(edges) <- c("X1", "Y1", "X2", "Y2")
    #for (i in seq_len(nrow(x))) {
    #  edges[i, ] <- c(nodes[which(nodes$AggSite == x[i, 1]), "East"],
    #                  nodes[which(nodes$AggSite == x[i, 1]), "North"],
    #                  nodes[which(nodes$AggSite == x[i, 2]), "East"],
    #                  nodes[which(nodes$AggSite == x[i, 2]), "North"])
    #}
    #xx=cbind(x,edges)
    #ls <- apply(xx, 1, function(x) 
    #{
    #  v <- as.numeric(x[c(3,4,5,6)])
    #  m <- matrix(v, nrow = 2)
    #  return(st_sfc(st_linestring(m), crs = st_crs(CatchLims.sf)))
    #})
    #
    #geometry = Reduce(c, ls)
    #knn_lines_sf=st_sf(geometry, agr="identity", ID=c(1:nrow(xx)))
    #knn_lines_sf$from = xx$V1
    #knn_lines_sf$to = xx$V2
    #colnames(xxx) <- c("a","geometry","from","to")
    #knn_lines_sf <- xxx
    
    knn_nb <- mat2listw(knn_am,style="B")
    knn_nb <- knn_nb$neighbours
    knn_lines_sf <- nb2lines(knn_nb, coords=coords, as_sf=T)
    st_crs(knn_lines_sf) = st_crs(CatchLims.sf)
    knn_path_CDmat <- distances(knn_ig_wt,  mode = "all")
    knn_netstats = net_stats(knn_ig,o,"knn",coords)
    
    if ("net_distmat" %in% outputs){net_distmat_list[[ia]] <- knn_path_CDmat ; ia=ia+1}
    if ("netstats_df" %in% outputs){netstats_df_list[[ib]]<-knn_netstats ; ib=ib+1}
    if ("igraph" %in% outputs){igraph_list[[ic]]<-knn_ig ; ic=ic+1}
    if ("igraph_wt" %in% outputs){igraph_wt_list[[id]]<-knn_ig_wt ; id=id+1}
    if ("lines_sf" %in% outputs){lines_sf_list[[ie]]<-knn_lines_sf ; ie=ie+1}
    if ("adjmat" %in% outputs){adjmat_list[[ih]]<-knn_am ; ih=ih+1}
    if ("adjmatCD" %in% outputs){adjmatCD_list[[ii]]<-knn_amCD ; ii=ii+1}
    if ("nb" %in% outputs){nb_list[[ij]]<-knn_nb ; ij=ij+1}
    
    if (ggplots == "MultiPlot"){
      gg_knn <- ggpltz(ggbasemap=ggbasemap, netstats_df=knn_netstats, 
                       CatchLims.sf=CatchLims.sf,poly=Catch.sf, s=c(0.75,1.5),
                       subtit=F,#subtitle = "K Nearest Neighbors", 
                       lines_sf = knn_lines_sf, geomet="geometry")
      MultiPlotList[[j]] <- gg_knn+theme(legend.position="none",
                                         plot.margin = unit(c(0,0,0,0), "cm"))
      j=j+1
    }
    if (mode == "single" & ggplots == "PlotResults"){
      ggout <- ggpltz(ggbasemap=ggbasemap, netstats_df=knn_netstats, 
                      CatchLims.sf=CatchLims.sf,poly=Catch.sf, 
                      subtitle = "K Nearest Neighbors", 
                      lines_sf = knn_lines_sf, geomet="geometry")
      ggout <- ggout  + labs(title=ggPlotTitle) + 
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
    }
  }
  
  if ("maxdist" %in% method){
    #Maximum Distance Networks
    maxdist_ig <- network(event2dichot(CD.mat, method = "absolute",thresh = maxdist,leq = TRUE),
                          directed = FALSE, vertex.attr = list(rownames(CD.mat)), 
                          vertex.attrnames = list(c("name")))
    maxdist_ig <- asIgraph(maxdist_ig, amap = attrmap(newdf=data.frame(type = "vertex", 
                                                                       fromcls="network", fromattr="name",
                                                                       tocls= "igraph", toattr="name")))
    V(maxdist_ig)$name <- V(maxdist_ig)$vertex.names
    maxdist_am <- as.matrix(as_adj(maxdist_ig,type = "both",attr = NULL,names = TRUE))
    maxdist_am <- maxdist_am[o, o]
    maxdist_amCD <- CD.mat*maxdist_am
    maxdist_ig_wt <- graph_from_adjacency_matrix(maxdist_amCD,mode = "undirected",
                                                 weighted=T,add.colnames = NULL,
                                                 add.rownames = NULL)
    maxdist_nb <- mat2listw(maxdist_am,style="B")
    maxdist_nb <- maxdist_nb$neighbours
    maxdist_lines_sf <- nb2lines(maxdist_nb, coords=coords, as_sf=T)
    st_crs(maxdist_lines_sf) <- st_crs(CatchLims.sf)
    maxdist_path_CDmat <- distances(maxdist_ig_wt,  mode = "all")
    maxdist_netstats = net_stats(maxdist_ig,o,"maxdist",coords)
    
    if ("net_distmat" %in% outputs){net_distmat_list[[ia]] <- maxdist_path_CDmat ; ia=ia+1}
    if ("netstats_df" %in% outputs){netstats_df_list[[ib]]<-maxdist_netstats ; ib=ib+1}
    if ("igraph" %in% outputs){igraph_list[[ic]]<-maxdist_ig ; ic=ic+1}
    if ("igraph_wt" %in% outputs){igraph_wt_list[[id]]<-maxdist_ig_wt ; id=id+1}
    if ("lines_sf" %in% outputs){lines_sf_list[[ie]]<-maxdist_lines_sf ; ie=ie+1}
    if ("adjmat" %in% outputs){adjmat_list[[ih]]<-maxdist_am ; ih=ih+1}
    if ("adjmatCD" %in% outputs){adjmatCD_list[[ii]]<-maxdist_amCD ; ii=ii+1}
    if ("nb" %in% outputs){nb_list[[ij]]<-maxdist_nb ; ij=ij+1}
    
    if (ggplots == "MultiPlot"){
      gg_maxdist <- ggpltz(ggbasemap=ggbasemap, netstats_df=maxdist_netstats, 
                           CatchLims.sf=CatchLims.sf,poly=Catch.sf, s=c(0.75,1.5),
                           subtit=F,#subtitle = "Max Distance Net", 
                           lines_sf = maxdist_lines_sf, geomet="geometry")
      MultiPlotList[[j]] <- gg_maxdist +theme(legend.position="none",
                                              plot.margin = unit(c(0,0,0,0), "cm"))
      j=j+1
    }
    if (mode == "single" & ggplots == "PlotResults"){
      ggout <- ggpltz(ggbasemap=ggbasemap, netstats_df=maxdist_netstats, 
                      CatchLims.sf=CatchLims.sf,poly=Catch.sf, 
                      subtitle = "Max Distance Net", 
                      lines_sf = maxdist_lines_sf, geomet="geometry")
      ggout <- ggout  + labs(title=ggPlotTitle) + 
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
    }
  }
  
  if (mode == "multicombine"){
    comb_am_wt <- Reduce('+', adjmat_list)
    comb_am <- as.matrix((comb_am_wt>0)+0)
    comb_amCD <- CD.mat*comb_am
    comb_ig_dwt <- graph_from_adjacency_matrix(comb_amCD,mode = "undirected",
                                               weighted=T,add.colnames = NULL,
                                               add.rownames = NULL)
    comb_ig_wt <- graph_from_adjacency_matrix(comb_am_wt,mode = "undirected",
                                              weighted=T,add.colnames = NULL,
                                              add.rownames = NULL)
    comb_ig <- graph_from_adjacency_matrix(comb_amCD,mode = "undirected",
                                           add.colnames = NULL,
                                           add.rownames = NULL)
    comb_ig = simplify(comb_ig)
    x=as.data.frame(as_edgelist(comb_ig_wt))
    x$wt =edge_attr(comb_ig_wt)$weight
    edges <- as.data.frame(matrix(NA, nrow(edgelist), 4))
    colnames(edges) <- c("X1", "Y1", "X2", "Y2")
    for (i in seq_len(nrow(x))) {
      edges[i, ] <- c(nodes[which(nodes$AggSite == x[i, 1]), "East"],
                      nodes[which(nodes$AggSite == x[i, 1]), "North"],
                      nodes[which(nodes$AggSite == x[i, 2]), "East"],
                      nodes[which(nodes$AggSite == x[i, 2]), "North"])
    }
    xx=cbind(x,edges)
    comb_nb <- mat2listw(comb_am,style="B")
    comb_nb <- comb_nb$neighbours
    comb_nb_wt <- mat2listw(comb_am_wt)
    wt <- comb_nb_wt$weights
    comb_nb_wt <- mat2listw(comb_am_wt)
    comb_nb_wt <- comb_nb_wt$neighbours
    comb_lines_sf <- nb2lines(comb_nb, coords=coords, as_sf=T)
    st_crs(comb_lines_sf) <- st_crs(CatchLims.sf)
    comb_lines_wt_sf <- nb2lines(comb_nb_wt, coords=coords, as_sf=T)
    st_crs(comb_lines_wt_sf) <- st_crs(CatchLims.sf)
    comb_lines_wt_sf$wt <- unlist(wt)
    comb_path_CDmat <- distances(comb_ig_dwt,  mode = "all")
    comb_netstats = net_stats(comb_ig,o,"comb",coords)
    comb_netstats_wt = net_stats(comb_ig_wt,o,"comb",coords)
    
    if (combine == "unweighted"){
      if ("net_distmat" %in% outputs){net_distmat_list[[ia]] <- comb_path_CDmat ; ia=ia+1}
      if ("netstats_df" %in% outputs){netstats_df_list[[ib]]<-comb_netstats ; ib=ib+1}
      if ("igraph" %in% outputs){igraph_list[[ic]]<-comb_ig ; ic=ic+1}
      if ("igraph_wt" %in% outputs){igraph_wt_list[[id]]<-comb_ig_wt ; id=id+1}
      if ("lines_sf" %in% outputs){lines_sf_list[[ie]]<-comb_lines_sf ; ie=ie+1}
      if ("adjmat" %in% outputs){adjmat_list[[ih]]<-comb_am ; ih=ih+1}
      if ("adjmatCD" %in% outputs){adjmatCD_list[[ii]]<-comb_amCD ; ii=ii+1}
      if ("nb" %in% outputs){nb_list[[ij]]<-comb_nb ; ij=ij+1}
      if (ggplots == "PlotResults"){
        ggout <- ggpltz(ggbasemap=ggbasemap, netstats_df=comb_netstats, 
                        CatchLims.sf=CatchLims.sf,poly=Catch.sf, 
                        subtitle = "Combined Metrics Unweighted Net", 
                        lines_sf = comb_lines_sf, geomet="geometry")
        ggout <- ggout + labs(title=ggPlotTitle) + 
          theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
      }
      
    }
    if (combine == "weighted"){
      if ("net_distmat" %in% outputs){net_distmat_list[[ia]] <- comb_path_CDmat ; ia=ia+1}
      if ("netstats_df" %in% outputs){netstats_df_list[[ib]]<-comb_netstats_wt ; ib=ib+1}
      if ("igraph" %in% outputs){igraph_list[[ic]]<-comb_ig_wt ; ic=ic+1}
      if ("igraph_wt" %in% outputs){igraph_wt_list[[id]]<-comb_ig_dwt ; id=id+1}
      if ("lines_sf" %in% outputs){lines_sf_list[[ie]]<-comb_lines_wt_sf ; ie=ie+1}
      if ("adjmat" %in% outputs){adjmat_list[[ih]]<-comb_am_wt ; ih=ih+1}
      if ("adjmatCD" %in% outputs){adjmatCD_list[[ii]]<-comb_amCD ; ii=ii+1}
      if ("nb" %in% outputs){nb_list[[ij]]<-comb_nb_wt ; ij=ij+1}
      if (ggplots == "PlotResults"){
        ggout <- ggpltz(ggbasemap=ggbasemap, netstats_df=comb_netstats_wt, 
                        CatchLims.sf=CatchLims.sf,poly=Catch.sf, 
                        subtitle = "Combined Metrics Weighted Net", 
                        lines_sf = xx, geomet="geometry",wt=T)
        ggout <- ggout  + labs(title=ggPlotTitle) + 
          theme(plot.title = element_text(hjust = 0.5, face="bold", size=16))
      }
    }
    
    
  }
  
  #outputs list of lists
  l=1
  if ("net_distmat" %in% outputs){outlist[[l]] <- net_distmat_list ; l=l+1}
  if ("netstats_df" %in% outputs){outlist[[l]]<-netstats_df_list ; l=l+1}
  if ("igraph" %in% outputs){outlist[[l]]<-igraph_list ; l=l+1}
  if ("igraph_wt" %in% outputs){outlist[[l]]<-igraph_wt_list ; l=l+1}
  if ("lines_sf" %in% outputs){outlist[[l]]<-lines_sf_list ; l=l+1}
  if ("adjmat" %in% outputs){outlist[[l]]<-adjmat_list  ; l=l+1}
  if ("adjmatCD" %in% outputs){outlist[[l]]<-adjmatCD_list ; l=l+1}
  if ("nb" %in% outputs){outlist[[l]]<-nb_list ; l=l+1}
  
  
  if (ggplots == "MultiPlot"){
    
    ggx = length(MultiPlotList)
    
    if (ggx == 1){ nr = 1 ; nc = 1
    ggplotgrid = plot_grid(MultiPlotList[[1]],align="hv", nrow = nr, ncol = nc) }
    if (ggx == 2){ nr = 1 ; nc = 2
    ggplotgrid = plot_grid(MultiPlotList[[1]], MultiPlotList[[2]], align="hv", nrow = nr, ncol = nc)}
    if (ggx == 3){ nr = 2 ; nc = 2
    ggplotgrid = plot_grid(MultiPlotList[[1]], MultiPlotList[[2]], MultiPlotList[[3]],
                           align="hv", nrow = nr, ncol = nc)}
    if (ggx == 4){ nr = 2 ; nc = 2
    ggplotgrid = plot_grid(MultiPlotList[[1]], MultiPlotList[[2]], MultiPlotList[[3]],
                           MultiPlotList[[4]],align="hv", nrow = nr, ncol = nc)}
    if (ggx == 5){ nr = 2 ; nc = 3
    ggplotgrid = plot_grid(MultiPlotList[[1]], MultiPlotList[[2]], MultiPlotList[[3]],
                           MultiPlotList[[4]], MultiPlotList[[5]],align="hv", nrow = nr, ncol = nc)}
    if (ggx == 6){ nr = 2 ; nc = 3
    ggplotgrid = plot_grid(MultiPlotList[[1]], MultiPlotList[[2]], MultiPlotList[[3]],
                           MultiPlotList[[4]], MultiPlotList[[5]], MultiPlotList[[6]],align="hv", nrow = nr, ncol = nc)}
    if (ggx == 7){nr = 3 ; nc = 3
    lllabz = c("Gabriel", "RNN", "Deldir", "SOI", "Catch", "KNN", "MaxDist")
    ggplotgrid = plot_grid(MultiPlotList[[1]], MultiPlotList[[2]], MultiPlotList[[3]],
                           MultiPlotList[[4]], MultiPlotList[[5]], MultiPlotList[[6]], MultiPlotList[[7]],
                           labels = lllabz[1:ggx],align="hv", nrow = nr, ncol = nc)}
    
    title <- ggdraw() + draw_label(ggPlotTitle,fontface = 'bold',hjust = 0.5, size=20 )
    
    ggout <- plot_grid(title, ggplotgrid,ncol = 1,rel_heights = c(0.05,0.95)) + 
      theme(plot.background = element_rect(fill = "white", colour = NA))
  }
  x<-length(outlist)+1
  outlist[[x]]<-ggout
  return(outlist)
  
}








