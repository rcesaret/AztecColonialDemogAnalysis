#net_stats2.R
#Rudolf Cesaretti
#7/6/2022

require(igraph)
require(scales)
require(centiserve)
require(CINNA)


#ig=Tij_ig
#ID=rownames(d_mat2)
#prefix="Rad_"
#direc=T
#weights=T


net_stats2 <- function(ig,ID,prefix="",direc=T, weights=T){
  out <- data.frame(ID = ID)
  rownames(out) <- out$ID
  if (weights==T){
    if (direc==T){
      igraph::graph.strength(ig, mode = "all")
      out$centr_indeg <- igraph::graph.strength(ig, mode = "in")
      out$centr_outdeg <- igraph::graph.strength(ig, mode = "out")
    }
    out$centr_totdeg <- igraph::graph.strength(ig, mode = "all")
  }
  if (weights==F){
    if (direc==T){
      out$centr_indeg <- igraph::degree(ig, mode = "in")
      out$centr_outdeg <- igraph::degree(ig, mode = "out")
    }
    out$centr_totdeg <- igraph::degree(ig, mode = "total")
  }
  if (weights==F){
    out$centr_botnk <- centiserve::bottleneck(ig, mode = "all")
    out$centr_clstrk <- centiserve::clusterrank(ig, directed = direc)
    out$centr_difus <- centiserve::diffusion.degree(ig, mode = "all")
    out$centr_dmnc <- centiserve::dmnc(ig, mode = "all")
    out$centr_semiloc <- centiserve::semilocal(ig, mode = "all")
    out$centr_entropy <- centiserve::entropy(ig, mode = "all")
  }
  if (direc==T){
    igraph::graph.strength(ig, mode = "all")
    out$centr_baryin <- centiserve::barycenter(ig, mode = "in")
    out$centr_baryout <- centiserve::barycenter(ig, mode = "out")
  }
  out$centr_barytot <- centiserve::barycenter(ig, mode = "all")
  if (direc==T){
    igraph::graph.strength(ig, mode = "all")
    out$centr_linin <- centiserve::lincent(ig, mode = "in")
    out$centr_linout <- centiserve::lincent(ig, mode = "out")
  }
  out$centr_lintot <- centiserve::lincent(ig, mode = "all")
  out$centr_btw <- igraph::betweenness(ig)
  out$centr_eig <- igraph::eigen_centrality(ig, directed=direc, scale = F)$vector
  if (direc==T){
    out$centr_inclos <- igraph::closeness(ig, mode = "in")
    out$centr_outclos <- igraph::closeness(ig, mode = "out")
  }
  out$centr_totclos <- igraph::closeness(ig, mode = "total")
  out$centr_hrmo <- igraph::harmonic_centrality(ig)
  out$centr_hub <- igraph::hub_score(ig, scale = F)$vector
  out$centr_auth <- igraph::authority_score(ig, scale = F)$vector
  out$centr_Kleinauth <- as.numeric(unlist(CINNA::calculate_centralities(ig, include = "Kleinberg's authority centrality scores")))
  out$centr_pgrk <- igraph::page_rank(ig)$vector
  out$centr_alph <- igraph::alpha_centrality(ig)
  
  
  out2 <- out[,-1]
  for (i in ncol(out2)){
    out2[,i] <- scales::rescale(out2[,i])
  }
  coln <- colnames(out2)
  out$centr_rs_avg <- out2 %>% rowwise() %>% mutate(
      centr_rs_avg = mean(!!!syms(coln), na.rm=T)) %>% ungroup() %>% pull(centr_rs_avg)
  
  if (weights==F){
    out$trans_loc <- igraph::transitivity(ig, type = "localundirected", isolates = "zero")
    out$trans_locavg <- igraph::transitivity(ig, type = "localaverageundirected", isolates = "zero")
    out$density <- igraph::edge_density(ig)
    ut$connectiv <- igraph::edge_connectivity(ig)
    out$trans_glob <- igraph::transitivity(ig, type = "globalundirected", isolates = "zero")
    if (direc==T){
      out$incoreness <- igraph::coreness(ig, mode="in")
      out$outcoreness <- igraph::coreness(ig, mode="out")
    }
    out$totcoreness <- igraph::coreness(ig, mode="all")
  }
    if (direc==F){out$diversity <- igraph::diversity(ig)}
  out$constraint <- igraph::constraint(ig)
  if (direc==T){
    out$avgdist_in <- centiserve::averagedis(ig, mode = "in")
    out$avgdist_out <- centiserve::averagedis(ig, mode = "out")
  }
  out$avgdist_tot <- centiserve::averagedis(ig, mode = "all")
  
  
  if (weights==F){
    out3 <- out2[,1]
    if (direc==T){
      out3$cntrlz_indeg <- out$cntrlz_indeg <- igraph::centr_degree(ig, mode="in")$centralization
      out3$cntrlz_outdeg <- out$cntrlz_outdeg <- igraph::centr_degree(ig, mode="out")$centralization
    }
    out3$cntrlz_totdeg <- out$cntrlz_totdeg <- igraph::centr_degree(ig, mode="total")$centralization
    out3$cntrlz_btw <- out$cntrlz_btw <- igraph::centr_betw(ig, directed=direc)$centralization
    out3$cntrlz_eig <- out$cntrlz_eig <- igraph::centr_eigen(ig, directed=direc, scale = F)$centralization
    if (direc==T){
      out3$cntrlz_inclo <- out$cntrlz_inclo <- igraph::centr_clo(ig, mode="in")$centralization
      out3$cntrlz_outclo <- out$cntrlz_outclo <- igraph::centr_clo(ig, mode="out")$centralization
    }
    out3$cntrlz_totclo <- out$cntrlz_totclo <- igraph::centr_clo(ig, mode="total")$centralization
    out3$cntrlz_hub <- out$cntrlz_hub <- igraph::hub_score(ig, scale = F)$value
    out3$cntrlz_auth <- out$cntrlz_auth <- igraph::authority_score(ig, scale = F)$value
  
    for (i in ncol(out3)){
      out3[,i] <- scales::rescale(out3[,i])
    }
    coln <- colnames(out3)
    out$centr_rs_avg <- out3 %>% rowwise() %>% mutate(
      centr_rs_avg = mean(!!!syms(coln), na.rm=T)) %>% ungroup() %>% pull(centr_rs_avg)
  }
  
  if (direc==F){
    out$fstgrd_com <- igraph::membership(cluster_fast_greedy(ig))
    out$fstgrd_mod <- igraph::modularity(ig,out$fstgrd_com)
    out$leiden_com <- igraph::membership(cluster_leiden(ig))
    out$leiden_mod <- igraph::modularity(ig, out$leiden_com)
    out$louvain_com <- igraph::membership(cluster_louvain(ig))
    out$louvain_mod <- igraph::modularity(ig, out$louvain_com)
  }
  
  #out$edgbtw_com <- igraph::membership(cluster_edge_betweenness(ig, directed=direc))
  #out$edgbtw_mod <- igraph::modularity(ig, out$edgbtw_com, directed=direc)
  #out$wlktrp_com <- membership(cluster_walktrap(ig))
  #out$wlktrp_mod <- igraph::modularity(ig, out$wlktrp_com, directed=direc)
  
  cn <- paste0(prefix,colnames(out))
  colnames(out) <- c("ID",cn[2:length(cn)])
  
  return(out)
}
