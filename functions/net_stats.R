#net_stats.R
#Rudolf Cesaretti
#7/6/2022

require(igraph)

net_stats <- function(ig,ID,graph,coords){
  out <- data.frame(ID = ID)
  rownames(out) <- out$ID
  out$East <- as.numeric(coords[,1])
  out$North <- as.numeric(coords[,2])
  out$graph <- graph
  out$centr_deg <- igraph::degree(ig)
  out$centr_btw <- igraph::betweenness(ig)
  out$centr_eig <- igraph::eigen_centrality(ig, scale = F)$vector
  out$centr_clos <- igraph::closeness(ig)
  out$centr_hrmo <- igraph::harmonic_centrality(ig)
  out$centr_hub <- igraph::hub_score(ig, scale = F)$vector
  #out$centr_pwr <- igraph::power_centrality(ig, rescale = F)
  out$centr_auth <- igraph::authority_score(ig, scale = F)$vector
  out$centr_pgrk <- igraph::page_rank(ig)$vector
  #out$centr_alph <- igraph::alpha_centrality(ig)
  
  out$centr_deg_n <- igraph::degree(ig, normalized = T)
  out$centr_btw_n <- igraph::betweenness(ig, normalized = T)
  out$centr_eig_n <- igraph::eigen_centrality(ig, scale = T)$vector
  out$centr_clos_n <- igraph::closeness(ig, normalized = T)
  out$centr_hrmo_n <- igraph::harmonic_centrality(ig, normalized = T)
  out$centr_hub_n <- igraph::hub_score(ig, scale = T)$vector
  #out$centr_pwr_n <- igraph::power_centrality(ig, rescale = T)
  out$centr_auth_n <- igraph::authority_score(ig, scale = T)$vector
  out <- out %>% rowwise() %>% mutate(
    centr_avg = mean(c(centr_deg_n,centr_btw_n,centr_eig_n, centr_clos_n, 
                       centr_hrmo_n, centr_hub_n, centr_auth_n), 
                     na.rm=T)) %>% ungroup()
  
  out$trans_loc <- igraph::transitivity(ig, type = "localundirected", isolates = "zero")
  out$trans_locavg <- igraph::transitivity(ig, type = "localaverageundirected", isolates = "zero")
  
  out$density <- igraph::edge_density(ig)
  out$connectiv <- igraph::edge_connectivity(ig)
  out$trans_glob <- igraph::transitivity(ig, type = "globalundirected", isolates = "zero")
  
  out$cntrlz_deg <- igraph::centr_degree(ig)$centralization
  out$cntrlz_btw <- igraph::centr_betw(ig)$centralization
  out$cntrlz_eig <- igraph::centr_eigen(ig, scale = F)$centralization
  out$cntrlz_clo <- igraph::centr_clo(ig)$centralization
  out$cntrlz_hub <- igraph::hub_score(ig, scale = F)$value
  out$cntrlz_auth <- igraph::authority_score(ig, scale = F)$value
  
  out$cntrlz_deg_n <- igraph::centr_degree(ig, normalized = T)$centralization
  out$cntrlz_btw_n <- igraph::centr_betw(ig, normalized = T)$centralization
  out$cntrlz_eig_n <- igraph::centr_eigen(ig, scale = T)$centralization
  out$cntrlz_clo_n <- igraph::centr_clo(ig, normalized = T)$centralization
  out$cntrlz_hub_n <- igraph::hub_score(ig, scale = T)$value
  out$cntrlz_auth_n <- igraph::authority_score(ig, scale = T)$value
  out <- out %>% rowwise() %>% mutate(
    cntrlz_avg = mean(c(cntrlz_deg_n,cntrlz_btw_n,cntrlz_eig_n, cntrlz_clo_n, 
                        cntrlz_hub_n, cntrlz_auth_n), 
                      na.rm=T)) %>% ungroup()
  
  #out$fstgrd_mod <- igraph::modularity(ig, membership(cluster_fast_greedy(ig)))
  #out$leiden_mod <- igraph::modularity(ig, membership(cluster_leiden(ig)))
  #out$louvain_mod <- igraph::modularity(ig, membership(cluster_louvain(ig)))
  #out$edgbtw_mod <- igraph::modularity(ig, membership(cluster_edge_betweenness(ig)))
  #out$wlktrp_mod <- igraph::modularity(ig, membership(cluster_walktrap(ig)))
  return(out)
}
