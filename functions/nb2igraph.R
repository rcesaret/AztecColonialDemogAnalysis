#### nb2igraph.R
#### Rudolf Cesaretti, 7/6/2022

#### "igraph2nb" converts an igraph class object to spdep neighbor class object

pak <- c("igraph", "spdep")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)

###############################################################
#########################  nb2igraph  #########################
###############################################################

nb2igraph <- function(nb, stylenb="B", weights=TRUE, nbwts=NULL, matwts = matrix(1, length(nb), length(nb)), type="undirected") {
  adjmat=nb2mat(nb, style=stylenb, glist=nbwts)
  adjmat=adjmat*matwts
  igr = graph_from_adjacency_matrix(adjmat,mode = type,weighted = weights,diag = T,add.colnames = NA,add.rownames = NULL)
  return(igr)
}
