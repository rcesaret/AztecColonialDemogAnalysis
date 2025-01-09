#### igraph2nb.R
#### Rudolf Cesaretti, 7/6/2022

#### "igraph2nb" converts an igraph class object to spdep neighbor class object

pak <- c("igraph", "ade4", "spdep")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)

###############################################################
#########################  igraph2nb  #########################
###############################################################

igraph2nb <- function(ig) {
  return(neig2nb(neig(edges=get.edgelist(ig))))
}
