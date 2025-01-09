#### SpInt_EstFlows.R
#### Rudolf Cesaretti, 6/31/2022

#### "SpInt_EstFlows" 
#### 
#### rihll_wilson
#### 
#### https://book.archnetworks.net/spatialinteraction

pak <- c("rgdal", "sp", "sf", "GISTools", "lwgeom", "tidyverse", "tidyr", "data.table", "zoo")
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
library(scales)
library(pracma)
library(modelsummary)
library(cowplot)


###############################################################
##########################  SpInt_EstFlows  #########################
###############################################################









# Set model parameters and initial variable states
Oi <- rep(1, nrow(dat))
Wj <- rep(1, nrow(dat))
alpha <- 1.05
beta <- 0.1
eps <- 0.01
K <- 1

# Define distance among points in kilometers. Because our
# points are in geographic coordinates we use the distm
# function. See the section on Spatial Networks for more.
library(geosphere)
d <- as.matrix(distm(dat[, c(7, 6)])) / 1000

# Dj is ineitial set as a vector of 1s like Wj
Dj <- Wj

# Create objects for keeping track of the number
# of iterations and the conditions required to end
# the loop.
end_condition <- 1
iter <- 0

# Define the deterrence function as a exponential
det <- exp(-beta * d)

# Create while loop that will continue to iterate Tij
# until 10,000 iterations or until the end_condition
# object is less than the threshold indicated.
while (!(end_condition < 1e-5) & iter < 10000) {
  # Set Wj to Dj
  Wj <- Dj
  # Calculate Tij as indicated above
  Tij <-
    apply(det * Oi %o% Wj ^ alpha, 2, '/',
          (Wj ^ alpha %*% det))
  # Calculate change in W using equation above
  delta_W <- eps * (colSums(Tij) - (K * Wj))
  
  # Calculate new Dj
  Dj <- delta_W + Wj
  
  # Add to iterator and check for end conditions
  iter <- iter  + 1
  end_condition <- sum((Dj - Wj) ^ 2)
}

hist(Wj, breaks = 15)








grav_mod <- function(attract, B, d) {
  res <- matrix(0, length(attract), length(attract))
  
  for (i in seq_len(length(attract))) {
    for (j in seq_len(length(attract))) {
      res[i, j] <-
        attract[i] * attract[j] * exp(-B * d[i,j])
    }
  }
  diag(res) <- 0
  return(res)
}





###############################
###### GRAVITY FUNCTIONS ######



#' Calculate Tij
#' 
#' Calculate Tij
#' 
#' @param Oi Vector of the ouputs of i
#' @param Wj Vector of the inputs of j
#' @param fcij matrix giving the deterrence function
#' @param alpha scaling factor for attractivity
#' 
#' @author Martin Hinz <martin.hinz@iaw.unibe.ch>
calculate_tij <- function(Oi, Wj, fcij, alpha) {
  
  Wj <- Wj^alpha
  fcij <- as.matrix(fcij)
  
  Tij <- apply((fcij * Oi %o% Wj), 2, `/`, t(Wj %*% t(fcij)))
  
  rownames(Tij) <- rownames(fcij)
  
  return(Tij)
}



#' Calculates the Rihll-Wilson-Model
#' 
#' Calculates the Rihll-Wilson-Model
#' 
#' @param Oi Outflows originating from i
#' @param Wj Attractiveness
#' @param fcij Deterrence function (has to be calculated from the distance
#'        cost matrix "cij" first)
#' @param alpha The scaling factor of attractivity#' 
#' @param eps .
#' @param maxrun maximum number of iterations
#' @return a list with the elements:
#' - a vector containing Ai
#' - a vector containing Ij
#' 
#' @export rihll_wilson
rihll_wilson <- function(Oi, Wj, fcij, alpha, eps = 1e-6, maxrun = 1000){
  
  iter <- 0
  
  Dj <- Wj
  
  epsilon <- 1
  
  while(!(epsilon<eps)&!(iter>maxrun)){
    
    Wj <- Dj
    
    Tij <- calculate_tij(Oi,Wj,fcij, alpha)
    
    Dj <- colSums(Tij, na.rm = TRUE)
    
    iter <- iter+1
    epsilon <- sum((Dj - Wj)^2)
  }
  
  names(Dj) <- rownames(fcij)
  rval <- list(Inputs = Dj,
               Tij = Tij,
               nb.iter = iter)
  return(rval)
}

####################################################################################

#' Simple Gravity Model
#' 
#' Simple Gravity Model
#' Mass = ability of entities to spread flows = ex : Mass or Population
#' 
#' @param Mi vector of mass for every site
#' @param fcij deterrence function from the distance matrix (has to be produced from a distance matrix before)
#' @param k adjustment variable
#' 
#' @return A matrix with the modeled flows from i to j
#' 
#' @author Clara Filet <clara.filet@gmail.com>
#' 
#' @export simple_gravity
simple_gravity <- function(Mi, fcij, k=1) {
  Mmatrix <- Mi%o%Mi
  Tij <- k*((Mmatrix)/fcij) 
  return(Tij)
}

##################################
###### DETERRENCE FUNCTIONS ######

#' Calculates the deterrence function
#' 
#' Calculates the deterrence function used for the different gravity based models
#' 
#' @param cij a matrix containing the cost matrix
#' @param beta the distance decay factor
#' @param type the family of deterrence function to use
#' @param alpha shape parameter 1 for ariadne type deterrence function
#' @param gamma shape parameter 2 for ariadne type deterrence function
#' 
#' @details Power is generally used for simple gravity models, exponential and ariadne can be used for the Rihll-Wilson model.
#' \code{alpha} and \code{gamma} are only relevant for the ariadne type model and default to the values alpha=4 and gamma=1 like in Evans&Rivers 2017.
#' 
#' 
#' @references Evans Tim S., Rivers Ray J., Was Thebes Necessary? Contingency in Spatial Modeling. Frontiers in Digital Humanities 4, 2017, \url{https://doi.org/10.3389/fdigh.2017.00008}
#' 
#' @return a matrix containing the deterrence function
#' 
#' @author Daniel Knitter <knitter@geographie.uni-kiel.de>
#' @author Martin Hinz <martin.hinz@iaw.unibe.ch>
#' @author Benedikt Grammer	<benedikt.grammer@univie.ac.at>
#' @author Kai Radloff <kai.radloff@hu-berlin.de>
#' @author Loren V. Cowin	<Lcowin@gshdl.uni-kiel.de>
#' @author Clara Filet <clara.filet@gmail.com>
#' 
#' @export deterrence_function
deterrence_function <- function(cij,beta,type = 'exponential', alpha=4, gamma=1) {
  if (type=='exponential'){
    fcij <- exp(-beta * cij)  
  }
  else if (type == 'negpower'){
    fcij <- cij^-beta 
  }
  else if (type == 'power'){
    fcij <- cij^beta 
  }
  else if (type == 'ariadne'){
    fcij <- (1 + (cij^beta)^alpha)^-gamma
  }
  else {
    stop ("Sorry, I do not know this kind of deterrence function")
  }
  return(fcij)
}