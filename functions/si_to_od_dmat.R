#### si_to_od_dmat.R
#### Rudolf Cesaretti, 6/31/2022

#### "si_to_od_dmat" 
#### 
#### 
#### 

require(sf)
require(od)
require(geodist)
require(simodels)
require(tidyverse)


#origins=st_as_sf(Pts_List[[1]])
#destinations=st_as_sf(Pts_List[[1]])
#dmat = CD.mats[[1]]
#out_crs = st_crs(origins)
#max_dist = Inf
#intrazonal = TRUE

si_to_od_dmat = function(origins, 
                         destinations, 
                         dmat,
                         out_crs = st_crs(origins),
                         max_dist = Inf, 
                         intrazonal = TRUE 
                    ) {
  st_crs(origins)
  if(identical(origins, destinations)) {
    od_df = od::points_to_od(origins)
  } else {
    od_df = od::points_to_od(origins, destinations)
  }
  
  dmat <- as.data.frame(dmat)
  dmat <- dmat %>% rownames_to_column(var = "O") %>% 
    tidyr::pivot_longer(!O,names_to="D",values_to ="distance")
  
  od_df <- od_df %>% left_join(dmat, by = c("O","D"))
    
  # Max dist check
  if(max(od_df$distance) > max_dist) {
    nrow_before = nrow(od_df)
    od_df = od_df[od_df$distance <= max_dist, ]
    nrow_after = nrow(od_df)
    pct_kept = round(nrow_after / nrow_before * 100)
    message(
      nrow_after,
      " OD pairs remaining after removing those with a distance greater than ", # nolint
      max_dist, " meters", ":\n",
      pct_kept, "% of all possible OD pairs"
    )
  }
  
  # Intrazonal check
  if(!intrazonal){
    od_df = dplyr::filter(od_df, distance > 0)
  }
  
  od_sfc = od::odc_to_sfc(od_df[3:6])
  sf::st_crs(od_sfc) = out_crs 
  od_df = od_df[-c(3:6)]
  
  # join origin attributes
  origins_to_join = sf::st_drop_geometry(origins)
  names(origins_to_join) = paste0("origin_", names(origins_to_join))
  names(origins_to_join)[1] = names(od_df)[1]
  od_dfj = dplyr::inner_join(od_df, origins_to_join, by = "O")
  # join destination attributes
  destinations_to_join = sf::st_drop_geometry(destinations)
  names(destinations_to_join) = paste0("destination_", names(destinations_to_join)) # nolint
  names(destinations_to_join)[1] = names(od_df)[2]
  od_dfj = dplyr::inner_join(od_dfj, destinations_to_join, by = "D")
  # names(od_dfj)
  # create and return sf object
  sf::st_sf(od_dfj, geometry = od_sfc)
}

