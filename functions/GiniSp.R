#### GiniSp.R
#### Rudolf Cesaretti, 7/4/2022

#### "LQ" = Location Quotient
#### 
#### 
#### 
#### 

require(tidyverse)
require(ineq)

###############################################################
########################    GiniSp.knn    #####################
###############################################################


GiniSp.knn <- function(df, var, m, k){

  Global <- ineq::Gini(as.numeric(unlist(df[,var])))
  
  if (is.nan(Global)){out <- NA}else{
    
    x <- as.data.frame(m)
    
    out <- x %>% rownames_to_column(var = "Site1") %>% 
      tidyr::pivot_longer(!Site1,names_to="Site2",values_to ="CDist") %>% 
      #filter(CDist >= r)) %>% 
      group_by(Site1) %>% arrange(CDist, group_by=T) %>% slice(1:(k+1)) %>%
      left_join(select(df, AggSite, !!sym(var)), by = c("Site2" = "AggSite")) %>% 
      summarize(spGini = ineq::Gini(!!sym(var))/Global) %>% 
      arrange(match(Site1, df$AggSite)) %>% pull(spGini)
  }
    
  return(out)
  
}

###############################################################
######################    GiniSp.dist    ######################
###############################################################

GiniSp.dist <- function(df, var, m, r){

  Global <- ineq::Gini(as.numeric(unlist(df[,var])))
  
  if (is.nan(Global)){out <- NA}else{
    
    x <- as.data.frame(m)
    
    out <- x %>% rownames_to_column(var = "Site1") %>% 
      tidyr::pivot_longer(!Site1,names_to="Site2",values_to ="CDist") %>% 
      filter(CDist <= r) %>% group_by(Site1) %>% 
      #arrange(CDist, group_by=T) %>% slice(1:(k+1)) %>%
      left_join(select(df, AggSite, !!sym(var)), by = c("Site2" = "AggSite")) %>% 
      summarize(spGini = ineq::Gini(!!sym(var))/Global) %>% 
      arrange(match(Site1, df$AggSite)) %>% pull(spGini)
    
  }
  
  return(out)
  
}

