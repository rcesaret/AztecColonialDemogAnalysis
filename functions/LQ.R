#### LQ.R
#### Rudolf Cesaretti, 7/4/2022

#### "LQ" = Location Quotient
#### 
#### 
#### 
#### 

require(tidyverse)

###############################################################
########################       LQ       #######################
###############################################################

LQ <- function(df, e, E){
  
  Denom <- df %>% summarize(Denom = sum(!!sym(e), na.rm=T)/sum(!!sym(E), na.rm=T)) %>% pull(Denom)
  
  if (Denom == 0){out <- NA}else{
    
  out <- df %>% summarize(LQ = (!!sym(e)/!!sym(E))/Denom) %>% pull(LQ)
  
  }
  
  return(out)
  
}

###############################################################
########################    LQ.knn    #########################
###############################################################

FLQ.knn <- function(df, e, E, m, k){

  Denom <- df %>% summarize(Denom = sum(!!sym(e), na.rm=T)/sum(!!sym(E), na.rm=T)) %>% pull(Denom)
  
  if (Denom == 0){out <- NA}else{
    
    x <- as.data.frame(m)
    
    out <- x %>% rownames_to_column(var = "Site1") %>% 
      tidyr::pivot_longer(!Site1,names_to="Site2",values_to ="CDist") %>% 
      group_by(Site1) %>% arrange(CDist, group_by=T) %>% slice(1:(k+1)) %>%
      left_join(select(df, AggSite, !!sym(e), !!sym(E)), by = c("Site2" = "AggSite")) %>% 
      summarize(LQ = (sum(!!sym(e), na.rm=T)/sum(!!sym(E), na.rm=T)/Denom)) %>% 
      mutate(LQ=ifelse(is.nan(LQ),NA,LQ)) %>% 
      arrange(match(Site1, df$AggSite)) %>% pull(LQ)
    
  }
  
  return(out)
  
}

###############################################################
########################    LQ.dist    ########################
###############################################################

FLQ.dist <- function(df, e, E, m, r){
  
  Denom <- df %>% summarize(Denom = sum(!!sym(e), na.rm=T)/sum(!!sym(E), na.rm=T)) %>% pull(Denom)
  
  if (Denom == 0){out <- NA}else{
    
    x <- as.data.frame(m)
    
    out <- x %>% rownames_to_column(var = "Site1") %>% 
      tidyr::pivot_longer(!Site1,names_to="Site2",values_to ="CDist") %>% 
      filter(CDist <= r) %>% group_by(Site1) %>% 
      left_join(select(df, AggSite, !!sym(e), !!sym(E)), by = c("Site2" = "AggSite")) %>% 
      summarize(LQ = (sum(!!sym(e), na.rm=T)/sum(!!sym(E), na.rm=T)/Denom)) %>% 
      mutate(LQ=ifelse(is.nan(LQ),NA,LQ)) %>% 
      arrange(match(Site1, df$AggSite)) %>% pull(LQ)
    
  }
  
  return(out)
  
}
