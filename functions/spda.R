#### spda.R
#### Rudolf Cesaretti, 6/1/2022


#### "spda", which stands for "Settlement Survey Probability Density Analysis",
#### is adapted from Matt Peeples' "updf" function implementation of Ortman (2016).
#### The function gives a number of options for calculating occupational probabilities 
#### from artifact counts and artifact type ranges (see adjoining .Rmd for details)


pak <- c("compiler", "msm", "snowfall", "parallel", "doParallel", "tidyverse", "era", "zoo", "scales", "data.table")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)

# Multicore setup
cl <- makeCluster(detectCores())
registerDoParallel(cl)

#############################################################################
#################################  spda  ####################################
#############################################################################


spda <- function(
    site,        # vector of site names with one entry for every site/type/count    
    ID,          # vector of site IDs with one entry for every site/type/count  
    ph.sites,    # vector of site names for overlap and phase sites
    ph.periods,  # vector of period names for overlap and ceramic phase periods
    cer.type,    # vector of type names with one entry for every site/type/count
    ct,          # vector of counts with one entry for every site/type/count
    start,       # vector of type beginning dates with one entry for every site/type/count
    start.era,   # vector of eras corresponding to each of the end dates (e.g. BC, AD, CE, BCE)
    end,         # vector of type ending dates with one entry for every site/type/count
    end.era,     # vector of eras corresponding to each of the end dates (e.g. BC, AD, CE, BCE)
    pop.ests,    # vector of population estimates for every site/type/count
    pc.input = NULL,
    obs.input = NULL,
    chron=rep_len(1,length.out=length(cer.type)),
    interval=1,
    cutoff=0.05,
    min.period=25, 
    pc.method = c("input", "uniform", "tnormal", "beta"), 
    method = c("bayesian","mean.obs"),
    z = rep_len(2,length.out=length(cer.type)),
    alpha = rep_len(2,length.out=length(cer.type)),
    beta = rep_len(2,length.out=length(cer.type)),
    UrbanThresh = 1500){
  
  site2 <- site
  
  ## Remove any type with a count of 0 from all vectors
  trim <- which(ct>0)
  cer.type <- as.vector(cer.type[trim])
  site <- as.vector(site[1])
  start <- start[trim]
  start.era <- start.era[trim]
  end <- end[trim]
  end.era <- end.era[trim]
  ct <- ct[trim]
  chron <- chron[trim]
  alpha <- alpha[trim]
  beta <- beta[trim]
  z <- z[trim]
  
  
  ## define function for rounding by base defined in "interval" argument and round start and end dates for each type
  mround <- function(x,base){base*round(x/base)}
  start <- mround(start,interval)
  end <- mround(end,interval)
  
  ## convert temporal scale to "Holocene Era (HE)" so that all years are positive
  start = abs(start)
  end = abs(end)
  start.lab = start
  end.lab = end
  for (i in 1:length(start)){
    start[i] = as.numeric(yr_transform(yr(start[i], start.era[i]), "HE"))
    end[i] = as.numeric(yr_transform(yr(end[i], end.era[i]), "HE"))
  }
  
  ## find minimal ceramic intervals based on rounded type date overlaps
  df = data.frame(years = c(start,end), years.lab = c(start.lab,end.lab), era.lab = c(start.era,end.era))
  df = df[!duplicated(df), ]
  df <- df[order(df$years),]
  years = df$years
  years.lab = df$years.lab
  years.era.lab = df$era.lab
  period <- do.call(rbind,foreach(i=1:(length(years)-1)) %dopar% rbind(c(years[i],years[i+1])))
  period.lab <- do.call(rbind,foreach(i=1:(length(years.lab)-1)) %dopar% rbind(c(years.lab[i],years.era.lab[i],years.lab[i+1],years.era.lab[i+1])))
  
  # Ensure that no period is shorter than the specificed min.period argument and combine periods to satisfy this requirement
  if (length(period)>2) {
    per.short <- foreach(i=1:(nrow(period)-1)) %dopar% (period[i+1,1]-period[i,1])
    comb.per <- which(per.short<min.period)
    if (length(comb.per)>0){ 
      period <- period[-comb.per,]
      period.lab <- period.lab[-comb.per,]
      for (i in 1:(nrow(period)-1)) {
        period[i,2] <- period[i+1,1]
        period.lab[i,3] <- period.lab[i+1,1]
        period.lab[i,4] <- period.lab[i+1,2]}}
  }
  
  ## period (column) names for intervals in AD/BC/CE/BCE
  
  l1a <- paste(period.lab[which(period.lab[,2] == "BC" | period.lab[,2]=="CE" | period.lab[,2]=="BCE"),1],period.lab[which(period.lab[,2] == "BC" | period.lab[,2]=="CE" | period.lab[,2]=="BCE"),2],sep="")
  l1b <- paste(period.lab[which(period.lab[,2] == "AD"),2],period.lab[which(period.lab[,2] == "AD"),1],sep="")
  l2a <- paste(period.lab[which(period.lab[,4] == "BC" | period.lab[,4]=="CE" | period.lab[,4]=="BCE"),3],period.lab[which(period.lab[,4] == "BC" | period.lab[,4]=="CE" | period.lab[,4]=="BCE"),4],sep="")
  l2b <- paste(period.lab[which(period.lab[,4] == "AD"),4],period.lab[which(period.lab[,4] == "AD"),3],sep="")
  col.nam <- paste0(c(l1a,l1b),"-",c(l2a,l2b))
  
  
  ## Define lists of ceramic period lengths and site period lengths
  cer.per <- foreach(i=1:length(start)) %dopar% (seq((start[i]:end[i]))+(start[i]-1))
  per.len <- foreach(i=1:nrow(period)) %dopar% (seq(period[i,1]:period[i,2])+(period[i,1])-1)
  
  plen = unlist(lapply(per.len, FUN=function(x){(length(x)-1)}))
  pl = as.data.frame(period.lab)
  pl[,1] <- as.numeric(pl[,1])
  pl[,3] <- as.numeric(pl[,3])
  colnames(pl) = c("PeriodBegin", "PeriodBegin.era", "PeriodEnd", "PeriodEnd.era")
  pl <- pl %>% rowwise() %>% mutate(
    PeriodBegin = ifelse(PeriodBegin.era == "BC"  | PeriodBegin.era == "BCE", PeriodBegin * -1, PeriodBegin),
    PeriodEnd = ifelse(PeriodEnd.era == "BC"  | PeriodEnd.era == "BCE", PeriodEnd * -1, PeriodEnd))%>% ungroup() %>% mutate(
      MidLength = plen/2)
  pl <- pl %>% mutate(PeriodMidpoint = PeriodBegin + MidLength) %>% 
    rowwise() %>% mutate(
      PeriodMidpoint.era = ifelse(PeriodMidpoint < 0, "BC", "AD")) %>% ungroup() %>% mutate(
        AbsPeriodMidpoint = abs(PeriodMidpoint))
  
  ### User Input Popularity Curve Method
  if (pc.method == "input"){
    #reduce popularity curve to the rows and columns evident in the data
    if (length(cer.type) == 1){
      pc.input = pc.input[rownames(pc.input) %in% cer.type, ]
      pc.input <- as.matrix(sum(pc.input))
      rownames(pc.input) <- cer.type
      colnames(pc.input) <- col.nam
    } else {
      pc.input = pc.input[rownames(pc.input) %in% cer.type, ]
      pc.input = pc.input[,colSums(pc.input)>0 ]
      merge.cols = setdiff(colnames(pc.input),col.nam)
      if (length(merge.cols) > 0){
        if (length(merge.cols) == 2){
          merge.cols.nums = which( colnames(pc.input)==merge.cols)
          if (merge.cols.nums[1] == 1){
            pc.input[,1] = pc.input[,1] + pc.input[,2]
            pc.input = pc.input[,-2]
          } else {
            pc.input[,ncol(pc.input)] = pc.input[,ncol(pc.input)] + pc.input[,(ncol(pc.input)-1)]
            pc.input = pc.input[,-(ncol(pc.input)-1)]
          }
        } else {
          pc.input[,1] = pc.input[,1] + pc.input[,2]
          pc.input = pc.input[,-2]
          pc.input[,ncol(pc.input)] = pc.input[,ncol(pc.input)] + pc.input[,(ncol(pc.input)-1)]
          pc.input = pc.input[,-(ncol(pc.input)-1)]
        }
      }
    }
    colnames(pc.input) <- col.nam
    prob <- pc.input
    if (nrow(prob)==1) {
      prob <- as.matrix(prob)
      colnames(prob) <- col.nam}
    
    
  }
  
  ### Beta Popularity Curve Method   
  if (pc.method == "beta"){
    ## Calculate beta distribution probabilities dataframe via parallel function
    bpcalc <- function(cp,pl,a,b) {
      q.cer.per = seq(0,1,length.out=length(cp))
      q.per.len = q.cer.per[match(pl,cp)]
      out <- pbeta(max(q.per.len), a, b, lower.tail=TRUE) - pbeta(min(q.per.len), a, b, lower.tail=TRUE)
      if (is.na(out)) {out <- 0}
      if (out<0) {out <- 0}
      return(out)}
    
    prob <- foreach(i=1:length(cer.per),.combine='rbind') %:% foreach(j=1:length(per.len),.combine='c') %dopar% (bpcalc(cp=cer.per[[i]], pl=per.len[[j]], a=alpha[i], b=beta[i]))
    if (length(prob)==1) {prob <- as.matrix(prob)}
    colnames(prob) <- col.nam
  }
  
  ### Truncated Normal Popularity Curve Method
  if (pc.method == "tnormal"){
    ## Calculate beta distribution probabilities dataframe via parallel function
    tnpcalc <- function(cp,pl,zq) {
      q.cer.per = seq(-2,2,length.out=length(cp))
      q.per.len = q.cer.per[match(pl,cp)]
      if (is.na(max(q.per.len))) {out <- 0} else {
        out <- msm::ptnorm(max(q.per.len), 0, 1, lower=-zq, upper=zq, lower.tail=TRUE) - msm::ptnorm(min(q.per.len), 0, 1, lower=-zq, upper=zq, lower.tail=TRUE)
      }
      if (out<0) {out <- 0}
      return(out)}
    
    prob <- foreach(i=1:length(cer.per),.combine='rbind') %:% foreach(j=1:length(per.len),.combine='c') %dopar% (tnpcalc(cp=cer.per[[i]], pl=per.len[[j]], zq=z[i]))
    if (length(prob)==1) {prob <- as.matrix(prob)}
    colnames(prob) <- col.nam
  }     
  
  ### Uniform Popularity Curve Method
  if (pc.method == "uniform"){
    ## Create Uniform Distance dataframe via parallel function
    ucalc <- function(a,b) {
      out <- (length(intersect(a,b))-1)/(length(a)-1)
      if (out<0) {out <- 0}
      return(out)}
    
    prob <- foreach(a=cer.per,.combine='rbind') %:% foreach(b=per.len,.combine='c') %dopar% (ucalc(a,b))
    if (length(prob)==1) {prob <- as.matrix(prob)}
    colnames(prob) <- col.nam
  }
  
  ## Create dataframe of prior values and sum to calculate prior probabilities
  prior.tab <- prob*ct
  if (length(prior.tab)==1) {prior.tab <- as.matrix(prior.tab)}
  if (ncol(prior.tab)>1) {
    if (length(prior.tab[which(chron==1),])==ncol(prior.tab)) {prior <- prior.tab[which(chron==1),]/sum(prior.tab[which(chron==1),])} else
    {prior <- colSums(prior.tab[which(chron==1),])/sum(prior.tab[which(chron==1),])}} else {
      prior <- prob[1,]
    }
  prior[is.na(prior)] <- 0
  
  ## methods using standard normal distribution conditional 
  
  if (method == "bayesian"){
    
    ## Create dataframe of count probabilities
    pij <- sweep(prior.tab,2,colSums(prior.tab),'/')
    pij[is.na(pij)] <- 0
    
    
    ## Create dataframe of probabilities based on popularity curve
    pcij <- sweep(prob,2,colSums(prob),'/')
    pcij[is.na(pcij)] <- 0
    
    ## Create dataframe of standard deviations of uniform distribution probabilities
    sd.t <- sqrt(sweep((pcij*(1-pcij)),2,colSums(ceiling(prob)),'/'))
    
    ## Create dataframe of conditionals and calculate conditional probabilities
    cij <- dnorm(pij,pcij,sd.t)
    cij[which(pcij-sd.t==1)] <- dnorm(1,1,1) # For intervals with only one ceramic type, set to standard normal
    is.na(cij) <- !is.finite(cij) 
    
    # calculate mean conditional probability and rescale from 0-1
    conditional <- apply(cij,2,mean,na.rm=T)
    conditional[is.na(conditional)] <- 0
    if (sum(conditional)>0) {conditional <- conditional/sum(conditional)}
    
    ## Calculate posterior proportions and remove any generated NA
    posterior <- foreach(i=1:length(prior),.combine='c') %dopar% ((prior[i]*conditional[i]))
    posterior <- posterior/sum(posterior)
    posterior[is.na(posterior)] <- 0
    
    ## Deal with edge cases where sum of posterior probabilities = 0 or conditional = prior
    if(sum(posterior)==0) {posterior <- prior}
    if (identical(posterior,prior)) {posterior <- prior}
    
    # calculated beginning (lwr) and ending (upr) dates based on the first and last period with posterior probability about the selected threshold
    lwr <- period[min(which(posterior>cutoff*max(posterior))),1]
    upr <- period[max(which(posterior>cutoff*max(posterior))),2]
    
    # estimate population for time intervals by rescaling to max pop est for site
    pop = scales::rescale(posterior, to = c(min(pop.ests),max(pop.ests)))
    
    ## Create output dataframe
    out <- data.frame(AggID = ID, SubOccSeqLoc = site2, AggSite = ph.sites, Period = ph.periods, Interval = col.nam, 
                      PeriodBegin = pl$PeriodBegin, PeriodBegin.era =pl$PeriodBegin.era, PeriodMidpoint = pl$PeriodMidpoint, 
                      PeriodMidpoint.era = pl$PeriodMidpoint.era, PeriodEnd = pl$PeriodEnd, 
                      PeriodEnd.era = pl$PeriodEnd.era, PeriodLength = plen, Prior=prior, Conditional=conditional,
                      Posterior =posterior, Population = pop, Log_Population = log(pop), Assemb = colSums(prior.tab))
    
  }
  
  ## methods using empirical likelihoods (user input site assemblages for main periods and overlap sites)
  if (method == "mean.obs"){
    
    ## give likelihood input vector names for modelling periods
    obs <- obs.input
    names(obs) <- col.nam
    
    ## calculate likelihood probabilities from likelihood input
    obs[is.na(obs)] <- 0
    if (sum(obs)>0) {obs <- obs/sum(obs)}
    
    ## Calculate posterior proportions and remove any generated NA
    mean.obs <- foreach(i=1:length(prior),.combine='c') %dopar% (mean(c(prior[i],obs[i])))
    mean.obs <- mean.obs/sum(mean.obs)
    mean.obs[is.na(mean.obs)] <- 0
    
    # calculated beginning (lwr) and ending (upr) dates based on the first and last period with posterior probability about the selected threshold
    lwr <- period[min(which(mean.obs>cutoff*max(mean.obs))),1]
    upr <- period[max(which(mean.obs>cutoff*max(mean.obs))),2]
    
    # estimate population for time intervals by rescaling to max pop est for site
    pop = scales::rescale(mean.obs, to = c(min(pop.ests),max(pop.ests)))
    
    ## Create output dataframe
    out <- data.frame(AggID = ID, SubOccSeqLoc = site2, AggSite = ph.sites, Period = ph.periods, Interval = col.nam, 
                      PeriodBegin = pl$PeriodBegin, PeriodBegin.era =pl$PeriodBegin.era, PeriodMidpoint = pl$PeriodMidpoint, 
                      PeriodMidpoint.era = pl$PeriodMidpoint.era, PeriodEnd = pl$PeriodEnd, 
                      PeriodEnd.era = pl$PeriodEnd.era, PeriodLength = plen, Prior=prior, Observed=obs,
                      MeanOccuProb =mean.obs, Population = pop, Log_Population = log(pop), Assemb = colSums(prior.tab))
  }
  
  rownames(out) = ID
  # calculate exponential growth rates using Pert equation via approximation method of Kintigh & Peeples (2020)
  #PertTab <- data.frame(Period = c(out$Period)
  
  
  out$r12_Pert <- NA
  for (i in 1:(nrow(out))) {
    out$r12_Pert[i] <- ((out$Population[i+1]/out$Population[i])^(1/(abs((out$PeriodMidpoint[i]) - (out$PeriodMidpoint[i+1])))))-1
    ## If site is abandoned, calculate growth rate from Period Midpoint to EndPoint, from population to 1
    if (is.na(out$r12_Pert[i])){
      if (out$PeriodEnd[i] == 1520) {out$r12_Pert[i] <- NA} else {
        out$r12_Pert[i] <- ((1/out$Population[i])^(1/(abs((out$PeriodMidpoint[i]) - (out$PeriodEnd[i])))))-1
      }
    }
  }
  
  out <- out %>% mutate(Pct_deltaPop12 = lead(Population, n=1, default = 1) / Population,
                        Pct_deltaPop01 = Population / lag(Population, n=1, default = 1),
                        r01_Pert = lag(r12_Pert, n=1),
                        r23_Pert = lead(r12_Pert, n=1),
                        helper = lag(PeriodLength/2, default = 0)) %>% rowwise() %>%
                 mutate(rPert = mean(c(rep(r12_Pert,(PeriodLength/2)), rep(r01_Pert,helper)),na.rm=T)) %>% ungroup() %>%
                 mutate(rPert0 = lag(rPert, n=1),
                        rPert2 = lead(rPert, n=1)) %>% 
                 rowwise() %>% 
                 mutate(Pct_deltaPop12 = ifelse(PeriodEnd == 1520, NA, Pct_deltaPop12),
                        UrbanPop = ifelse(Population > 1500, Population - UrbanThresh, 0),
                        UrbanOccu = ifelse(Population > 1500, PeriodLength, 0),
                        helper = ifelse(Population > 1500, 1, 0)) %>%
                 ungroup() %>%
                 mutate(OccuTime = cumsum(PeriodLength),
                        OccuInertia = cumsum((PeriodLength*Population)))
  
  out$UrbOccuTime <- as.numeric(unlist(by(out$UrbanOccu, rleid(out$helper), cumsum)))
  out$UrbOccuInertia <- as.numeric(unlist(by((out$UrbanPop*out$PeriodLength), rleid(out$helper), cumsum)))
  out <- out %>% dplyr::select(-UrbanOccu, -helper) %>% 
                 mutate(Pct_UrbdeltaPop12 =  ifelse(lead(UrbanPop, n=1, default = 1) == 0, 1, lead(UrbanPop, n=1, default = 1))/UrbanPop,
                        Pct_UrbdeltaPop12 =  ifelse(Pct_UrbdeltaPop12 == 0, NA, Pct_UrbdeltaPop12),
                        Pct_UrbdeltaPop01 = UrbanPop / lag(UrbanPop, n=1, default = 0)) %>%
                 rowwise() %>% 
                 mutate(Pct_UrbdeltaPop12 = ifelse(PeriodEnd == 1520, NA, Pct_UrbdeltaPop12),
                        UP = ifelse(UrbanPop == 0, 1, UrbanPop)) %>% ungroup()
                 
  
  out$Urb_r12_Pert <- NA
  
  for (i in 1:(nrow(out))) {
    out$Urb_r12_Pert[i] <- ((out$UP[i+1]/out$UP[i])^(1/(abs((out$PeriodMidpoint[i]) - (out$PeriodMidpoint[i+1])))))-1
    ## If site is abandoned, calculate growth rate from Period Midpoint to EndPoint, from population to 1
    if (is.na(out$Urb_r12_Pert[i])){
      if (out$PeriodEnd[i] == 1520) {out$Urb_r12_Pert[i] <- NA} else {
        out$Urb_r12_Pert[i] <- ((1/out$UP[i])^(1/(abs((out$PeriodMidpoint[i]) - (out$PeriodEnd[i])))))-1
      }
    }
  }
  out <- out %>% mutate(Urb_r12_Pert = ifelse(Urb_r12_Pert == 0, NA, Urb_r12_Pert),
                        Urb_r01_Pert = lag(Urb_r12_Pert, n=1),
                        Urb_r23_Pert = lead(Urb_r12_Pert, n=1),
                        helper = lag(PeriodLength/2, default = 0)) %>% rowwise() %>%
                 mutate(Urb_rPert = mean(c(rep(Urb_r12_Pert,(PeriodLength/2)), rep(Urb_r01_Pert,helper)),na.rm=T)) %>% ungroup() %>%
                 mutate(Urb_rPert0 = lag(rPert, n=1),
                        Urb_rPert2 = lead(rPert, n=1)) %>% dplyr::select(-UP, -helper)
  
  out <- do.call(data.frame, lapply(out, function(x) replace(x, is.infinite(x), NA)))
  out <- do.call(data.frame, lapply(out, function(x) replace(x, is.nan(x), NA)))
                                                    
  return(out)
}


