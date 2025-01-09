#### spda.ggplot.R
#### Rudolf Cesaretti, 6/1/2022


#### "spda.ggplot" is a function for plotting the output 
#### of the spda function as a time series of population,
#### assemblages, and occupational probabilities for each period.
#### Has parameter options for method-specific inputs, plot 
#### aesthetics, whether to plot Pert pop trajectory, and whether
#### to use a log scale for population


pak <- c("tidyverse", "era", "zoo", "scales", "ggnewscale")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)


#############################################################################
#################################  spda.ggplot  #############################
#############################################################################

spda.ggplot <- function(data, ## dataframe of ChronApportion output for a single site
                        ## or a different (merged) dataframe that includes the 
                        ## exact same variables (w/ same names)
                        xbreaks,
                        xlabels,
                        method = c("bayesian","mean.obs"),
                        pc.method = c("beta","tnorm","uniform","input"), #which method was used in the ChronApportion function
                        legend_position = c(0.01,0.99), #"bottom", "top", "none", "left", "right"
                        Pop.col = "black", 
                        Pert.col ="orange2",
                        Prior.col = "red",
                        Cond.Obs.col = "blue",
                        Post.OccPrb.col = "green2",
                        Assemb.col = "purple", 
                        PopLogScale = FALSE,
                        Assemb.plot = FALSE,
                        Pert.plot = TRUE,
                        yaxtitle = "P(x)"){
  
  
  ### Save data as vectors and dataframes
  prior <- data$Prior
  if (method == "bayesian"){
    conditional <- data$Conditional
    posterior <- data$Posterior}
  if (method == "mean.obs"){
    Observed <- data$Observed
    OccuProb <- data$MeanOccuProb}
  
  phases <- data$Period
  assemb <- data$Assemb
  p <- data[,c("PeriodBegin", "PeriodBegin.era", "PeriodEnd", "PeriodEnd.era", 
               "PeriodMidpoint", "PeriodMidpoint.era", "Population", "Log_Population", "r12_Pert")]
  site <- data$SubOccSeqLoc[1]
  pop <- data$Population
  
  ### Calculate Pert equation population trajectory
  if (Pert.plot == TRUE){
    # create empty lists for forloop
    t <- list()
    N0 <- list()
    r12 <- list()
    t1 <- list()
    # Iterate the equation for input parameters
    for (i in 1:(nrow(p))) {
      if (i == nrow(p)){
        t[[i]] <- seq(p$PeriodMidpoint[i], (p$PeriodEnd[i]), by=0.1)
      } else{
        t[[i]] <- seq(p$PeriodMidpoint[i], (p$PeriodMidpoint[i+1]-0.1), by=0.1)
      }
      r12[[i]] <- rep(p$r12_Pert[i],length.out = length(t[[i]]))
      N0[[i]] <- rep(p$Population[i],length.out = length(t[[i]]))
      t1[[i]] <- seq(0, ((length(t[[i]])/10)-0.1), by = 0.1)
    }
    # convert outputs into dataframe
    Pert = data.frame(t = unlist(t), r12=unlist(r12), N0 = unlist(N0), t1=unlist(t1))
    # calculate pop trajectory
    Pert$N1 = Pert$N0*((1+Pert$r12)^Pert$t1)
    # if Log scale is desired, take the log of N1, and remove rows with pop = 0 (NA)
    if (PopLogScale == TRUE) {
      Pert = Pert %>% filter(N1 > 0)
      Pert$N1 = log(Pert$N1)
    }
  }
  
  ### Polygon calculation based on method and parameters
  g.per <- as.numeric(as.vector(t(as.matrix(p[,c(1,3)]))))
  g.per <- c(g.per[1],g.per,g.per[length(g.per)])
  
  g.pop <- as.vector(t(cbind(pop,pop)))
  g.pop <- c(0,g.pop,0)
  if (PopLogScale == TRUE){# if using log scale
    g.pop = log(g.pop) #take natural logarithm
    g.pop[is.infinite(g.pop)] <- 0#replace Inf and -Inf with 0 (bc the scale is transformed from [0,1] so it needs zero values)
    p$Population <- p$Log_Population #change the Population column in table p to Log_Population for plotting points in ggplot
  }
  
  if (method == "bayesian"){
    g.post <- as.vector(t(cbind(posterior,posterior)))
    g.post <- c(0,g.post,0)
    g.prior <- as.vector(t(cbind(prior,prior)))
    g.prior <- c(0,g.prior,0)
    g.con <- as.vector(t(cbind(conditional,conditional)))
    g.con <- c(0,g.con,0)
    if (Assemb.plot == TRUE){
      g.assemb <- as.vector(t(cbind(assemb,assemb)))
      g.assemb <- c(0,g.assemb,0)
      polydf <- data.frame(Years = g.per, Posterior = g.post, Prior = g.prior, Conditional = g.con, Population = g.pop, Assemb = g.assemb)
    } else {
      polydf <- data.frame(Years = g.per, Posterior = g.post, Prior = g.prior, Conditional = g.con, Population = g.pop)
    }
    # Custom titles based on method
    if (pc.method=="input"){TitleMethod <- "Bayesian & Custom Empirical Popularity Curves"}
    if (pc.method=="beta"){TitleMethod <- "Bayesian & Custom Beta Popularity Curves"}
    if (pc.method=="tnorm"){TitleMethod <- "Bayesian & Truncated Normal Popularity Curves"}
    if (pc.method=="uniform"){TitleMethod <- "Bayesian & Uniform Popularity Curves"}
    # Factors for rescaling the graph ( )
    # divide values by this; multiply second axis by it
    sc.pop = (max(polydf$Posterior))/(max(polydf$Population)-min(polydf$Population))
    sc.log.pop = (max(polydf$Posterior))/(max(polydf$Population)-(min(polydf$Population)))
    scale.pop <- ifelse(PopLogScale == TRUE, sc.log.pop, sc.pop)
    yscale.pr <- ceiling(max(c(polydf$Prior,polydf$Conditional,polydf$Posterior))*10)/10 
    if (Assemb.plot == TRUE){scale.assemb <- (max(c(polydf$Prior,polydf$Conditional,polydf$Posterior)))/(max(polydf$Assemb)-min(polydf$Assemb)) }
  }
  
  if (method == "mean.obs"){
    g.occuprob <- as.vector(t(cbind(OccuProb,OccuProb)))
    g.occuprob <- c(0,g.occuprob,0)
    g.prior <- as.vector(t(cbind(prior,prior)))
    g.prior <- c(0,g.prior,0)
    g.obs <- as.vector(t(cbind(Observed,Observed)))
    g.obs <- c(0,g.obs,0)
    if (Assemb.plot == TRUE){
      g.assemb <- as.vector(t(cbind(assemb,assemb)))
      g.assemb <- c(0,g.assemb,0)
      polydf <- data.frame(Years = g.per, MeanOccuProb = g.occuprob, Prior = g.prior, Observed = g.obs, Population = g.pop, Assemb = g.assemb)
    } else {
      polydf <- data.frame(Years = g.per, MeanOccuProb = g.occuprob, Prior = g.prior, Observed = g.obs, Population = g.pop)
    }
    # Custom titles based on method
    if (pc.method=="input"){TitleMethod <- "Mean Occu. Prob. & Custom Empirical Popularity Curves"}
    if (pc.method=="beta"){TitleMethod <- "Mean Occu. Prob. & Custom Beta Popularity Curves"}
    if (pc.method=="tnorm"){TitleMethod <- "Mean Occu. Prob. & Truncated Normal Popularity Curves"}
    if (pc.method=="uniform"){TitleMethod <- "Mean Occu. Prob. & Uniform Popularity Curves"}
    # Factors for rescaling the graph ( )
    # divide values by this; multiply second axis by it
    sc.pop = (max(polydf$MeanOccuProb))/(max(polydf$Population)-min(polydf$Population))
    sc.log.pop = (max(polydf$MeanOccuProb))/(max(polydf$Population)-(min(polydf$Population)))
    scale.pop <- ifelse(PopLogScale == TRUE, sc.log.pop, sc.pop)
    yscale.pr <- ceiling(max(c(polydf$Prior,polydf$Observed,polydf$MeanOccuProb))*10)/10 
    if (Assemb.plot == TRUE){scale.assemb <- (max(c(polydf$Prior,polydf$Observed,polydf$MeanOccuProb)))/(max(polydf$Assemb)-min(polydf$Assemb)) }
  }
  
  Pop.Axis.Name = ifelse(PopLogScale == TRUE, "Log Population", "Population")
  if (Assemb.plot == TRUE){
    Y_Axis_Title <- paste0(yaxtitle," & Assemb Size")
    YSize <- 12
  } else {
    Y_Axis_Title <- yaxtitle
    YSize <- 14
  }
  
  ### The actual plot
  ggp <- ggplot() + 
    ## Population Polygon
    geom_polygon(polydf, mapping=aes(x=Years, y=Population*scale.pop, fill = "Population"), alpha=0.7, size=0) + 
    ## Posterior/MeanOccuProb Polygon
    {if (method == "bayesian") geom_polygon(polydf, mapping=aes(x=Years, y=Posterior, fill = "Posterior"), color="black", alpha=0.5)}+
    {if (method == "mean.obs") geom_polygon(polydf, mapping=aes(x=Years, y=MeanOccuProb, fill = "Mean OccuProb"), color="black", alpha=0.5)}+
    ## polygon fill color scale
    {if(method == "bayesian") scale_fill_manual(values = c("Posterior" = "grey", "Population" = "black"))}+
    {if (method == "mean.obs") scale_fill_manual(values = c("Mean OccuProb" = "grey", "Population" = "black"))}+
    ## Population Line
    geom_line(polydf, mapping=aes(x=Years,y=Population*scale.pop, color="Population"),size=1.9)+
    ## Posterior/OccuProb Line
    {if(method == "bayesian") geom_line(polydf, mapping=aes(x=Years,y=Posterior, color="Posterior"),size=1)}+
    {if(method == "mean.obs") geom_line(polydf, mapping=aes(x=Years,y=MeanOccuProb, color="Mean OccuProb"),size=1)}+
    ## Pert Line
    {if(Pert.plot == TRUE) geom_line(Pert, mapping=aes(x=t,y=N1*scale.pop, color="Pert Pop Trajectory"),size=1.4)}+
    ## Pert Points
    {if(Pert.plot == TRUE) geom_point(p, mapping=aes(x=PeriodMidpoint, y=Population*scale.pop), shape=24, color="black", fill="orange2", size=3)}+
    ## Prior Line
    geom_line(polydf, mapping=aes(x=Years,y=Prior, color="Prior"),size=1)+
    ## Conditional/Likelihood Line
    {if(method == "bayesian") geom_line(polydf, mapping=aes(x=Years,y=Conditional, color="Conditional"),size=1, linetype="dashed")}+
    {if(method == "mean.obs") geom_line(polydf, mapping=aes(x=Years,y=Observed, color="Observed"),size=1, linetype="dashed")}+
    ## Assemb Line
    {if(Assemb.plot == TRUE) geom_line(polydf, mapping=aes(x=Years,y=Assemb*scale.assemb, color="Assemb"),size=1)}+
    ## line color scale
    #### Custom color scales implemented based on parameters and method
    {if(method == "bayesian" & Pert.plot == TRUE & Assemb.plot == TRUE) scale_color_manual(values = c("Prior" = Prior.col, "Conditional" = Cond.Obs.col, "Posterior" = Post.OccPrb.col, "Population" = Pop.col, "Assemb" = Assemb.col, "Pert Pop Trajectory" = Pert.col))}+
    {if(method == "bayesian" & Pert.plot == FALSE & Assemb.plot == TRUE) scale_color_manual(values = c("Prior" = Prior.col, "Conditional" = Cond.Obs.col, "Posterior" = Post.OccPrb.col, "Population" = Pop.col, "Assemb" = Assemb.col))}+
    {if(method == "bayesian" & Pert.plot == TRUE & Assemb.plot == FALSE) scale_color_manual(values = c("Prior" = Prior.col, "Conditional" = Cond.Obs.col, "Posterior" = Post.OccPrb.col, "Population" = Pop.col, "Pert Pop Trajectory" = Pert.col))}+
    {if(method == "bayesian" & Pert.plot == FALSE & Assemb.plot == FALSE) scale_color_manual(values = c("Prior" = Prior.col, "Conditional" = Cond.Obs.col, "Posterior" = Post.OccPrb.col, "Population" = Pop.col))}+
    {if(method == "mean.obs" & Pert.plot == TRUE & Assemb.plot == TRUE) scale_color_manual(values = c("Prior" = Prior.col, "Observed" = Cond.Obs.col, "Mean OccuProb" = Post.OccPrb.col, "Population" = Pop.col, "Assemb" = Assemb.col, "Pert Pop Trajectory" = Pert.col))}+
    {if(method == "mean.obs" & Pert.plot == FALSE & Assemb.plot == TRUE) scale_color_manual(values = c("Prior" = Prior.col, "Observed" = Cond.Obs.col, "Mean OccuProb" = Post.OccPrb.col, "Population" = Pop.col, "Assemb" = Assemb.col))}+
    {if(method == "mean.obs" & Pert.plot == TRUE & Assemb.plot == FALSE) scale_color_manual(values = c("Prior" = Prior.col, "Observed" = Cond.Obs.col, "Mean OccuProb" = Post.OccPrb.col, "Population" = Pop.col, "Pert Pop Trajectory" = Pert.col))}+
    {if(method == "mean.obs" & Pert.plot == FALSE & Assemb.plot == FALSE) scale_color_manual(values = c("Prior" = Prior.col, "Observed" = Cond.Obs.col, "Mean OccuProb" = Post.OccPrb.col, "Population" = Pop.col))}+
    labs(title=paste0(site,"   (",xlabels[1]," - ",xlabels[length(xlabels)],")"), subtitle=TitleMethod, x ="Years BC/AD", y = Y_Axis_Title, 
         caption = "**Population = posterior probability rescaled to the min and max population estimates from script #1")+
    theme_bw() +
    scale_y_continuous(name = Y_Axis_Title,   sec.axis = sec_axis(~./scale.pop, name=Pop.Axis.Name) )+
    scale_x_continuous(breaks = xbreaks,labels = xlabels) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face="bold", color="black"), 
          axis.text.y = element_text(face="bold", color="black"), 
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=YSize, face="bold"),
          legend.justification=c(0,1), legend.position=legend_position, legend.title=element_blank(), legend.key.size = unit(0.4,"line"),
          legend.spacing.y = unit(0.1,"line"),legend.key = element_rect(colour = "transparent", fill = "white"), 
          plot.title = element_text(hjust = 0.5, face="bold", size=16), plot.subtitle = element_text(hjust = 0.5, face="bold", size=14))
  
  return(ggp)
}


