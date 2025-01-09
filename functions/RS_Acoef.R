#### RS_Acoef.R
#### Rudolf Cesaretti, 6/10/2022

#### "RS_Acoef" 
#### 
#### 
#### 
#### 

pak <- c("rgdal", "sp", "sf", "GISTools", "lwgeom", "tidyverse", "tidyr", "data.table", "zoo", "scales")
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

  #z=lnx
  #ids=lnx_id
  #z = c(5000, 3000, 2000, 1500, 1000, 800, 500, 500, 500, 200)
  #ids = c(LETTERS[1:length(X)])
  #setup data

#z=Poly_List[[2]]$Population.s2
#ids=Poly_List[[2]]$AggSite
#plot_title = paste0(2,". ",nam[2])
#yaxis_title = "Log Population"
#Acoef_rescale = "Zipf"

###############################################################
##########################  RS_Acoef  #########################
###############################################################

RS_Acoef <- function(z, 
                     ids,
                     return.data = F,
                     plot.results = T,
                     plots = "Both",#c("RS","Acoef","Both"),
                     plot_title = "Rank-Size A-Coefficient",
                     xaxis_title = "Log Rank",
                     yaxis_title = "Log Size",
                     Acoef_rescale = "Zipf", # "OLS"
                     trunc = F,
                     trunc_at = NULL,
					 plot_zipf = T){
                     
  df <- data.frame(ID = ids, Z = z)
  if (trunc == T){df <- subset(df, Z >= trunc_at)}
  
  dfc <- df[complete.cases(df), ]
  dfc$LogZ <- log(dfc$Z)
  dfc$Rank <- rank(-dfc$LogZ, ties.method = "random")
  dfc$LogRank <- log(dfc$Rank)
  
  dfc$C <- max(dfc$LogZ)
  fit <- lm(LogZ ~ I(LogRank) - 1, data = dfc, offset = C)
  
  fitsum <- summary(fit)
  
  dfc$LogZ_ols.pred <- fit$fitted.values #values predicted by OLS line
  dfc$LogZ_zipf.pred <- (-1*dfc$LogRank) + fit$offset[1] #values predicted by Zipfs Law

  # rescale the log-log data
  # Drennan and Peterson's (2004) original method rescales against the Zipf's law line
  # log(Pop smallest settlement)-log(Pop largest settlement)=sqrt(2)
  #and log(Pop smallest settlement)-log(Rank largest settlement)=sqrt(2)
  if (Acoef_rescale == "Zipf"){
    dfc$LogZ_zipf.pred.resc <- scales::rescale(dfc$LogZ_zipf.pred, to = c(0, sqrt(2)))
    dfc$LogRank.resc <- scales::rescale(dfc$LogRank, to = c(0, sqrt(2)))
    
    fit2 <- lm(LogZ_zipf.pred.resc ~ LogZ_zipf.pred, data = dfc)
    dfc$LogZ.resc <- as.numeric(fit2$coefficients[2])*dfc$LogZ + as.numeric(fit2$coefficients[1])
    dfc$LogZ_ols.pred.resc <- as.numeric(fit2$coefficients[2])*dfc$LogZ_ols.pred + as.numeric(fit2$coefficients[1])
  }
  #alternatively, we could rescale against the OLS line
  if (Acoef_rescale == "OLS"){
    dfc$LogZ_ols.pred.resc <- scales::rescale(dfc$LogZ_ols.pred, to = c(0, sqrt(2)))
    dfc$LogRank.resc <- scales::rescale(dfc$LogRank, to = c(0, sqrt(2)))
    
    fit2 <- lm(LogZ_ols.pred.resc ~ LogZ_ols.pred, data = dfc)
    dfc$LogZ.resc <- as.numeric(fit2$coefficients[2])*dfc$LogZ + as.numeric(fit2$coefficients[1])
    dfc$LogZ_zipf.pred.resc <- as.numeric(fit2$coefficients[2])*dfc$LogZ_zipf.pred + as.numeric(fit2$coefficients[1])
  }

  #Compute the area above the Zipfs law diagonal and below the observed
  #rank-size curve (A1), and then the area below the diagonal and above the
  #empirical data (A2)
  
  #first, if there are duplicate values, we need to jitter them 
  #by a very small amount so that they are unique
  un = match(unique(dfc$LogZ.resc), dfc$LogZ.resc)
  vv = 1:length(dfc$LogZ.resc)
  vv = setdiff(vv,un)
  if (length(vv)>0){
    dfc$LogZ.resc[vv] <- seq(-0.0001,length(vv)*-0.0001, by=-0.0001)+dfc$LogZ.resc[vv]
  }
  
  #numerically fit functions for the dividing line (Zipf or OLS) and the observed data
  if (Acoef_rescale == "OLS"){line <- approxfun(dfc$LogRank.resc, dfc$LogZ_ols.pred.resc)}
  if (Acoef_rescale == "Zipf"){line <- approxfun(dfc$LogRank.resc, dfc$LogZ_zipf.pred.resc)}
  obs <- approxfun(dfc$LogRank.resc, dfc$LogZ.resc)
  
  # defining x range and dividing it to sections (for example n=500)
  i <- seq(0, sqrt(2), length.out=500)
  
  # calculating the distance between the density curves
  h1a <- obs(i)-line(i)
  h2a <- line(i)-obs(i)
  h1 <- ifelse(h1a < 0, 0, h1a)
  h2 <- ifelse(h2a < 0, 0, h2a)
  
  #and using the trapezoidal rule here to approximate the integrals
  # about the RS line
  A1 <- trapz(i, h1)  # for the regions where line>obs
  A2 <- trapz(i, h2)  # for the regions where obs>line
  
  #Finally, compute A as the difference A1 - A2
  A <- A1 - A2 
  
  #class_df <- data.frame(i=i,QuantileX=scales::rescale(i),
  #                                Position=h1a,
  #                                A1_Position=h1,
  #                                A2_Position=h2)
  #if(abs(class_df[1,c("Position")]) < 0.00001){class_df <- class_df[-1,]}
  #class_df <- class_df %>%  rowwise() %>% 
  #                          mutate(PosClassNum=ifelse(Position < 0, 2, 1)) %>%
  #                          ungroup() %>%
  #                          mutate(Set = data.table::rleid(PosClassNum)) %>%
  #                          group_by(Set) %>%
  #                          summarize(integA1 = trapz(i, A1_Position),
  #                                    integA2 = trapz(i, A2_Position),
  #                                    Integral = integA1+integA2,
  #                                    PosClass = paste0("A",mean(PosClassNum)),
  #                                    MaxPosition = max(abs(Position)),
  #                                    PctTotal = max(QuantileX)-min(QuantileX)) %>% ungroup()
  
  #classification <- "Unclassified"

  #if (class_df[1,c("PosClass")]=="A2" & class_df[1,c("PctTotal")] >= 0.7){classification <- "Primate"}
  #if (class_df[1,c("PosClass")]=="A1" & class_df[1,c("PctTotal")] >= 0.7){classification <- "Convex"}
  #if (class_df[1,c("PosClass")]=="A1" & class_df[2,c("PosClass")]=="A2" & sum(class_df[1:2,c("MaxPosition")]) >= 0.8 & classification == "Unclassified"){classification <- "Convex-Primate"}
  #if (class_df[1,c("PosClass")]=="A2" & class_df[2,c("PosClass")]=="A1" & sum(class_df[1:2,c("MaxPosition")]) >= 0.8 & classification == "Unclassified"){classification <- "Primo-Convex"}
  #if (A >= 0.15 & class_df[1,c("PosClass")]=="A1" & class_df[2,c("PosClass")]=="A2" & sum(class_df[,c("integA2")]) <= 0.15 & classification == "Unclassified"){classification <- "Double-Convex"}
  #if (A <= -0.15 & class_df[1,c("PosClass")]=="A2" & class_df[2,c("PosClass")]=="A1" & sum(class_df[,c("integA1")]) <= 0.15 & classification == "Unclassified"){classification <- "Double-Primate"}
  
  #if (class_df$PosClass[1] == "A1" & A >= -0.1 & A <= 0.1 & A1 < 0.15 & A2 < 0.15 & classification == "Unclassified"){classification <- "Snaking Log Normal - Convex Head"}
  #if (class_df$PosClass[1] == "A2" & A >= -0.1 & A <= 0.1 & A1 < 0.15 & A2 < 0.15 & classification == "Unclassified"){classification <- "Snaking Log Normal - Convex Head"}
  #if (A >= -0.03 & A <= 0.03 & A1 < 0.07 & A2 < 0.07 & max(class_df$MaxPosition) < 0.25){classification <- "Near Log Normal"}
  
  m <- c("Slope", "Slope_se", "Slope_t", "Slope_p", "Intercept", "R2", "n", 
         "Ymin", "Ymax", "A1", "A2", "A", 
         #"Classification", 
         "Rescale Line")
  v <- c(fitsum$coefficients[1], fitsum$coefficients[2], fitsum$coefficients[3], 
         fitsum$coefficients[4], fit$offset[1], fitsum$adj.r.squared, 
         length(z), range(dfc$LogZ)[1], range(dfc$LogZ)[2], 
         A1, A2, A, 
         #classification, 
         Acoef_rescale)
  model.df <- data.frame(Metric = m, Value = v)
  data.df <- dfc
  out_list <- list()
  j=1
  
  if (plot.results == T){
    
    if (plots == "RS" | plots == "Both"){
      
      x <- c(dfc$LogRank,dfc$LogRank,dfc$LogRank)
      y1 <- dfc$LogZ_zipf.pred
      y2 <- dfc$LogZ_ols.pred
      y3 <- dfc$LogZ
      l1 <- rep("Zipf's Law", length(dfc$LogRank))
      l2 <- rep("OLS Estimate", length(dfc$LogRank))
      l3 <- rep("Observed", length(dfc$LogRank))
	  lines <- data.frame(x = x, y = c(y1, y2, y3), lab = c(l1, l2, l3))
	  
	  if (plot_zipf == F){
		lines <- subset(lines, lab != "Zipf's Law")
		cols <- c("Observed" = "black", 
                "OLS Estimate" = "firebrick")
		y1 = max(dfc$LogZ)+0.2		
		y2 = min(dfc$LogZ)-0.2	
	  }else{
		cols <- c("Observed" = "black", 
                "OLS Estimate" = "firebrick", 
                "Zipf's Law" = "forestgreen")
	  }
	  
      
                      
      OLSlab = paste0("atop(bold(italic(LogPop) == ",
                      round(fit$offset[1],2),
                      round(fitsum$coefficients[1],2),
                      "*italic(LogRank)),OLS ~ Regression ~ bold(italic(R^2) == ", 
                      round(fitsum$adj.r.squared,2),
                      " ~~ italic(n) == ",
                      length(z),
                      " ))")
      
      labypos = max(dfc$LogZ)*0.33
      labxpos = max(dfc$LogRank)*0.33
	  

      RS_plt <- ggplot() + geom_line(data=lines, aes(x,y,color=factor(lab)), size=1.2) + 
        scale_colour_manual(values = cols)+
        guides(color = guide_legend(nrow = 3, byrow = TRUE))+
        geom_point(data=dfc, aes(x=LogRank, y=LogZ), shape=19, size=1.5, color="black") + 
        annotate("label", x = labxpos, y = labypos, label = OLSlab, parse = TRUE, color ="firebrick", size=4)+
        labs(x = xaxis_title, y = yaxis_title) +
        theme_bw()+
        theme(
          axis.text.x = element_text(color="black", size=12, face="bold"), 
          axis.text.y = element_text(color="black", size=12, face="bold"), 
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          legend.position="right", legend.title=element_blank(),
          legend.spacing.y = unit(0.5,"cm"),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.text = element_text(color="black", size=10, face="bold"))
		  
	if (plot_zipf == F){RS_plt = RS_plt + coord_cartesian(ylim = c(y2, y1))}
      
    }
    
    if (plots == "Acoef" | plots == "Both"){
      
      dat = data.frame(x=c(i,i), y=c(obs(i),line(i)), y2=c(pmax(obs(i),line(i)),pmin(obs(i),line(i))),
                       line = c(rep("obs",length(i)),rep("zipf",length(i))), h=c(h1,h1))
      
      dat2 = data.frame(x=i, obs=obs(i), line=line(i), h1=h1, h2=h2)
      
      dat2 <- dat2 %>% rowwise() %>% 
        mutate(ymin = min(obs,line), ymax= max(obs,line),
               group = ifelse(h1 >= h2, "A1", "A2")) %>% ungroup()
      
      dat2$id = ave(dat2$h2,dat2$group, FUN=seq_along)
      dat2$id2 = cumsum(c(1, abs(diff(dat2$id)) > 1))
      
      ######################
      xs = c(0, sqrt(2))
      beta = c(sqrt(2), -1)
      ys = cbind(1, xs) %*% beta
      
	  
     
	 
	  
	  
      x <- c(dfc$LogRank.resc,dfc$LogRank.resc,dfc$LogRank.resc)
      y1 <- dfc$LogZ_zipf.pred.resc
      y2 <- dfc$LogZ_ols.pred.resc
      y3 <- dfc$LogZ.resc
      l1 <- rep("Zipf's Law", length(dfc$LogRank.resc))
      l2 <- rep("OLS Estimate", length(dfc$LogRank))
      l3 <- rep("Observed", length(dfc$LogRank))
      lines2 <- data.frame(x = x, y = c(y1, y2, y3), lab = c(l1, l2, l3))
      
      
											
	if (plot_zipf == F){
		lines2 <- subset(lines2, lab != "Zipf's Law")
		if (Acoef_rescale == "OLS"){cols2 <- c("Observed" = "black", 
                                            "OLS Estimate" = "black")}
		if (Acoef_rescale == "Zipf"){cols2 <- c("Observed" = "black",
                                             "OLS Estimate" = "firebrick")}
		y1 = max(dfc$LogZ)+0.2		
		y2 = min(dfc$LogZ)-0.2	
	  }else{
		if (Acoef_rescale == "Zipf"){cols2 <- c("Observed" = "black",
                                             "Zipf's Law" = "black",
                                             "OLS Estimate" = "firebrick")}
        if (Acoef_rescale == "OLS"){cols2 <- c("Observed" = "black", 
                                            "OLS Estimate" = "black", 
                                            "Zipf's Law" = "forestgreen")}
	  }
      
      ###########################
      
      a1t= paste0("A1=",round(A1,3))
      a2t= paste0("A2=",round(A2,3))
      at = paste0("A = ",round(A,4))
      
      xaxt <- paste0("Rescaled ", xaxis_title, " [0,sqrt(2)]")
      yaxt <- paste0("Rescaled ", yaxis_title, " [0,sqrt(2)]")
      
      Acoef_plt <- ggplot() + 
        geom_ribbon(data = dat2 %>% filter(obs >= line),
                    aes(x=x, ymin = line, ymax = obs, fill="dodgerblue"), 
                    alpha=0.7) +
        geom_ribbon(data = dat2 %>% filter(obs <= line),
                    aes(x=x, ymin = obs, ymax = line, fill = "coral1"), 
                    alpha=0.7) +
        scale_fill_identity(name = at,
                            breaks = c("dodgerblue", "coral1"),#firebricl
                            labels = c(a1t, a2t),
                            guide = "legend")+
        geom_line(data=lines2, aes(x,y,color=factor(lab)), size=1) +  
        scale_colour_manual(values = cols2)+
        geom_point(data=dfc, aes(x=LogRank.resc,y=LogZ.resc), shape=19, size=1.7, color="black")+
        guides(fill = guide_legend(nrow = 2, byrow = TRUE))+
        labs(x = xaxt, y = yaxt, colour = "") +
        theme_bw()+
        theme(
          axis.text.x = element_text(color="black", size=12, face="bold"), 
          axis.text.y = element_text(color="black", size=12, face="bold"), 
          axis.title.x = element_text(color="black", size=14, face="bold"),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          legend.position="right", legend.spacing.y = unit(0.5,"cm"),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.title = element_text(color="black", size=14, face="bold"),
          legend.text = element_text(color="black", size=12, face="bold"))
		  
      if (plot_zipf == F){Acoef_plt = Acoef_plt + coord_cartesian(ylim = c(y2, y1))}
	  
    }
    if (plots == "RS"){
      out_plot <- RS_plt + labs(title = plot_title, subtitle = "Rank-Size Plot") + 
        theme(plot.title = element_text(hjust = 0.5, face="bold"), 
              plot.subtitle = element_text(hjust = 0.5, face="bold"))
      
      
    }
    if (plots == "Acoef"){
      out_plot <- Acoef_plt + labs(title = plot_title, subtitle = "A-Coefficient Plot") + 
        theme(plot.title = element_text(hjust = 0.5, face="bold"), 
              plot.subtitle = element_text(hjust = 0.5, face="bold"))
      
      
    }
    if (plots == "Both"){
      RS_plt <- RS_plt + labs(subtitle = "Rank-Size Plot") + 
        guides(color = guide_legend(nrow = 1))+
        theme(plot.subtitle = element_text(hjust = 0.5, face="bold", size=14),
              legend.position = "bottom")
      Acoef_plt <- Acoef_plt + labs(subtitle = "A-Coefficient Plot") + 
        guides(fill = guide_legend(nrow = 1),color="none")+
        theme(plot.subtitle = element_text(hjust = 0.5, face="bold", size=14),
              legend.position = "bottom")
      plot_row <- plot_grid(RS_plt, Acoef_plt, align = "h", nrow = 1)
      title <- ggdraw() + 
        draw_label(
          plot_title,
          fontface = 'bold',
          x = 0.5,
          hjust = 0.5,
          size=24
        ) +
        theme(
          # add margin on the left of the drawing canvas,
          # so title is aligned with left edge of first plot
          plot.margin = margin(0, 0, 0, 7)
        )
      out_plot <- plot_grid(
        title,
        plot_row,
        ncol = 1,
        # rel_heights values control vertical title margins
        rel_heights = c(0.1, 1)
      )
      
    }
    
    out_list[[j]] <- out_plot
    j=j+1
    
  }
  
  out_list[[j]] <- model.df
  j=j+1
  
  if (return.data == T){
    out_list[[j]] <- data.df
  }
  
  return(out_list)
  
}


