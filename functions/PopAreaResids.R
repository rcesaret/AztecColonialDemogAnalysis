#### PopAreaResids.R
#### Rudolf Cesaretti, 7/4/2022

#### "PopAreaResids" 
#### 
#### 
#### 
#### 

packages <- c("Matrix","tidyverse","tidyr","cowplot","stars","ggnewscale",
              "broom","zoo","lmtest","sandwich")

#, "data.table",  "mgcv","igraph", "ggrepel","ggridges", "movecost",  "datplot", "scales",

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# load packages
invisible(lapply(packages, library, character.only = TRUE))

rm(packages,installed_packages)

###############################################################
#######################  PopAreaResids  #######################
###############################################################

PopAreaResids <- function(dat,
                          ggplotit=F){
outlist <- list()
x = dat$Log_Population.s2
y = log(dat$Area_ha)
fm <- lm(y ~ x)
dat$Log_Area_ha <- y
tidy_fm <- tidy(fm)
coeftest <- coeftest(fm, vcov = vcovHC(fm, "HC1"))
tidy_coeftest <- tidy(coeftest)
tidy_coeftest$Period <- dat$Period[1]
tidy_coeftest$PeriodNum <- dat$PeriodNum[1]
tidy_coeftest$term[1] <- "alpha"
tidy_coeftest$term[2] <- "beta"
tidy_coeftest$param <- c(1, 2)
tidy_coeftest$n <- c(length(x), length(x))
tidy_coeftest$rsquared <- c(summary(fm)$r.squared, summary(fm)$r.squared)
outlist[[1]] <- tidy_coeftest

resids <- data.frame(AggID=dat$AggID, AggSite=dat$AggSite, Period= dat$Period[1], PeriodNum= dat$PeriodNum[1], LogPop=x, LogArea=y, PredictLogArea=(tidy_coeftest$estimate[1]+tidy_coeftest$estimate[2]*x))
resids$Resid = resids$LogArea - resids$PredictLogArea
outlist[[2]] <- resids
#df$res = fm$residuals
if (ggplotit == T){
  ggp = ggplot(dat, aes(x=Log_Population.s2, y=Log_Area_ha)) +
  geom_point(size=3, shape=21, color = "black", fill="grey") +
  geom_smooth(fullrange=TRUE, method =lm, se=TRUE, colour="red", size=1) +
  geom_abline(mapping = NULL, data = NULL, slope = tidy_coeftest$estimate[2], intercept = tidy_coeftest$estimate[1], na.rm = FALSE, show.legend = NA, colour="black", size=1) + 
  xlab("Log Population") +
  ylab("Log Area (ha)") +
  theme_bw() +
  theme(legend.justification=c(1,0), legend.position=c(0.9,0.1)) +
  theme(legend.box.background = element_rect(colour = "black"), 
        panel.border = element_rect(colour = "black", fill=NA), 
        legend.background = element_blank()) +
  theme(axis.text = element_text(size = 12, color="black")) +
  theme(axis.title.y = element_text(size = 14, face = "bold")) +
  theme(axis.title.x = element_text(size = 14, face = "bold")) +
  theme(legend.title = element_text(face="bold"))
  
  outlist[[3]] <- ggp
  
  
}
return(outlist)
}

