#### RegressionResultsTab.R
#### Rudolf Cesaretti, 11/9/2022

#### "RegressionResultsTab.R"
#### 
#### 
#### 
#### 

pak <- c("tidyverse", "tidyr", "nls.multstart", "nlme", "brms", "broom", 
         "ggrepel", "minpack.lm", "data.table", "zoo", "lmtest", "sandwich", 
         "nlstools", "MASS", "NSM3", "gridExtra", "ggnewscale", "cowplot")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)

###############################################################
###################  RegressionResultsTab  ####################
###############################################################


RegressionResultsTab <- function(model, 
                                 data, 
                                 x, 
                                 y, 
                                 model_type, 
                                 coef="both", 
                                 residplots=TRUE, 
                                 residhists=TRUE){
  

  if (model_type=="ols"){
    
    s = summary(model)
    model_results <- data.frame(matrix(ncol = 0, nrow = 1))
    model_results$formula <- paste(as.character(s$call$formula)[[2]], as.character(s$call$formula)[[1]], as.character(s$call$formula)[[3]])
    model_results$n <- length(model$resid)
    model_results$df <- s$df[2]
    model_results$R2 <- s$r.squared
    model_results$Adj_R2 <- s$adj.r.squared
    model_results$resid_se <- s$sigma

    Results <- list()
    model_results <- list(model_results)
    Results <- append(Results,model_results)
    names(Results) <- c("OLS Model Results")
    fm.resid <- model$residuals
    
    if (coef=="normal" | coef=="both"){
      model_coef <- tidy(model)
      Cov <- vcov(model)
      tt <-qt(c(0.025,0.975),summary(model)$df[2])
      se <- sqrt(diag(Cov))
      ci <-coef(model) + se %o% tt
      model_coef$ci0.025 <- c(exp(ci[1,1]), ci[2,1])
      model_coef$ci0.975 <- c(exp(ci[1,2]), ci[2,2])
      model_coef$EstMethod <- "standard"
      model_coef <- as.data.frame(model_coef)
    }
    
    # sandwich results
    if (coef=="robust" | coef=="both"){
      coeftest <- coeftest(model, vcov = sandwich)
      model_coef_robust <- tidy(coeftest)
      Cov <- sandwich(model)
      tt <-qt(c(0.025,0.975),summary(model)$df[2])
      se <- sqrt(diag(Cov))
      ci <-coef(model) + se %o% tt
      model_coef_robust$ci0.025 <- c(exp(ci[1,1]), ci[2,1])
      model_coef_robust$ci0.975 <- c(exp(ci[1,2]), ci[2,2])
      model_coef_robust$EstMethod <- "robust"
      model_coef_robust <- as.data.frame(model_coef_robust)
    }
    
    if (coef=="normal"){
      model_coef <- list(model_coef)
      nam <- names(Results)
      Results <- append(Results,model_coef)
      names(Results) <- c(nam,"Model Coefficients")
    }
    if (coef=="robust"){
      model_coef_robust <- list(model_coef_robust)
      nam <- names(Results)
      Results <- append(Results,model_coef_robust)
      names(Results) <- c(nam,"Robust Model Coefficients")
    }
    if (coef=="both"){
      model_coef_both = rbind(model_coef,model_coef_robust)
      model_coef_both <- list(model_coef_both)
      nam <- names(Results)
      Results <- append(Results,model_coef_both)
      names(Results) <- c(nam,"Model Coefficients")
    }
    
    
    shap = shapiro.test(fm.resid)
    shap.t = tidy(shap)
    shap.t$alternative <- NA
    
    fit<-fitdistr(fm.resid,"normal")$estimate
    ksr = ks.test(fm.resid, "pnorm",fit[1],fit[2])
    
    ksr.t = tidy(ksr)
    
    resid_normal_test = rbind(shap.t,ksr.t) 
    resid_normal_test = as.data.frame(resid_normal_test) 
    resid_normal_test <- list(resid_normal_test)
    nam <- names(Results)
    Results <- append(Results,resid_normal_test)
    names(Results) <- c(nam,"Residuals Normality Tests")
    
    
    if (residplots==TRUE){
      opar <- par() 
      par(mfrow=c(2,2), mar = c(2, 2, 2, 2))
      plot(model)
      par(opar)
    }
    
    if (residhists==TRUE){
      Resid1 <- hist(fm.resid, plot=F)
      cts = max(Resid1$density) / max(Resid1$counts)
      freq <- pretty(Resid1$counts)
      freq2 = freq * cts
      pos <- pretty(Resid1$density)
      #freq <- round(pos * length(fm.resid) * with(Resid1, breaks[2] - breaks[1]))
      
      graphics:::plot.histogram(Resid1, freq = F, col="gray70", family="sans", main="Standardized Residuals vs. Gaussian Distribution", cex.main=1, 
                                font.main=2, font.lab=2, xlab = "Residuals", ylab="Frequency",
                                border="black", yaxt='n')#, ylim=c(0,1)
      Axis(side = 2, at = freq2, labels = freq)
      Axis(side = 4, at = pos, labels = pos)
      mtext("Density           ", side = 4, line = 3, family="sans", font=2)
      polygon(density(fm.resid), col = rgb(0, 0, 1, 0.3))
      fit<-fitdistr(fm.resid,"normal")$estimate
      xfit<-seq((min(fm.resid)-1),(max(fm.resid)+1),length=50) 
      yfit<-dnorm(xfit,fit[1],fit[2])
      lines(xfit, yfit, col="red", lwd=2)
      
      ecdf.ks.CI(fm.resid, main= "ECDF of KS Test of Standardized Residuals with Gaussian 95% CI", cex.main=1, font.main=2, font.lab=2, xlab="Stdz. Residuals", ylab="ECDF(Stdz. Residuals)")
      fit<-fitdistr(fm.resid,"normal")$estimate
      Norm3 <- rnorm(10000, fit[1],fit[2])
      lines(ecdf(Norm3), do.points = FALSE, verticals=T, lwd=2, col="blue")
    }
    
   }
    
  
  if (model_type=="nls"){
    
    # overall results
    s = summary(model)
    n = length(model$m$resid())
    RSS <- sum(model$m$resid()^2)
    TSS <- sum((y - mean(y))^2)
    R2 = 1 - (RSS/TSS)
    resid_se = s$sigma
    model_results <- data.frame(matrix(ncol = 0, nrow = 1))
    model_results$formula <- paste(summary(model)$formula[2], summary(model)$formula[1], summary(model)$formula[3])
    model_results$n <- n
    model_results$df <- s$df[2]
    model_results$R2 <- R2
    model_results$RSS <- RSS
    model_results$TSS <- TSS
    model_results$resid_se <- resid_se
    model_results$ConvergToler <- s$convInfo$finTol
    model_results$ConvergIter <- s$convInfo$finIter
    
    Results <- list()
    model_results <- list(model_results)
    Results <- append(Results,model_results)
    names(Results) <- c("NLS Model Results")
    model_resid1 <- nlsResiduals(model)
    fm.resid <- model_resid1$resi2[,2]
    #
    
    # Regular results
    if (coef=="normal" | coef=="both"){
      model_coef <- tidy(model)
      Cov <- vcov(model)
      tt <-qt(c(0.025,0.975),summary(model)$df[2])
      se <- sqrt(diag(Cov))
      ci <-coef(model) + se %o% tt
      model_coef$ci0.025 <- c(exp(ci[1,1]), ci[2,1])
      model_coef$ci0.975 <- c(exp(ci[1,2]), ci[2,2])
      model_coef$EstMethod <- "standard"
      model_coef <- as.data.frame(model_coef)
    }
    
    # sandwich results
    if (coef=="robust" | coef=="both"){
      coeftest <- coeftest(model, vcov = sandwich)
      model_coef_robust <- tidy(coeftest)
      Cov <- sandwich(model)
      tt <-qt(c(0.025,0.975),summary(model)$df[2])
      se <- sqrt(diag(Cov))
      ci <-coef(model) + se %o% tt
      model_coef_robust$ci0.025 <- c(exp(ci[1,1]), ci[2,1])
      model_coef_robust$ci0.975 <- c(exp(ci[1,2]), ci[2,2])
      model_coef_robust$EstMethod <- "robust"
      model_coef_robust <- as.data.frame(model_coef_robust)
    }
    
    if (coef=="normal"){
      model_coef <- list(model_coef)
      nam <- names(Results)
      Results <- append(Results,model_coef)
      names(Results) <- c(nam,"Model Coefficients")
    }
    if (coef=="robust"){
      model_coef_robust <- list(model_coef_robust)
      nam <- names(Results)
      Results <- append(Results,model_coef_robust)
      names(Results) <- c(nam,"Robust Model Coefficients")
    }
    if (coef=="both"){
      model_coef_both = rbind(model_coef,model_coef_robust)
      model_coef_both <- list(model_coef_both)
      nam <- names(Results)
      Results <- append(Results,model_coef_both)
      names(Results) <- c(nam,"Model Coefficients")
    }
    
    
    
    
    shap = shapiro.test(fm.resid)
    shap.t = tidy(shap)
    shap.t$alternative <- NA
    
    fit<-fitdistr(fm.resid,"normal")$estimate
    ksr = ks.test(fm.resid, "pnorm",fit[1],fit[2])
    
    ksr.t = tidy(ksr)
    
    resid_normal_test = rbind(shap.t,ksr.t) 
    resid_normal_test = as.data.frame(resid_normal_test) 
    resid_normal_test <- list(resid_normal_test)
    nam <- names(Results)
    Results <- append(Results,resid_normal_test)
    names(Results) <- c(nam,"Residuals Normality Tests")
    
    
    if (residplots==TRUE){
      model_resid1 <- nlsResiduals(model)
      plot(model_resid1)
    }
    

    
    if (residhists==TRUE){
      Resid1 <- hist(fm.resid, plot=F)
      cts = max(Resid1$density) / max(Resid1$counts)
      freq <- pretty(Resid1$counts)
      freq2 = freq * cts
      pos <- pretty(Resid1$density)
      #freq <- round(pos * length(fm.resid) * with(Resid1, breaks[2] - breaks[1]))
      
      graphics:::plot.histogram(Resid1, freq = F, col="gray70", family="sans", main="Standardized Residuals vs. Gaussian Distribution", cex.main=1, 
                                font.main=2, font.lab=2, xlab = "Residuals", ylab="Frequency",
                                border="black", yaxt='n')#, ylim=c(0,1)
      Axis(side = 2, at = freq2, labels = freq)#
      Axis(side = 4, at = pos, labels = pos)
      mtext("Density           ", side = 4, line = 3, family="sans", font=2)
      polygon(density(fm.resid), col = rgb(0, 0, 1, 0.3))
      fit<-fitdistr(fm.resid,"normal")$estimate
      xfit<-seq((min(fm.resid)-1),(max(fm.resid)+1),length=50) 
      yfit<-dnorm(xfit,fit[1],fit[2])
      lines(xfit, yfit, col="red", lwd=2)
      
      ecdf.ks.CI(fm.resid, main= "ECDF of KS Test of Standardized Residuals with Gaussian 95% CI", cex.main=1, font.main=2, font.lab=2, xlab="Stdz. Residuals", ylab="ECDF(Stdz. Residuals)")
      fit<-fitdistr(fm.resid,"normal")$estimate
      Norm3 <- rnorm(10000, fit[1],fit[2])
      lines(ecdf(Norm3), do.points = FALSE, verticals=T, lwd=2, col="blue")
    }
    
  }
  
  return(Results)
  
}

