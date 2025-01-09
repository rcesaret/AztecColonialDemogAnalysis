#### KDProb.R
#### Rudolf Cesaretti, 6/1/2022

#### "KDProb" is a function that calculates the probability under the curve for 
####  any probability distribution (e.g. empirical PDFs, non-theoretical posteriors, 
####  kernel density dists, etc.) via numerical approximation of the PDF integral

pak <- c("tidyverse", "era", "zoo", "scales")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)

#############################################################################
################################  KDProb  ###################################
#############################################################################


KDProb <- function(A, #a vector of x-axis observations
                   A.wt = NULL, #a vector of weights the same length as A
                   B, #a vector of x-axis observations
                   B.wt = NULL, #a vector of weights the same length as B
                   cuts = NULL, # vector of cutpoints for probability estimates; must be supplied if output = "prob" or "both"
                   xlims = NULL, # a vector length 2 with the x-axis extrema to be used
                   nscale = 1000000, #precision in PDF estimates
                   output = c("prob","dens","both"), #desired outputs
                   graph=F, #whether to return a graph of the results
                   lab = NA,
                   mov.avg=F,
                   mov.avg.win = 10,
                   mov.avg.rep = 10
){
  
  # find the extrema of the two distributions to provide the left and right-most 
  # points of the grid at which the density is to be estimated
  if (is.null(xlims)) {
    PDFA = density(A, weights = A.wt)
    PDFB = density(B, weights = B.wt)
  } else {
    PDFA = density(A, weights = A.wt, from=xlims[1], to=xlims[2])
    PDFB = density(B, weights = B.wt, from=xlims[1], to=xlims[2])
  }
  
  extrema <- c(max(PDFA$x),max(PDFB$x), min(PDFA$x), min(PDFB$x))
  t0 = min(extrema)
  t1 = max(extrema)
  
  
  #compute kernel density estimates for the two distributions fitted to the same min and max; 
  # both now have the same set of 512 x values (years/time)
  PDFA = density(A, weights = A.wt, from=t0, to=t1)
  PDFB = density(B, weights = B.wt, from=t0, to=t1)
  
  #save this range as its own vector
  X = PDFA$x
  
  #save pdfs of the two priors
  densA <- PDFA$y
  densB <- PDFB$y
  
  densAB = densA*densB
  
  if (mov.avg == TRUE){
    for (i in 1:mov.avg.rep){
      densAB = rollapply(densAB, width = mov.avg.win, function(...) {mean(...)}, partial = TRUE)
    }
  }
  
  f <- function(i)splinefun(X, densAB)(i) 
  
  int=integrate(f,t0,t1)
  
  densAB = densAB/int$value
  
  if (output == "dens" | output == "both"){
    
    df = data.frame(x=X,y=round(densAB*nscale*10))
    df = df %>% slice(rep(seq_len(n()), y))
    PDFAB <- density(df$x, from = t0, to = t1, n = 1024)
    
    if (output == "both"){
      out <- list()
      out[[1]] <- PDFAB
    }
    if (output == "dens"){
      out <- PDFAB
    }
  }
  
  if (output == "prob" | output == "both"){
    tmp = rep(NA,(length(cuts)-1))
    Pr <- data.frame(lab = tmp, Xmin = tmp, Xmax = tmp, Prob1 = tmp, Prob2 = tmp, Diff = tmp, ProbMean = tmp)
    
    for (i in 1:(length(cuts)-1)){
      #probability estimate #1
      est1 <- sum(densAB[X > cuts[i] & X < cuts[i+1]])/sum(densAB)
      
      #probability estimate #2
      df = data.frame(x=X,y=round(densAB*nscale))
      df = df %>% slice(rep(seq_len(n()), y))
      dens <- density(df$x, from = cuts[i], to = cuts[i+1], n = 1024)
      est2 <- with(dens, sum(y * diff(x)[1]))
      Pr[i,1] <- lab
      Pr[i,2] <- cuts[i]
      Pr[i,3] <- cuts[i+1]
      Pr[i,4] <- est1
      Pr[i,5] <- est2
      Pr[i,6] <- est1 - est2
      Pr[i,7] <- mean(est1,est2,na.rm=T)
    }
    if (output == "both"){
      out[[2]] <- Pr
    }
    if (output == "prob"){
      out <- Pr
    }
  }
  
  if (graph == TRUE){
    plot(X,densAB,type="l", col="red", lwd=2, ylab='PDF',xlab='Years BC/AD')
    lines(X, densA, type = "l", col = "green", lwd=1.6)
    lines(X, densB , type = "l", col = "blue", lwd=1.6)
    legend('top',bty='n',cex=1.75,inset=-0.2,ncol=3,c('A','B','AB'),col=c('green','blue','red'),lty=c(1,1,1),lwd=c(3,3,5),xpd=T)
  }
  
  return(out)
}
