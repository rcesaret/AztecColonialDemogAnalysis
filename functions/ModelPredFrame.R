#### ModelPredFrame.R
#### Rudolf Cesaretti, 11/9/2022

#### "ModelPredFrame.R"
#### 
#### 
#### 
#### 

pak <- c("nls.multstart", "nlme", "nlstools")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)

###############################################################
######################  ModelPredFrame  #######################
###############################################################


ModelPredFrame <- function(object, 
                           data,
                           x, 
                           y,
                           type = 'predict',
                           transformation = "none",#c("logx","logy", "loglog"),
                           logx = F,
                           logy = F,
                           loglog = F,
                           unlog_after = F,
                           unlog_after_y = F,
                           unlog_after_x = F,
                           xmin = NULL,
                           xmax = NULL,
                           n=1024){

  if (is.null(xmin)){xmin=min(data[,x],na.rm=T)}
  if (is.null(xmax)){xmax=max(data[,x],na.rm=T)}

  q = as.data.frame(seq(from=xmin, to=xmax, length.out = n))
  names(q) = x
  
  if (type == 'predict'){
    z = predict(object, newdata = q)
  }
  
  if (type == 'coef'){
    
    if (transformation == "none"){
      z = (object[2] * q[,1]) + object[1]
    }
    
    #provide log-linear model coefs, return unlogged predict
    if (transformation == "loglog"){
      z = (q[,1]^object[2]) * exp(object[1])
    }
    
    #provide log(x) model coefs, return unlogged predict
    if (transformation == "logx"){
      z = (object[2] * q[,1]) + object[1]
      q[,1] = exp(q[,1])
    }
    
    #provide log(y) model coefs, return unlogged predict
    if (transformation == "logy"){
      z = exp((object[2] * q[,1]) + object[1])
    }
    
  }
  if (unlog_after == T){
    z = exp(z)
    q[,1] <- exp(q[,1])
  }
  if (unlog_after_x == T){
    q[,1] <- exp(q[,1])
  }
  if (unlog_after_y == T){
    z = exp(z)
  }
  
  q$z <- z
  names(q) = c(x,paste0(y,".p"))
  
  return(q)
}
