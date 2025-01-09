#### Est_Beta.R
#### Rudolf Cesaretti, 7/17/2022

#### "Est_Beta" from radiation models
#### 
#### 
#### 
#### 


Est_Beta <- function(df,
                     d_mat,
                     type, ##"rad_ext_prob_a0.5","rad_ext_prob_a2"
                     crsc=26914
){
  var = "Population.s2"
  xc=df$East
  yc=df$North
  per = df$Period[1]
  pnum = df$PeriodNum[1]
  
  if (type == "rad_orig_prob"){
    rad_orig_prob = Radiation_EstFlows(df=df,var=var,d_mat=d_mat,xcoords=xc,ycoords=yc,
                                       crs_coords=crsc,extended=F,scale="none", 
                                       commuters="none",outputs ="sflines")
    q <- rad_orig_prob[[1]]
  }
  if (type == "rad_ext_prob_a0.5"){
    rad_ext_prob_a0.5 = Radiation_EstFlows(df=df,var=var,d_mat=d_mat,xcoords=xc,ycoords=yc,
                                           crs_coords=crsc, extended=T, scale="none", alpha=0.5,
                                           commuters="none",outputs ="sflines")
    q <- rad_ext_prob_a0.5[[1]]
  }
  if (type == "rad_ext_prob_a2"){
    rad_ext_prob_a2 = Radiation_EstFlows(df=df,var=var,d_mat=d_mat,xcoords=xc,ycoords=yc,
                                         crs_coords=crsc, extended=T, scale="none", alpha=2,
                                         commuters="none",outputs ="sflines")
    q <- rad_ext_prob_a2[[1]]
  }
  
  
  df2 <- q[q$P > 0, ]
  df2$Prob <- df2$P
  df2$Log_Prob <- log(df2$P)
  df2$C_ij <- df2$r_ij
  df2$Log_C_ij <- log(df2$r_ij)
  
  ################
  
  gam_model <- gam(Log_Prob ~ s(Log_C_ij), data = df2, method = "REML")
  gam_model$method <- "GAM REML"
  gam_model$k <- gam_model$rank
  gam_model$sp <- as.numeric(gam_model$sp)
  sum_gam_model <- summary(gam_model)
  power_model <- glm(Log_Prob ~ Log_C_ij, data = df2)
  power_model$k <- NA
  power_model$sp <- NA
  power_model$method <- "log-linear OLS"
  sum_power_model <- summary(power_model)
  exp_model <- glm(Prob ~ C_ij, family = gaussian(link = 'log'), data = df2)
  exp_model$k <- NA
  exp_model$sp <- NA
  exp_model$method <- "log-link OLS"
  sum_exp_model <- summary(exp_model)
  
  glance_custom.lm <- function(x, ...) {
    k <- x$k
    l <- x$sp
    formula <- paste(formula(x)[2],formula(x)[1],formula(x)[3])
    method <- x$method
    family <- x$family$family
    out <- data.frame("k" = k, "lambda" = l, "Formula" = formula, "Method" = method, "Family" = family)
    return(out)
  }
  
  models <- list("Exponential" = exp_model,
                 "Power" = power_model,
                 "GAM" = gam_model)
  gm <- tibble::tribble(
    ~raw,        ~clean,          ~fmt,
    "k",         "k",             0,
    "lambda",    "lambda",       3,
    "nobs",      "n obs.",        0,
    "r.squared", "R2", 3,
    "Formula",   "Formula",       0,
    "Method",    "Method",        0,
    "Family",    "Family",        0)
  #modelsummary(models)
  ms = modelsummary(models, stars = TRUE, gof_map = gm, output='flextable')
  msr=flextable::as_raster(ms)
  # bb = ms$body$dataset
  
  nn=500
  Log_x = seq(min(df2$Log_C_ij),max(df2$Log_C_ij),length.out = nn)
  x = exp(Log_x)
  
  GAM_y = as.numeric(predict.gam(gam_model,data.frame(Log_C_ij=Log_x)))
  GAM_se = (predict.gam(gam_model,data.frame(Log_C_ij=Log_x),se.fit	=T)$se.fit)
  GAM_CIh =  GAM_y + (qnorm(0.975)*GAM_se)
  GAM_CIl =  GAM_y - (qnorm(0.975)*GAM_se)
  
  Power_y = as.numeric(predict(power_model,data.frame(Log_C_ij=Log_x)))
  Power_se = (predict(power_model,data.frame(Log_C_ij=Log_x),se.fit	=T)$se.fit)
  Power_CIh =  Power_y + (qnorm(0.975)*Power_se)
  Power_CIl =  Power_y - (qnorm(0.975)*Power_se)
  
  EXP_y = (as.numeric(predict(exp_model,data.frame(C_ij=x))))
  EXP_se = (predict(exp_model,data.frame(C_ij=x),se.fit	=T)$se.fit)
  EXP_CIh =  EXP_y + (qnorm(0.975)*EXP_se)
  EXP_CIl =  EXP_y - (qnorm(0.975)*EXP_se)
  
  
  pltmod = data.frame(Log_X = c(Log_x,Log_x,Log_x),
                      Model = c(rep("Exponential",nn),rep("Power",nn),rep("GAM",nn)),
                      Predict = c(EXP_y, Power_y, GAM_y),
                      CI95_h = c(EXP_CIh,Power_CIh,GAM_CIh),
                      CI95_l = c(EXP_CIl,Power_CIl,GAM_CIl),
                      Rad = Model = c(rep(type,nn),rep(type,nn),rep(type,nn)))
  
  p2 <- ggplot() + 
    theme_void() + 
    annotation_custom(grid::rasterGrob(msr), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
    theme(plot.margin = margin(0, 0.25, 0, 0, "in"))
  
  p1=ggplot() + geom_point(data=df2, aes(y=Log_Prob, x=Log_C_ij)) +
    geom_line(data=pltmod, aes(x=Log_X, y=Predict, color=Model), size=1.2)+
    geom_ribbon(data=pltmod, aes(x = Log_X, ymin = CI95_l, ymax = CI95_h, fill=Model), alpha = 0.4)+
    scale_x_continuous(labels = math_format(e^.x))+
    scale_y_continuous(labels = math_format(e^.x), limits = c((min(df2$Log_Prob)-0.2), (max(df2$Log_Prob)+0.2)))+
    labs(x="Log Cost-Distance in Hours (C_ij)", y = "Log Flow Conditional Probability")+
    theme_bw()+
    theme(plot.margin = margin(0.25, 0, 0.25, 0.25, "in"),
          axis.title.x = element_text(color="black", size=12, face="bold"),
          axis.title.y = element_text(color="black", size=12, face="bold"),
          axis.text.y = element_text(color="black", size=12),
          axis.text.x = element_text(color="black", size=12),
          #legend.justification=c(0.5,0.5), 
          #legend.position=c(0.25,0.25), 
          legend.position="right", 
          #legend.margin = margin(1, 4, 1, 4),
          #legend.spacing.y = unit(0.1,"line"),legend.spacing.x = unit(0.55,"line"), 
          legend.text=element_text(size=11),legend.title=element_text(size=14)
          #legend.key = element_rect(colour = "transparent", fill = "white"), 
          #legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')
    )
  
  legend <- get_legend(p1)
  p1 <- p1 + theme(legend.position='none')
  
  right.col <- plot_grid(legend, p2, ncol = 1, rel_heights = c(0.2, 0.8))
  #right.col
  plot_row <- plot_grid(p1, right.col, ncol = 2, rel_widths = c(1.25,0.75))
  #plot_row
  title <- ggdraw() + draw_label(paste0("Est Dist Decay Functions for Radiation Model: ",pnum,". ",per),fontface = 'bold',hjust = 0.5, vjust=0.4, size=16 )
  
  out_plt <- plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.04,0.96)) + 
    theme(plot.background = element_rect(fill = "white", colour = NA))
  #out_plt
  nnnn=length(exp_model$residuals)
  
  modsum = data.frame(Period = c(per,per,per),
                      PeriodNum = c(pnum,pnum,pnum),
                      Model = c("Exponential","Power","GAM"),
                      R2 = c(with(summary(exp_model), 1 - deviance/null.deviance),
                             with(summary(power_model), 1 - deviance/null.deviance),
                             sum_gam_model$r.sq),
                      Intercept = c(as.numeric(exp_model$coefficients[1]),
                                    as.numeric(power_model$coefficients[1]), NA),
                      Beta = c(as.numeric(exp_model$coefficients[2]),
                               as.numeric(power_model$coefficients[2]), NA),
                      n = c(nnnn,nnnn,nnnn),
                      Rad = c(type,type,type),
                      Formula = c(paste(formula(exp_model)[2],formula(exp_model)[1],formula(exp_model)[3]),
                                  paste(formula(power_model)[2],formula(power_model)[1],formula(power_model)[3]),
                                  paste(formula(power_model)[2],formula(power_model)[1],formula(power_model)[3])))
  
  out_list = list()
  
  out_list[[1]] <- out_plt
  out_list[[2]] <- modsum
  out_list[[3]] <- pltmod
  
  return(out_list)
}