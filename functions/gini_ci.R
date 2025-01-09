#### gini_ci.R
#### Rudolf Cesaretti, 12/10/2023

#### "gini_ci" 
#### 
#### 
#### 
#### 

pak <- c("REAT", "boot")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)

###############################################################
##########################  gini_ci  ##########################
###############################################################
#Data = Tecalli_Mace_per_Pilli$NumMacehuales
#n_samples = 1000
#replacement = TRUE
#ci_level = 0.95
#Name = "Series1"
#na_rm = TRUE
#Weights = NULL
#CoefNorm = FALSE

gini_ci <- function(Data, # a numeric vector
                    n_samples = 1000, # integer
                    replacement = TRUE, # boolean; sampling with or without replacement
                    ci_level = 0.95, # confidence level to report; values between (0,1)
                    Name = "Series1",
                    na_rm = TRUE, #logical argument that indicates whether NA values should be excluded before computing results
                    Weights = NULL, #A numeric vector containing the weighting data (e.g. size of income classes when calculating a Gini coefficient for aggregated income data)
                    CoefNorm = FALSE # logical argument that indicates if the function output is the non-standardized or the standardized Gini coefficient (default: coefnorm = FALSE, that means the non-standardized Gini coefficient is returned)
                    ) {
  # Function to compute Gini index using the REAT::gini function
  compute_gini <- function(data, indices, CN = CoefNorm, W = Weights, NA_RM = na_rm) {
    d <- data[indices] #allows boot to select sample
    gini_value <- REAT::gini(x = d, lc = FALSE, na.rm = NA_RM, coefnorm = CN, weighting = W)
    return(gini_value)
  }
  
  # Set the seed for reproducibility
  #set.seed(123)
  
  # Perform bootstrap resampling
  #boot_result <- boot::boot(data = Data, statistic = compute_gini, R = n_samples, sim = "ordinary", 
 #                     ran.gen = NULL, mle = NULL, stype = "i", strata = NULL,
  #                    weights = NULL, subset = NULL, na.rm = na_rm, parallel = "no")
  
  boot_result <- boot::boot(data = Data, statistic = compute_gini, R = n_samples, CN = CoefNorm, W = Weights, NA_RM = na_rm)
  
  estimate = boot_result[["t0"]]
  #bias:
  Bias = mean(boot_result$t)-boot_result$t0
  #se: 
  SE = sd(boot_result$t)
  # number of observations
  
  n = length(boot_result[["data"]])
  # Calculate confidence intervals
  ci <- boot::boot.ci(boot_result, type = "perc", conf = ci_level)
  
  # Extract upper and lower bounds of the confidence interval
  lower_bound <- ci[["percent"]][4]
  upper_bound <- ci[["percent"]][5]
  
  name_low = paste0("CI_",ci_level,"_Low")
  name_high = paste0("CI_",ci_level,"_High")
  
  Out = data.frame(Name, n, n_samples, estimate, lower_bound, upper_bound, Bias, SE)
  colnames(Out) = c("Name", "n_obs", "n_samples", "Gini_obs", name_low, name_high, "Bias", "SE")
  
  return(Out)
}
















