library(lubridate)
library(R2jags)

### Uncomment the following lines according to the scenario to be estimated
load("R/parameters_linear_scen1.RData")
load("R/parameters_exp_scen1.RData")
#load("R/parameters_linear_scen2.RData")
#load("R/parameters_exp_scen2.RData")
#load("R/parameters_linear_scen3.RData")
#load("R/parameters_exp_scen3.RData")

### Estimate date of total reporting
estimate_date <- function(q0, alpha, cp, prob=0.9999, mode)
{
  q <- vector()
  cp <- round(cp)
  if (mode == "exp")
  {
    for (i in (cp+1):10000)
    {
      q[i-cp] <- 1-(1-q0)*exp(-alpha*(i-cp))
    }
    x <- which(q>prob)[1]
    if (is.na(x)) stop("Probability not reached in the considered period.")
    est_date <- as.Date("2010-01-01") %m+% weeks(cp) %m+% weeks(ceiling(x))
  }else{
    if (mode == "linear")
    {
      est_date <- as.Date("2010-01-01") %m+% weeks(ceiling(alpha))
    }else{
      stop("Wrong mode.")
    }
  }
  return(est_date)
}

### LINEAR MODE ###
### BASE 25
GBV_mcmc_25 <- as.mcmc(fit_GBV_25)
est_date_25 <- estimate_date(q0=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="q0"],
                             alpha=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="alpha"],
                             cp=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"], 
                             mode="linear")

### BASE 26
GBV_mcmc_26 <- as.mcmc(fit_GBV_26)
est_date_26 <- estimate_date(q0=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="q0"],
                             alpha=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="alpha"],
                             cp=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"], 
                             mode="linear")

### BASE 27
GBV_mcmc_27 <- as.mcmc(fit_GBV_27)
est_date_27 <- estimate_date(q0=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="q0"],
                             alpha=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="alpha"],
                             cp=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"], 
                             mode="linear")

### BASE 31
GBV_mcmc_31 <- as.mcmc(fit_GBV_31)
est_date_31 <- estimate_date(q0=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="q0"],
                             alpha=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="alpha"],
                             cp=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"], 
                             mode="linear")

### BASE 34
GBV_mcmc_34 <- as.mcmc(fit_GBV_34)
est_date_34 <- estimate_date(q0=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="q0"],
                             alpha=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="alpha"],
                             cp=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"], 
                             mode="linear")

### BASE 35
GBV_mcmc_35 <- as.mcmc(fit_GBV_35)
est_date_35 <- estimate_date(q0=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="q0"],
                             alpha=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="alpha"],
                             cp=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"], 
                             mode="linear")

### EXPONENTIAL MODE ###
### BASE 25
GBV_mcmc_25 <- as.mcmc(fit_GBV_25_exp)
est_date_25 <- estimate_date(q0=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="q0"],
                             alpha=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="alpha"],
                             cp=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"], 
                             mode="exp", prob=0.99)

### BASE 26
GBV_mcmc_26 <- as.mcmc(fit_GBV_26_exp)
est_date_26 <- estimate_date(q0=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="q0"],
                             alpha=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="alpha"],
                             cp=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"], 
                             mode="exp", prob=0.99)

### BASE 27
GBV_mcmc_27 <- as.mcmc(fit_GBV_27_exp)
est_date_27 <- estimate_date(q0=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="q0"],
                             alpha=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="alpha"],
                             cp=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"], 
                             mode="exp", prob=0.99)

### BASE 31
GBV_mcmc_31 <- as.mcmc(fit_GBV_31_exp)
est_date_31 <- estimate_date(q0=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="q0"],
                             alpha=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="alpha"],
                             cp=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"], 
                             mode="exp", prob=0.99)

### BASE 34
GBV_mcmc_34 <- as.mcmc(fit_GBV_34_exp)
est_date_34 <- estimate_date(q0=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="q0"],
                             alpha=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="alpha"],
                             cp=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"], 
                             mode="exp", prob=0.99)

### BASE 35
GBV_mcmc_35 <- as.mcmc(fit_GBV_35_exp)
est_date_35 <- estimate_date(q0=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="q0"],
                             alpha=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="alpha"],
                             cp=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"], 
                             mode="exp", prob=0.99)
