#### Model diagnostics
library(readxl)
library(dplyr)
library(MMWRweek)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(R2jags)
load("R/parameters_linear_scen1.RData")
load("R/parameters_exp_scen1.RData")
# load("R/parameters_linear_scen2.RData")
# load("R/parameters_exp_scen2.RData")
# load("R/parameters_linear_scen3.RData")
# load("R/parameters_exp_scen3.RData")

### LINEAR q_t GROWTH
### BASE 25
GBV_mcmc_25 <- as.mcmc(fit_GBV_25)
plot(GBV_mcmc_25)

### BASE 26
GBV_mcmc_26 <- as.mcmc(fit_GBV_26)
plot(GBV_mcmc_26)

### BASE 27
GBV_mcmc_27 <- as.mcmc(fit_GBV_27)
plot(GBV_mcmc_27)

### BASE 31
GBV_mcmc_31 <- as.mcmc(fit_GBV_31)
plot(GBV_mcmc_31)

### BASE 34
GBV_mcmc_34 <- as.mcmc(fit_GBV_34)
plot(GBV_mcmc_34)

### BASE 35
GBV_mcmc_35 <- as.mcmc(fit_GBV_35)
plot(GBV_mcmc_35)

### EXPONENTIAL q_t GROWTH
### BASE 25
GBV_mcmc_25 <- as.mcmc(fit_GBV_25_exp)
plot(GBV_mcmc_25)

### BASE 26
GBV_mcmc_26 <- as.mcmc(fit_GBV_26_exp)
plot(GBV_mcmc_26)

### BASE 27
GBV_mcmc_27 <- as.mcmc(fit_GBV_27_exp)
plot(GBV_mcmc_27)

### BASE 31
GBV_mcmc_31 <- as.mcmc(fit_GBV_31_exp)
plot(GBV_mcmc_31)

### BASE 34
GBV_mcmc_34 <- as.mcmc(fit_GBV_34_exp)
plot(GBV_mcmc_34)

### BASE 35
GBV_mcmc_35 <- as.mcmc(fit_GBV_35_exp)
plot(GBV_mcmc_35)
