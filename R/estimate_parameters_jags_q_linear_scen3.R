library(readxl)
library(dplyr)
library(MMWRweek)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(R2jags)
library(doParallel)

ncores <- detectCores() - 1  
registerDoParallel(cores=ncores)  
cl <- makeCluster(ncores, type="FORK")  

load("R/prop_survey_scen3.RData")
seguiments     <- read_xls("Data/3554569_DIAGNOSTICS.xls")
dades <- seguiments
dades$PR_DDE <- as.Date(dades$PR_DDE, format="%Y-%m-%d")
dades <- dades[dades$SEXE=="D", ]

### DATA PREPROCESSING
dades$Week <- MMWRweek(dades$PR_DDE)[, 2]
dades$Year <- as.numeric(substr(dades$PR_DDE, 1, 4))
dades$Month <- as.numeric(substr(dades$PR_DDE, 6, 7))
dades$cas <- 1
dades$AGE <- year(dades$PR_DDE)-year(dades$DAT_NAIX)
dades <- dades[dades$AGE >= 16, ] ### Only women 16 years old or older

### Only confirmed cases
dades <- dades[grepl("T74.", dades$PR_COD_PS), ]

data_ts_Week <- dades %>% group_by(BASE, Year, Week) %>% summarise(Cases=sum(cas))
data_ts_Week <- data_ts_Week %>% complete(Week=1:53, fill=list(Cases=0))

### Indicator for covid-19 confinement
I <- c(rep(0, 541), rep(1, 14), rep(0, 81))

# JAGS to estimate parameters
jags_GBV <- function() {
  for (i in 1:N) {
    lambda0[i] <- ifelse(i < cp, lambda1, lambda2)
    thin[i] <- q0*step(cp-i)+(q0+(i-cp)/(1/(1-q0)*(alpha-cp)))*step(i-cp) #linear
    z[i] ~ dpois(lambda0[i])
    x[i] ~ dpois(lambda + Ind[i]*beta)
    y[i] ~ dpois(thin[i]*x[i])
  }  
  # priors
  lambda1 ~ dgamma(1, 1)
  lambda2 ~ dgamma(3, 3)
  q0 ~ dbeta(1, 1)
  lambda ~ dnorm(lambda_mean, 100)
  beta ~ dgamma(2, 1)
  cp ~ dunif(0, N)
  alpha ~ dunif(N, 5000) #linear
}

#### BASE 25
lambda_mean <- prop_survey[names(prop_survey)=="BASE 25"]
datos <- list(y=data_ts_Week[data_ts_Week$BASE==25, ]$Cases, z=data_ts_Week[data_ts_Week$BASE==25, ]$Cases, 
              N=length(data_ts_Week[data_ts_Week$BASE==25, ]$Cases), Ind=I, lambda_mean=lambda_mean)
init_values <- function(){
  list(q0=runif(1), lambda=rnorm(1,lambda_mean,0.1), beta=runif(1,2,1), cp=runif(1,0,length(data_ts_Week[data_ts_Week$BASE==25, ]$Cases)),
       alpha=runif(1, length(data_ts_Week[data_ts_Week$BASE==25, ]$Cases), 5000), #linear
       .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
}
params <- c("q0", "lambda", "beta", "cp", "alpha")
clusterExport(cl, list("datos", "params", "init_values"))
fit_GBV_25 <- jags.parallel(data = datos, inits = NULL, parameters.to.save = params, model.file = jags_GBV,
                            n.chains = 5, n.iter = 50000, n.burnin = 5000, n.thin = 10, DIC = T)

#### BASE 26
lambda_mean <- prop_survey[names(prop_survey)=="BASE 26"]
datos <- list(y=data_ts_Week[data_ts_Week$BASE==26, ]$Cases, z=data_ts_Week[data_ts_Week$BASE==26, ]$Cases, 
              N=length(data_ts_Week[data_ts_Week$BASE==26, ]$Cases), Ind=I, lambda_mean=lambda_mean)
init_values <- function(){
  list(q0=runif(1), lambda=rnorm(1,lambda_mean,0.1), beta=runif(1,2,1), cp=runif(1,0,length(data_ts_Week[data_ts_Week$BASE==26, ]$Cases)),
       alpha=runif(1, length(data_ts_Week[data_ts_Week$BASE==26, ]$Cases), 5000), #linear
       .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
}
params <- c("q0", "lambda", "beta", "cp", "alpha")
clusterExport(cl, list("datos", "params", "init_values"))
fit_GBV_26 <- jags.parallel(data = datos, inits = NULL, parameters.to.save = params, model.file = jags_GBV,
                            n.chains = 5, n.iter = 50000, n.burnin = 5000, n.thin = 10, DIC = T)

#### BASE 27
lambda_mean <- prop_survey[names(prop_survey)=="BASE 27"]
datos <- list(y=data_ts_Week[data_ts_Week$BASE==27, ]$Cases, z=data_ts_Week[data_ts_Week$BASE==27, ]$Cases, 
              N=length(data_ts_Week[data_ts_Week$BASE==27, ]$Cases), Ind=I, lambda_mean=lambda_mean)
init_values <- function(){
  list(q0=runif(1), lambda=rnorm(1,lambda_mean,0.1), beta=runif(1,2,1), cp=runif(1,0,length(data_ts_Week[data_ts_Week$BASE==27, ]$Cases)),
       alpha=runif(1, length(data_ts_Week[data_ts_Week$BASE==27, ]$Cases), 5000), #linear
       .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
}
params <- c("q0", "lambda", "beta", "cp", "alpha")
clusterExport(cl, list("datos", "params", "init_values"))
fit_GBV_27 <- jags.parallel(data = datos, inits = NULL, parameters.to.save = params, model.file = jags_GBV,
                            n.chains = 5, n.iter = 50000, n.burnin = 5000, n.thin = 10, DIC = T)

#### BASE 31
lambda_mean <- prop_survey[names(prop_survey)=="BASE 31"]
datos <- list(y=data_ts_Week[data_ts_Week$BASE==31, ]$Cases, z=data_ts_Week[data_ts_Week$BASE==31, ]$Cases, 
              N=length(data_ts_Week[data_ts_Week$BASE==31, ]$Cases), Ind=I, lambda_mean=lambda_mean)
init_values <- function(){
  list(q0=runif(1), lambda=rnorm(1,lambda_mean,0.1), beta=runif(1,2,1), cp=runif(1,0,length(data_ts_Week[data_ts_Week$BASE==31, ]$Cases)),
       alpha=runif(1, length(data_ts_Week[data_ts_Week$BASE==31, ]$Cases), 5000), #linear
       .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
}
params <- c("q0", "lambda", "beta", "cp", "alpha")
clusterExport(cl, list("datos", "params", "init_values"))
fit_GBV_31 <- jags.parallel(data = datos, inits = NULL, parameters.to.save = params, model.file = jags_GBV,
                            n.chains = 5, n.iter = 50000, n.burnin = 5000, n.thin = 10, DIC = T)

#### BASE 34
lambda_mean <- prop_survey[names(prop_survey)=="BASE 34"]
datos <- list(y=data_ts_Week[data_ts_Week$BASE==34, ]$Cases, z=data_ts_Week[data_ts_Week$BASE==34, ]$Cases, 
              N=length(data_ts_Week[data_ts_Week$BASE==34, ]$Cases), Ind=I, lambda_mean=lambda_mean)
init_values <- function(){
  list(q0=runif(1), lambda=rnorm(1,lambda_mean,0.1), beta=runif(1,2,1), cp=runif(1,0,length(data_ts_Week[data_ts_Week$BASE==34, ]$Cases)),
       alpha=runif(1, length(data_ts_Week[data_ts_Week$BASE==34, ]$Cases), 5000), #linear
       .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
}
params <- c("q0", "lambda", "beta", "cp", "alpha")
clusterExport(cl, list("datos", "params", "init_values"))
fit_GBV_34 <- jags.parallel(data = datos, inits = NULL, parameters.to.save = params, model.file = jags_GBV,
                            n.chains = 5, n.iter = 50000, n.burnin = 5000, n.thin = 10, DIC = T)

#### BASE 35
lambda_mean <- prop_survey[names(prop_survey)=="BASE 35"]
datos <- list(y=data_ts_Week[data_ts_Week$BASE==35, ]$Cases, z=data_ts_Week[data_ts_Week$BASE==35, ]$Cases, 
              N=length(data_ts_Week[data_ts_Week$BASE==35, ]$Cases), Ind=I, lambda_mean=lambda_mean)
init_values <- function(){
  list(q0=runif(1), lambda=rnorm(1,lambda_mean,0.1), beta=runif(1,2,1), cp=runif(1,0,length(data_ts_Week[data_ts_Week$BASE==35, ]$Cases)),
       alpha=runif(1, length(data_ts_Week[data_ts_Week$BASE==35, ]$Cases), 5000), #linear
       .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
}
params <- c("q0", "lambda", "beta", "cp", "alpha")
clusterExport(cl, list("datos", "params", "init_values"))
fit_GBV_35 <- jags.parallel(data = datos, inits = NULL, parameters.to.save = params, model.file = jags_GBV,
                            n.chains = 5, n.iter = 50000, n.burnin = 5000, n.thin = 10, DIC = T)

save(list=c("fit_GBV_25", "fit_GBV_26", "fit_GBV_27", "fit_GBV_31", "fit_GBV_34", "fit_GBV_35"),
     file="R/parameters_linear_scen3.RData")
