library(dplyr)
library(tidyr)
library(tidyverse)
library(tibble)
library(R2jags)

### Files in the zip files should be uncompressed first.

### LINEAR q_t

## Load files
file_names <- list.files(path = "Results/Sim qt linear", pattern = "*.RData", full.names = T)

#List mcmc objects
all_data <- mcmc.list(lapply(file_names, function(dat) mcmc(get(load(dat)))))

# Join all data
df<- do.call(rbind.data.frame, all_data) 
df$parameter <- row.names(df)
rownames(df)<- NULL
df$parameter <- gsub("[^A-Za-z]+", "", df$parameter)

#Filter "deviance"
df<- subset(df,parameter!='deviance' )

#Create fail interval columns
df$fail_left <-  ifelse(df$real_value > df$`2.5%`,0,1) 
df$fail_right<-  ifelse(df$real_value < df$`97.5%`, 0, 1)


## Results summary

results<- df %>% group_by(parameter) %>% summarise(
  n = n(),
  tot_fail_left = sum(fail_left),
  tot_fail_right = sum(fail_right), 
  perc_left = 100 * tot_fail_left/n,
  perc_right = 100 * tot_fail_right/n,
  coverage = 100 - perc_left - perc_right, 
  mean_long = mean(`97.5%` - `2.5%`),
  min_long = min(`97.5%` - `2.5%`),
  max_long = max(`97.5%` - `2.5%`),
  mean_bias= 100 * mean(abs(real_value-`50%`))/n
)

results

df %>% group_by(parameter) %>% summarise(unique(real_value))

### EXPONENTIAL q_t

## Load files
file_names <- list.files(path = "Results/Sim qt exp", pattern = "*.RData", full.names = T)

#List mcmc objects
all_data <- mcmc.list(lapply(file_names, function(dat) mcmc(get(load(dat)))))

# Join all data
df<- do.call(rbind.data.frame, all_data) 
df$parameter <- row.names(df)
rownames(df)<- NULL
df$parameter <- gsub("[^A-Za-z]+", "", df$parameter)

#Filter "deviance"
df<- subset(df,parameter!='deviance' )

#Create fail interval columns
df$fail_left <-  ifelse(df$real_value > df$`2.5%`,0,1) 
df$fail_right<-  ifelse(df$real_value < df$`97.5%`, 0, 1)


## Results summary

results<- df %>% group_by(parameter) %>% summarise(
  n = n(),
  tot_fail_left = sum(fail_left),
  tot_fail_right = sum(fail_right), 
  perc_left = 100 * tot_fail_left/n,
  perc_right = 100 * tot_fail_right/n,
  coverage = 100 - perc_left - perc_right, 
  mean_long = mean(`97.5%` - `2.5%`),
  min_long = min(`97.5%` - `2.5%`),
  max_long = max(`97.5%` - `2.5%`),
  mean_bias= 100 * mean(abs(real_value-`50%`))/n
)

results

df %>% group_by(parameter) %>% summarise(unique(real_value))
