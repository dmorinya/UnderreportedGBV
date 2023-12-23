library(R2jags)
library(dplyr)
library(tidyr)
library(doSNOW)

#cl <- makeCluster(5, type="SOCK")
#registerDoSNOW(cl)  

simulate_data <- function(N, q0, lambda, beta, cp, alpha, I, mode) 
{
  x <- rpois(N, lambda+I*beta)
  y <- vector()
  for (i in 1:cp) 
  {
    y[i] <- rbinom(1, x[i], q0)
  }  
  if (mode=="exp")
  {
    q <- 1-(1-q0)*exp(-alpha*seq(1, (N-cp), 1)) #exponential
  }else{
    if (mode=="linear")
    {
      q <- q0+seq(1, round(alpha)-round(cp), 1)/((1/(1-q0))*(round(alpha)-round(cp))) #linear
    }else{
      stop("Wrong mode.")
    }
  }
  for (i in (cp+1):N)
  {
    y[i] <- rbinom(1, x[i], q[i-cp])
  }
  return(y)
}

I <- c(rep(0, 541), rep(1, 14), rep(0, 445)) ### Mandatory confinment indicator

# JAGS to estimate parameters
jags_GBV_sim <- function() {
  for (i in 1:N) {
    lambda0[i] <- ifelse(i < cp, lambda1, lambda2)
    thin[i] <- ifelse(mode==0, q0*step(cp-i)+(1-(1-q0)*exp(-alpha_exp*(i-cp)))*step(i-cp),
                                  q0*step(cp-i)+(q0+(i-cp)/(1/(1-q0)*(alpha_lin-cp)))*step(i-cp))
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
  alpha_lin ~ dunif(N, 5000) #linear
  alpha_exp ~ dgamma(1, 3) #exponential
}

###################################################
############ For parameter variation ##############
###################################################

# Define a list of parameter values
### Uncomment next statements to run the simulation in parallel
gs <- list(q0 = c(.1, .3, .5, .75, .9),
           lambda = c(5, 7, 10),
           beta = c(0.5, 2, 5),
           cp = c(100, 500, 900),
           alpha_lin = c(1200, 1500, 2000),
           alpha_exp = c(0.05, 0.1, 0.5)) %>% 
  expand.grid()

simulation_linear <- function(i)
{  
  print(paste0("simulation_", i))
  #cl2 <- makeCluster(5, type="SOCK")
  #registerDoSNOW(cl2)  
  for (j in 1:100)
  {
    print(paste0("repetition_", j))
    y <- simulate_data(1000, gs$q0[i], gs$lambda[i], gs$beta[i], gs$cp[i], gs$alpha_lin[i], I, "linear") #N, q0, lambda, beta, cp, alpha
    lambda_mean <- gs$lambda[i]
    datos <- list(y=y, N=length(y), Ind=I, lambda_mean=lambda_mean, mode=1)
    init_values <- function(){
      list(q0=runif(1), lambda=rnorm(1, lambda_mean, 0.1), beta=rgamma(1, 2, 1), cp=runif(1, 0, length(y)), 
           alpha_lin=runif(1, length(y), 5000), .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
    }
    params <- c("q0", "lambda", "beta", "cp", "alpha_lin")
    #clusterExport(cl=cl2, list("datos", "lambda_mean", "init_values", "params")) 
    fit_GBV_sim <- jags.parallel(data = datos, inits = init_values, parameters.to.save = params, model.file = jags_GBV_sim,
                                 n.chains = 5, n.iter = 5000, n.burnin = 1000, n.thin = 10, DIC = T)
    lm1_mcmc <- as.data.frame(summary(as.mcmc(fit_GBV_sim))$quantiles)
    
    #Add real parameters
    lm1_mcmc$real_value <- NULL
    lm1_mcmc[1, "real_value"] <- gs[i, 5]
    lm1_mcmc[2, "real_value"] <- gs[i, 3]
    lm1_mcmc[3, "real_value"] <- gs[i, 4]
    lm1_mcmc[5, "real_value"] <- gs[i, 2]
    lm1_mcmc[6, "real_value"] <- gs[i, 1]
    
    sim_name <- paste0("sim_n_", i)
    save(lm1_mcmc, file = paste0("Results/Sim qt linear/", sim_name, "_", j, ".RData"))
    #stopCluster(cl2)
  }
}

simulation_exp <- function(i)
{  
  print(paste0("simulation_", i))
  #cl2 <- makeCluster(5, type="SOCK")
  #registerDoSNOW(cl2)  
  for (j in 1:100)
  {
    print(paste0("repetition_", j))
    y <- simulate_data(1000, gs$q0[i], gs$lambda[i], gs$beta[i], gs$cp[i], gs$alpha_exp[i], I, "exp") #N, q0, lambda, beta, cp, alpha
    lambda_mean <- gs$lambda[i]
    datos <- list(y=y, N=length(y), Ind=I, lambda_mean=lambda_mean, mode=0)
    init_values <- function(){
      list(q0=runif(1), lambda=rnorm(1, lambda_mean, 0.1), beta=rgamma(1, 2, 1), cp=runif(1, 0, length(y)), 
           alpha_lin=runif(1, length(y), 5000), .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
    }
    params <- c("q0", "lambda", "beta", "cp", "alpha_exp")
    #clusterExport(cl=cl2, list("datos", "lambda_mean", "init_values", "params")) 
    fit_GBV_sim <- jags.parallel(data = datos, inits = init_values, parameters.to.save = params, model.file = jags_GBV_sim,
                                 n.chains = 5, n.iter = 5000, n.burnin = 1000, n.thin = 10, DIC = T)
    lm1_mcmc <- as.data.frame(summary(as.mcmc(fit_GBV_sim))$quantiles)
    
    #Add real parameters
    lm1_mcmc$real_value <- NULL
    lm1_mcmc[1, "real_value"] <- gs[i, 5]
    lm1_mcmc[2, "real_value"] <- gs[i, 3]
    lm1_mcmc[3, "real_value"] <- gs[i, 4]
    lm1_mcmc[5, "real_value"] <- gs[i, 2]
    lm1_mcmc[6, "real_value"] <- gs[i, 1]
    
    sim_name <- paste0("sim_n_", i)
    save(lm1_mcmc, file = paste0("Results/Sim qt exp/", sim_name, "_", j, ".RData"))
    #stopCluster(cl2)
  }
}
system.time(foreach(i=1:nrow(gs), combine=rbind, .packages=c("R2jags", "doSNOW")) %do% simulation_exp(i))
#stopCluster(cl)

system.time(foreach(i=1:nrow(gs), combine=rbind, .packages=c("R2jags", "doSNOW")) %do% simulation_linear(i))
#stopCluster(cl)
