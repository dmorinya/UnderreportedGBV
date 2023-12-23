source("R/reconstruct.R")

#### GLOBAL PERIOD
sum(data_ts_Week$Cases[data_ts_Week$BASE==25])+sum(data_ts_Week$Cases[data_ts_Week$BASE==26])+
  sum(data_ts_Week$Cases[data_ts_Week$BASE==27])+sum(data_ts_Week$Cases[data_ts_Week$BASE==31])+
  sum(data_ts_Week$Cases[data_ts_Week$BASE==34])+sum(data_ts_Week$Cases[data_ts_Week$BASE==35])


sum(rec_25)+sum(rec_26)+sum(rec_27)+sum(rec_31)+sum(rec_34)+sum(rec_35)

#### PRE TRAINING
sum(data_ts_Week$Cases[data_ts_Week$BASE==25][1:summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"]])+
  sum(data_ts_Week$Cases[data_ts_Week$BASE==26][1:summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"]])+
  sum(data_ts_Week$Cases[data_ts_Week$BASE==27][1:summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"]])+
  sum(data_ts_Week$Cases[data_ts_Week$BASE==31][1:summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"]])+
  sum(data_ts_Week$Cases[data_ts_Week$BASE==34][1:summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"]])+
  sum(data_ts_Week$Cases[data_ts_Week$BASE==35][1:summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"]])

sum(rec_25[1:summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"]])+
  sum(rec_26[1:summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"]])+
  sum(rec_27[1:summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"]])+
  sum(rec_31[1:summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"]])+
  sum(rec_34[1:summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"]])+
  sum(rec_35[1:summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"]])

#### POST TRAINING
sum(data_ts_Week$Cases[data_ts_Week$BASE==25][(summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"]+1):length(data_ts_Week$Cases[data_ts_Week$BASE==25])])+
  sum(data_ts_Week$Cases[data_ts_Week$BASE==26][(summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"]+1):length(data_ts_Week$Cases[data_ts_Week$BASE==26])])+
  sum(data_ts_Week$Cases[data_ts_Week$BASE==27][(summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"]+1):length(data_ts_Week$Cases[data_ts_Week$BASE==27])])+
  sum(data_ts_Week$Cases[data_ts_Week$BASE==31][(summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"]+1):length(data_ts_Week$Cases[data_ts_Week$BASE==31])])+
  sum(data_ts_Week$Cases[data_ts_Week$BASE==34][(summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"]+1):length(data_ts_Week$Cases[data_ts_Week$BASE==34])])+
  sum(data_ts_Week$Cases[data_ts_Week$BASE==35][(summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"]+1):length(data_ts_Week$Cases[data_ts_Week$BASE==35])])
  
sum(rec_25[(summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"]+1):length(rec_25)])+
  sum(rec_26[(summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"]+1):length(rec_26)])+
  sum(rec_27[(summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"]+1):length(rec_27)])+
  sum(rec_31[(summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"]+1):length(rec_31)])+
  sum(rec_34[(summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"]+1):length(rec_25)])+
  sum(rec_35[(summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"]+1):length(rec_35)])
  