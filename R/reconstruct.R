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

x_rec <- function(par, y, I, mode) ## par[1]=q0; par[2]=lambda; par[3]=beta; par[4]=alpha; par[5]=cp
{
  n <- length(y)
  x <- vector()
  if (mode=="exp")
  {
    q <- c(rep(par[1], round(par[5])), 1-(1-par[1])*exp(-par[4]*seq(1, n-round(par[5]), 1))) #exponential
  }else{
    if (mode=="linear")
    {
      q <- c(rep(par[1], round(par[5])), par[1]+seq(1, round(par[4])-round(par[5]), 1)/((1/(1-par[1]))*(round(par[4])-round(par[5])))) #linear
    }else{
      stop("Wrong mode.")
    }
  }
  for (i in 1:n)
  {
    prob_ant <- 0
    prob <- dbinom(y[i], y[i], q[i])*dpois(y[i], par[2]+I[i]*par[3])
    x[i] <- y[i]
    r <- y[i]+1
    while (r < 300 & prob > prob_ant)
    {
      prob_ant <- prob
      prob <- dbinom(y[i], r, q[i])*dpois(r, par[2]+I[i]*par[3])
      if (prob > prob_ant) x[i] <- r
      r <- r+1
    }
  }
  return(x)
}

### LINEAR q_t GROWTH (FIGURE SCALES SHOULD BE REVISED ACCORDING TO EACH SCENARIO)
### BASE 25
GBV_mcmc_25 <- as.mcmc(fit_GBV_25)
rec_25 <- x_rec(par=c(q0=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="q0"],
                      lambda=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="lambda"],
                      beta=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="beta"],
                      alpha=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="alpha"],
                      cp=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"]),
                y=data_ts_Week[data_ts_Week$BASE==25, ]$Cases,
                I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="linear")

df1 <- data_ts_Week[data_ts_Week$BASE==25, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_25))
df2 <- data.frame(BASE=rep("25", length(rec_25)), Week=rep(seq(1:53), length(rec_25)/53), Cases=rec_25, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_25)),
                  Series=rep("Reconstructed", length(rec_25)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p1 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() + ylab("") +
  ylim(0, 30)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") + xlab("") + ggtitle("SUBAREA A")

### BASE 26
GBV_mcmc_26 <- as.mcmc(fit_GBV_26)
rec_26 <- x_rec(par=c(q0=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="q0"],
                      lambda=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="lambda"],
                      beta=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="beta"],
                      alpha=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="alpha"],
                      cp=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"]),
                y=data_ts_Week[data_ts_Week$BASE==26, ]$Cases,
                I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="linear")

df1 <- data_ts_Week[data_ts_Week$BASE==26, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_26))
df2 <- data.frame(BASE=rep("26", length(rec_26)), Week=rep(seq(1:53), length(rec_26)/53), Cases=rec_26, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_26)),
                  Series=rep("Reconstructed", length(rec_26)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p2 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() + ylab("") +
  ylim(0, 30)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") + xlab("") + ggtitle("SUBAREA B")

### BASE 27
GBV_mcmc_27 <- as.mcmc(fit_GBV_27)
rec_27 <- x_rec(par=c(q0=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="q0"],
                      lambda=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="lambda"],
                      beta=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="beta"],
                      alpha=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="alpha"],
                      cp=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"]),
                y=data_ts_Week[data_ts_Week$BASE==27, ]$Cases,
                I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="linear")

df1 <- data_ts_Week[data_ts_Week$BASE==27, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_27))
df2 <- data.frame(BASE=rep("27", length(rec_27)), Week=rep(seq(1:53), length(rec_27)/53), Cases=rec_27, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_27)),
                  Series=rep("Reconstructed", length(rec_27)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p3 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() + ylab("") +
  ylim(0, 30)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") + xlab("") + ggtitle("SUBAREA C")

### BASE 31
GBV_mcmc_31 <- as.mcmc(fit_GBV_31)
rec_31 <- x_rec(par=c(q0=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="q0"],
                      lambda=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="lambda"],
                      beta=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="beta"],
                      alpha=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="alpha"],
                      cp=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"]),
                y=data_ts_Week[data_ts_Week$BASE==31, ]$Cases,
                I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="linear")

df1 <- data_ts_Week[data_ts_Week$BASE==31, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_31))
df2 <- data.frame(BASE=rep("31", length(rec_31)), Week=rep(seq(1:53), length(rec_31)/53), Cases=rec_31, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_31)),
                  Series=rep("Reconstructed", length(rec_31)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p4 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() + ylab("") +
  ylim(0, 55)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") + xlab("") + ggtitle("SUBAREA D")

### BASE 34
GBV_mcmc_34 <- as.mcmc(fit_GBV_34)
rec_34 <- x_rec(par=c(q0=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="q0"],
                      lambda=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="lambda"],
                      beta=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="beta"],
                      alpha=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="alpha"],
                      cp=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"]),
                y=data_ts_Week[data_ts_Week$BASE==34, ]$Cases,
                I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="linear")

df1 <- data_ts_Week[data_ts_Week$BASE==34, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_34))
df2 <- data.frame(BASE=rep("34", length(rec_34)), Week=rep(seq(1:53), length(rec_34)/53), Cases=rec_34, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_34)),
                  Series=rep("Reconstructed", length(rec_34)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p5 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line()  + ylab("") +
  ylim(0, 15)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") + xlab("") + ggtitle("SUBAREA E")

### BASE 35
GBV_mcmc_35 <- as.mcmc(fit_GBV_35)
rec_35 <- x_rec(par=c(q0=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="q0"],
                      lambda=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="lambda"],
                      beta=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="beta"],
                      alpha=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="alpha"],
                      cp=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"]),
                y=data_ts_Week[data_ts_Week$BASE==35, ]$Cases,
                I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="linear")

df1 <- data_ts_Week[data_ts_Week$BASE==35, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_35))
df2 <- data.frame(BASE=rep("35", length(rec_35)), Week=rep(seq(1:53), length(rec_35)/53), Cases=rec_35, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_35)),
                  Series=rep("Reconstructed", length(rec_35)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p6 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() + ylab("") + 
  ylim(0, 55)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") + xlab("") + ggtitle("SUBAREA F")

ggarrange(p1, p2, p3, p4, p5, p6,
          ncol = 3, nrow = 2)

### EXPONENTIAL q_t GROWTH
### BASE 25
GBV_mcmc_25 <- as.mcmc(fit_GBV_25_exp)
rec_25 <- x_rec(par=c(q0=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="q0"],
                      lambda=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="lambda"],
                      beta=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="beta"],
                      alpha=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="alpha"],
                      cp=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"]),
                y=data_ts_Week[data_ts_Week$BASE==25, ]$Cases,
                I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="exp")
p25 <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==25])/sum(rec_25)*100, 2)
p25_before <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==25]
                        [1:round(summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"])])/
                      sum(rec_25[1:round(summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"])])*100, 2)
p25_after <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==25]
                        [(round(summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"])+1):length(rec_25)])/
                      sum(rec_25[(round(summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"])+1):length(rec_25)])*100, 2)
df1 <- data_ts_Week[data_ts_Week$BASE==25, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_25))
df2 <- data.frame(BASE=rep("25", length(rec_25)), Week=rep(seq(1:53), length(rec_25)/53), Cases=rec_25, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_25)),
                  Series=rep("Reconstructed", length(rec_25)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p1 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() + ylab("") +
  ylim(0, 30)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") + xlab("") + ggtitle("SUBAREA A")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

### BASE 26
GBV_mcmc_26 <- as.mcmc(fit_GBV_26_exp)
rec_26 <- x_rec(par=c(q0=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="q0"],
                      lambda=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="lambda"],
                      beta=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="beta"],
                      alpha=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="alpha"],
                      cp=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"]),
                y=data_ts_Week[data_ts_Week$BASE==26, ]$Cases,
                I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="exp")
p26 <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==26])/sum(rec_26)*100, 2)
p26_before <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==26]
                        [1:round(summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"])])/
                      sum(rec_26[1:round(summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"])])*100, 2)
p26_after <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==26]
                       [(round(summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"])+1):length(rec_26)])/
                     sum(rec_26[(round(summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"])+1):length(rec_26)])*100, 2)
df1 <- data_ts_Week[data_ts_Week$BASE==26, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_26))
df2 <- data.frame(BASE=rep("26", length(rec_26)), Week=rep(seq(1:53), length(rec_26)/53), Cases=rec_26, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_26)),
                  Series=rep("Reconstructed", length(rec_26)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p2 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() + ylab("") +
  ylim(0, 30)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") + xlab("") + ggtitle("SUBAREA B")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

### BASE 27
GBV_mcmc_27 <- as.mcmc(fit_GBV_27_exp)
rec_27 <- x_rec(par=c(q0=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="q0"],
                      lambda=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="lambda"],
                      beta=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="beta"],
                      alpha=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="alpha"],
                      cp=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"]),
                y=data_ts_Week[data_ts_Week$BASE==27, ]$Cases,
                I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="exp")
p27 <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==27])/sum(rec_27)*100, 2)
p27_before <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==27]
                        [1:round(summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"])])/
                      sum(rec_27[1:round(summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"])])*100, 2)
p27_after <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==27]
                       [(round(summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"])+1):length(rec_27)])/
                     sum(rec_27[(round(summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"])+1):length(rec_27)])*100, 2)
df1 <- data_ts_Week[data_ts_Week$BASE==27, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_27))
df2 <- data.frame(BASE=rep("27", length(rec_27)), Week=rep(seq(1:53), length(rec_27)/53), Cases=rec_27, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_27)),
                  Series=rep("Reconstructed", length(rec_27)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p3 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() + ylab("") +
  ylim(0, 30)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") + xlab("") + ggtitle("SUBAREA C")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

### BASE 31
GBV_mcmc_31 <- as.mcmc(fit_GBV_31_exp)
rec_31 <- x_rec(par=c(q0=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="q0"],
                      lambda=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="lambda"],
                      beta=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="beta"],
                      alpha=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="alpha"],
                      cp=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"]),
                y=data_ts_Week[data_ts_Week$BASE==31, ]$Cases,
                I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="exp")
p31 <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==31])/sum(rec_31)*100, 2)
p31_before <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==31]
                        [1:round(summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"])])/
                      sum(rec_31[1:round(summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"])])*100, 2)
p31_after <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==31]
                       [(round(summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"])+1):length(rec_31)])/
                     sum(rec_31[(round(summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"])+1):length(rec_31)])*100, 2)
df1 <- data_ts_Week[data_ts_Week$BASE==31, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_31))
df2 <- data.frame(BASE=rep("31", length(rec_31)), Week=rep(seq(1:53), length(rec_31)/53), Cases=rec_31, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_31)),
                  Series=rep("Reconstructed", length(rec_31)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p4 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() + ylab("") +
  ylim(0, 50)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") + xlab("") + ggtitle("SUBAREA D")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

### BASE 34
GBV_mcmc_34 <- as.mcmc(fit_GBV_34_exp)
rec_34 <- x_rec(par=c(q0=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="q0"],
                      lambda=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="lambda"],
                      beta=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="beta"],
                      alpha=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="alpha"],
                      cp=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"]),
                y=data_ts_Week[data_ts_Week$BASE==34, ]$Cases,
                I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="exp")
p34 <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==34])/sum(rec_34)*100, 2)
p34_before <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==34]
                        [1:round(summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"])])/
                      sum(rec_34[1:round(summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"])])*100, 2)
p34_after <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==34]
                       [(round(summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"])+1):length(rec_34)])/
                     sum(rec_34[(round(summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"])+1):length(rec_34)])*100, 2)
df1 <- data_ts_Week[data_ts_Week$BASE==34, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_34))
df2 <- data.frame(BASE=rep("34", length(rec_34)), Week=rep(seq(1:53), length(rec_34)/53), Cases=rec_34, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_34)),
                  Series=rep("Reconstructed", length(rec_34)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p5 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line()  + ylab("") +
  ylim(0, 15)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") + xlab("") + ggtitle("SUBAREA E")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

### BASE 35
GBV_mcmc_35 <- as.mcmc(fit_GBV_35_exp)
rec_35 <- x_rec(par=c(q0=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="q0"],
                      lambda=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="lambda"],
                      beta=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="beta"],
                      alpha=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="alpha"],
                      cp=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"]),
                y=data_ts_Week[data_ts_Week$BASE==35, ]$Cases,
                I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="exp")
p35 <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==35])/sum(rec_35)*100, 2)
p35_before <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==35]
                        [1:round(summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"])])/
                      sum(rec_35[1:round(summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"])])*100, 2)
p35_after <- round(sum(data_ts_Week$Cases[data_ts_Week$BASE==35]
                       [(round(summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"])+1):length(rec_35)])/
                     sum(rec_35[(round(summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"])+1):length(rec_35)])*100, 2)
df1 <- data_ts_Week[data_ts_Week$BASE==35, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_35))
df2 <- data.frame(BASE=rep("35", length(rec_35)), Week=rep(seq(1:53), length(rec_35)/53), Cases=rec_35, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(rec_35)),
                  Series=rep("Reconstructed", length(rec_35)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p6 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() + ylab("") + 
  ylim(0, 50)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1), legend.position = "none") + xlab("") + ggtitle("SUBAREA F")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

ggarrange(p1, p2, p3, p4, p5, p6,
          ncol = 3, nrow = 2)

### Percentage of registered cases by subarea
data_plot <- data.frame(Subarea=c("A", "B", "C", "D", "E", "F"), 
                        Value=c(p25, p26, p27, p31, p34, p35),
                        Before=c(p25_before, p26_before, p27_before,
                                 p31_before, p34_before, p35_before),
                        After=c(p25_after, p26_after, p27_after,
                                p31_after, p34_after, p35_after))

data_plot %>% ggplot(aes(x=Subarea, y=Value)) + ylim(c(0, 100)) +
  geom_point()+geom_ribbon(aes(1:nrow(data_plot),ymax=After,ymin=Before),alpha=0.1) +
  xlab("Subarea")+ylab("Percentage of reported cases") + scale_x_discrete(labels=data_plot$Subarea) +
  theme(axis.text.x = element_text(angle = 25))
