### Simulate the process for each subarea 
library(readxl)
library(dplyr)
library(MMWRweek)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(lubridate)
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

simulate_data <- function(q0, lambda, beta, cp, alpha, I, mode) 
{
  x <- rpois(636, lambda+I*beta)
  y <- vector()
  for (i in 1:cp) 
  {
    y[i] <- rbinom(1, x[i], q0)
  }  
  if (mode=="exp")
  {
    q <- 1-(1-q0)*exp(-alpha*seq(1, (636-cp), 1)) #exponential
  }else{
    if (mode=="linear")
    {
      q <- q0+seq(1, round(alpha)-round(cp), 1)/((1/(1-q0))*(round(alpha)-round(cp))) #linear
    }else{
      stop("Wrong mode.")
    }
  }
  for (i in (cp+1):636)
  {
    y[i] <- rbinom(1, x[i], q[i-cp])
  }
  return(y)
}

### LINEAR q_t GROWTH
### BASE 25
set.seed(232023)
GBV_mcmc_25 <- as.mcmc(fit_GBV_25)
y_sim_25 <- simulate_data(q0=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="q0"],
                          lambda=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="lambda"],
                          beta=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="beta"],
                          cp=round(summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"]),
                          alpha=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="alpha"],
                          I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="linear")
df1 <- data_ts_Week[data_ts_Week$BASE==25, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_25))
df2 <- data.frame(BASE=rep("25", length(y_sim_25)), Week=rep(seq(1:53), length(y_sim_25)/53), Cases=y_sim_25, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_25)),
                  Series=rep("Simulated", length(y_sim_25)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p1 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() +
  ylim(0, 20)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA A")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

### BASE 26
set.seed(232023)
GBV_mcmc_26 <- as.mcmc(fit_GBV_26)
y_sim_26 <- simulate_data(q0=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="q0"],
                          lambda=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="lambda"],
                          beta=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="beta"],
                          cp=round(summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"]),
                          alpha=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="alpha"],
                          I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="linear")
df1 <- data_ts_Week[data_ts_Week$BASE==26, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_26))
df2 <- data.frame(BASE=rep("26", length(y_sim_26)), Week=rep(seq(1:53), length(y_sim_26)/53), Cases=y_sim_26, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_26)),
                  Series=rep("Simulated", length(y_sim_26)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p2 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() +
  ylim(0, 20)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA B")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

### BASE 27
set.seed(232023)
GBV_mcmc_27 <- as.mcmc(fit_GBV_27)
y_sim_27 <- simulate_data(q0=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="q0"],
                          lambda=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="lambda"],
                          beta=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="beta"],
                          cp=round(summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"]),
                          alpha=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="alpha"],
                          I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="linear")
df1 <- data_ts_Week[data_ts_Week$BASE==27, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_27))
df2 <- data.frame(BASE=rep("27", length(y_sim_27)), Week=rep(seq(1:53), length(y_sim_27)/53), Cases=y_sim_27, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_27)),
                  Series=rep("Simulated", length(y_sim_27)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p3 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() +
  ylim(0, 20)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA C")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

### BASE 31
set.seed(232023)
GBV_mcmc_31 <- as.mcmc(fit_GBV_31)
y_sim_31 <- simulate_data(q0=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="q0"],
                          lambda=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="lambda"],
                          beta=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="beta"],
                          cp=round(summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"]),
                          alpha=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="alpha"],
                          I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="linear")
df1 <- data_ts_Week[data_ts_Week$BASE==31, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_31))
df2 <- data.frame(BASE=rep("31", length(y_sim_31)), Week=rep(seq(1:53), length(y_sim_31)/53), Cases=y_sim_31, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_31)),
                  Series=rep("Simulated", length(y_sim_31)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p4 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() +
  ylim(0, 20)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA D")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

### BASE 34
set.seed(232023)
GBV_mcmc_34 <- as.mcmc(fit_GBV_34)
y_sim_34 <- simulate_data(q0=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="q0"],
                          lambda=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="lambda"],
                          beta=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="beta"],
                          cp=round(summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"]),
                          alpha=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="alpha"],
                          I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="linear")
df1 <- data_ts_Week[data_ts_Week$BASE==34, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_34))
df2 <- data.frame(BASE=rep("34", length(y_sim_34)), Week=rep(seq(1:53), length(y_sim_34)/53), Cases=y_sim_34, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_34)),
                  Series=rep("Simulated", length(y_sim_34)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p5 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() +
  ylim(0, 20)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA E")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

### BASE 35
set.seed(232023)
GBV_mcmc_35 <- as.mcmc(fit_GBV_35)
y_sim_35 <- simulate_data(q0=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="q0"],
                          lambda=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="lambda"],
                          beta=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="beta"],
                          cp=round(summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"]),
                          alpha=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="alpha"],
                          I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="linear")
df1 <- data_ts_Week[data_ts_Week$BASE==35, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_35))
df2 <- data.frame(BASE=rep("35", length(y_sim_35)), Week=rep(seq(1:53), length(y_sim_35)/53), Cases=y_sim_35, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_35)),
                  Series=rep("Simulated", length(y_sim_35)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p6 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() +
  ylim(0, 20)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA F")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

ggarrange(p1, p2, p3, p4, p5, p6,
          ncol = 3, nrow = 2)

### EXPONENTIAL q_t GROWTH
### BASE 25
set.seed(232023)
GBV_mcmc_25 <- as.mcmc(fit_GBV_25_exp)
y_sim_25 <- simulate_data(q0=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="q0"],
                          lambda=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="lambda"],
                          beta=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="beta"],
                          cp=round(summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"]),
                          alpha=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="alpha"],
                          I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="exp")
df1 <- data_ts_Week[data_ts_Week$BASE==25, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_25))
df2 <- data.frame(BASE=rep("25", length(y_sim_25)), Week=rep(seq(1:53), length(y_sim_25)/53), Cases=y_sim_25, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_25)),
                  Series=rep("Simulated", length(y_sim_25)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p1 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() +
  ylim(0, 20)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA A")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

### BASE 26
set.seed(232023)
GBV_mcmc_26 <- as.mcmc(fit_GBV_26_exp)
y_sim_26 <- simulate_data(q0=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="q0"],
                          lambda=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="lambda"],
                          beta=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="beta"],
                          cp=round(summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"]),
                          alpha=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="alpha"],
                          I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="exp")
df1 <- data_ts_Week[data_ts_Week$BASE==26, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_26))
df2 <- data.frame(BASE=rep("26", length(y_sim_26)), Week=rep(seq(1:53), length(y_sim_26)/53), Cases=y_sim_26, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_26)),
                  Series=rep("Simulated", length(y_sim_26)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p2 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() +
  ylim(0, 20)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA B")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

### BASE 27
set.seed(232023)
GBV_mcmc_27 <- as.mcmc(fit_GBV_27_exp)
y_sim_27 <- simulate_data(q0=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="q0"],
                          lambda=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="lambda"],
                          beta=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="beta"],
                          cp=round(summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"]),
                          alpha=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="alpha"],
                          I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="exp")
df1 <- data_ts_Week[data_ts_Week$BASE==27, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_27))
df2 <- data.frame(BASE=rep("27", length(y_sim_27)), Week=rep(seq(1:53), length(y_sim_27)/53), Cases=y_sim_27, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_27)),
                  Series=rep("Simulated", length(y_sim_27)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p3 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() +
  ylim(0, 20)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA C")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

### BASE 31
set.seed(232023)
GBV_mcmc_31 <- as.mcmc(fit_GBV_31_exp)
y_sim_31 <- simulate_data(q0=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="q0"],
                          lambda=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="lambda"],
                          beta=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="beta"],
                          cp=round(summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"]),
                          alpha=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="alpha"],
                          I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="exp")
df1 <- data_ts_Week[data_ts_Week$BASE==31, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_31))
df2 <- data.frame(BASE=rep("31", length(y_sim_31)), Week=rep(seq(1:53), length(y_sim_31)/53), Cases=y_sim_31, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_31)),
                  Series=rep("Simulated", length(y_sim_31)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p4 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() +
  ylim(0, 20)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA D")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

### BASE 34
set.seed(232023)
GBV_mcmc_34 <- as.mcmc(fit_GBV_34_exp)
y_sim_34 <- simulate_data(q0=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="q0"],
                          lambda=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="lambda"],
                          beta=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="beta"],
                          cp=round(summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"]),
                          alpha=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="alpha"],
                          I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="exp")
df1 <- data_ts_Week[data_ts_Week$BASE==34, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_34))
df2 <- data.frame(BASE=rep("34", length(y_sim_34)), Week=rep(seq(1:53), length(y_sim_34)/53), Cases=y_sim_34, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_34)),
                  Series=rep("Simulated", length(y_sim_34)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p5 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() +
  ylim(0, 20)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA E")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)

### BASE 35
set.seed(232023)
GBV_mcmc_35 <- as.mcmc(fit_GBV_35_exp)
y_sim_35 <- simulate_data(q0=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="q0"],
                          lambda=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="lambda"],
                          beta=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="beta"],
                          cp=round(summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"]),
                          alpha=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="alpha"],
                          I=c(rep(0, 541), rep(1, 14), rep(0, 81)), mode="exp")
df1 <- data_ts_Week[data_ts_Week$BASE==35, ]
df1$date <- seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_35))
df2 <- data.frame(BASE=rep("35", length(y_sim_35)), Week=rep(seq(1:53), length(y_sim_35)/53), Cases=y_sim_35, 
                  date=seq.Date(from=as.Date("2010-01-01", format="%Y-%m-%d"), by="week", length.out=length(y_sim_35)),
                  Series=rep("Simulated", length(y_sim_35)))
df2$Year <- year(df2$date)
df2 <- df2[, c(1, 6, 2, 3, 4, 5)]
df1$Series <- "Registered"
df <- rbind(df1, df2)
p6 <- df %>% 
  ggplot(aes(x=date, y=Cases, col=Series)) + geom_line() +
  ylim(0, 20)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA F")+
  geom_vline(xintercept=df$date[round(summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"])], linetype=1, color="black", linewidth=1) +
  geom_vline(xintercept=df$date[df$date=="2020-03-13"], linetype=1, color="red", linewidth=1)+
  geom_vline(xintercept=df$date[df$date=="2020-06-26"], linetype=1, color="red", linewidth=1)+
  annotate("rect", xmin = df$date[df$date=="2020-03-13"], xmax = df$date[df$date=="2020-06-26"], ymin = -Inf, ymax = Inf,
           alpha = .2)
ggarrange(p1, p2, p3, p4, p5, p6,
          ncol = 3, nrow = 2)