#### Evolution of q (impact of the training activity)
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

q_evolution <- function(q0, alpha, cp, mode)
{
  q <- vector()
  q[1:cp] <- q0
  if (mode == "linear")
  {
    q[(cp+1):alpha] <- q0+seq(1, round(alpha)-round(cp), 1)/((1/(1-q0))*(round(alpha)-round(cp)))
  }else{
    if (mode == "exp")
    {
      for (i in (cp+1):10000)
      {
        q[i] <- 1-(1-q0)*exp(-alpha*(i-cp))
      }
    }else{
      stop("Wrong mode.")
    }
  }
  return(q)
}

### LINEAR q_t GROWTH
### BASE 25
GBV_mcmc_25 <- as.mcmc(fit_GBV_25)
q_evolution_25 <- q_evolution(q0=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="q0"],
                              alpha=ceiling(summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="alpha"]),
                              cp=ceiling(summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"]),
                              mode="linear")
df <- data.frame(date=seq.Date(as.Date("2010-01-01"), by="week", length.out = length(q_evolution_25)), q=q_evolution_25)
p1 <- df %>% 
  ggplot(aes(x=date, y=q)) + geom_line() +
  ylim(0, 1)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA A")

### BASE 26
GBV_mcmc_26 <- as.mcmc(fit_GBV_26)
q_evolution_26 <- q_evolution(q0=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="q0"],
                              alpha=ceiling(summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="alpha"]),
                              cp=ceiling(summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"]),
                              mode="linear")
df <- data.frame(date=seq.Date(as.Date("2010-01-01"), by="week", length.out = length(q_evolution_26)), q=q_evolution_26)
p2 <- df %>% 
  ggplot(aes(x=date, y=q)) + geom_line() +
  ylim(0, 1)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA B")

### BASE 27
GBV_mcmc_27 <- as.mcmc(fit_GBV_27)
q_evolution_27 <- q_evolution(q0=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="q0"],
                              alpha=ceiling(summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="alpha"]),
                              cp=ceiling(summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"]),
                              mode="linear")
df <- data.frame(date=seq.Date(as.Date("2010-01-01"), by="week", length.out = length(q_evolution_27)), q=q_evolution_27)
p3 <- df %>% 
  ggplot(aes(x=date, y=q)) + geom_line() +
  ylim(0, 1)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA C")

### BASE 31
GBV_mcmc_31 <- as.mcmc(fit_GBV_31)
q_evolution_31 <- q_evolution(q0=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="q0"],
                              alpha=ceiling(summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="alpha"]),
                              cp=ceiling(summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"]),
                              mode="linear")
df <- data.frame(date=seq.Date(as.Date("2010-01-01"), by="week", length.out = length(q_evolution_31)), q=q_evolution_31)
p4 <- df %>% 
  ggplot(aes(x=date, y=q)) + geom_line() +
  ylim(0, 1)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA D")

### BASE 34
GBV_mcmc_34 <- as.mcmc(fit_GBV_34)
q_evolution_34 <- q_evolution(q0=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="q0"],
                              alpha=ceiling(summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="alpha"]),
                              cp=ceiling(summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"]),
                              mode="linear")
df <- data.frame(date=seq.Date(as.Date("2010-01-01"), by="week", length.out = length(q_evolution_34)), q=q_evolution_34)
p5 <- df %>% 
  ggplot(aes(x=date, y=q)) + geom_line() +
  ylim(0, 1)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA E")

### BASE 35
GBV_mcmc_35 <- as.mcmc(fit_GBV_35)
q_evolution_35 <- q_evolution(q0=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="q0"],
                              alpha=ceiling(summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="alpha"]),
                              cp=ceiling(summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"]),
                              mode="linear")
df <- data.frame(date=seq.Date(as.Date("2010-01-01"), by="week", length.out = length(q_evolution_35)), q=q_evolution_35)
p6 <- df %>% 
  ggplot(aes(x=date, y=q)) + geom_line() +
  ylim(0, 1)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA F")

ggarrange(p1, p2, p3, p4, p5, p6,
          ncol = 3, nrow = 2)

### EXPONENTIAL q_t GROWTH
### BASE 25
GBV_mcmc_25 <- as.mcmc(fit_GBV_25_exp)
q_evolution_25 <- q_evolution(q0=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="q0"],
                              alpha=summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="alpha"],
                              cp=ceiling(summary(GBV_mcmc_25)$quantiles[,3][rownames(summary(GBV_mcmc_25)$quantiles)=="cp"]),
                              mode="exp")
q_evolution_25 <- q_evolution_25[1:which(q_evolution_25>0.99)[1]]
df <- data.frame(date=seq.Date(as.Date("2010-01-01"), by="week", length.out = length(q_evolution_25)), q=q_evolution_25)
p1 <- df %>% 
  ggplot(aes(x=date, y=q)) + geom_line() +
  ylim(0, 1)+scale_x_date(date_labels = "%m-%Y", date_breaks = "3 year")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA A")

### BASE 26
GBV_mcmc_26 <- as.mcmc(fit_GBV_26_exp)
q_evolution_26 <- q_evolution(q0=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="q0"],
                              alpha=summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="alpha"],
                              cp=ceiling(summary(GBV_mcmc_26)$quantiles[,3][rownames(summary(GBV_mcmc_26)$quantiles)=="cp"]),
                              mode="exp")
q_evolution_26 <- q_evolution_26[1:which(q_evolution_26>0.99)[1]]
df <- data.frame(date=seq.Date(as.Date("2010-01-01"), by="week", length.out = length(q_evolution_26)), q=q_evolution_26)
p2 <- df %>% 
  ggplot(aes(x=date, y=q)) + geom_line() +
  ylim(0, 1)+scale_x_date(date_labels = "%m-%Y", date_breaks = "3 year")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA B")

### BASE 27
GBV_mcmc_27 <- as.mcmc(fit_GBV_27_exp)
q_evolution_27 <- q_evolution(q0=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="q0"],
                              alpha=summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="alpha"],
                              cp=ceiling(summary(GBV_mcmc_27)$quantiles[,3][rownames(summary(GBV_mcmc_27)$quantiles)=="cp"]),
                              mode="exp")
q_evolution_27 <- q_evolution_27[1:which(q_evolution_27>0.99)[1]]
df <- data.frame(date=seq.Date(as.Date("2010-01-01"), by="week", length.out = length(q_evolution_27)), q=q_evolution_27)
p3 <- df %>% 
  ggplot(aes(x=date, y=q)) + geom_line() +
  ylim(0, 1)+scale_x_date(date_labels = "%m-%Y", date_breaks = "5 year")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA C")

### BASE 31
GBV_mcmc_31 <- as.mcmc(fit_GBV_31_exp)
q_evolution_31 <- q_evolution(q0=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="q0"],
                              alpha=summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="alpha"],
                              cp=ceiling(summary(GBV_mcmc_31)$quantiles[,3][rownames(summary(GBV_mcmc_31)$quantiles)=="cp"]),
                              mode="exp")
q_evolution_31 <- q_evolution_31[1:which(q_evolution_31>0.99)[1]]
df <- data.frame(date=seq.Date(as.Date("2010-01-01"), by="week", length.out = length(q_evolution_31)), q=q_evolution_31)
p4 <- df %>% 
  ggplot(aes(x=date, y=q)) + geom_line() +
  ylim(0, 1)+scale_x_date(date_labels = "%m-%Y", date_breaks = "3 year")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA D")

### BASE 34
GBV_mcmc_34 <- as.mcmc(fit_GBV_34_exp)
q_evolution_34 <- q_evolution(q0=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="q0"],
                              alpha=summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="alpha"],
                              cp=ceiling(summary(GBV_mcmc_34)$quantiles[,3][rownames(summary(GBV_mcmc_34)$quantiles)=="cp"]),
                              mode="exp")
q_evolution_34 <- q_evolution_34[1:which(q_evolution_34>0.9999)[1]]
df <- data.frame(date=seq.Date(as.Date("2010-01-01"), by="week", length.out = length(q_evolution_34)), q=q_evolution_34)
p5 <- df %>% 
  ggplot(aes(x=date, y=q)) + geom_line() +
  ylim(0, 1)+scale_x_date(date_labels = "%m-%Y", date_breaks = "6 month")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA E")

### BASE 35
GBV_mcmc_35 <- as.mcmc(fit_GBV_35_exp)
q_evolution_35 <- q_evolution(q0=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="q0"],
                              alpha=summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="alpha"],
                              cp=ceiling(summary(GBV_mcmc_35)$quantiles[,3][rownames(summary(GBV_mcmc_35)$quantiles)=="cp"]),
                              mode="exp")
q_evolution_35 <- q_evolution_35[1:which(q_evolution_35>0.99)[1]]
df <- data.frame(date=seq.Date(as.Date("2010-01-01"), by="week", length.out = length(q_evolution_35)), q=q_evolution_35)
p6 <- df %>% 
  ggplot(aes(x=date, y=q)) + geom_line() +
  ylim(0, 1)+scale_x_date(date_labels = "%m-%Y", date_breaks = "3 year")+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=60, hjust=1)) + xlab("") + ggtitle("SUBAREA F")

ggarrange(p1, p2, p3, p4, p5, p6,
          ncol = 3, nrow = 2)