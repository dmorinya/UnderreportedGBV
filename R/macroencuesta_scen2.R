library(dplyr)

#### GENDER BASED VIOLENCE MACRO SURVEY (SPANISH GOVERNMENT)
dades <- read.table("Data/macroencuesta2019.csv", header=T, sep=",", dec=",")

### BARCELONA
dades <- dades[dades$PROV==8, c("VFISICA_P", "VSEXUAL_P", 
                                "VFISICA_EXP", "VSEXUAL_EXP",
                                "VFISICA", "V_SEXUAL", "M1P8",
                                "M2P8", "M3P2K")]

### PHYSICAL OR SEXUAL GBV
dades$VIOLENCIA <- ifelse(dades$VFISICA_P==0 & dades$VFISICA_EXP==0 &
                            dades$VFISICA==0 & dades$VSEXUAL_P==0 &
                            dades$VSEXUAL_EXP==0 & dades$V_SEXUAL==0, 0, 1)
dades$VIOLENCIA <- factor(dades$VIOLENCIA, levels=c(0, 1), labels=c("NO", "SI"))
p1 <- table(dades$VIOLENCIA)[2]/sum(!is.na(dades$VIOLENCIA)) ### PERCENTATGE DE DONES QUE HAN PATIT VIOLÈNCIA EN ALGUN MOMENT (fisica o sexual)

### USAGE OF PRIMARY CARE SYSTEM
p2 <- 0.5

### NUMBER OF WOMEN ASSIGNED TO EACH AREA
assignades25 <- readxl::read_xls("Data/ASSIG_VG.xls", sheet="25")
assignades26 <- readxl::read_xls("Data/ASSIG_VG.xls", sheet="26")
assignades27 <- readxl::read_xls("Data/ASSIG_VG.xls", sheet="27")
assignades31 <- readxl::read_xls("Data/ASSIG_VG.xls", sheet="31")
assignades34 <- readxl::read_xls("Data/ASSIG_VG.xls", sheet="34")
assignades35 <- readxl::read_xls("Data/ASSIG_VG.xls", sheet="35")

assignades25 <- assignades25[assignades25$EDAT>=16,]
assignades26 <- assignades26[assignades26$EDAT>=16,]
assignades27 <- assignades27[assignades27$EDAT>=16,]
assignades31 <- assignades31[assignades31$EDAT>=16,]
assignades34 <- assignades34[assignades34$EDAT>=16,]
assignades35 <- assignades35[assignades35$EDAT>=16,]

assignades25 <- assignades25 %>% group_by(BASE) %>% summarise(Assignades=sum(ASSIG))
assignades26 <- assignades26 %>% group_by(BASE) %>% summarise(Assignades=sum(ASSIG))
assignades27 <- assignades27 %>% group_by(BASE) %>% summarise(Assignades=sum(ASSIG))
assignades31 <- assignades31 %>% group_by(BASE) %>% summarise(Assignades=sum(ASSIG))
assignades34 <- assignades34 %>% group_by(BASE) %>% summarise(Assignades=sum(ASSIG))
assignades35 <- assignades35 %>% group_by(BASE) %>% summarise(Assignades=sum(ASSIG))

n_25 <- p1*p2*assignades25$Assignades/(53*12)
n_26 <- p1*p2*assignades26$Assignades/(53*12)
n_27 <- p1*p2*assignades27$Assignades/(53*12)
n_31 <- p1*p2*assignades31$Assignades/(53*12)
n_34 <- p1*p2*assignades34$Assignades/(53*12)
n_35 <- p1*p2*assignades35$Assignades/(53*12)

prop_survey <- c(n_25, n_26, n_27, n_31, n_34, n_35)

names(prop_survey) <- c("BASE 25", "BASE 26", "BASE 27", "BASE 31", "BASE 34", "BASE 35")
save(list=c("prop_survey"), file="R/prop_survey_scen2.RData")