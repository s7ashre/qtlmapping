rm(list = ls())
#setwd('D:/HHU/poya')
setwd('/home/shrestha/qggp/marvin/seedsizeQTL/h2estimate')
#setwd("/home/shrestha/mounts/projects/marvin/seedsizeQTL/h2estimate")
#setwd('~/qggp/marvin/TGW')
dir()
source("main_function_aem_h2.R")
library(Matrix)
library(lme4)
library(emmeans)

dat <- read.csv('phenoWidth_final.csv')
dat <- dat[,c(1,2,4:9,11,3,10,12:14)]
#dat <- dat[dat$Population == 'DRR_Parents',]
#load("alldata.check.RData")

traitnames <- colnames(dat)[10:14]
j <- 1

cat("trait:", traitnames[j], "\n")


formu.aem.origin <- as.formula(paste(traitnames[j], "~ (1|Location) + Genotype + (1|Genotype:Location)"))
formu.h2.origin <- as.formula(paste(traitnames[j], " ~ (1|Location) + (1|Genotype) + (1|Genotype:Location)"))

## fit model: G-fixed
fit.aem.origin <- lmer(formu.aem.origin, data = dat)

## fit model: G-random
fit.h2.origin <- lmer(formu.h2.origin, data = dat)
## my function
out.origin <- AEM.fun(fit = fit.aem.origin, select.name = "Genotype")
# obtain AEM 
AEM.my.origin <- out.origin$aem
colnames(AEM.my.origin)[1] <- "Genotype"
#print(AEM.my.origin)
write.csv(AEM.my.origin, file = paste('aem_', traitnames[j],'_poya.csv', sep = ""))
# obain error variance 
errorvar.my.origin <- out.origin$contrast.var.mean/2

## extract variance components
vcov.origin <- as.data.frame(VarCorr(fit.h2.origin))
vg.origin <- vcov.origin[vcov.origin$grp == "Genotype", "vcov"]

# heritability
h2 <- vg.origin/(vg.origin + errorvar.my.origin)
print(h2)
write.csv(h2, file = paste("h2_",traitnames[j],"_poya.csv", sep = ""))
save(fit.aem.origin, fit.h2.origin, out.origin, AEM.my.origin, h2, file = paste("fitmodel_output_",traitnames[j], ".RData", sep = ""))
