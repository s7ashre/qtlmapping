setwd('/home/shrestha/mounts/project/shrestha/marvin/seedsizeQTL/experiments/qtlmapping/Width')
library(qtl)
datposition <- read.csv('/home/shrestha/mounts/hpc_home/experiments/qtlmapping/Barley_MarkerPositions3Versions.csv')
datWidth <- readRDS('Widthcross_new.RDS')
#######HvDRR02#######
hvdrr <- datWidth$HvDRR02 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) # 134inds and 93.9% genotyped
png('summary_Width02.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width02.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-107818   2 167.3 14.64 0.0000
#JHI-Hv50k-2016-335087   5 199.9  3.63 0.0170
#JHI-Hv50k-2016-402874   6  79.8  4.17 0.0065
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 5, 6), c(167.3, 199.9, 79.8), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-231131   4  45 3.43 0.0255
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 4, 5, 6), c(167.3, 45, 199.9, 79.8), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 2@167.3   2 167.301     2
#Q2  4@45.0   4  45.001     2
#Q3 5@203.2   5 203.172     2
#Q4  6@79.0   6  79.000     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1 + Q2 + Q3 + Q4)) #adds if there are QTL effects
pt(abs(3.937), 125, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width02.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr2
find.marker(hvdrrprob, 2, 167.301) #JHI-Hv50k-2016-107818
png('effetplot_ch2_Width02.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-107818")
dev.off()
png('effetplotsca_ch2_Width02.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-107818") #Scatter plot
dev.off()

#effect plot marker chr4
find.marker(hvdrrprob, 4, 45.001) #JHI-Hv50k-2016-231133
png('effetplot_ch4_Width02.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-231133')
dev.off()
png('effetplotsca_ch4_Width20.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-231133') #Scatter plot
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 203.172) #BOPA2_12_30930
png('effetplot_ch5_Width02.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'BOPA2_12_30930')
dev.off()
png('effetplotsca_ch5_Width02.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'BOPA2_12_30930') #Scatter plot
dev.off()

#effect plot marker chr6
find.marker(hvdrrprob, 6, 79) #JHI-Hv50k-2016-402874
png('effetplot_ch6_Width02.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-402874')
dev.off()
png('effetplotsca_ch6_Width02.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-402874') #Scatter plot
dev.off()
#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-107351   2 166.3059 16.02840
#JHI-Hv50k-2016-107818   2 167.3005 17.56922
#JHI-Hv50k-2016-109301   2 170.7037 14.71128
datposition[(which(datposition == 'JHI-Hv50k-2016-107351')),9] #580550249
datposition[(which(datposition == 'JHI-Hv50k-2016-107818')),9] #582526776
datposition[(which(datposition == 'JHI-Hv50k-2016-109301')),9] #591197589

hvdrrCI <- lodint(hvdrr_rqtl, chr = 4, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#SCRI_RS_98443           4 38.04715 1.114544
#JHI-Hv50k-2016-231131   4 45.00114 3.300993
#JHI-Hv50k-2016-231677   4 49.26122 1.676463
datposition[(which(datposition == 'SCRI_RS_98443')),9] #14583570
datposition[(which(datposition == 'JHI-Hv50k-2016-231133')),9] #17408045
datposition[(which(datposition == 'JHI-Hv50k-2016-231677')),9] #19780549

hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 3, expandtomarkers = TRUE)
hvdrrCI
#SCRI_RS_130992          5 198.1170 4.891383
#BOPA2_12_30930          5 203.1718 7.224296
#JHI-Hv50k-2016-338178   5 215.8993 4.689457
datposition[(which(datposition == 'SCRI_RS_130992')),9] #534439825
datposition[(which(datposition == 'BOPA2_12_30930')),9] #537284215
datposition[(which(datposition == 'JHI-Hv50k-2016-338178')),9] #542769265

hvdrrCI <- lodint(hvdrr_rqtl, chr = 6, qtl.index = 4, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-382621   6 56.66057 2.083564
#c6.loc79                6 79.00000 3.595754
#JHI-Hv50k-2016-407207   6 84.78945 2.069384
datposition[(which(datposition == 'JHI-Hv50k-2016-382621')),9] #33392186
datposition[(which(datposition == 'JHI-Hv50k-2016-402874')),9] #399721803
datposition[(which(datposition == 'JHI-Hv50k-2016-407207')),9] #456127431

#######HvDRR03#######
hvdrr <- datWidth$HvDRR03 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#hvdrrpheno <- hvdrr$pheno
#hvdrr$pheno[83, 1] <- NA
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #65 inds and 94.9% genotyped
png('summary_Width03.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width03.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#c2.loc149   2 149 8.01    0
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2), c(149), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 2@149.0   2 149     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1)) #adds if there are QTL effects
pt(abs(6.936), 62, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width03.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr2
find.marker(hvdrrprob, 2, 149) #JHI-Hv50k-2016-107384
png('effetplot_ch2_Width03.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-107384")
dev.off()
png('effetplotsca_ch2_Width03.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-107384") #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-106749   2 143.1319 6.435250
#c2.loc149               2 149.0000 8.007573
#JHI-Hv50k-2016-108011   2 151.1707 5.934100
datposition[(which(datposition == 'JHI-Hv50k-2016-106749')),9] #578660476
datposition[(which(datposition == 'JHI-Hv50k-2016-107384')),9] #580717220
datposition[(which(datposition == 'JHI-Hv50k-2016-108011')),9] #583148593


#######HvDRR04#######
hvdrr <- datWidth$HvDRR04 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) # 136 inds and 93.3% genotyped
png('summary_Width04.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width04.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr7
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#c7.loc192   7 192 14.5    0
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(7), c(192), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 7@192.0   7 192     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1)) #adds if there are QTL effects
pt(abs(-9.205), 133, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width04.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr7
find.marker(hvdrrprob, 7, 192) #JHI-Hv50k-2016-491158
png('effetplot_ch7_Width04.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-491158")
dev.off()
png('effetplotsca_ch7_Width04.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-491158") #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 7, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-491158   7 188.7255 12.89136 #SCRI_RS_124478 187.29cM
#c7.loc192               7 192.0000 14.47179
#BOPA2_12_30026          7 195.8720 12.42120
datposition[(which(datposition == 'SCRI_RS_124478')),9] #518220704
datposition[(which(datposition == 'JHI-Hv50k-2016-491158')),9] #520712431
datposition[(which(datposition == 'BOPA2_12_30026')),9] #542585738
hvdrr$geno$`7`

#######HvDRR07#######
hvdrr <- datWidth$HvDRR07 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #98 inds and 93.8% genotyped
png('summary_Width07.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width07.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr3, 5 and 7
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-205634   3 158 3.57 0.01700
#JHI-Hv50k-2016-285916   5  66 3.59 0.01625
#c7.loc142               7 142 4.12 0.00575
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(3, 5, 7), c(158, 66, 142), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 3@163.0   3 163.000     2
#Q2  5@63.2   5  63.215     2
#Q3 7@141.0   7 141.000     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2+Q3)) #adds if there are QTL effects
pt(abs(4.468), 91, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width07.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr3
find.marker(hvdrrprob, 3, 163) #JHI-Hv50k-2016-205634
png('effetplot_ch3_Width.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-205634")
dev.off()
png('effetplotsca_ch3_Width.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-205634") #Scatter plot
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 63.215) #JHI-Hv50k-2016-285820
png('effetplot_ch5_Width.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-285820')
dev.off()
png('effetplotsca_ch5_Width.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-285820') #Scatter plot
dev.off()

#effect plot marker chr7
find.marker(hvdrrprob, 7, 141) #JHI-Hv50k-2016-491383
png('effetplot_ch7_Width.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-491383')
dev.off()
png('effetplotsca_ch7_Width.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-491383') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 3, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#SCRI_RS_150944   3 152.5719 1.643225
#c3.loc163        3 163.0000 4.098412
#SCRI_RS_169325   3 187.1746 1.405280
datposition[(which(datposition == 'SCRI_RS_150944')),9] #568277774
datposition[(which(datposition == 'JHI-Hv50k-2016-205634')),9] #573139245
datposition[(which(datposition == 'SCRI_RS_169325')),9] #596671509

hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-283752   5 56.26116 3.476344
#JHI-Hv50k-2016-285820   5 63.21480 5.090461
#JHI-Hv50k-2016-307382   5 90.85661 3.028652
datposition[(which(datposition == 'JHI-Hv50k-2016-283752')),9] #16189986
datposition[(which(datposition == 'JHI-Hv50k-2016-285820')),9] #22394787
datposition[(which(datposition == 'JHI-Hv50k-2016-307382')),9] #435714162

hvdrrCI <- lodint(hvdrr_rqtl, chr = 7, qtl.index = 3, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-486580   7 135.0823 3.664461
#c7.loc141               7 141.0000 5.469315
#SCRI_RS_200020          7 145.6885 2.513518
datposition[(which(datposition == 'JHI-Hv50k-2016-486580')),9] #431927444
datposition[(which(datposition == 'JHI-Hv50k-2016-491383')),9] #522867563
datposition[(which(datposition == 'SCRI_RS_200020')),9] #554516185

#######HvDRR08#######
hvdrr <- datWidth$HvDRR08 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #95 inds and 92.7% genotyped
png('summary_Width08.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width08.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr3 and 5
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#c3.loc263        3 263 4.31 0.00475
#BOPA2_12_30062   5 335 3.40 0.02750
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(3, 5), c(263, 335), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 3@267.0   3 267.00     2
#Q2 5@343.1   5 343.07     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2)) #adds if there are QTL effects
pt(abs(-4.143), 90, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width08.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr3
find.marker(hvdrrprob, 3, 267) #JHI-Hv50k-2016-211367
png('effetplot_ch3_Width08.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-211367")
dev.off()
png('effetplotsca_ch3_Width08.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-211367") #Scatter plot
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 343) #JHI-Hv50k-2016-346486
png('effetplot_ch5_Width08.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-346486')
dev.off()
png('effetplotsca_ch5_Width08.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-346486') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 3, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#SCRI_RS_165334          3 224.5493 2.599322
#c3.loc267               3 267.0000 4.116693
#JHI-Hv50k-2016-217798   3 318.1903 2.554132
datposition[(which(datposition == 'SCRI_RS_165334')),9] #572326171
datposition[(which(datposition == 'JHI-Hv50k-2016-211367')),9] #593614464
datposition[(which(datposition == 'JHI-Hv50k-2016-217798')),9] #606802787

hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#SCRI_RS_166296          5 157.3690 1.718131
#JHI-Hv50k-2016-346486   5 343.0668 3.528937
#JHI-Hv50k-2016-349937   5 354.2024 1.639561
datposition[(which(datposition == 'SCRI_RS_166296')),9] #464855148
datposition[(which(datposition == 'JHI-Hv50k-2016-346486')),9] #559900040
datposition[(which(datposition == 'JHI-Hv50k-2016-349937')),9] #569856578


#######HvDRR09#######
hvdrr <- datWidth$HvDRR09 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #100 inds and 32.2% genotyped
png('summary_Width.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width09.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #no significynt peaks
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)

#######HvDRR10#######
hvdrr <- datWidth$HvDRR10 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #89 inds and 93.9% genotyped
png('summary_Width10.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width10.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2, 3 and 5
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#c2.loc191               2 191 3.26 0.0345
#JHI-Hv50k-2016-205741   3 125 5.39 0.0000
#c5.loc131               5 131 3.16 0.0423
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 3, 5), c(191, 125, 131), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 2@189.8   2 189.838     2
#Q2 3@124.8   3 124.757     2
#Q3  5@91.4   5  91.414     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2+Q3)) #adds if there are QTL effects
pt(abs(-4.546), 82, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr2
find.marker(hvdrrprob, 2, 189.838) #
png('effetplot_ch2_Width10.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-130413")
dev.off()
png('effetplotsca_ch2_Width10.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-130413") #Scatter plot
dev.off()

#effect plot marker chr3
find.marker(hvdrrprob, 3, 124.757) #JHI-Hv50k-2016-205741
png('effetplot_ch3_Width10.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-205741')
dev.off()
png('effetplotsca_ch3_Width10.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-205741') #Scatter plot
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 91.414) #JHI-Hv50k-2016-309463
png('effetplot_ch5_Width10.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-309463')
dev.off()
png('effetplotsca_ch5_Width10.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-309463') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#BOPA1_ABC03181-1-1-164   2  90.62116 4.039316
#JHI-Hv50k-2016-130260    2 189.83845 5.634983
#JHI-Hv50k-2016-132236    2 197.89142 3.283853
datposition[(which(datposition == 'BOPA1_ABC03181-1-1-164')),9] #510838917
datposition[(which(datposition == 'JHI-Hv50k-2016-130413')),9] #644544848
datposition[(which(datposition == 'JHI-Hv50k-2016-132236')),9] #648297526

hvdrrCI <- lodint(hvdrr_rqtl, chr = 3, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-205209   3 121.4699 6.331474
#JHI-Hv50k-2016-205741   3 124.7567 8.372046
#JHI-Hv50k-2016-206235   3 128.1700 6.538151
datposition[(which(datposition == 'JHI-Hv50k-2016-205209')),9] #571533368
datposition[(which(datposition == 'JHI-Hv50k-2016-205741')),9] #573504167
datposition[(which(datposition == 'JHI-Hv50k-2016-206235')),9] #576274814

hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 3, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-297745   5  69.35613 2.547964
#JHI-Hv50k-2016-309463   5  91.41398 4.205889
#SCRI_RS_45011           5 132.03557 2.640240
datposition[(which(datposition == 'JHI-Hv50k-2016-297745')),9] #284026804
datposition[(which(datposition == 'JHI-Hv50k-2016-309463')),9] #453224475
datposition[(which(datposition == 'SCRI_RS_45011')),9] #508830904

#######HvDRR11#######
hvdrr <- datWidth$HvDRR11 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #97 inds and 93,9% genotyped
png('summary_Width11.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width11.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#c2.loc69   2  69 3.64 0.011
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2), c(69), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-282637   5 53.1 3.33 0.024
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 5), c(69, 53.1), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 2@107.6   2 107.564     2
#Q2  5@53.1   5  53.088     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2)) #adds if there are QTL effects
pt(abs(-4.235), 92, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width11.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr2
find.marker(hvdrrprob, 2, 107.56) #JHI-Hv50k-2016-111319
png('effetplot_ch2_Width11.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-111319")
dev.off()
png('effetplotsca_ch2_Width11.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-111319") #Scatter plot
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 53.09) #JHI-Hv50k-2016-282637
png('effetplot_ch5_Width11.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-282637')
dev.off()
png('effetplotsca_ch5_Width11.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-282637') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-80736    2  52.45301 2.489854
#JHI-Hv50k-2016-111319   2 107.56360 5.072400
#JHI-Hv50k-2016-118919   2 132.51647 3.006793
datposition[(which(datposition == 'JHI-Hv50k-2016-80736')),9] #56169213
datposition[(which(datposition == 'JHI-Hv50k-2016-111319')),9] #604464234
datposition[(which(datposition == 'JHI-Hv50k-2016-118919')),9] #624036873

hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-281033   5 36.58245 2.102492
#JHI-Hv50k-2016-282637   5 53.08788 3.677794
#JHI-Hv50k-2016-284304   5 59.99111 2.140782
datposition[(which(datposition == 'JHI-Hv50k-2016-281033')),9] #8549699
datposition[(which(datposition == 'JHI-Hv50k-2016-282637')),9] #11280864
datposition[(which(datposition == 'JHI-Hv50k-2016-284304')),9] #17583787

#######HvDRR12#######
hvdrr <- datWidth$HvDRR12 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #66 inds and 93,8% genotyped
png('summary_Width12.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width12.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2, 5
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#SCRI_RS_153811          2 203.5 3.15 0.0328
#JHI-Hv50k-2016-297876   5  60.3 2.94 0.0485
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 5), c(203, 60.3), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 2@203.5   2 203.489     2
#Q2  5@60.3   5  60.276     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2)) #adds if there are QTL effects
pt(abs(3.839), 61, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width12.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr2
find.marker(hvdrrprob, 2, 203.489) #SCRI_RS_153811
png('effetplot_ch2_Width12.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "SCRI_RS_153811")
dev.off()
png('effetplotsca_ch2_Width12.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "SCRI_RS_153811") #Scatter plot
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 60.276) #JHI-Hv50k-2016-297876
png('effetplot_ch5_Width12.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-297876')
dev.off()
png('effetplotsca_ch5_Width12.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-297876') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-74500    2  52.02371 0.9896354
#SCRI_RS_153811          2 203.48883 3.2174956
#JHI-Hv50k-2016-126155   2 210.81597 1.9399680
datposition[(which(datposition == 'JHI-Hv50k-2016-74500')),9] #26383663
datposition[(which(datposition == 'SCRI_RS_153811')),9] #633281761
datposition[(which(datposition == 'JHI-Hv50k-2016-126155')),9] #636135071

hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-290020   5 53.89448 1.042803
#JHI-Hv50k-2016-297876   5 60.27645 3.013279
#JHI-Hv50k-2016-312375   5 83.20775 1.411044
datposition[(which(datposition == 'JHI-Hv50k-2016-290020')),9] #41926504
datposition[(which(datposition == 'JHI-Hv50k-2016-297876')),9] #402745588
datposition[(which(datposition == 'JHI-Hv50k-2016-312375')),9] #480377334


#######HvDRR13#######
hvdrr <- datWidth$HvDRR13 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #68 inds and 94.3% genotyped
png('summary_Width13.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width13.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr3 and 5
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#c3.loc56                3  56 3.35 0.028
#JHI-Hv50k-2016-330574   5 200 3.74 0.012
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(3, 5), c(56, 200), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1  3@54.8   3  54.833     2
#Q2 5@201.0   5 201.000     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2)) #adds if there are QTL effects
pt(abs(-4.037), 63, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width13.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr3
find.marker(hvdrrprob, 3, 54.833) #JHI-Hv50k-2016-161342
png('effetplot_ch3_Width13.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-161342")
dev.off()
png('effetplotsca_ch3_Width13.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-161342") #Scatter plot
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 201) #JHI-Hv50k-2016-330574
png('effetplot_ch5_Width13.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-330574')
dev.off()
png('effetplotsca_ch5_Width13.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-330574') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 3, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-158683   3 40.29438 1.353433
#JHI-Hv50k-2016-161342   3 54.83304 3.303428
#JHI-Hv50k-2016-182700   3 72.01723 1.627586
datposition[(which(datposition == 'JHI-Hv50k-2016-158683')),9] #17984177
datposition[(which(datposition == 'JHI-Hv50k-2016-161342')),9] #26096119
datposition[(which(datposition == 'JHI-Hv50k-2016-182700')),9] #438612841

hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-326205   5 195.5698 1.508282
#c5.loc201               5 201.0000 3.713137
#JHI-Hv50k-2016-335125   5 206.1578 2.130289
datposition[(which(datposition == 'JHI-Hv50k-2016-326205')),9] #519836337
datposition[(which(datposition == 'JHI-Hv50k-2016-330574')),9] #527791267
datposition[(which(datposition == 'JHI-Hv50k-2016-335125')),9] #536096044

#######HvDRR14#######
hvdrr <- datWidth$HvDRR14 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #76 inds and 94.1% genotyped
png('summary_Width14.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width14.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#BOPA2_12_31270   2 151 3.41 0.0333
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2), c(151), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 2@151.4   2 151.39     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1)) #adds if there are QTL effects
pt(abs(-4.123), 73, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width14.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr2
find.marker(hvdrrprob, 2, 151.39) #BOPA2_12_31270
png('effetplot_ch2_Width14.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "BOPA2_12_31270")
dev.off()
png('effetplotsca_ch2_Width14.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "BOPA2_12_31270") #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-113683   2 145.8905 1.803690
#BOPA2_12_31270          2 151.3886 3.413177
#JHI-Hv50k-2016-131327   2 197.5784 1.858816
datposition[(which(datposition == 'JHI-Hv50k-2016-113683')),9] #610874966
datposition[(which(datposition == 'BOPA2_12_31270')),9] #615216417
datposition[(which(datposition == 'JHI-Hv50k-2016-131327')),9] #646476109


#######HvDRR15#######
hvdrr <- datWidth$HvDRR15 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #47 inds and 94.5% genotyped
png('summary_Width15.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width15.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #no significant peaks
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no peak above the threshold

#######HvDRR16#######
hvdrr <- datWidth$HvDRR16 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #86 inds and 93.4% genotyped
png('summary_Width16.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width16.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-129807   2 264 7.77    0
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2), c(264), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-383725   6 114 4.52 0.0025
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 6), c(264, 114), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 2@264.2   2 264.19     2
#Q2 6@114.1   6 114.10     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2)) #adds if there are QTL effects
pt(abs(4.753), 81, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width16.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr2
find.marker(hvdrrprob, 2, 264.19) #JHI-Hv50k-2016-129807
png('effetplot_ch2_Width16.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-129807")
dev.off()
png('effetplotsca_ch2_Width16.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-129807") #Scatter plot
dev.off()

#effect plot marker chr6
find.marker(hvdrrprob, 6, 114.10) #JHI-Hv50k-2016-383725
png('effetplot_ch6_Width16.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-383725')
dev.off()
png('effetplotsca_ch6_Width16.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-383725') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-127383   2 260.0660 7.050639
#JHI-Hv50k-2016-129807   2 264.1915 9.712822
#JHI-Hv50k-2016-131048   2 268.8867 7.913705
datposition[(which(datposition == 'JHI-Hv50k-2016-127383')),9] #638507553
datposition[(which(datposition == 'JHI-Hv50k-2016-129807')),9] #643472880
datposition[(which(datposition == 'JHI-Hv50k-2016-131048')),9] #645931108

hvdrrCI <- lodint(hvdrr_rqtl, chr = 6, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-381532   6  91.49433 2.972981
#JHI-Hv50k-2016-383725   6 114.09922 4.496184
#JHI-Hv50k-2016-401443   6 181.89793 2.374837
datposition[(which(datposition == 'JHI-Hv50k-2016-381830')),9] #30272196
datposition[(which(datposition == 'JHI-Hv50k-2016-383725')),9] #35724828
datposition[(which(datposition == 'JHI-Hv50k-2016-401443')),9] #384025372
hvdrr$geno$`6`

#######HvDRR17#######
hvdrr <- datWidth$HvDRR17 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrrpheno <- hvdrr$pheno
hvdrr$pheno[92, 1] <- NA
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #112 inds and 93.2% genotyped
png('summary_Width17.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width17.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr6
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-386716   6 79.5 7.15    0

#JHI-Hv50k-2016-386716   6 79.5 7.15    0
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(6), c(79.5), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-77411   2 74.9 3.09 0.0467
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 6), c(74.9, 79.5), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-79235   2 91.2 3.73 0.0105
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 2, 6), c(74.9, 91.2, 79.5), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 2@74.9   2 74.928     2
#Q2 2@91.2   2 91.199     2
#Q3 6@79.5   6 79.453     2
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1 + Q2 + Q3)) #adds if there are QTL effects
pt(abs(-4.232), 106, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width17.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr6
find.marker(hvdrrprob, 6, 79.453) #JHI-Hv50k-2016-386716
png('effetplot_ch6_Width17.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-386716')
dev.off()
png('effetplotsca_ch6_Width17.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-386716') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 6, qtl.index = 3, expandtomarkers = TRUE)
hvdrrCI
#Q1
#JHI-Hv50k-2016-77146   2 73.47605 4.440726
#JHI-Hv50k-2016-77411   2 74.92822 6.786169
#SCRI_RS_182408         2 77.35658 5.025526
datposition[(which(datposition == 'JHI-Hv50k-2016-77146')),9] #37885471
datposition[(which(datposition == 'JHI-Hv50k-2016-77411')),9] #38921714
datposition[(which(datposition == 'SCRI_RS_182408')),9] #39806820
#Q2
#BOPA1_3709-716         2 82.77277 1.791628
#JHI-Hv50k-2016-79235   2 91.19880 3.733645
#JHI-Hv50k-2016-79711   2 92.59751 1.860826
datposition[(which(datposition == 'BOPA1_3709-716')),9] #41368107
datposition[(which(datposition == 'JHI-Hv50k-2016-79228')),9] #46146872
datposition[(which(datposition == 'JHI-Hv50k-2016-79711')),9] #48453575
#Q3
#JHI-Hv50k-2016-384956   6 78.34840 7.370811
#JHI-Hv50k-2016-386716   6 79.45288 9.065601
#JHI-Hv50k-2016-390826   6 97.49002 7.513090
datposition[(which(datposition == 'JHI-Hv50k-2016-384956')),9] #40544757
datposition[(which(datposition == 'JHI-Hv50k-2016-386716')),9] #50074479
datposition[(which(datposition == 'JHI-Hv50k-2016-390826')),9] #117015712
hvdrr$geno$`6`
#######HvDRR18#######
hvdrr <- datWidth$HvDRR18 #extracting the single population
summary(hvdrr$pheno)
hvdrrpheno <- hvdrr$pheno
hvdrr$pheno[71, 1] <- NA
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
plot(hvdrr)
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #78 inds and 93.3% genotyped
png('summary_Width18.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width18.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr6
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-416642   6 89.2 3.38 0.022
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(6), c(89.2), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-308333   5 84.4 3.35 0.025
hvdrrqtl <- makeqtl(hvdrrprob, c(5, 6), c(84.4, 89.2), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 5@84.4   5 84.437     2
#Q2 6@96.0   6 96.000     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2)) #adds if there are QTL effects
pt(abs(-4.19), 76, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos18_Width.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 84.437) #JHI-Hv50k-2016-308333
png('effetplot_ch5_Width18.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-308333")
dev.off()
png('effetplotsca_ch5_Width18.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-308333") #Scatter plot
dev.off()
#effect plot marker chr6
find.marker(hvdrrprob, 6, 96) #JHI-Hv50k-2016-417206
png('effetplot_ch6_Width18.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-417206")
dev.off()
png('effetplotsca_ch6_Width18.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-417206") #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#BOPA1_7504-485          5  82.54658 1.969903
#JHI-Hv50k-2016-308333   5  84.43703 3.916562
#JHI-Hv50k-2016-310285   5 102.49252 2.289108
datposition[(which(datposition == 'BOPA1_7504-485')),9] #437264627
datposition[(which(datposition == 'JHI-Hv50k-2016-308333')),9] #440797169
datposition[(which(datposition == 'JHI-Hv50k-2016-310285')),9] #463049204
hvdrr$geno$`5`
hvdrrCI <- lodint(hvdrr_rqtl, chr = 6, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-416242   6  85.67668 3.371821
#c6.loc96                6  96.00000 5.085756
#SCRI_RS_95857           6 103.29936 3.033650
datposition[(which(datposition == 'JHI-Hv50k-2016-416242')),9] # 528697061
datposition[(which(datposition == 'JHI-Hv50k-2016-417206')),9] #534304885
datposition[(which(datposition == 'SCRI_RS_95857')),9] #542095015
datposition[(which(datposition == 'JHI-Hv50k-2016-420301')),9]
#######HvDRR19#######
hvdrr <- datWidth$HvDRR19 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #83 inds and 93.6% genotyped
png('summary_Width.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width19.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr5
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-287068   5  56 4.68 0.0025
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(5), c(56), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-28745   1 110 3.39 0.0257
hvdrrqtl <- makeqtl(hvdrrprob, c(1, 5), c(110, 56), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 1@110.5   1 110.456     2
#Q2  5@56.0   5  56.008     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2)) #adds if there are QTL effects
pt(abs(4.070), 78, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width19.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr1
find.marker(hvdrrprob, 1, 110.456) #
png('effetplot_ch1_Width19.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-28745")
dev.off()
png('effetplotsca_ch1_Width19.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-28745") #Scatter plot
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 56.008) #JHI-Hv50k-2016-287068
png('effetplot_ch5_Width19.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-287068')
dev.off()
png('effetplotsca_ch5_Width19.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-287068') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 1, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#SCRI_RS_114047         1  74.6073 1.703084
#JHI-Hv50k-2016-28745   1 110.4559 3.391320
#JHI-Hv50k-2016-34535   1 130.6576 1.748760
datposition[(which(datposition == 'SCRI_RS_114047')),9] #29745896
datposition[(which(datposition == 'JHI-Hv50k-2016-28745')),9] #381831869
datposition[(which(datposition == 'JHI-Hv50k-2016-34535')),9] #427348709

hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-286412   5 52.82293 4.575278
#JHI-Hv50k-2016-287068   5 56.00849 6.145094
#SCRI_RS_13395           5 64.38447 4.288382
datposition[(which(datposition == 'JHI-Hv50k-2016-286412')),9] #25133724
datposition[(which(datposition == 'JHI-Hv50k-2016-287068')),9] #27766729
datposition[(which(datposition == 'SCRI_RS_13395')),9] #324719977
hvdrr$geno$`5`

#######HvDRR20#######
hvdrr <- datWidth$HvDRR20 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) # 98 inds and 93.9% genotyped
png('summary_Width20.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width20.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr6 and 7
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-426707   6 161 6.83    0
#JHI-Hv50k-2016-486204   7 118 8.74    0
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(6, 7), c(161, 118), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 6@160.5   6 160.47     2
#Q2 7@118.0   7 118.00     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2)) #adds if there are QTL effects
pt(abs(5.347), 93, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width20.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr6
find.marker(hvdrrprob, 6, 160.47) #JHI-Hv50k-2016-426633
png('effetplot_ch6_Width20.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-426633")
dev.off()
png('effetplotsca_ch6_Width20.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-426633") #Scatter plot
dev.off()

#effect plot marker chr7
find.marker(hvdrrprob, 7, 118) #JHI-Hv50k-2016-486204
png('effetplot_ch7_Width20.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-486204')
dev.off()
png('effetplotsca_ch7_Width20.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-486204') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 6, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-424364   6 155.6959 3.817608
#JHI-Hv50k-2016-426633   6 160.4652 5.598483
#JHI-Hv50k-2016-427559   6 162.3319 4.077363
datposition[(which(datposition == 'JHI-Hv50k-2016-424364')),9] #554318838
datposition[(which(datposition == 'JHI-Hv50k-2016-426633')),9] #558656624
datposition[(which(datposition == 'JHI-Hv50k-2016-427559')),9] #561778610

hvdrrCI <- lodint(hvdrr_rqtl, chr = 7, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-474021   7 106.9027 6.035435
#c7.loc118               7 118.0000 7.540476
#SCRI_RS_4562            7 142.1958 5.877797
datposition[(which(datposition == 'JHI-Hv50k-2016-474021')),9] #110765114
datposition[(which(datposition == 'JHI-Hv50k-2016-486204')),9] #428068719
datposition[(which(datposition == 'SCRI_RS_4562')),9] #551468516

hvdrrCI <- lodint(hvdrr_rqtl, chr = , qtl.index = 3, expandtomarkers = TRUE)
hvdrrCI
#
datposition[(which(datposition == '')),9] #
datposition[(which(datposition == '')),9] #
datposition[(which(datposition == '')),9] #
#######HvDRR21#######
hvdrr <- datWidth$HvDRR21 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrrpheno <- hvdrr$pheno
hvdrr$pheno[40, 1] <- NA
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #69 inds and 93.9% genotyped
png('summary_Width21.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width21.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr3
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-23499   1 74.8 3.52 0.02325
#c3.loc82               3 82.0 6.19 0.00025
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(1, 3), c(74.8, 82), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 1@77.0   1  77     2
#Q2 3@86.0   3  86     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2)) #adds if there are QTL effects
pt(abs(-4.624), 64, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width21.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr1
find.marker(hvdrrprob, 1, 77) #BOPA1_7090-1260
png('effetplot_ch1_Width21.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "BOPA1_7090-1260")
dev.off()
png('effetplotsca_ch1_Width21.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "BOPA1_7090-1260") #Scatter plot
dev.off()

#effect plot marker chr3
find.marker(hvdrrprob, 3, 86) #JHI-Hv50k-2016-184503
png('effetplot_ch3_Width21.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-184503')
dev.off()
png('effetplotsca_ch3_Width21.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-184503') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 1, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#BOPA2_12_10314           1 66.72813 2.510091
#c1.loc77                 1 77.00000 4.204988
#BOPA1_ABC10636-1-4-285   1 84.21575 2.003348
datposition[(which(datposition == 'BOPA2_12_10314')),9] #43706821
datposition[(which(datposition == 'BOPA1_7090-1260')),9] #308262183
datposition[(which(datposition == 'BOPA1_ABC10636-1-4-285')),9] #370302846

hvdrrCI <- lodint(hvdrr_rqtl, chr = 3, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-182720   3 80.63626 4.358956
#c3.loc86                3 86.00000 6.126046
#JHI-Hv50k-2016-186165   3 88.04747 4.371497
datposition[(which(datposition == 'JHI-Hv50k-2016-182720')),9] #438802542
datposition[(which(datposition == 'JHI-Hv50k-2016-184503')),9] #454933755
datposition[(which(datposition == 'JHI-Hv50k-2016-186165')),9] #467665750


#######HvDRR22#######
hvdrr <- datWidth$HvDRR22 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #90 inds and 93.7% genotyped
png('summary_Width22.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width22.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #NO SIGNIFICANT PEAKS
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#There were no LOD peaks above the threshold.

#######HvDRR23#######
hvdrr <- datWidth$HvDRR23 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #83 inds and 93.5% genotyped
png('summary_Width23.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width23.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr7
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#c7.loc155   7 155 13.9    0
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(7), c(155), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 7@155.0   7 155     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1)) #adds if there are QTL effects
pt(abs(9.689), 80, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width23.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr7
find.marker(hvdrrprob, 7, 155) #BOPA2_12_20685
png('effetplot_ch7_Width23.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "BOPA2_12_20685")
dev.off()
png('effetplotsca_ch7_Width23.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "BOPA2_12_20685") #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 7, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-491256   7 150.2046 11.56790
#c7.loc155               7 155.0000 13.87163
#BOPA2_12_11437          7 158.1702 11.21532
datposition[(which(datposition == 'JHI-Hv50k-2016-491256')),9] #521128569
datposition[(which(datposition == 'BOPA2_12_20685')),9] #532125027
datposition[(which(datposition == 'BOPA2_12_11437')),9] #536163387

#######HvDRR24#######
hvdrr <- datWidth$HvDRR24 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrrpheno <- hvdrr$pheno
hvdrr$pheno[49, 1] <- NA
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #59 inds and 93.8% genotyped
png('summary_Width24.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width24.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-33097    1 101 4.35 0.0050
#JHI-Hv50k-2016-426471   6 172 3.53 0.0285
#JHI-Hv50k-2016-491439   7 164 5.09 0.0010

#JHI-Hv50k-2016-7757     1  23.1 3.84 0.0175
#JHI-Hv50k-2016-426471   6 172.2 3.69 0.0210
#c7.loc170               7 170.0 3.46 0.0302
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(1, 6, 7), c(23.1, 172.2, 170.0), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)

#hvdrrqtl <- makeqtl(hvdrrprob, c(1, 6, 7), c(23.1, 172.2, 164), what = "prob")
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 1@100.6   1 100.60     2
#Q2 6@169.0   6 169.00     2
#Q3 7@163.9   7 163.88     2
#Q1  1@23.1   1  23.126     2
#Q2 6@172.2   6 172.205     2
#Q3 7@182.0   7 182.000     2
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2+Q3)) #adds if there are QTL effects
pt(abs(4.336), 52, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width24.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr1
find.marker(hvdrrprob, 1, 23.126) #JHI-Hv50k-2016-33097 #JHI-Hv50k-2016-7757
#find.marker(hvdrrprob, 1, 23.1)
png('effetplot_ch1_Width24.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-33097")
dev.off()
png('effetplotsca_ch1_Width24.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-33097") #Scatter plot
dev.off()

#effect plot marker chr6
find.marker(hvdrrprob, 6, 172.205) #JHI-Hv50k-2016-425342 #JHI-Hv50k-2016-426471
png('effetplot_ch6_Width24.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-425342')
dev.off()
png('effetplotsca_ch6_Width24.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-425342') #Scatter plot
dev.off()

#effect plot marker chr7
find.marker(hvdrrprob, 7, 182) #JHI-Hv50k-2016-491472 #JHI-Hv50k-2016-493417 #JHI-Hv50k-2016-493417
png('effetplot_ch7_Width24.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-491472')
dev.off()
png('effetplotsca_ch7_Width.png24', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-491472') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 1, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-29153   1  88.01407 4.213469
#JHI-Hv50k-2016-33097   1 100.59695 6.088148
#JHI-Hv50k-2016-34594   1 105.94633 2.868839

#JHI-Hv50k-2016-170     1   4.801071 1.429461
#JHI-Hv50k-2016-7757    1  23.126441 4.303081
#JHI-Hv50k-2016-33109   1 101.298206 2.563810
datposition[(which(datposition == 'JHI-Hv50k-2016-170')),9] #384541598 #174905
datposition[(which(datposition == 'JHI-Hv50k-2016-7757')),9] #415112798 #6428040
datposition[(which(datposition == 'JHI-Hv50k-2016-33109')),9] #427830820 #415358170

hvdrrCI <- lodint(hvdrr_rqtl, chr = 6, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-423961   6 159.1876 1.856817
#c6.loc169               6 169.0000 3.766430
#JHI-Hv50k-2016-428430   6 180.0904 2.180511

#JHI-Hv50k-2016-424045   6 161.3444 1.930830
#JHI-Hv50k-2016-426471   6 172.2050 4.290478
#JHI-Hv50k-2016-429246   6 182.3696 2.554586
datposition[(which(datposition == 'JHI-Hv50k-2016-424045')),9] #553316571 #553429109
datposition[(which(datposition == 'JHI-Hv50k-2016-426471')),9] #556245192 #558316219
datposition[(which(datposition == 'JHI-Hv50k-2016-429246')),9] #563722787 #565768882

hvdrrCI <- lodint(hvdrr_rqtl, chr = 7, qtl.index = 3, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-489914   7 154.4482 2.401091
#JHI-Hv50k-2016-491472   7 163.8849 4.413913
#BOPA1_4791-1541         7 172.2029 2.115472

#BOPA2_12_30026          7 165.9623 1.420911
#c7.loc182               7 182.0000 3.221500
#JHI-Hv50k-2016-497609   7 192.8382 0.812670
datposition[(which(datposition == 'BOPA2_12_30026')),9] #504893609 #542585738
datposition[(which(datposition == 'JHI-Hv50k-2016-493417')),9] #529112900 #573028500
datposition[(which(datposition == 'JHI-Hv50k-2016-497609')),9] #571922809 #587308534

#######HvDRR25#######
hvdrr <- datWidth$HvDRR25 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrrpheno <- hvdrr$pheno
hvdrr$pheno[64, 1] <- NA
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #74 inds and 93,5% genotyped
png('summary_Width25.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width25.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2 and 7
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-72378   2  38.4 3.32 0.0435
#BOPA2_12_30368         7 214.0 4.52 0.0025

#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2,7), c(38.4,214.0), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-101118   2 112 3.39 0.0393
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 2, 7), c(38.4, 112, 214.0), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-159468   3 35.3 4.15 0.00825
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 2, 3, 7), c(38.4, 112, 35.3, 214.0), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#c4.loc36   4  36 3.51 0.0285
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 2, 3, 4, 7), c(38.4, 112, 35.3, 36, 214.0), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1  2@40.8   2  40.826     2
#Q2 2@112.2   2 112.182     2
#Q3  3@35.3   3  35.275     2
#Q4  4@37.0   4  37.000     2
#Q5 7@213.4   7 213.444     2

#Q1 7@214.0   7 214     2
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2+Q3+Q4+Q5)) #adds if there are QTL effects
pt(abs(-4.581), 63, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width25.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr2a
find.marker(hvdrrprob, 2, 40.826) #JHI-Hv50k-2016-73697
png('effetplot_ch2a_Width25.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-73697")
dev.off()
png('effetplotsca_ch2a_Width25.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-73697") #Scatter plot
dev.off()

#effect plot marker chr2b
find.marker(hvdrrprob, 2, 112.182) #JHI-Hv50k-2016-101118
png('effetplot_ch2b_Width25.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-101118')
dev.off()
png('effetplotsca_ch2b_Width25.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-101118') #Scatter plot
dev.off()

#effect plot marker chr3
find.marker(hvdrrprob, 3, 35.275) #JHI-Hv50k-2016-159441
png('effetplot_ch3_Width25.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-159441')
dev.off()
png('effetplotsca_ch3_Width25.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-159441') #Scatter plot
dev.off()

#effect plot marker chr4
find.marker(hvdrrprob, 4, 37.000) #JHI-Hv50k-2016-230655
png('effetplot_ch4_Width25.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-230655')
dev.off()
png('effetplotsca_ch4_Width25.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-230655') #Scatter plot
dev.off()

#effect plot marker chr7
find.marker(hvdrrprob, 7, 213.444) #JHI-Hv50k-2016-501148
png('effetplot_ch7_Width25.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-501148')
dev.off()
png('effetplotsca_ch7_Width25.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-501148') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#SCRI_RS_66401          2 28.78644 3.709067
#JHI-Hv50k-2016-73697   2 40.82567 7.730218
#JHI-Hv50k-2016-73822   2 47.02424 4.599390
datposition[(which(datposition == 'SCRI_RS_66401')),9] #19007448
datposition[(which(datposition == 'JHI-Hv50k-2016-73697')),9] #24211311
datposition[(which(datposition == 'JHI-Hv50k-2016-73822')),9] #24377464

hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-99471    2 106.4764 3.811385
#JHI-Hv50k-2016-101118   2 112.1825 6.261115
#JHI-Hv50k-2016-103059   2 115.7502 4.411259
datposition[(which(datposition == 'JHI-Hv50k-2016-99471')),9] #520894241
datposition[(which(datposition == 'JHI-Hv50k-2016-101118')),9] #547304012
datposition[(which(datposition == 'JHI-Hv50k-2016-103059')),9] #558758364

hvdrrCI <- lodint(hvdrr_rqtl, chr = 3, qtl.index = 3, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-157771   3 27.49980 4.424608
#JHI-Hv50k-2016-159468   3 35.27468 6.188982
#JHI-Hv50k-2016-159793   3 37.66443 4.456603
datposition[(which(datposition == 'JHI-Hv50k-2016-157771')),9] #15035502
datposition[(which(datposition == 'JHI-Hv50k-2016-159441')),9] #18887837
datposition[(which(datposition == 'JHI-Hv50k-2016-159793')),9] #19978985

hvdrrCI <- lodint(hvdrr_rqtl, chr = 4, qtl.index = 4, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-230329   4 29.63876 2.331707
#c4.loc37                4 37.00000 4.322402
#JHI-Hv50k-2016-231733   4 44.78992 2.716380
datposition[(which(datposition == 'JHI-Hv50k-2016-230329')),9] #11691148
datposition[(which(datposition == 'JHI-Hv50k-2016-230655')),9] #13737339
datposition[(which(datposition == 'JHI-Hv50k-2016-231733')),9] #19994864

hvdrrCI <- lodint(hvdrr_rqtl, chr = 7, qtl.index = 5, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-498414   7 200.1788 4.059011
#JHI-Hv50k-2016-501148   7 213.4437 6.244974
#JHI-Hv50k-2016-503621   7 219.4232 3.759620
datposition[(which(datposition == 'JHI-Hv50k-2016-498414')),9] #589772676
datposition[(which(datposition == 'JHI-Hv50k-2016-501148')),9] #598281223
datposition[(which(datposition == 'JHI-Hv50k-2016-503621')),9] #603011955

#######HvDRR26#######
hvdrr <- datWidth$HvDRR26 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #54 inds and 93,9% genotyped
png('summary_Width26.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width26.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #no signficant peaks
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#    There were no LOD peaks above the threshold.

#######HvDRR27#######
hvdrr <- datWidth$HvDRR27 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #93 inds and 93.5% genotyped
png('summary_Width27.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width27.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2 and 4
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-77256    2  68.3 5.19 0.00025
#c3.loc235               3 235.0 3.75 0.01375
#JHI-Hv50k-2016-231492   4  33.5 3.29 0.03525
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 3, 4), c(68.3, 235.0, 33.5), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-95965    2 119 3.67 0.0165
#JHI-Hv50k-2016-163953   3  97 3.51 0.0222
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 2, 3, 3, 4), c(68.3, 119, 97, 235.0, 33.5), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1  2@58.0   2  57.985     2
#Q2 2@117.0   2 117.016     2
#Q3  3@74.0   3  74.000     2
#Q4 3@229.0   3 229.000     2
#Q5  4@33.5   4  33.512     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2+Q3+Q4+Q5)) #adds if there are QTL effects
pt(abs(4.532), 82, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width27.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr2a
find.marker(hvdrrprob, 2, 57.985) #JHI-Hv50k-2016-75250
png('effetplot_ch2a_Width27.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-75250")
dev.off()
png('effetplotsca_ch2a_Width27.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-75250") #Scatter plot
dev.off()

#effect plot marker chr2b
find.marker(hvdrrprob, 2, 117.016) #JHI-Hv50k-2016-93359
png('effetplot_ch2b_Width27.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-93359')
dev.off()
png('effetplotsca_ch2b_Width27.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-93359') #Scatter plot
dev.off()

#effect plot marker chr3a
find.marker(hvdrrprob, 3, 74.000) #BOPA1_15141-257
png('effetplot_ch3a_Width27.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'BOPA1_15141-257')
dev.off()
png('effetplotsca_ch3a_Width27.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'BOPA1_15141-257') #Scatter plot
dev.off()

#effect plot marker chr3b
find.marker(hvdrrprob, 4, 33.512) #JHI-Hv50k-2016-231492
png('effetplot_ch4_Width27.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-231492')
dev.off()
png('effetplotsca_ch4_Width27.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-231492') #Scatter plot
dev.off()

#effect plot marker chr4
find.marker(hvdrrprob, 4, 229.000) #JHI-Hv50k-2016-209865
png('effetplot_ch3b_Width27.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-209865')
dev.off()
png('effetplotsca_ch3b_Width27.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-209865') #Scatter plot
dev.off()
#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#SCRI_RS_142851         2 8.696783e-05 2.609094
#JHI-Hv50k-2016-75250   2 5.798512e+01 4.386543
#JHI-Hv50k-2016-77411   2 6.968696e+01 2.372252
datposition[(which(datposition == 'SCRI_RS_142851')),9] #880829
datposition[(which(datposition == 'JHI-Hv50k-2016-75250')),9] #28816411
datposition[(which(datposition == 'JHI-Hv50k-2016-77411')),9] #38921714
hvdrr$geno$`2`

hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-87492   2 113.1019 4.407590
#JHI-Hv50k-2016-93359   2 117.0165 6.060633
#SCRI_RS_4802           2 123.1315 3.707340
datposition[(which(datposition == 'JHI-Hv50k-2016-87492')),9] #102484526
datposition[(which(datposition == 'JHI-Hv50k-2016-93359')),9] #394287896
datposition[(which(datposition == 'JHI-Hv50k-2016-99022')),9] #515681394 129.97cM

hvdrrCI <- lodint(hvdrr_rqtl, chr = 3, qtl.index = 3, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-159228   3  47.33977 3.234864
#c3.loc74                3  74.00000 4.946855
#JHI-Hv50k-2016-164723   3 100.11208 2.969163
datposition[(which(datposition == 'JHI-Hv50k-2016-159228')),9] #18392788
datposition[(which(datposition == 'BOPA1_15141-257')),9] #32422295
datposition[(which(datposition == 'JHI-Hv50k-2016-164723')),9] #45111515

hvdrrCI <- lodint(hvdrr_rqtl, chr = 3, qtl.index = 4, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-207551   3 220.4236 2.623965
#c3.loc229               3 229.0000 4.279878
#JHI-Hv50k-2016-212423   3 236.8873 2.169905
datposition[(which(datposition == 'JHI-Hv50k-2016-207551')),9] #582546764
datposition[(which(datposition == 'JHI-Hv50k-2016-209865')),9] #590032574
datposition[(which(datposition == 'JHI-Hv50k-2016-212423')),9] #594996922

hvdrrCI <- lodint(hvdrr_rqtl, chr = 4, qtl.index = 5, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-231183   4 32.57075 4.492105
#JHI-Hv50k-2016-231492   4 33.51234 6.095483
#JHI-Hv50k-2016-231961   4 37.51181 4.130193
datposition[(which(datposition == 'JHI-Hv50k-2016-231183')),9] #17611398
datposition[(which(datposition == 'JHI-Hv50k-2016-231492')),9] #18251048
datposition[(which(datposition == 'JHI-Hv50k-2016-231961')),9] #21238643

#######HvDRR28#######
hvdrr <- datWidth$HvDRR28 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #81 inds and 94% genotyped
png('summary_Width28.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width28.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-91390   2 113 4.31 0.003
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2), c(113), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-231202   4 33.8 4.36 0.003
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 4), c(113, 33.8), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 2@113.4   2 113.364     2
#Q2  4@33.8   4  33.828     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2)) #adds if there are QTL effects
pt(abs(4.627), 76, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width28.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr2
find.marker(hvdrrprob, 2, 113.364) #JHI-Hv50k-2016-90638
png('effetplot_ch2_Width28.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-90638")
dev.off()
png('effetplotsca_ch2_Width28.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-90638") #Scatter plot
dev.off()

#effect plot marker chr4
find.marker(hvdrrprob, 4, 33.828) #JHI-Hv50k-2016-231202
png('effetplot_ch4_Width28.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-231202')
dev.off()
png('effetplotsca_ch4_Width28.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-231202') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-87559   2 110.8022 4.312946
#JHI-Hv50k-2016-91390   2 113.3645 6.848530
#JHI-Hv50k-2016-95965   2 116.5968 4.647870
datposition[(which(datposition == 'JHI-Hv50k-2016-87559')),9] #103463080
datposition[(which(datposition == 'JHI-Hv50k-2016-90638')),9] #165428931
datposition[(which(datposition == 'JHI-Hv50k-2016-95965')),9] #469888490

hvdrrCI <- lodint(hvdrr_rqtl, chr = 4, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#SCRI_RS_202425          4 31.01320 2.424882
#JHI-Hv50k-2016-231202   4 33.82770 4.265221
#JHI-Hv50k-2016-231851   4 38.46079 2.088511
datposition[(which(datposition == 'SCRI_RS_202425')),9] #15560515
datposition[(which(datposition == 'JHI-Hv50k-2016-231202')),9] #17720939
datposition[(which(datposition == 'JHI-Hv50k-2016-231851')),9] #20592280

#######HvDRR29#######
hvdrr <- datWidth$HvDRR29 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #110 inds and 93.8% genotyped
png('summary_Width29.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width29.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#BOPA2_12_30897   2 146 14.9    0
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2), c(146), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-230927   4  39.2 4.65 0.00150
#c5.loc206               5 206.0 4.84 0.00125
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 4, 5), c(146, 39.2, 206.0), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-408360   6 102 5.18 5e-04
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 4, 5, 6), c(146, 39.2, 206.0, 102), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#SCRI_RS_160461   4 151 3.53 0.024
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 4, 4, 5, 6), c(146, 39.2, 151, 206.0, 102), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 2@146.5   2 146.54     2
#Q2  4@36.3   4  36.33     2
#Q3 4@150.8   4 150.84     2
#Q4 5@207.1   5 207.10     2
#Q5 6@101.8   6 101.81     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2+Q3+Q4+Q5)) #adds if there are QTL effects
pt(abs( 4.731), 99, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width29.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr2
find.marker(hvdrrprob, 2, 146.54) #JHI-Hv50k-2016-107351
png('effetplot_ch2_Width29.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-107351")
dev.off()
png('effetplotsca_ch2_Width29.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-107351") #Scatter plot
dev.off()

#effect plot marker chr4a
find.marker(hvdrrprob, 4, 36.33) #SCRI_RS_98443
png('effetplot_ch4a_Width29.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'SCRI_RS_98443')
dev.off()
png('effetplotsca_ch4a_Width29.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'SCRI_RS_98443') #Scatter plot
dev.off()

#effect plot marker chr4b
find.marker(hvdrrprob, 4, 150.84) #SCRI_RS_160461
png('effetplot_ch4b_Width29.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'SCRI_RS_160461')
dev.off()
png('effetplotsca_ch4b_Width29.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'SCRI_RS_160461') #Scatter plot
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 207.10) #JHI-Hv50k-2016-327190
png('effetplot_ch5_Width29.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-327190')
dev.off()
png('effetplotsca_ch5_Width29.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-327190') #Scatter plot
dev.off()

#effect plot marker chr6
find.marker(hvdrrprob, 6, 101.81) #JHI-Hv50k-2016-408298
png('effetplot_ch6_Width29.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-408298')
dev.off()
png('effetplotsca_ch6_Width29.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-408298') #Scatter plot
dev.off()
#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-106820   2 144.7666 21.14825
#JHI-Hv50k-2016-107351   2 146.5440 25.51781
#JHI-Hv50k-2016-108149   2 151.0455 16.72863
datposition[(which(datposition == 'JHI-Hv50k-2016-106820')),9] #579356193
datposition[(which(datposition == 'JHI-Hv50k-2016-107351')),9] #580550249
datposition[(which(datposition == 'JHI-Hv50k-2016-108149')),9] #583656135

hvdrrCI <- lodint(hvdrr_rqtl, chr = 4, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-230743   4 35.53144 6.932658
#SCRI_RS_98443           4 36.32958 8.612280
#JHI-Hv50k-2016-231184   4 40.79691 6.368356
datposition[(which(datposition == 'JHI-Hv50k-2016-230743')),9] #14129869
datposition[(which(datposition == 'SCRI_RS_98443')),9] #14583570
datposition[(which(datposition == 'JHI-Hv50k-2016-231184')),9] #17610769

hvdrrCI <- lodint(hvdrr_rqtl, chr = 4, qtl.index = 3, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-265904   4 145.6822 2.692773
#SCRI_RS_160461          4 150.8381 4.656257
#JHI-Hv50k-2016-268238   4 152.1080 2.981206
datposition[(which(datposition == 'JHI-Hv50k-2016-265904')),9] #601981565
datposition[(which(datposition == 'SCRI_RS_160461')),9] #604450209
datposition[(which(datposition == 'JHI-Hv50k-2016-268238')),9] #606352633

hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 4, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-318391   5 167.3430 3.485001
#JHI-Hv50k-2016-327190   5 207.1043 5.290078
#JHI-Hv50k-2016-328212   5 210.2117 2.522959
datposition[(which(datposition == 'JHI-Hv50k-2016-318391')),9] #499014161
datposition[(which(datposition == 'JHI-Hv50k-2016-327190')),9] #522269695
datposition[(which(datposition == 'JHI-Hv50k-2016-328212')),9] #524330188

hvdrrCI <- lodint(hvdrr_rqtl, chr = 6, qtl.index = 5, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-404600   6  96.98611 4.709100
#JHI-Hv50k-2016-408298   6 101.80997 6.550508
#JHI-Hv50k-2016-411991   6 105.28111 4.208065
datposition[(which(datposition == 'JHI-Hv50k-2016-404600')),9] #414346286
datposition[(which(datposition == 'JHI-Hv50k-2016-408298')),9] #470278903
datposition[(which(datposition == 'JHI-Hv50k-2016-411991')),9] #499718548

#######HvDRR30#######
hvdrr <- datWidth$HvDRR30 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #123 inds and 93.8% genotyped
png('summary_Width30.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width30.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#c2.loc188   2 188 24.2    0
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2), c(188), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 2@188.0   2 188     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1)) #adds if there are QTL effects
pt(abs(13.37), 120, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width30.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr2
find.marker(hvdrrprob, 2, 188) #JHI-Hv50k-2016-107351
png('effetplot_ch2_Width30.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-107351")
dev.off()
png('effetplotsca_ch2_Width30.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-107351") #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#SCRI_RS_221795   2 184.2607 22.11927
#c2.loc188        2 188.0000 24.22559
#SCRI_RS_88704    2 189.0364 21.09529
datposition[(which(datposition == 'SCRI_RS_221795')),9] #576550511
datposition[(which(datposition == 'JHI-Hv50k-2016-107351')),9] #580550249
datposition[(which(datposition == 'SCRI_RS_88704')),9] #582465837


#######HvDRR31#######
hvdrr <- datWidth$HvDRR31 #extracting the single population
testmap <- hvdrr$geno$`2`
testmarkername <- as.data.frame(testmap$map)
which(row.names(testmarkername) == 'JHI-Hv50k-2016-109469')
marker2drop <- markernames(hvdrr, chr = 2)[228:230]
testhvdrr <- drop.markers(hvdrr, marker2drop)
testhvdrr$geno$`2`
testhvdrrmap <- est.map(testhvdrr, error.prob = 0.001)
testhvdrrfinal <- replace.map(testhvdrr, testhvdrrmap)
plot(testhvdrrfinal)
summary(testhvdrrfinal$pheno)
hvdrrpheno <- testhvdrrfinal$pheno
testhvdrrfinal$pheno[76, 1] <- NA
testhvdrrfinal <- subset(testhvdrrfinal, ind=!is.na(pull.pheno(testhvdrrfinal,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
testhvdrrfinal$geno$'1'$data <- apply(testhvdrrfinal$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
testhvdrrfinal$geno$'2'$data <- apply(testhvdrrfinal$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
testhvdrrfinal$geno$'3'$data <- apply(testhvdrrfinal$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
testhvdrrfinal$geno$'4'$data <- apply(testhvdrrfinal$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
testhvdrrfinal$geno$'5'$data <- apply(testhvdrrfinal$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
testhvdrrfinal$geno$'6'$data <- apply(testhvdrrfinal$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
testhvdrrfinal$geno$'7'$data <- apply(testhvdrrfinal$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
testhvdrrfinal$geno$'1'$data <- apply(testhvdrrfinal$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
testhvdrrfinal$geno$'2'$data <- apply(testhvdrrfinal$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
testhvdrrfinal$geno$'3'$data <- apply(testhvdrrfinal$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
testhvdrrfinal$geno$'4'$data <- apply(testhvdrrfinal$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
testhvdrrfinal$geno$'5'$data <- apply(testhvdrrfinal$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
testhvdrrfinal$geno$'6'$data <- apply(testhvdrrfinal$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
testhvdrrfinal$geno$'7'$data <- apply(testhvdrrfinal$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(testhvdrrfinal)[1] <- 'riself'
summary(testhvdrrfinal) #128 inds and 94.4% genotyped
png('summary_Width31.png', width = 6, height = 6, units = 'in', res = 300)
plot(testhvdrrfinal)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(testhvdrrfinal, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width31.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-108254   2 219 3.55 0.0235
#c5.loc55                5  55 3.26 0.0490
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 5), c(219, 55), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-257076   4 146 3.2 0.0465
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 4, 5), c(219, 146, 55), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 2@217.0   2 217.00     2
#Q2 4@140.0   4 140.00     2
#Q3 5@151.3   5 151.33     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2+Q3)) #adds if there are QTL effects
pt(abs(-3.935), 121, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width31.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr2
find.marker(hvdrrprob, 2, 217) #JHI-Hv50k-2016-108928
png('effetplot_ch2_Width31.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-108928")
dev.off()
png('effetplotsca_ch2_Width31.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-108928") #Scatter plot
dev.off()

#effect plot marker chr4
find.marker(hvdrrprob, 4, 140) #JHI-Hv50k-2016-255545
png('effetplot_ch4_Width31.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-255545')
dev.off()
png('effetplotsca_ch4_Width31.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-255545') #Scatter plot
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 151.33) #JHI-Hv50k-2016-316287
png('effetplot_ch5_Width31.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-316287')
dev.off()
png('effetplotsca_ch5_Width31.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-316287') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#SCRI_RS_220533          2 209.7944 3.033061
#c2.loc217               2 217.0000 4.650897
#JHI-Hv50k-2016-108103   2 224.3177 2.239554
datposition[(which(datposition == 'SCRI_RS_220533')),9] #568573773
datposition[(which(datposition == 'JHI-Hv50k-2016-108928')),9] #587426328
datposition[(which(datposition == 'JHI-Hv50k-2016-108103')),9] #583438184

hvdrrCI <- lodint(hvdrr_rqtl, chr = 4, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-254984   4 136.8000 2.916937
#c4.loc140               4 140.0000 4.720757
#SCRI_RS_157650          4 147.0286 3.211396
datposition[(which(datposition == 'JHI-Hv50k-2016-254984')),9] #543746589
datposition[(which(datposition == 'JHI-Hv50k-2016-255545')),9] #549656721
datposition[(which(datposition == 'SCRI_RS_157650')),9] #558063586

hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 3, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-279669   5  48.84688 1.519097
#JHI-Hv50k-2016-316287   5 151.33014 3.271150
#BOPA2_12_30580          5 213.19577 1.445182
datposition[(which(datposition == 'JHI-Hv50k-2016-279669')),9] #5735259
datposition[(which(datposition == 'JHI-Hv50k-2016-316287')),9] #493026387
datposition[(which(datposition == 'BOPA2_12_30580')),9] #553460069

#######HvDRR32#######
hvdrr <- datWidth$HvDRR32 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #37 inds and 94.9% genotyped
png('summary_Width32.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width32.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #no significant peaks
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#There were no LOD peaks above the threshold.
#######HvDRR33#######
hvdrr <- datWidth$HvDRR33 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #79 inds and 93.2% genotyped
png('summary_Width33.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width33.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#BOPA2_12_30897   2 210 15.7    0
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2), c(210), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#BOPA2_12_31312   5 114 3.81 0.0147
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 5), c(210, 114), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-28380   1 123 3.28 0.0408
hvdrrqtl <- makeqtl(hvdrrprob, c(1, 2, 5), c(123, 210, 114), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 1@122.6   1 122.61     2
#Q2 2@209.8   2 209.82     2
#Q3 5@113.9   5 113.92     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2+Q3)) #adds if there are QTL effects
pt(abs(3.994), 72, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width33.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr1
find.marker(hvdrrprob, 1, 122.61) #JHI-Hv50k-2016-28380
png('effetplot_ch1_Width33.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-28380")
dev.off()
png('effetplotsca_ch1_Width33.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-28380") #Scatter plot
dev.off()

#effect plot marker chr2
find.marker(hvdrrprob, 2, 209.82) #JHI-Hv50k-2016-107729
png('effetplot_ch2_Width33.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-107729')
dev.off()
png('effetplotsca_ch2_Width33.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-107729') #Scatter plot
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 113.92) #SCRI_RS_103840
png('effetplot_ch5_Width33.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'SCRI_RS_103840')
dev.off()
png('effetplotsca_ch5_Width33.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'SCRI_RS_103840') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 1, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-16799   1  71.57569 1.115347
#JHI-Hv50k-2016-28380   1 122.61355 3.308257
#JHI-Hv50k-2016-30859   1 133.33797 1.646465
datposition[(which(datposition == 'JHI-Hv50k-2016-16799')),9] #26574375
datposition[(which(datposition == 'JHI-Hv50k-2016-28380')),9] #379746038
datposition[(which(datposition == 'JHI-Hv50k-2016-30859')),9] #397332767

hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-106954   2 206.3718 14.47253
#BOPA2_12_30897          2 209.8152 18.03004
#JHI-Hv50k-2016-108843   2 211.4576 15.83624
datposition[(which(datposition == 'JHI-Hv50k-2016-106954')),9] #579459694
datposition[(which(datposition == 'JHI-Hv50k-2016-107729')),9] #581638922
datposition[(which(datposition == 'JHI-Hv50k-2016-108843')),9] #586582763

hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 3, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-288872   5 109.6857 2.551751
#BOPA2_12_31312          5 113.9210 4.208850
#JHI-Hv50k-2016-301004   5 132.5898 2.686237
datposition[(which(datposition == 'JHI-Hv50k-2016-288872')),9] #34625394
datposition[(which(datposition == 'SCRI_RS_103840')),9] #44888824
datposition[(which(datposition == 'JHI-Hv50k-2016-301004')),9] #343047446

#######HvDRR34#######
hvdrr <- datWidth$HvDRR34 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #34 inds and 92.2% genotyped
png('summary_Width34.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width34.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #no significant peaks
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#There were no LOD peaks above the threshold

#######HvDRR35#######
hvdrr <- datWidth$HvDRR35 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #85 inds and 94.6% genotyped
png('summary_Width35.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width35.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr4 and 6
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-241208   4 95.9 4.42 0.004
#c6.loc17                6 17.0 4.15 0.007
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(4, 6), c(95.9, 17), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 4@95.9   4 95.94     2
#Q2 6@17.0   6 17.00     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2)) #adds if there are QTL effects
pt(abs(-4.459), 80, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width35.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr4
find.marker(hvdrrprob, 4, 95.94) #JHI-Hv50k-2016-241208
png('effetplot_ch4_Width35.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-241208")
dev.off()
png('effetplotsca_ch4_Width35.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-241208") #Scatter plot
dev.off()

#effect plot marker chr6
find.marker(hvdrrprob, 6, 17) #
png('effetplot_ch6_Width35.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-373096')
dev.off()
png('effetplotsca_ch6_Width35.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-373096') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 4, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-232606   4  60.11477 0.7322323
#JHI-Hv50k-2016-241208   4  95.93954 4.2760704
#BOPA1_7069-1149         4 114.14754 2.4095186
datposition[(which(datposition == 'JHI-Hv50k-2016-232606')),9] #24338016
datposition[(which(datposition == 'JHI-Hv50k-2016-241208')),9] #313407598
datposition[(which(datposition == 'BOPA1_7069-1149')),9] #507395834

hvdrrCI <- lodint(hvdrr_rqtl, chr = 6, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#BOPA1_2057-412          6  8.022318 1.890893
#c6.loc17                6 17.000000 4.007957
#JHI-Hv50k-2016-374677   6 21.649562 2.030088
datposition[(which(datposition == 'BOPA2_12_30651')),9] #7703894 7.47cM
datposition[(which(datposition == 'JHI-Hv50k-2016-373096')),9] #10463350
datposition[(which(datposition == 'JHI-Hv50k-2016-374677')),9] #12797251
hvdrr$geno$`6`

#######HvDRR36#######
hvdrr <- datWidth$HvDRR36 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#hvdrrpheno <- hvdrr$pheno
#hvdrr$pheno[83, 1] <- NA
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #58 inds and 94% genotyped
png('summary_Width36.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width36.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #no significant peak
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
# There were no LOD peaks above the threshold

#######HvDRR37#######
hvdrr <- datWidth$HvDRR37 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #63 inds and 94% genotyped
png('summary_Width37.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width37.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #no significant peaks
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
# There were no LOD peaks above the threshold.

#######HvDRR38#######
hvdrr <- datWidth$HvDRR38 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #71 inds and 94.5% genotyped
png('summary_Width38.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width38.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr4 and 5
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-227817   4 9.6 3.58 0.0275
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(4), c(9.6), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-296636   5 91.5 4.11 0.00925
hvdrrqtl <- makeqtl(hvdrrprob, c(4, 5), c(9.6, 91.5), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#c6.loc38   6  38 5.02 5e-04
hvdrrqtl <- makeqtl(hvdrrprob, c(4, 5, 6), c(9.6, 91.5, 38), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1  4@9.6   4  9.5958     2
#Q2 5@91.5   5 91.4584     2
#Q3 6@38.0   6 38.0000     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2+Q3)) #adds if there are QTL effects
pt(abs(-4.580), 64, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr4
find.marker(hvdrrprob, 4, 9.5958) #JHI-Hv50k-2016-227817
png('effetplot_ch4_Width38.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-227817")
dev.off()
png('effetplotsca_ch4_Width38.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-227817") #Scatter plot
dev.off()

#effect plot marker chr
find.marker(hvdrrprob, 5, 91.4584) #SCRI_RS_208177
png('effetplot_ch5_Width38.png', wiSCRI_RS_20817dth = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'SCRI_RS_208177')
dev.off()
png('effetplotsca_ch5_Width38.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'SCRI_RS_208177') #Scatter plot
dev.off()

#effect plot marker chr6
find.marker(hvdrrprob, 6, 38) #JHI-Hv50k-2016-376938
png('effetplot_ch6_Width38.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-376938')
dev.off()
png('effetplotsca_ch6_Width38.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'JHI-Hv50k-2016-376938') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 4, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-226108   4  0.000000 2.874918
#JHI-Hv50k-2016-227817   4  9.595772 4.198906
#JHI-Hv50k-2016-228512   4 22.915908 1.012265
datposition[(which(datposition == 'JHI-Hv50k-2016-226108')),9] #1670480
datposition[(which(datposition == 'JHI-Hv50k-2016-227817')),9] #5062419
datposition[(which(datposition == 'JHI-Hv50k-2016-228512')),9] #8100758

hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#SCRI_RS_170957          5 79.49554 5.333852
#JHI-Hv50k-2016-296636   5 91.45843 7.342655
#JHI-Hv50k-2016-298727   5 95.24981 4.845292
datposition[(which(datposition == 'SCRI_RS_170957')),9] #67641368
datposition[(which(datposition == 'SCRI_RS_208177')),9] #210965360
datposition[(which(datposition == 'JHI-Hv50k-2016-298727')),9] #314030053
hvdrr$geno$`5`
hvdrrCI <- lodint(hvdrr_rqtl, chr = 6, qtl.index = 3, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-375661   6 34.16766 3.030407
#c6.loc38                6 38.00000 5.016336
#SCRI_RS_189878          6 45.80602 2.823729
datposition[(which(datposition == 'JHI-Hv50k-2016-375661')),9] #14174001
datposition[(which(datposition == 'JHI-Hv50k-2016-376938')),9] #15517920
datposition[(which(datposition == 'SCRI_RS_189878')),9] #19163949

#######HvDRR39#######
hvdrr <- datWidth$HvDRR39 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #105 inds and 94.1% genotyped
png('summary_Width39.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width39.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-123897   2 226 5.08 0.00075
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2), c(226), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 2@225.8   2 225.81     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1)) #adds if there are QTL effects
pt(abs(-5.071), 102, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width39.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr2
find.marker(hvdrrprob, 2, 225.81) #JHI-Hv50k-2016-123572
png('effetplot_ch2_Width39.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-123572")
dev.off()
png('effetplotsca_ch2_Width39.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-123572") #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#BOPA1_5541-418          2 221.6221 3.399326
#JHI-Hv50k-2016-123897   2 225.8096 5.082067
#JHI-Hv50k-2016-129048   2 237.9591 3.557508
datposition[(which(datposition == 'BOPA1_5541-418')),9] #629624878
datposition[(which(datposition == 'JHI-Hv50k-2016-123572')),9] #631645054
datposition[(which(datposition == 'JHI-Hv50k-2016-129048')),9] #641669723

#######HvDRR40#######
hvdrr <- datWidth$HvDRR40 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrrpheno <- hvdrr$pheno
hvdrr$pheno[49, 1] <- NA
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #73 inds and 94% genotyped
png('summary_Width40.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width40.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #no significant peaks
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#    There were no LOD peaks above the threshold.

#######HvDRR41#######
hvdrr <- datWidth$HvDRR41 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrrpheno <- hvdrr$pheno
hvdrr$pheno[68, 1] <- NA
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #89 inds and 93.6% genotyped
png('summary_Width41.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width41.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr5
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-299658   5 80.3 3.8 0.00875

#JHI-Hv50k-2016-295329   5 77.2 3.34 0.014
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(5), c(77.2), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-195505   3 125 3.38 0.0262
hvdrrqtl <- makeqtl(hvdrrprob, c(3, 5), c(125, 80.3), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 3@125.3   3 125.34     2
#Q2  5@80.0   5  80.00     2

#Q1 5@77.2   5 77.178     2
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1)) #adds if there are QTL effects
pt(abs(-4.121), 84, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width41.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr3
find.marker(hvdrrprob, 3, 125.34) #JHI-Hv50k-2016-195505
png('effetplot_ch3_Width41.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-195505")
dev.off()
png('effetplotsca_ch3_Width41.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-195505") #Scatter plot
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 77.178) #SCRI_RS_203128 #BOPA1_6035-654
png('effetplot_ch5_Width41.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'SCRI_RS_203128')
dev.off()
png('effetplotsca_ch5_Width41.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'SCRI_RS_203128') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-190756   3 105.8054 1.962982
#JHI-Hv50k-2016-195505   3 125.3380 3.482407
#SCRI_RS_229593          3 156.0883 1.669677
datposition[(which(datposition == 'JHI-Hv50k-2016-190756')),9] #501186565
datposition[(which(datposition == 'JHI-Hv50k-2016-195505')),9] #529586832
datposition[(which(datposition == 'SCRI_RS_229593')),9] #550686097

hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#SCRI_RS_68142           5 69.45208 2.316641
#c5.loc80                5 80.00000 4.167717
#JHI-Hv50k-2016-307382   5 94.50152 2.275345

#JHI-Hv50k-2016-284256   5 52.80957 1.627308
#JHI-Hv50k-2016-295329   5 77.17806 3.341179
#JHI-Hv50k-2016-299042   5 79.80222 1.714297
datposition[(which(datposition == 'JHI-Hv50k-2016-284256')),9] #64968257 #17429653
datposition[(which(datposition == 'BOPA1_6035-654')),9] #322624458 #201636247
datposition[(which(datposition == 'JHI-Hv50k-2016-299042')),9] #435714162 #324725679

#######HvDRR42#######
hvdrr <- datWidth$HvDRR42 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #71 inds and 94.3% genotyped
png('summary_Width42.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width42.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr5 and 7
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#c7.loc293   7 293 4.41 0.00625
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(7), c(293), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#SCRI_RS_210928   5 456 3.68 0.02700
#c6.loc295        6 295 5.24 0.00075
hvdrrqtl <- makeqtl(hvdrrprob, c(5, 6, 7), c(456, 295, 293), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 5@456.3   5 456.34     2
#Q2 6@296.0   6 296.00     2
#Q3 7@294.0   7 294.00     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2+Q3)) #adds if there are QTL effects
pt(abs(4.610), 64, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width42.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 456.34) #SCRI_RS_210928
png('effetplot_ch5_Width42.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "SCRI_RS_210928")
dev.off()
png('effetplotsca_ch5_Width42.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "SCRI_RS_210928") #Scatter plot
dev.off()

#effect plot marker chr6
find.marker(hvdrrprob, 6, 296) #BOPA2_12_31053
png('effetplot_ch6_Width42.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'BOPA2_12_31053')
dev.off()
png('effetplotsca_ch6_Width42.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'BOPA2_12_31053') #Scatter plot
dev.off()

#effect plot marker chr7
find.marker(hvdrrprob, 7, 294) #BOPA1_11619-618
png('effetplot_ch7_Width42.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'BOPA1_11619-618')
dev.off()
png('effetplotsca_ch7_Width42.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'BOPA1_11619-618') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-337292   5 445.9212 2.371806
#SCRI_RS_210928          5 456.3395 4.247915
#JHI-Hv50k-2016-339153   5 478.2443 2.583183
datposition[(which(datposition == 'JHI-Hv50k-2016-337292')),9] #539884781
datposition[(which(datposition == 'SCRI_RS_210928')),9] #545172151
datposition[(which(datposition == 'JHI-Hv50k-2016-339153')),9] #545540141

hvdrrCI <- lodint(hvdrr_rqtl, chr = 6, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-424666   6 285.3680 3.401457
#c6.loc296               6 296.0000 5.709064
#JHI-Hv50k-2016-428413   6 298.8638 3.560152
datposition[(which(datposition == 'JHI-Hv50k-2016-424666')),9] #555104524
datposition[(which(datposition == 'BOPA2_12_31053')),9] #562145907
datposition[(which(datposition == 'JHI-Hv50k-2016-428413')),9] #563721815

hvdrrCI <- lodint(hvdrr_rqtl, chr = 7, qtl.index = 3, expandtomarkers = TRUE)
hvdrrCI
#BOPA2_12_31199   7 286.1861 5.968735
#c7.loc294        7 294.0000 8.296331
#SCRI_RS_134640   7 300.3667 6.307146
datposition[(which(datposition == 'BOPA2_12_31199')),9] #507475301
datposition[(which(datposition == 'BOPA1_11619-618')),9] #520120674
datposition[(which(datposition == 'SCRI_RS_134640')),9] #555814735
hvdrr$geno$`7`

#######HvDRR43#######
hvdrr <- datWidth$HvDRR43 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrrpheno <- hvdrr$pheno
hvdrr$pheno[c(12, 32), 1] <- NA
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #99 inds and 91.5% genotyped
png('summary_Width43.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width43.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr4
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-268870   4 131 3.38 0.0192
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(4), c(131), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#BOPA2_12_31007   6 80.8 3.22 0.028
hvdrrqtl <- makeqtl(hvdrrprob, c(4, 6), c(131, 80.8), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 4@131.1   4 131.07     2

#Q1 4@131.1   4 131.066     2
#Q2  6@80.8   6  80.793     2
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2)) #adds if there are QTL effects
pt(abs(4.063), 96, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width43.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr4
find.marker(hvdrrprob, 6, 80.793) #JHI-Hv50k-2016-268870
png('effetplot_ch4_Width43.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-268870")
dev.off()
png('effetplotsca_ch4_Width43.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-268870") #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 6, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-267991    4 125.7294 1.475787
#JHI-Hv50k-2016-268870    4 131.0658 3.378304
#JHI-Hv50k-2016-268870    4 131.0658 3.378304

#JHI-Hv50k-2016-268606    4 128.0588 3.003648
#JHI-Hv50k-2016-268870    4 131.0658 4.821718
#JHI-Hv50k-2016-268870    4 131.0658 4.821718
datposition[(which(datposition == 'JHI-Hv50k-2016-268606')),9] #606001172 #607247320
datposition[(which(datposition == 'JHI-Hv50k-2016-268870')),9] #609284484
datposition[(which(datposition == 'JHI-Hv50k-2016-268870')),9] #609284484
#SCRI_RS_189254   6  75.10525 1.619047
#BOPA2_12_31007   6  80.79291 3.228206
#BOPA1_3038-541   6 115.71902 1.101551
datposition[(which(datposition == 'SCRI_RS_189254')),9] #68265315
datposition[(which(datposition == 'BOPA2_12_31007')),9] #129652885
datposition[(which(datposition == 'BOPA1_3038-541')),9] #525277321
hvdrr$geno$`6`

#######HvDRR44#######
hvdrr <- datWidth$HvDRR44 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #87 inds and 92.9% genotyped
png('summary_Width44.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width44.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr6 and 7
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-429466   6 201 4.26 0.00375
#c7.loc182               7 182 7.06 0.00000
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(6, 7), c(201, 182), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 6@190.0   6 190     2
#Q2 7@182.0   7 182     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2)) #adds if there are QTL effects
pt(abs(5.684), 82, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width44.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr6
find.marker(hvdrrprob, 6, 190) #JHI-Hv50k-2016-427116
png('effetplot_ch6_Width44.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-427116")
dev.off()
png('effetplotsca_ch6_Width44.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-427116") #Scatter plot
dev.off()

#effect plot marker chr7
find.marker(hvdrrprob, 7, 182) #BOPA2_12_20685
png('effetplot_ch7_Width44.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'BOPA2_12_20685')
dev.off()
png('effetplotsca_ch7_Width44.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'BOPA2_12_20685') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 6, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#SCRI_RS_149269   6 183.7816 3.632086
#c6.loc190        6 190.0000 6.147960
#SCRI_RS_143994   6 203.8472 4.098914
datposition[(which(datposition == 'SCRI_RS_149269')),9] #556896720
datposition[(which(datposition == 'JHI-Hv50k-2016-427116')),9] #560750847
datposition[(which(datposition == 'SCRI_RS_143994')),9] #569415594

hvdrrCI <- lodint(hvdrr_rqtl, chr = 7, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#BOPA1_3140-491          7 177.9211 7.672288
#c7.loc182               7 182.0000 9.317727
#JHI-Hv50k-2016-492576   7 204.5096 3.666432
datposition[(which(datposition == 'BOPA1_3140-491')),9] #522384232
datposition[(which(datposition == 'BOPA2_12_20685')),9] #
datposition[(which(datposition == 'JHI-Hv50k-2016-492576')),9] #569510428

hvdrrCI <- lodint(hvdrr_rqtl, chr = , qtl.index = 3, expandtomarkers = TRUE)
hvdrrCI
#
datposition[(which(datposition == '')),9] #
datposition[(which(datposition == '')),9] #
datposition[(which(datposition == '')),9] #
#######HvDRR45#######
hvdrr <- datWidth$HvDRR45 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrrpheno <- hvdrr$pheno
hvdrr$pheno[c(7, 11), 1] <- NA
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #70 inds and 94.1% genotyped
png('summary_Width45.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width45.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr3
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#c2.loc97    2  97 2.97 0.0415
#c3.loc167   3 167 3.00 0.0382
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 3), c(97, 167), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-336856   5 199 3.15 0.0253
hvdrrqtl <- makeqtl(hvdrrprob, c(2, 3, 5), c(97, 167,199), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1  2@91.2   2  91.248     2
#Q2 3@165.3   3 165.315     2
#Q3 5@198.0   5 198.000     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2+Q3)) #adds if there are QTL effects
pt(abs(-3.969), 93, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width47.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 198.000) #JHI-Hv50k-2016-334878
png('effetplot_ch_Width47.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-334878")
dev.off()


#effect plot marker chr3
find.marker(hvdrrprob, 3, 165.315) #JHI-Hv50k-2016-101074
png('effetplot_ch6_Width47.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-101074')
dev.off()


#effect plot marker chr2
find.marker(hvdrrprob, 2, 91.248) #JHI-Hv50k-2016-168303
png('effetplot_ch6_Width47.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'JHI-Hv50k-2016-168303')
dev.off()


#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-73712    2  37.44744 0.1109606
#JHI-Hv50k-2016-101074   2  91.24834 1.7758674
#SCRI_RS_138848          2 256.96051 0.6896479
datposition[(which(datposition == 'JHI-Hv50k-2016-73712')),9] #24210792
datposition[(which(datposition == 'JHI-Hv50k-2016-101074')),9] #547071743
datposition[(which(datposition == 'SCRI_RS_138848')),9] #674152163

hvdrrCI <- lodint(hvdrr_rqtl, chr = 3, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-194821   3 127.5014 1.720594
#JHI-Hv50k-2016-202418   3 165.3155 3.428552
#JHI-Hv50k-2016-203630   3 176.2443 1.581754
datposition[(which(datposition == 'JHI-Hv50k-2016-194821')),9] #521184244
datposition[(which(datposition == 'JHI-Hv50k-2016-202418')),9] #561475004
datposition[(which(datposition == 'JHI-Hv50k-2016-203630')),9] #565272859

hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 3, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-324435   5 160.0968 1.721111
#c5.loc198               5 198.0000 3.259507
#JHI-Hv50k-2016-337874   5 204.1690 1.237111
datposition[(which(datposition == 'JHI-Hv50k-2016-324435')),9] #516935567
datposition[(which(datposition == 'JHI-Hv50k-2016-334878')),9] #535950569
datposition[(which(datposition == 'JHI-Hv50k-2016-337874')),9] #541448361
#######HvDRR46#######
hvdrr <- datWidth$HvDRR46 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #31 inds and 93.7% genotyped
png('summary_Width46.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width46.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr7
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
# There were no LOD peaks above the threshold

#######HvDRR47#######
hvdrr <- datWidth$HvDRR47 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrrpheno <- hvdrr$pheno
hvdrr$pheno[67, 1] <- NA
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #98 inds and 93% genotyped
png('summary_Width47.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width47.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr5 and 6
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#c5.loc256               5 256 5.31 0.00025
#JHI-Hv50k-2016-432894   6 251 3.84 0.01275

#c5.loc259               5 259 5.26 0.00075
#JHI-Hv50k-2016-432894   6 251 3.93 0.00700
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(5, 6), c(259, 251), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-330440   5 257 3.49 0.0203
hvdrrqtl <- makeqtl(hvdrrprob, c(5, 5, 6), c(257, 259, 251), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 5@261.4   5 261.37     2
#Q2 6@110.7   6 110.69     2

#Q1 5@257.0   5 256.98     2
#Q2 5@258.3   5 258.25     2
#Q3 6@110.7   6 110.69     2
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1+Q2+Q3)) #adds if there are QTL effects
pt(abs(-3.969), 93, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width47.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr5
find.marker(hvdrrprob, 5, 258.25) #BOPA2_12_30590
png('effetplot_ch_Width47.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-330440")
dev.off()
png('effetplotsca_ch2_Width47.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-330440") #Scatter plot
dev.off()
hvdrr$geno$`5`

#effect plot marker chr6
find.marker(hvdrrprob, 6, 110.69) #SCRI_RS_239642
png('effetplot_ch6_Width47.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = 'SCRI_RS_239642')
dev.off()
png('effetplotsca_ch6_Width47.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, 'SCRI_RS_239642') #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 5, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#SCRI_RS_126419          5 254.2333 0.2821166
#JHI-Hv50k-2016-330440   5 256.9767 4.8966784
#JHI-Hv50k-2016-330609   5 257.8306 0.1442926
datposition[(which(datposition == 'JHI-Hv50k-2016-329411')),9] #526186606
datposition[(which(datposition == 'JHI-Hv50k-2016-330440')),9] #527600633
datposition[(which(datposition == 'JHI-Hv50k-2016-330609')),9] #527785291

#JHI-Hv50k-2016-330447   5 257.4006 0.1703462
#JHI-Hv50k-2016-330490   5 258.2516 6.7669639
#SCRI_RS_90041           5 260.9464 2.6323930
datposition[(which(datposition == 'JHI-Hv50k-2016-330447')),9] #527607276
datposition[(which(datposition == 'JHI-Hv50k-2016-330490')),9] #527666352
datposition[(which(datposition == 'SCRI_RS_90041')),9] #528994103

hvdrrCI <- lodint(hvdrr_rqtl, chr = 6, qtl.index = 2, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-380943   6  56.99803 1.612240
#SCRI_RS_138529          6 110.69215 3.264579
#JHI-Hv50k-2016-433886   6 252.36583 2.239515
datposition[(which(datposition == 'JHI-Hv50k-2016-380943')),9] #29083597
datposition[(which(datposition == 'SCRI_RS_239642')),9] #356194651
datposition[(which(datposition == 'JHI-Hv50k-2016-433886')),9] #572623101

#######HvDRR48#######
hvdrr <- datWidth$HvDRR48 #extracting the single population
summary(hvdrr$pheno)
plot(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #remove the lines with missing phenotype
#set heterozygote to missing
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 2, NA, x)
})

#replace 3 with 2
hvdrr$geno$'1'$data <- apply(hvdrr$geno$'1'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'2'$data <- apply(hvdrr$geno$'2'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'3'$data <- apply(hvdrr$geno$'3'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'4'$data <- apply(hvdrr$geno$'4'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'5'$data <- apply(hvdrr$geno$'5'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'6'$data <- apply(hvdrr$geno$'6'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
hvdrr$geno$'7'$data <- apply(hvdrr$geno$'7'$data, 2, function(x){
  ifelse(x == 3, 2, x)
})
#indicate that the lines are RILS
class(hvdrr)[1] <- 'riself'
summary(hvdrr) #82 inds and 93.3% genotyped
png('summary_Width48.png', width = 6, height = 6, units = 'in', res = 300)
plot(hvdrr)
dev.off()
#Initial qtl mapping
hvdrrprob <- calc.genoprob(hvdrr, step = 1, err = 0.01, map.function = "haldane") #calculate conditional qtl probability
hvdrrout <- scanone(hvdrrprob, method = 'hk') #genomewide scanning of qtl
png('lod_Width48.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrrout, ylab = 'LOD score') #peak on chr2
dev.off()
#permutation test
set.seed(523938) #seed to make the next step repeatable
hvdrrperm <- scanone(hvdrrprob, method = "hk", n.perm = 4000)
summary(hvdrrout, perms =hvdrrperm, alpha = 0.05, pvalues = TRUE)
#JHI-Hv50k-2016-107258   2 168 12.1    0
#check for additional qtl
hvdrrqtl <- makeqtl(hvdrrprob, c(2), c(168), what = "prob")
hvdrrout.c <- addqtl(hvdrrprob,qtl =hvdrrqtl, method = "hk")
summary(hvdrrout.c, perms = hvdrrperm, alpha = 0.05, pvalues = TRUE)
#no lod peak above the threshold
#refine the estimated location of QTL
hvdrr_rqtl <- refineqtl(hvdrrprob, qtl = hvdrrqtl, method = "hk", verbose = FALSE)
hvdrr_rqtl
#Q1 2@168.3   2 168.32     2

summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method="hk"), pvalues=FALSE) #summarises all qtls in a list and shows LODs and how much they add to the variance
summary(fitqtl(hvdrrprob, qtl = hvdrr_rqtl, method = "hk", get.ests=TRUE, formula=Width~ Q1)) #adds if there are QTL effects
pt(abs(8.843), 79, lower.tail = FALSE) # t-test value degree of freedom #p < 0.01
png('qtlpos_Width48.png', width = 6, height = 4, units = 'in', res = 300)
plot(hvdrr_rqtl)
dev.off()

#effect plot marker chr2
find.marker(hvdrrprob, 2, 168.32) #JHI-Hv50k-2016-107258
png('effetplot_ch2_Width48.png', width = 4, height = 4, units = 'in', res = 300)
effectplot(hvdrrprob, mname1 = "JHI-Hv50k-2016-107258")
dev.off()
png('effetplotsca_ch2_Width48.png', width = 4, height = 4, units = 'in', res = 300)
plotPXG(hvdrrprob, "JHI-Hv50k-2016-107258") #Scatter plot
dev.off()

#Confidence interval
hvdrrCI <- lodint(hvdrr_rqtl, chr = 2, qtl.index = 1, expandtomarkers = TRUE)
hvdrrCI
#JHI-Hv50k-2016-106749   2 166.9115  9.676210
#JHI-Hv50k-2016-107258   2 168.3227 12.141332
#JHI-Hv50k-2016-109641   2 175.3968  8.972447
datposition[(which(datposition == 'JHI-Hv50k-2016-106749')),9] #578660476
datposition[(which(datposition == 'JHI-Hv50k-2016-107258')),9] #579781504
datposition[(which(datposition == 'JHI-Hv50k-2016-109641')),9] #597915659

