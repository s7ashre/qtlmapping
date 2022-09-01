library(qtl)
setwd('~/mounts/project/marvin/seedsizeQTL/experiments/qtlmapping/Width')
dat <- readRDS('Widthcross_new.RDS')
df1 <- read.csv('/1data/Asis/epiqtl/GW/epiqtl.gw.lod3.csv')
df1.ls <- levels(df1$populaiton)
#i <- 1

for (i in 1:length(df1.ls)){
  hvdrr <- dat[[df1.ls[i]]]
  print(summary(hvdrr))
  hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width')))
  
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
  
  class(hvdrr)[1] <- 'riself'
  
  hvdrrprob2 <- calc.genoprob(hvdrr, step = 2.5, err = 0.01, map.function = "haldane")
  hvdrrout2 <- scantwo(hvdrrprob2, verbose=FALSE, clean.output = T)
  save(hvdrrout2, file = paste('gw.out2.',df1.ls[i],'.RData', sep = ''))
}

#####HvDRR03####
load("operm.HvDRR03.r1_2.5.RData")
load("operm.HvDRR03.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR03.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
# (1000 permutations)
#   full  fv1  int  add  av1  one
#5% 6.62 5.06 4.15 5.31 2.88 3.01
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 4.15, ]
#####HvDRR04####
load("operm.HvDRR04.r1_2.5.RData")
load("operm.HvDRR04.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR04.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
# (1000 permutations)
#   full  fv1  int  add  av1  one
#5% 6.87 5.05 4.13 5.47 2.87 3.12
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 4.13, ]
#####HvDRR07####
load("operm.HvDRR07.r1_2.5.RData")
load("operm.HvDRR07.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR07.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
# (1000 permutations)
#   full fv1  int  add  av1  one
#5% 6.52 4.7 3.98 5.15 2.86 3.09
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 3.98, ]
#####HvDRR09####
load("operm.HvDRR09.r1_2.5.RData")
load("operm.HvDRR09.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR09.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
# (1000 permutations)
#   full  fv1  int  add  av1  one
#5% 6.22 4.57 3.84 5.07 2.86 2.92
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 3.84, ]
#####HvDRR12####
load("operm.HvDRR12.r1_2.5.RData")
load("operm.HvDRR12.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR12.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
# (1000 permutations)
#   full  fv1  int  add  av1  one
#5% 6.56 5.01 3.88 5.24 3.56 3.01
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 3.88, ]
#chr1 chr2 pos1f pos2f lod.full pval  lod.fv1 pval  lod.int  pval pos1a pos2a  lod.add  pval lod.av1  pval
#   2    5  97.5  52.5 10.09207    0 7.246651    0 4.419132 0.017   205  57.5 5.672939 0.018 2.82752 0.243
hvdrr <- dat$HvDRR12
summary(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width')))
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
class(hvdrr)[1] <- 'riself'
mar <- find.marker(hvdrr, c('2','5'), c(97.5, 52.5))
mar
png('d12.gw.h2_5.epiqtl.png', height = 5, width = 10, units = 'in', res = 300)
par(mfrow=c(1,2))
plotPXG(hvdrr, marker=mar)
effectplot(hvdrr, mname1="JHI-Hv50k-2016-90688", mname2="JHI-Hv50k-2016-289452", ylim=range(pull.pheno(hvdrr,1)), main = NULL)
dev.off()



#####HvDRR15####
load("operm.HvDRR15.r1_2.5.RData")
load("operm.HvDRR15.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR15.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
# (1000 permutations)
#   full  fv1  int  add  av1  one
#5%  6.1 4.63 3.91 5.07 2.95 2.88
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 3.91, ]

#####HvDRR20####
load("operm.HvDRR20.r1_2.5.RData")
load("operm.HvDRR20.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR20.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
# (1000 permutations)
#   full  fv1  int  add  av1  one
#5% 6.47 4.81 3.83 5.21 2.89 3.09
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 3.83, ]


#####HvDRR25####
load("operm.HvDRR25.r1_2.5.RData")
load("operm.HvDRR25.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR25.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
# (1000 permutations)
#   full fv1  int  add av1  one
#5% 17.3  16 8.61 16.3  15 2.03
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 8.61, ]
#####HvDRR27####
load("operm.HvDRR27.r1_2.5.RData")
load("operm.HvDRR27.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR27.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
# (1000 permutations)
#   full  fv1  int  add  av1  one
#5% 6.76 4.89 3.98 5.43 3.23 3.14
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 3.98, ]
#####HvDRR30####
load("operm.HvDRR30.r1_2.5.RData")
load("operm.HvDRR30.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR30.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
# (1000 permutations)
#   full  fv1  int add  av1  one
#5%  6.9 5.13 4.27 5.3 2.72 3.22
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 4.27, ]
#####HvDRR32####
load("operm.HvDRR32.r1_2.5.RData")
load("operm.HvDRR32.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR32.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
# (1000 permutations)
#   full  fv1  int  add  av1  one
#5% 7.18 5.32 4.26 5.56 3.46 3.08
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 4.26, ]
#####HvDRR34####
load("operm.HvDRR34.r1_2.5.RData")
load("operm.HvDRR34.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR34.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
# (1000 permutations)
#   full  fv1  int  add  av1  one
#5% 6.67 5.01 4.01 5.74 3.23 3.11
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 4.01, ]





#####HvDRR37####
load("operm.HvDRR37.r1_2.5.RData")
load("operm.HvDRR37.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR37.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
#   full  fv1  int  add  av1  one
#5% 6.28 4.65 3.61 5.53 3.53 2.99
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 3.61, ]
#####HvDRR41####
load("operm.HvDRR41.r1_2.5.RData")
load("operm.HvDRR41.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR41.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
#   full  fv1  int  add  av1  one
#5% 9.86 8.67 3.96 8.52 7.21 2.66
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 3.61, ]
# chr1 chr2 pos1f pos2f lod.full  pval  lod.fv1  pval  lod.int  pval pos1a pos2a  lod.add pval  lod.av1 pval
#   6    6   135   140 9.076808 0.225 7.399745 0.629 4.146768 0.032    30  32.5 4.930041    1 3.252977    1
hvdrr <- dat$HvDRR41
summary(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width')))
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
class(hvdrr)[1] <- 'riself'
mar <- find.marker(hvdrr, c('6','6'), c(135, 140))
mar
png('d41.gw.h6_6.epiqtl.png', height = 5, width = 10, units = 'in', res = 300)
par(mfrow=c(1,2))
plotPXG(hvdrr, marker=mar)
effectplot(hvdrr, mname1="JHI-Hv50k-2016-419702", mname2="BOPA2_12_31115", ylim=range(pull.pheno(hvdrr,1)), main = NULL)
dev.off()
#####HvDRR43####
load("operm.HvDRR43.r1_2.5.RData")
load("operm.HvDRR43.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR43.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
#   full  fv1  int add  av1  one
#5% 10.5 9.58 6.67 7.2 5.97 2.76
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 6.67, ]
#####HvDRR44####
load("operm.HvDRR44.r1_2.5.RData")
load("operm.HvDRR44.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR44.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
#   full  fv1 int  add  av1  one
#5% 6.49 4.85   4 5.12 2.89 3.05
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 4, ]
#####HvDRR45####
load("operm.HvDRR45.r1_2.5.RData")
load("operm.HvDRR45.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR45.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
#   full  fv1  int  add  av1  one
#5% 8.64 7.45 4.49 7.36 6.27 2.845
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 4.49, ]
#   chr1 chr2 pos1f pos2f lod.full  pval  lod.fv1  pval  lod.int  pval pos1a pos2a  lod.add pval  lod.av1 pval
#    4    5   180   195 7.805074 0.244 5.680194 0.709 4.526035 0.048  82.5   195 3.279039    1 1.154159    1
hvdrr <- dat$HvDRR45
summary(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width')))
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
class(hvdrr)[1] <- 'riself'
mar <- find.marker(hvdrr, c('2','5'), c(180, 195))
mar
png('d45.gw.h4_5.epiqtl.png', height = 5, width = 10, units = 'in', res = 300)
par(mfrow=c(1,2))
plotPXG(hvdrr, marker=mar)
effectplot(hvdrr, mname1="JHI-Hv50k-2016-117094", mname2="JHI-Hv50k-2016-333733", ylim=range(pull.pheno(hvdrr,1)), main = NULL)
dev.off()



#####HvDRR46####
load("operm.HvDRR46.r1_2.5.RData")
load("operm.HvDRR46.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR46.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
#   full  fv1  int  add av1  one
#5% 7.69 5.97 5.01 6.14 3.2 3.65
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 5.01, ]


#####HvDRR48####
load("operm.HvDRR48.r1_2.5.RData")
load("operm.HvDRR48.r2_2.5.RData")
out.perm1k <- c(out.perm2a, out.perm2b)
load("gw.out2.HvDRR48.RData")
summary(out.perm1k, 0.05)
#Width (1000 permutations)
#   full  fv1 int  add  av1  one
#5% 6.55 4.94 4.2 5.44 2.73 3.18
df.hvdrro2 <- as.data.frame(summary(hvdrrout2, perms=out.perm1k, alphas=c(1, 1, 0, 1, 1), pvalues=TRUE))
df.hvdrro2[df.hvdrro2$lod.int >= 4.2, ]
#   chr1 chr2 pos1f pos2f lod.full pval  lod.fv1  pval  lod.int  pval pos1a pos2a  lod.add pval  lod.av1 pval
#    2    4   170    65 16.82425    0 5.416803 0.016 4.383658 0.029   170    50 12.44059    0 1.033145    1
hvdrr <- dat$HvDRR48
summary(hvdrr)
hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width')))
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
class(hvdrr)[1] <- 'riself'
  
mar <- find.marker(hvdrr, c('2','4'), c(170, 65))
mar
png('d48.gw.h2_4.epiqtl.png', height = 5, width = 10, units = 'in', res = 300)
par(mfrow=c(1,2))
plotPXG(hvdrr, marker=mar)
effectplot(hvdrr, mname1="JHI-Hv50k-2016-108147", mname2="JHI-Hv50k-2016-231722", ylim=range(pull.pheno(hvdrr,1)), main = NULL)
dev.off()
