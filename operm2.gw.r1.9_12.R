setwd('/gpfs/project/shrestha/marvin/seedsizeQTL/experiments/qtlmapping')
#setwd('/home/shrestha/mounts/project/shrestha/marvin/seedsizeQTL/experiments/qtlmapping/')
library(qtl)
dat <- readRDS('Width/Widthcross_new.RDS') #change
df1 <- read.csv('Width/epiqtl.gw.lod3.csv') #change
df1.ls <- levels(df1$populaiton)
# 20 levels
#i <- 1

for (i in 9:12){
  hvdrr <- dat[[df1.ls[i]]]
  print(summary(hvdrr))
  hvdrr <- subset(hvdrr, ind=!is.na(pull.pheno(hvdrr,'Width'))) #change
  
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
  set.seed(85843000) #change
  out.perm2a <- scantwo(hvdrrprob2, n.perm = 500, clean.output=TRUE)
  save(out.perm2a, file = paste('Width/','','operm.',df1.ls[i],'.r1_2.5.RData', sep = '')) #change
}

