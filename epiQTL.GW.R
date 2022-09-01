setwd('/gpfs/project/shrestha/marvin/seedsizeQTL/experiments/qtlmapping')
#setwd('~/mounts/project/shrestha/marvin/seedsizeQTL/experiments/qtlmapping')
#install.packages("qtl")
library(qtl)
dat <- readRDS('Width/Widthcross_new.RDS')
#i <- 2
df1 <- NULL
for ( i in 1:length(dat)){
  hvdrr <- dat[[i]]
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
  summary(hvdrr)
  hvdrrprob2 <- calc.genoprob(hvdrr, step = 2.5, err = 0.01, map.function = "haldane")
  hvdrrout2 <- scantwo(hvdrrprob2, verbose=FALSE)
  summary(hvdrrout2)
  df.s2sum <- data.frame(summary(hvdrrout2))
  df.s2sum$populaiton <- names(dat[i])
  df1 <- rbind(df1, df.s2sum)
} 
write.csv(df1, 'Width/epiQTL.gw2.5.csv', row.names = F)