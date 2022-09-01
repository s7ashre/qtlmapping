setwd('/home/shrestha/Documents/dailywork')
dat <- read.csv('koln2019_final.csv')
test_dup <- dat[(duplicated(dat[,1:2]) | duplicated(dat[,1:2])), ]
dat <- dat[!(duplicated(dat[,1:2]) | duplicated(dat[,1:2])), ]
fieldplan2019 <- read.csv('Field_Plan_Cologne33K_2019_CORRECTED.csv')
fieldplan2019 <- as.matrix(fieldplan2019)
fieldplan2019 <- fieldplan2019[,-c(1:8, 84:ncol(fieldplan2019))]
fieldplan2019 <- fieldplan2019[-c(1:2, 65:nrow(fieldplan2019)),]
rownames(fieldplan2019) <- nrow(fieldplan2019):1
colnames(fieldplan2019) <- 1:ncol(fieldplan2019)
#melting the matrix
library(reshape2)
fieldplan2019_melt <- melt(fieldplan2019)
colnames(fieldplan2019_melt) <- c("Range", "Row", "Lotcode")
write.csv(fieldplan2019_melt, 'fieldmap_koln2019.csv', row.names = F)
pheno_mvg <- merge(fieldplan2019_melt, dat, by = c('Range', 'Row'), all.x = TRUE)

####crossvalidation if the merging was correct#####
test_pheno_mvg <- na.omit(pheno_mvg)
str(test_pheno_mvg)
test_pheno_mvg$Lotcode.x <- as.character(test_pheno_mvg$Lotcode.x)
test_pheno_mvg$Lotcode.y <- as.character(test_pheno_mvg$Lotcode.y)
#identical(test_pheno_mvg[['Lotcode.x']],test_pheno_mvg[['Lotcode.y']])
identical(test_pheno_mvg$Lotcode.x,test_pheno_mvg$Lotcode.y)
all(test_pheno_mvg$Lotcode.x == test_pheno_mvg$Lotcode.y)


koln_2019_rows <- rep(1:62, each = 75)
koln_2019_columns <- rep(1:75, 62)
#######Area######
Unlisted_pheno_plan <- as.vector(pheno_mvg$Area)
datArea <- dat[, -c(3:4, 6:7)]
###1:3
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:3, 1:3, 1:3, 1:3), layers = 1:3) 
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjArea_3 <- merge(datArea, info, by = c('Range', 'Row'), all.x =TRUE)
adjArea_3 <- adjArea_3[,-c(12:13)]
colnames(adjArea_3)[10:11] <- c('adjArea_3', 'obsArea')
###1:5
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:5, 1:5, 1:5, 1:5), layers = 1:5)
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjArea_5 <- merge(datArea, info, by = c('Range', 'Row'), all.x =TRUE)
adjArea_5 <- adjArea_5[,-c(4:9, 11:13)]
colnames(adjArea_5)[3:4] <- c('Areacheck_5','adjArea_5')
adjArea <- merge(adjArea_3, adjArea_5, by = c('Range', 'Row'), all.x =TRUE )
####1:10
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:10, 1:10, 1:10, 1:10), layers = 1:10) 
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjArea_10 <- merge(datArea, info, by = c('Range', 'Row'), all.x =TRUE)
adjArea_10 <- adjArea_10[,-c(4:9, 11:13)]
colnames(adjArea_10)[3:4] <- c('Areacheck_10','adjArea_10')
adjArea <- merge(adjArea, adjArea_10, by = c('Range', 'Row'), all.x =TRUE )
###1:15
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:15, 1:15, 1:15, 1:15), layers = 1:15) 
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjArea_15 <- merge(datArea, info, by = c('Range', 'Row'), all.x =TRUE)
adjArea_15 <- adjArea_15[,-c(4:9, 11:13)]
colnames(adjArea_15)[3:4] <- c('Areacheck_15','adjArea_15')
adjArea <- merge(adjArea, adjArea_15, by = c('Range', 'Row'), all.x =TRUE )
write.csv(adjArea, 'adjArea_movgrad_koln2019.csv', row.names = FALSE)

#######Width######
Unlisted_pheno_plan <- as.vector(pheno_mvg$Width)
datWidth <- dat[, -c(3:5, 7)]
###1:3
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:3, 1:3, 1:3, 1:3), layers = 1:3) 
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjWidth_3 <- merge(datWidth, info, by = c('Range', 'Row'), all.x =TRUE)
adjWidth_3 <- adjWidth_3[,-c(12:13)]
colnames(adjWidth_3)[10:11] <- c('adjWidth_3', 'obsWidth')
###1:5
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:5, 1:5, 1:5, 1:5), layers = 1:5)
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjWidth_5 <- merge(datWidth, info, by = c('Range', 'Row'), all.x =TRUE)
adjWidth_5 <- adjWidth_5[,-c(4:9, 11:13)]
colnames(adjWidth_5)[3:4] <- c('Widthcheck_5','adjWidth_5')
adjWidth <- merge(adjWidth_3, adjWidth_5, by = c('Range', 'Row'), all.x =TRUE )
####1:10
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:10, 1:10, 1:10, 1:10), layers = 1:10) 
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjWidth_10 <- merge(datWidth, info, by = c('Range', 'Row'), all.x =TRUE)
adjWidth_10 <- adjWidth_10[,-c(4:9, 11:13)]
colnames(adjWidth_10)[3:4] <- c('Widthcheck_10','adjWidth_10')
adjWidth <- merge(adjWidth, adjWidth_10, by = c('Range', 'Row'), all.x =TRUE )
###1:15
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:15, 1:15, 1:15, 1:15), layers = 1:15) 
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjWidth_15 <- merge(datWidth, info, by = c('Range', 'Row'), all.x =TRUE)
adjWidth_15 <- adjWidth_15[,-c(4:9, 11:13)]
colnames(adjWidth_15)[3:4] <- c('Widthcheck_15','adjWidth_15')
adjWidth <- merge(adjWidth, adjWidth_15, by = c('Range', 'Row'), all.x =TRUE )
write.csv(adjWidth, 'adjWidth_movgrad_koln2019.csv', row.names = FALSE)

#######Length######
Unlisted_pheno_plan <- as.vector(pheno_mvg$Length)
datLength <- dat[, -c(3:6)]
###1:3
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:3, 1:3, 1:3, 1:3), layers = 1:3) 
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjLength_3 <- merge(datLength, info, by = c('Range', 'Row'), all.x =TRUE)
adjLength_3 <- adjLength_3[,-c(12:13)]
colnames(adjLength_3)[10:11] <- c('adjLength_3', 'obsLength')
###1:5
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:5, 1:5, 1:5, 1:5), layers = 1:5)
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjLength_5 <- merge(datLength, info, by = c('Range', 'Row'), all.x =TRUE)
adjLength_5 <- adjLength_5[,-c(4:9, 11:13)]
colnames(adjLength_5)[3:4] <- c('Lengthcheck_5','adjLength_5')
adjLength <- merge(adjLength_3, adjLength_5, by = c('Range', 'Row'), all.x =TRUE )
####1:10
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:10, 1:10, 1:10, 1:10), layers = 1:10) 
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjLength_10 <- merge(datLength, info, by = c('Range', 'Row'), all.x =TRUE)
adjLength_10 <- adjLength_10[,-c(4:9, 11:13)]
colnames(adjLength_10)[3:4] <- c('Lengthcheck_10','adjLength_10')
adjLength <- merge(adjLength, adjLength_10, by = c('Range', 'Row'), all.x =TRUE )
###1:15
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:15, 1:15, 1:15, 1:15), layers = 1:15) 
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjLength_15 <- merge(datLength, info, by = c('Range', 'Row'), all.x =TRUE)
adjLength_15 <- adjLength_15[,-c(4:9, 11:13)]
colnames(adjLength_15)[3:4] <- c('Lengthcheck_15','adjLength_15')
adjLength <- merge(adjLength, adjLength_15, by = c('Range', 'Row'), all.x =TRUE )
write.csv(adjLength, 'adjLength_movgrad_koln2019.csv', row.names = FALSE)

#######TGW######
Unlisted_pheno_plan <- as.vector(pheno_mvg$TGW)
datTGW <- dat[, -c(3, 5:7)]
###1:3
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:3, 1:3, 1:3, 1:3), layers = 1:3) 
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjTGW_3 <- merge(datTGW, info, by = c('Range', 'Row'), all.x =TRUE)
adjTGW_3 <- adjTGW_3[,-c(12:13)]
colnames(adjTGW_3)[10:11] <- c('adjTGW_3', 'obsTGW')
###1:5
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:5, 1:5, 1:5, 1:5), layers = 1:5)
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjTGW_5 <- merge(datTGW, info, by = c('Range', 'Row'), all.x =TRUE)
adjTGW_5 <- adjTGW_5[,-c(4:9, 11:13)]
colnames(adjTGW_5)[3:4] <- c('TGWcheck_5','adjTGW_5')
adjTGW <- merge(adjTGW_3, adjTGW_5, by = c('Range', 'Row'), all.x =TRUE )
####1:10
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:10, 1:10, 1:10, 1:10), layers = 1:10) 
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjTGW_10 <- merge(datTGW, info, by = c('Range', 'Row'), all.x =TRUE)
adjTGW_10 <- adjTGW_10[,-c(4:9, 11:13)]
colnames(adjTGW_10)[3:4] <- c('TGWcheck_10','adjTGW_10')
adjTGW <- merge(adjTGW, adjTGW_10, by = c('Range', 'Row'), all.x =TRUE )
###1:15
test <- movingGrid(rows = koln_2019_rows, columns = koln_2019_columns, obsPhe = Unlisted_pheno_plan,
                   shapeCross = list(1:15, 1:15, 1:15, 1:15), layers = 1:15) 
info <- entryData(test)
colnames(info)[1:2] <- c('Range', 'Row')
adjTGW_15 <- merge(datTGW, info, by = c('Range', 'Row'), all.x =TRUE)
adjTGW_15 <- adjTGW_15[,-c(4:9, 11:13)]
colnames(adjTGW_15)[3:4] <- c('TGWcheck_15','adjTGW_15')
adjTGW <- merge(adjTGW, adjTGW_15, by = c('Range', 'Row'), all.x =TRUE )
write.csv(adjTGW, 'adjTGW_movgrad_koln2019.csv', row.names = FALSE)

#####correlationtest####
library(PerformanceAnalytics)
datplot <- read.csv('adjArea_movgrad_koln2019.csv')
mydata <- datplot[, c(3, 10, 13, 15, 17)]
chart.Correlation(mydata, histogram = TRUE, pch = 19)

datplot <- read.csv('adjWidth_movgrad_koln2019.csv')
mydata <- datplot[, c(3, 10, 13, 15, 17)]
chart.Correlation(mydata, histogram = TRUE, pch = 19)

datplot <- read.csv('adjLength_movgrad_koln2019.csv')
mydata <- datplot[, c(3, 10, 13, 15, 17)]
chart.Correlation(mydata, histogram = TRUE, pch = 19)

datplot <- read.csv('adjTGW_movgrad_koln2019.csv')
mydata <- datplot[, c(3, 10, 13, 15, 17)]
chart.Correlation(mydata, histogram = TRUE, pch = 19)
