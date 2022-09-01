
### This is a modification of the aug.rowcol function of the package plant breeding
### For a detailed documentation - check the original documentation
### install.packages("plantbreeding", repos="http://R-Forge.R-project.org")
### library(plantbreeding)
### help(plantbreeding)
### it performs a correction of field heterogeneity in a augmented row column design
### compared to to the original version of the function it has no restriction any more on the number of included checks and their replication number
### B.Stich - 1.12.2018

aug.rowcol.BST <- function (dataframe, rows, columns, genotypes, yield, minimumrepcheck){

    dataframe <- data.frame(rows = dataframe[, rows], columns = dataframe[, columns], genotypes = dataframe[, genotypes], yield = dataframe[,yield])
    dataframe$rows = factor(dataframe$rows)
    dataframe$columns = factor(dataframe$columns)
    dataframe$genotypes = factor(dataframe$genotypes)
    geno <- data.frame(table(dataframe$genotypes))
    #gr = as.vector(geno$Var1[geno$Freq == length(levels(dataframe$rows))])
    gr = as.vector(geno$Var1[geno$Freq >= minimumrepcheck])
    mydf1 = dataframe[dataframe$genotypes %in% gr, ]
    # remove NA
    mydf1 <- mydf1[!is.na(mydf1[,4]),]
    model1 = lm(yield ~ genotypes + rows + columns, data = mydf1)
    amodel1 <- anova(model1)
    model1.res = resid(model1)
    dev.new(width = 8, height = 8)
    plot(mydf1$yield, model1.res, ylab = "Residuals", xlab = "Trait", 
        main = "residual plot of checks")
    text(mydf1$yield, model1.res, labels = mydf1$genotypes, pos = 3, 
        cex = 0.75, col = "red")
    abline(0, 0, col = "blue4")
    Df <- amodel1$Df
    "Sum Sq" <- amodel1$"Sum Sq"
    "Mean Sq" <- amodel1$"Mean Sq"
    "F value" <- amodel1$"F value"
    "Pr(>F)" <- amodel1$"Pr(>F)"
    amod2 <- data.frame(Df = amodel1$Df, `Sum Sq` = amodel1$"Sum Sq", 
        `Mean Sq` = amodel1$"Mean Sq", `F value` = amodel1$"F value", 
        `Pr(>F)` = amodel1$"Pr(>F)", check.names = FALSE)
    rownames(amod2) <- c("Genotypes", "rows", "columns", "Residual")
    class(amod2) <- c("anova", "data.frame")
    print(amod2)
    rowdf = data.frame(aggregate(yield ~ rows, data = mydf1, 
        mean))
    rowdf$yield <- rowdf$yield - mean(rowdf$yield)
    columndf = data.frame(aggregate(yield ~ columns, data = mydf1, 
        mean))
    columndf$yield <- columndf$yield - mean(columndf$yield)
    pr = as.vector(geno$Var1[geno$Freq < length(levels(dataframe$rows))])
    newaugdata1 <- data.frame(dataframe[dataframe$genotypes %in% 
        gr, ])
    row_match <- match(newaugdata1$rows, rowdf$rows)
    col_match <- match(newaugdata1$columns, columndf$columns)
    mydf1.adj = transform(newaugdata1, yield.adj = yield - rowdf[row_match, 
        "yield"] - columndf[col_match, "yield"])
    cat("Phenotypes and adjusted values : ", "\n\n")
    print(mydf1.adj)
    plot.new()
    plot(mydf1.adj$yield, mydf1.adj$yield.adj, xlab = "yield", 
        ylab = " yield adjusted")
    abline(a = 0, b = 1, col = "red", lty = 2)
    plot.new()
    dev.new(width = 8, height = 8)
    plot(mydf1.adj$yield, mydf1.adj$yield.adj, xlab = "yield", 
        ylab = " yield adjusted", asp = 1, pch = 16, col = "lightseagreen")
    a = 0
    b = 1
    abline(a, b, col = "red", lty = 1)
    segments(mydf1.adj$yield, mydf1.adj$yield.adj, mydf1.adj$yield, 
        mydf1.adj$yield, col = "blue1", lty = 2)
    text(mydf1.adj$yield, mydf1.adj$yield.adj, pos = 3, cex = 0.5)
    grid(NULL, NULL, lty = 6, col = "cornsilk2")
    cat("Standard error of  different comparisions", "\n\n")
    MSE = amodel1$"Mean Sq"[4]
    SEcheck = (2 * MSE/(amodel1$Df[2] + 1))^0.5
    cat("Difference between check means: ", SEcheck, "\n\n")
    SEwithin = (2 * MSE + (2 * MSE)/(amodel1$Df[2] + 1))^0.5
    cat("Difference adjusted yield of two varities in same row or column : ", 
        SEwithin, "\n\n")
    SEdiff = (2 * MSE + (4 * MSE)/(amodel1$Df[2] + 1))^0.5
    cat("Difference between two varieties in different rows or blocks: ", 
        SEdiff, "\n\n")
    SEgcheck = (MSE + ((3 * MSE)/(amodel1$Df[2] + 1)) - ((2 * 
        MSE)/((amodel1$Df[2] + 1)^2)))^0.5
    cat("Difference between two varieties and a check mean: ", 
        SEgcheck, "\n\n")
    results = list(ANOVA = amod2, Adjustment = mydf1.adj, se_check = SEcheck, 
        se_within = SEwithin, se_diff = SEdiff, se_geno_check = SEgcheck)
    return(results)
}

### Helper function of gobacktoF3
### datamatrix is the csv file of the database
### returns the homozygosity of the lotcode
### B.Stich - 1.12.2018
### mod 23.1.19
 search.homozyg <- function(lotcode,datamatrix){
   collot <- c(1:ncol(datamatrix))[is.element(colnames(datamatrix),"Lotcode")]
   colhomo <- c(1:ncol(datamatrix))[is.element(colnames(datamatrix),"Homozygosity")]
   a <- is.element(datamatrix[,collot],lotcode)
   return( datamatrix[a,colhomo])
 }

### Helper function of gobacktoF3
### datamatrix is the csv file of the database
### returns the lot and genotype code of P1
### B.Stich - 1.12.2018
### mod 23.1.19
 search.P1 <- function(lotcode,datamatrix){
   collot <- c(1:ncol(datamatrix))[is.element(colnames(datamatrix),"Lotcode")]
   colgenoP1 <- c(1:ncol(datamatrix))[is.element(colnames(datamatrix),"GenotypecodeP1")]
   collotP1 <- c(1:ncol(datamatrix))[is.element(colnames(datamatrix),"LotcodeP1")]
   a <- is.element(datamatrix[,2],lotcode)
   list(genoP1=as.character(datamatrix[a,13]),lotP1=as.character(datamatrix[a,14]))
 }

### Helper function of gobacktoF3
### datamatrix is the csv file of the database
### returns the production location of the lotcode
### B.Stich - 1.12.2018
### mod 23.1.19
 search.prodloc <- function(lotcode,datamatrix){
   collot <- c(1:ncol(datamatrix))[is.element(colnames(datamatrix),"Lotcode")]
   colprodloc <- c(1:ncol(datamatrix))[is.element(colnames(datamatrix),"Productionlocation")]
   a <- is.element(datamatrix[,collot],lotcode)
   as.character(datamatrix[a,colprodloc ])
 }

### Helper function of gobacktoF3
### datamatrix is the csv file of the database
### returns the genotype code of the lotcode
### B.Stich - 1.12.2018
### mod 23.1.19
 search.geno <- function(lotcode,datamatrix){
   collot <- c(1:ncol(datamatrix))[is.element(colnames(datamatrix),"Lotcode")]
   colgeno <- c(1:ncol(datamatrix))[is.element(colnames(datamatrix),"Genotypecode")]
   a <- is.element(datamatrix[,collot],lotcode)
   as.character(datamatrix[a,colgeno])
 }

### Helper function of gobacktoF3
### datamatrix is the csv file of the database
### returns the population name of the genotype code
### B.Stich - 1.12.2018
### mod 23.1.19
 search.pop <- function(genocode,datamatrix){
   colgeno <- c(1:ncol(datamatrix))[is.element(colnames(datamatrix),"Genotypecode")]
   colpop <- c(1:ncol(datamatrix))[is.element(colnames(datamatrix),"Population")]
   a <- is.element(datamatrix[,colgeno],genocode)
   as.character(datamatrix[a,colpop])
 }

### this function is usefull to analyse phenotypic data of different inbreeding generations of a selfing species
### datamatrix is the csv file of the database
### it returns the lotcode of the Lot that is the F3 ancestor of the searched Lotcode
### B.Stich - 1.12.2018
 gobacktoF3 <- function(Lotcode,datamatrix){
  name <- Lotcode
  for(j in 1:1000){
     prodloc <- search.prodloc(lotcode=name,datamatrix)
     if(prodloc!="ET"){
     homozyg <- search.homozyg(lotcode=name,datamatrix)
       if(homozyg>3){
         res <- search.P1(lotcode=name,datamatrix)
         name <- res$lotP1
       }
       else{
         namegeno <- search.geno(lotcode=name,datamatrix)
         namegeno
         break()
       }
     }
     else break()
   }
   namegeno <- search.geno(lotcode=name,datamatrix)
   namegeno
 }


