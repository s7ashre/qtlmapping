setwd('/home/shrestha/experiments/qtlmapping/mpp_marvin/mppAsis')
library(mppR)
library(ggplot2)
library(farver)

mppData <- readRDS("/home/shrestha/experiments/qtlmapping/mpp_marvin/mppTGW.RDS")

mppData <- QC.mppData(mppData = mppData, n.lim = 15, MAF.pop.lim = 0.05, MAF.cr.miss = TRUE, mk.miss = 0.1,
                        gen.miss = 0.25, verbose = TRUE)

mppData <- IBS.mppData(mppData = mppData)

mppData <- IBD.mppData(mppData = mppData, het.miss.par = TRUE, type = "RIL", type.mating = "selfing")

#c è la parte mancante di par clu

QTL_proc <- mpp_proc(pop.name = "HvDRR", trait.name = "TGW", trait = "TGW", mppData = mppData,
                     Q.eff = "par", plot.gen.eff = TRUE, N.cim = 1, thre.cof = 3, win.cof = 20, window = 20,
                     thre.QTL = 3, win.QTL = 20, CI = TRUE, drop = 1.5, verbose = FALSE, output.loc = tempdir())

QTL_proc

MQE <- MQE_proc(pop.name = "HvDRR", trait.name = "TGW", mppData = mppData, Q.eff = c("par","biall","cr"), window = 20, 
                verbose = FALSE, output.loc = tempdir())

MQE

SIM <- mpp_SIM(mppData = mppData, Q.eff = "par")

SIM

cofactors <- QTL_select(Qprof = SIM)

cofactors

CIM <- mpp_CIM(mppData = mppData, Q.eff = "par", cofactors = cofactors, plot.gen.eff = TRUE)

CIM

QTL <- QTL_select(Qprof = CIM)

QTL

gen.eff <- QTL_gen_effects(mppData = mppData, QTL = QTL, Q.eff = "par")

gen.eff

summary(gen.eff, QTL = 1)

plot(x = CIM, QTL = QTL, type = "l")
png('QTLplot_parental.png', width = 20, height = 12, units = 'in', res = 300)
Plot_result <- plot(x = CIM, QTL = QTL, type = "l")
dev.off()

Plot_result <- plot(x = CIM, QTL = QTL, type = "l")
save(Plot_result, file = "Parental_Plot_result.RData")

#plot(x = CIM, gen.eff= TRUE, mppData = mppData, QTL = QTL, Q.eff = "par", main = "QTL genetic effect plot")

png('plot_parenteffect.png', width = 20, height = 12, units = 'in', res = 300)
Plot_parents_effect <- plot(x = CIM, gen.eff= TRUE, mppData = mppData, QTL = QTL, Q.eff = "par", main = "QTL genetic effect plot")
dev.off()

Plot_parents_effect <- plot(x = CIM, gen.eff= TRUE, mppData = mppData, QTL = QTL, Q.eff = "par", main = "QTL genetic effect plot")
save(Plot_parents_effect, file = "Plot_parents_effect.RData")


#mpp_CV(pop.name = "HvDRR", trait.name = "TGW", mppData = mppData, Q.eff = "par", her = 0.4, Rep = 1,
       #k = 3, verbose = FALSE, output.loc = tempdir())
CV <- mpp_CV(pop.name = "HvDRR", trait.name = "TGW", mppData = mppData, Q.eff = "par", her = 0.4, Rep = 1,
             k = 3, verbose = FALSE, output.loc = tempdir('/home/shrestha/experiments/qtlmapping/mpp_marvin'))
save(CV, file = "Parental_CV.RData")

mpp_CV(pop.name = "HvDRR", trait.name = "TGW", mppData = mppData, Q.eff = "par", her = 0.4, Rep = 1,
             k = 3, verbose = FALSE, output.loc = tempdir('/home/shrestha/experiments/qtlmapping/mpp_marvin'))


