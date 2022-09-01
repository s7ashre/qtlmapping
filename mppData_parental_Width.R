#setwd('/home/shrestha/mounts/hpc_home/experiments/qtlmapping/mpp_marvin/mppAsis/QTL_effect_parental/Width')
setwd('/home/shrestha/experiments/qtlmapping/mpp_marvin/mppAsis/QTL_effect_parental/Width')
library(mppR)
library(ggplot2)
library(farver)
library(parallel)
my.loc <- '/home/shrestha/experiments/qtlmapping/mpp_marvin/mppAsis/QTL_effect_parental/Width'

#mppData <- readRDS("/home/shrestha/mounts/hpc_home/experiments/qtlmapping/mpp_marvin/mppWidth.RDS")
mppData <- readRDS("/home/shrestha/experiments/qtlmapping/mpp_marvin/mppWidth.RDS")
mppData <- QC.mppData(mppData = mppData, n.lim = 15, MAF.pop.lim = 0.05, MAF.cr.miss = TRUE, mk.miss = 0.1,
                        gen.miss = 0.25, verbose = TRUE, n.cores = 4)
mppData <- IBS.mppData(mppData = mppData)
mppData <- IBD.mppData(mppData = mppData, het.miss.par = TRUE, type = "RIL", type.mating = "selfing")

#c è la parte mancante di par clu

#QTL_proc <- mpp_proc(pop.name = "HvDRR", trait.name = "Width", trait = "Width", mppData = mppData,
                    # Q.eff = "par", plot.gen.eff = TRUE, N.cim = 1, thre.cof = 3, win.cof = 20, window = 20,
                     #thre.QTL = 3, win.QTL = 20, CI = TRUE, drop = 1.5, verbose = FALSE, n.cores = 4, output.loc = my.loc)
#QTL_proc
#save(QTL_proc, file = 'QTL_proc_Width.RData')

#MQE <- MQE_proc(pop.name = "HvDRR", trait.name = "Width", mppData = mppData, Q.eff = c('par'), window = 20, n.cores = 4,
                #verbose = FALSE, output.loc = my.loc)
#MQE

SIM <- mpp_SIM(mppData = mppData, Q.eff = "par", n.cores = 4)
#SIM

cofactors <- QTL_select(Qprof = SIM)
#cofactors

CIM <- mpp_CIM(mppData = mppData, Q.eff = "par", cofactors = cofactors, plot.gen.eff = TRUE, n.cores = 4)
#CIM
write.csv(CIM, 'CIM_width.csv', row.names = F)
save(CIM, file = 'CIM_width.RData')

QTL <- QTL_select(Qprof = CIM)
#QTL
write.csv(QTL, 'QTL_width.csv', row.names = F)
gen.eff <- QTL_gen_effects(mppData = mppData, QTL = QTL, Q.eff = "par")
#gen.eff
#summary(gen.eff, QTL = 1)
save(gen.eff, file = 'gen.eff_width.RData')

#plot(x = CIM, QTL = QTL, type = "l")
#png('QTLplot_parental.png', width = 20, height = 12, units = 'in', res = 300)
#Plot_result <- plot(x = CIM, QTL = QTL, type = "l")
#dev.off()

#Plot_result <- plot(x = CIM, QTL = QTL, type = "l")
#save(Plot_result, file = "Parental_Plot_result_Width.RData")

#plot(x = CIM, gen.eff= TRUE, mppData = mppData, QTL = QTL, Q.eff = "par", main = "QTL genetic effect plot")
#png('plot_parenteffect.png', width = 20, height = 12, units = 'in', res = 300)
#Plot_parents_effect <- plot(x = CIM, gen.eff= TRUE, mppData = mppData, QTL = QTL, Q.eff = "par", main = "QTL genetic effect plot")
#dev.off()

#Plot_parents_effect <- plot(x = CIM, gen.eff= TRUE, mppData = mppData, QTL = QTL, Q.eff = "par", main = "QTL genetic effect plot")
#save(Plot_parents_effect, file = "Plot_parents_effect_Width.RData")

#mpp_CV(pop.name = "HvDRR", trait.name = "Width", mppData = mppData, Q.eff = "par", her = 0.4, Rep = 1,
       #k = 3, verbose = FALSE, output.loc = tempdir())
CV <- mpp_CV(pop.name = "HvDRR", trait.name = "Width", mppData = mppData, Q.eff = "par", her = 0.4, Rep = 1,
             k = 3, verbose = FALSE, output.loc = my.loc)
save(CV, file = "Parental_CV_Width.RData")

#Perm <- mpp_perm(mppData, trait = 'Width', Q.eff = "par", N = 1000, q.val = 0.95,
         #verbose = TRUE, n.cores = 4)
#save(Perm, file = 'Perm_Width.RData')



