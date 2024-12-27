# Chuxuan Li 03/08/2023
# Visualization of networks with heatmaps

library(WGCNA)
setwd("~/NVME/scARC_Duan_018/case_control_DE_analysis/WGCNA")

# GABA 0hr ####
load("./TOMs/GABA_0hr_TOM-block.1.RData")
load("./module_colors_geneTrees_for_plot/GABA_0hr.RData")
TOM <- as.matrix(TOM)
sizeGrWindow(9, 9)
TOMplot(TOM, dendro = geneTree, Colors = moduleColors, 
        main = "Network heatmap - GABA 0hr (All genes)")
# length(moduleColors) #17028
# dim(TOM)[[1]]

# GABA 1hr ####
load("./TOMs/GABA_1hr_TOM-block.1.RData")
load("./module_colors_geneTrees_for_plot/GABA_1hr.RData")
TOM <- as.matrix(TOM)
sizeGrWindow(9, 9)
TOMplot(TOM, dendro = geneTree, Colors = moduleColors, 
        main = "Network heatmap - GABA 1hr (All genes)")

# GABA 6hr ####
load("./TOMs/GABA_6hr_TOM-block.1.RData")
load("./module_colors_geneTrees_for_plot/GABA_6hr.RData")
TOM <- as.matrix(TOM)
sizeGrWindow(9, 9)
TOMplot(TOM, dendro = geneTree, Colors = moduleColors, 
        main = "Network heatmap - GABA 6hr (All genes)")

# nmglut 0hr ####
load("./TOMs/nmglut_0hr_TOM-block.1.RData")
load("./module_colors_geneTrees_for_plot/nmglut_0hr.RData")
TOM <- as.matrix(TOM)
sizeGrWindow(9, 9)
TOMplot(TOM, dendro = geneTree, Colors = moduleColors, 
        main = "Network heatmap - nmglut 0hr (All genes)")

# nmglut 1hr ####
load("./TOMs/nmglut_1hr_TOM-block.1.RData")
load("./module_colors_geneTrees_for_plot/nmglut_1hr.RData")
TOM <- as.matrix(TOM)
sizeGrWindow(9, 9)
TOMplot(TOM, dendro = geneTree, Colors = moduleColors, 
        main = "Network heatmap - nmglut 1hr (All genes)")

# nmglut 6hr ####
load("./TOMs/nmglut_6hr_TOM-block.1.RData")
load("./module_colors_geneTrees_for_plot/nmglut_6hr.RData")
TOM <- as.matrix(TOM)
sizeGrWindow(9, 9)
TOMplot(TOM, dendro = geneTree, Colors = moduleColors, 
        main = "Network heatmap - nmglut 6hr (All genes)")