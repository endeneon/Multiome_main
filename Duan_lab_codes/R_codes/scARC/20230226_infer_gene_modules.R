# Chuxuan Li 02/27/2023
# Get gene modules

# init ####
library(WGCNA)
options(stringsAsFactors = FALSE)
setwd("~/NVME/scARC_Duan_018/case_control_DE_analysis/WGCNA")
load("./count_matrices_filtered_genes_removed_outlier_samples.RData")
load("./trait_dfs.RData")
names(trait_df_list) <- names(cleaned_lst)
enableWGCNAThreads(nThreads = 24)

# GABA_0hr ####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=25, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(cleaned_lst[[1]], powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.60,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
setwd("./TOMs/")
# one step network construction and module detection

#dS = 1
# 0     1     2     3     4     5     6     7     8 
# 43 10267  4107  1092   552   522   276    99    70 

#dS = 3
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
# 32 6847 4793 1552 1108  696  436  371  364  212  203  152  127   70   65 

#dS = 4
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13 
# 333 9805 2992 1450  616  398  369  362  198  137  120  115   68   65 


net = blockwiseModules(cleaned_lst[[1]], power = 14, #soft thresholding power
                       deepSplit = 1, # sensitivity to cluster splitting
                       TOMType = "unsigned", # unsigned network
                       minModuleSize = 50,
                       #reassignThreshold = 0.1, 
                       maxBlockSize = ncol(cleaned_lst[[1]]), # how large the largest block of splitted data is
                       mergeCutHeight = 0.05, #threshold for merging modules 
                       numericLabels = T, pamRespectsDendro = F, 
                       saveTOMs = T,
                       saveTOMFileBase = paste0(names(cleaned_lst)[1], "_TOM"), 
                       verbose = 3)
table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs0 = moduleEigengenes(cleaned_lst[[1]], moduleColors)$eigengenes #total 19 modules, each has 1 value for each sample
MEs = orderMEs(MEs0)
geneTree = net$dendrograms[[1]]
save(moduleColors, geneTree, MEs, moduleLabels,
     file = "../module_colors_geneTrees_for_plot/GABA_0hr.RData")

# examine module correlation to trait
moduleTraitCor = cor(MEs, trait_df_list[[1]], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(cleaned_lst[[1]]))

# plot the correlation between eigengenes (summary of modules) and traits
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(", # round up to 2 significant figures
                    signif(moduleTraitPvalue, 1), ")", sep = "") # round up to 1 significant figure
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(trait_df_list[[1]]),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# output genes
genes = colnames(cleaned_lst[[1]])
# Choose interesting modules
table(moduleColors)
intModules = unique(moduleColors)
dir.create(path = paste0("../module_genelists/", names(cleaned_lst)[1]))
for (module in intModules) {
  modGenes = (moduleColors==module) # get boolean vector for genes in a module
  modLLIDs = genes[modGenes] # select genes in the module, match with entrez ID
  # Write them into a file
  fileName = paste("../module_genelists/", 
                   names(cleaned_lst)[1],
                   "/Genes_in_module_", module, ".txt", sep="")
  write.table(as.data.frame(modLLIDs), file = fileName, quote = F, 
              row.names = FALSE, col.names = FALSE) # each file is gene IDs in a module
}


# GABA_1hr ####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=25, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(cleaned_lst[[2]], powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.60,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# one step network construction and module detection
net = blockwiseModules(cleaned_lst[[2]], power = 14, #soft thresholding power
                       deepSplit = 3, # sensitivity to cluster splitting
                       TOMType = "unsigned", # unsigned network
                       minModuleSize = 50,
                       #reassignThreshold = 0.1, 
                       maxBlockSize = ncol(cleaned_lst[[2]]), # how large the largest block of splitted data is
                       mergeCutHeight = 0.05, #threshold for merging modules 
                       numericLabels = T, pamRespectsDendro = F, saveTOMs = T,
                       saveTOMFileBase = paste0(names(cleaned_lst)[2], "_TOM"), 
                       verbose = 3)
table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs0 = moduleEigengenes(cleaned_lst[[2]], moduleColors)$eigengenes #total 19 modules, each has 1 value for each sample
MEs = orderMEs(MEs0)
geneTree = net$dendrograms[[1]]
save(moduleColors, geneTree, MEs, moduleLabels,
     file = "../module_colors_geneTrees_for_plot/GABA_1hr.RData")

# examine module correlation to trait
moduleTraitCor = cor(MEs, trait_df_list[[2]], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(cleaned_lst[[2]]))

# plot the correlation between eigengenes (summary of modules) and traits
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(", # round up to 2 significant figures
                    signif(moduleTraitPvalue, 1), ")", sep = "") # round up to 1 significant figure
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(trait_df_list[[2]]),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# output genes 
genes = colnames(cleaned_lst[[2]])
# Choose interesting modules
table(moduleColors)
intModules = unique(moduleColors)
dir.create(path = paste0("../module_genelists/", names(cleaned_lst)[2]))
for (module in intModules) {
  modGenes = (moduleColors==module) # get boolean vector for genes in a module
  modLLIDs = genes[modGenes] # select genes in the module, match with entrez ID
  # Write them into a file
  fileName = paste("../module_genelists/", 
                   names(cleaned_lst)[2],
                   "/Genes_in_module_", module, ".txt", sep="")
  write.table(as.data.frame(modLLIDs), file = fileName, quote = F, 
              row.names = FALSE, col.names = FALSE) # each file is gene IDs in a module
}

# GABA_6hr ####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=25, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(cleaned_lst[[3]], powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.70,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# one step network construction and module detection
net = blockwiseModules(cleaned_lst[[3]], power = 14, #soft thresholding power
                       deepSplit = 2, # sensitivity to cluster splitting
                       TOMType = "unsigned", # unsigned network
                       minModuleSize = 50,
                       #reassignThreshold = 0.1, 
                       maxBlockSize = ncol(cleaned_lst[[3]]), # how large the largest block of splitted data is
                       mergeCutHeight = 0.05, #threshold for merging modules 
                       numericLabels = T, pamRespectsDendro = F, saveTOMs = T,
                       saveTOMFileBase = paste0(names(cleaned_lst)[3], "_TOM"), 
                       verbose = 3)
table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs0 = moduleEigengenes(cleaned_lst[[3]], moduleColors)$eigengenes #total 19 modules, each has 1 value for each sample
MEs = orderMEs(MEs0)
table(moduleColors)
geneTree = net$dendrograms[[1]]
save(moduleColors, geneTree, MEs, moduleLabels,
     file = "../module_colors_geneTrees_for_plot/GABA_6hr.RData")

# examine module correlation to trait
moduleTraitCor = cor(MEs, trait_df_list[[3]], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(cleaned_lst[[3]]))

# plot the correlation between eigengenes (summary of modules) and traits
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(", # round up to 2 significant figures
                    signif(moduleTraitPvalue, 1), ")", sep = "") # round up to 1 significant figure
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(trait_df_list[[3]]),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# output genes
genes = colnames(cleaned_lst[[3]])
# Choose interesting modules
intModules = unique(moduleColors)
dir.create(path = paste0("../module_genelists/", names(cleaned_lst)[3]))
for (module in intModules) {
  modGenes = (moduleColors==module) # get boolean vector for genes in a module
  modLLIDs = genes[modGenes] # select genes in the module, match with entrez ID
  # Write them into a file
  fileName = paste("../module_genelists/", 
                   names(cleaned_lst)[3],
                   "/Genes_in_module_", module, ".txt", sep="")
  write.table(as.data.frame(modLLIDs), file = fileName, quote = F, 
              row.names = FALSE, col.names = FALSE) # each file is gene IDs in a module
}

# nmglut_0hr ####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=25, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(cleaned_lst$nmglut_0hr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.70,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# one step network construction and module detection
net = blockwiseModules(cleaned_lst$nmglut_0hr, power = 20, #soft thresholding power
                       deepSplit = 3, # sensitivity to cluster splitting
                       TOMType = "unsigned", # unsigned network
                       minModuleSize = 50,
                       #reassignThreshold = 0.1, 
                       maxBlockSize = ncol(cleaned_lst$nmglut_0hr), # how large the largest block of splitted data is
                       mergeCutHeight = 0.05, #threshold for merging modules 
                       numericLabels = T, pamRespectsDendro = F, saveTOMs = T,
                       saveTOMFileBase = paste0(names(cleaned_lst)[4], "_TOM"), 
                       verbose = 3)
table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs0 = moduleEigengenes(cleaned_lst$nmglut_0hr, moduleColors)$eigengenes #total 19 modules, each has 1 value for each sample
MEs = orderMEs(MEs0)
geneTree = net$dendrograms[[1]]
save(moduleColors, geneTree, MEs, moduleLabels,
     file = "../module_colors_geneTrees_for_plot/nmglut_0hr.RData")

# examine module correlation to trait
moduleTraitCor = cor(MEs, trait_df_list$nmglut_0hr, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(cleaned_lst$nmglut_0hr))

# plot the correlation between eigengenes (summary of modules) and traits
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(", # round up to 2 significant figures
                    signif(moduleTraitPvalue, 1), ")", sep = "") # round up to 1 significant figure
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(trait_df_list$nmglut_0hr),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# output genes
genes = colnames(cleaned_lst$nmglut_0hr)
# Choose interesting modules
table(moduleColors)
intModules = unique(moduleColors)
dir.create(path = paste0("../module_genelists/", names(cleaned_lst)[4]))
for (module in intModules) {
  modGenes = (moduleColors==module) # get boolean vector for genes in a module
  modLLIDs = genes[modGenes] # select genes in the module, match with entrez ID
  # Write them into a file
  fileName = paste("../module_genelists/", 
                   names(cleaned_lst)[4],
                   "/Genes_in_module_", module, ".txt", sep="")
  write.table(as.data.frame(modLLIDs), file = fileName, quote = F, 
              row.names = FALSE, col.names = FALSE) # each file is gene IDs in a module
}

# nmglut_1hr ####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=50, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(cleaned_lst$nmglut_1hr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.75,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# one step network construction and module detection
net = blockwiseModules(cleaned_lst$nmglut_1hr, power = 30, #soft thresholding power
                       deepSplit = 4, # sensitivity to cluster splitting
                       TOMType = "unsigned", # unsigned network
                       minModuleSize = 50,
                       #reassignThreshold = 0.1, 
                       maxBlockSize = ncol(cleaned_lst$nmglut_1hr), # how large the largest block of splitted data is
                       mergeCutHeight = 0.05, #threshold for merging modules 
                       numericLabels = T, pamRespectsDendro = F, saveTOMs = T,
                       saveTOMFileBase = paste0(names(cleaned_lst)[5], "_TOM"), 
                       verbose = 3)
table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs0 = moduleEigengenes(cleaned_lst$nmglut_1hr, moduleColors)$eigengenes #total 19 modules, each has 1 value for each sample
MEs = orderMEs(MEs0)
geneTree = net$dendrograms[[1]]
save(moduleColors, geneTree, MEs, moduleLabels,
     file = "../module_colors_geneTrees_for_plot/nmglut_1hr.RData")

# examine module correlation to trait
moduleTraitCor = cor(MEs, trait_df_list$nmglut_1hr, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(cleaned_lst$nmglut_1hr))

# plot the correlation between eigengenes (summary of modules) and traits
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(", # round up to 2 significant figures
                    signif(moduleTraitPvalue, 1), ")", sep = "") # round up to 1 significant figure
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(trait_df_list$nmglut_1hr),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# output genes 
genes = colnames(cleaned_lst$nmglut_1hr)
# Choose interesting modules
table(moduleColors)
intModules = unique(moduleColors)
dir.create(path = paste0("../module_genelists/", names(cleaned_lst)[5]))
for (module in intModules) {
  modGenes = (moduleColors==module) # get boolean vector for genes in a module
  modLLIDs = genes[modGenes] # select genes in the module, match with entrez ID
  # Write them into a file
  fileName = paste("../module_genelists/", 
                   names(cleaned_lst)[5],
                   "/Genes_in_module_", module, ".txt", sep="")
  write.table(as.data.frame(modLLIDs), file = fileName, quote = F, 
              row.names = FALSE, col.names = FALSE) # each file is gene IDs in a module
}

# nmglut_6hr ####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=25, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(cleaned_lst$nmglut_6hr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# one step network construction and module detection
net = blockwiseModules(cleaned_lst$nmglut_6hr, power = 16, #soft thresholding power
                       deepSplit = 2, # sensitivity to cluster splitting
                       TOMType = "unsigned", # unsigned network
                       minModuleSize = 50,
                       #reassignThreshold = 0.1, 
                       maxBlockSize = ncol(cleaned_lst$nmglut_6hr), # how large the largest block of splitted data is
                       mergeCutHeight = 0.05, #threshold for merging modules 
                       numericLabels = T, pamRespectsDendro = F, saveTOMs = T,
                       saveTOMFileBase = paste0(names(cleaned_lst)[6], "_TOM"), 
                       verbose = 3)
table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs0 = moduleEigengenes(cleaned_lst$nmglut_6hr, moduleColors)$eigengenes #total 19 modules, each has 1 value for each sample
MEs = orderMEs(MEs0)
geneTree = net$dendrograms[[1]]
save(moduleColors, geneTree, MEs, moduleLabels,
     file = "../module_colors_geneTrees_for_plot/nmglut_6hr.RData")

# examine module correlation to trait 
moduleTraitCor = cor(MEs, trait_df_list$nmglut_6hr, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(cleaned_lst$nmglut_6hr))

# plot the correlation between eigengenes (summary of modules) and traits
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(", # round up to 2 significant figures
                    signif(moduleTraitPvalue, 1), ")", sep = "") # round up to 1 significant figure
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(trait_df_list$nmglut_6hr),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# output genes
genes = colnames(cleaned_lst$nmglut_6hr)
# Choose interesting modules
table(moduleColors)
intModules = unique(moduleColors)
dir.create(path = paste0("../module_genelists/", names(cleaned_lst)[6]))
for (module in intModules) {
  modGenes = (moduleColors==module) # get boolean vector for genes in a module
  modLLIDs = genes[modGenes] # select genes in the module, match with entrez ID
  # Write them into a file
  fileName = paste("../module_genelists/", 
                   names(cleaned_lst)[6],
                   "/Genes_in_module_", module, ".txt", sep="")
  write.table(as.data.frame(modLLIDs), file = fileName, quote = F, 
              row.names = FALSE, col.names = FALSE) # each file is gene IDs in a module
}

# npglut_0hr ####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(cleaned_lst$npglut_0hr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.60,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# one step network construction and module detection
net = blockwiseModules(cleaned_lst$npglut_0hr, power = 26, #soft thresholding power
                       deepSplit = 3, # sensitivity to cluster splitting
                       TOMType = "unsigned", # unsigned network
                       minModuleSize = 50,
                       #reassignThreshold = 0.1, 
                       maxBlockSize = ncol(cleaned_lst$npglut_0hr), # how large the largest block of splitted data is
                       mergeCutHeight = 0.05, #threshold for merging modules 
                       numericLabels = T, pamRespectsDendro = F, saveTOMs = T,
                       saveTOMFileBase = paste0(names(cleaned_lst)[7], "_TOM"), 
                       verbose = 3)
table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs0 = moduleEigengenes(cleaned_lst$npglut_0hr, moduleColors)$eigengenes #total 19 modules, each has 1 value for each sample
MEs = orderMEs(MEs0)
geneTree = net$dendrograms[[1]]
save(moduleColors, geneTree, MEs, moduleLabels,
     file = "../module_colors_geneTrees_for_plot/npglut_0hr.RData")

# examine module correlation to trait
moduleTraitCor = cor(MEs, trait_df_list$npglut_0hr, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(cleaned_lst$npglut_0hr))

# plot the correlation between eigengenes (summary of modules) and traits
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(", # round up to 2 significant figures
                    signif(moduleTraitPvalue, 1), ")", sep = "") # round up to 1 significant figure
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(trait_df_list$npglut_0hr),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# output genes
genes = colnames(cleaned_lst$npglut_0hr)
# Choose interesting modules
table(moduleColors)
intModules = unique(moduleColors)
dir.create(path = paste0("../module_genelists/", names(cleaned_lst)[7]))
for (module in intModules) {
  modGenes = (moduleColors==module) # get boolean vector for genes in a module
  modLLIDs = genes[modGenes] # select genes in the module, match with entrez ID
  # Write them into a file
  fileName = paste("../module_genelists/", 
                   names(cleaned_lst)[7],
                   "/Genes_in_module_", module, ".txt", sep="")
  write.table(as.data.frame(modLLIDs), file = fileName, quote = F, 
              row.names = FALSE, col.names = FALSE) # each file is gene IDs in a module
}

# npglut_1hr ####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(cleaned_lst$npglut_1hr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.6,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# one step network construction and module detection
net = blockwiseModules(cleaned_lst$npglut_1hr, power = 30, #soft thresholding power
                       deepSplit = 4, # sensitivity to cluster splitting
                       TOMType = "unsigned", # unsigned network
                       minModuleSize = 50,
                       #reassignThreshold = 0.1, 
                       maxBlockSize = ncol(cleaned_lst$npglut_1hr), # how large the largest block of splitted data is
                       mergeCutHeight = 0.05, #threshold for merging modules 
                       numericLabels = T, pamRespectsDendro = F, saveTOMs = T,
                       saveTOMFileBase = paste0(names(cleaned_lst)[8], "_TOM"), 
                       verbose = 3)
table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs0 = moduleEigengenes(cleaned_lst$npglut_1hr, moduleColors)$eigengenes #total 19 modules, each has 1 value for each sample
MEs = orderMEs(MEs0)
geneTree = net$dendrograms[[1]]
save(moduleColors, geneTree, MEs, moduleLabels,
     file = "../module_colors_geneTrees_for_plot/npglut_1hr.RData")

# examine module correlation to trait
moduleTraitCor = cor(MEs, trait_df_list$npglut_1hr, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(cleaned_lst$npglut_1hr))

# plot the correlation between eigengenes (summary of modules) and traits
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(", # round up to 2 significant figures
                    signif(moduleTraitPvalue, 1), ")", sep = "") # round up to 1 significant figure
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(trait_df_list$npglut_1hr),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# output genes 
genes = colnames(cleaned_lst$npglut_1hr)
# Choose interesting modules
table(moduleColors)
intModules = unique(moduleColors)
dir.create(path = paste0("../module_genelists/", names(cleaned_lst)[8]))
for (module in intModules) {
  modGenes = (moduleColors==module) # get boolean vector for genes in a module
  modLLIDs = genes[modGenes] # select genes in the module, match with entrez ID
  # Write them into a file
  fileName = paste("../module_genelists/", 
                   names(cleaned_lst)[8],
                   "/Genes_in_module_", module, ".txt", sep="")
  write.table(as.data.frame(modLLIDs), file = fileName, quote = F, 
              row.names = FALSE, col.names = FALSE) # each file is gene IDs in a module
}

# npglut_6hr ####
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(cleaned_lst$npglut_6hr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.7,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# one step network construction and module detection
net = blockwiseModules(cleaned_lst$npglut_6hr, power = 30, #soft thresholding power
                       deepSplit = 4, # sensitivity to cluster splitting
                       TOMType = "unsigned", # unsigned network
                       minModuleSize = 50,
                       #reassignThreshold = 0.1, 
                       maxBlockSize = ncol(cleaned_lst$npglut_6hr), # how large the largest block of splitted data is
                       mergeCutHeight = 0.05, #threshold for merging modules 
                       numericLabels = T, pamRespectsDendro = F, saveTOMs = T,
                       saveTOMFileBase = paste0(names(cleaned_lst)[9], "_TOM"), 
                       verbose = 3)
table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs0 = moduleEigengenes(cleaned_lst$npglut_6hr, moduleColors)$eigengenes #total 19 modules, each has 1 value for each sample
MEs = orderMEs(MEs0)
geneTree = net$dendrograms[[1]]
save(moduleColors, geneTree, MEs, moduleLabels,
     file = "../module_colors_geneTrees_for_plot/npglut_6hr.RData")

# examine module correlation to trait 
moduleTraitCor = cor(MEs, trait_df_list$npglut_6hr, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(cleaned_lst$npglut_6hr))

# plot the correlation between eigengenes (summary of modules) and traits
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(", # round up to 2 significant figures
                    signif(moduleTraitPvalue, 1), ")", sep = "") # round up to 1 significant figure
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(trait_df_list$npglut_6hr),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# output genes
genes = colnames(cleaned_lst$npglut_6hr)
# Choose interesting modules
table(moduleColors)
intModules = unique(moduleColors)
dir.create(path = paste0("../module_genelists/", names(cleaned_lst)[9]))
for (module in intModules) {
  modGenes = (moduleColors==module) # get boolean vector for genes in a module
  modLLIDs = genes[modGenes] # select genes in the module, match with entrez ID
  # Write them into a file
  fileName = paste("../module_genelists/", 
                   names(cleaned_lst)[9],
                   "/Genes_in_module_", module, ".txt", sep="")
  write.table(as.data.frame(modLLIDs), file = fileName, quote = F, 
              row.names = FALSE, col.names = FALSE) # each file is gene IDs in a module
}
