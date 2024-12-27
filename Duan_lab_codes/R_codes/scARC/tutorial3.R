# Chuxuan Li 02/27/2023

# init ###
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the expression and trait data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "FemaleLiver-02-networkConstruction-auto.RData")
lnames

# quantifying correlation between modules and traits ####
# Define numbers of genes and samples
nGenes = ncol(datExpr) #3600
nSamples = nrow(datExpr) #134
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes #total 19 modules, each has 1 value for each sample
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# plot the correlation between eigengenes (summary of modules) and traits
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(", # round up to 2 significant figures
                    signif(moduleTraitPvalue, 1), ")", sep = "") # round up to 1 significant figure
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# Gene relationship to trait and modules ####
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$weight_g) # this is one of the covariables
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3) # cut off the "ME" in front of the color names

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p")) # this is corr between individual genes and modules
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p")) # corr between genes and 1 trait
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

module = "brown" # select 1 module to look at
column = match(module, modNames)
moduleGenes = moduleColors==module # moduleColors is basically the module ident of each gene

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), # gene x module
                   abs(geneTraitSignificance[moduleGenes, 1]), # gene x trait
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
names(datExpr)
names(datExpr)[moduleColors=="brown"] # list the genes in this module

# output genes with high significance in correlation with trait in a module ####
annot = read.csv(file = "GeneAnnotation.csv")
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot], # match gene symbol with probe#
                       LocusLinkID = annot$LocusLinkID[probes2annot], # this is the Entrez code
                       moduleColor = moduleColors,
                       geneTraitSignificance, # this is all genes' corr with 1 trait
                       GSPvalue) # this is the significance of genes' corr with trait
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p"))) # modules with more positive values are ranked top 
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]) # add 2 cols from most pos module to most neg corr with the trait
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight)) # 2 orders
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")
