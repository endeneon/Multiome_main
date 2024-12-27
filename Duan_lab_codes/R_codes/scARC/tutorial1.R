# Chuxuan Li 02/23/2023
# Tutorials for WGCNA

# init ####
library(WGCNA)
setwd("~/NVME/scARC_Duan_018/case_control_DE_analysis/WGCNA/FemaleLiver-Data")
options(stringsAsFactors = F)

femData = read.csv("./LiverFemale3600.csv")
dim(femData)
names(femData)

# clean data ####
# extract expression data, transpose to row = sample, col = gene
datExpr0 = as.data.frame(t(femData[, -c(1:8)])) 
names(datExpr0) = femData$substanceBXH
rownames(datExpr0) = names(femData)[-c(1:8)]

# remove genes with past-threshold missing entries or zero variance
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# hierarchical clustering and cleaning of the samples ####
sampleTree = hclust(dist(datExpr0), method = "average")
plot(sampleTree, main = "sample clustering to detect outliers", sub = "", 
     xlab = "")
abline(h = 15, col = "red") # one sample seems to be the outlier here
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10) # retain branches with at least size 10
table(clust) # cut 1 out, keep 134 (clust is just a numerical vector)
datExpr = datExpr0[clust == 1, ] # remove the one row corresponding to sample 221
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# load traits ####
traitData = read.csv("ClinicalTraits.csv");
dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
allTraits = traitData[, -c(31, 16)];
allTraits = allTraits[, c(2, 11:36) ];
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$Mice);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();
# Re-cluster samples ####
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
save(datExpr, datTraits, file = "FemaleLiver-01-dataInput.RData")
