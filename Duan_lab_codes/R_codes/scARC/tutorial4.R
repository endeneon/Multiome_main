# Chuxuan Li 03/01/2023
# WGCNA tutorial - Interfacing network analysis with other data such as 
#functional annotation and gene ontology

# init ####
library(WGCNA)
getwd()
options(stringsAsFactors = FALSE)
# Load the expression and trait data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "FemaleLiver-02-networkConstruction-auto.RData")
lnames

# write out gene lists for external analysis ####
annot = read.csv(file = "GeneAnnotation.csv")
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH) #sort probes by their order in the annot df
# Get the corresponding Locuis Link IDs (entrez ID)
allLLIDs = annot$LocusLinkID[probes2annot]
# Choose interesting modules
intModules = c("brown", "red", "salmon")
for (module in intModules) {
  modGenes = (moduleColors==module) # get boolean vector for genes in a module
  modLLIDs = allLLIDs[modGenes] # select genes in the module, match with entrez ID
  # Write them into a file
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="")
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE) # each file is gene IDs in a module
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-all.txt", sep="")
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)

# do GO term enrichment ####
GOenr = GOenrichmentAnalysis(labels = moduleColors, # can be vector or matrix of labels (colors)
                             entrezCodes = allLLIDs, 
                             organism = "mouse", 
                             ontologies = c("BP", "CC", "MF"), 
                             backgroundType = "givenInGO", # this takes all genes in any GO categories as bg
                             nBestP = 10)

tab = GOenr$bestPTerms[[4]]$enrichment # the 4th list is the one with all 3 GO term categories
names(tab)

write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)

keepCols = c(1, 2, 5, 6, 7, 12, 13)
screenTab = tab[, keepCols]
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4)
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name")
rownames(screenTab) = NULL
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab


