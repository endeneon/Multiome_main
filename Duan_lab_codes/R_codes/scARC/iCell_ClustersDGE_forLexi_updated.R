##DIFFERENTIAL GENE EXPRESSION WITH MONOCLE

library(Seurat)
library(Matrix)
library(methods)
library(ggplot2)
library(cowplot)
library(reticulate)

rm(list=ls())
setwd("") # DIRECTORY WHERE R OBJECT IS
load("FILENAME.Rda") # THIS IS WHAT's DOWNLOADED FROM GEO

# save csv of cluster assignments and raw data:

# cluster information
write.table((iCell.integrated@meta.data$integrated_snn_res.0.5),
	file = "FILENAME_clusters.csv",
	sep = ',',
	row.names = rownames(iCell.integrated@meta.data),
	col.names = NA)

# raw data
tmpData = iCell.integrated@assays$RNA@counts
tmpData = as.matrix(tmpData)
tmpData = as.data.frame(tmpData)
write.table(tmpData, file = "FILENAME_raw-data.csv", sep = ',', row.names = rownames(tmpData), col.names = NA)

#===================================================================================================

library(devtools)
library(dplyr)
library(useful)
library(ggplot2)
library(stringr)
library(Matrix)
library(xlsx)
library(reshape)
library(qlcMatrix)
library(monocle)
library(DESeq2)
library(limma)
packageVersion("monocle") #‘2.6.4’

rm(list=ls())
setwd("FILEPATH") # to where you saved _clusters.csv and _raw-data.csv

#load cluster info
clusters = read.table('FILENAME_clusters.csv', header = T, fill = T, row.names = 1, sep = ',')
cellNames = rownames(clusters)

#load dataset
data = read.table('FILENAME_raw-data.csv', header = T, fill = T, row.names = 1, sep =',')
colnames(data) = cellNames
geneNames = rownames(data)

iCell <- substr(cellNames, start = 1, stop = 6)
iCell[grep('0hr',iCell)] = '0hr'
iCell[grep('1hr',iCell)] = '1hr'
iCell[grep('2hr',iCell)] = '2hr'
iCell[grep('4hr',iCell)] = '4hr'
iCell <- as.data.frame(iCell, stringsAsFactors=FALSE)
rownames(iCell) = cellNames

#CHANGE THE SAVE DIRECTORY
dir.create("SAVE-FOLDER")
setwd("PATH/SAVE-FOLDER")
for (timept in c("1hr", "2hr", "4hr")) {

  iCell.subset <- subset(iCell, iCell=='0hr' | iCell == timept)
  names <- rownames(iCell.subset)
  clusters.subset <- subset(clusters, rownames(clusters) %in% names)
  clusters.subset$x = factor(clusters.subset$x)
  clusters.subset$cellType = factor(iCell.subset$iCell)

  #for each cluster
  for(i in as.numeric(levels(clusters.subset$x))) {
    
    #gets all the time points for a cluster i; factor with two levels
    cellType = clusters.subset$cellType[clusters.subset$x == i]
    #gets the cellNames of each cell in the cluster; character vector of cell names
    cellNameiCells = rownames(clusters.subset[clusters.subset$x == i,])
    #looks for those cells in the full dataset; dataframe (rows are genes, columns are cells)
    dataiCells = data[,cellNameiCells]
    #dataframe of each cell with its time point
    cellT = data.frame(stim = cellType, row.names = cellNameiCells)
    pd = new('AnnotatedDataFrame', data = cellT)
    geneType = data.frame(type = rep('protein_coding', length(geneNames)), gene_short_name = geneNames, row.names = geneNames)
    fd = new('AnnotatedDataFrame', data = geneType)

    neurData = newCellDataSet(as(as.matrix(dataiCells), 'sparseMatrix'), phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
    neurData = estimateSizeFactors(neurData)
    neurData = estimateDispersions(neurData)


    pData(neurData)$Total_mRNAs = colSums(as.matrix(exprs(neurData)))
    fData(neurData)$num_cells_expressed = rowSums(as.matrix(exprs(neurData)>0)) #
    fData(neurData)$perc_cells_expressed =rowSums(as.matrix(exprs(neurData))>0)/ncol(as.matrix(exprs(neurData))) #
    expressed_genes = row.names(subset(fData(neurData), perc_cells_expressed >= 0.01)) #CAN INTRODUCE CUTOFF
    print(length(expressed_genes))
    
    diff_test_res = differentialGeneTest(neurData[expressed_genes,], fullModelFormulaStr = '~stim', cores = 10) #run parallel
    sig_genes = subset(diff_test_res, qval < 0.05)
    #Monocle calculates qval as such: p.adjust(subset(diff_test_res, status == 'OK')[, 'pval'], method="BH")
    #Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57, 289–300. http://www.jstor.org/stable/2346101.
    #same as FDR in R

    #calculate fold change
    mean0hr = rowMeans(as.matrix(exprs(neurData[expressed_genes,pData(neurData)$stim == '0hr'])))
    mean_timept = rowMeans(as.matrix(exprs(neurData[expressed_genes,pData(neurData)$stim == timept])))
    FC = (mean_timept/mean0hr)
    log2FC = log2(mean_timept/mean0hr)
    FC_sig = FC[rownames(sig_genes)]
    log2FC_sig = log2FC[rownames(sig_genes)]
    sig_genes$FC = FC_sig
    sig_genes$log2FC = log2FC_sig


    sig_genes = sig_genes[order(-sig_genes$FC),] #CAN INTRODUCE CUTOFF
    foldChange = 1.5
    sig_genes_inducible = sig_genes[sig_genes$FC > foldChange,]
    # sig_genes_decreasing = sig_genes[sig_genes$FC < foldChange,] #not a good measure of decreasing
    sig_genes_decreasing = sig_genes[sig_genes$FC < 1,]

    write.csv(sig_genes_inducible,file = paste('Cluster_',i,'_',timept,'_vs_0hr_inducible.csv', sep = ''))
    write.csv(sig_genes_decreasing,file = paste('Cluster_',i,'_',timept,'_vs_0hr_decreasing.csv', sep = '')) #sep for file name, not file format
    write.csv(sig_genes,file = paste('Cluster_',i,'_',timept,'_vs_0hr_all.csv', sep = ''))

  }
}


#===================================================================================================

##REDO with no cutoffs so we can see all the genes and their expression profiles

library(devtools)
library(dplyr)
library(useful)
library(ggplot2)
library(stringr)
library(Matrix)
library(xlsx)
library(reshape)
library(qlcMatrix)
library(monocle)
library(DESeq2)
library(limma)
packageVersion("monocle") #‘2.6.4’

rm(list=ls())
setwd("FILEPATH") # to where you saved _clusters.csv and _raw-data.csv

#load cluster info
clusters = read.table('FILENAME_clusters.csv', header = T, fill = T, row.names = 1, sep = ',')
cellNames = rownames(clusters)

#load dataset
data = read.table('FILENAME_raw-data.csv', header = T, fill = T, row.names = 1, sep =',')
colnames(data) = cellNames

geneNames = rownames(data)

iCell <- substr(cellNames, start = 1, stop = 6)
iCell[grep('0hr',iCell)] = '0hr'
iCell[grep('1hr',iCell)] = '1hr'
iCell[grep('2hr',iCell)] = '2hr'
iCell[grep('4hr',iCell)] = '4hr'
iCell <- as.data.frame(iCell, stringsAsFactors=FALSE)
rownames(iCell) = cellNames

#CHANGE THE SAVE DIRECTORY
dir.create("SAVE-FOLDER")
setwd("PATH/SAVE-FOLDER")
for (timept in c("1hr", "2hr", "4hr")) {

  iCell.subset <- subset(iCell, iCell=='0hr' | iCell == timept)
  names <- rownames(iCell.subset)
  clusters.subset <- subset(clusters, rownames(clusters) %in% names)
  clusters.subset$x = factor(clusters.subset$x)
  clusters.subset$cellType = factor(iCell.subset$iCell)

  #for each cluster
  for(i in as.numeric(levels(clusters.subset$x))) {
    
    #gets all the time points for a cluster i; factor with two levels
    cellType = clusters.subset$cellType[clusters.subset$x == i]
    #gets the cellNames of each cell in the cluster; character vector of cell names
    cellNameiCells = rownames(clusters.subset[clusters.subset$x == i,])
    #looks for those cells in the full dataset; dataframe (rows are genes, columns are cells)
    dataiCells = data[,cellNameiCells]
    #dataframe of each cell with its time point
    cellT = data.frame(stim = cellType, row.names = cellNameiCells)
    pd = new('AnnotatedDataFrame', data = cellT)
    geneType = data.frame(type = rep('protein_coding', length(geneNames)), gene_short_name = geneNames, row.names = geneNames)
    fd = new('AnnotatedDataFrame', data = geneType)

    neurData = newCellDataSet(as(as.matrix(dataiCells), 'sparseMatrix'), phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
    neurData = estimateSizeFactors(neurData)
    neurData = estimateDispersions(neurData)


    pData(neurData)$Total_mRNAs = colSums(as.matrix(exprs(neurData)))
    fData(neurData)$num_cells_expressed = rowSums(as.matrix(exprs(neurData)>0)) #per cluster!
    fData(neurData)$perc_cells_expressed =rowSums(as.matrix(exprs(neurData))>0)/ncol(as.matrix(exprs(neurData))) #per cluster!
    expressed_genes = row.names(subset(fData(neurData), perc_cells_expressed >= 0)) #CAN INTRODUCE CUTOFF
    print(length(expressed_genes))
    
    diff_test_res = differentialGeneTest(neurData[expressed_genes,], fullModelFormulaStr = '~stim', cores = 10) #run parallel
    sig_genes = subset(diff_test_res, qval <= 1)

    #calculate fold change
    mean0hr = rowMeans(as.matrix(exprs(neurData[expressed_genes,pData(neurData)$stim == '0hr'])))
    mean_timept = rowMeans(as.matrix(exprs(neurData[expressed_genes,pData(neurData)$stim == timept])))
    FC = log2(mean_timept/mean0hr)
    FC_sig = FC[rownames(sig_genes)]
    sig_genes$FC = FC_sig


    sig_genes = sig_genes[order(-sig_genes$FC),] #CAN INTRODUCE CUTOFF
    # foldChange = 0
    # sig_genes_inducible = sig_genes[sig_genes$FC > foldChange,]
    # sig_genes_decreasing = sig_genes[sig_genes$FC < foldChange,]

    # write.csv(sig_genes_inducible,file = paste('Cluster_',i,'_',timept,'_vs_0hr_inducible.csv', sep = ''))
    # write.csv(sig_genes_decreasing,file = paste('Cluster_',i,'_',timept,'_vs_0hr_decreasing.csv', sep = '')) #sep for file name, not file format
    write.csv(sig_genes,file = paste('Cluster_',i,'_',timept,'_vs_0hr_all.csv', sep = ''))

  }
}



#--------------------------------------

#looking at the DGE data


setwd("PATH/SAVE-FOLDER")

#change this vector
x = c("MMP1")

name <- c()
timepoint <- c()
cluster <- c()
log2FC <- c()
FC <- c()
pval <- c()
q.val <- c()
num.cells.exp = c()
per.cells.exp = c()

for (g in x) {

  for (timept in c("1hr", "2hr", "4hr")) {
    for(i in 0:14) {
      file.name = paste('Cluster_',i,'_',timept,'_vs_0hr_all.csv', sep = '')
      var.name = paste('Cluster_',i,'_',timept, sep = '')
      assign(var.name, read.csv(file.name, stringsAsFactors = FALSE, row.names = 1))
      if (g %in% get(var.name)$gene_short_name) {
        name <- c(name, g)
        timepoint <- c(timepoint, timept)
        cluster <- c(cluster, i)
        log2FC <- c(log2FC, get(var.name)[g,]$FC)
        FC <- c(FC, (2^get(var.name)[g,]$FC))
        pval <- c(pval, get(var.name)[g,]$pval)
        q.val <- c(q.val, get(var.name)[g,]$qval)
        num.cells.exp = c(num.cells.exp, get(var.name)[g,]$num_cells_expressed)
        per.cells.exp = c(per.cells.exp, get(var.name)[g,]$perc_cells_expressed)

      }
    }
  }
}

table <- cbind(name, timepoint, cluster, log2FC, FC, pval, q.val, num.cells.exp, per.cells.exp)
print(table)



#===================================================================================================
#2019-05-08
#Making geom_point ClusterDGE plots with "percent of cells in the cluster expressing the gene" as the size of the point

library(ggplot2)
library(plyr)
library(scales)
library(zoo)
library(reshape)
rm(list=ls())
setwd("PATH/SAVE-FOLDER")

all_FC <- c()
all_percell <- c()

cat <- matrix(nrow=15, ncol=0)
cat2 <- matrix(nrow=15, ncol=0)
genes <- c()

for (g in c("FOS","LINC00473","BDNF", "MMP1", "ZNF331", "NPTX2", "CACNB2", "DHCR7", "CIB2")) {
#for (g in c("SIK1","DHCR7","CIB2","NTRK2", "UBE3C", "SOX5", "FMR1", "KCNQ3", "PTCHD1", "SHANK3")) {
#for (g in c("EGR1","ETF1","CACNB2","PPP1R13B", "MAN2A1", "SETD8", "NAGA", "ZSWIM6", "CLP1", "ATP2A2")) {
#for (g in c("SST","TAC1","EBF1","NKX2.1", "SCGN", "ROBO2", "GABRB3", "SLC6A17"))  {
#for (g in c("XIRP1","LINC00602","SIK1","DHCR7", "HCN1", "EGR1", "CACNB2", "CYP2D6", "CNTN4"))  {
#for (g in c("AHI1","CIB2","CNR1","DHCR7", "FMR1", "GRID1", "NFIX", "NTRK2", "PER2", "PPM1D", "PRICKLE2", "SIK1", "SOX5", "SLITRK5", "UBE3C", "WASF1")) {
#for (g in c("TAF5", "PITPNM2", "FURIN","PPP1R13B", "RFTN2", "NAGA", "LUZP2", "EGR1", "ETF1", "MAN2A1")) {
#for (g in c("LINC01226", "ARL4D", "DUSP1", "KLF10", "MSN", "NTRK2", "RND3", "VIM", "ZNF184")) {
#for (g in c("INA", "TAF5", "IMMP2L", "ZNF804A", "CACNB2", "NAB2", "R3HDM2", "STAT6", "CREB3L1", "DGKZ", "ACTR5", "PPP1R16B", "MARS2", "RFTN2", "EP300", "GNL3", "SMIM4", "STAG1", "GRM3", "ANP32E", "SNAP91", "SRPK2", "CLCN3", "NAGA", "CLP1", "LUZP2", "PTN", "ZNF536", "HSPA9", "KDM3B", "GPM6A", "MMP16", "ENKD1", "NFATC3", "PLA2G15", "RANBP10", "SLC7A6", "MYO15A", "MAN2A1", "TMTC1", "PRR12")) {
#for (g in c("AGAP2", "ANK2", "ARID1B", "BAZ2B", "BTAF1", "CACNB2", "CAMK2B", "CCNK", "CDK13", "CIB2", "CNOT3", "DDX3X", "DOCK8", "DPP10", "EEF1A2", "ELP4", "EP300", "FOXG1", "HIVEP3", "ILF2", "KDM6B", "KMT2A", "MBOAT7", "MED13", "MED13L", "NFIX", "OXTR", "PCDH19", "PER2", "PRICKLE1", "PRICKLE2", "PRR12", "RLIM", "RORA", "SHANK2", "SHANK3", "SLC12A5", "SLC6A1", "SLC9A6", "SMAD4", "SNX14", "SLITRK5", "STAG1", "STXBP1", "TBC1D23", "TET2", "TRIO", "TRIP12", "TTN", "UBR5", "USP15", "USP7", "ZBTB20", "ZNF804A")) {
    cat(g, "\n")
    #stores Fold Change
    mat <- matrix(nrow = 15, ncol = 3)
    #stores percentage of cells in the cluster that express the gene
    mat2 <- matrix(nrow = 15, ncol = 3)
    one <- paste(g, "1hr")
    two <- paste(g, "2hr")
    four <- paste(g, "4hr")
    colnames(mat) <- c(one, two, four)
    rownames(mat) <- c(0:14)
    colnames(mat2) <- c(one, two, four)
    rownames(mat2) <- c(0:14)

    for (timept in c("1hr", "2hr", "4hr")) {
        for(i in 0:14) {
          file.name = paste('Cluster_',i,'_',timept,'_vs_0hr_inducible.csv', sep = '')
          # file.name = paste('Cluster_',i,'_',timept,'_vs_0hr_all.csv', sep = '')
          var.name = paste('Cluster_',i,'_',timept, sep = '')
          assign(var.name, read.csv(file.name, stringsAsFactors = FALSE, row.names = 1))
          if (g %in% get(var.name)$gene_short_name) {
            mat[(i+1), paste(g,timept)] <- as.numeric(get(var.name)[g,]$log2FC)
            mat2[(i+1), paste(g,timept)] <- as.numeric(get(var.name)[g,]$perc_cells_expressed)
            } else {
                mat[(i+1), paste(g,timept)] <- 0
                mat2[(i+1), paste(g,timept)] <- 0
            }
            all_percell <- c(all_percell, get(var.name)$perc_cells_expressed)
            all_FC <- c(all_FC, get(var.name)$log2FC)
        }
    }
    cat <- cbind(cat, mat)
    cat2 <- cbind(cat2, mat2)
}

all_FC <- all_FC[is.finite(all_FC)]
range_FC <- max(range(all_FC)) - min(range(all_FC))
range_percell <- max(range(all_percell)) - min(range(all_percell))

dat <- melt(cat)
dat2 <- melt(cat2)

dat$percell <- dat2$value * 100
dat$percell[dat$percell == 0] <- NA
dat$X1=as.factor(dat$X1)

dat$timept <- NA
dat$timept[grep("1hr", dat$X2)] <- "1hr"
dat$timept[grep("2hr", dat$X2)] <- "2hr"
dat$timept[grep("4hr", dat$X2)] <- "4hr"

dat$cluster <- NA
for (i in 0:14) {
  dat$cluster[grep(i, dat$X1)] <- paste("Cluster ", i)
}
dat$cluster <- as.factor(dat$cluster)
dat$cluster=factor(dat$cluster, levels=levels(dat$cluster)[rev(c(1:2,8:15,3:7))])


dat$gene <- NA
#for (g in c("AGAP2", "ANK2", "ARID1B", "BAZ2B", "BTAF1", "CACNB2", "CAMK2B", "CCNK", "CDK13", "CIB2", "CNOT3", "DDX3X", "DOCK8", "DPP10", "EEF1A2", "ELP4", "EP300", "FOXG1", "HIVEP3", "ILF2", "KDM6B", "KMT2A", "MBOAT7", "MED13", "MED13L", "NFIX", "OXTR", "PCDH19", "PER2", "PRICKLE1", "PRICKLE2", "PRR12", "RLIM", "RORA", "SHANK2", "SHANK3", "SLC12A5", "SLC6A1", "SLC9A6", "SMAD4", "SNX14", "SLITRK5", "STAG1", "STXBP1", "TBC1D23", "TET2", "TRIO", "TRIP12", "TTN", "UBR5", "USP15", "USP7", "ZBTB20", "ZNF804A")) {
#for (g in c("INA", "TAF5", "IMMP2L", "ZNF804A", "CACNB2", "NAB2", "R3HDM2", "STAT6", "CREB3L1", "DGKZ", "ACTR5", "PPP1R16B", "MARS2", "RFTN2", "EP300", "GNL3", "SMIM4", "STAG1", "GRM3", "ANP32E", "SNAP91", "SRPK2", "CLCN3", "NAGA", "CLP1", "LUZP2", "PTN", "ZNF536", "HSPA9", "KDM3B", "GPM6A", "MMP16", "ENKD1", "NFATC3", "PLA2G15", "RANBP10", "SLC7A6", "MYO15A", "MAN2A1", "TMTC1", "PRR12")) {
for (g in c("FOS","LINC00473","BDNF", "MMP1", "ZNF331", "NPTX2", "CACNB2", "DHCR7", "CIB2")) {
#for (g in c("SIK1","DHCR7","CIB2","NTRK2", "UBE3C", "SOX5", "FMR1", "KCNQ3", "PTCHD1", "SHANK3")) {
#for (g in c("EGR1","ETF1","CACNB2","PPP1R13B", "MAN2A1", "SETD8", "NAGA", "ZSWIM6", "CLP1", "ATP2A2")) {
#for (g in c("XIRP1","LINC00602","SIK1","DHCR7", "HCN1", "EGR1", "CACNB2", "CYP2D6", "CNTN4")) {
#for (g in c("AHI1","CIB2","CNR1","DHCR7", "FMR1", "GRID1", "NFIX", "NTRK2", "PER2", "PPM1D", "PRICKLE2", "SIK1", "SOX5", "SLITRK5", "UBE3C", "WASF1")) {
#for (g in c("TAF5", "PITPNM2", "FURIN","PPP1R13B", "RFTN2", "NAGA", "LUZP2", "EGR1", "ETF1", "MAN2A1")) {
#for (g in c("LINC01226", "ARL4D", "DUSP1", "KLF10", "MSN", "NTRK2", "RND3", "VIM", "ZNF184")) {
  dat$gene[grep(g, dat$X2)] <- paste(g)
}

dat$gene <- as.factor(dat$gene)
#dat$gene=factor(dat$gene, levels=levels(dat$gene)[c(1:2,5,3:4,6:7)])

#dat$value = as.numeric(dat$value)
dat$value = as.numeric(dat$value)
#dat$value <- 2^(as.numeric(dat$value))
dat$value[is.infinite(dat$value)] <- max(range(all_FC))
dat$value[dat$value <= 0] <- NA
dat$percell[is.na(dat$value)] <- NA

DGEplot <- ggplot(dat, aes(timept, cluster)) + 
geom_point(aes(size = percell, color = value)) + 
scale_colour_gradient(low="goldenrod1", high = "red") +
facet_grid(rows = NULL, vars(dat$gene), drop =F) + 
#scale_fill_manual(labels = c("No", "Yes"), values = c("gray95","black")) +
labs(x="Timepoint",
 y="Cluster",
 title = "Differential Gene Expression Analysis", 
 subtitle="Fold-change cutoff of 1.5x"
   #fill = "Inducible?"
    ) +
scale_size_continuous(range = c(3,12)) +
theme(panel.grid.major =   element_line(colour = "gray90",size=0.1),
  panel.grid.major.y = element_line(colour = "gray90",size=0.1),
  panel.grid.minor.y = element_blank(),
  panel.background = element_rect(fill = "gray95"),
  plot.background = element_blank(),
  #strip.text = element_text(size = 25)
  #panel.spacing.x = unit(4, "lines")

  #strip.background = element_blank(),
)
DGEplot

ggsave("FILENAME_geompoint1.pdf", plot = DGEplot, useDingbats=FALSE, width = 15, height = 8, units = "in")
ggsave("FILENAME_Schiz108.pdf", plot = DGEplot, useDingbats=FALSE, width = 15, height = 8, units = "in")
ggsave("FILENAME_SFARI.pdf", plot = DGEplot, useDingbats=FALSE, width = 15, height = 8, units = "in")


#===================================================================================================

#Privately induced genes (only induced in 1 cluster)

rm(list=ls())
library(plyr)
library(dplyr)
library(tidyr)
setwd("PATH/SAVE-FOLDER")

combined.genes <- c()

for(i in 0:14) {
  for (timept in c("1hr", "2hr", "4hr")) {
    file.name = paste('Cluster_',i,'_',timept,'_vs_0hr_inducible.csv', sep = '')
    var.name = paste('Cluster_',i,'_',timept, sep = '')
    assign(var.name, read.csv(file.name, stringsAsFactors = FALSE, row.names = 1))
  }

  induc.name = paste('Cluster_',i,'_inducible', sep = '')
  onehr = get(paste('Cluster_',i,'_1hr', sep = ''))
  twohr = get(paste('Cluster_',i,'_2hr', sep = ''))
  fourhr = get(paste('Cluster_',i,'_4hr', sep = ''))
  full.list = join_all(list(onehr, twohr, fourhr), type = "full")
  assign(induc.name, full.list)

  #get the inducible genes for each cluster
  combined.genes <- c(combined.genes, unique(get(induc.name)$gene_short_name)) #7640 genes
}

#count how many clusters a particular gene is inducible in
combined.genes.table <- as.data.frame(table(combined.genes))
#only get the genes which are inducible in 1 cluster
private.gene.names <- as.character(combined.genes.table[combined.genes.table$Freq==1,]$combined.genes) #1810 genes

#go back and figure out what cluster the privately inducible gene is
private.gene.clusters <- c()
for(g in private.gene.names) {
  cat(g, "\n")
  for(i in 0:14) {
    cluster <- paste('Cluster_',i,'_inducible', sep = '')
    if (g %in% get(cluster)$gene_short_name) {
      cat(cluster,"\n")
      private.gene.clusters <- c(private.gene.clusters, cluster)
    }  
  }
}

length(private.gene.clusters) #1810

#table of genes and what cluster they are inducible in
private.table <- as.data.frame(cbind(private.gene.names, private.gene.clusters), stringsAsFactors=F) #1810

# cluster1.private <- private.table[private.table$cluster=="Cluster_1_inducible",]$private
# cluster2.private <- private.table[private.table$cluster=="Cluster_2_inducible",]$private
# intersect(cluster1.private, cluster2.private) #character(0)

#now go back and figure out what time points the privately inducible gene was in
gene <- c()
cluster <- c()
for(g in private.gene.names) {
  cat(g, "\n")
  for(i in 0:14) {
    for (timept in c("1hr", "2hr", "4hr")) {
      cluster.timept <- paste('Cluster_',i,'_',timept, sep = '')
      if (g %in% get(cluster.timept)$gene_short_name) {
        cat(cluster.timept,"\n")
        gene <- c(gene, g)
        cluster <- c(cluster, cluster.timept)
      }
    }
  }
}

#table of privately induced genes and which cluster and timepoint they were inducible in
private.table.mat <- cbind(gene, cluster)
private.table <- as.data.frame(cbind(gene, cluster), stringsAsFactors=F) #1964 GOOD should be bigger number

#go back and extract out FC and associated values for the privately inducible gene in its cluster & timepoint
numbercell <- c()
percentexp <- c()
pvalue <- c()
qvalue <- c()
foldchange <- c()
for (i in 1:length(rownames(private.table))) {
  FC <- get(private.table[i,]$cluster)[private.table[i,]$gene,]$FC
  pval <- get(private.table[i,]$cluster)[private.table[i,]$gene,]$pval
  qval <- get(private.table[i,]$cluster)[private.table[i,]$gene,]$qval
  percell <- get(private.table[i,]$cluster)[private.table[i,]$gene,]$perc_cells_expressed
  numcell <- get(private.table[i,]$cluster)[private.table[i,]$gene,]$num_cells_expressed
  foldchange <- c(foldchange, FC)
  pvalue <- c(pvalue, pval)
  qvalue <- c(qvalue, qval)
  numbercell <- c(numbercell, numcell)
  percentexp <- c(percentexp, percell)

}

private.table <- as.data.frame(cbind(gene, cluster, foldchange, pvalue, qvalue, numbercell, percentexp), stringsAsFactors=F)
#ordered by fold change
private.table.ordered <- private.table[order(private.table$foldchange, decreasing = TRUE),]

write.table(private.table, file = "FILENAME_private-genes.csv", sep = ',', col.names = T)
