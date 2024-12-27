# Chuxuan Li 03/15/2022
# Compare the results from DESeq2 and Monocle differential expression analysis of
#the 5-line RNAseq data

rm(list = ls())

library(readr)
library(stringr)

library(topGO)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Seurat)


# load pseudobulk data ####
setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_pseudobulk/DE_results/filtered_by_basemean/significant")
pseudo_path_list <- list.files(path = "./")
pseudo_df_list <- vector(mode = "list", length = length(pseudo_path_list))
for (i in 1:length(pseudo_path_list)){
  pseudo_df_list[[i]] <- read.csv(pseudo_path_list[i], 
                                  col.names = c("gene", "baseMean", "log2FC", "lfcSE", "stat", "pValue", "padj", "qvalue"))
}

# load monocle data ####
setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v5_monocle")
monocle_path_list <- list.files("./5line_output/Significant/", full.names = T, recursive = F)
monocle_df_list <- vector(mode = "list", length = length(monocle_path_list))
for (i in 1:length(monocle_path_list)){
  monocle_df_list[[i]] <- read.csv(monocle_path_list[i])
}
# vec <- as.character(ENSG_EntrezID_hs$`NCBI gene (formerly Entrezgene) ID`)[!is.na(ENSG_EntrezID_hs$`NCBI gene (formerly Entrezgene) ID`)] 
# test <- annFUN.org(whichOnto = "BP", feasibleGenes = vec, 
#                    mapping = "org.Hs.eg.db", ID = "entrez")
# vec <- ENSG_EntrezID_hs$`Gene stable ID`
# sum(is.na(vec))
# test <- annFUN.org(whichOnto = "BP", feasibleGenes = vec, 
#            mapping = "org.Hs.eg.db", ID = "ensembl") 
test <- annFUN.org(whichOnto = "BP", feasibleGenes = monocle_df_list[[2]]$gene_short_name, 
                   mapping = "org.Hs.eg.db", ID = "symbol") 

write.table(monocle_df_list[[4]]$gene_short_name, file = "GABA_6v0_up.txt", quote = F, sep = "\t", row.names = F,
            col.names = F)
# common DEGs ####


# GO term analysis ####
genelist <- pseudo_df_list[[4]]$pValue
names(genelist) <- pseudo_df_list[[4]]$gene
func <- function(x){return(genelist)}
fulllist <- res_0v6_list$GABA$pvalue
names(fulllist) <- res_0v6_list$GABA$gene
fulllist <- fulllist[!is.na(fulllist)]
testGO <- new("topGOdata",
              ontology = "BP", 
              allGenes = fulllist, 
              geneSel = func,
              nodeSize = 10, 
              annot = annFUN.org,
              #whichOnto = "BP", 
              #feasibleGenes = monocle_df_list[[1]]$gene_short_name, 
              mapping = "org.Hs.eg.db", ID = "symbol")
testres <- runTest(testGO, algorithm = "classic", statistic = "ks")
res <- GenTable(testGO,
                classicKS = testres,
                topNodes = 20)

