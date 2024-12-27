# Chuxuan Li 03/03/2022
# create gene activity matrix using 5000 variable genes from RNAseq analysis v3,
#ATACseq object from ATAcseq analysis v2 (new_peak_set_after_motif.RData),
#then plot gene activity and gene expression for 5000 genes by time

library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)

library(gplots)
set.seed(1234)
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/gene_activity_and_expression_matrices.RData")
# get Seurat object with only peaks assay

zero_ATAC <- subset(obj_complete, time.ident == "0hr")
one_ATAC <- subset(obj_complete, time.ident == "1hr")
six_ATAC <- subset(obj_complete, time.ident == "6hr")
zero_ATAC <- zero_ATAC@assays$peaks
zero_ATAC <- CreateSeuratObject(
  counts = zero_ATAC,
  assay = "peaks"
)
one_ATAC <- one_ATAC@assays$peaks
one_ATAC <- CreateSeuratObject(
  counts = one_ATAC,
  assay = "peaks"
)
six_ATAC <- six_ATAC@assays$peaks
six_ATAC <- CreateSeuratObject(
  counts = six_ATAC,
  assay = "peaks"
)
zero_gact <- GeneActivity(zero_ATAC, 
                          features = human_only@assays$integrated@data@Dimnames[[1]],
                          extend.upstream = 1000,
                          extend.downstream = 0)
one_gact <- GeneActivity(one_ATAC, 
                          features = human_only@assays$integrated@data@Dimnames[[1]],
                          extend.upstream = 1000,
                          extend.downstream = 0)
six_gact <- GeneActivity(six_ATAC, 
                          features = human_only@assays$integrated@data@Dimnames[[1]],
                          extend.upstream = 1000,
                          extend.downstream = 0)

# gene expression matrix
zero_gexp <- subset(human_only, time.ident == "0hr")@assays$integrated@data
one_gexp <- subset(human_only, time.ident == "1hr")@assays$integrated@data
six_gexp <- subset(human_only, time.ident == "6hr")@assays$integrated@data


# plot heatmap
# get alll gene expression data
library(gplots)
library(colorspace)
gene_activity_matrix <- 
  as.matrix(t(data.frame(rbind(rowSums(zero_gact),
                               rowSums(one_gact),
                               rowSums(six_gact)),
                         stringsAsFactors = F)))
rownames(gene_activity_matrix)

heatmap.2(log1p(gene_activity_matrix),
          Rowv = T,
          Colv = F,
          dendrogram = "row",
          col = diverge_hcl(n = 12, 
                            h = c(255, 330),
                            l = c(40, 90)),
          trace = "none",
          labRow = "",
          labCol = "",
          scale = "row")

gene_expression_matrix <- 
  as.matrix(t(data.frame(rbind(rowSums(zero_gexp),
                               rowSums(one_gexp),
                               rowSums(six_gexp)),
                         stringsAsFactors = F)))

heatmap.2(log1p(gene_expression_matrix + 147),
          Rowv = T,
          Colv = F,
          dendrogram = "none",
          col = diverge_hcl(n = 12, 
                            h = c(130, 43),
                            l = c(70, 90),
                            c = 100),
          trace = "none",
          labRow = "",
          labCol = "",
          scale = "row",
          na.rm = T)
sum(gene_expression_matrix < 0)
min(gene_expression_matrix)

distCor <- function(x) as.dist(1-cor(x))
zClust <- function(x, scale= "row" , zlim=c(-3,3), method="average") {
  if (scale == "row") z <- t(scale(t(x)))
  # if (scale == "row") z <- x
  if (scale == "col") z <- scale(x)
  z <- pmin(pmax(z, zlim[1]), zlim[2])
  hcl_row <- hclust(distCor(t(z)), method = method)
  hcl_col <- hclust(distCor(z), method = method)
  return(list(data = z, Rowv = as.dendrogram(hcl_row), Colv = as.dendrogram(hcl_col)))
}



heatmap.2(log1p(z$data),
          Rowv = z$Rowv,
          Colv = F,
          dendrogram = "none",
          col = diverge_hcl(n = 12, 
                            h = c(130, 43),
                            l = c(70, 90),
                            c = 100),
          trace = "none",
          labRow = "",
          labCol = "",
          scale = "row",
          na.rm = T)

library(ggplot2)


sum(is.na(gene_activity_matrix))
sum(is.infinite(gene_activity_matrix))
sum(is.nan(gene_activity_matrix))
# gene activity matrix
z <- zClust(log1p(gene_activity_matrix))

heatmap.2(z$data,
          Rowv = z$Rowv,
          Colv = F,
          dendrogram = "none",
          col = diverge_hcl(n = 12, 
                            h = c(255, 330),
                            l = c(40, 90)),
          trace = "none",
          labRow = "",
          labCol = c("0hr", "1hr", "6hr"),
          scale = "row",
          main = "Gene Activity Matrix") 
  ggtitle("Gene Activity Matrix")

  
  z <- zClust(log1p(gene_expression_matrix + 147))
  heatmap.2(z$data,
            Rowv = z$Rowv,
            Colv = F,
            dendrogram = "none",
            col = diverge_hcl(n = 12, 
                              h = c(130, 43),
                              l = c(70, 90),
                              c = 100),
            trace = "none",
            labRow = "",
            labCol = c("0hr", "1hr", "6hr"),
            scale = "row",
            main = "Gene Expression Matrix") 
