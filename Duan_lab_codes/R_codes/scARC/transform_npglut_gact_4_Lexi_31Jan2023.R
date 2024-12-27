# Siwei 31 Jan 2024
# Normalize npglut GeneActivity assay for Alexi's DEG analysis

# init ####
library(Seurat)
library(Signac)

# library(harmony)

library(future)

library(Matrix)

library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)

library(stringr)
library(ggplot2)
# library(scales)
library(RColorBrewer)
# library(plyranges)
# library(gplots)
# library(grDevices)
# library(viridis)
# library(colorspace)

plan("multisession", workers = 4)
set.seed(42)
options(future.globals.maxSize = 229496729600)
options(Seurat.object.assay.version = "v5")

# temp function

# load data ####
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/npglut_geneActivity_list_all_times.RData")


####
# reinstall Matrix 1.6.1?
for (i in 1:length(npglut_geneActivity)) {
  print(i)
  temp_obj <- npglut_geneActivity[[i]]
  temp_obj[['RNA']] <- NULL
  temp_obj <-
    UpdateSeuratObject(temp_obj)
  temp_obj[["gact5"]] <-
    as(object = temp_obj[["gact5"]],
       Class = "Assay5")
  
  temp_obj[["gact5"]] <-
    split(temp_obj[["gact5"]],
          f = temp_obj$group.ident)
  
    # split(npglut_geneActivity[[i]],
    #       f = npglut_geneActivity[[i]]$group.ident)
  temp_obj <-
    SCTransform(temp_obj)
  temp_obj <-
    RunPCA(temp_obj)
  temp_obj <-
    RunUMAP(temp_obj,
            dims = 1:30)
  
  temp_obj <-
    IntegrateLayers(temp_obj,
                    method = CCAIntegration,
                    normalization.method = "SCT",
                    verbose = T)
  temp_obj <-
    FindNeighbors(temp_obj,
                  reduction = "integrated.dr",
                  dims = 1:30)
  temp_obj <-
    FindClusters(temp_obj,
                 resolution = 0.6)
  temp_obj <-
    RunUMAP(temp_obj,
            dims = 1:30,
            reduction = "integrated.dr")
  npglut_geneActivity[[i]] <-
    temp_obj
}

temp_obj@assays$gact5@layers$counts@Dimnames[[1]] <- rownames(temp_obj)
temp_obj@assays$gact5@layers$counts@Dimnames[[2]] <- colnames(temp_obj)
rownames(temp_obj@assays$gact5@layers$counts)

temp_obj@assays$gact5@layers$data@Dimnames[[1]] <- rownames(temp_obj)
temp_obj@assays$gact5@layers$data@Dimnames[[2]] <- colnames(temp_obj)

rownames(temp_obj@assays$gact5@layers$data)

rownames(temp_obj@assays$gact5@cells)
rownames(temp_obj@assays$gact5@features)
