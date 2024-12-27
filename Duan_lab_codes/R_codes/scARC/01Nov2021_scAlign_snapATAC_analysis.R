# Chuxuan Li 11/1/2021
# using scAlign and snapATAC to analyze the group 2, group 8 ATACseq data

library(scAlign)
library(Seurat)
library(SingleCellExperiment)
library(gridExtra)
library(PMA)
library(future)
set.seed(1886)
# set threads and parallelization
plan("multisession", workers = 3)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)

# read in separate raw files

# # Read files, separate the ATACseq data from the .h5 matrix
# g_2_0_raw <-
#   Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only//libraries_2_0/outs/filtered_feature_bc_matrix.h5")
# # load ATACseq information
# g_2_0_raw_pks <-
#   CreateSeuratObject(counts = g_2_0_raw$Peaks,
#                      project = "g_2_0_raw_pks")
# 
# g_2_1_raw <-
#   Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_2_1/outs/filtered_feature_bc_matrix.h5")
# # load ATACseq information
# g_2_1_raw_pks <-
#   CreateSeuratObject(counts = g_2_1_raw$Peaks,
#                      project = "g_2_1_raw_pks")
# 
# g_2_6_raw <-
#   Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_2_6/outs/filtered_feature_bc_matrix.h5")
# # load ATACseq information
# g_2_6_raw_pks <-
#   CreateSeuratObject(counts = g_2_6_raw$Peaks,
#                      project = "g_2_6_raw_pks")
# 
# 
# g_8_0_raw <-
#   Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_8_0/outs/filtered_feature_bc_matrix.h5")
# # load ATACseq information
# g_8_0_raw_pks <-
#   CreateSeuratObject(counts = g_8_0_raw$Peaks,
#                      project = "g_8_0_raw_pks")
# 
# g_8_1_raw <-
#   Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_8_1/outs/filtered_feature_bc_matrix.h5")
# # load ATACseq information
# g_8_1_raw_pks <-
#   CreateSeuratObject(counts = g_8_1_raw$Peaks,
#                      project = "g_8_1_raw_pks")
# 
# g_8_6_raw <-
#   Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_8_6/outs/filtered_feature_bc_matrix.h5")
# # load ATACseq information
# g_8_6_raw_pks <-
#   CreateSeuratObject(counts = g_8_6_raw$Peaks,
#                      project = "g_8_6_raw_pks")


# load the Signac-merged object
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/codes/aggregated_signac_object.RData")
# Then split the object
obj_list = SplitObject(aggr_signac, split.by = "group_time.ident")


## Independently preprocess each datasets
for (i in 1:length(obj_list)) {
  obj_list[[i]] <- NormalizeData(obj_list[[i]])
}
new_list <- obj_list
for (i in 1:length(obj_list)) {
  print(i)
  new_list[[i]] <- ScaleData(obj_list[[i]],
                             do.scale = T,
                             do.center = T,
                             verbose = T)
}

for (i in 1:length(new_list)) {
  new_list[[i]] <- FindVariableFeatures(new_list[[i]], 
                                        nfeatures = 30000)
}

# Extract common set of genes across all datasets
varfeat <- lapply(new_list, 
       function(seurat.obj){
         VariableFeatures(seurat.obj)})
sub_list <- list(obj_list[[1]], obj_list[[2]])
genes.use = Reduce(intersect, 
                   lapply(new_list, 
                          function(seurat.obj){
                            VariableFeatures(seurat.obj)
                            }))

# Convert each Seurat object into an SCE object, being sure to retain Seurat's metadata in SCE's colData field
for (i in 1:length(obj_list)) {
  
}
sce_list = lapply(new_list,
                  function(seurat.obj){
                    SingleCellExperiment(assays = list(counts = seurat.obj@assays$peaks@counts[genes.use, ],
                                                                          logcounts = seurat.obj@assays$peaks@data[genes.use, ],
                                                                          scale.data = seurat.obj@assays$peaks@scale.data[genes.use, ]),
                                                            colData = seurat.obj@meta.data)
                    })

# Create the combined scAlign object from list of SCE(s).
scAlignobj = scAlignCreateObject(sce.objects = sce_list,
                                 genes.use = genes.use,
                                 cca.reduce = T,
                                 ccs.compute = 20,
                                 project.name = "scAlign_ATACseq")

save("new_list", file = "list_of_processed_split_ATACseq_obj.RData")
# scAlignobj = scAlignCreateObject(sce.objects = sce_list,
#                                  genes.use = genes.use,
#                                  pca.reduce = T,
#                                  pcs.compute = 20,
#                                  project.name = "scAlign_ATACseq")

# Run scAlign with CCA results as input to the encoder (alignment).
scAlignobj_combined = scAlignMulti(scAlignobj,
                                   options=scAlignOptions(steps=15000,
                                                          log.every=5000,
                                                          batch.size=300,
                                                          perplexity=30,
                                                          norm=TRUE,
                                                          batch.norm.layer=FALSE,
                                                          architecture="large",  ## 3 layer neural network
                                                          num.dim=64),            ## Number of latent dimensions
                                   encoder.data="MultiCCA",
                                   decoder.data="scale.data",
                                   supervised='none',
                                   run.encoder=TRUE,
                                   run.decoder=TRUE,
                                   log.results=TRUE,
                                   log.dir=file.path('./tmp'),
                                   device="GPU")
