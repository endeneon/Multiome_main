load("final_lst_before_integration.RData")
# prepare reciprocal PCA integration ####
# find anchors
features <- SelectIntegrationFeatures(object.list = final_lst_for_integration,
                                      nfeatures = 3000,
                                      fvf.nfeatures = 3000)
final_lst_for_integration <- PrepSCTIntegration(final_lst_for_integration, anchor.features = features)
final_lst_for_integration <- lapply(X = final_lst_for_integration, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})

anchors <- FindIntegrationAnchors(object.list = final_lst_for_integration,
                                  anchor.features = features,
                                  reference = c(1, 2, 3),
                                  reduction = "rpca",
                                  normalization.method = "SCT",
                                  scale = F, 
                                  dims = 1:50)
#Found 9036 anchors
#all.genes <- transformed_lst[[1]]@assays$RNA@counts@Dimnames[[1]]

# integrate ####
integrated <- IntegrateData(anchorset = anchors,
                            verbose = T)
save(final_lst_for_integration, file = "final_lst_prepped_4_SCT_integration.RData")
save(integrated, file = "integrated_obj_nfeature_8000.RData")
rm(final_lst_for_integration)
