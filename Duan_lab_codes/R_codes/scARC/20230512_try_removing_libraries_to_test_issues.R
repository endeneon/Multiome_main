cleaned_lst <- cleaned_lst[11:length(cleaned_lst)]
features <- SelectIntegrationFeatures(object.list = cleaned_lst,
                                      nfeatures = 3000,
                                      fvf.nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = cleaned_lst,
                                  anchor.features = features,
                                  reference = c(1, 2, 3),
                                  reduction = "rpca",
                                  normalization.method = "SCT",
                                  scale = F, 
                                  dims = 1:50)
integrated <- IntegrateData(anchorset = anchors,
                            verbose = T)
