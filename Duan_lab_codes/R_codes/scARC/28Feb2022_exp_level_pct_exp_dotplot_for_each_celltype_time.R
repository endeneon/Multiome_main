# Chuxuan Li 28 Sept 2021
# look at the counts per cell

# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
library(readr)
library(ggplot2)
library(future)
library(viridis)
library(RColorBrewer)

features <- c("FOS", "NPAS4", "EGR1")

unique(filtered_obj$type.time.for.plot)
DefaultAssay(filtered_obj)

# pure_human_sub <- subset(pure_human, subset = spec.cell.type != "others")

DPnew(filtered_obj, 
      features = features,
      assay = "RNA",
      cols = c("red3", "white"),
      group.by = "type.time.for.plot") +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom"
        ) +
  coord_flip()
