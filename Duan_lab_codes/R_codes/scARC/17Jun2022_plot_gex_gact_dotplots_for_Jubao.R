# Chuxuan Li 06/17/2022
# plot gene expression and gene activity score dotplots with 5-line data

# init ####
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(future)
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes/ATACseq_5line_objonly_labeled.RData")

plan("multisession", workers = 2)
options(future.globals.maxSize = 2097374182400)

# calculate gene activity  ####
gene.activities <- GeneActivity(obj_complete)
obj_complete[['GACT']] <- CreateAssayObject(counts = gene.activities)
obj_complete <- NormalizeData(
  object = obj_complete,
  assay = 'GACT',
  normalization.method = 'LogNormalize'
)
unique(obj_complete$fine.cell.type)
obj_complete <- subset(obj_complete, fine.cell.type %in% c("NEFM+/CUX2- glut",
                                                           "NEFM-/CUX2+ glut",
                                                           "GABA"))
sub <- unique(obj_complete$celltype.time.ident)
sub <- sub[!is.na(sub)]
timextype_objs <- vector("list", length(sub))

for (i in 1:length(sub)) {
  timextype_objs[[i]] <- subset(obj_complete, celltype.time.ident == sub[i])
}

# generate plot df ####
df_to_plot <- data.frame(time = rep_len("", length(timextype_objs)),
                         cell.type = rep_len("", length(timextype_objs)),
                         gex = rep_len(0, length(timextype_objs)),
                         pct.exp = rep_len(0, length(timextype_objs)),
                         gact = rep_len(0, length(timextype_objs)), 
                         pct.gact = rep_len(0, length(timextype_objs)), stringsAsFactors = F)

for (i in 1:length(timextype_objs)) {
  obj <- timextype_objs[[i]]
  t <- unique(obj$time.ident)
  c <- as.vector(unique(obj$fine.cell.type))
  if (c == "NEFM+/CUX2- glut") {
    c <- "NEFM+ glut"
  } else if (c == "NEFM-/CUX2+ glut") {
    c <- "NEFM- glut"
  }
  print(paste(t, c))
  df_to_plot$time[i] <- t
  df_to_plot$cell.type[i] <- c
  DefaultAssay(obj) <- "SCT"
  print(ncol(obj))
  df_to_plot$gex[i] <- sum(expm1(obj@assays$SCT@data[rownames(obj) == "SP4", ])) /
    ncol(obj)
  df_to_plot$pct.exp[i] <- sum(obj@assays$SCT@data[rownames(obj) == "SP4", ] != 0) /
    ncol(obj)
  DefaultAssay(obj) <- "GACT"
  df_to_plot$gact[i] <- sum(expm1(obj@assays$GACT@data[rownames(obj) == "SP4", ])) /
    ncol(obj)
  df_to_plot$pct.gact[i] <- sum(obj@assays$GACT@data[rownames(obj) == "SP4", ] != 0) /
    ncol(obj)
}

# plot ####

ggplot(df_to_plot,
       aes(x = as.factor(time),
           y = as.factor(cell.type),
           size = pct.exp,
           fill = gex)) +
  geom_point(shape = 21) +
  scale_size(range = c(5, 9)) +
  ylab("cell type") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "lightsteelblue2", "steelblue4")) +
  theme_bw() +
  labs(size = "pct. cell expressed", fill = "avg. gene expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        text = element_text(size = 10)) +
  ggtitle("FOS gene expression - 5 line")


ggplot(df_to_plot,
       aes(x = as.factor(time),
           y = as.factor(cell.type),
           size = pct.gact,
           fill = gact)) +
  geom_point(shape = 21) +
  scale_size(range = c(5, 9)) +
  ylab("cell type") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "lightsalmon", "salmon4")) +
  theme_bw() +
  labs(size = "pct. cell with\nnonzero score", fill = "avg. gene activity score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        text = element_text(size = 10)) +
  ggtitle("FOS gene activity score - 5 line")

# check ####
hist(as.vector(timextype_objs[[7]]@assays$SCT@data[rownames(obj) == "SP4", ]),
     breaks = 10, main = "GABA 0hr", xlab = "count")

obj_complete$cell.line.ident <- NA
obj_complete$cell.line.ident[str_sub(obj_complete@assays$RNA@counts@Dimnames[[2]], end = -3L) %in% 
                                            c(g_2_0_bc_lst[[1]], g_2_1_bc_lst[[1]], g_2_6_bc_lst[[1]])] <- "CD_27"
obj_complete$cell.line.ident[str_sub(obj_complete@assays$RNA@counts@Dimnames[[2]], end = -3L) %in% 
                                            c(g_2_0_bc_lst[[2]], g_2_1_bc_lst[[2]], g_2_6_bc_lst[[2]])] <- "CD_54"

obj_complete$cell.line.ident[str_sub(obj_complete@assays$RNA@counts@Dimnames[[2]], end = -3L) %in% 
                                            c(g_8_0_bc_lst[[1]], g_8_1_bc_lst[[1]], g_8_6_bc_lst[[1]])] <- "CD_08"
obj_complete$cell.line.ident[str_sub(obj_complete@assays$RNA@counts@Dimnames[[2]], end = -3L) %in% 
                                            c(g_8_0_bc_lst[[2]], g_8_1_bc_lst[[2]], g_8_6_bc_lst[[2]])] <- "CD_25"
obj_complete$cell.line.ident[str_sub(obj_complete@assays$RNA@counts@Dimnames[[2]], end = -3L) %in% 
                                            c(g_8_0_bc_lst[[3]], g_8_1_bc_lst[[3]], g_8_6_bc_lst[[3]])] <- "CD_26"
sum(is.na(obj_complete$cell.line.ident))
unique(obj_complete$cell.line.ident)
save(obj_complete, file = "5line_multiomic_obj_complete_added_geneactivities.RData")
g2 <- subset(obj_complete, cell.line.ident %in% c("CD_54", "CD_27"))
g8 <- subset(obj_complete, cell.line.ident %in% c("CD_26", "CD_25", "CD_08"))
sub <- unique(g8$celltype.time.ident)
sub <- sub[!is.na(sub)]
timextype_objs <- vector("list", length(sub))

for (i in 1:length(sub)) {
  timextype_objs[[i]] <- subset(g8, celltype.time.ident == sub[i])
}

# generate plot df ####
df_to_plot <- data.frame(time = rep_len("", length(timextype_objs)),
                         cell.type = rep_len("", length(timextype_objs)),
                         gex = rep_len(0, length(timextype_objs)),
                         pct.exp = rep_len(0, length(timextype_objs)),
                         gact = rep_len(0, length(timextype_objs)), 
                         pct.gact = rep_len(0, length(timextype_objs)), stringsAsFactors = F)

for (i in 1:length(timextype_objs)) {
  obj <- timextype_objs[[i]]
  t <- unique(obj$time.ident)
  c <- as.vector(unique(obj$fine.cell.type))
  if (c == "NEFM+/CUX2- glut") {
    c <- "npglut"
  } else if (c == "NEFM-/CUX2+ glut") {
    c <- "nmglut"
  }
  print(paste(t, c))
  df_to_plot$time[i] <- t
  df_to_plot$cell.type[i] <- c
  DefaultAssay(obj) <- "SCT"
  print(ncol(obj))
  df_to_plot$gex[i] <- sum(obj@assays$SCT@data[rownames(obj) == "SP4", ]) /
    ncol(obj)
  df_to_plot$pct.exp[i] <- sum(obj@assays$SCT@data[rownames(obj) == "SP4", ] != 0) /
    ncol(obj)
  DefaultAssay(obj) <- "GACT"
  df_to_plot$gact[i] <- sum(obj@assays$GACT@data[rownames(obj) == "SP4", ]) /
    ncol(obj)
  df_to_plot$pct.gact[i] <- sum(obj@assays$GACT@data[rownames(obj) == "SP4", ] != 0) /
    ncol(obj)
}

# plot ####

ggplot(df_to_plot,
       aes(x = as.factor(time),
           y = as.factor(cell.type),
           size = pct.exp,
           fill = gex)) +
  geom_point(shape = 21) +
  scale_size(range = c(5, 9)) +
  ylab("cell type") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "lightsteelblue2", "steelblue4")) +
  theme_bw() +
  labs(size = "pct. cell expressed", fill = "avg. gene expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        text = element_text(size = 10)) +
  ggtitle("SP4 gene expression - 5 line - group8")


ggplot(df_to_plot,
       aes(x = as.factor(time),
           y = as.factor(cell.type),
           size = pct.gact,
           fill = gact)) +
  geom_point(shape = 21) +
  scale_size(range = c(5, 9)) +
  ylab("cell type") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "lightsalmon", "salmon4")) +
  theme_bw() +
  labs(size = "pct. cell with\nnonzero score", fill = "avg. gene activity score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        text = element_text(size = 10)) +
  ggtitle("SP4 gene activity score - 5 line - group8")
