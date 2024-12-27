# Chuxuan Li 06/16/2022
# Plot dot plot for gene expression + gene activity for 18line

# init ####
load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/multiomic_obj.RData")
library(Seurat)
library(ggplot2)
library(Signac)
library(RColorBrewer)
library(future)
plan("multisession", workers = 2)
options(future.globals.maxSize = 2097374182400)

# calculate gene activity  ####
gene.activities <- GeneActivity(multiomic_obj)
multiomic_obj[['GACT']] <- CreateAssayObject(counts = gene.activities)
multiomic_obj <- NormalizeData(
  object = multiomic_obj,
  assay = 'GACT',
  normalization.method = 'LogNormalize'
)

sub <- unique(multiomic_obj$timextype.ident)
timextype_objs <- vector("list", length(types)*length(times))

for (i in 1:length(sub)) {
  timextype_objs[[i]] <- subset(multiomic_obj, timextype.ident == sub[i])
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
  c <- as.vector(unique(obj$cell.type))
  c <- str_replace_all(c, "_", " ")
  c <- str_replace(c, " pos", "+")
  c <- str_replace(c, " neg", "-")
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
  ggtitle("SP4 gene expression")


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
  ggtitle("SP4 gene activity score")
