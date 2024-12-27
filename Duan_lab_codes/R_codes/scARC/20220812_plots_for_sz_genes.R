# Chuxuan Li 08/11/2022
# plot vlnplots, dimplot, dotplot for a set of genes from Jubao

genes <- c("NAB2", "GRIN2A", "IMMP2L", "SNAP91", 
           "PTPRK", "SP4", "RORB", "TCF4", "DIP2A")

# vlnplot for expression in 3 cell types, each has 3 time points ####
load("~/NVME/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/18line_pseudobulk_limma_post_normfactor_voom_dfs.RData")
lexprs <- as.data.frame(v_GABA$E[rownames(v_GABA$E) %in% genes, ])
lexprs$gene <- rownames(lexprs)
lexprs$cell.type <- "GABA"
lappend <- as.data.frame(v_nmglut$E[rownames(v_nmglut$E) %in% genes, ])
lappend$gene <- rownames(lappend)
lappend$cell.type <- "nmglut"
lexprs <- rbind(lexprs, lappend)
lappend <- as.data.frame(v_npglut$E[rownames(v_npglut$E) %in% genes, ])
lappend$gene <- rownames(lappend)
lappend$cell.type <- "npglut"
lexprs <- rbind(lexprs, lappend)

cnames <- (colnames(lexprs))[1:54]
times <- str_extract(cnames, "[0|1|6]hr")
lines <- str_extract(cnames, "CD_[0-9]+")

ldf <- matrix(nrow = 54*3*9, ncol = 5, dimnames = list(rep_len("", 54*3*9), 
                                                     c("Gene", "cell.line", 
                                                       "Time", "cell.type", 
                                                       "Expression")))
ldf <- as.data.frame(ldf)
ldf$Time <- rep(times, by = 3)
ldf$cell.type <- rep(c("GABA", "nmglut", "npglut"), each = 3)
ldf$cell.line <- rep(lines, each = 27)
ldf$Gene <- rep(rep(genes, each = 9), by = 18)
for (i in 1:nrow(ldf)) {
  g <- ldf$Gene[i]
  l <- ldf$cell.line[i]
  ti <- ldf$Time[i]
  ty <- ldf$cell.type[i]
  print(paste(g, l, ti, ty, sep = "-"))
  selectcol <- str_detect(colnames(lexprs), l) & str_detect(colnames(lexprs), ti)
  selectrow <- lexprs$gene == g & lexprs$cell.type == ty
  ldf$Expression[i] <- lexprs[selectrow, selectcol]
}

ldf$Expression <- as.numeric(ldf$Expression)
#ldf$Expression <- exp(ldf$Expression)

for (i in 1:length(genes)) {
  g = genes[i]
  print(g)
  subdf <- ldf[ldf$Gene == g, ] 
  filename <- paste0("./pseudobulk_DE/res_use_combatseq_mat/limma_voom_res/adj_cpm_9in18_res/limma_cnfved_exp_jubao_genes_vlnplot/",
                    g, "_limma_normalized_exp_vlnplot.pdf")
  pdf(file = filename)
  p <- ggplot(data = subdf, aes(x = Time, y = Expression)) + 
    geom_violin(aes(fill = Time), alpha = 0.5) +
    facet_grid(cols = vars(cell.type), scales = "free_y") +
    geom_dotplot(stackdir = "center", binaxis = "y", #binwidth = 0.03,
                 dotsize = 0.3) +
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    ggtitle(paste(g, "Average expression by cell type and time", sep = " - ")) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    NoLegend() 
  print(p)
  dev.off()
  filename <- paste0("./pseudobulk_DE/res_use_combatseq_mat/limma_voom_res/adj_cpm_9in18_res/limma_cnfved_exp_jubao_genes_vlnplot/",
                    g, "_limma_normalized_exp_vlnplot.png")
  png(filename = filename, width = 700, height = 550)
  print(p)
  dev.off()
}

# expression overlaying on dimplot ####
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000.RData")
library(Seurat)
integrated_labeled$cell.type[integrated_labeled$cell.type %in% c("GABA", "SEMA3E_pos_GABA", "SST_pos_GABA")] <- "GABA"
scobj <- subset(integrated_labeled, cell.type != "unknown")

unique(scobj$timexline.ident.forplot)


types <- unique(scobj$cell.type)
times <- unique(scobj$time.ident)
for (i in types) {
  for (j in times) {
    scobj$timextype.ident.forplot[scobj$cell.type == i & scobj$time.ident == j] <- paste(i, j, sep = ' ')
  }
}
unique(scobj$timextype.ident.forplot)
for (j in 1:length(types)) {
  subobj <- subset(scobj, cell.type == types[j])
  for (i in 1:length(genes)) {
    g = genes[i]
    print(g)
    filename <- paste0("./integrated_3000vargene_no_rat_gene_plots/featplots_jubao_genes/",
                       types[j], "_", g, "_featplot.pdf")
    pdf(file = filename, width = 27)
    p <- FeaturePlot(object = subobj, features = genes[i], slot = "data", cols = c("grey", "steelblue"),
                     split.by = "time.ident", keep.scale = "all") & 
      theme(legend.position = "right")
    
    print(p)
    dev.off()
  }
  
}

# Dotplot single cell expression by cell type and time point ####
library(dplyr)
library(patchwork)
library(ggplot2)
library(stringr)

unique(scobj$timextype.ident.forplot)
DotPlot(scobj, assay = "SCT", features = genes, cols = c("white", "red3"),
        group.by = "timextype.ident.forplot") + 
  coord_flip() +
  xlab("Genes") +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(), text = element_text(size = 10))

unique(scobj$timexline.ident.forplot)
genes <-
  c("BDNF",
    "VGF",
    "FOS",
    "NPAS4")

scobj$timexline.ident.forplot_v2 <-
  str_replace_all(string = scobj$timexline.ident.forplot,
                  pattern = '_',
                  replacement = '-')


DotPlot(scobj, 
        assay = "SCT", 
        features = genes, 
        cols = c("white", "red3"),
        group.by = "timexline.ident.forplot_v2") + 
  coord_flip() +
  xlab("Genes") +
  theme(axis.text.x = element_text(hjust = 0, 
                                   vjust = 1, angle = 315, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(), 
        text = element_text(size = 10),
        legend.position = "right")




# link peak to gene plots ####
load("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/multiomic_obj.RData")
load("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/type_obj_with_link_list.RData")
library(Signac)
idents_list <- unique(multiomic_obj$timextype.ident)
idents_list <- sort(idents_list)
idents_list
for (i in 1:length(obj_with_link_list)){
  for (j in genes) {
    print(j)
    filename <- paste0("../../Duan_project_024_ATAC/link_peak_to_gene_jubao_genes/", 
                       j, "_", idents_list[i], 
                     "_linkPeak_plot.pdf")
    print(filename)
    Idents(obj_with_link_list[[i]]) <- "timextype.ident"
    pdf(file = filename, height = 6, width = 12)
    p <- CoveragePlot(
      object = obj_with_link_list[[i]],
      region = j,
      features = j,
      expression.assay = "SCT", 
      extend.upstream = 100000,
      extend.downstream = 100000
    )
    print(p)
    dev.off()
  }
}



CoveragePlot(
  object = obj_with_link_list[[4]],
  region = "SNAP91",
  features = "SNAP91",
  expression.assay = "SCT", 
  extend.upstream = 10000, group.by = "time.ident", ymax = 500, 
  extend.downstream = 150000#, region.highlight = StringToGRanges("chr16-10199933-10200432")
)

CoveragePlot(
  object = type_obj_with_link_list[[1]],
  region = "SNAP91",
  features = "SNAP91",
  expression.assay = "SCT", ymax = 200,  
  extend.upstream = 100000, group.by = "time.ident", 
  extend.downstream = 120000, region.highlight = StringToGRanges("chr6-83656000-83658200")
)
