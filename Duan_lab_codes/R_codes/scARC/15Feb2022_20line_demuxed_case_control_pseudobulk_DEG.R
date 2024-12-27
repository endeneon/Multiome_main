# Chuxuan Li 02/15/2022
# using pseudobulk on demuxed_integrated_labeled_obj.RData, compute DEG for
#case vs. control, with samples being each of the 20 cell lines * 3 time points = 
#60 samples.

# init####
library(Seurat)
library(stringr)
library(readr)
library(ggplot2)
library(RColorBrewer)

load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/demux_integrated_labeled_obj.RData")


# make df ####
GABA <- subset(integrated_labeled, cell.type == "GABA")
glut <- subset(integrated_labeled, cell.type %in% c("NEFM_pos_glut"))

line_time_names <- rep_len("", length(unique(GABA$cell.line.ident)) * length(unique(GABA$time.ident)))
i = 0
for(line in unique(GABA$cell.line.ident)){
  for(time in unique(GABA$time.ident)){
    i = i + 1
    print(i)
    name <- paste(line, time, sep = "_")
    print(name)
    line_time_names[i] <- name 
    if (i == 1){
      print("first sample")
      subobj <- subset(GABA, subset = cell.line.ident == line &
                         time.ident == time)
      GABA_df <- as.array(rowSums(subobj@assays$RNA@counts))
    } else {
      print("the rest")
      subobj <- subset(GABA, subset = cell.line.ident == line &
                         time.ident == time)
      c_to_bind <- as.array(rowSums(subobj@assays$RNA@counts))
      GABA_df <- cbind(GABA_df, c_to_bind)
    }
  }
}

colnames(GABA_df) <- line_time_names

i = 0
for(line in unique(glut$cell.line.ident)){
  for(time in unique(glut$time.ident)){
    i = i + 1
    print(i)
    name <- paste(line, time, sep = "_")
    print(name)
    line_time_names[i] <- name 
    if (i == 1){
      print("first sample")
      subobj <- subset(glut, subset = cell.line.ident == line &
                         time.ident == time)
      glut_df <- as.array(rowSums(subobj@assays$RNA@counts))
    } else {
      print("the rest")
      subobj <- subset(glut, subset = cell.line.ident == line &
                         time.ident == time)
      c_to_bind <- as.array(rowSums(subobj@assays$RNA@counts))
      glut_df <- cbind(glut_df, c_to_bind)
    }
  }
}
colnames(glut_df) <- line_time_names

# edgeR ####

library(edgeR)

# make colData df
coldata_GABA <- data.frame(time = str_extract(string = colnames(GABA_pseudobulk),
                                                   pattern = "[0-6]hr"),
                           line = str_extract(string = colnames(GABA_pseudobulk),
                                              pattern = "CD_[0-9][0-9]"),
                           condition = str_extract(string = colnames(GABA_pseudobulk),
                                              pattern = "CD_[0-9][0-9]"))
rownames(coldata_GABA) <- colnames(GABA_pseudobulk)

# check
all(rownames(coldata_GABA) %in% colnames(GABA_pseudobulk)) # TRUE
all(rownames(coldata_GABA) == colnames(GABA_pseudobulk)) #TRUE


Glut_DGEList <- DGEList(counts = as.matrix(glut_df),
                        genes = rownames(glut_df), 
                        samples = colnames(glut_df),
                        remove.zeros = T)

Glut_DGEList <- calcNormFactors(Glut_DGEList)
Glut_DGEList <- estimateDisp(Glut_DGEList)

hist(cpm(Glut_DGEList), breaks = 10000, xlim = c(0, 10))
sum(cpm(Glut_DGEList) < 1)

# filter DGEList
Glut_DGEList <- Glut_DGEList[rowSums(cpm(Glut_DGEList) > 1) >= 10, ,
                             keep.lib.sizes = F]
# Glut_DGEList$genes

Glut_design_matrix <- coldata_glut
Glut_design_matrix$condition <- as.factor(Glut_design_matrix$condition)
Glut_design_matrix$type <- as.factor(Glut_design_matrix$type)

Glut_design <- model.matrix(~ 0 + 
                              Glut_design_matrix$condition +
                              Glut_design_matrix$type)

Glut_fit <- glmQLFit(Glut_DGEList,
                     design = Glut_design)
Glut_1vs0 <- glmQLFTest(Glut_fit, 
                        contrast = c(-1, 1, 0, 0, 0, 0, 0))
Glut_6vs0 <- glmQLFTest(Glut_fit, 
                        contrast = c(-1, 0, 1, 0, 0, 0, 0))
Glut_6vs1 <- glmQLFTest(Glut_fit, 
                        contrast = c(0, -1, 1, 0, 0, 0, 0))

Glut_1vs0 <- Glut_1vs0$table
Glut_6vs0 <- Glut_6vs0$table
Glut_6vs1 <- Glut_6vs1$table


Glut_1vs0[rownames(Glut_1vs0) %in% genelist, ]
Glut_6vs0[rownames(Glut_6vs0) %in% genelist, ]
Glut_6vs1[rownames(Glut_6vs1) %in% genelist, ]



GABA_DGEList <- DGEList(counts = as.matrix(GABA_df),
                        genes = rownames(GABA_df), 
                        samples = colnames(GABA_df),
                        remove.zeros = T)

GABA_DGEList <- calcNormFactors(GABA_DGEList)
GABA_DGEList <- estimateDisp(GABA_DGEList)

hist(cpm(GABA_DGEList), breaks = 10000, xlim = c(0, 10))
sum(cpm(GABA_DGEList) < 1)

# filter DGEList
GABA_DGEList <- GABA_DGEList[rowSums(cpm(GABA_DGEList) > 1) >= 5, ,
                             keep.lib.sizes = F]
# GABA_DGEList$genes

GABA_design_matrix <- coldata_GABA
GABA_design_matrix$condition <- as.factor(GABA_design_matrix$condition)
GABA_design_matrix$type <- as.factor(GABA_design_matrix$type)

GABA_design <- model.matrix(~ 0 + 
                              GABA_design_matrix$condition +
                              GABA_design_matrix$type)

GABA_fit <- glmQLFit(GABA_DGEList,
                     design = GABA_design)
GABA_1vs0 <- glmQLFTest(GABA_fit, 
                        contrast = c(-1, 1, 0, 0, 0, 0, 0))
GABA_6vs0 <- glmQLFTest(GABA_fit, 
                        contrast = c(-1, 0, 1, 0, 0, 0, 0))
GABA_6vs1 <- glmQLFTest(GABA_fit, 
                        contrast = c(0, -1, 1, 0, 0, 0, 0))

GABA_1vs0 <- GABA_1vs0$table
GABA_6vs0 <- GABA_6vs0$table
GABA_6vs1 <- GABA_6vs1$table


GABA_1vs0[rownames(GABA_1vs0) %in% genelist, ]
GABA_6vs0[rownames(GABA_6vs0) %in% genelist, ]
GABA_6vs1[rownames(GABA_6vs1) %in% genelist, ]


# bar graph ####
genelist <- genelist[genelist != c("GAPDH")]
df_to_plot <- rbind(res_glut_0v1[genelist, ],
                    res_glut_0v6[genelist, ],
                    res_GABA_0v1[genelist, ],
                    res_GABA_0v6[genelist, ])
df_to_plot$time <- c(rep_len("0v1", 4), rep_len("0v6", 4), rep_len("0v1", 4), rep_len("0v6", 4))                    
df_to_plot$cell.type <- c(rep_len("glut", 8), rep_len("GABA", 8))
df_to_plot$gene.name <- c(genelist, genelist)
df_to_plot <- as.data.frame(df_to_plot)

ggplot(df_to_plot, aes(x = cell.type, 
                       y = log2FoldChange,
                       color = cell.type,
                       fill = time,
                       group = time,
                       ymax = log2FoldChange - 1/2*lfcSE, 
                       ymin = log2FoldChange + 1/2*lfcSE
)) + 
  xlab("") +
  ylab("log2(FC)") +
  geom_col(position = position_dodge(0.6),
           width = 0.5,
           color = "black") +
  geom_errorbar(color = "black",
                position = position_dodge(0.6),
                width = 0.5) +
  facet_grid(cols = vars(gene.name)) +
  scale_fill_manual(values = c("indianred", "steelblue")) +
  theme_bw() +
  theme(legend.title = element_blank())
