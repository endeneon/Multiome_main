# Chuxuan Li 02/21/2022
# 1) check number of DEG (per cell type) in each of the 5 cell types in the 20-line
#dataset obtained from different min.pct cutoffs and plot the numbers
# 2) check the effect on DEGs when expression in one cell population is 0 using
#library 39 whose 0hr has no expression of NPAS4


# init ####
library(ggplot2)
library(Seurat)
library(stringr)
library(ggrepel)
library(RColorBrewer)

# load data 
#load("DEG_pct0-0.15.RData")
#load("DEG_5line_0.03.RData)
#load("DEG_5line_0.02.RData")
#load("DEG_5line_0.01.RData")
#load("DEG_5line_0.04-6.RData")
#load("DEG_5line_0.1.RData")
#load("./5line_codes_and_rdata/renormalized_DEG_pct_0.05.RData")
#deg_list <- list(dfs_pct0, dfs_pct0.05, dfs_pct0.1, dfs_pct0.15)
#deg_list <- list(DEG_5line_0.01, DEG_5line_0.02, DEG_5line_0.03)
#deg_list <- list(dfs_pct0.01, dfs_pct0.02, dfs_pct0.03)
#deg_list <- list(DEG_5line_0.04, DEG_5line_0.05, DEG_5line_0.06)
deg_list <- list(DEG_5line_0.06)

# evaluate min.pct values ####
# initialize dfs to store gene counts
num_types <- length(unique(deg_list[[1]][[1]]$cell_type))
dfs_to_plot <- vector(mode = "list", length = length(deg_list))
col_length <- num_types * 4

for (i in 1:length(dfs_to_plot)){
  dfs_to_plot[[i]] <- data.frame(cell.type = rep_len(NA, col_length),
                                 times = rep_len(NA, col_length), # 2 time groups
                                 n.genes = rep_len(0, col_length),
                                 filtered.status = rep_len(NA, col_length) # either total or filtered counts
  )
}
# separate 1v0, 6v0
times <- c("1v0", "6v0")
fstatus <- c("q-value > 0.05", "q-value < 0.05")
i = 0
for (df in deg_list){
  i = i + 1
  for (j in 1:2){
    subdf <- df[[j]]
    types <- unique(subdf$cell_type)
    for (k in 1:length(types)){
      type <- types[k] 
      print(type)
      for (l in 1:2){
        print(j * 8 + k * 2 + l - 10)
        dfs_to_plot[[i]]$cell.type[j * 8 + k * 2 + l - 10] <- type
        dfs_to_plot[[i]]$times[j * 8 + k * 2 + l - 10] <- times[j]
        if (l == 1){
          dfs_to_plot[[i]]$n.genes[j * 8 + k * 2 + l - 10] <- sum(subdf$cell_type == type &
                                                                     subdf$p_val_adj > 0.05)
        } else {
          dfs_to_plot[[i]]$n.genes[j * 8 + k * 2 + l - 10] <- sum(subdf$cell_type == type &
                                                                     subdf$p_val_adj < 0.05)
        }
        dfs_to_plot[[i]]$filtered.status[j * 8 + k * 2 + l - 10] <- fstatus[l]
      }
    }
  }
}

print(dfs_to_plot[[1]]$cell.type)
print(dfs_to_plot[[1]]$n.genes[dfs_to_plot[[1]]$filtered.status == "q-value < 0.05"])
print(dfs_to_plot[[2]]$n.genes[dfs_to_plot[[2]]$filtered.status == "q-value < 0.05"])
# print(dfs_to_plot[[3]]$n.genes[dfs_to_plot[[3]]$filtered.status == "q-value < 0.05"])

titles <- c("min.pct = 0.06")
filenames <- gsub("=", "_", gsub("\\.", "_", gsub(" ", "", titles)))

i = 0
for (df in dfs_to_plot){
  i = i + 1
  print(titles[i])
  jpeg(paste0("../5line_plots/renormalized_plots/DEG_plots/", filenames[i], "_DEG_counts.jpeg"))
  p <- ggplot(df, aes(x = cell.type,
                      y = n.genes,
                      fill = filtered.status)) +
    geom_col(position = "stack") +
    facet_grid(cols = vars(times)) + 
    scale_fill_manual(values = brewer.pal(3, "Paired")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) + 
    xlab("Cell Type") +
    ylab("Number of DEGs") +
    labs(fill = element_blank()) +
    ggtitle(titles[i])
  print(p)
  dev.off()
}

# test using log2FC cutoff too ####

# test GO term enrichment on the subset with log2FC > 0.5 and min.pct > 0.05
DEG_test_list_0v1_pos <- DEG_5line_0.05[[1]][DEG_5line_0.05[[1]]$avg_diff > 0.5 & 
  DEG_5line_0.05[[1]]$p_val_adj < 0.05, ]
DEG_test_list_0v1_neg <- DEG_5line_0.05[[1]][DEG_5line_0.05[[1]]$avg_diff < -0.5 & 
                                               DEG_5line_0.05[[1]]$p_val_adj < 0.05, ]

DEG_test_list_0v6_pos <- DEG_5line_0.05[[2]][DEG_5line_0.05[[2]]$avg_diff > 0.5 & 
                                           DEG_5line_0.05[[2]]$p_val_adj < 0.05, ]
DEG_test_list_0v6_neg <- DEG_5line_0.05[[2]][DEG_5line_0.05[[2]]$avg_diff< -0.5 & 
                                               DEG_5line_0.05[[2]]$p_val_adj < 0.05, ]


sum(DEG_5line_0.01[[2]]$avg_diff > 0)
sum(DEG_5line_0.01[[2]]$avg_diff < 0)

# try log2FC threshold 0.1
DEG_test_list_0v1_pos <- DEG_5line_0.05[[1]][DEG_5line_0.05[[1]]$avg_diff > 0.1 & 
                                               DEG_5line_0.05[[1]]$p_val_adj < 0.05, ]
DEG_test_list_0v1_neg <- DEG_5line_0.05[[1]][DEG_5line_0.05[[1]]$avg_diff < -0.1 & 
                                               DEG_5line_0.05[[1]]$p_val_adj < 0.05, ]

DEG_test_list_0v6_pos <- DEG_5line_0.05[[2]][DEG_5line_0.05[[2]]$avg_diff > 0.1 & 
                                               DEG_5line_0.05[[2]]$p_val_adj < 0.05, ]
DEG_test_list_0v6_neg <- DEG_5line_0.05[[2]][DEG_5line_0.05[[2]]$avg_diff< -0.1 & 
                                               DEG_5line_0.05[[2]]$p_val_adj < 0.05, ]

# no log2FC threshold
DEG_test_list_0v1_pos <- DEG_5line_0.05[[1]][DEG_5line_0.05[[1]]$avg_diff > 0 & 
                                               DEG_5line_0.05[[1]]$p_val_adj < 0.05, ]
DEG_test_list_0v1_neg <- DEG_5line_0.05[[1]][DEG_5line_0.05[[1]]$avg_diff < 0 & 
                                               DEG_5line_0.05[[1]]$p_val_adj < 0.05, ]

DEG_test_list_0v6_pos <- DEG_5line_0.05[[2]][DEG_5line_0.05[[2]]$avg_diff > 0 & 
                                               DEG_5line_0.05[[2]]$p_val_adj < 0.05, ]
DEG_test_list_0v6_neg <- DEG_5line_0.05[[2]][DEG_5line_0.05[[2]]$avg_diff< 0 & 
                                               DEG_5line_0.05[[2]]$p_val_adj < 0.05, ]

# separate cell type
types <- c("NEFM_pos_glut", "NEFM_neg_glut", "GABA", "NPC")
p1 <- vector(mode = "list", length = length(types))
p6 <- vector(mode = "list", length = length(types))
n1 <- vector(mode = "list", length = length(types))
n6 <- vector(mode = "list", length = length(types))
for (i in 1:length(types)){
  p1[[i]] <- DEG_test_list_0v1_pos[DEG_test_list_0v1_pos$cell_type == types[i], ]
  n1[[i]] <- DEG_test_list_0v1_neg[DEG_test_list_0v1_neg$cell_type == types[i], ]
  p6[[i]] <- DEG_test_list_0v6_pos[DEG_test_list_0v6_pos$cell_type == types[i], ]
  n6[[i]] <- DEG_test_list_0v6_neg[DEG_test_list_0v6_neg$cell_type == types[i], ]
  if (nrow(p1[[i]]) > 2000){
    p1[[i]] <- p1[[i]][1:2000, ]
  }
  if (nrow(n1[[i]]) > 2000){
    n1[[i]] <- n1[[i]][1:2000, ]
  }
  if (nrow(p6[[i]]) > 2000){
    p6[[i]] <- p6[[i]][1:2000, ]
  }
  if (nrow(n6[[i]]) > 2000){
    n6[[i]] <- n6[[i]][1:2000, ]
  }
  
  write.table(p1[[i]]$gene_symbol, 
              file = paste0(types[i], "_1v0_pos.txt"), 
              quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(n1[[i]]$gene_symbol, 
              file = paste0(types[i], "_1v0_neg.txt"), 
              quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(p6[[i]]$gene_symbol, 
              file = paste0(types[i], "_6v0_pos.txt"), 
              quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(n6[[i]]$gene_symbol, 
              file = paste0(types[i], "_6v0_neg.txt"), 
              quote = F, sep = "\t", row.names = F, col.names = F)
}



# GO term analysis ####
txtlist <- list.files(".", pattern = "*.txt")
plotlist <- vector(mode = "list", length = 2*length(txtlist))
for(i in 1:length(txtlist)){
  df <- read.delim(txtlist[[i]])
  col_names <- c("term", "REF", "#", "expected", "over_under", 
                 "fold_enrichment", "raw_p_value", "FDR")
  colnames(df) <- col_names
  df$term <- gsub(pattern = '\\(GO:[0-9]+\\)', replacement = "", x = df$term)
  pos <- df[df$over_under == "+", ]
  neg <- df[df$over_under == "-", ]
  pos <- pos[1:10, ]
  pos$FDR4label <- paste0("FDR = ", pos$FDR)
  neg <- neg[1:10, ]
  neg$FDR4label <- paste0("FDR = ", neg$FDR)
  
  title <- str_extract(txtlist[[i]], "[a-z]+_[1,6]v0_[a-z]+")
  title_pos <- paste0(title, "_pos")
  title_neg <- paste0(title, "_neg")
  print(title)

  plotlist[[2*i - 1]] <- ggplot(pos, 
              aes(x = term,
                  y = as.numeric(fold_enrichment))) +
    geom_col(fill = "skyblue") + 
    geom_text(aes(label = FDR4label), 
              position = position_dodge(width = 0.9), 
              hjust = 1,
              vjust = 0.5) + 
    xlab("GO term") +
    ylab("Fold Enrichment") + 
    ggtitle(title_pos) + 
    coord_flip()
  plotlist[[2*i]] <- ggplot(neg, 
                  aes(x = term,
                      y = as.numeric(fold_enrichment))) +
    geom_col(fill = "skyblue") + 
    geom_text(aes(label = FDR4label), 
              position = position_dodge(width = 0.9), 
              hjust = 1,
              vjust = 0.5) + 
    xlab("GO term") +
    ylab("Fold Enrichment") + 
    ggtitle(title_neg) + 
    coord_flip()
}

plotlist[[16]]$data$term[[5]] <- 
  "detection of chemical stimulus\ninvolved in sensory perception of smell "
plotlist[[16]]


# plot expression of marker genes ####
df_1v0 <- DEG_5line_0.05[[1]][DEG_5line_0.05[[1]]$gene_symbol %in% c("FOS", "EGR1", "BDNF", "VGF"), ]
df_6v0 <- DEG_5line_0.05[[2]][DEG_5line_0.05[[2]]$gene_symbol %in% c("FOS", "EGR1", "BDNF", "VGF"), ]

df_to_plot <- rbind(df_1v0,
                    df_6v0)
df_to_plot$source <- c(rep_len("1v0hr", length.out = 16),
                       rep_len("6v0hr", length.out = 16))

ggplot(df_to_plot,
       aes(x = as.factor(source),
           y = as.factor(cell_type),
           size = pct.1 * 100,
           fill = 2^avg_diff * pct.1/pct.2)) +
  geom_point(shape = 21) +
  ylab("cluster") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "darkred")) +
  facet_grid(cols = vars(gene_symbol)) + 
  theme_bw() 


# volcano plots ####
DEG_5line_0.05[[2]]$significance <- "nonsignificant"
DEG_5line_0.05[[2]]$significance[DEG_5line_0.05[[2]]$p_val_adj < 0.05 & DEG_5line_0.05[[2]]$avg_diff > 0] <- "up"
DEG_5line_0.05[[2]]$significance[DEG_5line_0.05[[2]]$p_val_adj < 0.05 & DEG_5line_0.05[[2]]$avg_diff < 0] <- "down"
unique(DEG_5line_0.05[[2]]$significance)
DEG_5line_0.05[[2]]$significance <- factor(DEG_5line_0.05[[2]]$significance, levels = c("up", "nonsignificant", "down"))
DEG_5line_0.05[[2]]$neg_log_pval <- (0 - log2(DEG_5line_0.05[[2]]$p_val))
DEG_5line_0.05[[2]]$labelling <- ""
for (i in c("BDNF", "FOS", "EGR1")){
  DEG_5line_0.05[[2]]$labelling[DEG_5line_0.05[[2]]$gene_symbol %in% i] <- i
}
for (t in unique(DEG_5line_0.05[[2]]$cell_type)){
  print(t)
  jpeg(paste0(t, "_6v0_volcano_plot.jpeg"))
  p <- ggplot(data = DEG_5line_0.05[[2]][DEG_5line_0.05[[2]]$cell_type == t, ],
              aes(x = avg_diff, 
                  y = neg_log_pval, 
                  color = significance,
                  label = labelling)) + 
    geom_point(size = 0.2) +
    scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
    theme_minimal() +
    geom_text_repel(box.padding = unit(0.05, 'lines'),
                    min.segment.length = 0,
                    force = 2,
                    max.overlaps = 8000,
                    force_pull = 0.5,
                    show.legend = F) +
    ggtitle(paste0("DEG at 0hr for ", t, " cells"))
  print(p)
  dev.off()
}



# test if expression in one population is 0 ####
DefaultAssay(integrated_labeled)
Idents(integrated_labeled) <- "time.ident"

# use scale.data first
df_1vs0 <- FindMarkers(object = integrated_labeled,
                       slot = "scale.data",
                       ident.1 = "1hr",
                       ident.2 = "0hr",
                       min.pct = 0.0,
                       logfc.threshold = 0.0,
                       test.use = "MAST",
                       features = c("NPAS4", "FOS"))
df_6vs0 <- FindMarkers(object = integrated_labeled,
                       slot = "scale.data",
                       ident.1 = "6hr",
                       ident.2 = "0hr",
                       min.pct = 0.0,
                       logfc.threshold = 0.0,
                       test.use = "MAST",
                       features = c("NPAS4", "FOS"))

# use RNA assay without scaling
DefaultAssay(integrated_labeled) <- "RNA"
df_1vs0 <- FindMarkers(object = integrated_labeled,
                       ident.1 = "1hr",
                       ident.2 = "0hr",
                       min.pct = 0.0,
                       logfc.threshold = 0.0,
                       test.use = "MAST",
                       features = c("NPAS4", "FOS"))
df_6vs0 <- FindMarkers(object = integrated_labeled,
                       ident.1 = "6hr",
                       ident.2 = "0hr",
                       min.pct = 0.0,
                       logfc.threshold = 0.0,
                       test.use = "MAST",
                       features = c("NPAS4", "FOS"))
