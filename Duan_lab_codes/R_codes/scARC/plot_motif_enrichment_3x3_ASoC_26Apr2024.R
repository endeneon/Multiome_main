# Siwei 26 Apr 2024
# plot motif enrichment - 3x3 ASoC, 
# Use hierachical clustering;
# Use top 20?

# init ####
{
  library(readxl)
  
  library(stringr)
  library(reshape2)
  
  library(ggplot2)
  library(gplots)
  
  library(RColorBrewer)
  library(viridis)
  
}

# load data ####
df_raw <-
  read_excel("ASoC_SNP_motif_enrich_all.xlsx")

df_2_plot <-
  df_raw
df_2_plot <-
  df_2_plot[!duplicated(df_2_plot$`Motif Name`), ]

motif_names <-
  df_2_plot$`Motif Name`

df_2_plot$`Motif Name` <- NULL
rownames(df_2_plot) <- motif_names
df_2_plot <-
  log1p(df_2_plot)
rownames(df_2_plot) <-
  motif_names

max(df_2_plot)
min(df_2_plot)
colour_scale <-
  c(seq(min(df_2_plot),
        log1p(5),
        length = 10),
    seq(log1p(5),
        max(df_2_plot),
        length = 100))

grey_viridis <-
  c(rep_len("#888888",
            length.out = 10),
    viridis::magma(100))

# rownames(df_2_plot) <-
#   str_split(string = rownames(df_2_plot),
#             pattern = '/',
#             simplify = T)[, 1]

# dev.off()
heatmap.2(as.matrix(df_2_plot),
          Colv = F,
          symm = F,
          symkey = F,
          symbreaks = F,
          scale = "none",
          col = grey_viridis,
          dendrogram = "row",
          margins = c(8, 12),
          density.info = "density",
          key.title = 'Transformed -log10P',
          cexRow = 0.6,
          keysize = 1,
          lmat = rbind(c(0, 4),
                       c(2, 1),
                       c(0, 3)),
          lhei = c(1, 4, 0.1),
          trace = "none")
dev.off()

data(mtcars)
x  <- as.matrix(mtcars)
dev.off()

hv <- 
  heatmap.2(x, 
            col = bluered, 
            scale = "column", 
            trace = "none")

whiteBin <- 
  unlist(hv$colorTable[hv$colorTable[, "color"] == "green4",
                       1:2])
rbind(whiteBin[1] * hv$colSDs + hv$colMeans,
      whiteBin[2] * hv$colSDs + hv$colMeans )
