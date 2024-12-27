# Chuxuan Li 12/15/2022
# Plot 025 psuedobulk da peak results

# init ####
library(readr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(graphics)
library(stringr)

set.seed(2022)

# load annotated bed files ####
pathlist <- list.files(path = "./dapeaks_limma/unfiltered_bed/annotated/", pattern = "full",
                      full.names = T)
reslist <- vector("list", length(pathlist))
for (i in 1:length(pathlist)) {
  reslist[[i]] <- read_delim(pathlist[i], delim = "\t", col_names = F, skip = 1, escape_double = F)
}
for (i in 1:length(reslist)){
  reslist[[i]] <- reslist[[i]][reslist[[i]]$X18 %in% c("FOS", "NPAS4", "BDNF", "VGF"), ]
}
reslist[[2]]$X18

# extract peaks for response genes ####
onegenelist <- vector("list", length(reslist))
for (i in 1:length(reslist)) {
  genes.exist <- unique(reslist[[i]]$X18)
  if (length(genes.exist) != 0){
    for (j in 1:length(genes.exist)) {
      tempres <- reslist[[i]][reslist[[i]]$X18 == genes.exist[j], ]
      tempres <- tempres[order(tempres$X1, decreasing = T), ]
      tempres <- tempres[1, ] #pick the peak with highest pvalue
      tempres[, 19] <- paste(tempres$X4, tempres$X5 - 1, tempres$X6, sep = "-")
      tempres <- tempres[, 18:19]
      if (j == 1) {
        onegenetemp <- tempres
      } else {
        onegenetemp <- rbind(onegenetemp, tempres)
      }
    }
    onegenelist[[i]] <- as.data.frame(onegenetemp)
  }
}
names(onegenelist) <- str_extract(pathlist, "[A-Za-z]+_[1|6]v0hr_[a-z]+")
genes.use.1v0 <- onegenelist[seq(1, 6, 2)]
genes.use.6v0 <- onegenelist[seq(2, 6, 2)]

# plot volcano plots ####
names_1v0 <- c("GABA_1v0hr", "nmglut_1v0hr", "npglut_1v0hr")
names_6v0 <- c("GABA_6v0hr", "nmglut_6v0hr", "npglut_6v0hr")
colors <- c("steelblue3", "grey", "indianred3")
da_peaks_by_time_list_1v0 <- vector("list", length(fit_all))
da_peaks_by_time_list_6v0 <- vector("list", length(fit_all))
for (i in 1:length(fit_all)){
  da_peaks_by_time_list_1v0[[i]] <- topTable(fit_all[[i]], coef = "onevszero", p.value = Inf, 
                                             sort.by = "P", number = Inf)
  da_peaks_by_time_list_6v0[[i]] <- topTable(fit_all[[i]], coef = "sixvszero", p.value = Inf, 
                                             sort.by = "P", number = Inf)
}

for (i in 1:3){
  name1 <- names_1v0[i]
  print(name1)
  da_peaks_by_time_list_1v0[[i]]$peak <- rownames(da_peaks_by_time_list_1v0[[i]])
  da_peaks_by_time_list_1v0[[i]]$genes <- NULL
  df <- da_peaks_by_time_list_1v0[[i]][, c("peak", "AveExpr", "logFC", "P.Value", "adj.P.Val")]
  df$significance <- "nonsig"
  df$significance[df$adj.P.Val < 0.05 & df$logFC > 0] <- "pos"
  df$significance[df$adj.P.Val < 0.05 & df$logFC < 0] <- "neg"
  df$neg_log_p_val <- (0 - log10(df$P.Value))
  df$labelling <- ""
  for (j in 1:length(genes.use.1v0[[i]])){
    df$labelling[df$peak %in% genes.use.1v0[[i]]$X19] <- genes.use.1v0[[i]]$X18
  }
  print(unique(df$labelling))
  
  # plot
  file_name <- paste0("./dapeaks_limma/", name1, ".png")
  png(file_name)
  p <- ggplot(data = df,
              aes(x = logFC, 
                  y = neg_log_p_val, 
                  label = labelling,
                  color = significance)) + 
    geom_point() +
    scale_color_manual(values = colors) +
    theme_minimal() +
    ggtitle(str_replace_all(name1, "_", " ")) +
    geom_text_repel(box.padding = unit(0.05, 'lines'),
                    min.segment.length = 0,
                    force = 2,
                    force_pull = 10,
                    max.overlaps = 100000) +
    xlim(-5, 5) 
  print(p)
  dev.off()
  
  name6 <- names_6v0[i]
  print(name6)
  da_peaks_by_time_list_6v0[[i]]$peak <- rownames(da_peaks_by_time_list_6v0[[i]])
  da_peaks_by_time_list_6v0[[i]]$genes <- NULL
  df <- da_peaks_by_time_list_6v0[[i]][, c("peak", "AveExpr", "logFC", "P.Value", "adj.P.Val")]
  df$significance <- "nonsig"
  df$significance[df$adj.P.Val < 0.05 & df$logFC > 0] <- "pos"
  df$significance[df$adj.P.Val < 0.05 & df$logFC < 0] <- "neg"
  df$neg_log_p_val <- (0 - log10(df$P.Value))
  df$labelling <- ""
  for (j in 1:length(genes.use.6v0[[i]])){
    df$labelling[df$peak %in% genes.use.6v0[[i]]$X19] <- genes.use.6v0[[i]]$X18
  }
  print(unique(df$labelling))
  file_name <- paste0("./dapeaks_limma/", name6, ".png")
  png(file_name)
  p <- ggplot(data = df,
              aes(x = logFC, 
                  y = neg_log_p_val, 
                  label = labelling,
                  color = significance)) + 
    geom_point() +
    scale_color_manual(values = colors) +
    theme_minimal() +
    ggtitle(str_replace_all(name6, "_", " ")) +
    geom_text_repel(box.padding = unit(0.05, 'lines'),
                    min.segment.length = 0,
                    force = 2,
                    force_pull = 10,
                    max.overlaps = 100000) +
    xlim(-5, 5) 
  print(p)
  dev.off()
}

# load annotation stats ####
pathlist <- list.files(path = "./dapeaks_limma/unfiltered_bed/stats/", pattern = "full",
                      full.names = T)
statlist <- vector("list", length(pathlist))
for (i in 1:length(pathlist)) {
  statlist[[i]] <- read_delim("dapeaks_limma/unfiltered_bed/stats/GABA_1v0hr_fullstats.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
}
