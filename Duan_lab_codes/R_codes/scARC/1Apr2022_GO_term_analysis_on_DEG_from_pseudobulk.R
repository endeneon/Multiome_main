# Chuxuan Li 04/01/2022
# GO term analysis on the DEG lists obtained from pseudobulk analysis on 5-line
#and 18-line combined data

# init ####
library(stringr)
library(readr)
library(ggplot2)
library(RColorBrewer)

# filter and output gene lists ####
pathlist <- list.files(path = "./pseudobulk_DE/geneonly", pattern = ".txt", full.names = T)
filelist <- vector("list", length(pathlist))
for(i in 1:length(pathlist)){
  filelist[[i]] <- read_csv(pathlist[i])
}

