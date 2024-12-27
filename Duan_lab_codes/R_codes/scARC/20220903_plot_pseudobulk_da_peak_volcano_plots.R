library(readr)
pathlist <-list.files(path = "./da_peaks_limma/bed/annotated/", pattern = ".tsv", full.names = T)
reslist <- vector("list", length(pathlist))
for (i in 1:length(pathlist)) {
  reslist[[i]] <- read_delim(pathlist[i], delim = "\t", col_names = F, skip = 1, escape_double = F)
}
reslist[[1]]$X18
for (i in 1:length(reslist)){
  reslist[[i]] <- reslist[[i]][reslist[[i]]$X18 %in% c("FOS", "NPAS4", "BDNF"), ]
}
for (i in 1:length(reslist)) {
  print(unique(reslist[[i]]$X18))
}
reslist[[2]]
