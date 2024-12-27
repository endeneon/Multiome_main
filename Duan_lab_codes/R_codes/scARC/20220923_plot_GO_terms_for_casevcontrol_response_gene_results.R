# Chuxuan Li 09/23/2022
# Plot GO terms from quantile normalized, limma then lm fitted, filtered by p-value 
#differential response analysis results; GO terms also filtered by p-value only

# init ####
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# read files
paths <- list.files(path = "case_v_control_DE/qnres/GOterms/", full.names = T)
names <- str_extract(paths, "[A-Za-z]+_[1|6]v[0|1|6]")
txts <- vector("list", length(paths))
for (i in 1:length(paths)) {
  txts[[i]] <- read_delim(paths[i], delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)
  colnames(txts[[i]]) <- 
    c("GO_biological_process", "REF", "n_genes", "expected", "over/under", 
      "Fold_Enrichment", "p_value", "FDR")
  txts[[i]]$GO_biological_process <- str_remove(txts[[i]]$GO_biological_process,
                                                " \\(GO:[0-9]+\\)")
}
names(txts) <- names

dfs_1v0 <- txts[str_detect(names, "1v0")]
dfs_6v0 <- txts[str_detect(names, "6v0")]
dfs_6v1 <- txts[str_detect(names, "6v1")]

# count significant ####
types1 <- names(dfs_1v0)
types6 <- names(dfs_6v0)
types61 <- names(dfs_6v1)
deg_counts1 <- array(dim = c(3, 4), 
                     dimnames = list(types1, c("positive", "negative", "nonsig", "total"))) 
deg_counts6 <- array(dim = c(3, 4), 
                     dimnames = list(types6, c("positive", "negative", "nonsig", "total"))) 
deg_counts61 <- array(dim = c(3, 4), 
                      dimnames = list(types61, c("positive", "negative", "nonsig", "total"))) 

for (i in 1:length(dfs_1v0)) {
  deg_counts1[i,1] <- sum(dfs_1v0[[i]]$FDR < 0.05 & dfs_1v0[[i]]$`over/under` == "+")
  deg_counts1[i,2] <- sum(dfs_1v0[[i]]$FDR < 0.05 & dfs_1v0[[i]]$`over/under` == "-")
  deg_counts1[i,3] <- sum(dfs_1v0[[i]]$FDR > 0.05)
  deg_counts1[i,4] <- nrow(dfs_1v0[[i]])
  deg_counts6[i,1] <- sum(dfs_6v0[[i]]$FDR < 0.05 & dfs_6v0[[i]]$`over/under` == "+")
  deg_counts6[i,2] <- sum(dfs_6v0[[i]]$FDR < 0.05 & dfs_6v0[[i]]$`over/under` == "-")
  deg_counts6[i,3] <- sum(dfs_6v0[[i]]$FDR > 0.05)
  deg_counts6[i,4] <- nrow(dfs_6v0[[i]])
  deg_counts61[i,1] <- sum(dfs_6v1[[i]]$FDR < 0.05 & dfs_6v1[[i]]$`over/under` == "+")
  deg_counts61[i,2] <- sum(dfs_6v1[[i]]$FDR < 0.05 & dfs_6v1[[i]]$`over/under` == "-")
  deg_counts61[i,3] <- sum(dfs_6v1[[i]]$FDR > 0.05)
  deg_counts61[i,4] <- nrow(dfs_6v1[[i]])
} 
deg_counts <- rbind(deg_counts1, deg_counts6, deg_counts61)
write.table(deg_counts, 
            file = "./case_v_control_DE/qnres/GOterms/summary_counts_by_padj.csv", 
            quote = F, sep = ",", row.names = T, col.names = T)

# find common terms ####
top10 <- vector("list", length(dfs_6v0))
for (i in 1:length(top10)) {
  top10[[i]] <- arrange(dfs_6v0[[i]], p_value)[1:10, ]
  if (i == 1) {
    common_terms <- top10[[i]]$GO_biological_process
  } else {
    common_terms <- union(common_terms, top10[[i]]$GO_biological_process)
  }
}
names(top10) <- names(dfs_6v0)
for (i in 1:length(top10)) {
  top10[[i]]$cell.type <- names(top10)[i]
  top10[[i]]$neglogp <- (-1) * log2(top10[[i]]$p_value)
  top10[[i]]$Fold_Enrichment <- as.numeric(top10[[i]]$Fold_Enrichment)
  if (i == 1) {
    big_df1 <- top10[[i]]
  } else {
    big_df1 <- rbind(big_df1, top10[[i]])
  }
}
# plot ####
ggplot(big_df1, aes(x = GO_biological_process, y = cell.type, 
                    size = neglogp, color = Fold_Enrichment)) +
  geom_point()+
  scale_color_gradientn(colours = c("white", "red", "red2", "red3", "darkred"),
                        na.value = "transparent") +
  labs(size = "-log2(p-value)", color = "Fold\nEnrichment", x = "GO Terms", y = "Cell Type") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 7)) +
  ggtitle("Top 10 GO terms from each cell type\nEnriched in differential response genes at 6v0hr")
