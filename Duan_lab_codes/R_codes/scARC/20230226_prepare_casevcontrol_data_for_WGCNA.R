# Chuxuan Li 02/26/2022
# Prepare case vs control count data for WGCNA analysis

# init ####
library(WGCNA)
library(stringr)
library(edgeR)

options(stringsAsFactors = F)
load("~/NVME/scARC_Duan_018/case_control_DE_analysis/018-029_by_type_and_time_mat_lst_bytime.RData")

# clean data ####
cleaned_lst <- mat_lst_bytime
gene_count_summary <- matrix(nrow = 2, ncol = length(mat_lst_bytime), 
                            dimnames = list(c("total", "cpm > 1 in at least half of the samples"),
                                            names(mat_lst_bytime)))
for (i in 1:length(mat_lst_bytime)) {
  mat <- mat_lst_bytime[[i]]
  cpm <- cpm(mat)
  passfilter <- rowSums(cpm > 1) > 0.5 * ncol(cpm)
  gene_count_summary[1, i] <- nrow(cpm)
  gene_count_summary[2, i] <- sum(passfilter)
  mat <- mat[passfilter, ]
  mat_transposed <- t(mat)
  cleaned_lst[[i]] <- mat_transposed
}

# hierarchical clustering and cleaning of the samples ####
# GABA
sampleTree <- hclust(dist(cleaned_lst[[1]]), method = "average")
plot(sampleTree, main = "sample clustering to detect outliers - GABA 0hr", sub = "", 
     xlab = "")
abline(h = 3.6e+05, col = "red") 
clust <- cutreeStatic(sampleTree, cutHeight = 3.6e+05) 
table(clust) 
cleaned_lst[[1]] <- cleaned_lst[[1]][clust != 0, ] # remove 6 samples

sampleTree <- hclust(dist(cleaned_lst[[2]]), method = "average")
plot(sampleTree, main = "sample clustering to detect outliers - GABA 1hr", sub = "", 
     xlab = "")
abline(h = 3e+05, col = "red") 
clust <- cutreeStatic(sampleTree, cutHeight = 3e+05) 
table(clust) 
cleaned_lst[[2]] <- cleaned_lst[[2]][clust != 0, ] # remove 3 samples

sampleTree <- hclust(dist(cleaned_lst[[3]]), method = "average")
plot(sampleTree, main = "sample clustering to detect outliers - GABA 6hr", sub = "", 
     xlab = "")
abline(h = 6e+05, col = "red") 
clust <- cutreeStatic(sampleTree, cutHeight = 6e+05) 
table(clust) 
cleaned_lst[[3]] <- cleaned_lst[[3]][clust != 0, ] # remove 1 sample


# nmglut
sampleTree <- hclust(dist(cleaned_lst[[4]]), method = "average")
plot(sampleTree, main = "sample clustering to detect outliers - NEFM- glut 0hr", sub = "", 
     xlab = "")
abline(h = 3e+05, col = "red") 
clust <- cutreeStatic(sampleTree, cutHeight = 3e+05) 
table(clust) 
cleaned_lst[[4]] <- cleaned_lst[[4]][clust != 0, ] # remove 2 samples

sampleTree <- hclust(dist(cleaned_lst[[5]]), method = "average")
plot(sampleTree, main = "sample clustering to detect outliers - NEFM- glut 1hr", sub = "", 
     xlab = "")
abline(h = 3e+05, col = "red") 
clust <- cutreeStatic(sampleTree, cutHeight = 3e+05) 
table(clust) 
cleaned_lst[[5]] <- cleaned_lst[[5]][clust != 0, ] # remove 5 samples

sampleTree <- hclust(dist(cleaned_lst[[6]]), method = "average")
plot(sampleTree, main = "sample clustering to detect outliers - NEFM- glut 6hr", sub = "", 
     xlab = "")
abline(h = 250000, col = "red") 
clust <- cutreeStatic(sampleTree, cutHeight = 250000) 
table(clust) 
cleaned_lst[[6]] <- cleaned_lst[[6]][clust != 0, ] # remove 2 samples


# npglut
sampleTree <- hclust(dist(cleaned_lst[[7]]), method = "average")
plot(sampleTree, main = "sample clustering to detect outliers - NEFM+ glut 0hr", sub = "", 
     xlab = "")
abline(h = 1e+06, col = "red") 
clust <- cutreeStatic(sampleTree, cutHeight = 1e+06) 
table(clust) 
cleaned_lst[[7]] <- cleaned_lst[[7]][clust != 0, ] # remove 1 sample

sampleTree <- hclust(dist(cleaned_lst[[8]]), method = "average")
plot(sampleTree, main = "sample clustering to detect outliers - NEFM+ glut 1hr", sub = "", 
     xlab = "")
abline(h = 6e+05, col = "red") 
clust <- cutreeStatic(sampleTree, cutHeight = 6e+05) 
table(clust) 
cleaned_lst[[8]] <- cleaned_lst[[8]][clust != 0, ] # remove 4 samples

sampleTree <- hclust(dist(cleaned_lst[[9]]), method = "average")
plot(sampleTree, main = "sample clustering to detect outliers - NEFM+ glut 6hr", sub = "", 
     xlab = "")
abline(h = 5.5e+05, col = "red") 
clust <- cutreeStatic(sampleTree, cutHeight = 5.5e+05) 
table(clust) 
cleaned_lst[[9]] <- cleaned_lst[[9]][clust != 0, ] # remove 1 sample


save(cleaned_lst, file = "count_matrices_filtered_genes_removed_outlier_samples.RData")

# prepare metadata df ####
load("~/NVME/scARC_Duan_018/case_control_DE_analysis/018-029_covar_table_final.RData")
trait_df <- covar_table_final[order(covar_table_final$cell_line), ]
trait_df_final <- trait_df[, -c(1, 2)]
trait_df_final$sex <- str_replace(trait_df_final$sex, "F", "1")
trait_df_final$sex <- str_replace(trait_df_final$sex, "M", "2")
trait_df_final$aff <- str_replace(trait_df_final$aff, "case", "1")
trait_df_final$aff <- str_replace(trait_df_final$aff, "control", "0")
trait_df_final <- lapply(trait_df_final, as.numeric)
trait_df_final <- as.data.frame(trait_df_final)
rownames(trait_df_final) <- paste(trait_df$cell_line, trait_df$time, sep = "_")
collectGarbage()

# visualize metadata + expression data ####
trait_df_list <- vector("list", length(cleaned_lst))
trait_colors_list <- vector("list", length(cleaned_lst))
for (i in 1:length(trait_df_list)) {
  sampleTree2 <- hclust(dist(cleaned_lst[[i]]), method = "average")
  trait_df_list[[i]] <- trait_df_final[rownames(trait_df_final) %in% rownames(cleaned_lst[[i]]), ]
  trait_colors_list[[i]] <- numbers2colors(trait_df_list[[i]], signed = F)
  png(paste0("./sample_dendro_covar_heatmap/dendroheatmap_", names(cleaned_lst)[i],
            ".png"), width = 1200, height = 500)
  p <- plotDendroAndColors(sampleTree2, trait_colors_list[[i]], 
                           groupLabels = colnames(trait_df_list[[i]]),
                      main = paste0("sample dendrogram and covariables heatmap - ",
                                    names(cleaned_lst)[i]))
  print(p)
  dev.off()
}

save(trait_df_list, file = "trait_dfs.RData")
save(trait_colors_list, file = "trait_colors_list.RData")
