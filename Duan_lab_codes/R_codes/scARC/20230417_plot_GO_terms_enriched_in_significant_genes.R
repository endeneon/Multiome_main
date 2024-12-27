# Chuxuan Li 04/17/2023
# plot top GO terms from case vs control DEs

# init ####
library(ggplot2)
library(RColorBrewer)
library(scales)
library(readr)
library(stringr)

up_txtlist <- list.files("./MAST_scDE/DE_results/GO_terms_on_significant_genes/",
                         pattern = "upregulated", full.names = T)
down_txtlist <- list.files("./MAST_scDE/DE_results/GO_terms_on_significant_genes/",
                         pattern = "downregulated", full.names = T)
up_names <- str_remove(str_remove(up_txtlist, "_in_case.txt"),
                       "\\./MAST_scDE/DE_results/GO_terms_on_significant_genes//")
down_names <- str_remove(str_remove(down_txtlist, "_in_case.txt"),
                         "\\./MAST_scDE/DE_results/GO_terms_on_significant_genes//")
up_dfs <- vector("list", length(up_txtlist))
names(up_dfs) <- up_names
down_dfs <- vector("list", length(down_txtlist))
names(down_dfs) <- down_names

for (i in 1:length(up_txtlist)) {
  up_df <- read_delim(up_txtlist[i], 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE, skip = 11)
  colnames(up_df) <- c("GOterm", "reflist", "actual", "expected", "over_under",
                       "fold_enrichment", "pvalue", "FDR")
  up_df$GOterm <- str_remove(str_remove(up_df$GOterm, " \\(GO:[0-9]+\\)"), " \\([A-Z]+\\)")
  up_df$neglogp <- -log10(up_df$pvalue)
  up_df <- up_df[order(up_df$fold_enrichment, decreasing = T), ]
  up_dfs[[i]] <- up_df
  down_df <- read_delim(down_txtlist[i], 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE, skip = 11)
  colnames(down_df) <- c("GOterm", "reflist", "actual", "expected", "over_under",
                       "fold_enrichment", "pvalue", "FDR")
  down_df$GOterm <- str_remove(str_remove(down_df$GOterm, " \\(GO:[0-9]+\\)"), 
                                " \\([A-Z]+\\)")
  down_df$neglogp <- -log10(down_df$pvalue)
  down_df <- down_df[order(down_df$fold_enrichment, decreasing = T), ]
  down_dfs[[i]] <- down_df
}

# plot ####
for (i in 1:length(up_dfs)) {
  up_df <- up_dfs[[i]][1:10, ]
  orders <- rev(up_df$GOterm)
  up_df$GOterm <- factor(up_df$GOterm, levels = orders)
  up_df$fold_enrichment <- as.numeric(up_df$fold_enrichment)
  png(paste0("./MAST_scDE/DE_results/GO_terms_bargraphs/", 
             names(up_dfs)[i], "_top10_GOterms_bargraph.png"), width = 900,
      height = 400)
  p <- ggplot(up_df, mapping = aes(x = fold_enrichment, y = GOterm)) +
    geom_col(aes(fill = neglogp), position = "dodge") +
    scale_fill_gradient(low = "red", high = "blue") +
    ggtitle(paste0("top 10 GO terms enriched in\n", 
                   str_replace_all(names(up_dfs)[i], "_", " "), " genes")) +
    theme(text = element_text(size = 10), axis.text.y = element_text(size = 12))
  print(p)
  dev.off()
  
  down_df <- down_dfs[[i]][1:10, ]
  orders <- rev(down_df$GOterm)
  down_df$GOterm <- factor(down_df$GOterm, levels = orders)
  down_df$fold_enrichment <- as.numeric(down_df$fold_enrichment)
  png(paste0("./MAST_scDE/DE_results/GO_terms_bargraphs/", 
             names(down_dfs)[i], "_top10_GOterms_bargraph.png"), width = 900,
      height = 400)
  p <- ggplot(down_df, mapping = aes(x = fold_enrichment, y = GOterm)) +
    geom_col(aes(fill = neglogp), position = "dodge") +
    scale_fill_gradient(low = "red", high = "blue") +
    ggtitle(paste0("top 10 GO terms enriched in\n", 
                   str_replace_all(names(down_dfs)[i], "_", " "), " genes")) +
    theme(text = element_text(size = 10), axis.text.y = element_text(size = 12))
  print(p)
  dev.off()
}

