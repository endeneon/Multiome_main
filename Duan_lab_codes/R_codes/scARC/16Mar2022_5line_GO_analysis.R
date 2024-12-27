# Chuxuan Li 03/16/2022
# Do GO term analysis on pseudobulk-5line DEGs

# init ####
library(readr)
library(stringr)
library(DESeq2)
library(ggplot2)

# load pseudobulk data ####

# load and count full list
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_pseudobulk/data_frames_for_pseudobulk_raw.RData")

df_list <- list(GABA_pseudobulk, nmglut_pseudobulk,
                npglut_pseudobulk, NPC_pseudobulk)
types <- c("GABA", "nmglut",
           "npglut", "NPC")
res_0v1_list <- vector(mode = "list", length = length(df_list))
names(res_0v1_list) <- types
res_0v6_list <- vector(mode = "list", length = length(df_list))
names(res_0v6_list) <- types

for (i in 1:length(df_list)){
  # make colData df
  coldata <- data.frame(condition = str_extract(string = colnames(df_list[[i]]),
                                                pattern = "[0-6]hr"),
                        type = str_extract(string = colnames(df_list[[i]]),
                                           pattern = "CD_[0-9][0-9]"))
  rownames(coldata) <- colnames(df_list[[i]])
  
  # check
  print(all(rownames(coldata) %in% colnames(df_list[[i]]))) # TRUE
  print(all(rownames(coldata) == colnames(df_list[[i]]))) #TRUE
  
  # make deseq2 obj
  print(i)
  dds <- DESeqDataSetFromMatrix(countData = df_list[[i]],
                                colData = coldata,
                                design = ~ type + condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "1hr", "0hr"))
  res <- as.data.frame(res)
  res <- res[res$baseMean > 10, ]
  res$q_value <- p.adjust(p = res$pvalue, method = "BH")
  res$gene <- rownames(res)
  res_0v1_list[[i]] <- res
  
  res <- results(dds, contrast = c("condition", "6hr", "0hr"))
  res <- as.data.frame(res)
  res <- res[res$baseMean > 10, ]
  res$q_value <- p.adjust(p = res$pvalue, method = "BH")
  res$gene <- rownames(res)
  res_0v6_list[[i]] <- res
}

length(res_0v1_list[[4]]$gene)

setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_pseudobulk/5line_DE_results/filtered_by_basemean/gene_only_lists_after_basemean_filtering_only")
for (i in 1:length(res_0v1_list)){
  filename <- paste0(types[i], "_1v0_genelist.txt")
  print(filename)
  write.table(res_0v1_list[[i]]$gene, file = filename, 
              quote = F, sep = ",", col.names = F, row.names = F)
  
  filename <- paste0(types[i], "_6v0_genelist.txt")
  print(filename)
  write.table(res_0v6_list[[i]]$gene, file = filename, 
              quote = F, sep = ",", col.names = F, row.names = F)
}

# after filtering and separating up/down
setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_pseudobulk/DE_results/filtered_by_basemean/significant")
pseudo_path_list <- list.files(path = "./")
pseudo_df_list <- vector(mode = "list", length = length(pseudo_path_list))
for (i in 1:length(pseudo_path_list)){
  pseudo_df_list[[i]] <- read.csv(pseudo_path_list[i], 
                                  col.names = c("gene", "baseMean", "log2FC", "lfcSE", "stat", "pValue", "padj", "qvalue"))
}

# after GO enrichment on PANTHER ####
setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_pseudobulk/5line_GO_results/background_genes_updated")# load data
plist <- list.files("./")
dflist <- vector(mode = "list", length = length(plist))
for (i in 1:length(plist)){
  df <- read_delim(plist[i], col_names = T)
  colnames(df) <- c("GOterm", "ref", "#", "expected", "over_under", "fold_enrichment",    
                    "p_value", "FDR")
  dflist[[i]] <- df
}

termlist <- c("axonal transport (GO:0098930)",
              "regulation of neuronal synaptic plasticity (GO:0048168)",
              "axo-dendritic transport (GO:0008088)",
              "signal release from synapse (GO:0099643)",
              "regulation of neuron migration (GO:2001222)",
              "neuron projection organization (GO:0106027)",
              "neuron projection morphogenesis (GO:0048812)",	
              "positive regulation of axon extension (GO:0045773)",
              "synaptic vesicle cycle (GO:0099504)",
              "synapse assembly (GO:0007416)",
              "dendrite morphogenesis (GO:0048813)",
              "neurotransmitter transport (GO:0006836)",
              "regulation of axon extension (GO:0030516)",
              "neural tube development (GO:0021915)",
              "regulation of neuron projection development (GO:0010975)",
              "axonogenesis (GO:0007409)",
              "nervous system development (GO:0007399)",
              "axon development (GO:0061564)",	
              "neuron differentiation (GO:0030182)",	
              "brain development (GO:0007420)",
              "central nervous system development (GO:0007417)",
              "neurogenesis (GO:0022008)",	
              "regulation of neuron projection development (GO:0010975)",
              "neuron migration (GO:0001764)",	
              "synaptic vesicle cycle (GO:0099504)",	
              "synaptic vesicle exocytosis (GO:0016079)",
              "neurotransmitter secretion (GO:0007269)",
              "regulation of neurotransmitter secretion (GO:0046928)",	
              "regulation of neurotransmitter transport (GO:0051588)",
              "regulation of synaptic plasticity (GO:0048167)",
              "regulation of synapse organization (GO:0050807)",
              "regulation of synapse structure or activity (GO:0050803)",
              "positive regulation of neuron projection development (GO:0010976)",
              "regulation of axonogenesis (GO:0050770))")

filteredlist <- vector(mode = "list", length = length(dflist))
for (i in 1:length(dflist)){
  filteredlist[[i]] <- dflist[[i]][dflist[[i]]$GOterm %in% termlist, ] 
}
unique_terms <- c()
for (i in 1:length(filteredlist)){
  unique_terms <- c(unique_terms, filteredlist[[i]]$GOterm)
}
unique_terms <- unique(unique_terms)
unique_terms
df_to_plot <- data.frame(GO_term = rep_len(unique_terms, length(dflist)),
                         neg_log_pval = rep_len(NA, length(unique_terms) * length(dflist)),
                         fold_enrichment = rep_len(NA, length(unique_terms) * length(dflist)),
                         type_time = rep_len("", length(unique_terms) * length(dflist)))
for (i in 1:length(dflist)){
  df <- dflist[[i]]
  for (j in 1:length(unique_terms)){
    print(unique_terms[j])
    df_to_plot$GO_term[length(unique_terms) * i - length(unique_terms) + j] <- 
      str_remove(unique_terms[j], " \\(GO:[0-9]+\\)")
    if (unique_terms[j] %in% df$GOterm){
      subdf <- df[df$GOterm == unique_terms[j], ]
      df_to_plot$neg_log_pval[length(unique_terms) * i - length(unique_terms) + j] <- 
        0 - log10(subdf$p_value)
      df_to_plot$fold_enrichment[length(unique_terms) * i - length(unique_terms) + j] <- 
        as.numeric(subdf$fold_enrichment)
    }
    df_to_plot$type_time[length(unique_terms) * i - length(unique_terms) + j] <- 
      str_remove(string = plist[i], pattern = ".txt")
  }
}
# df_to_plot <- df_to_plot[c(sequence(nvec = rep_len(6, 16), from = seq(1, 96, 12)),
#                            sequence(nvec = rep_len(6, 16), from = seq(7, 96, 12))), ]
# df_to_plot$num <- factor(df_to_plot$type_time, levels = unique(df_to_plot$type_time))
# df_to_plot <- df_to_plot[sequence(nvec = rep_len(6, 16), from = seq(7, 96, 12)), ]

ggplot(data = df_to_plot, mapping = aes(x = type_time,
                                        y = as.factor(GO_term),
                                        fill = neg_log_pval,
                                        size = fold_enrichment)) +
  geom_point(shape = 21) +
  scale_fill_gradientn(colours = c("white", "royalblue3")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab("Cell type/time points") +
  ylab("GO term") + 
  labs(fill = "-log(p-value)",
       size = "fold enrichment") 
