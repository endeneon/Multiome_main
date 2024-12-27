# Chuxuan Li 04/18/2022
# Find the correlation coefficient between the expression of BDNF and accessibility
#of the peak 

# init ####
library(Signac)
library(Seurat)
library(stringr)
load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/multiomic_obj.RData")

DefaultAssay(multiomic_obj) <- "peaks"

glut <- subset(multiomic_obj, RNA.cell.type %in% c("NEFM_neg_glut", "NEFM_pos_glut"))
nmglut <- subset(multiomic_obj, RNA.cell.type %in% c("NEFM_neg_glut"))
npglut <- subset(multiomic_obj, RNA.cell.type %in% c("NEFM_pos_glut"))

# extract the two datasets ####
RNA <- glut@assays$SCT@data[rownames(glut@assays$SCT) == "BDNF", ]
glut@assays$peaks[str_detect(rownames(glut@assays$peaks), "chr11-277"), ]
ATAC <- glut@assays$peaks@data[rownames(glut@assays$peaks) == "chr11-27770213-27771025", ]
left1_ATAC <- glut@assays$peaks@data[rownames(glut@assays$peaks) == "chr11-27759385-27759809", ]
left2_ATAC <- glut@assays$peaks@data[rownames(glut@assays$peaks) == "chr11-27740518-27741389", ]
r1_ATAC <- glut@assays$peaks@data[rownames(glut@assays$peaks) == "chr11-27796008-27797123", ]

res <- cor.test(log1p(RNA), log1p(ATAC)) 
res_left1 <- cor.test(log1p(RNA), log1p(left1_ATAC))
res_left2 <- cor.test(log1p(RNA), log1p(left2_ATAC))
res_r1 <- cor.test(log1p(RNA), log1p(r1_ATAC))

# get a list of peaks 5e05 up and downstream of the peak of interest ####
RNA <- glut@assays$SCT@data[rownames(glut@assays$SCT) == "BDNF", ]
27770213 - 5e05
27770213 + 5e05
# start: chr11-27246928-27247249
# end: chr11-28295404-28296254
ATAC <- glut@assays$peaks@data[str_detect(rownames(glut@assays$peaks), 
                                          "chr11-2[7|8][0-9]{6}"), ]
rownames(ATAC)
df <- data.frame(p_value = rep_len(0, nrow(ATAC)),
                 cor = rep_len(0, nrow(ATAC)))
for (i in 1:nrow(ATAC)) {
  res <- cor.test(log1p(RNA), log1p(ATAC[i, ]))
  df$p_value[i] <- res$p.value
  df$cor[i] <- res$estimate
}
hist(-log10(df$p_value), breaks = 100)
hist(df$cor, breaks = 200)

# separate data into 3 time points and do calculation again 
# nmglut ####
unique(multiomic_obj$timextype.ident)
nmglut_0hr <- subset(multiomic_obj, timextype.ident == "NEFM_neg_glut_0hr")
nmglut_1hr <- subset(multiomic_obj, timextype.ident == "NEFM_neg_glut_1hr")
nmglut_6hr <- subset(multiomic_obj, timextype.ident == "NEFM_neg_glut_6hr")
nmglut_df_list <- vector(mode = "list", 3L)
nmglut_obj_list <- list(nmglut_0hr, nmglut_1hr, nmglut_6hr)
for (i in 1:length(nmglut_df_list)){
  RNA <- nmglut_obj_list[[i]]@assays$SCT@data[rownames(nmglut_obj_list[[i]]@assays$SCT) == "BDNF", ]
  ATAC <- nmglut_obj_list[[i]]@assays$peaks@data[str_detect(rownames(nmglut_obj_list[[i]]@assays$peaks), 
                                            "chr11-2[7|8][0-9]{6}"), ]
  df <- data.frame(p_value = rep_len(0, nrow(ATAC)),
                   cor = rep_len(0, nrow(ATAC)),
                   peak = rownames(ATAC),
                   zscore = rep_len(0, nrow(ATAC)))
  
  for (j in 1:nrow(ATAC)) {
    res <- cor.test(log1p(RNA), log1p(ATAC[j, ]))
    df$p_value[j] <- res$p.value
    df$cor[j] <- res$estimate
    df$zscore[j] <- res$statistic
  }
  df$q_value <- p.adjust(df$p_value, method = "BH")
  #df <- df[df$q_value < 0.05, ]
  #df <- df[df$q_value < 0.1, ]
  df <- df[df$p_value < 0.05, ]
  
  nmglut_df_list[[i]] <- df
  es <- (as.numeric(str_remove_all(str_extract(df$peak, "-[0-9]+-"), "-")) + 
                       as.numeric(str_remove(str_extract(df$peak, "[0-9]+$"), "-"))) / 2
  ss <- rep_len(27721989, nrow(df))
  for (k in 1:length(es)){
    e <- es[k]
    s <- ss[k]
    print(paste0("start: ", s, "end: ", e))
    if (e < s) {
      ss[k] <- e
      es[k] <- 27721989
      print("switched")
    }
  }
  gr <- GRanges(seqnames = str_extract(df$peak, "chr11"), 
                strand = "+",
                ranges = IRanges(start = ss,
                                 end = es))
  elementMetadata(gr)[["score"]] <- df$cor
  #elementMetadata(gr)[["pvalue"]] <- df$q_value
  elementMetadata(gr)[["pvalue"]] <- df$p_value
  elementMetadata(gr)[["gene"]] <- "BDNF"
  elementMetadata(gr)[["peak"]] <- df$peak
  elementMetadata(gr)[["zscore"]] <- df$zscore
  Links(nmglut_obj_list[[i]]) <- gr
  Annotation(nmglut_obj_list[[i]]) <- annotation_ranges
}

for (i in 1:length(nmglut_obj_list)){
  nmglut_obj_list[[i]]$timextype.ident.for.plot <- str_replace_all(nmglut_obj_list[[i]]$timextype.ident,
                                                               "_", " ")
  nmglut_obj_list[[i]]$timextype.ident.for.plot <- str_replace(nmglut_obj_list[[i]]$timextype.ident.for.plot,
                                                               " neg", "-")
  
  Idents(nmglut_obj_list[[i]]) <- "timextype.ident.for.plot"
  name <- unique(nmglut_obj_list[[i]]$timextype.ident)
  filename <- paste0("./link_peak_to_gene_plots/BDNF_alt_links_with_gex_p0.05/", 
                     name, "_link_peak_to_BDNF_plot.pdf")
  pdf(filename, width = 8, height = 4)
  p <- CoveragePlot(nmglut_obj_list[[i]], region = "BDNF", 
                    expression.assay = "SCT", 
                    extend.downstream = 100000, 
                    extend.upstream = 20000, 
                    links = T,
                    region.highlight = GRanges(seqnames = "chr11", 
                                               strand = "+",
                                               ranges = IRanges(start = 27763500, 
                                                                end = 27777500))) +
    theme(text = element_text(size = 12))
  print(p)
  dev.off()
  
  filename <- paste0("./link_peak_to_gene_plots/BDNF_alt_links_with_gex_p0.05/", 
                     name, "exprs_plot.pdf")
  pdf(filename, width = 1, height = 2)
  expr_plot <- ExpressionPlot(
    object = nmglut_obj_list[[i]],
    features = c("BDNF"),
    assay = "SCT"
  )
  print(expr_plot)
  dev.off()
}


# npglut ####
npglut_0hr <- subset(multiomic_obj, timextype.ident == "NEFM_pos_glut_0hr")
npglut_1hr <- subset(multiomic_obj, timextype.ident == "NEFM_pos_glut_1hr")
npglut_6hr <- subset(multiomic_obj, timextype.ident == "NEFM_pos_glut_6hr")
npglut_df_list <- vector(mode = "list", 3L)
npglut_obj_list <- list(npglut_0hr, npglut_1hr, npglut_6hr)
for (i in 1:length(npglut_df_list)){
  RNA <- npglut_obj_list[[i]]@assays$SCT@data[rownames(npglut_obj_list[[i]]@assays$SCT) == "BDNF", ]
  ATAC <- npglut_obj_list[[i]]@assays$peaks@data[str_detect(rownames(npglut_obj_list[[i]]@assays$peaks), 
                                                            "chr11-2[7|8][0-9]{6}"), ]
  df <- data.frame(p_value = rep_len(0, nrow(ATAC)),
                   cor = rep_len(0, nrow(ATAC)),
                   peak = rownames(ATAC),
                   zscore = rep_len(0, nrow(ATAC)))
  
  for (j in 1:nrow(ATAC)) {
    res <- cor.test(log1p(RNA), log1p(ATAC[j, ]))
    df$p_value[j] <- res$p.value
    df$cor[j] <- res$estimate
    df$zscore[j] <- res$statistic
  }
  
  df <- df[!is.na(df$p_value), ]
  df$q_value <- p.adjust(df$p_value, method = "BH")
  #df <- df[df$q_value < 0.05, ]
  #df <- df[df$q_value < 0.1, ]
  df <- df[df$p_value < 0.05, ]
  npglut_df_list[[i]] <- df
  es <- (as.numeric(str_remove_all(str_extract(df$peak, "-[0-9]+-"), "-")) + 
           as.numeric(str_remove(str_extract(df$peak, "[0-9]+$"), "-"))) / 2
  ss <- rep_len(27721989, nrow(df))
  for (k in 1:length(es)){
    e <- es[k]
    s <- ss[k]
    print(paste0("start: ", s, "end: ", e))
    if (e < s) {
      ss[k] <- e
      es[k] <- 27721989
      print("switched")
    }
  }
  gr <- GRanges(seqnames = str_extract(df$peak, "chr11"), 
                strand = "+",
                ranges = IRanges(start = ss,
                                 end = es))
  elementMetadata(gr)[["score"]] <- df$cor
  #elementMetadata(gr)[["pvalue"]] <- df$q_value
  elementMetadata(gr)[["pvalue"]] <- df$p_value  
  elementMetadata(gr)[["gene"]] <- "BDNF"
  elementMetadata(gr)[["peak"]] <- df$peak
  elementMetadata(gr)[["zscore"]] <- df$zscore
  Links(npglut_obj_list[[i]]) <- gr
  Annotation(npglut_obj_list[[i]]) <- annotation_ranges
  
}


for (i in 1:length(npglut_obj_list)){
  npglut_obj_list[[i]]$timextype.ident.for.plot <- str_replace_all(npglut_obj_list[[i]]$timextype.ident,
                                                                   "_", " ")
  npglut_obj_list[[i]]$timextype.ident.for.plot <- str_replace(npglut_obj_list[[i]]$timextype.ident.for.plot,
                                                               " pos", "+")
  
  Idents(npglut_obj_list[[i]]) <- "timextype.ident.for.plot"
  name <- unique(npglut_obj_list[[i]]$timextype.ident)
  filename <- paste0("./link_peak_to_gene_plots/BDNF_alt_links_with_gex_p0.05/", 
                     name, "_link_peak_to_BDNF_plot.pdf")
  pdf(filename, width = 8, height = 4)
  p <- CoveragePlot(npglut_obj_list[[i]], region = "BDNF", 
                    expression.assay = "SCT", 
                    extend.downstream = 100000, 
                    extend.upstream = 20000, 
                    links = T,
                    region.highlight = GRanges(seqnames = "chr11", 
                                               strand = "+",
                                               ranges = IRanges(start = 27763500, 
                                                                end = 27777500))) +
    theme(text = element_text(size = 12))
  print(p)
  dev.off()
  
  filename <- paste0("./link_peak_to_gene_plots/BDNF_alt_links_with_gex_p0.05/", 
                     name, "exprs_plot.pdf")
  pdf(filename, width = 1, height = 2)
  expr_plot <- ExpressionPlot(
    object = npglut_obj_list[[i]],
    features = c("BDNF"),
    assay = "SCT"
  )
  print(expr_plot)
  dev.off()
}

# GABA ####
GABA_0hr <- subset(multiomic_obj, timextype.ident == "GABA_0hr")
GABA_1hr <- subset(multiomic_obj, timextype.ident == "GABA_1hr")
GABA_6hr <- subset(multiomic_obj, timextype.ident == "GABA_6hr")
GABA_df_list <- vector(mode = "list", 3L)
GABA_obj_list <- list(GABA_0hr, GABA_1hr, GABA_6hr)
for (i in 1:length(GABA_df_list)){
  RNA <- GABA_obj_list[[i]]@assays$SCT@data[rownames(GABA_obj_list[[i]]@assays$SCT) == "BDNF", ]
  ATAC <- GABA_obj_list[[i]]@assays$peaks@data[str_detect(rownames(GABA_obj_list[[i]]@assays$peaks), 
                                                            "chr11-2[7|8][0-9]{6}"), ]
  df <- data.frame(p_value = rep_len(0, nrow(ATAC)),
                   cor = rep_len(0, nrow(ATAC)),
                   peak = rownames(ATAC),
                   zscore = rep_len(0, nrow(ATAC)))
  
  for (j in 1:nrow(ATAC)) {
    res <- cor.test(log1p(RNA), log1p(ATAC[j, ]))
    df$p_value[j] <- res$p.value
    df$cor[j] <- res$estimate
    df$zscore[j] <- res$statistic
  }
  
  df <- df[!is.na(df$p_value), ]
  df$q_value <- p.adjust(df$p_value, method = "BH")
  #df <- df[df$q_value < 0.05, ]
  #df <- df[df$q_value < 0.1, ]
  df <- df[df$p_value < 0.05, ]
  GABA_df_list[[i]] <- df
  es <- (as.numeric(str_remove_all(str_extract(df$peak, "-[0-9]+-"), "-")) + 
           as.numeric(str_remove(str_extract(df$peak, "[0-9]+$"), "-"))) / 2
  ss <- rep_len(27721989, nrow(df))
  for (k in 1:length(es)){
    e <- es[k]
    s <- ss[k]
    print(paste0("start: ", s, "end: ", e))
    if (e < s) {
      ss[k] <- e
      es[k] <- 27721989
      print("switched")
    }
  }
  gr <- GRanges(seqnames = str_extract(df$peak, "chr11"), 
                strand = "+",
                ranges = IRanges(start = ss,
                                 end = es))
  elementMetadata(gr)[["score"]] <- df$cor
  #elementMetadata(gr)[["pvalue"]] <- df$q_value
  elementMetadata(gr)[["pvalue"]] <- df$p_value
  elementMetadata(gr)[["gene"]] <- "BDNF"
  elementMetadata(gr)[["peak"]] <- df$peak
  elementMetadata(gr)[["zscore"]] <- df$zscore
  Links(GABA_obj_list[[i]]) <- gr
  Annotation(GABA_obj_list[[i]]) <- annotation_ranges
  
}


for (i in 1:length(GABA_obj_list)){
  GABA_obj_list[[i]]$timextype.ident.for.plot <- str_replace_all(GABA_obj_list[[i]]$timextype.ident,
                                                                   "_", " ")
  GABA_obj_list[[i]]$timextype.ident.for.plot <- str_replace(GABA_obj_list[[i]]$timextype.ident.for.plot,
                                                               " pos", "+")
  
  Idents(GABA_obj_list[[i]]) <- "timextype.ident.for.plot"
  name <- unique(GABA_obj_list[[i]]$timextype.ident)
  filename <- paste0("./link_peak_to_gene_plots/BDNF_alt_links_with_gex_p0.05/", 
                     name, "_link_peak_to_BDNF_plot.pdf")
  
  pdf(filename, width = 8, height = 4)
  p <- CoveragePlot(GABA_obj_list[[i]], region = "BDNF", 
                    expression.assay = "SCT", 
                    extend.downstream = 100000, 
                    extend.upstream = 20000, 
                    links = T,
                    region.highlight = GRanges(seqnames = "chr11", 
                                               strand = "+",
                                               ranges = IRanges(start = 27763500, 
                                                                end = 27777500))) +
    theme(text = element_text(size = 12))
  print(p)
  dev.off()
  
  filename <- paste0("./link_peak_to_gene_plots/BDNF_alt_links_with_gex_p0.05/", 
                     name, "exprs_plot.pdf")
  pdf(filename, width = 1, height = 2)
  expr_plot <- ExpressionPlot(
    object = GABA_obj_list[[i]],
    features = c("BDNF"),
    assay = "SCT"
  )
  print(expr_plot)
  dev.off()
}


# remove lines 02, 05, 06 and calculate again ####
unique(multiomic_obj$cell.line.ident)
removed_lines <- subset(multiomic_obj, 
                        cell.line.ident %in% c("CD_62", "CD_50", "CD_56", "CD_34",
                                               "CD_35", "CD_64", "CD_51", "CD_33",
                                               "CD_66", "CD_63", "CD_65", "CD_60",
                                               "CD_55", "CD_04", "CD_52"))

# separate object by cell line and type ####
lines <-unique(multiomic_obj$cell.line.ident)
types <- unique(multiomic_obj$cell.type)
times <- unique(multiomic_obj$time.ident)
for (line in lines){
  for (type in types){
    for (time in times){
      print(paste(line, type, time, sep = "; "))
      tempobj <- subset(multiomic_obj, cell.line.ident == line)
      tempobj <- subset(tempobj, cell.type == type)
      tempobj <- subset(tempobj, time.ident == time)
      
      Annotation(tempobj) <- annotation_ranges
      
      RNA <- tempobj@assays$SCT@data[rownames(tempobj@assays$SCT) == "BDNF", ]
      ATAC <- tempobj@assays$peaks@data[str_detect(rownames(tempobj@assays$peaks), 
                                                   "chr11-2[7|8][0-9]{6}"), ]
      df <- data.frame(p_value = rep_len(0, nrow(ATAC)),
                       cor = rep_len(0, nrow(ATAC)),
                       peak = rownames(ATAC),
                       zscore = rep_len(0, nrow(ATAC)))
      
      for (j in 1:nrow(ATAC)) {
        res <- cor.test(log1p(RNA), log1p(ATAC[j, ]))
        df$p_value[j] <- res$p.value
        df$cor[j] <- res$estimate
        df$zscore[j] <- res$statistic
      }
      
      df <- df[!is.na(df$p_value), ]
      df$q_value <- p.adjust(df$p_value, method = "BH")
      #df <- df[df$q_value < 0.05, ]
      #df <- df[df$q_value < 0.1, ]
      df <- df[df$p_value < 0.05, ]
      es <- (as.numeric(str_remove_all(str_extract(df$peak, "-[0-9]+-"), "-")) + 
               as.numeric(str_remove(str_extract(df$peak, "[0-9]+$"), "-"))) / 2
      ss <- rep_len(27721989, nrow(df))
      for (k in 1:length(es)){
        e <- es[k]
        s <- ss[k]
        print(paste0("start: ", s, "end: ", e))
        if (e < s) {
          ss[k] <- e
          es[k] <- 27721989
          print("switched")
        }
      }
      gr <- GRanges(seqnames = str_extract(df$peak, "chr11"), 
                    strand = "+",
                    ranges = IRanges(start = ss,
                                     end = es))
      elementMetadata(gr)[["score"]] <- df$cor
      #elementMetadata(gr)[["pvalue"]] <- df$q_value
      elementMetadata(gr)[["pvalue"]] <- df$p_value
      elementMetadata(gr)[["gene"]] <- "BDNF"
      elementMetadata(gr)[["peak"]] <- df$peak
      elementMetadata(gr)[["zscore"]] <- df$zscore
      
      Links(tempobj) <- gr
      
      tempobj$timextype.ident.for.plot <- str_replace_all(tempobj$timextype.ident,
                                                          "_", " ")
      tempobj$timextype.ident.for.plot <- str_replace(tempobj$timextype.ident.for.plot,
                                                      " pos", "+")
      
      tempobj$timextype.ident.for.plot <- str_replace(tempobj$timextype.ident.for.plot,
                                                      " neg", "-")
      
      Idents(tempobj) <- "timextype.ident.for.plot"
      filename <- paste0("./link_peak_to_gene_plots/BDNF_alt_links_with_gex_p0.05/", 
                         paste(line, type, time, sep = "_"), "_link_peak_to_BDNF_plot.pdf")
      
      print(filename)
      pdf(filename, width = 10, height = 5)
      p <- CoveragePlot(tempobj, region = "BDNF", 
                        expression.assay = "SCT", 
                        extend.downstream = 100000, 
                        extend.upstream = 20000, 
                        links = T,
                        region.highlight = GRanges(seqnames = "chr11", 
                                                   strand = "+",
                                                   ranges = IRanges(start = 27763500, 
                                                                    end = 27777500))) +
        theme(text = element_text(size = 12))
      print(p)
      dev.off()
    }
  }
}
