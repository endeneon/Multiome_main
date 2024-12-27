# Siwei 25 Sept 2023
# process Duan-024 (18), Duan-025 (20), Duan-029 (16+8) lines
# Sept 2023

# init
{
  library(readr)
  library(vcfR)
  library(stringr)
  library(ggplot2)
  
  library(parallel)
  library(future)
  
  library(MASS)
  
  library(RColorBrewer)
  library(grDevices)
  
  library(ggvenn)
}

# param
plan("multisession", workers = 12)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)
set.seed(42)

# load the table from vcf (the vcf has been prefiltered to include DP >= 20 only)
df_vcf_list <-
  dir(path = ".",
      pattern = "*25Sept2023\\.vcf$",
      recursive = T)
df_vcf_list

sample_names <-
  str_split(string = df_vcf_list,
            pattern = "scARC_60_lines_merge_all_VQSR_",
            simplify = T)[, 2]
sample_names <-
  str_remove(string = sample_names,
             pattern = "_25Sept2023\\.vcf$")

master_vcf_list <-
  vector(mode = "list",
         length = length(sample_names))
names(master_vcf_list) <-
  sample_names

names(df_vcf_list) <-
  sample_names
# names(master_vcf_list) <- sample_names

# recalc vcr of 74 lines

func_calc_vcr <-
  function(vcr_readin) {
    return_list <-
      vector(mode = "list",
             length = 3L)
    print(names(vcr_readin))
    
    df_raw <-
      read.vcfR(file = vcr_readin,
                limit = 2e9,
                verbose = T)
    
    cell_type <- names(vcr_readin)
    # str_split(string = names(vcr_readin),
    #           pattern = "_",
    #           simplify = T)[, 2]
    
    df_fix <-
      as.data.frame(df_raw@fix)
    
    sum_gt_AD <-
      extract.gt(x = df_raw,
                 element = "AD",
                 as.numeric = F,
                 convertNA = T)
    line_count <- ncol(sum_gt_AD)
    
    ncalled_counts <-
      extract.gt(x = df_raw,
                 element = "AD",
                 as.numeric = F,
                 convertNA = T)
    ncalled_counts <-
      str_split(string = ncalled_counts,
                pattern = ",",
                simplify = T)[, 1]
    ncalled_counts <- as.numeric(ncalled_counts)
    ncalled_counts <-
      matrix(data = ncalled_counts,
             ncol = ncol(sum_gt_AD))
    ncalled_counts <-
      as.data.frame(ncalled_counts)
    ncalled_counts <-
      rowSums(!is.na(ncalled_counts))
    
    sum_gt_AD[is.na(sum_gt_AD)] <- "0,0"
    
    sum_REF <-
      str_split(string = sum_gt_AD,
                pattern = ",",
                simplify = T)[, 1]
    sum_REF <- as.numeric(sum_REF)
    sum_REF <-
      matrix(data = sum_REF,
             ncol = ncol(sum_gt_AD))
    sum_REF <-
      rowSums(sum_REF,
              na.rm = T)
    
    sum_ALT <-
      str_split(string = sum_gt_AD,
                pattern = ",",
                simplify = T)[, 2]
    sum_ALT <- as.numeric(sum_ALT)
    sum_ALT <-
      matrix(data = sum_ALT,
             ncol = ncol(sum_gt_AD))
    sum_ALT <-
      rowSums(sum_ALT,
              na.rm = T)
    
    
    
    df_assembled <-
      as.data.frame(cbind(df_fix,
                          data.frame(NCALLED = ncalled_counts),
                          data.frame(REF_N = sum_REF),
                          data.frame(ALT_N = sum_ALT)))
    
    df_assembled <-
      df_assembled[str_detect(string = df_assembled$ID,
                              pattern = "^rs.*"),
      ]
    df_assembled <-
      df_assembled[(df_assembled$REF_N > 2 & df_assembled$ALT_N > 2), ]
    df_assembled$DP <-
      (df_assembled$REF_N + df_assembled$ALT_N)
    df_assembled <-
      df_assembled[df_assembled$DP > 29, ]
    df_assembled <-
      df_assembled[!is.na(df_assembled$DP), ]
    # 
    
    print(ncol(df_assembled)) # make sure column counts match df_assembled$pVal calculation
    print(nrow(df_assembled)) # check remaining SNP counts
    # print(head(df_assembled))
    # break
    clust_df_assembled <- makeCluster(type = "FORK", 20)
    clusterExport(cl = clust_df_assembled, 
                  varlist = "df_assembled",
                  envir = environment())
    df_assembled$pVal <-
      parApply(cl = clust_df_assembled,
               X = df_assembled,
               MARGIN = 1,
               FUN = function(x)(binom.test(x = as.numeric(x[10]),
                                            n = as.numeric(x[12]),
                                            p = 0.5,
                                            alternative = "t")$p.value))
    stopCluster(clust_df_assembled)
    rm(clust_df_assembled)
    
    df_assembled$FDR <-
      p.adjust(p = df_assembled$pVal,
               method = "fdr")
    df_assembled <-
      df_assembled[order(df_assembled$pVal), ]
    
    df_assembled$`-logFDR` <-
      0 - log10(df_assembled$FDR)
    df_assembled$REF_ratio <-
      df_assembled$REF_N / df_assembled$DP
    df_assembled <-
      df_assembled[order(df_assembled$pVal), ]
    
    return_list[[1]] <- df_assembled
    
    
    return_list[[2]] <- ggplot(df_assembled,
                               aes(x = (REF_N/DP),
                                   y = (0 - log10(pVal)))) +
      geom_point(size = 0.2,
                 aes(colour = factor(ifelse(FDR < 0.05,
                                            "FDR < 0.05",
                                            "FDR > 0.05")))) +
      scale_color_manual(values = c("red", "black")) +
      labs(colour = "FDR value") +
      ggtitle(paste(names(vcr_readin), 
                    " DP >= 30, minAllele >= 2,",
                    line_count, "lines;\n",
                    # "SNP lists from 23 Sept results, removed all NCALLED >= 14\n",
                    nrow(df_assembled),
                    "SNPs, of",
                    sum(df_assembled$FDR < 0.05),
                    "FDR < 0.05;\n",
                    "FDR < 0.05 & REF/DP < 0.5 =",
                    sum((df_assembled$REF_N/df_assembled$DP < 0.5) &
                          (df_assembled$FDR < 0.05)),
                    ";\n FDR < 0.05 & REF/DP > 0.5 =",
                    sum((df_assembled$REF_N/df_assembled$DP > 0.5) &
                          (df_assembled$FDR < 0.05))
                    # , ";\n",
                    # "rs2027349:", 
                    # "REF =", df_assembled[df_assembled$ID == "rs2027349", ]$REF_N,
                    # ", ALT =", df_assembled[df_assembled$ID == "rs2027349", ]$ALT_N,
                    # ", DP =", df_assembled[df_assembled$ID == "rs2027349", ]$DP,
                    # " P =", formatC(df_assembled[df_assembled$ID == "rs2027349", ]$pVal,
                    #                 digits = 3),
                    # ", FDR =", formatC(df_assembled[df_assembled$ID == "rs2027349", ]$FDR,
                    #                    digits = 3)
              )) +
      ylab("-log10P") +
      ylim(0, 100) +
      theme_classic() 
    png(filename = paste0("plot_",
                          names(vcr_readin),
                          '.png'),
        width = 800,
        height = 800, 
        units = "px")
    print(return_list[[2]])
    dev.off()
    
    df_assembled_FDR_005 <-
      df_assembled[df_assembled$FDR < 0.05, ]
    
    df_avinput <-
      df_assembled_FDR_005[, c(1,2,2,4,5,3,6:ncol(df_assembled_FDR_005))]
    # write.table(df_avinput,
    #             file = "DN_39_lines_005_21Aug2023.avinput",
    #             quote = F, sep = "\t",
    #             row.names = F, col.names = T)
    return_list[[3]] <-
      df_avinput
    
    return(return_list)
  }

## Run #####
for (i in 1:length(df_vcf_list)) {
  print(paste0("sample is ",
               df_vcf_list[i]))
  master_vcf_list[[i]] <-
    func_calc_vcr(df_vcf_list[i])
}


# this function count the depths of each vcf file,
# of both gross and ones passed current filter 
# (DP >= 30, minREF/ALT > 2)

function_calc_depth <-
  function(vcr_readin) {
    return_list <-
      vector(mode = "list",
             length = 2L)
    names(return_list) <-
      c("gross_depth",
        "filtered_depth")
    
    print(names(vcr_readin))
    
    df_raw <-
      read.vcfR(file = vcr_readin,
                limit = 2e9,
                verbose = T)
    
    cell_type <- names(vcr_readin)
    # str_split(string = names(vcr_readin),
    #           pattern = "_",
    #           simplify = T)[, 2]
    df_fix <-
      as.data.frame(df_raw@fix)
    
    sum_gt_AD <-
      extract.gt(x = df_raw,
                 element = "AD",
                 as.numeric = F,
                 convertNA = T)
    line_count <- ncol(sum_gt_AD)
    
    ncalled_counts <-
      extract.gt(x = df_raw,
                 element = "AD",
                 as.numeric = F,
                 convertNA = T)
    ncalled_counts <-
      str_split(string = ncalled_counts,
                pattern = ",",
                simplify = T)[, 1]
    ncalled_counts <- as.numeric(ncalled_counts)
    ncalled_counts <-
      matrix(data = ncalled_counts,
             ncol = ncol(sum_gt_AD))
    ncalled_counts <-
      as.data.frame(ncalled_counts)
    ncalled_counts <-
      rowSums(!is.na(ncalled_counts))
    
    sum_gt_AD[is.na(sum_gt_AD)] <- "0,0"
    
    sum_REF <-
      str_split(string = sum_gt_AD,
                pattern = ",",
                simplify = T)[, 1]
    sum_REF <- as.numeric(sum_REF)
    sum_REF <-
      matrix(data = sum_REF,
             ncol = ncol(sum_gt_AD))
    sum_REF <-
      rowSums(sum_REF,
              na.rm = T)
    
    sum_ALT <-
      str_split(string = sum_gt_AD,
                pattern = ",",
                simplify = T)[, 2]
    sum_ALT <- as.numeric(sum_ALT)
    sum_ALT <-
      matrix(data = sum_ALT,
             ncol = ncol(sum_gt_AD))
    sum_ALT <-
      rowSums(sum_ALT,
              na.rm = T)
    
    return_list[[1]] <-
      sum(sum_REF, 
          sum_ALT)
    
    df_assembled <-
      as.data.frame(cbind(df_fix,
                          data.frame(NCALLED = ncalled_counts),
                          data.frame(REF_N = sum_REF),
                          data.frame(ALT_N = sum_ALT)))
    
    df_assembled <-
      df_assembled[str_detect(string = df_assembled$ID,
                              pattern = "^rs.*"),
      ]
    df_assembled <-
      df_assembled[(df_assembled$REF_N > 2 & df_assembled$ALT_N > 2), ]
    df_assembled$DP <-
      (df_assembled$REF_N + df_assembled$ALT_N)
    df_assembled <-
      df_assembled[df_assembled$DP > 29, ]
    df_assembled <-
      df_assembled[!is.na(df_assembled$DP), ]
    
    return_list[[2]] <-
      sum(df_assembled$DP)
    
    return(return_list)
  }

for (i in 1:length(df_vcf_list)) {
  print(paste0("the current file is ",
               df_vcf_list[i]))
  results_from_depth_calc <-
    function_calc_depth(vcr_readin = df_vcf_list[i])
  line_to_write <-
    paste("the gross depth of",
          names(df_vcf_list)[i],
          "is",
          results_from_depth_calc[[1]],
          ";\nthe sum depth of SNPs passed filter is",
          results_from_depth_calc[[2]])
  cat(line_to_write,
      file = "Duan_018_024_025_029_76_lines_depth.txt",
      sep = "\t",
      fill = F,
      append = T)
}

# re-run 6hr
for (i in 7:length(df_vcf_list)) {
  print(paste0("sample is ",
               df_vcf_list[i]))
  master_vcf_list[[i]] <-
    func_calc_vcr(df_vcf_list[i])
}


saveRDS(object = master_vcf_list,
        file = "Duan_Project_018_024_025_029_ASoC_74_lines_11Oct2023.RDs")
 # save.image("Duan_Project_025_all_ASoC_12Sept2023.RData")

## make venn plot of each cell type #####
# > sample_names
# [1] "0hr_GABA"   "0hr_nmglut" "0hr_npglut" "1hr_GABA"   "1hr_nmglut" "1hr_npglut"
# [7] "6hr_GABA"   "6hr_nmglut" "6hr_npglut"

### GABA #####
list_venn_2_plot <-
  vector(mode = "list",
         length = 3L)
names(list_venn_2_plot) <-
  c("hr_1",
    "hr_6",
    "hr_0")

i <- 1
list_venn_2_plot[[1]] <-
  master_vcf_list[[i + 3]][[3]]$ID
list_venn_2_plot[[2]] <-
  master_vcf_list[[i + 6]][[3]]$ID
list_venn_2_plot[[3]] <-
  master_vcf_list[[i]][[3]]$ID

ggvenn(data = list_venn_2_plot,
       fill_color = brewer.pal(n = 3, 
                               name = "Dark2"),
       stroke_size = 0.5,
       set_name_size = 6,
       text_size = 4) +
  ggtitle(paste("GABA,\ntotal FDR < 0.05 SNP count (unique) =",
                length(unique(c(list_venn_2_plot[[1]],
                                list_venn_2_plot[[2]],
                                list_venn_2_plot[[3]])))))

i <- 2
list_venn_2_plot[[1]] <-
  master_vcf_list[[i + 3]][[3]]$ID
list_venn_2_plot[[2]] <-
  master_vcf_list[[i + 6]][[3]]$ID
list_venn_2_plot[[3]] <-
  master_vcf_list[[i]][[3]]$ID

ggvenn(data = list_venn_2_plot,
       fill_color = brewer.pal(n = 3, 
                               name = "Dark2"),
       stroke_size = 0.5,
       set_name_size = 6,
       text_size = 4) +
  ggtitle(paste("NEFM- Glut,\ntotal FDR < 0.05 SNP count (unique) =",
                length(unique(c(list_venn_2_plot[[1]],
                                list_venn_2_plot[[2]],
                                list_venn_2_plot[[3]])))))


i <- 3
list_venn_2_plot[[1]] <-
  master_vcf_list[[i + 3]][[3]]$ID
list_venn_2_plot[[2]] <-
  master_vcf_list[[i + 6]][[3]]$ID
list_venn_2_plot[[3]] <-
  master_vcf_list[[i]][[3]]$ID

ggvenn(data = list_venn_2_plot,
       fill_color = brewer.pal(n = 3, 
                               name = "Dark2"),
       stroke_size = 0.5,
       set_name_size = 6,
       text_size = 4) +
  ggtitle(paste("NEFM+ Glut,\ntotal FDR < 0.05 SNP count (unique) =",
                length(unique(c(list_venn_2_plot[[1]],
                                list_venn_2_plot[[2]],
                                list_venn_2_plot[[3]])))))


## make volcano plots from master vcf lists
for (i in 1:length(master_vcf_list)) {
  df_2_plot <- 
    master_vcf_list[[i]][[1]]
  volcano_plot_to_print <-
    ggplot(df_2_plot,
           aes(x = (REF_N/DP),
             y = (0 - log10(pVal)))) +
    geom_point(size = 0.2,
               aes(colour = factor(ifelse(FDR < 0.05,
                                          "FDR < 0.05",
                                          "FDR > 0.05")))) +
    scale_color_manual(values = c("red", "black")) +
    labs(colour = "FDR value") +
    ggtitle(paste(names(master_vcf_list)[i], 
                  " DP >= 30, minAllele >= 2,\n",
                  # line_count, "lines;\n",
                  # "SNP lists from 23 Sept results, removed all NCALLED >= 14\n",
                  nrow(df_2_plot),
                  "SNPs, of",
                  sum(df_2_plot$FDR < 0.05),
                  "FDR < 0.05;\n",
                  "FDR < 0.05 & REF/DP < 0.5 =",
                  sum((df_2_plot$REF_N/df_2_plot$DP < 0.5) &
                        (df_2_plot$FDR < 0.05)),
                  ";\n FDR < 0.05 & REF/DP > 0.5 =",
                  sum((df_2_plot$REF_N/df_2_plot$DP > 0.5) &
                        (df_2_plot$FDR < 0.05))
    )) +
    ylab("-log10P") +
    ylim(0, 100) +
    theme_classic() 

  png(filename = paste0("plot_",
                          names(master_vcf_list)[i],
                          '.png'),
      width = 400,
      height = 400,
      units = "px")
  print(volcano_plot_to_print)
  dev.off()
}

dir.create("results_76lines_17Oct2023")
## write out summary tables
for (i in 1:length(master_vcf_list)) {
  print(names(master_vcf_list)[i])
  df_2_writeout <- 
    master_vcf_list[[i]][[1]]
  df_2_writeout <-
    df_2_writeout[order(df_2_writeout$pVal), ]
  df_2_writeout$INFO <- NULL
  df_2_writeout$FILTER <- NULL
  df_2_writeout$`-logFDR` <- NULL
  write.table(df_2_writeout,
              file = paste0("results_76lines_17Oct2023/",
                            names(master_vcf_list)[i],
                            "_SNP_2_30_full.tsv"),
              quote = F, sep = "\t",
              row.names = F, col.names = T)
  
}

df_2_plot <- 
  master_vcf_list[[1]][[1]]
ggplot(df_2_plot,
       aes(x = (REF_N/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.2,
             aes(colour = factor(ifelse(FDR < 0.05,
                                        "FDR < 0.05",
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle(paste(names(master_vcf_list)[i], 
                " DP >= 30, minAllele >= 2,",
                # line_count, "lines;\n",
                # "SNP lists from 23 Sept results, removed all NCALLED >= 14\n",
                nrow(df_2_plot),
                "SNPs, of",
                sum(df_2_plot$FDR < 0.05),
                "FDR < 0.05;\n",
                "FDR < 0.05 & REF/DP < 0.5 =",
                sum((df_2_plot$REF_N/df_2_plot$DP < 0.5) &
                      (df_2_plot$FDR < 0.05)),
                ";\n FDR < 0.05 & REF/DP > 0.5 =",
                sum((df_2_plot$REF_N/df_2_plot$DP > 0.5) &
                      (df_2_plot$FDR < 0.05))
  )) +
  ylab("-log10P") +
  ylim(0, 100) +
  theme_classic() 
