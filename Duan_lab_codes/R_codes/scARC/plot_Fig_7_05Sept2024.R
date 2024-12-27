# Siwei 5 Sept 2024
# process all 100 (actually ~ 94-97 lines depending on cell type) lines
# Remake the plots in pdf for publication
# Dec 2023

# init
{
  library(readr)
  library(vcfR)
  library(stringr)
  library(ggplot2)
  library(scales)
  
  library(ggstatsplot)
  library(ggside)
  
  library(parallel)
  library(future)
  
  library(MASS)
  
  library(RColorBrewer)
  library(grDevices)
  
  library(ggvenn)
  
  library(gridExtra)
  library(cowplot)
}

# param
plan("multisession", workers = 4)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)
set.seed(42)

load("~/NVME/scARC_Duan_018/R_ASoC/Sum_100_lines_21Dec2023.RData")


# print(master_vcf_list[[1]][[2]])
# 
# packageurl <- "http://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.4.tar.gz"
# install.packages(packageurl,
#                  repos = NULL,
#                  type = "source")

# Need to regenerate plot
# recalc vcr of 100 lines #####

func_calc_vcr <-
  function(vcr_readin) {
    return_list <-
      vector(mode = "list",
             length = 2L)
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
    
    
    return_list[[2]] <- 
    ggplot(df_assembled,
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
    # png(filename = paste0("plot_",
    #                       names(vcr_readin),
    #                       '.png'),
    #     width = 800,
    #     height = 800, 
        # units = "px")
    # print(return_list[[2]])
    # dev.off()
    
    # df_assembled_FDR_005 <-
    #   df_assembled[df_assembled$FDR < 0.05, ]
    
    # df_avinput <-
    #   df_assembled_FDR_005[, c(1,2,2,4,5,3,6:ncol(df_assembled_FDR_005))]
    # write.table(df_avinput,
    #             file = "DN_39_lines_005_21Aug2023.avinput",
    #             quote = F, sep = "\t",
    #             row.names = F, col.names = T)
    # return_list[[3]] <-
    #   df_avinput
    
    return(return_list)
  }

# Run #####

master_plot_list <-
  vector(mode = "list",
         length = length(df_vcf_list))
names(master_plot_list) <-
  sample_names

master_plot_df_list <-
  vector(mode = "list",
         length = length(df_vcf_list))
names(master_plot_df_list) <-
  sample_names

for (i in 1:length(df_vcf_list)) {
  print(paste0("The current sample is ",
               df_vcf_list[i]))
  returned_list <-
    func_calc_vcr(paste0("./merge_60_lines_024_025_029_25Sept2023/output/",
                         df_vcf_list[i]))
  master_plot_df_list[[i]] <-
    returned_list[[1]]
  master_plot_list[[i]] <-
    returned_list[[2]]
}

save(list = c("master_plot_list",
              "master_plot_df_list"),
     file = "regenerated_volcano_w_ggplot_06Sept2024.RData")

func_make_volcano_plot_only <-
  function(df_assembled) {
    returned_plot <-
      ggplot(df_assembled,
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
  }


print(master_plot_list[[2]])

for (i in 1:length(master_plot_df_list)) {
  print(i)
  df_assembled <-
    master_plot_df_list[[i]]
  
  vcr_readin <-
    paste0("./merge_60_lines_024_025_029_25Sept2023/output/",
           df_vcf_list[i])
  df_raw <-
    read.vcfR(file = vcr_readin,
              limit = 2e9,
              verbose = T)
  
  sum_gt_AD <-
    extract.gt(x = df_raw,
               element = "AD",
               as.numeric = F,
               convertNA = T)
  line_count <- ncol(sum_gt_AD)
  
  master_plot_list[[i]] <-
    ggplot(df_assembled,
           aes(x = (REF_N/DP),
               y = (0 - log10(pVal)))) +
    geom_point(size = 0.2,
               aes(colour = factor(ifelse(FDR < 0.05,
                                          "FDR < 0.05",
                                          "FDR > 0.05")))) +
    scale_color_manual(values = c("red", "black")) +
    labs(colour = "FDR value") +
    ggtitle(paste(names(master_plot_df_list)[i], 
                  # " DP >= 30, minAllele >= 2,",
                  "\n",
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
    theme_classic()  +
    theme(legend.position = "none",
          plot.title = element_text(size = 8))
}

print(master_plot_list[[2]])

plot_grid(master_plot_list[[1]], 
          master_plot_list[[2]],
          master_plot_list[[3]],
          master_plot_list[[4]],
          master_plot_list[[5]],
          master_plot_list[[6]],
          master_plot_list[[7]],
          master_plot_list[[8]],
          master_plot_list[[9]],
          nrow = 3)

save(list = c("master_plot_list",
              "master_plot_df_list"),
     file = "regenerated_volcano_w_ggplot_06Sept2024.RData")
