# Siwei 12 Sept 2023
# Process all Duan-024 ()
# Sept 2023

# init
{
  library(readr)
  library(vcfR)
  library(stringr)
  library(ggplot2)
  
  library(parallel)
  
  library(MASS)
  
  library(RColorBrewer)
  library(grDevices)
}

# load the table from vcf (the vcf has been prefiltered to include DP >= 20 only)
df_vcf_list <-
  dir(path = "no_VQSR_Duan_025_12Sept2023",
      pattern = "*.vcf")

sample_names <-
  str_remove(string = df_vcf_list,
             pattern = "Duan_Project_024_")
sample_names <-
  str_remove(string = sample_names,
             pattern = "_trimmed_merged_SNP_VF\\.vcf$")

master_vcf_list <-
  vector(mode = "list",
         length = length(sample_names))
names(master_vcf_list) <
  sample_names

names(df_vcf_list) <-
  sample_names



# make a function returns the all SNPs, FDR < 0.05 SNPs, and volcano plots

func_calc_vcr <-
  function(vcr_readin) {
    return_list <-
      vector(mode = "list",
             length = 3L)
    # print(names(vcr_readin))
    
    df_raw <-
      read.vcfR(file = paste0("no_VQSR_Duan_025_12Sept2023",
                             "/",
                             vcr_readin),
                limit = 2e9,
                verbose = T)
    
    cell_type <- 
      str_split(string = names(vcr_readin),
                pattern = "_",
                simplify = T)[, 3]
    
    df_fix <-
      as.data.frame(df_raw@fix)
    
    sum_gt_AD <-
      extract.gt(x = df_raw,
                 element = "AD",
                 as.numeric = F,
                 convertNA = T)

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
      df_assembled[df_assembled$DP > 19, ]
    df_assembled <-
      df_assembled[!is.na(df_assembled$DP), ]
    # 

    print(ncol(df_assembled))
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
      geom_point(size = 0.5,
                 aes(colour = factor(ifelse(FDR < 0.05,
                                            "FDR < 0.05",
                                            "FDR > 0.05")))) +
      scale_color_manual(values = c("red", "black")) +
      labs(colour = "FDR value") +
      ggtitle(paste(names(vcr_readin), 
                    " DP >= 20, minAllele >= 2, 20 lines;\n",
                    # "SNP lists from 23 Sept results, removed all NCALLED >= 14\n",
                    nrow(df_assembled),
                    "SNPs, of",
                    sum(df_assembled$FDR < 0.05),
                    "FDR < 0.05;\n",
                    "FDR < 0.05 & REF/DP < 0.5 = ",
                    sum((df_assembled$REF_N/df_assembled$DP < 0.5) &
                          (df_assembled$FDR < 0.05)),
                    ";\n FDR < 0.05 & REF/DP > 0.5 = ", "",
                    sum((df_assembled$REF_N/df_assembled$DP > 0.5) &
                          (df_assembled$FDR < 0.05)), ";\n")) +
      ylab("-log10P") +
      ylim(0, 100) +
      theme_classic()
    
    print(return_list[[2]])
    
    df_assembled_FDR_005 <-
      df_assembled[df_assembled$FDR < 0.05, ]
    

    # save.image("DN_39_lines_21Aug2023.RData")
    
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

#####
for (i in 1:length(df_vcf_list)) {
  print(paste0("sample is ",
               df_vcf_list[i]))
  master_vcf_list[[i]] <-
    func_calc_vcr(df_vcf_list[i])
}
names(master_vcf_list) <- sample_names
save.image("Duan_Project_025_all_ASoC_12Sept2023.RData")

