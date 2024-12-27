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
}

# param
plan("multisession", workers = 12)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

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
  line_to_write <-
    paste("the gross depth of",
          names(df_vcf_list)[i],
          "is",
          function_calc_depth[[1]],
          ";\nthe sum depth of SNPs passed filter is",
          function_calc_depth[[2]])
  cat(line_to_write,
      file = "Duan_024_025_029_70_lines_depth.txt",
      sep = "\t",
      fill = F,
      append = T)
}

# re-run 6hr
# for (i in 7:length(df_vcf_list)) {
#   print(paste0("sample is ",
#                df_vcf_list[i]))
#   master_vcf_list[[i]] <-
#     func_calc_vcr(df_vcf_list[i])
# }


saveRDS(object = master_vcf_list,
        file = "Duan_Project_024_025_029_ASoC_70_lines_04Oct2023.RDs")
 # save.image("Duan_Project_025_all_ASoC_12Sept2023.RData")

