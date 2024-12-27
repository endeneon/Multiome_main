# Siwei 03 Dec 2022
# Calculate the direct count output results from bcftools mpileup
# from 29 Dec 2022 as input
# use vcfR


# init
library(readr)
library(vcfR)
library(stringr)
library(ggplot2)

library(parallel)

library(MASS)

library(RColorBrewer)
library(grDevices)


### load data

# ## mbq=10, mmq=10, -B
# vcf_files <-
#   list.files(path = "~/NVME/scARC_Duan_024_raw_fastq/GRCh38_output/hr_0/known_genotypes_neg_pt1/vcf_output_direct_count_mgs_29Dec2022",
#              pattern = ".*vcf$",
#              full.names = T)

# ## mbq=10, mmq=10, NO -B
# vcf_files <-
#   list.files(path = "~/NVME/scARC_Duan_024_raw_fastq/GRCh38_output/hr_0/known_genotypes_neg_pt1/vcf_output_direct_count_mmq10_mbq10_no_B_05Jan2023",
#              pattern = ".*vcf$",
#              full.names = T)

## mbq=5, mmq=20, -B
vcf_files <-
  list.files(path = "~/NVME/scARC_Duan_024_raw_fastq/GRCh38_output/hr_0/known_genotypes_neg_pt1/vcf_output_direct_count_mm20_mbq5_B_06Jan2023",
             pattern = ".*vcf$",
             full.names = T)

# ## mbq=10, mmq=20, -B
# vcf_files <-
#   list.files(path = "~/NVME/scARC_Duan_024_raw_fastq/GRCh38_output/hr_0/known_genotypes_neg_pt1/vcf_output_direct_count_mmq_20_mbq_10_B_05Jan2023",
#              pattern = ".*vcf$",
#              full.names = T)

## dedup only, no WASP
# vcf_files <-
#   list.files(path = "~/NVME/scARC_Duan_024_raw_fastq/GRCh38_output/hr_0_dedup_only/known_genotypes_bwa_by_barcode_hg38_WASP_deduped_only_BAMs/vcf_output",
#              pattern = ".*vcf$",
#              full.names = T)

# ## mbq=1, mmq=20, -B
# vcf_files <-
#   list.files(path = "~/NVME/scARC_Duan_024_raw_fastq/GRCh38_output/hr_0/known_genotypes_neg_pt1/vcf_output_direct_count_mmq_20_mbq_0_B_13Jan2023",
#              pattern = ".*vcf$",
#              full.names = T)

## mbq=0, mmq=20, -B
vcf_files <-
  list.files(path = "~/NVME/scARC_Duan_024_raw_fastq/GRCh38_output/hr_0/known_genotypes_neg_pt1/vcf_output_direct_count_mmq_20_mbq_0_B_16Jan2023",
             pattern = ".*vcf$",
             full.names = T)

## mbq=1, mmq=20, -B
vcf_files <-
  list.files(path = "~/NVME/scARC_Duan_024_raw_fastq/GRCh38_output/hr_0/known_genotypes_neg_pt1/vcf_output_direct_count_mmq_20_mbq_1_B_16Jan2023",
             pattern = ".*vcf$",
             full.names = T)

## make a list to hold the contents of all vcf files
vcf_list <-
  vector(mode = "list",
         length = length(vcf_files))

## readin vcf data one by one
## if used within a loop (local variable), the loop counter (i) does not need
## to be pre-defined
for (i in 1:length(vcf_list)) {
  print(i)
  vcf_list[[i]] <-
    read.vcfR(file = vcf_files[i],
              limit = 5e+8,
              verbose = T)
}

vcf_names <-
  str_split(string = vcf_files,
            pattern = "/",
            simplify = T)[, 10]

# vcf_names <-
#   str_remove(string = vcf_names,
#              pattern = "_gluthg38only_bwa_barcode_known_genotypes_WASPed_28Dec_direct_counting.vcf")
vcf_names <-
  str_remove(string = vcf_names,
             pattern = "_gluthg38only_bwa_barcode_known_genotypes_WASPed_28Dec_direct_counting.vcf")
# vcf_names <-
#   str_remove(string = vcf_names,
#              pattern = "_gluthg38only_bwa_barcode_known_genotypes_WASP_deduped_only_05Jan2023_direct_counting_dedup_only.vcf")

vcf_names
names(vcf_list) <-
  vcf_names


# make this function to get the REF/ALT value of each SNPs per sample and
# assign a uuid (CHR+POS) for each SNP

get_DP_w_uuid <-
  function(x, sample_name) {
    uuid <- str_c(x@fix[, 1],
                  x@fix[, 2],
                  sep = "_")
    REF_name <- str_c(sample_name,
                      "REF_C",
                      sep = "_")
    ALT_name <- str_c(sample_name,
                      "ALT_C",
                      sep = "_")
    DP4_df <-
      extract_info_tidy(x = x,
                        info_fields = "DP4",
                        info_types = c(DP4 = NULL))
    
    DP4_df_split <-
      as.data.frame(str_split(string = DP4_df$DP4,
                              pattern = ",",
                              simplify = T))
    DP4_df_split <-
      as.data.frame(sapply(DP4_df_split,
                           FUN = as.numeric))
    
    DP4_df_return <-
      data.frame(uuid = uuid,
                 REF_C = DP4_df_split$V1 + DP4_df_split$V2, # REF counts
                 ALT_C = DP4_df_split$V3 + DP4_df_split$V4, # ALT counts
                 stringsAsFactors = F)
    colnames(DP4_df_return) <-
      c("uuid", REF_name, ALT_name)
    
    return(DP4_df_return)
  }

# make another list to hold the returning dfs of REF/ALT counts from the function
# a multithread version through mcmapply
# vcf_indiv_count_results <-
#   vector(mode = "list",
#          length = length(vcf_files))
# names(vcf_indiv_count_results) <-
#   vcf_names

# for (i in 1:length(vcf_list)) {
#   print(vcf_names[i])
#   vcf_indiv_count_results[[i]] <-
#     get_DP_w_uuid(x = vcf_list[[i]],
#                   sample_name = vcf_names[i])
# }

# try a multithread version through mcmapply
# # make another list to hold the returning dfs of REF/ALT counts from the function
# vcf_indiv_count_results_para <-
#   vector(mode = "list",
#          length = length(vcf_files))
# names(vcf_indiv_count_results_para) <-
#   vcf_names

# rm(vcf_indiv_count_results_para)
vcf_indiv_count_results <-
  mcmapply(FUN = get_DP_w_uuid, 
           x = vcf_list[1:length(vcf_list)],
           sample_name = vcf_names[1:length(vcf_names)], # if you want to apply vectorised arguments SEQUENTIALLY you have to provide them in-line
                                                         # MoreArgs will apply the WHOLE vector to EVERY x rather than one-by-one
           SIMPLIFY = F,
           USE.NAMES = F,
           mc.preschedule = F,
           mc.cleanup = T,
           mc.cores = 9L)
names(vcf_indiv_count_results) <-
  vcf_names

# merge the 18 individual SNP lists into one main list using the uuid
master_count_df <- data.frame()
for (i in 1:length(vcf_list)) {
  print(vcf_names[i])
  if (i == 1) {
    master_count_df <-
      vcf_indiv_count_results[[1]]
  } else {
    y_to_merge <-
      vcf_indiv_count_results[[i]]
    y_to_merge <-
      y_to_merge[!duplicated(y_to_merge$uuid), ]
    master_count_df <-
      merge(x = master_count_df,
            y = y_to_merge,
            by = "uuid",
            all = T,
            sort = T)
  }
  master_count_df <-
    master_count_df[!duplicated(master_count_df$uuid), ]
}


## load all MGS intervals from the master file
mgs_all_intervals <-
  read.vcfR(file = "~/Data/Databases/GWAS/mgs_new_imputation_guo/output_vcf/Duan_Project_024_samples.vcf",
            limit = 1e+8,
            verbose = T)
mgs_all_intervals <-
  as.data.frame(mgs_all_intervals@fix)
mgs_all_intervals$rsID <-
  str_split(string = mgs_all_intervals$ID,
            pattern = "\\:",
            simplify = T)[, 1]
mgs_all_intervals$uuid <-
  str_c(mgs_all_intervals$CHROM,
        mgs_all_intervals$POS,
        sep = "_")
mgs_all_intervals <-
  mgs_all_intervals[!duplicated(mgs_all_intervals$uuid), ]

mgs_master_table <-
  merge(mgs_all_intervals,
        master_count_df,
        by = "uuid")

mgs_master_table_counts <-
  mgs_master_table[, 11:ncol(mgs_master_table)]
mgs_master_table_meta <-
  mgs_master_table[, 1:10]

mgs_master_table_ref_count <-
  mgs_master_table_counts[, c(2 * (1:18) - 1)]
mgs_master_table_alt_count <-
  mgs_master_table_counts[, c(2 * (1:18))]


mgs_master_table_meta$REF_C <-
  rowSums(mgs_master_table_ref_count, na.rm = T)
mgs_master_table_meta$ALT_C <-
  rowSums(mgs_master_table_alt_count, na.rm = T)
mgs_master_table_meta$DP <-
  mgs_master_table_meta$REF_C +
  mgs_master_table_meta$ALT_C

mgs_master_table_meta <-
  mgs_master_table_meta[mgs_master_table_meta$DP > 19, ]
mgs_master_table_meta <-
  mgs_master_table_meta[(mgs_master_table_meta$REF_C > 1) &
                          (mgs_master_table_meta$ALT_C > 1), ]

mgs_master_table_meta <-
  mgs_master_table_meta[str_detect(string = mgs_master_table_meta$rsID,
                          pattern = "^rs.*"), ]  

clust_df_plot_use <- makeCluster(type = "FORK", 6)
clusterExport(clust_df_plot_use, "mgs_master_table_meta")
mgs_master_table_meta$pVal <- 
  parApply(cl = clust_df_plot_use, 
           X = mgs_master_table_meta, 
           MARGIN = 1,
           FUN = function(x)(binom.test(x = as.numeric(x[11]),
                                        n = as.numeric(x[ncol(mgs_master_table_meta)]),
                                        p = 0.5,
                                        alternative = "t")$p.value))
stopCluster(clust_df_plot_use)
rm(clust_df_plot_use)

mgs_master_table_meta$FDR <- 
  p.adjust(p = mgs_master_table_meta$pVal,
           method = "fdr")
mgs_master_table_meta <-
  mgs_master_table_meta[order(mgs_master_table_meta$pVal), ]

ggplot(mgs_master_table_meta,
       aes(x = (REF_C/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle(paste(#"",
    "bwa_align_count_WASP_calib_dedup_bcftools;\n",
    "mbq=1, mmq=20, -B;\n",
    "NEFM_neg_0hr, DP >= 20, minAllele >= 2, 18 lines;\n",
    # "SNP lists from 23 Sept results, removed all NCALLED >= 14\n",
    nrow(mgs_master_table_meta),
    "SNPs, of",
    sum(mgs_master_table_meta$FDR < 0.05),
    "FDR < 0.05;\n",
    "FDR < 0.05 & REF/DP < 0.5 = ", 
    sum((mgs_master_table_meta$REF_C/mgs_master_table_meta$DP < 0.5) & 
          (mgs_master_table_meta$FDR < 0.05)),
    ";\n FDR < 0.05 & REF/DP > 0.5 = ", "",
    sum((mgs_master_table_meta$REF_C/mgs_master_table_meta$DP > 0.5) & 
          (mgs_master_table_meta$FDR < 0.05)), ";\n",
    "rs2027349: REF_C =",
    mgs_master_table_meta[mgs_master_table_meta$rsID %in% "rs2027349", ]$REF_C,
    ", DP = ",
    mgs_master_table_meta[mgs_master_table_meta$rsID %in% "rs2027349", ]$DP, 
    ",\n Pval =",
    format(signif(mgs_master_table_meta[mgs_master_table_meta$rsID %in% "rs2027349", ]$pVal, digits = 2),
           nsmall = 2), 
    ", FDR =",
    format(signif(mgs_master_table_meta[mgs_master_table_meta$rsID %in% "rs2027349", ]$FDR, digits = 2),
           nsmall = 2))) +
  ylab("-log10P") +
  ylim(0, 75) +
  theme_classic()

## !! NOTE: it appears that the REF/DP balance can be changed by adjusting the
## !! -min-MQ, -q, default 0 and -min-BQ, -Q, default 13 values in bcftools mpileup
## !! (corresponding values of --minimum-mapping-quality, default 20 
## [HCMappingQualityFilter] and
## !! --min-base-quality-score -mbq, default 10 in GATK HaplotypeCaller)

## !! Use -q 29 and -Q 13 resulted in 32260 FDR < 0.05, of REF/DP < 0.5 = 12397
## rs2027349, REF_C = 92, DP = 227 (~ 100 less, too high min-MQ?)
## However gross REF/DP looks OK
# > sum(mgs_master_table_meta$REF_C/mgs_master_table_meta$DP < 0.5)
# [1] 149239
# > sum(mgs_master_table_meta$REF_C/mgs_master_table_meta$DP > 0.5)
# [1] 149605
sum(mgs_master_table_meta$REF_C/mgs_master_table_meta$DP < 0.5)
sum(mgs_master_table_meta$REF_C/mgs_master_table_meta$DP > 0.5)


save(mgs_master_table_meta,
     file = "mgs_master_table_meta_mbq_0_mmq_10.RData")

save(mgs_master_table_meta,
     file = "mgs_master_table_meta_mbq_5_mmq_10.RData")

# save(mgs_master_table_meta,
#      file = "mgs_master_table_meta_mbq_10_mmq_10.RData")

# save(mgs_master_table_meta,
#      file = "mgs_master_table_meta_mbq_10_mmq_20_B.RData")

# save(mgs_master_table_meta,
#      file = "mgs_master_table_meta_mbq_10_mmq_10_B_no_WASP.RData")

# save(mgs_master_table_meta,
#      file = "mgs_master_table_meta_mbq_10_mmq_10_no_B.RData")

# test_df <-
#   extract_info_tidy(x = vcf_list$`0hr_CD_02_NEFM_neg`,
#                     info_fields = "DP4",
#                     info_types = c(DP4 = NULL))
# test_df_split <-
#   as.data.frame(str_split(string = test_df$DP4,
#                           pattern = ",",
#                           simplify = T))
# test_df_split <-
#   as.data.frame(sapply(test_df_split,
#                        FUN = as.numeric))
# 
# vcf_list$`0hr_CD_02_NEFM_neg`@fix[, 1]
