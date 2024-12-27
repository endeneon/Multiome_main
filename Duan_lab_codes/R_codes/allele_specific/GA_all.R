# Siwei 21 Jun 2023
# Process all GAs (FASTQs trimmed (2 batches, 1 re-seqed)
# + previous ones from 2019)
# note the sample names are complex, need processing
# Jun 2023

# init
library(readr)
library(vcfR)
library(stringr)
library(ggplot2)

library(parallel)

library(MASS)

library(RColorBrewer)
library(grDevices)

# load the table from vcf #####
# (the vcf has been prefiltered to include DP >= 20 only)
df_raw <-
  read.vcfR(file = "Aug2023_vcfs/GA_21Aug2023_trimmed_merged_SNP_VF.vcf",
            limit = 2e9,
            verbose = T)

df_fix <-
  as.data.frame(df_raw@fix)

sum_gt_AD <-
  extract.gt(x = df_raw,
             element = "AD",
             as.numeric = F,
             convertNA = T)
# sum(!is.na(sum_gt_AD))

# extract individual data columns #####
## get NCALLED, start from sum_gt_AD #####
ncalled_counts <- sum_gt_AD
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

## get REF #####
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

# get ALT #####
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

## assemble df #####
df_assembled <-
  as.data.frame(cbind(df_fix,
                      data.frame(NCALLED = ncalled_counts),
                      data.frame(REF_N = sum_REF),
                      data.frame(ALT_N = sum_ALT)))

### remove any SNP that does not have rsID #####
df_assembled <-
  df_assembled[str_detect(string = df_assembled$ID,
                          pattern = "^rs.*"),
  ]
### require both REF_N >= 3 and ALT_N >= 3 #####
df_assembled <-
  df_assembled[(df_assembled$REF_N > 2 & df_assembled$ALT_N > 2), ]
df_assembled$DP <-
  (df_assembled$REF_N + df_assembled$ALT_N)
### require DP >= 40 #####
# df_assembled <-
#   df_assembled[df_assembled$DP > 25, ]
df_assembled <-
  df_assembled[!is.na(df_assembled$DP), ]

df_assembled <-
  df_assembled[df_assembled$DP > 39, ]

## calc and check REF ratio #####
df_assembled$REF_ratio <-
  df_assembled$REF_N / df_assembled$DP
sum(df_assembled$REF_ratio > 0.5) / nrow(df_assembled)
sum(df_assembled$REF_ratio < 0.5) / nrow(df_assembled)

# calc binom P value and FDR #####
# save.image("GA_30_lines_21Jun2023.RData")

clust_df_assembled <- makeCluster(type = "FORK", 20)
clusterExport(clust_df_assembled, "df_assembled")
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
# df_assembled$REF_ratio <-
#   df_assembled$REF_N / df_assembled$DP
df_assembled[df_assembled$ID %in% "rs2027349", ]
df_assembled[df_assembled$ID %in% "rs1532278", ]
# write output #####

write.table(df_assembled,
            file = "GA_30_lines_all_SNPs_25Sept2023.txt",
            quote = F, sep = "\t",
            row.names = F, col.names = T)


# df_assembled_plot <-
#   df_assembled[sample(x = nrow(df_assembled),
#                       size = 10000), ]

df_assembled_FDR_005 <-
  df_assembled[df_assembled$FDR < 0.05, ]

write.table(df_assembled_FDR_005,
            file = "GA_30_lines_005_22Aug2023.txt",
            quote = F, sep = "\t",
            row.names = F, col.names = T)

df_avinput <-
  df_assembled_FDR_005[, c(1,2,2,4,5,3,6:ncol(df_assembled_FDR_005))]
write.table(df_avinput,
            file = "GA_30_lines_005_22Aug2023.avinput",
            quote = F, sep = "\t",
            row.names = F, col.names = T)

df_assembled <-
  readRDS("GA_30_lines_22Aug2023.RDs")

ggplot(df_assembled,
       aes(x = (REF_N/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05,
                                        "FDR < 0.05",
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle(paste("GABA, DP >= 40, minAllele >= 2, 30 lines;\n",
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
# make violin plot (note: 10000 subsetted) #####

saveRDS(object = df_assembled,
        file = "GA_30_lines_22Aug2023.RDs")
