# Siwei 11 May 2023
# Process all Asts (FASTQs trimmed)
# May 2023

# init
library(readr)
library(vcfR)
library(stringr)
library(ggplot2)

library(parallel)

library(MASS)

library(RColorBrewer)
library(grDevices)

# load the table from vcf (the vcf has been prefiltered to include DP >= 20 only)
df_raw <-
  read.vcfR(file = "../Hanwen_ATACseq_test/Ast/WASPed_vcfs/output/Ast_26May2023_trimmed_merged_SNP_VF.vcf",
            limit = 1e9,
            verbose = T)

df_fix <-
  as.data.frame(df_raw@fix)


sum_gt_DP <-
  extract.gt(x = df_raw,
             element = "DP",
             as.numeric = T,
             convertNA = T)
sum_gt_DP[is.na(sum_gt_DP)] <- 0
# rowSums(sum_gt_DP)

sum_gt_AD <-
  extract.gt(x = df_raw,
             element = "AD",
             as.numeric = F,
             convertNA = T)
# sum(!is.na(sum_gt_AD))
# get NCALLED
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

# binom.test(x = df_assembled[[9]][1],
#            n = df_assembled[[11]][1])
#
# binom.test(x = as.numeric(unlist(df_assembled[1, ])[9]),
#            n = as.numeric(unlist(df_assembled[1, ])[11]))
#
# sum(is.na(df_assembled$DP))

clust_df_assembled <- makeCluster(type = "FORK", 10)
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
df_assembled$REF_ratio <-
  df_assembled$REF_N / df_assembled$DP

write.table(df_assembled,
            file = "Astrocytes_18_lines_all_SNPs_30May2023.txt",
            quote = F, sep = "\t",
            row.names = F, col.names = T)


ggplot(df_assembled,
       aes(x = (REF_N/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05,
                                        "FDR < 0.05",
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle(paste("Astrocytes, DP >= 20, minAllele >= 2, 18 lines;\n",
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
  # ylim(0, 150) +
  theme_classic()


df_assembled_FDR_005 <-
  df_assembled[df_assembled$FDR < 0.05, ]

write.table(df_assembled_FDR_005,
            file = "Astrocytes_18_lines_005_30May2023.txt",
            quote = F, sep = "\t",
            row.names = F, col.names = T)

save.image("Astrocytes_18_lines_30May2023.RData")

df_avinput <-
  df_assembled_FDR_005[, c(1,2,2,4,5,3,6:ncol(df_assembled_FDR_005))]
write.table(df_avinput,
            file = "Astrocytes_18_lines_005_30May2023.avinput",
            quote = F, sep = "\t",
            row.names = F, col.names = T)

### test
df_test_assembled <-
  df_assembled[1:50, ]


df_test_assembled$pVal <-
  apply(X = df_test_assembled,
        MARGIN = 1,
        FUN = function(x)(binom.test(x = as.numeric(x[9]),
                                     n = as.numeric(x[11]),
                                     p = 0.5,
                                     alternative = "t")$p.value))
