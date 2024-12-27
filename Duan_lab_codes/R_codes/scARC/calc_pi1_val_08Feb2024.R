# Siwei 08 Feb 2024
# process all 100 (actually ~ 94-97 lines depending on cell type) lines
# Calculate Storey's pi1 value using Yifan's code

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
  
  library(dplyr)
}

# param ####
plan("multisession", workers = 4)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)
set.seed(42)

# load data ####
load("~/NVME/scARC_Duan_018/R_ASoC/Sum_100_lines_21Dec2023.RData")

# Aux functions by Yifan ####
get_pi1 <- 
  function(pval, lambda = 0.5) {
  pi1 <- 
    (length(pval) - sum(pval > lambda) / (1 - lambda)) / length(pval)
  return(pi1)
  }

pi1_A_in_B <- 
  function(A.full, B.full) {
  A.sig <- 
    dplyr::filter(A.full, FDR < 0.05)$ID
  A.sig_in_B <- 
    dplyr::filter(B.full, ID %in% A.sig)
  return(get_pi1(A.sig_in_B$pVal) * nrow(A.sig_in_B) / length(A.sig))
}

## extract the lists w all SNPs ([[1]])
vcf_lst_4_get_pi1 <-
  vector(mode = "list",
         length = length(master_vcf_list))

for (i in 1:length(vcf_lst_4_get_pi1)) {
  print(i)
  vcf_lst_4_get_pi1[[i]] <-
    master_vcf_list[[i]][[1]]
}
names(vcf_lst_4_get_pi1) <-
  names(master_vcf_list)

pairwise_pi1 <-
  data.frame(diag(length(vcf_lst_4_get_pi1)),
             row.names = names(vcf_lst_4_get_pi1))
colnames(pairwise_pi1) <-
  names(vcf_lst_4_get_pi1) # for df, names() == colnames()

## run a nested loop that makes pairwise comparison ####
for (i in 1:(length(vcf_lst_4_get_pi1) - 1)) {
  for (j in (i + 1):length(vcf_lst_4_get_pi1)) {
    print(paste0("A =",
                 names(vcf_lst_4_get_pi1)[i],
                 ", ",
                 "B =",
                 names(vcf_lst_4_get_pi1)[j]))
    pairwise_pi1[i, j] <-
      pi1_A_in_B(A.full = vcf_lst_4_get_pi1[[i]],
                 B.full = vcf_lst_4_get_pi1[[j]])
    pairwise_pi1[j, i] <-
      pi1_A_in_B(A.full = vcf_lst_4_get_pi1[[j]],
                 B.full = vcf_lst_4_get_pi1[[i]])
  }
}


pairwise_pi1.plot <- 
  expand.grid(dimnames(pairwise_pi1))
pairwise_pi1.plot$value <- 
  unlist(unname(pairwise_pi1))
colnames(pairwise_pi1.plot) <- 
  c('Leading','Matched','pi1')

pairwise_pi1.plot$Leading <-
  factor(pairwise_pi1.plot$Leading,
         levels = c("0hr_GABA",
                    "1hr_GABA",
                    "6hr_GABA",
                    "0hr_nmglut",
                    "1hr_nmglut",
                    "6hr_nmglut",
                    "0hr_npglut",
                    "1hr_npglut",
                    "6hr_npglut"))

pairwise_pi1.plot$Matched <-
  factor(pairwise_pi1.plot$Matched,
         levels = c("0hr_GABA",
                    "1hr_GABA",
                    "6hr_GABA",
                    "0hr_nmglut",
                    "1hr_nmglut",
                    "6hr_nmglut",
                    "0hr_npglut",
                    "1hr_npglut",
                    "6hr_npglut"))

ggplot(data = pairwise_pi1.plot, 
       aes(x = Leading, 
           y = Matched,
           fill = pi1)) + 
  geom_tile(aes(fill = pi1)) +
  scale_fill_gradient(high = 'firebrick3', 
                      low = 'white',# low = "deepskyblue4", high = "lightskyblue1",
                      name = expression(pi*'1'),
                      limits = c(0.4, 1),
                      na.value = "white") + 
  theme_minimal() + 
  xlab('Leading Cell Type') + 
  ylab('Matched Cell Type') + 
  theme(axis.title = element_text(size = 12,
                                  face = "bold"), 
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        legend.title = element_text(size = 12,
                                    face = "bold"), 
        legend.text = element_text(size = 12))

sample_order_by_type <-
  c("0hr_GABA",
    "1hr_GABA",
    "6hr_GABA",
    "0hr_nmglut",
    "1hr_nmglut",
    "6hr_nmglut",
    "0hr_npglut",
    "1hr_npglut",
    "6hr_npglut")

df_writeout_pi1 <-
  pairwise_pi1[match(x = sample_order_by_type,
                     table = rownames(pairwise_pi1)),
               match(x = sample_order_by_type,
                     table = colnames(pairwise_pi1))]
write.table(df_writeout_pi1,
            file = "pi1_3x3_matrix.txt",
            quote = F, sep = "\t")

## extract the lists w ASoC SNPs only ([[3]])
vcf_lst_ASoC_only <-
  vector(mode = "list",
         length = length(master_vcf_list))

for (i in 1:length(vcf_lst_ASoC_only)) {
  print(i)
  vcf_lst_ASoC_only[[i]] <-
    master_vcf_list[[i]][[3]]
}
names(vcf_lst_ASoC_only) <-
  names(master_vcf_list)

npglut_0hr_specific_SNPs <-
  vcf_lst_ASoC_only[[3]]
length(vcf_lst_ASoC_only[['1hr_npglut']]$ID %in% npglut_0hr_specific_SNPs$ID)
length(npglut_0hr_specific_SNPs$ID %in% vcf_lst_ASoC_only[['1hr_npglut']]$ID)

npglut_0hr_specific_SNPs <-
  npglut_0hr_specific_SNPs[!(npglut_0hr_specific_SNPs$ID %in% vcf_lst_ASoC_only[['1hr_npglut']]$ID), ]
npglut_0hr_specific_SNPs <-
  npglut_0hr_specific_SNPs[!(npglut_0hr_specific_SNPs$ID %in% vcf_lst_ASoC_only[['6hr_npglut']]$ID), ]

save(npglut_0hr_specific_SNPs,
     file = "npglut_0hr_specific_SNP_list.RData")

for (i in 1:length(peaks_uncombined)) {
  print(paste(i, 
              unique(peaks_uncombined[[i]]$ident)))
}
# unique(peaks_uncombined[[1]]$ident)
npglut_0hr_specific_peaks <-
  peaks_uncombined[[8]]
save(npglut_0hr_specific_peaks,
     file = "npglut_0hr_specific_peaks_list.RData")
