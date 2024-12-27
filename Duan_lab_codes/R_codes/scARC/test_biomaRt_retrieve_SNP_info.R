# Siwei 12 Oct 2023
# Test biomaRt to retrive SNP by rsid 
# (then with the index rsid there will be no need to liftOver between genomes)

# hg19 sometimes returns 2 results for 1 rsid, may not use this way...

# init #####
library(biomaRt)

# define SNP list and the ENSEMBL mart that will be used
SNP_list <- "rs2027349"

listEnsemblArchives()
hg19_ensembl_mart <-
  useEnsembl(biomart = "snp",
             dataset = "hsapiens_snp",
             version = "GRCh37") # use hg19
hg38_ensembl_mart <-
  useEnsembl(biomart = "snp",
             dataset = "hsapiens_snp",
             version = "110") # use hg38

listFilters(mart = hg19_ensembl_mart)
listAttributes(hg38_ensembl_mart)
getBM(attributes = c("allele",
                     "snp",
                     "allele_1",
                     "minor_allele",
                     "mapweight"),
      values = SNP_list,
      filters = "snp_filter",
      mart = hg38_ensembl_mart,
      uniqueRows = T,
      bmHeader = T)

getBM(attributes = c("allele",
                     "snp",
                     "allele_1",
                     "minor_allele",
                     "mapweight"),
      values = SNP_list,
      filters = "snp_filter",
      mart = hg19_ensembl_mart,
      uniqueRows = T,
      bmHeader = T)
