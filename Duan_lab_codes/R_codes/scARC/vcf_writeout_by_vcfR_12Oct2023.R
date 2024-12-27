# Siwei 12 Oct 2023
# test writeout vcf using vcfR package

# init #####
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

master_vcf_list <-
  readRDS(file = "Duan_Project_018_024_025_029_ASoC_74_lines_11Oct2023.RDs")

## test if ggvenn can return intersected lists ######
## can, but not very efficient
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

list_to_return <-
  ggvenn(data = list_venn_2_plot,
       fill_color = brewer.pal(n = 3, 
                               name = "Dark2"),
       stroke_size = 0.5,
       set_name_size = 6,
       text_size = 4, 
       show_elements = T) +
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

rm(list_to_return)
rm(list_venn_2_plot)

## test the "traditional" method #####
# for each cell type, 8 lists will be required:
# 3 lists of each time point;
# 3 lists unique to each time point;
# 1 list that is the union of at least 1 time point;
# 1 list that is the shared SNPS of all 3 time points

### define function #####
calc_lst_index <-
  function(i) {
    lst_index <-
      vector(mode = "list",
             length = 8L)
    names(lst_index) <-
      c("0hr",
        "1hr",
        "6hr",
        "unique_0hr",
        "unique_1hr",
        "unique_6hr",
        "union_atleast_1_time",
        "shared_3_times")
    
    # i <- 1
    lst_index[[1]] <-
      master_vcf_list[[i]][[1]][, 1:8]
    lst_index[[1]]$QUAL <- 100
    lst_index[[1]]$INFO <- "AF=0"
    lst_index[[1]]$INFO[lst_index[[1]]$ID %in% 
                          master_vcf_list[[i]][[3]]$ID] <-
      "AF=1"
    # sum(lst_index[[1]]$ID %in%
    #       master_vcf_list[[i]][[3]]$ID)
    
    lst_index[[2]] <-
      master_vcf_list[[i + 3]][[1]][, 1:8]
    lst_index[[2]]$QUAL <- 100
    lst_index[[2]]$INFO <- "AF=0"
    lst_index[[2]]$INFO[lst_index[[2]]$ID %in% 
                          master_vcf_list[[i + 3]][[3]]$ID] <-
      "AF=1"
    
    lst_index[[3]] <-
      master_vcf_list[[i + 6]][[1]][, 1:8]
    lst_index[[3]]$QUAL <- 100
    lst_index[[3]]$INFO <- "AF=0"
    lst_index[[3]]$INFO[lst_index[[3]]$ID %in% 
                          master_vcf_list[[i + 6]][[3]]$ID] <-
      "AF=1"
    
    
    lst_index[[4]] <-
      master_vcf_list[[i]][[1]][, 1:8]
    lst_index[[4]]$QUAL <- 100
    lst_index[[4]]$INFO <- "AF=0"
    lst_index[[4]]$INFO[lst_index[[4]]$ID %in%
                          master_vcf_list[[i]][[3]]$ID 
                        [!(master_vcf_list[[i]][[3]]$ID %in% 
                             unique(c(master_vcf_list[[i + 3]][[3]]$ID,
                                      master_vcf_list[[i + 6]][[3]]$ID)))]] <-
      "AF=1"
    sum(lst_index[[4]]$INFO == "AF=1")
    
    lst_index[[5]] <-
      master_vcf_list[[i + 3]][[1]][, 1:8]
    lst_index[[5]]$QUAL <- 100
    lst_index[[5]]$INFO <- "AF=0"
    lst_index[[5]]$INFO[lst_index[[5]]$ID %in%
                          master_vcf_list[[i + 3]][[3]]$ID 
                        [!(master_vcf_list[[i + 3]][[3]]$ID %in% 
                             unique(c(master_vcf_list[[i]][[3]]$ID,
                                      master_vcf_list[[i + 6]][[3]]$ID)))]] <-
      "AF=1"
    sum(lst_index[[5]]$INFO == "AF=1")
    
    lst_index[[6]] <-
      master_vcf_list[[i + 6]][[1]][, 1:8]
    lst_index[[6]]$QUAL <- 100
    lst_index[[6]]$INFO <- "AF=0"
    lst_index[[6]]$INFO[lst_index[[6]]$ID %in%
                          master_vcf_list[[i + 6]][[3]]$ID 
                        [!(master_vcf_list[[i + 6]][[3]]$ID %in% 
                             unique(c(master_vcf_list[[i]][[3]]$ID,
                                      master_vcf_list[[i + 3]][[3]]$ID)))]] <-
      "AF=1"
    sum(lst_index[[6]]$INFO == "AF=1")
    
    
    lst_index[[7]] <-
      as.data.frame(rbind(lst_index[[1]],
                          lst_index[[2]],
                          lst_index[[3]]))
    lst_index[[7]]$QUAL <- 100
    lst_index[[7]]$INFO <- "AF=0"
    lst_index[[7]] <-
      lst_index[[7]][!duplicated(lst_index[[7]]$ID), ]
    lst_index[[7]]$INFO[lst_index[[7]]$ID %in%
                          unique(c(master_vcf_list[[i]][[3]]$ID,
                                   master_vcf_list[[i + 3]][[3]]$ID,
                                   master_vcf_list[[i + 6]][[3]]$ID))] <-
      "AF=1"
    sum(lst_index[[7]]$INFO == "AF=1")
    
    
    lst_index[[8]] <-
      as.data.frame(rbind(lst_index[[1]],
                          lst_index[[2]],
                          lst_index[[3]]))
    lst_index[[8]]$QUAL <- 100
    lst_index[[8]]$INFO <- "AF=0"
    lst_index[[8]] <-
      lst_index[[8]][!duplicated(lst_index[[8]]$ID), ]
    lst_index[[8]]$INFO[lst_index[[8]]$ID %in%
                          master_vcf_list[[i]][[3]]$ID 
                        [master_vcf_list[[i]][[3]]$ID %in%
                            master_vcf_list[[i + 3]][[3]]$ID
                          [master_vcf_list[[i + 3]][[3]]$ID %in%
                              master_vcf_list[[i + 6]][[3]]$ID]]] <-
      "AF=1"
    sum(lst_index[[8]]$INFO == "AF=1")
    
    return(lst_index)
  }

write_vcf_data <-
  function(output_dir = "",
           vcf_prefix = "cell_type",
           write_list = write_list) {
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    
    for (k in 1:length(write_list)) {
      df_writeout <-
        write_list[[k]]
      colnames(df_writeout)[1] <- '#CHROM'
      write_path <-
        paste0(output_dir,
               '/',
               vcf_prefix,
               '_',
               names(write_list)[k],
               '.txt')
      print(write_path)
      write.table(df_writeout,
                  file = write_path,
                  quote = F, sep = "\t",
                  row.names = F, col.names = T)
    }
    return(0)
  }

# if (!dir.exists("VCF_data_part")) {
#   dir.create("VCF_data_part")
# }

names(master_vcf_list)
# > names(master_vcf_list)
# [1] "0hr_GABA"   "0hr_nmglut" "0hr_npglut" "1hr_GABA"   "1hr_nmglut" "1hr_npglut" "6hr_GABA"  
# [8] "6hr_nmglut" "6hr_npglut"
### GABA #####
list_to_write <-
  calc_lst_index(i = 1)
names(list_to_write)

write_vcf_data(output_dir = "VCF_data_part",
               vcf_prefix = "GABA",
               write_list = list_to_write)


### nmglut #####
list_to_write <-
  calc_lst_index(i = 2)
names(list_to_write)

write_vcf_data(output_dir = "VCF_data_part",
               vcf_prefix = "nmglut",
               write_list = list_to_write)

### npglut #####
list_to_write <-
  calc_lst_index(i = 3)
names(list_to_write)

write_vcf_data(output_dir = "VCF_data_part",
               vcf_prefix = "npglut",
               write_list = list_to_write)
