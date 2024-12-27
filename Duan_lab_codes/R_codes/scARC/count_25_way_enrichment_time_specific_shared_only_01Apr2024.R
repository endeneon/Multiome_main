# Siwei 27 Mar 2023
# count ASoC enrichment in 25-way regions
# generate necessary bed6 files for bash counting

# init
{
  library(readr)
  # library(vcfR)
  library(stringr)
  library(ggplot2)
  
  library(scales)
  
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

load("~/NVME/scARC_Duan_018/R_ASoC/Sum_100_lines_21Dec2023.RData")

ASoC_intersect_sum_dir <-
  "ASoC_25_ways_sum_dir"
# dir.create(path = ASoC_intersect_sum_dir,
#            recursive = T)

## make a list that holds all ASoC SNPs
names(master_vcf_list)
cell_type_list <-
  c("GABA", "nmglut", "npglut")
cell_time_list <-
  c("hr_0", "hr_1", "hr_6")

for (i in 1:length(cell_type_list)) {
  all_ASoC_SNPs <-
    vector(mode = "list",
           length = length(cell_time_list))
  names(all_ASoC_SNPs) <-
    cell_time_list
  for (j in 1:length(cell_time_list)) {
    # print(paste(i, j))
    all_ASoC_SNPs[[j]] <-
      master_vcf_list[[3 * (j - 1) + i]][[3]]
    print(3 * (j - 1) + i)
  }
  
  ### 0 hr output #####
  {
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[['hr_0']][!(all_ASoC_SNPs[['hr_0']]$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[!(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
    print(paste("length of ASoC list is",
                nrow(ASoC_stage_specific_output)))
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- '.'
    
    write.table(table_writeout,
                file = paste0(ASoC_intersect_sum_dir,
                              "/",
                              cell_type_list[i],
                              "_",
                              "0hr_specific.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }

  ### 1 hr output #####
  {
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[['hr_1']][!(all_ASoC_SNPs[['hr_1']]$ID %in% all_ASoC_SNPs[['hr_0']]$ID), ]
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[!(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
    print(paste("length of ASoC list is",
                nrow(ASoC_stage_specific_output)))
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- '.'
    
    write.table(table_writeout,
                file = paste0(ASoC_intersect_sum_dir,
                              "/",
                              cell_type_list[i],
                              "_",
                              "1hr_specific.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
  
  ### 6 hr output #####
  {
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[['hr_6']][!(all_ASoC_SNPs[['hr_6']]$ID %in% all_ASoC_SNPs[['hr_0']]$ID), ]
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[!(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
    print(paste("length of ASoC list is",
                nrow(ASoC_stage_specific_output)))
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- '.'
    
    write.table(table_writeout,
                file = paste0(ASoC_intersect_sum_dir,
                              "/",
                              cell_type_list[i],
                              "_",
                              "6hr_specific.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
  
  ### Common all 3 output
  {
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[['hr_0']][(all_ASoC_SNPs[['hr_0']]$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
    print(paste("length of ASoC list is",
                nrow(ASoC_stage_specific_output)))
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- '.'
    
    write.table(table_writeout,
                file = paste0(ASoC_intersect_sum_dir,
                              "/",
                              cell_type_list[i],
                              "_",
                              "shared_in_all.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
  
  ### At least one output
  {
    ASoC_stage_specific_output <-
      do.call(what = rbind,
              args = all_ASoC_SNPs)
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[!duplicated(ASoC_stage_specific_output$ID), ]
    #   all_ASoC_SNPs[['hr_0']][(all_ASoC_SNPs[['hr_0']]$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
    # ASoC_stage_specific_output <-
    #   ASoC_stage_specific_output[(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
    print(paste("length of ASoC list is",
                nrow(ASoC_stage_specific_output)))
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- '.'
    
    write.table(table_writeout,
                file = paste0(ASoC_intersect_sum_dir,
                              "/",
                              cell_type_list[i],
                              "_",
                              "at_least_one.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
}
  
# make separate, dedicated lists for SNPs that a common/at least one #####
# of each time
## shared, common in all 3 cell types of the same stage #####
cell_type_list <-
  c("sharedByAllCellTypes")
cell_time_list <-
  c("hr_0", "hr_1", "hr_6")

for (i in 1:length(cell_type_list)) {
  all_ASoC_SNPs <-
    vector(mode = "list",
           length = length(cell_time_list))
  names(all_ASoC_SNPs) <-
    cell_time_list
  for (j in 1:length(cell_time_list)) {
    # print(paste(i, j))
    df_GABA <-
      master_vcf_list[[3 * (j - 1) + 1]][[3]]
    df_nmglut <-
      master_vcf_list[[3 * (j - 1) + 2]][[3]]
    df_npglut <-
      master_vcf_list[[3 * (j - 1) + 3]][[3]]
    
    common_rsID <-
      df_GABA$ID[(df_GABA$ID %in% df_nmglut$ID)]
    common_rsID <-
      common_rsID[common_rsID %in% df_npglut$ID]
    
    df_all <-
      df_GABA[df_GABA$ID %in% common_rsID, ]
    
    all_ASoC_SNPs[[j]] <-
      df_all
    print(3 * (j - 1 + i))
  }
  
  ### 0 hr output #####
  {
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[['hr_0']][!(all_ASoC_SNPs[['hr_0']]$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[!(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
    print(paste("length of ASoC list is",
                nrow(ASoC_stage_specific_output)))
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- '.'
    
    write.table(table_writeout,
                file = paste0(ASoC_intersect_sum_dir,
                              "/",
                              cell_type_list[i],
                              "_",
                              "0hr_specific.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
  
  ### 1 hr output #####
  {
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[['hr_1']][!(all_ASoC_SNPs[['hr_1']]$ID %in% all_ASoC_SNPs[['hr_0']]$ID), ]
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[!(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
    print(paste("length of ASoC list is",
                nrow(ASoC_stage_specific_output)))
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- '.'
    
    write.table(table_writeout,
                file = paste0(ASoC_intersect_sum_dir,
                              "/",
                              cell_type_list[i],
                              "_",
                              "1hr_specific.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
  
  ### 6 hr output #####
  {
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[['hr_6']][!(all_ASoC_SNPs[['hr_6']]$ID %in% all_ASoC_SNPs[['hr_0']]$ID), ]
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[!(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
    print(paste("length of ASoC list is",
                nrow(ASoC_stage_specific_output)))
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- '.'
    
    write.table(table_writeout,
                file = paste0(ASoC_intersect_sum_dir,
                              "/",
                              cell_type_list[i],
                              "_",
                              "6hr_specific.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
  
  ### Common all 3 output
  {
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[['hr_0']][(all_ASoC_SNPs[['hr_0']]$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
    print(paste("length of ASoC list is",
                nrow(ASoC_stage_specific_output)))
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- '.'
    
    write.table(table_writeout,
                file = paste0(ASoC_intersect_sum_dir,
                              "/",
                              cell_type_list[i],
                              "_",
                              "shared_in_all.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
  
  ### At least one output
  {
    ASoC_stage_specific_output <-
      do.call(what = rbind,
              args = all_ASoC_SNPs)
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[!duplicated(ASoC_stage_specific_output$ID), ]
    #   all_ASoC_SNPs[['hr_0']][(all_ASoC_SNPs[['hr_0']]$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
    # ASoC_stage_specific_output <-
    #   ASoC_stage_specific_output[(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
    print(paste("length of ASoC list is",
                nrow(ASoC_stage_specific_output)))
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- '.'
    
    write.table(table_writeout,
                file = paste0(ASoC_intersect_sum_dir,
                              "/",
                              cell_type_list[i],
                              "_",
                              "at_least_one.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
}

## SNP in at least one of the 3 cell types #####
cell_type_list <-
  c("atLeastInOneCellType")
cell_time_list <-
  c("hr_0", "hr_1", "hr_6")

for (i in 1:length(cell_type_list)) {
  all_ASoC_SNPs <-
    vector(mode = "list",
           length = length(cell_time_list))
  names(all_ASoC_SNPs) <-
    cell_time_list
  for (j in 1:length(cell_time_list)) {
    # print(paste(i, j))
    df_GABA <-
      master_vcf_list[[3 * (j - 1) + 1]][[3]]
    df_nmglut <-
      master_vcf_list[[3 * (j - 1) + 2]][[3]]
    df_npglut <-
      master_vcf_list[[3 * (j - 1) + 3]][[3]]
    
    df_all <-
      as.data.frame(rbind(df_GABA,
                          df_nmglut,
                          df_npglut))
    df_all <-
      df_all[!duplicated(df_all$ID), ]
    
    all_ASoC_SNPs[[j]] <-
      df_all
    print(3 * (j - 1 + i))
  }
  
  ### 0 hr output #####
  {
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[['hr_0']][!(all_ASoC_SNPs[['hr_0']]$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[!(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
    print(paste("length of ASoC list is",
                nrow(ASoC_stage_specific_output)))
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- '.'
    
    write.table(table_writeout,
                file = paste0(ASoC_intersect_sum_dir,
                              "/",
                              cell_type_list[i],
                              "_",
                              "0hr_specific.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
  
  ### 1 hr output #####
  {
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[['hr_1']][!(all_ASoC_SNPs[['hr_1']]$ID %in% all_ASoC_SNPs[['hr_0']]$ID), ]
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[!(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
    print(paste("length of ASoC list is",
                nrow(ASoC_stage_specific_output)))
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- '.'
    
    write.table(table_writeout,
                file = paste0(ASoC_intersect_sum_dir,
                              "/",
                              cell_type_list[i],
                              "_",
                              "1hr_specific.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
  
  ### 6 hr output #####
  {
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[['hr_6']][!(all_ASoC_SNPs[['hr_6']]$ID %in% all_ASoC_SNPs[['hr_0']]$ID), ]
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[!(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
    print(paste("length of ASoC list is",
                nrow(ASoC_stage_specific_output)))
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- '.'
    
    write.table(table_writeout,
                file = paste0(ASoC_intersect_sum_dir,
                              "/",
                              cell_type_list[i],
                              "_",
                              "6hr_specific.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
  
  ### Common all 3 output
  {
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[['hr_0']][(all_ASoC_SNPs[['hr_0']]$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
    print(paste("length of ASoC list is",
                nrow(ASoC_stage_specific_output)))
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- '.'
    
    write.table(table_writeout,
                file = paste0(ASoC_intersect_sum_dir,
                              "/",
                              cell_type_list[i],
                              "_",
                              "shared_in_all.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
  
  ### At least one output
  {
    ASoC_stage_specific_output <-
      do.call(what = rbind,
              args = all_ASoC_SNPs)
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[!duplicated(ASoC_stage_specific_output$ID), ]
    #   all_ASoC_SNPs[['hr_0']][(all_ASoC_SNPs[['hr_0']]$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
    # ASoC_stage_specific_output <-
    #   ASoC_stage_specific_output[(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
    print(paste("length of ASoC list is",
                nrow(ASoC_stage_specific_output)))
    table_writeout <-
      ASoC_stage_specific_output[, c(1, 2, 2, 6)]
    table_writeout[[2]] <-
      as.numeric(table_writeout[[2]]) - 1
    table_writeout$SCORE <- 100
    table_writeout$STRAND <- '.'
    
    write.table(table_writeout,
                file = paste0(ASoC_intersect_sum_dir,
                              "/",
                              cell_type_list[i],
                              "_",
                              "at_least_one.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
}
#####


# length("a" %in% c("a", "b","c"))
# length(c("a", "b","c") %in% "a")

hg38_genome_size <- 3234830000

df_raw_data <-
  read.table(file = "ASoC_25_ways_sum_dir/intersect_table.txt",
             header = F, sep = "\t",
             quote = "")
colnames(df_raw_data) <-
  c("Sample_SNP", "SNP_coverage", "Interval_region", "Interval_coverage", "Overlap_coverage")
# df_raw_data <-
#   as.data.frame(rbind(df_raw_data,
#                       c("Z_Enh_all",
#                         )))

df_2_plot <-
  df_raw_data
df_2_plot$Enrichment <-
  (df_2_plot$Overlap_coverage / hg38_genome_size) / 
  ((df_2_plot$Interval_coverage / hg38_genome_size) * (df_2_plot$SNP_coverage / hg38_genome_size))
df_2_plot$log10_enrichment <-
  log10(df_2_plot$Enrichment)
df_2_plot$log10_enrichment[is.infinite(df_2_plot$log10_enrichment)] <- 0

df_2_plot$Cell_type <-
  str_split(df_2_plot$Sample_SNP,
            pattern = "_",
            n = 2,
            simplify = T)[ , 1]
df_2_plot$Cell_time <-
  str_split(df_2_plot$Sample_SNP,
            pattern = "_",
            n = 2,
            simplify = T)[ , 2]

df_2_plot$Region <-
  str_split(df_2_plot$Interval_region,
            pattern = "_",
            n = 2,
            simplify = T)[, 2]
df_2_plot$Region <-
  str_replace_all(df_2_plot$Region,
             pattern = "^12F",
             replacement = "Me_F")
df_2_plot$Region <-
  str_replace_all(df_2_plot$Region,
              pattern = "^12O",
              replacement = "Me_O")

df_2_plot$Cell_type <-
  factor(df_2_plot$Cell_type,
         levels = c("GABA",
                    "nmglut",
                    "npglut",
                    "atLeastInOneCellType",
                    "sharedByAllCellTypes"))

# exclude atLeastInOneCellType
df_2_plot <-
  df_2_plot[!(df_2_plot$Cell_type %in% "atLeastInOneCellType"), ]
df_2_plot <-
  df_2_plot[!(df_2_plot$Cell_time %in% "at_least_one"), ]

ggplot(df_2_plot,
       aes(x = Region,
           y = log10_enrichment)) +
  geom_bar(aes(fill = factor(ifelse(test = as.numeric(factor(df_2_plot$Region)) %% 2,
                                    yes = 0,
                                    no = 1)),
               colour = factor(as.numeric(factor(df_2_plot$Cell_type)))),
           stat = "identity",
           linewidth = 0.2) +
  geom_text(aes(label = df_2_plot$Overlap_coverage),
            angle = 270,
            inherit.aes = T,
            vjust = 0.5,
            hjust = 0,
            colour = "black") +
  scale_fill_manual(values = c("pink", "lightblue"),
                    # values = brewer.pal(n = 3,
                    #                     name = "Dark2")[c(1, 2)],
                    guide = "none") +
  scale_colour_manual(values = brewer.pal(n = 5,
                                          name = "Dark2"),
                      guide = "none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 270,
                                   vjust = 0.5,
                                   hjust = 0,
                                   colour = "black")) +
  facet_wrap(Cell_type ~ Cell_time)

df_2_plot$GREAT_p <-
  unlist(apply(X = as.matrix(df_2_plot), 
               MARGIN = 1, 
               FUN = function(x) {
                 0 - log10(binom.test(x = as.numeric(x[5]),
                                      n = as.numeric(x[2]),
                                      p = as.numeric(x[4]) / hg38_genome_size,
                                      alternative = "g")$p.value)
               }))

  # unlist(apply(X = as.matrix(df_2_plot), 
  #              MARGIN = 1, 
  #              FUN = function(x) {
  #                return(x[2])
  #              }))
df_2_plot$GREAT_p[is.infinite(df_2_plot$GREAT_p)] <- 500
  
ggplot(df_2_plot,
       aes(x = Region,
           y = GREAT_p)) +
  geom_bar(aes(fill = factor(ifelse(test = as.numeric(factor(df_2_plot$Region)) %% 2,
                                    yes = 0,
                                    no = 1)),
               colour = factor(as.numeric(factor(df_2_plot$Cell_type)))),
           stat = "identity",
           linewidth = 0.2) +
  geom_text(aes(label = df_2_plot$Overlap_coverage),
            angle = 270,
            inherit.aes = T,
            vjust = 0.5,
            hjust = 0,
            colour = "black") +
  scale_fill_manual(values = c("pink", "lightblue"),
                    # values = brewer.pal(n = 3,
                    #                     name = "Dark2")[c(1, 2)],
                    guide = "none") +
  scale_colour_manual(values = brewer.pal(n = 5,
                                          name = "Dark2"),
                      guide = "none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 270,
                                   vjust = 0.5,
                                   hjust = 0,
                                   colour = "black")) +
  # facet_grid(Cell_type ~ Cell_time)
  facet_wrap(~ Cell_time)



#### test #####
common_rsID <-
  df_GABA$ID[(df_GABA$ID %in% df_nmglut$ID)]
common_rsID <-
  common_rsID[common_rsID %in% df_npglut$ID]

## plot Enh vs Prom in cell-type-specific (3 timepoints) /shared_by_all #####
## make non-SNP lists of atLeastInOncCellType
## shared, common in all 3 cell types of the same stage #####
df_raw_data <-
  read.table(file = "ASoC_25_ways_sum_dir/intersect_table_EnhvsProm.tsv",
             header = F, sep = "\t",
             quote = "")
colnames(df_raw_data) <-
  c("Sample_SNP", "SNP_coverage", "Interval_region", "Interval_coverage", "Overlap_coverage")

df_2_plot <- df_raw_data
df_2_plot$Percentage <-
  df_2_plot$Overlap_coverage / df_2_plot$SNP_coverage

df_2_plot$Cell_type <-
  str_split(df_2_plot$Sample_SNP,
            pattern = "_",
            n = 2,
            simplify = T)[ , 1]
df_2_plot$Cell_time <-
  str_split(df_2_plot$Sample_SNP,
            pattern = "_",
            n = 2,
            simplify = T)[ , 2]

df_2_plot$Region <-
  str_split(df_2_plot$Interval_region,
            pattern = "_",
            n = 2,
            simplify = T)[, 2]
df_2_plot$Region <-
  str_replace_all(df_2_plot$Region,
                  pattern = "^12F",
                  replacement = "Me_F")
df_2_plot$Region <-
  str_replace_all(df_2_plot$Region,
                  pattern = "^12O",
                  replacement = "Me_O")

df_2_plot$Region[df_2_plot$Region %in% "zEnh_all"] <- "Enhancers"
df_2_plot$Region[df_2_plot$Region %in% "zProm"] <- "Promoters"


unique(df_2_plot$Cell_type)
df_2_plot <-
  df_2_plot[!(df_2_plot$Cell_type == "atLeastInOneCellType"), ]
df_2_plot <-
  df_2_plot[!(df_2_plot$Cell_type == "sharedByAllCellTypes"), ]
# df_2_plot <-
#   df_2_plot[(df_2_plot$Cell_type == "npglut"), ]
df_2_plot <-
  df_2_plot[!(df_2_plot$Cell_time == "at_least_one"), ]

# df_2_plot

ggplot(df_2_plot,
       aes(x = Cell_time,
           y = Percentage)) +
  geom_bar(aes(fill = Region,
               colour = Cell_type),
           stat = "identity",
           position = "dodge",
           # position = position_dodge(width = 1),
           # width = 1,
           linewidth = 0) +
  geom_text(aes(label = paste(df_2_plot$Overlap_coverage,
                              df_2_plot$SNP_coverage,
                              sep = ' / ')),
            angle = 270,
            inherit.aes = T,
            vjust = 0.5,
            hjust = 0.5,
            position = position_dodge2(width = 0.9),
            colour = "black") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Set1"))  +
  scale_colour_manual(values = brewer.pal(n = 5,
                                          name = "Set2"),
                      guide = "none") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 0.8),
                     labels = scales::percent) +
  labs(x = "") +
  theme_classic() +
  theme(axis.text.x.bottom = element_text(angle = 315,
                                   vjust = 0.5,
                                   hjust = 0,
                                   colour = "black")) +
  facet_grid(. ~ Cell_type)

# df_2_plot <-
#   df_2_plot[(df_2_plot$Cell_type == "atLeastInOneCellType"), ]
# df_2_plot <-
#   df_2_plot[(df_2_plot$Cell_type == "npglut"), ]
df_2_plot <-
  df_2_plot[!(df_2_plot$Cell_time == "at_least_one"), ]
df_2_plot <-
  df_2_plot[(df_2_plot$Cell_type == "atLeastInOneCellType"), ]
# df_2_plot

ggplot(df_2_plot,
       aes(x = Cell_time,
           y = Percentage)) +
  geom_bar(aes(fill = Region,
               colour = Cell_type),
           stat = "identity",
           position = "dodge",
           # position = position_dodge(width = 1),
           # width = 1,
           linewidth = 0) +
  geom_text(aes(label = paste(df_2_plot$Overlap_coverage,
                              df_2_plot$SNP_coverage,
                              sep = ' / ')),
            angle = 270,
            inherit.aes = T,
            vjust = 0.5,
            hjust = 0.75,
            position = position_dodge2(width = 0.9),
            colour = "black") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Set1"))  +
  scale_colour_manual(values = brewer.pal(n = 5,
                                          name = "Set2"),
                      guide = "none") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 0.8),
                     labels = scales::percent) +
  labs(x = "") +
  theme_classic() +
  theme(axis.text.x.bottom = element_text(angle = 315,
                                          vjust = 0.5,
                                          hjust = 0,
                                          colour = "black")) +
  facet_grid(. ~ Cell_type)
