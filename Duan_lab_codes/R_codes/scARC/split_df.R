# Siwei 08 Apr 2024
# Split excel cells and convert to row-wise

library(readxl)
library(stringr)

df_raw <-
  read_excel("split_rs.xlsx",
             sheet = 1)

list_split <-
  mapply(str_split, 
         df_raw$elements,
         MoreArgs = list(pattern = ' '),
         SIMPLIFY = T)

final_df <-
  lapply(1:length(list_split), function(x) {
    return(unlist(str_split(list_split[[x]],
                            pattern = ' ',
                            simplify = T)))
  })


load("~/NVME/scARC_Duan_018/R_ASoC/Sum_100_lines_21Dec2023.RData")

ASoC_intersect_sum_dir <-
  "ASoC_25_ways_sum_dir"
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
  
  # ### 0 hr output #####
  # {
  #   ASoC_stage_specific_output <-
  #     all_ASoC_SNPs[['hr_0']][!(all_ASoC_SNPs[['hr_0']]$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
  #   ASoC_stage_specific_output <-
  #     ASoC_stage_specific_output[!(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
  #   print(paste("length of ASoC list is",
  #               nrow(ASoC_stage_specific_output)))
  #   table_writeout <-
  #     ASoC_stage_specific_output[, c(1, 2, 2, 6)]
  #   table_writeout[[2]] <-
  #     as.numeric(table_writeout[[2]]) - 1
  #   table_writeout$SCORE <- 100
  #   table_writeout$STRAND <- '.'
  #   
  #   write.table(table_writeout,
  #               file = paste0(ASoC_intersect_sum_dir,
  #                             "/",
  #                             cell_type_list[i],
  #                             "_",
  #                             "0hr_specific.txt"),
  #               quote = F, sep = "\t",
  #               row.names = F, col.names = F)
  # }
  # 
  # ### 1 hr output #####
  # {
  #   ASoC_stage_specific_output <-
  #     all_ASoC_SNPs[['hr_1']][!(all_ASoC_SNPs[['hr_1']]$ID %in% all_ASoC_SNPs[['hr_0']]$ID), ]
  #   ASoC_stage_specific_output <-
  #     ASoC_stage_specific_output[!(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
  #   print(paste("length of ASoC list is",
  #               nrow(ASoC_stage_specific_output)))
  #   table_writeout <-
  #     ASoC_stage_specific_output[, c(1, 2, 2, 6)]
  #   table_writeout[[2]] <-
  #     as.numeric(table_writeout[[2]]) - 1
  #   table_writeout$SCORE <- 100
  #   table_writeout$STRAND <- '.'
  #   
  #   write.table(table_writeout,
  #               file = paste0(ASoC_intersect_sum_dir,
  #                             "/",
  #                             cell_type_list[i],
  #                             "_",
  #                             "1hr_specific.txt"),
  #               quote = F, sep = "\t",
  #               row.names = F, col.names = F)
  # }
  # 
  # ### 6 hr output #####
  # {
  #   ASoC_stage_specific_output <-
  #     all_ASoC_SNPs[['hr_6']][!(all_ASoC_SNPs[['hr_6']]$ID %in% all_ASoC_SNPs[['hr_0']]$ID), ]
  #   ASoC_stage_specific_output <-
  #     ASoC_stage_specific_output[!(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
  #   print(paste("length of ASoC list is",
  #               nrow(ASoC_stage_specific_output)))
  #   table_writeout <-
  #     ASoC_stage_specific_output[, c(1, 2, 2, 6)]
  #   table_writeout[[2]] <-
  #     as.numeric(table_writeout[[2]]) - 1
  #   table_writeout$SCORE <- 100
  #   table_writeout$STRAND <- '.'
  #   
  #   write.table(table_writeout,
  #               file = paste0(ASoC_intersect_sum_dir,
  #                             "/",
  #                             cell_type_list[i],
  #                             "_",
  #                             "6hr_specific.txt"),
  #               quote = F, sep = "\t",
  #               row.names = F, col.names = F)
  # }
  # 
  # ### Common all 3 output
  # {
  #   ASoC_stage_specific_output <-
  #     all_ASoC_SNPs[['hr_0']][(all_ASoC_SNPs[['hr_0']]$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
  #   ASoC_stage_specific_output <-
  #     ASoC_stage_specific_output[(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
  #   print(paste("length of ASoC list is",
  #               nrow(ASoC_stage_specific_output)))
  #   table_writeout <-
  #     ASoC_stage_specific_output[, c(1, 2, 2, 6)]
  #   table_writeout[[2]] <-
  #     as.numeric(table_writeout[[2]]) - 1
  #   table_writeout$SCORE <- 100
  #   table_writeout$STRAND <- '.'
  #   
  #   write.table(table_writeout,
  #               file = paste0(ASoC_intersect_sum_dir,
  #                             "/",
  #                             cell_type_list[i],
  #                             "_",
  #                             "shared_in_all.txt"),
  #               quote = F, sep = "\t",
  #               row.names = F, col.names = F)
  # }
  
  ### 0 & 6 hr only
  {
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[['hr_0']][(all_ASoC_SNPs[['hr_0']]$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
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
                              "0_and_6_hr_exclusive.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
  
  ### 0 & 1 hr only
  {
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[['hr_0']][(all_ASoC_SNPs[['hr_0']]$ID %in% all_ASoC_SNPs[['hr_1']]$ID), ]
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
                              "0_and_1_hr_exclusive.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
  ### 1 & 6 hr only
  {
    ASoC_stage_specific_output <-
      all_ASoC_SNPs[['hr_1']][(all_ASoC_SNPs[['hr_1']]$ID %in% all_ASoC_SNPs[['hr_6']]$ID), ]
    ASoC_stage_specific_output <-
      ASoC_stage_specific_output[!(ASoC_stage_specific_output$ID %in% all_ASoC_SNPs[['hr_0']]$ID), ]
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
                              "1_and_6_hr_exclusive.txt"),
                quote = F, sep = "\t",
                row.names = F, col.names = F)
  }
 
}
