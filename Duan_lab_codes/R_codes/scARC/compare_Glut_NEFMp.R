## Siwei 10 Feb 2022
## compare Glut with GABA


##### NEFMp_CUXm #####
## make 0 hr_NEFMp group
df_1hr_NEFMp <- read_delim("raw_import/NEFMp__1hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_1hr_NEFMp <- df_1hr_NEFMp[!is.na(df_1hr_NEFMp$REF), ]

df_1hr_NEFMp$REF_C <- apply(df_1hr_NEFMp, 1, 
                    function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12) + 2]))))
df_1hr_NEFMp$ALT_C <- apply(df_1hr_NEFMp, 1, 
                    function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12) + 3]))))

df_1hr_NEFMp <- df_1hr_NEFMp[df_1hr_NEFMp$REF_C > 1, ]
df_1hr_NEFMp <- df_1hr_NEFMp[df_1hr_NEFMp$ALT_C > 1, ]
df_1hr_NEFMp$DP_C <- df_1hr_NEFMp$REF_C + df_1hr_NEFMp$ALT_C

df_1hr_NEFMp <- df_1hr_NEFMp[df_1hr_NEFMp$DP_C > 19, ]


df_1hr_NEFMp$pVal <- apply(df_1hr_NEFMp, 1, 
                     function(x)(binom.test(x = as.numeric(x[16]),
                                            n = sum(as.numeric(x[16]), 
                                                    as.numeric(x[17])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_1hr_NEFMp$FDR <- p.adjust(p = df_1hr_NEFMp$pVal,
                       method = "fdr")

sum(df_1hr_NEFMp$FDR < 0.05)

##### NEFMp_CUXm #####
## make 6 hr_NEFMp group
df_6hr_NEFMp <- read_delim("raw_import/NEFMp__6hr.table", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

df_6hr_NEFMp <- df_6hr_NEFMp[!is.na(df_6hr_NEFMp$REF), ]

df_6hr_NEFMp$REF_C <- apply(df_6hr_NEFMp, 1, 
                            function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12) + 2]))))
df_6hr_NEFMp$ALT_C <- apply(df_6hr_NEFMp, 1, 
                            function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12) + 3]))))

df_6hr_NEFMp <- df_6hr_NEFMp[df_6hr_NEFMp$REF_C > 1, ]
df_6hr_NEFMp <- df_6hr_NEFMp[df_6hr_NEFMp$ALT_C > 1, ]
df_6hr_NEFMp$DP_C <- df_6hr_NEFMp$REF_C + df_6hr_NEFMp$ALT_C

df_6hr_NEFMp <- df_6hr_NEFMp[df_6hr_NEFMp$DP_C > 19, ]


df_6hr_NEFMp$pVal <- apply(df_6hr_NEFMp, 1, 
                           function(x)(binom.test(x = as.numeric(x[17]),
                                                  n = sum(as.numeric(x[17]), 
                                                          as.numeric(x[16])),
                                                  p = 0.5,
                                                  alternative = "t")$p.value))
df_6hr_NEFMp$FDR <- p.adjust(p = df_6hr_NEFMp$pVal,
                             method = "fdr")


sum(df_6hr_NEFMp$FDR < 0.05)

sum(df_0hr$ID[df_0hr$FDR < 0.05] %in% 
      df_0hr_NEFMp$ID[df_0hr_NEFMp$FDR < 0.05]) # 11
length(df_0hr$ID[df_0hr$FDR < 0.05] %in% 
      df_0hr_NEFMp$ID[df_0hr_NEFMp$FDR < 0.05]) # 11
length(df_0hr$ID[df_0hr$FDR < 0.05])
length(df_0hr_NEFMp$ID[df_0hr_NEFMp$FDR < 0.05])

df_0hr$ID[df_0hr$FDR < 0.05][df_0hr$ID[df_0hr$FDR < 0.05] %in% 
             df_0hr_NEFMp$ID[df_0hr_NEFMp$FDR < 0.05]]

(df_0hr$ID[df_0hr$FDR < 0.05] %in% 
    df_0hr_NEFMp$ID[df_0hr_NEFMp$FDR < 0.05])


sum(df_1hr$ID[df_1hr$FDR < 0.05] %in% 
      df_1hr_NEFMp$ID[df_1hr_NEFMp$FDR < 0.05]) # 0

sum(df_6hr$ID[df_6hr$FDR < 0.05] %in% 
      df_6hr_NEFMp$ID[df_6hr_NEFMp$FDR < 0.05]) # 6

df_6hr$ID[df_6hr$FDR < 0.05][df_6hr$ID[df_6hr$FDR < 0.05] %in% 
      df_6hr_NEFMp$ID[df_6hr_NEFMp$FDR < 0.05]] # 6

sum(df_0hr_NEFMp$ID[df_0hr_NEFMp$FDR < 0.05] %in% 
      df_1hr_NEFMp$ID[df_1hr_NEFMp$FDR < 0.05]) # 8
df_0hr_NEFMp$ID[df_0hr_NEFMp$FDR < 0.05][(df_0hr_NEFMp$ID[df_0hr_NEFMp$FDR < 0.05] %in% 
    df_1hr_NEFMp$ID[df_1hr_NEFMp$FDR < 0.05])]

sum(df_0hr_NEFMp$ID[df_0hr_NEFMp$FDR < 0.05] %in% 
      df_6hr_NEFMp$ID[df_6hr_NEFMp$FDR < 0.05]) # 7
df_0hr_NEFMp$ID[df_0hr_NEFMp$FDR < 0.05][(df_0hr_NEFMp$ID[df_0hr_NEFMp$FDR < 0.05] %in% 
    df_6hr_NEFMp$ID[df_6hr_NEFMp$FDR < 0.05])]

sum(df_1hr_NEFMp$ID[df_1hr_NEFMp$FDR < 0.05] %in% 
      df_6hr_NEFMp$ID[df_6hr_NEFMp$FDR < 0.05]) # 9
df_1hr_NEFMp$ID[df_1hr_NEFMp$FDR < 0.05][(df_1hr_NEFMp$ID[df_1hr_NEFMp$FDR < 0.05] %in% 
    df_6hr_NEFMp$ID[df_6hr_NEFMp$FDR < 0.05])]

sum(df_0hr_NEFMp$ID[df_0hr_NEFMp$FDR < 0.05] %in% 
      df_1hr_NEFMp$ID[df_1hr_NEFMp$FDR < 0.05] %in%
      df_6hr_NEFMp$ID[df_6hr_NEFMp$FDR < 0.05]) # 0

max(df_0hr$DP)
max(df_1hr$DP)
max(df_6hr$DP)

max(df_0hr_NEFMp$DP_C)
max(df_1hr_NEFMp$DP_C)
max(df_6hr_NEFMp$DP_C)

median(df_0hr_NEFMp$DP_C)
median(df_1hr_NEFMp$DP_C)
median(df_6hr_NEFMp$DP_C)

temp_writeout_table <- df_0hr_NEFMp
temp_writeout_table <- df_1hr_NEFMp
temp_writeout_table <- df_6hr_NEFMp
output_table <- data.frame(cbind(temp_writeout_table$CHROM,
                                 temp_writeout_table$POS,
                                 temp_writeout_table$POS,
                                 temp_writeout_table$REF,
                                 temp_writeout_table$ALT,
                                 temp_writeout_table$ID,
                                 temp_writeout_table$REF_C,
                                 temp_writeout_table$ALT_C,
                                 temp_writeout_table$DP_C,
                                 temp_writeout_table$pVal,
                                 temp_writeout_table$FDR))

write.table(output_table,
            file = "txt_4_annovar/NEFMp_1hr.avinput",
            quote = F,
            sep = "\t", 
            col.names = F, row.names = F)


##### NEMFm_CUXp #####
## make 0 hr_NEMFm group
df_0hr_NEMFm <- read_delim("raw_import/NEMFm__0hr.table", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

df_0hr_NEMFm <- df_0hr_NEMFm[!is.na(df_0hr_NEMFm$REF), ]

df_0hr_NEMFm$REF_C <- apply(df_0hr_NEMFm, 1, 
                            function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12) + 2]))))
df_0hr_NEMFm$ALT_C <- apply(df_0hr_NEMFm, 1, 
                            function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12) + 3]))))

df_0hr_NEMFm <- df_0hr_NEMFm[df_0hr_NEMFm$REF_C > 1, ]
df_0hr_NEMFm <- df_0hr_NEMFm[df_0hr_NEMFm$ALT_C > 1, ]
df_0hr_NEMFm$DP_C <- df_0hr_NEMFm$REF_C + df_0hr_NEMFm$ALT_C

sum(df_0hr_NEMFm$DP_C > 9)
df_0hr_NEMFm <- df_0hr_NEMFm[df_0hr_NEMFm$DP_C > 9, ]


df_0hr_NEMFm$pVal <- apply(df_0hr_NEMFm, 1, 
                           function(x)(binom.test(x = as.numeric(x[16]),
                                                  n = sum(as.numeric(x[16]), 
                                                          as.numeric(x[17])),
                                                  p = 0.5,
                                                  alternative = "t")$p.value))
df_0hr_NEMFm$FDR <- p.adjust(p = df_0hr_NEMFm$pVal,
                             method = "fdr")

sum(df_0hr_NEMFm$FDR < 0.05)


###
df_1hr_NEMFm <- read_delim("raw_import/NEMFm__1hr.table", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

df_1hr_NEMFm <- df_1hr_NEMFm[!is.na(df_1hr_NEMFm$REF), ]

df_1hr_NEMFm$REF_C <- apply(df_1hr_NEMFm, 1, 
                            function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12) + 2]))))
df_1hr_NEMFm$ALT_C <- apply(df_1hr_NEMFm, 1, 
                            function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12) + 3]))))

df_1hr_NEMFm <- df_1hr_NEMFm[df_1hr_NEMFm$REF_C > 1, ]
df_1hr_NEMFm <- df_1hr_NEMFm[df_1hr_NEMFm$ALT_C > 1, ]
df_1hr_NEMFm$DP_C <- df_1hr_NEMFm$REF_C + df_1hr_NEMFm$ALT_C

df_1hr_NEMFm <- df_1hr_NEMFm[df_1hr_NEMFm$DP_C > 9, ]


df_1hr_NEMFm$pVal <- apply(df_1hr_NEMFm, 1, 
                           function(x)(binom.test(x = as.numeric(x[16]),
                                                  n = sum(as.numeric(x[16]), 
                                                          as.numeric(x[17])),
                                                  p = 0.5,
                                                  alternative = "t")$p.value))
df_1hr_NEMFm$FDR <- p.adjust(p = df_1hr_NEMFm$pVal,
                             method = "fdr")

sum(df_1hr_NEMFm$FDR < 0.05)

##### NEFMm_CUXp #####
## make 6 hr_NEMFm group
df_6hr_NEMFm <- read_delim("raw_import/NEMFm__6hr.table", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

df_6hr_NEMFm <- df_6hr_NEMFm[!is.na(df_6hr_NEMFm$REF), ]

df_6hr_NEMFm$REF_C <- apply(df_6hr_NEMFm, 1, 
                            function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12) + 2]))))
df_6hr_NEMFm$ALT_C <- apply(df_6hr_NEMFm, 1, 
                            function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12) + 3]))))

df_6hr_NEMFm <- df_6hr_NEMFm[df_6hr_NEMFm$REF_C > 1, ]
df_6hr_NEMFm <- df_6hr_NEMFm[df_6hr_NEMFm$ALT_C > 1, ]
df_6hr_NEMFm$DP_C <- df_6hr_NEMFm$REF_C + df_6hr_NEMFm$ALT_C

df_6hr_NEMFm <- df_6hr_NEMFm[df_6hr_NEMFm$DP_C > 9, ]


df_6hr_NEMFm$pVal <- apply(df_6hr_NEMFm, 1, 
                           function(x)(binom.test(x = as.numeric(x[17]),
                                                  n = sum(as.numeric(x[17]), 
                                                          as.numeric(x[16])),
                                                  p = 0.5,
                                                  alternative = "t")$p.value))
df_6hr_NEMFm$FDR <- p.adjust(p = df_6hr_NEMFm$pVal,
                             method = "fdr")


sum(df_6hr_NEMFm$FDR < 0.05)

sum(df_0hr$ID[df_0hr$FDR < 0.05] %in% 
      df_0hr_NEMFm$ID[df_0hr_NEMFm$FDR < 0.05]) # 11

sum(df_1hr$ID[df_1hr$FDR < 0.05] %in% 
      df_1hr_NEMFm$ID[df_1hr_NEMFm$FDR < 0.05]) # 0

sum(df_6hr$ID[df_6hr$FDR < 0.05] %in% 
      df_6hr_NEMFm$ID[df_6hr_NEMFm$FDR < 0.05]) # 6

sum(df_0hr_NEMFm$ID[df_0hr_NEMFm$FDR < 0.05] %in% 
      df_1hr_NEMFm$ID[df_1hr_NEMFm$FDR < 0.05]) # 8

sum(df_0hr_NEMFm$ID[df_0hr_NEMFm$FDR < 0.05] %in% 
      df_6hr_NEMFm$ID[df_6hr_NEMFm$FDR < 0.05]) # 7

sum(df_1hr_NEMFm$ID[df_1hr_NEMFm$FDR < 0.05] %in% 
      df_6hr_NEMFm$ID[df_6hr_NEMFm$FDR < 0.05]) # 9

sum(df_0hr_NEMFm$ID[df_0hr_NEMFm$FDR < 0.05] %in% 
      df_1hr_NEMFm$ID[df_1hr_NEMFm$FDR < 0.05] %in%
      df_6hr_NEMFm$ID[df_6hr_NEMFm$FDR < 0.05]) # 0

max(df_0hr_NEMFm$DP_C)
max(df_1hr_NEMFm$DP_C)
max(df_6hr_NEMFm$DP_C)

median(df_0hr_NEMFm$DP_C)
median(df_1hr_NEMFm$DP_C)
median(df_6hr_NEMFm$DP_C)

temp_writeout_table <- df_0hr_NEMFm
temp_writeout_table <- df_1hr_NEMFm
temp_writeout_table <- df_6hr_NEMFm
output_table <- data.frame(cbind(temp_writeout_table$CHROM,
                                 temp_writeout_table$POS,
                                 temp_writeout_table$POS,
                                 temp_writeout_table$REF,
                                 temp_writeout_table$ALT,
                                 temp_writeout_table$ID,
                                 temp_writeout_table$REF_C,
                                 temp_writeout_table$ALT_C,
                                 temp_writeout_table$DP_C,
                                 temp_writeout_table$pVal,
                                 temp_writeout_table$FDR))
# output_table <- output_table[order(output_table$X8)]

write.table(output_table,
            file = "txt_4_annovar/NEFMm_1hr.avinput",
            quote = F,
            sep = "\t", 
            col.names = F, row.names = F)
