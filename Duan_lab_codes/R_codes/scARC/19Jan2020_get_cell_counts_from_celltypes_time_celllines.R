# Chuxuan Li 01/19/2022
# count number of cells in each cell type at each time point from each cell type

library(Seurat)

# remove NAs
integrated_renamed_1$cell.line.ident[is.na(integrated_renamed_1$cell.line.ident)] <- "unmatched"

# get lists of variable values
types <- unique(integrated_renamed_1$fine.cell.type)
times <- unique(integrated_renamed_1$time.ident.mod)
lines <- unique(integrated_renamed_1$cell.line.ident)
total.len <- length(types) * length(times) * length(lines)

df <- data.frame(cell.type = rep("", total.len),
                 time.point = rep("", total.len),
                 cell.line = rep("", total.len),
                 counts = rep(0, total.len))


i = 0
for (c in types){
  print(paste0("cell type: ", c))
  for (t in times){
    print(paste0("time: ", t))
    for (l in lines){
      print(paste0("cell line: ", l))
      i = i + 1
      count <- sum(integrated_renamed_1$fine.cell.type == c &
                     integrated_renamed_1$time.ident.mod == t &
                     integrated_renamed_1$cell.line.ident == l, na.rm = T)
      print(count)
      df[i, 1] <- c
      df[i, 2] <- t
      df[i, 3] <- l
      df[i, 4] <- count
    }
  }
}

# output df
write.table(df, file = "RNAseq_cell_count.csv",
            quote = F, sep = ",", col.names = c("cell.type", "time.point", "cell.line", "counts"))

# plot it
ggplot(data = df, aes(x = cell.line, 
                      y = cell.type,
                      fill = time.point,
                      size = counts,
                      group = time.point)) +
  scale_fill_manual(values = brewer.pal(3, "Set2")) +
  geom_point(shape = 21,
             na.rm = TRUE,
             position = position_dodge(width = 0.9)) +
  scale_size_area() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  coord_flip()
