for (i in 1:length(barcode_list)){
  lines <- unique(barcode_list[[i]]$line)
  #print(lines)
  print(unique(nomouselist$`70305-6`$orig.ident))
  
  nomouselist$`70305-6`$cell.line.ident <- "unmatched"
  for (j in 1:length(lines)){
    line_spec_barcodes <- barcode_list[[i]]$barcode[barcode_list[[i]]$line == lines[j]]
    nomouselist[[i]]$cell.line.ident[nomouselist$`70305-6`@assays$RNA@counts@Dimnames[[2]] %in% 
                                       line_spec_barcodes] <- lines[j]
  }
  print(paste0("\n70305-6 with barcodes from ", names(barcode_list)[i], ": "))
  print(sum(nomouselist$`70305-6`$cell.line.ident != "unmatched"))
}

group_70305_gex_1 <- read_delim("~/Data/FASTQ/Duan_Project_030/barcode_demux_output/gex_output_70305/group_70305_gex.1.best", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)
group_70305_gex_6 <- read_delim("~/Data/FASTQ/Duan_Project_030/barcode_demux_output/gex_output_70305/group_70305_gex.2.best", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)
h5file <- Read10X_h5(filename = h5list[12])
cat("h5: ", str_extract(string = h5list[12],
                        pattern = "[0-9]+-[0-6]"))
obj_70305_6 <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                          project = str_extract(string = h5list[12],
                                                pattern = "[0-9]+-[0-6]"))
h5file <- Read10X_h5(filename = h5list[11])
cat("h5: ", str_extract(string = h5list[11],
                        pattern = "[0-9]+-[0-6]"))
obj_70305_1 <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                                  project = str_extract(string = h5list[11],
                                                        pattern = "[0-9]+-[0-6]"))
sum(obj_70305_1@assays$RNA@counts@Dimnames[[2]] %in% group_70305_gex_1$BARCODE) #209
sum(obj_70305_6@assays$RNA@counts@Dimnames[[2]] %in% group_70305_gex_1$BARCODE) #13020

sum(obj_70305_1@assays$RNA@counts@Dimnames[[2]] %in% group_70305_gex_6$BARCODE) #11256
sum(obj_70305_6@assays$RNA@counts@Dimnames[[2]] %in% group_70305_gex_6$BARCODE) #209

# 70305-0
group_70305_gex_0 <- read_delim("~/Data/FASTQ/Duan_Project_030/barcode_demux_output/gex_output_70305/group_70305_gex.0.best", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)
h5file <- Read10X_h5(filename = h5list[10])
cat("h5: ", str_extract(string = h5list[10],
                        pattern = "[0-9]+-[0-6]"))
obj_70305_0 <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                                  project = str_extract(string = h5list[10],
                                                        pattern = "[0-9]+-[0-6]"))
sum(obj_70305_0@assays$RNA@counts@Dimnames[[2]] %in% group_70305_gex_0$BARCODE) #10123

# use output barcode 
sum(obj_70305_6@assays$RNA@counts@Dimnames[[2]] %in% barcode_list$`70305_1hr`$barcode)
sum(obj_70305_1@assays$RNA@counts@Dimnames[[2]] %in% barcode_list$`70305_6hr`$barcode)


# recount 
for (i in 1:length(objlist)){
  obj <- objlist[[i]]
  cell_counts$total[i] <- ncol(obj)
  cell_counts$mouse[i] <- sum(obj$mouse.ident == "mouse")
  cell_counts$nonmouse[i] <- sum(obj$mouse.ident == "human")
}
