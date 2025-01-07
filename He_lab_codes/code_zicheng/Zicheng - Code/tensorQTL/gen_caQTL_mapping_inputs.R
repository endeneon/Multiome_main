caqtl_df = data.frame(time = rep(c("0hr","1hr","6hr"),3),
		                            cell_type = c(rep("npglut",3),rep("GABA",3),rep("nmglut",3)))

caqtl_df$ncovar = c(4,5,8,8,6,7,5,9,6)
write.table(caqtl_df,"caQTL_mapping_inputs.txt",col.names = F,quote = F,sep = " ",row.names = F)

