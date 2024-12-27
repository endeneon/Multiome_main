library(EnsDb.Hsapiens.v86)
library(biomaRt)

edb <- EnsDb.Hsapiens.v86
ensgenes <- transcriptsBy(edb, by = "gene", filter = SeqNameFilter(as.character(1:22)))
ensembl_hs38 <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
View(listAttributes(ensembl_hs38))
transcript_lookup_table <- getBM(attributes = c('ensembl_transcript_id', 
                                                'external_gene_name', 
                                                'gene_biotype'),
                                 mart = ensembl_hs38)

query_df <- data.frame(transcript_id = ens_use@elementMetadata$tx_id)

assembled_df <- merge(query_df,
                      transcript_lookup_table,
                      by.x = "transcript_id",
                      by.y = "ensembl_transcript_id",
                      all.x = T)

sum(is.na(assembled_df$external_gene_name))
sorted_df <- assembled_df[order(match(assembled_df$transcript_id, query_df$transcript_id)), ]

ens_use <- ensgenes@unlistData
# ens_use@elementMetadata[[7]] <- NULL
ens_use@elementMetadata$gene_name <- sorted_df$external_gene_name
ens_use@elementMetadata$gene_biotype <- sorted_df$gene_biotype
ens_use@elementMetadata$type <- "gene"

# change chromosome name
ens_use@seqnames@values <- paste0("chr", ens_use@seqnames@values)
ens_use@seqinfo@genome <- rep("hg38", 22)
sum(is.na(ens_use@elementMetadata$gene_biotype))
ens_use@elementMetadata$gene_biotype[is.na(ens_use@elementMetadata$gene_biotype)] <- "unknown"

# test
ens_use$gene_name[100]


save("ens_use", file = "EnsDb_UCSC_hg38_annotation_new.RData")
