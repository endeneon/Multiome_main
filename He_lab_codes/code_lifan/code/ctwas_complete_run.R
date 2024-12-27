library(ctwas)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

if (length(args) <7) {
  stop(" 7 arguments:
       * zscore file name
       * ld R directory
       * weight
       * config file name (.R)
       * out expr z file name
       * out file name
       * outputdir", call.=FALSE)
}

print(args)

outputdir <- args[7]

dir.create(outputdir, showWarnings=F, recursive=T)

z_snp_stem <- unlist(strsplit(rev(unlist(strsplit(args[1], "/")))[1],"[.]"))[1]
z_snp_outfile <- paste0(outputdir, "/", z_snp_stem, ".RDS")

if (file.exists(z_snp_outfile)){
  z_snp <- readRDS(z_snp_outfile)
} else {
  z_snp <- VariantAnnotation::readVcf(args[1])
  z_snp <- as.data.frame(gwasvcf::vcf_to_tibble(z_snp))
  z_snp$Z <- z_snp$ES/z_snp$SE
  z_snp <- z_snp[,c("rsid", "ALT", "REF", "Z", "SS")]
  colnames(z_snp) <- c("id", "A1", "A2", "z", "ss")
  z_snp <- z_snp[!(z_snp$id %in% z_snp$id[duplicated(z_snp$id)]),] #drop multiallelic variants (id not unique)
  saveRDS(z_snp, file=z_snp_outfile)
}

ld_R_dir <- args[2]
weight <- args[3]
weight <- unlist(strsplit(weight, ";"))

outname.e <- args[5]
outname <- args[6]

source(args[4]) # config

#preharmonize snp z score
if (file.exists(paste0(outputdir, "/", outname, "_z_snp.Rd"))){
  load(file = paste0(outputdir, "/", outname, "_z_snp.Rd"))
} else {
  res <- preharmonize_z_ld(z_snp=z_snp, 
                   ld_R_dir=ld_R_dir, 
                   outputdir=outputdir,
                   outname=outname,
                   harmonize_z=T, 
                   strand_ambig_action_z="none")
  z_snp <- res$z_snp
  save(z_snp, file = paste0(outputdir, "/", outname, "_z_snp.Rd"))
  rm(res)
}

#impute gene z-scores for both sets of prediction weights by chromosome
for (i in 1:22){
  if (!file.exists(paste0(outputdir, "/", outname.e, "_chr", i, ".exprqc.Rd"))){
    res <- impute_expr_z(z_snp = z_snp,
                         weight = weight,
                         ld_R_dir = ld_R_dir,
                         outputdir = outputdir,
                         outname = outname.e,
                         harmonize_z = F,
                         harmonize_wgt = F,
                         strand_ambig_action_wgt="none",
                         ncore=8, 
                         chrom=i,
                         scale_by_ld_variance=T)
  }
}

thin <- 0.1
max_snp_region <- 20000
ncore <- 6

#additional options for the multigroup branch
fine_map <- F
ncore_LDR <- 6

#prior variance structure
group_prior_var_structure <- "shared"

#combine the imputed gene z-scores
ld_exprfs <- paste0(outputdir, "/", outname, "_chr", 1:22, ".expr.gz")
z_gene <- list()
for (i in 1:22){
  load(paste0(outputdir, "/", outname, "_chr", i, ".exprqc.Rd"))
  z_gene[[i]] <- z_gene_chr
}
z_gene <- do.call(rbind, z_gene)
rownames(z_gene) <- NULL
save(z_gene, file = paste0(outputdir, "/", outname, "_z_gene.Rd"))

z_gene$type <- sapply(z_gene$id, function(x){paste(unlist(strsplit(unlist(strsplit(x, "[|]"))[2],"_")), collapse="_") })
#z_gene$QTLtype <- sapply(z_gene$id, function(x){paste(unlist(strsplit(unlist(strsplit(x, "[|]"))[2],"_"))[2], collapse="_") })

## Estimate group prior
if (file.exists(paste0(outputdir, outname, ".s2.susieIrssres.Rd"))){
  load(paste0(outputdir, outname, ".s2.susieIrssres.Rd"))
  
  group_prior <- group_prior_rec[,ncol(group_prior_rec)]
  #group_prior["SNP"] <- group_prior["SNP"]*thin # Obsolete
  group_prior_var <- group_prior_var_rec[,ncol(group_prior_var_rec)]
  rm(group_prior_rec, group_prior_var_rec)
} else {
  parameters <- ctwas_rss(z_gene = z_gene, 
                          z_snp = z_snp, 
                          ld_exprfs = ld_exprfs, 
                          ld_R_dir = ld_R_dir, 
                          ld_regions = "EUR",
                          ld_regions_version = "b38",
                          outputdir = outputdir, 
                          outname = outname,
                          thin = thin,
                          max_snp_region = max_snp_region,
                          ncore = ncore,
                          ncore_LDR = ncore_LDR,
                          group_prior_var_structure = "shared",
                          merge = FALSE, # Only top 1 QTL was used
                          fine_map = F)
  group_prior <- parameters$group_prior
  #group_prior["SNP"] <- group_prior["SNP"]*thin # Obsolete
  
  group_prior_var <- parameters$group_prior_var
}

## Run cTWAS
ctwas_rss(z_gene = z_gene, 
          z_snp = z_snp, 
          ld_exprfs = ld_exprfs, 
          ld_R_dir = ld_R_dir, 
          ld_regions = "EUR",
          ld_regions_version = "b38",
          outputdir = outputdir, 
          outname = paste0(outname, "_cover9"),
          thin = thin,
          max_snp_region = max_snp_region,
          ncore = ncore,
          ncore_LDR = ncore_LDR,
          group_prior = group_prior, 
          group_prior_var = group_prior_var, 
          estimate_group_prior = F, 
          estimate_group_prior_var = F,
          rerun_gene_PIP=0.2,
          merge=F, # We only used 1 top QTL
          L=5,
          fine_map = TRUE,
          coverage=0.9)
