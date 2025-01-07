time=$1
cell_type=$2
i=$3


plink_prefix_path=data/plink_bed/CIRM_v1.2_MGS_merged_100cells_hg38_rsq0.7_maf0.05
expression_bed=data/count_matrix/TMM/${time}__${cell_type}.bed.gz
prefix=output/caQTL_mapping/Cis_${time}__${cell_type}__Geno3_Pheno${i}_25kb
covariates_file=data/covariate/TMM/${cell_type}/${time}/Geno3_Pheno${i}_PCA.txt

# run tensorQTL with permutations, report top association for each peak
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
--mode cis --seed 123456 \
--covariates ${covariates_file} \
--window 25000

# run tensorQTL without, report associations for all cis SNPs
nominal_prefix=output/caQTL_mapping/Cis_Nominal_${time}__${cell_type}__Geno3_Pheno${i}_25kb

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${nominal_prefix} \
	--mode cis_nominal \
	--covariates ${covariates_file} \
	 --window 25000
