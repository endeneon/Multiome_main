time=$1
cell_type=$2
i=$3


plink_prefix_path=data/plink_bed/CIRM_v1.2_MGS_merged_100cells_hg38_rsq0.7_maf0.05
expression_bed=data/count_matrix/TMM/${time}__${cell_type}.bed.gz
prefix=output/covar_eval/${time}__${cell_type}__Geno3_Pheno${i}_Empirical
covariates_file=data/covariate/TMM/${cell_type}/${time}/Geno3_Pheno${i}_PCA.txt

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
--mode cis --seed 123456 \
--covariates ${covariates_file} \
--window 25000 \
--permutations 1000
