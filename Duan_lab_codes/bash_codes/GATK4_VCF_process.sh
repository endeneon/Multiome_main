#!/bin/bash

# use for ghiprd03

# 05 Jul 2022
# Siwei upgrade to dbsnp154

# 11 May 2021
# Siwei rewrite in GATK4

# 11 Jun 2020

vcf_suffix="_scATAC.vcf"
gatk4="/data/zhangs3/software/gatk-4.1.8.1/gatk"

ref_path="/data/zhangs3/Databases"
ref_genome="/data/zhangs3/Databases/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
dbsnp="/data/zhangs3/Databases/dbsnp_154.hg38.vcf"

mkdir -p vcf_output

date > log.txt

rm *.vcf
rm *.vcf.idx
rm *.recal
rm *.tranches
rm *.R


for EACHFILE in *.bam
do
	rm *.vcf
	rm *.vcf.idx
	rm *.recal
	rm *.tranches
	rm *.R


	echo $EACHFILE
	echo $EACHFILE >> log.txt
	date >> log.txt
############## call variants ##############
	
	samtools index -@ 20 $EACHFILE
	
	$gatk4 --java-options "-Xmx100g" \
		HaplotypeCaller \
		--native-pair-hmm-threads 36 \
		-R $ref_genome \
		--dbsnp $dbsnp \
		-I $EACHFILE \
		-O raw_variants.vcf
	
#		--minimum-mapping-quality 30 \

############variants filtering##########
############SNP#########################
############Remove previous calibration files######

############ SNP Recalibration ###########

        $gatk4 --java-options "-Xmx100g" \
	 	VariantRecalibrator \
		-R $ref_genome \
		-V raw_variants.vcf \
		--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $ref_path/hapmap_3.3.hg38.vcf \
		--resource:omni,known=false,training=true,truth=true,prior=12.0 $ref_path/1000G_omni2.5.hg38.vcf \
		--resource:1000G,known=false,training=true,truth=false,prior=10.0 $ref_path/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf \
		--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
		--max-gaussians 2 \
		-an DP -an QD -an FS -an SOR -an MQ -an ReadPosRankSum \
		-mode SNP \
		-AS \
		-tranche 100.0 -tranche 99.0 -tranche 95.0 -tranche 50.0 \
		-O raw_recalibrate_SNP.recal \
		--tranches-file raw_recalibrate_SNP.tranches \
		--rscript-file raw_recalibrate_SNP_plots.R  

########### apply SNP recalib###########

        $gatk4 --java-options "-Xmx100g" \
		ApplyVQSR \
		-R $ref_genome \
		-V raw_variants.vcf \
		-mode SNP \
		-AS \
		--truth-sensitivity-filter-level 99.0 \
		--recal-file raw_recalibrate_SNP.recal \
		--tranches-file raw_recalibrate_SNP.tranches \
		-O variants_calibrated.vcf

########## filter SNPs located in the hg38 blacklisted region #######

	$gatk4 --java-options "-Xmx100g" \
		VariantFiltration \
		-R $ref_genome \
		-V variants_calibrated.vcf \
		--mask $ref_path/hg38_blacklisted_regions.bed --mask-name "blacklisted_regions" \
		-O variants_calibrated_blacklisted.vcf

########## select SNP for output from autosomes only #######

	$gatk4 --java-options "-Xmx100g" \
		SelectVariants \
		-R $ref_genome \
		-V variants_calibrated_blacklisted.vcf \
		--select-type-to-include SNP \
		--restrict-alleles-to BIALLELIC \
		-L $ref_path/autosomes.list \
		-O vcf_output/${EACHFILE/%.bam/}$vcf_suffix

done

