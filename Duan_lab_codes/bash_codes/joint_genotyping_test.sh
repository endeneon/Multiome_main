#!/bin/bash

# Siwei 19 Jan 2023
# Direct count for 20 SNP sites designed for ABE/CBE multiplex editing
# Use CombineGVCFs and GenotypeGVCFs to jointly genotype all gvcf files produced from the previous step

vcf_suffix="_12Nov_use_known_genotyping.g.vcf"
gatk4="/home/zhangs3/Data/Tools/gatk-4.2.6.1/gatk"

ref_path="/home/zhangs3/Data/Databases/Genomes/hg38"
ref_genome="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
dbsnp="/home/zhangs3/Data/Databases/Genomes/hg38/dbsnp_154.hg38.vcf"


mkdir -p output

ls *.g.vcf > gvcf_input.list

$gatk4 \
	CombineGVCFs \
	-R $ref_genome \
	-D $dbsnp \
	--call-genotypes true \
	-V gvcf_input.list \
	-O output/cohort.g.vcf

$gatk4 \
	GenotypeGVCFs \
	-R $ref_genome \
	-D $dbsnp \
	-V output/cohort.g.vcf \
	--call-genotypes true \
	-all-sites true \
	-O output/genotyped.vcf


