#!/bin/bash


# 24 May 2022
# 07 May 2021 Siwei

## use conda aligners environment

## Set these environment vars to point to
## your local installation

gatk37_path="/home/zhangs3/Data/Tools/GATK37"
gatk36_path="/home/zhangs3/Data/Tools/gatk_36/opt/gatk-3.6"
gatk4_path="/home/zhangs3/Data/Tools/gatk-4.1.8.1"

jre_8="/home/zhangs3/Data/Tools/jre1.8.0_291/bin"
openjdk_16="/home/zhangs3/Data/Tools/jdk-16.0.1/bin/java"


# hg38_index="/home/zhangs3/Data/Databases/Genomes/hg38/GRCh38p13.primary_assembly.genome.fa"
# hg38_index="/home/zhangs3/Data/Databases/Genomes/hg38/hg38_p12_chr_only.fasta"
hg38_index="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"

mkdir -p output

## list all .vcf files
#shopt -s nullglob
#vcf_array=(*.vcf)
#echo "${vcf_array[@]}"

## send all .vcf files to a .list file
ls *.vcf > vcf_list.list

## use picard to merge all vcfs
## note that picard is located under the gatk4 folder
#$openjdk_16 -Xmx200G -jar $gatk4_path/picard2_25_4.jar MergeVcfs \
#	-I vcf_list.list \
#	-O raw_merged_vcfs.vcf

## !! Note!! Use GATK3.6 and jre 1.8 CombineVariants to merge all VCFs
## The new Picard from GATK4 only accepts GVCF files therefore 
## cannot merge vcfs with different sample identities

$jre_8/java -Xmx200G -jar $gatk36_path/GenomeAnalysisTK.jar \
	-T CombineVariants \
	-R $hg38_index \
	--variant vcf_list.list \
	-genotypeMergeOptions UNIQUIFY \
	-o output/raw_merged_vcfs.vcf

## select variants that only includes SNP
$gatk4_path/gatk SelectVariants \
	-R $hg38_index \
	-V output/raw_merged_vcfs.vcf \
	--select-type-to-include SNP \
	--restrict-alleles-to BIALLELIC \
	-XL /home/zhangs3/Data/Databases/Genomes/hg38/hg38_blacklisted_regions.bed \
	-O output/GABA_merged_SNP.vcf

# convert FILTER column to pass
cat output/GABA_merged_SNP.vcf \
	| grep "^#" \
	> output/header.temp
cat output/GABA_merged_SNP.vcf \
	| grep -v "^#" \
	| awk -F '\t' 'BEGIN {OFS=FS}{ $7="PASS" ; print }' \
	> output/body.temp

cat output/header.temp output/body.temp \
	> output/GABA_merged_SNP_24May2022.vcf



## select variants with DP >= 10
#$gatk4_path/gatk VariantFiltration \
#        -R $hg38_index \
#        -V output/GABA_merged_SNP.vcf \
#	--filter-expression "DP < 10.0" \
#        --filter-name "DP_filter" \
#        -O output/GABA_merged_SNP_24May2022.vcf

## export vcfs to table
$gatk4_path/gatk VariantsToTable \
	-V output/GABA_merged_SNP_24May2022.vcf \
	-F CHROM -F POS -F ID -F REF -F ALT -F NCALLED \
	-GF AD \
	-O output/GABA_scATAC_0hr_merged_SNP_20Jun2022.txt

## convert the txt title line
cat output/GABA_scATAC_0hr_merged_SNP_20Jun2022.txt \
	| sed 's/NA/0,0/g' \
	| sed 's/,/\t/g' \
        | sed 's/variant\.AD/variant\tAD/g' \
        | sed 's/variant.\.AD/variant\tAD/g' \
        | sed 's/variant..\.AD/variant\tAD/g' \
	> output/GABA_scATAC_0hr_merged_SNP_20Jun2022_4_R.txt

