#!/bin/bash

# Siwei 05 Jul 2022
# Use the GATK pipeline to convert 10x subset-bam generated per-line/time/type
# bam files for
# 1) realign with bwa
# 2) calibrate by WASP
# and prepare for vcf calls

# use dbSNP v154

vcf_suffix="_scATAC.vcf"

jre_18="/home/zhangs3/Data/Tools/jre1.8.0_291/bin/java"

gatk4="/home/zhangs3/Data/Tools/gatk-4.1.8.1/gatk"
picard="/home/zhangs3/Data/Tools/gatk-4.1.8.1/picard2_25_4.jar"

ref_path="/home/zhangs3/Data/Databases/Genomes/hg38"
ref_genome="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
dbsnp="/home/zhangs3/Data/Databases/Genomes/hg38/dbsnp_154.hg38.vcf"

RevertSam_temp="/home/zhangs3/NVME/package_temp/RevertSam_temp"

mkdir -p vcf_output
mkdir -p $RevertSam_temp

date > log.txt

rm *.vcf
rm *.vcf.idx
rm *.recal
rm *.tranches
rm *.R


for eachfile in *.bam
do
        rm *.vcf
        rm *.vcf.idx
        rm *.recal
        rm *.tranches
        rm *.R

	$jre_18 -Xmx100G -jar \
		$picard \
		RevertSam \
		I=$eachfile \
		O=ubam_revertsam.bam \
		TMP_DIR=$RevertSam_temp

	

