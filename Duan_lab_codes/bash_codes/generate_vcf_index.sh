#!/bin/bash


# 24 May 2022
# 07 May 2021 Siwei

## use conda aligners environment

## Set these environment vars to point to
## your local installation

gatk4_path="/home/zhangs3/Data/Tools/gatk-4.1.8.1"

jre_8="/home/zhangs3/Data/Tools/jre1.8.0_291/bin"
openjdk_16="/home/zhangs3/Data/Tools/jdk-16.0.1/bin/java"

hg38_index="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"

for eachfile in *.vcf
do
	
	$gatk4_path/gatk \
		IndexFeatureFile \
        	-I $eachfile 
done
