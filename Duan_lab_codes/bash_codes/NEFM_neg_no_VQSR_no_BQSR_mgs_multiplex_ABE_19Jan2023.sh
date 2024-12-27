#!/bin/bash

# Siwei 19 Jan 2023
# Direct count for 20 SNP sites designed for ABE/CBE multiplex editing

# Disable BQSR, use the Phred value directly from CellRanger

# Siwei 08 Aug 2022
# Remove all filters, set minMappingQual=0
# Do not use VQSR recalibration

# Siwei 21 Jul 2022
# update SNP to dbsnp v154
# move all intermediate files to temp dir

# 11 May 2021
# Siwei rewrite in GATK4

# 11 Jun 2020

# source /data/Tools/cellranger-arc-2.0.0/sourceme.bash

vcf_suffix="_12Nov_use_known_genotyping.vcf"
gatk4="/home/zhangs3/Data/Tools/gatk-4.2.6.1/gatk"

ref_path="/home/zhangs3/Data/Databases/Genomes/hg38"
ref_genome="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
dbsnp="/home/zhangs3/Data/Databases/Genomes/hg38/dbsnp_154.hg38.vcf"

# use a /tmp location to host temp files
temp_dir="/home/zhangs3/NVME/package_temp/GATK/hr_0_NEFM_neg_pt2"
working_temp=$temp_dir

mkdir -p vcf_output

date > log.txt

# rm *.vcf
# rm *.vcf.idx
# rm *.recal
# rm *.tranches
# rm *.R
rm -r $temp_dir

for eachfile in *NEFM_neg*WASPed.bam
do

############## locate interval file #######
#        interval=$(echo $eachfile | \
#                sed 's/.*\(CD_[0-9][0-9]\).*/\1/g')
#        interval_file="interval_bed/"$interval"_het_only_from_genotyping.bed"
	interval_file="../../ABE_multiplex_20_SNPs_hg38.bed"

	echo $interval_file

#	rm *.vcf
#	rm *.vcf.idx
#	rm *.recal
#	rm *.tranches
#	rm *.R

	rm -r $temp_dir
	mkdir -p $temp_dir

	echo $eachfile
	echo $eachfile >> log.txt
	date >> log.txt
############## call variants ##############
	
	samtools index -@ 20 $eachfile

#        $gatk4 \
#                BaseRecalibrator \
#                -I $eachfile \
#                -R $ref_genome \
#                --known-sites $dbsnp \
#                --bqsr-baq-gap-open-penalty 30 \
#                -XL $ref_path/hg38_blacklisted_regions.bed \
#		-XL $interval_file \
#                -O $working_temp/BQSR_recal_data.table
#
#        ### apply BQSR tables
#        $gatk4 \
#                ApplyBQSR \
#                -R $ref_genome \
#                -I $eachfile \
#                --bqsr-recal-file $working_temp/BQSR_recal_data.table \
#                -O $working_temp/merged_bam_BQSR_calibrated.bam
#
#	samtools index -@ 20 $working_temp/merged_bam_BQSR_calibrated.bam
	
	$gatk4 --java-options "-Xmx150g" \
		HaplotypeCaller \
		--native-pair-hmm-threads 16 \
		-R $ref_genome \
		--dbsnp $dbsnp \
		-I $eachfile \
		-L $interval_file \
                -XL $ref_path/hg38_blacklisted_regions.bed \
		--output-mode EMIT_ALL_CONFIDENT_SITES \
		-ERC BP_RESOLUTION \
		-O $temp_dir/raw_variants.vcf
	
#		--minimum-mapping-quality 30 \

############variants filtering##########
############SNP#########################
############Remove previous calibration files######


	$gatk4 \
		SelectVariants \
		-R $ref_genome \
		-V $temp_dir/raw_variants.vcf \
		-L $ref_path/autosomes.list \
		--exclude-filtered true \
		-O vcf_output/${eachfile/%.bam/}$vcf_suffix

#	cat $temp_dir/header.txt $temp_dir/body.txt \
#		> vcf_output/${eachfile/%.bam/}$vcf_suffix

	rm -r $temp_dir

done

