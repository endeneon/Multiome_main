#!/bin/bash

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

vcf_suffix="_28Dec_direct_counting.vcf"
gatk4="/home/zhangs3/Data/Tools/gatk-4.2.6.1/gatk"
picard_2254="/home/zhangs3/Data/Tools/gatk-4.1.8.1/picard2_25_4.jar"


extract_bam_coverage="/home/zhangs3/Data/Tools/jvarkit/dist/biostar78285.jar"

ref_path="/home/zhangs3/Data/Databases/Genomes/hg38"
ref_genome="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
dbsnp="/home/zhangs3/Data/Databases/Genomes/hg38/dbsnp_154.hg38.vcf"
blacklisted_regions="/home/zhangs3/Data/Databases/Genomes/hg38/hg38_blacklisted_regions.bed.gz"


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
        interval=$(echo $eachfile | \
                sed 's/.*\(CD_[0-9][0-9]\).*/\1/g')
	interval_file="../../interval_bed/"$interval"_het_only_from_genotyping.bed.bgzip"
	interval_file_bed="../../interval_bed/"$interval"_het_only_from_genotyping.bed"
#  interval_file="../../mgs_608001_all_SNP_sites_hg38_sorted.bed"

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
	date
		
	samtools index -@ 20 $eachfile
#
#        $gatk4 \
#                BaseRecalibrator \
#                -I $eachfile \
#                -R $ref_genome \
#                --known-sites $dbsnp \
#                --bqsr-baq-gap-open-penalty 30 \
#                -XL $ref_path/hg38_blacklisted_regions.bed \
#		-XL $interval_file_bed \
#                -O $working_temp/BQSR_recal_data.table
#
#        ### apply BQSR tables
#        $gatk4 \
#                ApplyBQSR \
#                -R $ref_genome \
#                -I $eachfile \
#                --bqsr-recal-file $working_temp/BQSR_recal_data.table \
#                -O ${eachfile/%.bam/_BQSRed.bam}

#	samtools index -@ 20 ${eachfile/%.bam/_BQSRed.bam}

	bcftools mpileup \
		-O v \
		-f $ref_genome \
		-R $interval_file \
		-q 20 \
		-B \
		-Q 1 \
		-d 2000 \
		-I \
		-o $temp_dir/bcftools_mpileup_output.vcf \
		--threads 16 \
		$eachfile

	### sort vcf
	java -jar $picard_2254 \
		SortVcf \
		-R $ref_genome \
		-I $temp_dir/bcftools_mpileup_output.vcf \
		-O $temp_dir/bcftools_mpileup_output_sorted.vcf.gz


#	tabix \
#		-f \
#		-p vcf \
#		$temp_dir/bcftools_mpileup_output_sorted.vcf.gz

	bcftools call \
		-O z \
		-c \
		-R $interval_file \
		-V indels \
		--threads 16 \
		-o $temp_dir/bcftools_called.vcf.gz \
		$temp_dir/bcftools_mpileup_output_sorted.vcf.gz

	$gatk4 \
		IndexFeatureFile \
		-I $temp_dir/bcftools_called.vcf.gz

#	bcftools filter \
#		-i 'TYPE="snp"' \
#		-O z \
#		-o $temp_dir/bcftools_called_filtered.vcf.gz \
#		$temp_dir/bcftools_called.vcf.gz
#
#	$gatk4 \
#		IndexFeatureFile \
#		-I $temp_dir/bcftools_called_filtered.vcf.gz

	bcftools view \
		-O v \
		-T ^$blacklisted_regions \
		-o vcf_output/${eachfile/%.bam/}$vcf_suffix \
		--threads 16 \
		$temp_dir/bcftools_called.vcf.gz

	
#                -o vcf_output/${eachfile/%.bam/}$vcf_suffix \

	rm -r $temp_dir

done

