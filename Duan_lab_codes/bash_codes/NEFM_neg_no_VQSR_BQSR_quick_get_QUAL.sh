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

vcf_suffix="_19Dec_QUAL.vcf"
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

for eachfile in 0hr_CD_60_NEFM_neg*.bam
do

############## locate interval file #######
        interval=$(echo $eachfile | \
                sed 's/.*\(CD_[0-9][0-9]\).*/\1/g')
        interval_file="interval_bed/"$interval"_het_only_from_genotyping.bed"
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

        $gatk4 \
                BaseRecalibrator \
                -I $eachfile \
                -R $ref_genome \
                --known-sites $dbsnp \
                --bqsr-baq-gap-open-penalty 30 \
                -XL $ref_path/hg38_blacklisted_regions.bed \
		-XL $interval_file \
                -O $working_temp/BQSR_recal_data.table

        ### apply BQSR tables
        $gatk4 \
                ApplyBQSR \
                -R $ref_genome \
                -I $eachfile \
                --bqsr-recal-file $working_temp/BQSR_recal_data.table \
                -O $working_temp/merged_bam_BQSR_calibrated.bam

	samtools index -@ 20 $working_temp/merged_bam_BQSR_calibrated.bam
	
	$gatk4 --java-options "-Xmx150g" \
		HaplotypeCaller \
		--native-pair-hmm-threads 24 \
		-R $ref_genome \
		--dbsnp $dbsnp \
		-I $working_temp/merged_bam_BQSR_calibrated.bam \
                -XL $ref_path/hg38_blacklisted_regions.bed \
		-O $temp_dir/raw_variants.vcf
	
#		--minimum-mapping-quality 30 \

############variants filtering##########
############SNP#########################
############Remove previous calibration files######

############ SNP Recalibration ###########
## use 8 processes

       $gatk4 \
       	VariantRecalibrator \
      	-R $ref_genome \
      	-V $temp_dir/raw_variants.vcf \
      	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $ref_path/hapmap_3.3.hg38.vcf \
      	--resource:omni,known=false,training=true,truth=true,prior=12.0 $ref_path/1000G_omni2.5.hg38.vcf \
      	--resource:1000G,known=false,training=true,truth=false,prior=10.0 $ref_path/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf \
      	--resource:dbsnp,known=true,training=false,truth=false,prior=7.0 $dbsnp \
      	--max-gaussians 4 \
      	-an DP -an QD -an FS -an SOR -an MQ -an ReadPosRankSum \
      	-mode SNP \
      	-AS \
      	-tranche 100.0 -tranche 99.5 -tranche 99.0 -tranche 50.0 \
      	-O $temp_dir/raw_recalibrate_SNP.recal \
      	--tranches-file $temp_dir/raw_recalibrate_SNP.tranches \
      	--rscript-file $temp_dir/raw_recalibrate_SNP_plots.R 

######### apply SNP recalib###########

       $gatk4 \
      	ApplyVQSR \
      	-R $ref_genome \
      	-V $temp_dir/raw_variants.vcf \
      	-mode SNP \
      	-AS \
	--truth-sensitivity-filter-level 50.0 \
      	--recal-file $temp_dir/raw_recalibrate_SNP.recal \
      	--tranches-file $temp_dir/raw_recalibrate_SNP.tranches \
      	-O vcf_output/${eachfile/%.bam/}$vcf_suffix

#		> vcf_output/${eachfile/%.bam/}$vcf_suffix

#	rm -r $temp_dir

done

