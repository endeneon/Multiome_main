#!/bin/bash

# Siwei 14 Nov 2022

## !! make softlink of BED4 intervals in ./ or just use $bed4_path variable

gatk4="/home/zhangs3/Data/Tools/gatk-4.2.6.1/gatk"
samtools="/home/zhangs3/Data/Anaconda3-envs/aligners/bin/samtools"

bwa_ref="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"

bed4_path="/home/zhangs3/Data/FASTQ/Duan_Project_024/hybrid_output/per_library_vcf/cleaned_vcf/output/het_vcf/bed4"
SNP_prefix="/home/zhangs3/Data/FASTQ/Duan_Project_024/hybrid_output/per_library_vcf/cleaned_vcf/output/het_vcf/output"

ref_path="/home/zhangs3/Data/Databases/Genomes/hg38"
ref_genome="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"

nonaligned_bam_output_dir="nonaligned_bams"
mkdir -p $nonaligned_bam_output_dir
vcf_output_dir="vcf_output"
mkdir -p $vcf_output_dir

temp_dir="/home/zhangs3/NVME/package_temp/GATK_calib_temp_0hr_pt1"
mkdir -p $temp_dir


eachfile="0hr_CD_52_NEFM_neg_gluthg38onlyintersected.bam"

	samtools index -@ 20 $eachfile

        ############## locate interval file #######
	        interval=$(echo $eachfile | \
			                sed 's/.*\(CD_[0-9][0-9]\).*/\1/g')
	        interval_file="interval_bed"/$interval"_het_only_from_genotyping.bed"
	        echo $interval_file

		samtools view \
			-b -h -@ 16 \
			--write-index \
			-L $interval_file \
			-U $nonaligned_bam_output_dir/non_intersected_reads.bam \
			-o $nonaligned_bam_output_dir/intersected_reads.bam \
			$eachfile
	
	rm -r $temp_dir

#        ## run GATK HaplotypeCaller
#               $gatk4 \
#                        HaplotypeCaller \
#                        --native-pair-hmm-threads 32 \
#                        -R $ref_genome \
#                        -I $eachfile \
#                        -L $interval_file \
#                        -XL $ref_path/hg38_blacklisted_regions.bed \
#                        -O $vcf_output_dir/${eachfile/%.bam/.vcf} \
#                        -bamout $bam_output_dir/${eachfile/%.bam/.remap_paired.bam} \
#                        -DF NotDuplicateReadFilter

