#!/bin/bash

# Siwei 15 Aug 2023

# Merge demultiplexed bam file at individual (sample) level (merge 0/1/6 hrs);
# Generate an array contains all individual names and iterate over;
# Will need fast add readgroup and dedup methods (does not need to be precise)
# 	based on samtools;
# Keep chr1-22 only;
# call SNPs use GATK 4261

# set error trap
# exit when any command fails
set -e
# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

#####
gatk4="/home/zhangs3/Data/Tools/gatk-4.2.6.1/gatk"
picard_2232="/home/zhangs3/Data/Tools/gatk-4.2.6.1/picard.jar"

ref_path="/home/zhangs3/Data/Databases/Genomes/hg38"
interval_file="/home/zhangs3/Data/Databases/Genomes/hg38/GRCh38_10x_somatic_chr_only.bed"
blacklisted_regions_file="/home/zhangs3/Data/Databases/Genomes/hg38/hg38_blacklisted_regions.bed"
ref_genome="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
dbsnp="/home/zhangs3/Data/Databases/Genomes/hg38/dbsnp_154.hg38.vcf"
#####
output_dir="output_vcf_4_WASP"
mkdir -p $output_dir
##### magic, do not change!!!
indiv_list=( $(ls *.bam \
	| sed 's/.*CD_\([0-9]\{1,\}\)_.*/\1/p' \
	| sort \
	| uniq) )

echo ${indiv_list[@]}

for each_indiv in ${indiv_list[@]}
do
	echo $each_indiv
	
	bam_list=( *"CD_"$each_indiv*".bam" )
	echo ${bam_list[@]}

	temp_dir="/home/zhangs3/NVME/package_temp/temp_scATAC_CD_"$each_indiv
	mkdir -p $temp_dir

	##### merge reads and sort
	samtools merge \
		-l 5 -@ 20 -f \
		--no-PG \
		$temp_dir"/merged_unsorted.bam" \
		${bam_list[@]}

	samtools sort \
		-l 9 -m 5G -@ 20 \
		-o $temp_dir"/merged_sorted.bam" \
		$temp_dir"/merged_unsorted.bam"
	samtools index \
		-@ 20 \
		$temp_dir"/merged_sorted.bam"

	##### Add VN: tag to header (as required by AddOrReplaceReadGroups)
        samtools view -H \
               $temp_dir"/merged_sorted.bam" \
               | sed 's/^@HD/@HD\tVN:1\.5/g' \
               > $temp_dir"/new_header.sam"

        ##### replace the bam header
        samtools reheader \
                $temp_dir"/new_header.sam" \
                $temp_dir"/merged_sorted.bam" \
                | samtools sort \
                -@ 16 -m 5G -l 9 \
                -o $temp_dir"/merged_sorted_header.bam"
        samtools index \
		-@ 20 \
		$temp_dir"/merged_sorted_header.bam"

        ##### dedup
        java -Xmx100g -jar $picard_2232 \
                MarkDuplicates \
                -I $temp_dir"/merged_sorted_header.bam" \
                -O $temp_dir"/merged_sorted_dedup_header.bam" \
                -M $temp_dir"/metrics.txt" \
                --REMOVE_DUPLICATES true \
		--REMOVE_SEQUENCING_DUPLICATES true \
                --ASSUME_SORTED true
	samtools index \
		-@ 20 \
		$temp_dir"/merged_sorted_dedup_header.bam"

	##### Strip all duplicate flags (0x400)	
	samtools view \
		-h -S -@ 10 \
               --remove-flags 0x400 \
               $temp_dir"/merged_sorted_dedup_header.bam" \
               | samtools sort \
               -@ 16 -m 5G -l 9 \
               -o $temp_dir"/merged_sorted_dedup_header_nodup.bam"
	samtools index \
		-@ 20 \
		$temp_dir"/merged_sorted_dedup_header_nodup.bam"

	##### Re-add readgroup
        samtools addreplacerg \
                -w \
                -r "ID:CD-$each_indiv\tLB:lib1\tPL:illumina\tSM:CD-$each_indiv\tPU:unit1" \
                -@ 16 \
                $temp_dir"/merged_sorted_dedup_header_nodup.bam" \
                | samtools sort \
                -l 9 -m 10G -@ 16 \
		-o $temp_dir"/GATK_call_source.bam"
	samtools index \
		-@ 20 \
		$temp_dir"/GATK_call_source.bam"

	##### Call SNP by GATK
        $gatk4 --java-options "-Xmx150g" \
                HaplotypeCaller \
                --native-pair-hmm-threads 16 \
                -R $ref_genome \
                --dbsnp $dbsnp \
                -I $temp_dir"/GATK_call_source.bam" \
                -L $interval_file \
                -XL $blacklisted_regions_file \
                -O $temp_dir"/raw_variants.vcf"

	##### select SNP for output from autosomes only
        $gatk4 --java-options "-Xmx150g" \
                SelectVariants \
                -R $ref_genome \
                -V $temp_dir"/raw_variants.vcf" \
                --select-type-to-include SNP \
                --restrict-alleles-to BIALLELIC \
                -L $ref_path"/autosomes.list" \
                --exclude-filtered true \
                -O $temp_dir"/pre_output.vcf"

	##### extract het only
        cat $temp_dir"/pre_output.vcf" \
                | grep "^#" \
                > $temp_dir"/header.txt"

        cat $temp_dir"/pre_output.vcf" \
                | grep -v "^#" \
                | grep "0/1" \
                > $temp_dir"/body.txt"

        cat $temp_dir"/header.txt" \
                $temp_dir"/body.txt" \
                > $output_dir"/vcf_4_WASP_CD-"$each_indiv".vcf"

	##### cleanup
        rm -r $temp_dir

done

