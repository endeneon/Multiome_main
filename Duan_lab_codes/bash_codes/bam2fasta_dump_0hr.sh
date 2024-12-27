#!/bin/bash

# Siwei 01 Jun 2022

# Split bam files by time, line, and cell type
# process 0hr, 1hr, 6hr files separately



out_path="bam_dump_4_fasta_04Jul2022_hg38/hr_0"
temp_path="/home/zhangs3/NVME/package_temp/samtools_temp/hr_0"

shopt -s nullglob

hg38_autosome_bed="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/GRCh38_chr_1_22.bed"
bam_suffix='.bam'

# load bam list of different times
readarray -t bam_input_list < all_groups/bam_list_0hr.list

# make the array of group names, need to have the same order
# as in the bam input list
group_name_prefix="/home/zhangs3/Data/FASTQ/Duan_Project_024/hybrid_output/txt_barcodes_for_ASoC/hr_0/"
group_names=('group_51')

# for each bam file in the bam_input_list,
# use the corresponding element in the group_names

mkdir -p $out_path
mkdir -p $temp_path

for ((i=0; i<${#bam_input_list[@]}; i++ ))
do
	echo "${bam_input_list[$i]}"
	echo "${group_names[$i]}"

	current_working_group=$group_name_prefix${group_names[$i]}
#	ls $current_working_group

	for eachfile in $current_working_group/*.txt
	do
#		echo $eachfile
                eachfile_basename="$(basename -- $eachfile)"
                echo $eachfile_basename

                ~/Data/Tools/subset-bam/subset-bam_linux \
	                --bam ${bam_input_list[$i]}\
                        --cell-barcodes $eachfile \
                        --out-bam $temp_path/to_filter.bam \
                        --cores 10

                samtools sort \
                        -l 9 -m 8G -@ 10 \
                        $temp_path/to_filter.bam \
                | samtools view \
                        -L $hg38_autosome_bed \
                        -h -b -u \
                        -@ 10 \
		| samtools sort \
			-l 9 -m 8G -@ 10 \
			-o $out_path/${eachfile_basename/%_barcodes.txt/.bam}

		rm $temp_path/*
	done

done

rm -r $temp_path


