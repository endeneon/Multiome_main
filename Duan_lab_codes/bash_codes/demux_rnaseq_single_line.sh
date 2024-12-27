#!/bin/bash

# Siwei 05 March 2022

mkdir -p output_demux_bams

shopt -s nullglob

input_folder_list_by_time=( *0_input_lib_GRCh38 )
barcode_time_prefix=('0hr')
# barcodes_list_by_time=( $( ls common_barcodes_*/*.tsv ) )

hg38_autosome_bed="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/GRCh38_chr_1_22.bed"
bam_suffix='hg38only.bam'

mkdir -p "gex_output_Jeff"

# echo "${barcodes_list_by_time[@]}"

for (( i=0; i <${#input_folder_list_by_time[@]}; i++ ))
do	
#	echo ${input_folder_list_by_time[@]}
	echo ${input_folder_list_by_time[$i]}
#	echo "${barcodes_list_by_time[$i]}"
#	echo ${individual_list[@]}
#	echo "###############"
#	cat ${barcodes_list_by_time[i]} | cut -f 2 | sort | uniq | sed 's/\s//g' > indiv.list
#	readarray -t individual_list < indiv.list
	individual_list=( $( ls all_barcodes/${barcode_time_prefix[$i]}*pos*barcodes.txt ) )

	for (( j=0; j<${#individual_list[@]}; j++ ))
	do
		echo ${individual_list[@]}
		echo ${individual_list[$j]}


		output_bam_name=$( basename -- ${individual_list[$j]} )
		output_bam_name="gex_output_Jeff"/${output_bam_name/%_barcodes.txt/}$bam_suffix
		echo $output_bam_name

		echo "*********"
		cat ${individual_list[$j]} | head -n 10

		rm to_filter.bam
		rm to_filter_intersected.bam

		~/Data/Tools/subset-bam/subset-bam_linux \
			--bam ${input_folder_list_by_time[$i]}/outs/gex_possorted_bam.bam \
			--cell-barcodes ${individual_list[$j]} \
			--out-bam to_filter.bam \
			--cores 20
		
		samtools index -@ 20 to_filter.bam

		samtools sort \
			-l 9 -m 5G -@ 20 \
			to_filter.bam \
		| samtools view \
			-L $hg38_autosome_bed \
			-h -b -u \
			-@ 20 \
		| samtools sort \
			-l 9 -m 8G \
			-o $output_bam_name \
			-@ 20

		rm to_filter.bam

	done
done




#	individual_list=$( echo "${barcodes_list_by_time[$i]}" | cut -f 2 | sort | uniq | sed 's/\s/\n/g' )
