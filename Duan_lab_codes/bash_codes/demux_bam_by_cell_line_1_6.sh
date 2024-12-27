#!/bin/bash

# Siwei 03 July 2022
# Siwei 05 March 2022

mkdir -p output_demux_bams

shopt -s nullglob

input_folder_list_by_time=(*_lib_GRCh38)
barcodes_list_by_time=( common_barcodes_*/*.tsv )

hg38_autosome_bed="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/GRCh38_chr_1_22.bed"
bam_suffix='.bam'


# echo "${barcodes_list_by_time[@]}"

for (( i=1; i <${#input_folder_list_by_time[@]}; i++ ))
do	
#	echo ${input_folder_list_by_time[@]}
	echo ${input_folder_list_by_time[$i]}
	echo "${barcodes_list_by_time[$i]}"
#	echo ${individual_list[@]}
#	echo "###############"
	cat ${barcodes_list_by_time[i]} | cut -f 2 | sort | uniq | sed 's/\s//g' > indiv.list
	readarray -t individual_list < indiv.list
	for (( j=0; j<${#individual_list[@]}; j++ ))
	do
		echo ${individual_list[@]}
		echo ${individual_list[$j]}

#		output_path='output_demux_bams/'
#		output_path+=${input_folder_list_by_time[$i]/%_input_lib_hybrid/_}
#		output_path+=${individual_list[$j]}
#		output_path+=$bam_suffix
#		echo "$output_path"


		echo "output_demux_bams/${input_folder_list_by_time[$i]/%_input_lib_GRCh38/_}${individual_list[$j]}$bam_suffix"
		echo "*********"
		cat ${barcodes_list_by_time[$i]} | head -n 10
		cat ${barcodes_list_by_time[$i]} \
			| grep "${individual_list[$j]}" \
			| cut -f 1 \
			> barcodes_to_subset.tsv

		rm to_filter.bam
#		rm to_filter_intersected.bam

		~/Data/Tools/subset-bam/subset-bam_linux \
			--bam ${input_folder_list_by_time[$i]}/outs/atac_possorted_bam.bam \
			--cell-barcodes barcodes_to_subset.tsv \
			--out-bam to_filter.bam \
			--cores 10
		
		samtools index -@ 10 to_filter.bam

#		samtools sort \
#			-l 9 -m 10G -@ 10 \
#			to_filter.bam \
#		| samtools view \
#			-L $hg38_autosome_bed \
#			-h -b -u \
#			-@ 10 \
#			-o to_filter_intersected.bam 

#		samtools index -@ 10 to_filter_intersected.bam


		samtools sort \
			-l 9 -m 10G \
			-o output_demux_bams/${input_folder_list_by_time[$i]/%_input_lib_GRCh38/}_${individual_list[$j]}$bam_suffix \
			-@ 10 \
			to_filter.bam

#		rm to_filter_intersected.bam
		rm to_filter.bam

	done
done




#	individual_list=$( echo "${barcodes_list_by_time[$i]}" | cut -f 2 | sort | uniq | sed 's/\s/\n/g' )
