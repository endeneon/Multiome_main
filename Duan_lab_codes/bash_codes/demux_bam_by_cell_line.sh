#!/bin/bash

# Siwei 05 March 2022

#  mkdir -p output_demux_bams

# shopt -s nullglob

cat atac*.list | sort > input_bam_list.input
readarray -t input_bam_list < input_bam_list.input
# barcodes_list_by_type=('GABA' 'NEFM_neg_glut' 'NEFM_pos_glut')
barcodes_list_by_type=('NEFM_neg_glut')

hg38_autosome_bed="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/GRCh38_chr_1_22.bed"
bam_suffix='.bam'
timestamp='0hr'


# echo "${barcodes_list_by_type[@]}"

mkdir -p output_demux_bams
mkdir -p output_temp_bams



for (( i=0; i < ${#barcodes_list_by_type[@]}; i++ ))
do	
#	echo ${input_bam_list[@]}
#	echo ${input_bam_list[$i]}
#	echo "${#barcodes_list_by_type[@]}"
#	echo "${barcodes_list_by_type[2]}"

	rm output_temp_bams/*

#	echo "${barcodes_list_by_type[@]}"
	echo "${barcodes_list_by_type[$i]}"
#	echo ${individual_list[@]}
#	echo "###############"
#	cat ${barcodes_list_by_type[i]} | cut -f 2 | sort | uniq | sed 's/\s//g' > indiv.list
#	ls barcodes_by_library/${barcodes_list_by_type[$i]}*0_barcodes.txt | sort > indiv.list
#	readarray -t individual_list < indiv.list
	individual_list=( $( ls "barcodes_by_library/${barcodes_list_by_type[$i]}"*0_barcodes.txt | sort ) )

	echo ${individual_list[@]}

	for (( j=0; j<${#input_bam_list[@]}; j++ ))
	do
#		echo ${input_bam_list[@]}
		echo ${input_bam_list[$j]}

		bam_base=$(basename ${input_bam_list[$j]/%.bam/_temp})
		echo "output_temp_bams/$bam_base$j$bam_suffix"
		echo "*********"
		echo ${individual_list[$j]}

		
		echo "*********"
		final_output_var=$( echo ${individual_list[$j]} \
			| sed 's/.*\(\(GABA\|NEFM\).*\-[016]\).*/\1/' )
		echo $final_output_var

		cell_type=$( echo $final_output_var \
			| sed 's/.*\(\(GABA\|NEFM_neg_glut\|NEFM_pos_glut\)\).*/\1/' )
		cell_line=$( echo $final_output_var \
			| sed 's/.*\([0-9][0-9]\-[0-9]\).*/\1/' \
			| sed 's/\-/\t/g' \
			| awk '{print "CD_"$1}' )

		final_output_bam=$timestamp"_"$cell_line"_hg38"$cell_type$bam_suffix

		echo "Output file name is $final_output_bam"

		~/Data/Tools/subset-bam/subset-bam_linux \
			--bam ${input_bam_list[$j]} \
			--cell-barcodes "${individual_list[$j]}" \
			--out-bam output_temp_bams/to_filter.bam \
			--cores 20
		

	

        samtools sort \
	        -l 0 -m 10G -@ 10 \
        	output_temp_bams/to_filter.bam \
                | samtools view \
                -L $hg38_autosome_bed \
                -h -b -1 \
                -@ 10 \
                | samtools sort \
                -l 9 -m 10G -@ 10 \
                -o output_demux_bams/$final_output_bam

	rm output_temp_bams/*

	done

	echo "##########"
#	final_output_bam=$barcodes_list_by_type$timestamp$bam_suffix

#	samtools merge \
#		-l 9 -@ 10 -f \
#		output_demux_bams/merged_unsorted.bam \
#		output_temp_bams/*.bam



done




#	individual_list=$( echo "${barcodes_list_by_type[$i]}" | cut -f 2 | sort | uniq | sed 's/\s/\n/g' )
