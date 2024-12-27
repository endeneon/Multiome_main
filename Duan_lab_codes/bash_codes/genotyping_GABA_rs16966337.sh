#! /bin/bash
# Siwei 24 Oct 2023


#! /bin/bash
# Siwei 23 Oct 2023

# make 3x3 bam files around rs16966337 (+/- 1MB) for GViz plotting of ASoC

# find NVME/scARC_Duan_024_raw_fastq/GRCh38_output NVME/scARC_Duan_018/GRCh38_mapped_only Data/FASTQ/Duan_Project_029/GRCh38_only Data/FASTQ/Duan_Project_029_40201 Data/FASTQ/Duan_Project_025_17_46_full_run -type f -regex "^\0*.*WASPed\.bam" -regextype grep | grep "0hr" | awk -F '\/' '{print $NF}' | sort | uniq | wc -l

dir_to_find=('/home/zhangs3/NVME/scARC_Duan_024_raw_fastq/GRCh38_output' '/home/zhangs3/Data/FASTQ/Duan_025_WASPed_BAMs' '/home/zhangs3/NVME/scARC_Duan_018/GRCh38_mapped_only' '/home/zhangs3/Data/FASTQ/Duan_Project_029/GRCh38_only' '/home/zhangs3/Data/FASTQ/Duan_Project_029_40201' '/home/zhangs3/Data/FASTQ/Duan_Project_025_17_46_full_run')
time_stimulation=('0hr' '1hr' '6hr')
temp_folder="/home/zhangs3/NVME/package_temp/subset_rs16966337_temp"

bam_output_dir="/home/zhangs3/NVME/Bulk_ATACseq_ASoC/Gviz_code/plot_rs16966337"

[[ -f rs16966337_genotype_list.txt ]] &&
	rm rs16966337_genotype_list.txt


cell_type_to_find="GABA"
for (( i=0; i<${#time_stimulation[@]}; i++ ))
do
	printf "${time_stimulation[$i]}"

	echo "${dir_to_find[@]}"
	[[ -d $temp_folder ]] &&
		rm -r $temp_folder
	mkdir -p $temp_folder
	# find all xhr GABA cells
	## ! the -regex applies to the WHOLE string (i.e., ".*" needs to be placed at the BEGINNING of the whole path
	mapfile \
		-d $'\0' \
		all_bams_to_merge \
		< \
		<(find \
			"${dir_to_find[@]}" \
			\( -type f \
			-regex ".*"${time_stimulation[$i]}".*"${cell_type_to_find}".*WASPed\.bam$" \
			-regextype grep \
			-print0 \) )

#	echo ${all_bams_to_merge[@]}

	for ((k=0; k<${#all_bams_to_merge[@]}; k++ ))
	do
		echo $(basename -- ${all_bams_to_merge[$k]})"__"${cell_type_to_find}"__"${time_stimulation[$i]}

		### check if the BAM index file is there
		[[ ! -f ${all_bams_to_merge[$k]}.bai ]] &&
			samtools index -@ 20 \
				${all_bams_to_merge[$k]}
		
		####
		printf "$(basename -- ${all_bams_to_merge[$k]})\t${cell_type_to_find}\t${time_stimulation[$i]}\t" >> rs16966337_genotype_GABA_list.txt
		### genotype at the site rs16966337 (chr16:9788096)
		samtools mpileup \
			-r "chr16:9788096-9788096" \
		       	${all_bams_to_merge[$k]} \
			>> rs16966337_genotype_GABA_list.txt
	
	done	

done
rm -r $temp_dir
