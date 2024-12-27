#! /bin/bash
# Siwei 23 Oct 2023

# make 3x3 bam files around rs16966337 (+/- 1MB) for GViz plotting of ASoC

# find NVME/scARC_Duan_024_raw_fastq/GRCh38_output NVME/scARC_Duan_018/GRCh38_mapped_only Data/FASTQ/Duan_Project_029/GRCh38_only Data/FASTQ/Duan_Project_029_40201 Data/FASTQ/Duan_Project_025_17_46_full_run -type f -regex "^\0*.*WASPed\.bam" -regextype grep | grep "0hr" | awk -F '\/' '{print $NF}' | sort | uniq | wc -l

dir_to_find=('/home/zhangs3/NVME/scARC_Duan_024_raw_fastq/GRCh38_output' '/home/zhangs3/Data/FASTQ/Duan_025_WASPed_BAMs' '/home/zhangs3/NVME/scARC_Duan_018/GRCh38_mapped_only' '/home/zhangs3/Data/FASTQ/Duan_Project_029/GRCh38_only' '/home/zhangs3/Data/FASTQ/Duan_Project_029_40201' '/home/zhangs3/Data/FASTQ/Duan_Project_025_17_46_full_run')
time_stimulation=('0hr' '1hr' '6hr')
temp_folder="/home/zhangs3/NVME/package_temp/subset_rs16966337_temp"

bam_output_dir="/home/zhangs3/NVME/Bulk_ATACseq_ASoC/Gviz_code/plot_rs16966337"

# load het line lists
# readarray -t het_line_list_rs16966337 < "het_rs16966337_cell_lines_2_capture.txt"

cell_type_to_find="nmglut"
for (( i=0; i<${#time_stimulation[@]}; i++ ))
do
	printf "${time_stimulation[$i]}"

	echo "${dir_to_find[@]}"
	[[ -d $temp_folder ]] &&
		rm -r $temp_folder
	mkdir -p $temp_folder

	# load het line lists
	readarray -t het_line_list_rs16966337 < "ASoC_het_samples/"$time_stimulation"_"$cell_type_to_find"_het_sample_names.tsv"

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
			-print0 \) \
			-o \
			\( -type f \
			-regex .*"${time_stimulation[$i]}".*"NEFM_neg".*"WASPed.bam" \
			-regextype grep \
			-print0 \) )

#	echo ${all_bams_to_merge[@]}

	for ((k=0; k<${#all_bams_to_merge[@]}; k++ ))
	do
		echo $(basename -- ${all_bams_to_merge[$k]})

#		[[ -d $temp_folder ]] &&
#			rm -r $temp_folder

		### check if the BAM index file is there
		[[ ! -f ${all_bams_to_merge[$k]}.bai ]] &&
			samtools index -@ 20 \
				${all_bams_to_merge[$k]}

		### check if the BAM is in the het list
		name_to_verify=$(echo $(basename -- ${all_bams_to_merge[$k]}) | sed 's/hr_/\t/g' | sed 's/_NEFM/\t/g' | sed 's/_nmglut/\t/g' | cut -f 2)
		echo "$name_to_verify"
		if [[ $(echo ${het_line_list_rs16966337[@]} | fgrep -w $name_to_verify) ]]
		then
			echo "True"
			### extract +/- 10 MB of rs16966337 (chr16:9788096), total length of chr16=90338345
			samtools view \
				-b \
				-h \
				-@ 30 \
				-o "$temp_folder/$(basename -- ${all_bams_to_merge[$k]})" \
				${all_bams_to_merge[$k]} \
				"chr16:1-20000000"
			samtools index \
				-@ 20 \
				"$temp_folder/$(basename -- ${all_bams_to_merge[$k]})"
		else
			echo "False"
		fi
	done	

	## create and process merged bam
	echo "Merging BAMs..."
	[[ -d "$temp_folder/merging" ]] &&
		rm -r "$temp_folder/merging"
	mkdir -p "$temp_folder/merging"

	samtools merge \
		-@ 40 \
		-o "$temp_folder/merging/merged_unsorted.bam" \
		-s 42 \
		--no-PG \
		-f \
		$temp_folder/*.bam

	samtools sort \
		-l 9 \
		-m 10G \
		-@ 40 \
		-o "$temp_folder/merging/merged_sorted.bam" \
		"$temp_folder/merging/merged_unsorted.bam"
	samtools index \
		-@ 20 \
		"$temp_folder/merging/merged_sorted.bam"

#	samtools view \
#		-b -h -@ 30 \
#		"$temp_folder/merging/merged_sorted.bam" \
#		"chr16:9688096-9888096" \
#		samtools sort \
#		-l 9 -m 5G -@ 40 \
#		-o "$temp_folder/merging/merged_sorted_rs16966337.bam"
#	samtools index \
#		-@ 20 \
#		"$temp_folder/merging/merged_sorted_rs16966337.bam"

	# subsampling by reformat.sh
	reformat.sh \
		ow=t \
		in="$temp_folder/merging/merged_sorted.bam" \
		out="$temp_folder/merging/merged_sorted_10M.bam" \
		sampleseed=42 \
		samplereadstarget=3000000 \
		mappedonly=t \
		pairedonly=t \
		primaryonly=t
	samtools index \
		-@ 20 \
		"$temp_folder/merging/merged_sorted_10M.bam"

	samtools view \
		-b -h -@ 30 \
		"$temp_folder/merging/merged_sorted_10M.bam" \
		"chr16:9688096-9888096" \
		| samtools sort \
		-l 9 -m 5G -@ 40 \
		-o "$bam_output_dir"/"subsampled_""${time_stimulation[$i]}"_"$cell_type_to_find".bam
	samtools index \
		-@ 20 \
		"$bam_output_dir"/"subsampled_""${time_stimulation[$i]}"_"$cell_type_to_find".bam
	
#		"$bam_output_dir"/"subsampled_""${time_stimulation[$i]}"_"$cell_type_to_find".bam
#	echo ${all_bams_to_merge[@]}
done


cell_type_to_find="npglut"
for (( i=0; i<${#time_stimulation[@]}; i++ ))
do
	printf "${time_stimulation[$i]}"

	echo "${dir_to_find[@]}"
	[[ -d $temp_folder ]] &&
		rm -r $temp_folder
	mkdir -p $temp_folder

	# load het line lists
	readarray -t het_line_list_rs16966337 < "ASoC_het_samples/"$time_stimulation"_"$cell_type_to_find"_het_sample_names.tsv"
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
			-print0 \) \
			-o \
			\( -type f \
			-regex .*"${time_stimulation[$i]}".*"NEFM_pos".*"WASPed.bam" \
			-regextype grep \
			-print0 \) )

#	echo ${all_bams_to_merge[@]}

	for ((k=0; k<${#all_bams_to_merge[@]}; k++ ))
	do
		echo $(basename -- ${all_bams_to_merge[$k]})

#		[[ -d $temp_folder ]] &&
#			rm -r $temp_folder

		### check if the BAM index file is there
		[[ ! -f ${all_bams_to_merge[$k]}.bai ]] &&
			samtools index -@ 20 \
				${all_bams_to_merge[$k]}

		### check if the BAM is in the het list
		name_to_verify=$(echo $(basename -- ${all_bams_to_merge[$k]}) | sed 's/hr_/\t/g' | sed 's/_NEFM/\t/g' | sed 's/_npglut/\t/g' | cut -f 2)
		echo "$name_to_verify"
		if [[ $(echo ${het_line_list_rs16966337[@]} | fgrep -w $name_to_verify) ]]
		then
			### extract +/- 10 MB of rs16966337 (chr16:9788096), total length of chr16=90338345
			samtools view \
				-b \
				-h \
				-@ 30 \
				-o "$temp_folder/$(basename -- ${all_bams_to_merge[$k]})" \
				${all_bams_to_merge[$k]} \
				"chr16:1-20000000"
			samtools index \
				-@ 20 \
				"$temp_folder/$(basename -- ${all_bams_to_merge[$k]})"
		fi
		done	

	## create and process merged bam
	[[ -d "$temp_folder/merging" ]] &&
		rm -r "$temp_folder/merging"

	mkdir -p "$temp_folder/merging"
	samtools merge \
		-@ 40 \
		-o "$temp_folder/merging/merged_unsorted.bam" \
		-s 42 \
		--no-PG \
		-f \
		$temp_folder/*.bam

	samtools sort \
		-l 9 \
		-m 10G \
		-@ 40 \
		-o "$temp_folder/merging/merged_sorted.bam" \
		"$temp_folder/merging/merged_unsorted.bam"
	samtools index \
		-@ 20 \
		"$temp_folder/merging/merged_sorted.bam"
	
	reformat.sh \
		ow=t \
		in="$temp_folder/merging/merged_sorted.bam" \
		out="$temp_folder/merging/merged_sorted_10M.bam" \
		sampleseed=42 \
		samplereadstarget=3000000 \
		mappedonly=t \
		pairedonly=t \
		primaryonly=t
	samtools index \
		-@ 20 \
		"$temp_folder/merging/merged_sorted_10M.bam"

	samtools view \
		-b -h -@ 30 \
		"$temp_folder/merging/merged_sorted_10M.bam" \
		"chr16:9688096-9888096" \
		| samtools sort \
		-l 9 -m 5G -@ 40 \
		-o "$bam_output_dir"/"subsampled_""${time_stimulation[$i]}"_"$cell_type_to_find".bam
	samtools index \
		-@ 20 \
		"$bam_output_dir"/"subsampled_""${time_stimulation[$i]}"_"$cell_type_to_find".bam
done

cell_type_to_find="GABA"
for (( i=0; i<${#time_stimulation[@]}; i++ ))
do
	printf "${time_stimulation[$i]}"

	echo "${dir_to_find[@]}"
	[[ -d $temp_folder ]] &&
		rm -r $temp_folder
	mkdir -p $temp_folder
	# find all xhr GABA cells
	# load het line lists
	readarray -t het_line_list_rs16966337 < "ASoC_het_samples/"$time_stimulation"_"$cell_type_to_find"_het_sample_names.tsv"
	## ! the -regex applies to the WHOLE string (i.e., ".*" needs to be placed at the BEGINNING of the whole path
	mapfile \
		-d $'\0' \
		all_bams_to_merge \
		< \
		<(find \
			"${dir_to_find[@]}" \
			-type f \
			-regex ".*"${time_stimulation[$i]}".*"${cell_type_to_find}".*WASPed\.bam$" \
			-regextype grep \
			-print0)

	for ((k=0; k<${#all_bams_to_merge[@]}; k++ ))
	do
		echo $(basename -- ${all_bams_to_merge[$k]})

#		[[ -d $temp_folder ]] &&
#			rm -r $temp_folder

		### check if the BAM index file is there
		[[ ! -f ${all_bams_to_merge[$k]}.bai ]] &&
			samtools index -@ 20 \
				${all_bams_to_merge[$k]}

		### check if the BAM is in the het list
		name_to_verify=$(echo $(basename -- ${all_bams_to_merge[$k]}) | sed 's/hr_/\t/g' | sed 's/_GABA/\t/g' | cut -f 2)
		echo "$name_to_verify"

		if [[ $(echo ${het_line_list_rs16966337[@]} | fgrep -w $name_to_verify) ]]
		then
		### extract +/- 10 MB of rs16966337 (chr16:9788096), total length of chr16=90338345
		samtools view \
			-b \
			-h \
			-@ 30 \
			-o "$temp_folder/$(basename -- ${all_bams_to_merge[$k]})" \
			${all_bams_to_merge[$k]} \
			"chr16:1-20000000"
		samtools index \
			-@ 20 \
			"$temp_folder/$(basename -- ${all_bams_to_merge[$k]})"
		fi
	done	

	## create and process merged bam
	[[ -d "$temp_folder/merging" ]] &&
		rm -r "$temp_folder/merging"
	mkdir -p "$temp_folder/merging"

	samtools merge \
		-@ 40 \
		-o "$temp_folder/merging/merged_unsorted.bam" \
		-s 42 \
		--no-PG \
		-f \
		$temp_folder/*.bam

	samtools sort \
		-l 9 \
		-m 10G \
		-@ 40 \
		-o "$temp_folder/merging/merged_sorted.bam" \
		"$temp_folder/merging/merged_unsorted.bam"
	samtools index -@ 20 \
		"$temp_folder/merging/merged_sorted.bam"

	reformat.sh \
		ow=t \
		in="$temp_folder/merging/merged_sorted.bam" \
		out="$temp_folder/merging/merged_sorted_10M.bam" \
		sampleseed=42 \
		samplereadstarget=3000000 \
		mappedonly=t \
		pairedonly=t \
		primaryonly=t
	samtools index \
		-@ 20 \
		"$temp_folder/merging/merged_sorted_10M.bam"

	samtools view \
		-b -h -@ 30 \
		"$temp_folder/merging/merged_sorted_10M.bam" \
		"chr16:9688096-9888096" \
		| samtools sort \
		-l 9 -m 5G -@ 40 \
		-o "$bam_output_dir"/"subsampled_""${time_stimulation[$i]}"_"$cell_type_to_find".bam
	samtools index \
		-@ 20 \
		"$bam_output_dir"/"subsampled_""${time_stimulation[$i]}"_"$cell_type_to_find".bam

#	echo ${all_bams_to_merge[@]}
done

rm -r $temp_folder

