#!/bin/bash

# Siwei 13 Nov 2022
# Extract reads intersected with a designated bed file 
# (will also lift its mate pair, so will output paired ends rather than 
# sometimes fragmented reads only)

# Designed for GATK-WASP calibration workflow on scARC files

samviewwithmate="/home/zhangs3/Data/Tools/jvarkit/dist/samviewwithmate.jar"

temp_dir="/home/zhangs3/NVME/package_temp/samviewwithmate_0hr_pt1"

mkdir -p bam_intersected
# mkdir -p $temp_dir

for eachfile in *.bam
do
	mkdir -p $temp_dir

	samtools index -@ 20 $eachfile

	############## locate interval file #######
        interval=$(echo $eachfile | \
                sed 's/.*\(CD_[0-9][0-9]\).*/\1/g')
        interval_file="interval_bed"/$interval"_het_only_from_genotyping.bed"
        echo $interval_file

       samtools sort \
	       -l 5 -m 5G -@ 16 \
	       -o $temp_dir/merged_sorted_dedup.bam \
               $eachfile
       samtools index -@ 20 $temp_dir/merged_sorted_dedup.bam


        ############### Add SO:tag to header (as required by AddOrReplaceReadGroups
        # need to extract sam header first
        samtools view -H \
              $temp_dir/merged_sorted_dedup.bam \
              | sed 's/^@HD/@HD\tVN:1\.5/g' \
              > $temp_dir/new_header.sam

        # replace the bam header
        samtools reheader \
                $temp_dir/new_header.sam \
                $temp_dir/merged_sorted_dedup.bam \
                | samtools sort \
                -@ 16 -m 5G -l 5 \
                -o $temp_dir/merged_sorted_dedup_header.bam

        samtools index -@ 20 $temp_dir/merged_sorted_dedup_header.bam

	# extract the reads intersected with the bed file
	# note that the output bam does not have a header so needs to be added back
	# from the extracted "new_header.sam"
	echo "extracting reads..."
	java -jar \
		$samviewwithmate \
		--bamcompression 9 \
		-o $temp_dir/extracted_intersected.bam \
		-b $interval_file \
		--samoutputformat BAM \
		$temp_dir/merged_sorted_dedup_header.bam

        # Add the bam header
        samtools reheader \
                $temp_dir/new_header.sam \
                $temp_dir/extracted_intersected.bam \
                | samtools sort \
                -@ 16 -m 5G -l 5 \
                -o bam_intersected/${eachfile/%.bam/intersected.bam}

#	${eachfile/%.bam/intersected.bam}
	
	rm -r $temp_dir
done


