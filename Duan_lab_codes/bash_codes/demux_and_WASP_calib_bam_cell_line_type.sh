#!/bin/bash

# Siwei 23 August 2023
# rewrite to
# 1) Directly subset from master atac_possorted files and extract bam
# 2) WASP-calibrate bam on-the-fly
# 3) write WASPed bams
# 4) for each library (4 lines), run 3 time points sequentially using for loop

# Siwei 08 March 2023
# Siwei 05 March 2022

# Set these environment vars to point to
# your local installation of WASP
source /home/zhangs3/Data/Tools/cellranger-arc-2.0.2/sourceme.bash
export OMP_NUM_THREADS=8

# assign constants
WASP="/home/zhangs3/Data/Tools/WASP"
picard_291="/home/zhangs3/Data/Tools/picard291.jar"
picard_2254="/home/zhangs3/Data/Tools/gatk-4.1.8.1/picard2_25_4.jar"
gatk4="/home/zhangs3/Data/Tools/gatk-4.1.8.1/gatk"
samtools="/home/zhangs3/Data/Anaconda3-envs/aligners/bin/samtools"
parallel="/home/zhangs3/Data/Anaconda3-envs/aligners/bin/parallel"

# in the examples below, the reference genome is assumed
# to be indexed for use with BWA (comes from Cell Ranger)
# use BWA reference from Cell Ranger_ref
bwa_ref="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
barcode_time_prefix=('0hr' '1hr' '6hr')
hg38_autosome_bed="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/GRCh38_chr_1_22.bed"
bam_suffix='_bwa_barcode_WASPed.bam'

# set error trap
# exit when any command fails
set -e
# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting

input_folder_list_by_time=( $( ls -d Duan_024_* ) )


mkdir -p output_demux_bams
shopt -s nullglob

for (( i=0; i <${#input_folder_list_by_time[@]}; i++ ))
do	
	echo ${input_folder_list_by_time[$i]}
	echo "###############"
	individual_list=( $( ls all_barcodes/${barcode_time_prefix[$i]}*barcodes.txt ) )

	for (( j=0; j<${#individual_list[@]}; j++ ))
	do
		temp_dir=$( basename -- ${individual_list[$j]} )
		temp_dir=${temp_dir/%_barcodes.txt/}
		temp_dir="/home/zhangs3/NVME/package_temp"/$temp_dir
                echo "temp_dir = "$temp_dir
		mkdir -p $temp_dir

		echo ${individual_list[$j]}

		output_bam_name=$( basename -- ${individual_list[$j]} )
		snp_dir=$output_bam_name
		output_bam_name="output_demux_bams"/${output_bam_name/%_barcodes.txt/}$bam_suffix
                echo "output_bam_name = "$output_bam_name

		echo "#########"
		echo "Subsetting atac_possorted_bam.bam by barcodes, this may take a long time..."
		cat ${individual_list[$j]} \
			| cut -d "_" -f 1 \
			| sort \
			| uniq \
			> $temp_dir/"trimmed_barcode_list.txt"

		~/Data/Tools/subset-bam/subset-bam_linux \
			--bam ${input_folder_list_by_time[$i]}/"outs/atac_possorted_bam.bam" \
			--cell-barcodes $temp_dir/"trimmed_barcode_list.txt" \
			--out-bam $temp_dir/"to_filter.bam" \
			--cores 20
		
		$samtools index \
			-@ 20 \
			$temp_dir/"to_filter.bam"

		$samtools sort \
			-l 9 -m 5G -@ 20 \
			$temp_dir/"to_filter.bam" \
		| $samtools view \
			-L $hg38_autosome_bed \
			-h -b -u \
			-@ 20 \
		| $samtools sort \
			-l 9 -m 8G \
			-o $temp_dir/"filtered_autosome_chr_only.bam" \
			-@ 20

		$samtools index \
			-@ 20 \
			$temp_dir/"filtered_autosome_chr_only.bam"

		rm $temp_dir/"to_filter.bam"

		##### Calibrate and dedup the output bam ($temp_dir/"filtered_autosome_chr_only.bam") by WASP
		
		# set WASP calibration SNP path
                snp_dir=$( echo $snp_dir \
                        | sed 's/CD_/CD-/g' \
                        | cut -d "_" -f 2 )
                snp_dir="/home/zhangs3/NVME/scARC_Duan_024_raw_fastq/GRCh38_output/WASP_calib_all_lines"/$snp_dir
                echo "snp_dir = "$snp_dir

		final_output_file_name_prefix=${snp_dir/%.bam/} #! need to deal with this line later

	        ################# Split the source file by barcodes, dedup at barcode-level bam files,
	        ################# and pool all deduped bam files in a temp folder
	        ################# Use 10x subset-bam to subset by barcodes

                ### make a WASP temp folder to hold all barcode-level bam files
		WASP_dir=$temp_dir/"WASP_temp"
		mkdir -p $WASP_dir

                mkdir -p $WASP_dir/"barcode_bams"
                mkdir -p $WASP_dir/"WASP_run"
                mkdir -p $WASP_dir/"deduped_barcode_bams"

                ### make a function for parallel to call WASP rmdup
                ### BAM calibration will also be performed at this stage
		# !! do not change this function! let it run !!
                call_WASP_rmdup () {

                        ### !! constants need to be re-assigned within functions !!

                        bwa_ref="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"
                        samtools="/home/zhangs3/Data/Anaconda3-envs/aligners/bin/samtools"

                        DATA_DIR=$2
                        WASP=$3
                        SNP_DIR=$4
                        eachfile_input=$5

                        echo $1 > $DATA_DIR'/barcode_bams/'$1'_barcode_2_subset.txt'

                        DATA_barcodes_dir=$DATA_DIR'/barcode_bams'
                        WASP_run_dir=$DATA_DIR/'WASP_run'
                        deduped_bam_dir=$DATA_DIR/'deduped_barcode_bams'

                        ### subset by 10x subset bam

                        ### !!! remove the previous subset-bam output bam file first !!!
                        ### !!! otherwise the subset-bam will refuse to overwrite !!!
                        ### use [[]] to test the presence of the target file first so the system does not
                        ### generate an stderr message if no bam (which should be the case) found

                        [[ -f $DATA_barcodes_dir/$1'_barcode_to_dedup_unsorted.bam' ]] && \
                                rm $DATA_barcodes_dir/$1'_barcode_to_dedup_unsorted.bam'

                        ~/Data/Tools/subset-bam/subset-bam_linux \
                                --bam $eachfile_input \
                                --cell-barcodes $DATA_barcodes_dir/$1'_barcode_2_subset.txt' \
                                --out-bam $DATA_barcodes_dir/$1'_barcode_to_dedup_unsorted.bam' \
                                --cores 1

                        $samtools sort \
                                -l 9 -@ 1 \
                                -o $WASP_run_dir/$1'_raw.bam' \
                                $DATA_barcodes_dir/$1'_barcode_to_dedup_unsorted.bam'

                        $samtools index $WASP_run_dir/$1'_raw.bam'

                        eachfile=$1'_raw.bam'

		        # Pull out reads that need to be remapped to check for bias
		        # Use the -p option for paired-end reads.
		        #
	                python $WASP/mapping/find_intersecting_snps.py \
	                        --is_paired_end \
	                        --is_sorted \
	                        --output_dir $WASP_run_dir/ \
	                        --snp_dir $SNP_DIR \
	                        $WASP_run_dir/$eachfile
		
		        # Remap the reads, using same the program and options as before.
		        # NOTE: If you use an option in the first mapping step that modifies the
		        # reads (e.g. the -5 read trimming option to bowtie2) you should omit this
		        # option during the second mapping step here (otherwise the reads will be modified
		        # twice)!
	                echo "map for paired ends"
	                bwa mem \
	                        -I 250.0,150.0 \
	                        -t 1 \
	                        $bwa_ref \
	                        $WASP_run_dir/${eachfile/%.bam/.remap.fq1.gz} \
	                        $WASP_run_dir/${eachfile/%.bam/.remap.fq2.gz} \
	                        | $samtools sort \
	                        -o $WASP_run_dir/${eachfile/%.bam/.remap_paired.bam}
	
	                $samtools index $WASP_run_dir/${eachfile/%.bam/.remap_paired.bam}

			## remap to same position
	                python $WASP/mapping/filter_remapped_reads.py \
	                        $WASP_run_dir/${eachfile/%.bam/.to.remap.bam} \
	                        $WASP_run_dir/${eachfile/%.bam/.remap_paired.bam} \
	                        $WASP_run_dir/${eachfile/%.bam/.remap_paired.keep.bam}
	
		        #
		        ################################### merge, sort, and dedup
		        #
		        ## Create a merged BAM containing [1] reads (PE+SE) that did
		        ## not need remapping [2] filtered remapped reads
		        #       echo "merge bam"
		
	               $samtools merge \
	                       -f \
	                       --no-PG \
	                        $WASP_run_dir/${eachfile/%.bam/.sim_pe_reads.merged.bam} \
	                        $WASP_run_dir/${eachfile/%.bam/.remap_paired.keep.bam} \
	                        $WASP_run_dir/${eachfile/%.bam/.keep.bam}
		
		
		        #################### Sort and index the bam file
		        ####### !! Retain Uniquely mapped reads only (-v XA:Z:, SA:Z:)before send
		                echo "filter and sort bam"
		                $samtools view \
		                        -h \
		                        $WASP_run_dir/${eachfile/%.bam/.sim_pe_reads.merged.bam} \
		                       | $samtools sort \
		                       -o $WASP_run_dir/${eachfile/%.bam/.sim_pe_reads.merged.sorted.bam}
	                $samtools index $WASP_run_dir/${eachfile/%.bam/.sim_pe_reads.merged.sorted.bam}

                        ### Dedup by WASP rmdup

                        python $WASP/mapping/rmdup_pe.py \
                                $WASP_run_dir/${eachfile/%.bam/.sim_pe_reads.merged.sorted.bam} \
                                $deduped_bam_dir/$1'.bam'

                        ### !! remove the subsetted files in $DATA_DIR !!
                        rm $DATA_barcodes_dir/$1*
                        rm $WASP_run_dir/$1*

                }

                # export the function
                export -f call_WASP_rmdup
                # ! need to take care of the line below as well !
		# input_barcode_file="../../../barcodes_by_time_4_WASP_rmdup/barcode_suffix_stripped/${eachfile_input/%_hg38only.bam/_barcodes.txt}"
		input_barcode_file=$temp_dir/"trimmed_barcode_list.txt"
                ### read in barcodes using parallel -a
                ### there is only one barcode per line
		echo "Running WASP calibration..."
                $parallel \
			-a $input_barcode_file \
			-j 24 \
                        call_WASP_rmdup {} $WASP_dir $WASP $snp_dir $temp_dir/"filtered_autosome_chr_only.bam"

		        #################### merge all barcode-level bam files into one

                $samtools merge \
                        -l 9 -@ 16 -f \
                        --no-PG \
                        $WASP_dir/merged_dedup.bam \
                        $WASP_dir/deduped_barcode_bams/*.bam

	        ################### remove barcodes bam files
	
                rm -r $WASP_dir/deduped_barcode_bams/

		echo "Merging and post-processing WASPed barcode-wise BAMs..."

                $samtools sort -l 9 -m 10G -@ 16 \
                        -o $WASP_dir/merged_sorted_dedup.bam \
                        $WASP_dir/merged_dedup.bam
                $samtools index -@ 20 $WASP_dir/merged_sorted_dedup.bam
	        #
	        ############### Add SO:tag to header (as required by AddOrReplaceReadGroups
                # need to extract sam header first
                $samtools view -H \
                        $WASP_dir/merged_sorted_dedup.bam \
                        | sed 's/^@HD/@HD\tVN:1\.5/g' \
                        > $WASP_dir/new_header.sam

                # replace the bam header
                $samtools reheader \
                        $WASP_dir/new_header.sam \
                        $WASP_dir/merged_sorted_dedup.bam \
                        | $samtools sort \
                        -@ 16 -m 5G -l 9 \
			-o $WASP_dir/merged_sorted_dedup_header.bam

		$samtools index -@ 20 $WASP_dir/merged_sorted_dedup_header.bam

	        ################### Strip all duplicate flags (0x400)
	
                $samtools view -h -S -@ 10 \
                        --remove-flags 0x400 \
			$WASP_dir/merged_sorted_dedup_header.bam \
                        | $samtools sort \
                        -@ 16 -m 5G -l 9 \
                        -o $WASP_dir/merged_sorted_dedup_header_nodup.bam
	
	        ################### Re-add readgroup
		rg_string=$( basename -- $final_output_file_name_prefix )
		rg_string=$( echo $rg_string \
			| cut -d "." -f 1 )
		echo "rg_string = "$rg_string

	        $samtools addreplacerg \
	                -w \
	                -r "ID:$rg_string\tLB:lib1\tPL:illumina\tSM:$rg_string\tPU:unit1" \
	                -@ 16 \
	                $WASP_dir/merged_sorted_dedup_header_nodup.bam \
	                | $samtools sort \
	                -l 9 -m 10G -@ 16 \
	                -o $output_bam_name
		
		###### clean up before exit loop
		rm -r $temp_dir
		echo '############################'
	done
done

unset OMP_NUM_THREADS
