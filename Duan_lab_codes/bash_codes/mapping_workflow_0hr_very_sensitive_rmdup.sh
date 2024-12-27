#!/bin/bash

# Siwei 10 May 2021
# customise code for genres01

# Siwei 20 Dec 2019
# updated WASP to 0.3.4 (30 Apr 2019)


# Set these environment vars to point to
# your local installation of WASP

export OMP_NUM_THREADS=1

WASP="/home/zhangs3/Data/Tools/WASP"
jre_8="/home/zhangs3/Data/Tools/jre1.8.0_291/bin/java"
picard_291="/home/zhangs3/Data/Tools/picard291.jar"


DATA_DIR="/home/zhangs3/NVME/package_temp/WASP_temp_0hr"
OUTPUT_FINAL="WASPed_BAMs"
# SNP_DIR="CD_50_WASP_SNP"

# These environment vars point to the reference genome and bowtie2.
# in the examples below, the reference genome is assumed
# to be indexed for use with bowtie2
INDEX="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"

mkdir -p WASPed_BAMs

# Start the loop
for EACHFILE in *.bam
do

	samtools index -@ 20 $EACHFILE

# set WASP calibration SNP path 	
#        SNP_DIR=$(echo $EACHFILE | \
#                sed 's/.\+\(CD_[0-9][0-9]\).*/\1/g')
#        SNP_DIR=../WASP_calib_source/$SNP_DIR"_WASP_SNP"
	SNP_DIR="/home/zhangs3/Data/FASTQ/Duan_Project_024/hybrid_output/bam_dump_4_fasta_01Jun2022/WASP_to_calibrate_by_line/WASP_calib_source/CD_51_WASP_SNP"
        echo $SNP_DIR

# pre-cleanup
	rm -r $DATA_DIR/

# initialise
	mkdir $DATA_DIR/
	echo $EACHFILE
	echo $DATA_DIR/${EACHFILE/%.bam/.fq.gz}
	date


	ls $DATA_DIR/

# Pull out reads that need to be remapped to check for bias
# Use the -p option for paired-end reads.

	python $WASP/mapping/find_intersecting_snps.py \
		--is_paired_end \
		--is_sorted \
		--output_dir $DATA_DIR/ \
		--snp_dir $SNP_DIR \
		$EACHFILE

# Remap the reads, using same the program and options as before.
# NOTE: If you use an option in the first mapping step that modifies the
# reads (e.g. the -5 read trimming option to bowtie2) you should omit this
# option during the second mapping step here (otherwise the reads will be modified
# twice)!         
	echo "map for paired ends"      
        bowtie2 \
		--very-sensitive-local \
		-p 40 \
		-t \
		-X 10000 \
                -x $INDEX \
		-1 $DATA_DIR/${EACHFILE/%.bam/.remap.fq1.gz} \
		-2 $DATA_DIR/${EACHFILE/%.bam/.remap.fq2.gz} \
                | samtools sort \
		-l 9 -@ 20 -m 10G \
		-o $DATA_DIR/${EACHFILE/%.bam/.remap_paired.bam}

# remap to same position
        python $WASP/mapping/filter_remapped_reads.py \
                $DATA_DIR/${EACHFILE/%.bam/.to.remap.bam} \
                $DATA_DIR/${EACHFILE/%.bam/.remap_paired.bam} \
                $DATA_DIR/${EACHFILE/%.bam/.remap_paired.keep.bam}
                ls -lh $DATA_DIR/

################### map and filter for unpaired ones
        
#        echo "map for unpaired ends"
#        bowtie2 -p 16 -X 2000 \
#		--mm --qc-filter \
#		--met 1 -t \
#                --sensitive --no-mixed --no-discordant \
#                -x $INDEX \
#		-U $DATA_DIR/${EACHFILE/%.bam/.remap.single.fq.gz} \
#                | samtools sort -l 9 -@ 20 -m 10G \
#		-o $DATA_DIR/${EACHFILE/%.bam/.remap_unpaired.bam}

# remap to same position
#	echo "remap"
#        python $WASP/mapping/filter_remapped_reads.py \
#                $DATA_DIR/${EACHFILE/%.bam/.to.remap.bam} \
#                $DATA_DIR/${EACHFILE/%.bam/.remap_unpaired.bam} \
#                $DATA_DIR/${EACHFILE/%.bam/.remap_unpaired.keep.bam}
#        ls -lh $DATA_DIR/

################################## merge, sort, and dedup

# Create a merged BAM containing [1] reads (PE+SE) that did
# not need remapping [2] filtered remapped reads
	echo "merge bam"
       samtools merge -l 9 -@ 10 -f \
                $DATA_DIR/sim_pe_reads.keep.merged.bam \
                $DATA_DIR/${EACHFILE/%.bam/.remap_paired.keep.bam} \
                $DATA_DIR/${EACHFILE/%.bam/.keep.bam}

#                $DATA_DIR/${EACHFILE/%.bam/.remap_unpaired.keep.bam} \
################### Sort and index the bam file
	echo "sort bam"
        samtools sort -l 9 -m 10G -@ 20 \
                -o $DATA_DIR/sim_pe_reads.keep.merged.sorted.bam \
                $DATA_DIR/sim_pe_reads.keep.merged.bam 

        samtools index -@ 20 $DATA_DIR/sim_pe_reads.keep.merged.sorted.bam
	ls -lh $DATA_DIR/

################### Dedup
#        $jre_8 -Xmx100g -jar $picard_291 MarkDuplicates \
#		I=$DATA_DIR/sim_pe_reads.keep.merged.sorted.bam \
#		O=$DATA_DIR/merged_sorted_dedup.bam \
#		M=metrics.txt \
#               REMOVE_DUPLICATES=True

################# Dedup by WASP rmdup
	python $WASP/mapping/rmdup_pe.py \
		$DATA_DIR/sim_pe_reads.keep.merged.sorted.bam \
		$DATA_DIR/merged_dedup.bam

	ls -lh $DATA_DIR/

        samtools sort -l 9 -m 10G -@ 20 \
                -o $DATA_DIR/merged_sorted_dedup.bam \
                $DATA_DIR/merged_dedup.bam


################## Re-add readgroup     
        $jre_8 -Xmx100g -jar $picard_291 AddOrReplaceReadGroups \
                        INPUT=$DATA_DIR/merged_sorted_dedup.bam \
                        OUTPUT=$OUTPUT_FINAL/${EACHFILE/%.bam/_WASPed.bam} \
                        RGID=sample_51_0 \
                        RGLB=lib1 \
                        RGPL=illumina \
                        RGPU=unit1 \
                        RGSM=sample_51_0


	ls -lh $DATA_DIR/
        ls -lh $OUTPUT_FINAL/

################################## clean up the temp folder

        rm -r $DATA_DIR/


done

unset OMP_NUM_THREADS
