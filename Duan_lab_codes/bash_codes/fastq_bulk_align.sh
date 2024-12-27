#!/bin/bash

# align fastq files dumped from demultiplexed 10x Cellranger ARC (ATAC)
# files as bulk
# Use bowtie2 as the aligner instead of bwa

WASP="/home/zhangs3/Data/Tools/WASP"
OUTPUT_FINAL="WASPed_BAMs"

picard="/home/zhangs3/Data/Tools/gatk-4.2.6.1/picard.jar"
ref_fasta="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"

# output_path="pre_bams"
# mkdir -p $output_path

temp_dir="/home/zhangs3/NVME/package_temp/dumped_fastq_align_hr0_NEFM_neg_pt1"
mkdir -p $OUTPUT_FINAL
# mkdir -p $temp_dir

# get the array of all sample-level directories
sample_directories=( $( ls -d *glut/) )

# rm -r $temp_dir

for ((i=0; i<${#sample_directories[@]}; i++))
do

        [[ -d $temp_dir ]] && rm -r $temp_dir/


	mkdir -p $temp_dir
	DATA_DIR=$temp_dir/"WASP_run_dir"
	mkdir -p $temp_dir/"to_merge"
	mkdir -p $DATA_DIR

	current_sample_name=${sample_directories[$i]}
	current_sample_name=${current_sample_name/%\//}
	echo "current sample name is: "$current_sample_name

        SNP_DIR=$(echo $current_sample_name | \
                sed 's/.\+\(CD_[0-9][0-9]\).*/\1/g')
        SNP_DIR="WASP_calib_known_genotypes"/$SNP_DIR"_WASP_SNP_fr_genotyping"

	echo $SNP_DIR
	echo ${sample_directories[$i]}



      # get the array of R1 reads below the current sample directory
      R1_fastqs=( $( find ${sample_directories[$i]} -type f -name "*R1*.fastq.gz") )

      for((j=0; j<${#R1_fastqs[@]}; j++))
      do

      	echo ${R1_fastqs[$j]}
      	R1_input=${R1_fastqs[$j]}

      	base_temp_name=$( basename -- ${R1_fastqs[$j]} )
      	base_temp_name=${base_temp_name/%.fastq.gz/.bam}
      	echo $base_temp_name

      	##### bowtie2 align
      	bowtie2 \
      		-p 20 -X 2000 \
      		--mm --qc-filter --met 1 -t \
                      --sensitive-local \
      		--no-mixed --no-discordant \
       	        -x $ref_fasta \
                      -1 $R1_input \
      		-2 ${R1_input/_R1_/_R3_} \
                      | samtools sort \
       	        -l 9 -m 5G -@ 16 \
                      -o $temp_dir/"to_merge"/$base_temp_name
      done

      ##### merge and sort all bam files under the to_merge dir
       samtools merge \
              -l 5 -@ 16 -f \
              --no-PG \
              $temp_dir/merged.bam \
             $temp_dir/"to_merge"/*.bam

       samtools sort -l 9 -m 8G -@ 16 \
               -o $temp_dir/merged_sorted.bam \
               $temp_dir/merged.bam
       samtools index -@ 20 $temp_dir/merged_sorted.bam

      # Pull out reads that need to be remapped to check for bias
      # Use the -p option for paired-end reads.

       python $WASP/mapping/find_intersecting_snps.py \
               --is_paired_end \
               --is_sorted \
               --output_dir $DATA_DIR/ \
               --snp_dir $SNP_DIR \
               $temp_dir/merged_sorted.bam

      # Remap the reads, using same the program and options as before.
       echo "map for paired ends"
       bowtie2 -p 20 -X 2000 \
      	--mm --qc-filter --met 1 -t \
               --sensitive-local \
      	--no-mixed --no-discordant \
               -x $ref_fasta \
               -1 $DATA_DIR/"merged_sorted.remap.fq1.gz" \
               -2 $DATA_DIR/"merged_sorted.remap.fq2.gz" \
               | samtools sort \
      	-l 9 -@ 20 -m 5G \
               -o $DATA_DIR/"merged_sorted.remap_paired.bam"
      samtools index -@ 20 $DATA_DIR/"merged_sorted.remap_paired.bam"

      # remap to same position
       python $WASP/mapping/filter_remapped_reads.py \
               $DATA_DIR/"merged_sorted.to.remap.bam" \
               $DATA_DIR/"merged_sorted.remap_paired.bam" \
               $DATA_DIR/"merged_sorted.remap_paired.keep.bam"

      # create the cleaned up bam by merging
       echo "merge bam"
       samtools merge \
	       -l 9 -@ 16 -f \
                $DATA_DIR/"sim_pe_reads.keep.merged.bam" \
                $DATA_DIR/"merged_sorted.remap_paired.keep.bam" \
                $DATA_DIR/"merged_sorted.keep.bam"
	
	samtools sort \
                -l 9 -@ 16 -m 5G \
                -o $DATA_DIR/"sim_pe_reads.keep.merged.sorted.bam" \
		$DATA_DIR/"sim_pe_reads.keep.merged.bam"

	################### Re-add readgroup
	final_output_file_name_prefix=$current_sample_name
	rg_string=$final_output_file_name_prefix

	samtools addreplacerg \
		-w \
		-r "ID:$rg_string\tLB:lib1\tPL:illumina\tSM:$rg_string\tPU:unit1" \
		-@ 16 \
		$DATA_DIR/"sim_pe_reads.keep.merged.sorted.bam" \
		| samtools sort \
		-l 9 -m 5G -@ 16 \
		-o $DATA_DIR/"sim_pe_reads.keep.merged.sorted.rg.bam"

	################### MarkDuplicate
        java -jar $picard MarkDuplicates \
                -I $DATA_DIR/"sim_pe_reads.keep.merged.sorted.rg.bam" \
                -O $DATA_DIR/"sim_pe_reads.keep.merged.sorted.rg.dupmarked.bam" \
		-M metrics.txt \
                --REMOVE_DUPLICATES false
	samtools index -@ 20 $DATA_DIR/"sim_pe_reads.keep.merged.sorted.rg.dupmarked.bam"

	### Dedup by WASP rmdup
        python $WASP/mapping/rmdup_pe.py \
               $DATA_DIR/"sim_pe_reads.keep.merged.sorted.rg.dupmarked.bam" \
               $DATA_DIR/"sim_pe_reads.keep.merged.sorted.rg.dupmarked.removed.bam"
       
       # final sort and output
       samtools sort \
               -l 9 -@ 16 -m 5G \
               -o $OUTPUT_FINAL/$final_output_file_name_prefix'_bulk_known_genotypes_WASPed.bam' \
               $DATA_DIR/"sim_pe_reads.keep.merged.sorted.rg.dupmarked.removed.bam"

	rm -r $temp_dir

done


