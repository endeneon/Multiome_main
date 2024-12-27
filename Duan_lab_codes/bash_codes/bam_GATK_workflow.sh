#!/bin/bash

# Siwei 14 Nov 2022

## !! make softlink of BED4 intervals in ./ or just use $bed4_path variable

WASP="/home/zhangs3/Data/Tools/WASP"
gatk4="/home/zhangs3/Data/Tools/gatk-4.2.6.1/gatk"
samtools="/home/zhangs3/Data/Anaconda3-envs/aligners/bin/samtools"

bwa_ref="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"

bed4_path="/home/zhangs3/Data/FASTQ/Duan_Project_024/hybrid_output/per_library_vcf/cleaned_vcf/output/het_vcf/bed4"
SNP_prefix="/home/zhangs3/Data/FASTQ/Duan_Project_024/hybrid_output/per_library_vcf/cleaned_vcf/output/het_vcf/output"

ref_path="/home/zhangs3/Data/Databases/Genomes/hg38"
ref_genome="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"

bam_output_dir="GATK_realigned_bams"
mkdir -p $bam_output_dir
vcf_output_dir="vcf_output"
mkdir -p $vcf_output_dir

temp_dir="/home/zhangs3/NVME/package_temp/GATK_calib_temp_0hr_pt1"
DATA_DIR=$temp_dir
# mkdir $temp_dir


eachfile="0hr_CD_52_NEFM_neg_gluthg38onlyintersected.bam"

	mkdir -p $DATA_DIR
	WASP_run_dir=$DATA_DIR/'WASP_run'
#	mkdir -p $WASP_run_dir

	samtools index -@ 20 $eachfile

        ############## locate interval file #######
        interval=$(echo $eachfile | \
                sed 's/.*\(CD_[0-9][0-9]\).*/\1/g')
        interval_file="interval_bed"/$interval"_het_only_from_genotyping.bed"
        echo $interval_file

        SNP_DIR=$(echo $eachfile | \
                sed 's/.\+\(CD_[0-9][0-9]\).*/\1/g')
        SNP_DIR=WASP_calib_known_genotypes/$SNP_DIR"_WASP_SNP_fr_genotyping"
        echo $SNP_DIR


	######## make a header file for later use
        # need to extract sam header first
        $samtools view -H \
                $eachfile \
                > $DATA_DIR/new_header.sam

        ## run GATK HaplotypeCaller
       $gatk4 \
                HaplotypeCaller \
                --native-pair-hmm-threads 16 \
                -R $ref_genome \
                -I $eachfile \
                -L $interval_file \
                -XL $ref_path/hg38_blacklisted_regions.bed \
                -O $DATA_DIR/"gatk_vcf_output.vcf" \
                -bamout $DATA_DIR/"GATK_writeout.bam" \
                -DF NotDuplicateReadFilter

       ####### get reads that are not intersected with interval_file sets
       $samtools view \
               -b -h -@ 16 \
               --write-index \
               -L $interval_file \
               -U $DATA_DIR/"non_intersected_reads.bam" \
               -o $DATA_DIR/"intersected_reads.bam" \
               $eachfile

       ####### merge the writeout file and non-intersected_reads file
       ####### use the file header directly from $eachfile
       $samtools merge \
	       -l 9 -@ 16 -f \
	       -h $eachfile \
	       --no-PG \
	       -o $DATA_DIR/"merged_GATK_for_WASP.bam" \
	       $DATA_DIR/"GATK_writeout.bam" \
	       $DATA_DIR/"non_intersected_reads.bam"

       ####### filter and sort bam
       $samtools sort \
	       -l 4 -m 5G -@ 16 \
	       -o $DATA_DIR/"merged_GATK_for_WASP_sorted.bam" \
	       $DATA_DIR/"merged_GATK_for_WASP.bam"
	$samtools index -@ 16 $DATA_DIR/"merged_GATK_for_WASP_sorted.bam"

        ###### Pull out reads that need to be remapped to check for bias
        # Use the -p option for paired-end reads.
        #
        python $WASP/mapping/find_intersecting_snps.py \
                --is_paired_end \
                --is_sorted \
                --output_dir $WASP_run_dir/ \
                --snp_dir $SNP_DIR \
                $DATA_DIR/"merged_GATK_for_WASP_sorted.bam"

        ##### Remap the reads, using same the program and options as before.
	# !! the output bams has to go a second round of GATK before sent for 
	# !! processing with filter_remapped_reads.py
        echo "map for paired ends"
        bwa mem \
                -I 250.0,150.0 \
                -t 16 \
                $bwa_ref \
                $WASP_run_dir/"merged_GATK_for_WASP_sorted.remap.fq1.gz" \
                $WASP_run_dir/"merged_GATK_for_WASP_sorted.remap.fq2.gz" \
                | $samtools sort \
		-l 4 -m 5G -@ 16 \
                -o $WASP_run_dir/"remapped_bam_to_GATK.bam"
        $samtools index $WASP_run_dir/"remapped_bam_to_GATK.bam"

	##### Send the bam file to GATK for 2nd round alignments
        $gatk4 \
                HaplotypeCaller \
                --native-pair-hmm-threads 16 \
                -R $ref_genome \
                -I $eachfile \
                -L $interval_file \
                -XL $ref_path/hg38_blacklisted_regions.bed \
                -O $WASP_run_dir/"gatk_vcf_output.vcf" \
                -bamout $WASP_run_dir/"GATK_writeout.bam"


       ####### get reads that are not intersected with interval_file sets
       $samtools view \
               -b -h -@ 16 \
               --write-index \
               -L $interval_file \
               -U $WASP_run_dir/"non_intersected_reads.bam" \
               -o $WASP_run_dir/"intersected_reads.bam" \
               $eachfile

       ####### merge the writeout file and non-intersected_reads file
       ####### use the file header directly from $eachfile
       $samtools merge \
	       -l 9 -@ 16 -f \
	       -h $eachfile \
	       --no-PG \
	       -o $WASP_run_dir/"merged_GATK_for_WASP.bam" \
	       $WASP_run_dir/"GATK_writeout.bam" \
	       $WASP_run_dir/"non_intersected_reads.bam"

       ####### filter and sort bam
       $samtools sort \
	       -l 4 -m 5G -@ 16 \
	       -o $WASP_run_dir/"merged_GATK_for_WASP_sorted.remap_paired.bam" \
	       $WASP_run_dir/"merged_GATK_for_WASP.bam"
	$samtools index -@ 16 $WASP_run_dir/"merged_GATK_for_WASP_sorted.remap_paired.bam"

#        $samtools index $WASP_run_dir/${eachfile/%.bam/.remap_paired.bam}

        ## remap to same position
        python $WASP/mapping/filter_remapped_reads.py \
                $WASP_run_dir/"merged_GATK_for_WASP_sorted.to.remap.bam" \
                $WASP_run_dir/"merged_GATK_for_WASP_sorted.remap_paired.bam" \
                $WASP_run_dir/"merged_GATK_for_WASP_sorted.remap_paired.keep.bam"

        ## Create a merged BAM containing [1] reads (PE+SE) that did
        ## not need remapping [2] filtered remapped reads
        #       echo "merge bam"

       $samtools merge \
               -l 5 -@ 16 -f \
               --no-PG \
               $WASP_run_dir/"merged_GATK_for_WASP_sorted.sim_pe_reads.merged.bam" \
               $WASP_run_dir/"merged_GATK_for_WASP_sorted.remap_paired.keep.bam" \
               $WASP_run_dir/"merged_GATK_for_WASP_sorted.keep.bam"

       $samtools view \
                -h \
                $WASP_run_dir/"merged_GATK_for_WASP_sorted.sim_pe_reads.merged.bam" \
               | $samtools sort \
	       -l 9 -m 5G -@ 16 \
               -o $bam_output_dir/${eachfile/%.bam/_GATK_calib.bam}


