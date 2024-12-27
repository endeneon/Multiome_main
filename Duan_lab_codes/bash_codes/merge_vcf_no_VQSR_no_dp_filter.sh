#!/bin/bash


# 24 May 2022
# 07 May 2021 Siwei

## use conda aligners environment

## Set these environment vars to point to
## your local installation

gatk37_path="/home/zhangs3/Data/Tools/GATK37"
gatk36_path="/home/zhangs3/Data/Tools/gatk_36/opt/gatk-3.6"
gatk4_path="/home/zhangs3/Data/Tools/gatk-4.1.8.1"

jre_8="/home/zhangs3/Data/Tools/jre1.8.0_291/bin"
openjdk_16="/home/zhangs3/Data/Tools/jdk-16.0.1/bin/java"

hg38_index="/home/zhangs3/Data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa"

rm -r output/

mkdir -p output

## list all .vcf files
#shopt -s nullglob
#vcf_array=(*.vcf)
#echo "${vcf_array[@]}"

## send all .vcf files to a .list file
ls *.vcf > vcf_list.list

## use picard to merge all vcfs
## note that picard is located under the gatk4 folder
#$openjdk_16 -Xmx200G -jar $gatk4_path/picard2_25_4.jar MergeVcfs \
#	-I vcf_list.list \
#	-O raw_merged_vcfs.vcf

## !! Note!! Use GATK3.6 and jre 1.8 CombineVariants to merge all VCFs
## The new Picard from GATK4 only accepts GVCF files therefore 
## cannot merge vcfs with different sample identities

## detect cell type
cell_type_hint="$(cat vcf_list.list \
	| head -n 1 \
	| sed 's/_/\t/g' \
	| sed 's/hg38only//g' \
	| cut -f 1-5 \
	| awk 'BEGIN { OFS = "_" } {print $0}' \
	| sed 's/\s/_/g')"
echo $cell_type_hint

if [[ $cell_type_hint == *"GABA_bwa"* ]]
then
	output_cell_type="GABA_"
elif [[ $cell_type_hint == *"NEFM_pos"* ]]
then
	output_cell_type="NEFM_pos_"
elif [[ $cell_type_hint == *"NEFM_neg"* ]]
then
	output_cell_type="NEFM_neg_"
else
	echo "Unknown cell type found!"
	exit 1
fi

echo $output_cell_type

# exit 1
if [[ $cell_type_hint == "0hr"* ]]
then
	output_cell_time="0hr"
elif [[ $cell_type_hint == "1hr"* ]]
then
	output_cell_time="1hr"
elif [[ $cell_type_hint == "6hr"* ]]
then
	output_cell_time="6hr"
else
	echo "Unknown cell timing found!"
	exit 1
fi

echo $output_cell_time

# exit 1
# output_cell_type="NEFM_pos_glut"

$jre_8/java -Xmx200G -jar $gatk36_path/GenomeAnalysisTK.jar \
	-T CombineVariants \
	-R $hg38_index \
	--variant vcf_list.list \
	-genotypeMergeOptions UNIQUIFY \
	-o output/raw_merged_vcfs.vcf

## select variants that only includes SNP
$gatk4_path/gatk SelectVariants \
	-R $hg38_index \
	-V output/raw_merged_vcfs.vcf \
	--select-type-to-include SNP \
	--restrict-alleles-to BIALLELIC \
	-XL /home/zhangs3/Data/Databases/Genomes/hg38/hg38_blacklisted_regions.bed \
	-O output/$output_cell_type$output_cell_time"_merged_SNP.vcf"

## select variants with DP >= 2
$gatk4_path/gatk VariantFiltration \
        -R $hg38_index \
        -V output/$output_cell_type$output_cell_time"_merged_SNP.vcf" \
        --filter-name "DP_filter" -filter "DP < 1.1" \
        -O output/$output_cell_type$output_cell_time"_merged_SNP_VF.vcf"

## export vcfs to table
$gatk4_path/gatk VariantsToTable \
	-V output/$output_cell_type$output_cell_time"_merged_SNP_VF.vcf" \
	-F CHROM -F POS -F ID -F REF -F ALT -F NCALLED \
	-GF AD \
	-O output/$output_cell_type$output_cell_time"_merged_SNP_VTT.txt"

## convert the txt title line
cat output/$output_cell_type$output_cell_time"_merged_SNP_VTT.txt" \
	| sed 's/NA/0,0/g' \
	| sed 's/,/\t/g' \
        | sed 's/variant\.AD/variant\tAD/g' \
        | sed 's/variant.\.AD/variant\tAD/g' \
        | sed 's/variant..\.AD/variant\tAD/g' \
	> output/$output_cell_type$output_cell_time"_scATAC_BQSR_16Dec2022_4_R.txt"

