#! /bin/bash
# Siwei 05 Sept 2023

# mkdir -p "../../../interval_bed"

for eachfile in *.vcf
do
	vcf_sample_name=$( echo $(basename -- $eachfile) \
			| sed 's/vcf_4_WASP_//g' \
			| sed 's/.vcf//g' \
			| sed 's/-/_/g' )
	vcf_sample_name=$vcf_sample_name"_het_only_from_genotyping.bed"
	echo $vcf_sample_name

#	echo $eachfile
	cat $eachfile \
		| grep -v "^#" \
		| grep "^chr" \
		| awk -F '\t' '{print $1"\t"$2-1"\t"$2"\t."}' \
		> "./"$vcf_sample_name
done
