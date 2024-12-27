#!/bin/sh

## 5 Jul 2017 Siwei
## modified on 2 Oct 2017 

## This script takes one vcf input file and generate WASP-compatible txt inputs of chr 1-22 (gz compressed)

INPUT_VCF=$1
OUTPUT_DIR=$2

k=1
char_chr='chr'
rm -r $OUTPUT_DIR

mkdir -p $OUTPUT_DIR

while [ $k -lt 23 ]
do
	echo $char_chr$k
	k_chr=$char_chr$k
	echo '$char_chr='$char_chr
	echo '$k_chr='$k_chr
	echo '$k='$k
	OUTPUT_FILE=$OUTPUT_DIR/$k_chr.snps.txt.gz
	cat $INPUT_VCF | awk -F '\t' -v awk_k=$k_chr '$1 == awk_k {print $2"\t"$4"\t"$5}' | gzip > $OUTPUT_FILE
	k=$[$k+1]
done

